#pragma once

#include <vector>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <barrier>
#include <cassert>

#include "road_network.h"

namespace chl::util {

// start new time measurement
void start_timer();
// returns time in seconds since last unconsumed start_timer call and consumes it
double stop_timer();

// sort vector and remove duplicate elements
template<typename T>
void make_set(std::vector<T> &v)
{
    size_t v_size = v.size();
    if (v_size == 0)
        return;
    std::sort(v.begin(), v.end());
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next] != v[last_distinct])
        {
            last_distinct++;
            std::swap(v[next], v[last_distinct]);
        }
    v.resize(last_distinct + 1);
}

template<typename T, class Compare>
void make_set_(std::vector<T> &v, Compare comp)
{
    size_t v_size = v.size();
    if (v_size < 2)
        return;
    std::sort(v.begin(), v.end(), comp);
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next] != v[last_distinct])
        {
            last_distinct++;
            std::swap(v[next], v[last_distinct]);
        }
    v.erase(v.begin() + (last_distinct + 1), v.end());
}

template<typename T, class Compare>
void make_set(std::vector<T> &v, Compare comp)
{
    size_t v_size = v.size();
    if (v_size < 2)
        return;
    std::sort(v.begin(), v.end(), comp);
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next].neighbor != v[last_distinct].neighbor)
        {
            last_distinct++;
	    road_network::CHNeighbor temp = v[next];
	    v[next] = v[last_distinct];
	    v[last_distinct] = temp;
        }
    v.erase(v.begin() + (last_distinct + 1), v.end());
}

// remove elements in set from v
// set must be sorted
template<typename T>
void remove_set(std::vector<T> &v, const std::vector<T> set)
{
    assert(is_sorted(set.cbegin(), set.cend()));
    if (v.empty() || set.empty())
        return;
    std::erase_if(v, [&set](T value) { return std::binary_search(set.cbegin(), set.cend(), value); });
}

struct Summary
{
    double min;
    double max;
    double avg;
    Summary operator*(double x) const;
};

template<typename T, class Map>
Summary summarize(const std::vector<T> &v, Map f)
{
    Summary summary = {};
    if (!v.empty())
        summary.min = f(v[0]);
    for (const T& e : v)
    {
        double x = f(e);
        summary.avg += x;
        if (x < summary.min)
            summary.min = x;
        if (x > summary.max)
            summary.max = x;
    }
    if (!v.empty())
        summary.avg /= v.size();
    return summary;
}

// compute total number of elements in vector of collections
template<typename T>
size_t size_sum(const std::vector<T> &v)
{
    size_t sum = 0;
    for (const T &x : v)
        sum += x.size();
    return sum;
}

// extract size values of vector of collection
template<typename T>
std::vector<size_t> sizes(const std::vector<T> &v)
{
    std::vector<size_t> s;
    for (const T &x : v)
        s.push_back(x.size());
    return s;
}

template<typename T>
T random(const std::vector<T> &v)
{
    assert(v.size() > 0);
    return v[rand() % v.size()];
}

// list of buckets for synchronized access; buckets can be iterated over in either ascending or descending order, and can be reset for re-use
template<typename T, size_t threads>
class par_bucket_list
{
    std::vector<std::vector<T>> buckets;
    size_t current_bucket = 0, next_in_bucket = 0;
    std::barrier<> sync_point;
    std::mutex m_mutex;
    bool is_empty = false, ascending;
    void on_complete()
    {
        assert(!is_empty);
        size_t last_bucket = (ascending && buckets.size() > 0) ? buckets.size() - 1 : 0;
        if (current_bucket == last_bucket)
        {
            is_empty = true;
            return;
        }
        next_in_bucket = 0;
        // advance to next non-empty bucket
        while (buckets[ascending ? ++current_bucket : --current_bucket].empty())
            if (current_bucket == last_bucket)
            {
                is_empty = true;
                return;
            }
    }
public:
    par_bucket_list(size_t max_bucket, bool ascending) : sync_point(threads), ascending(ascending)
    {
        // prevent re-allocation which could cause syncronization issues
        buckets.reserve(max_bucket + 1);
	ascending ? current_bucket = 0 : current_bucket = max_bucket;
    }
    void push(T value, size_t bucket)
    {
        if (buckets.size() <= bucket)
        {
            buckets.resize(bucket + 1);
        }
        buckets[bucket].push_back(value);
    }
    void print() {
	std::cout << current_bucket << std::endl;
        for(std::vector<uint32_t> vec: buckets) {
	    
	    for(uint32_t val: vec)
		std::cout << val << " ";
	    std::cout << std::endl;
	}
    }
    bool next(T& value, size_t thread)
    {
	
        while (true)
        {
            if (is_empty)
                return false;
            {
                // if there's a value available in current bucket, simply return it
                std::lock_guard<std::mutex> lock(m_mutex);
                if (next_in_bucket < buckets[current_bucket].size()) {
                    value = buckets[current_bucket][next_in_bucket++];
                    return true;
                }
            }
            sync_point.arrive_and_wait();
            // ugly workaround for barrier template fuckery
            if (thread == 0)
                on_complete();
            sync_point.arrive_and_wait();
        }
    }
    void reset(bool ascending)
    {
        this->ascending = ascending;
        current_bucket = (ascending || buckets.size() == 0) ? 0 : buckets.size() - 1;
        next_in_bucket = 0;
        is_empty = false;
    }
};

} // util

namespace std {

enum class ListFormat { plain, indexed };

void set_list_format(ListFormat format);
ListFormat get_list_format();

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
    if (v.empty())
        return os << "[]";
    if (get_list_format() == ListFormat::indexed)
    {
        os << "[0:" << v[0];
        for (size_t i = 1; i < v.size(); i++)
            os << ',' << i << ":" << v[i];
    }
    else
    {
        os << "[" << v[0];
        for (size_t i = 1; i < v.size(); i++)
            os << ',' << v[i];
    }
    return os << ']';
}

template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::pair<A,B> &p)
{
    return os << "(" << p.first << "," << p.second << ")";
}

std::ostream& operator<<(std::ostream& os, chl::util::Summary s);

} // std
