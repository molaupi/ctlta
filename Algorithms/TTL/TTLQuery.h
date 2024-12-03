#pragma once

#include "Algorithms/TTL/TruncatedTreeLabelling.h"

class TTLQuery {

public:
    TTLQuery(const TopologyCentricTreeHierarchy& hierarchy, const TruncatedTreeLabelling& ttl) : hierarchy(hierarchy), ttl(ttl) {}

    // Expects ranks in the underlying separator decomposition order as inputs.
    void run(const int32_t s, const int32_t t) {
        lastDistance = INFTY;
        lastMeetingHubIdx = INVALID_INDEX;
        lastS = s;
        lastT = t;
        for (uint32_t i = 0; i < hierarchy.getLowestCommonHub(s, t); ++i) {
            KASSERT(ttl.upHub(s,i) == INFTY || ttl.downHub(t, i) == INFTY || ttl.upHub(s, i) == ttl.downHub(t, i));
            const int32_t distForHubI = ttl.upDist(s, i) + ttl.downDist(t, i);
            if (distForHubI < lastDistance) {
                lastDistance = distForHubI;
                lastMeetingHubIdx = i;
            }
        }
    }

    // Returns the distance of the last query.
    int32_t getDistance() const {
        return lastDistance;
    }

    // Returns the CCH edges in the upward graph on the up segment of the up-down path (in reverse order to conform to
    // standard).
    const std::vector<int32_t>& getUpEdgePath() {
        lastUpPath.clear();
        auto v = lastS;
        int32_t e;
        while ((e = ttl.upPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            lastUpPath.push_back(e);
        }
        std::reverse(lastUpPath.begin(), lastUpPath.end());
        return lastUpPath;
    }

    // Returns the CCH edges in the upward graph on the down segment of the up-down path.
    const std::vector<int32_t>& getDownEdgePath() {
        lastDownPath.clear();
        auto v = lastT;
        int32_t e;
        while ((e = ttl.downPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            lastDownPath.push_back(e);
        }
        std::reverse(lastDownPath.begin(), lastDownPath.end());
        return lastDownPath;
    }


private:

    int32_t lastS;
    int32_t lastT;
    int32_t lastDistance;
    uint32_t lastMeetingHubIdx;
    std::vector<int32_t> lastUpPath;
    std::vector<int32_t> lastDownPath;

    const TopologyCentricTreeHierarchy& hierarchy;
    const TruncatedTreeLabelling& ttl;

};

