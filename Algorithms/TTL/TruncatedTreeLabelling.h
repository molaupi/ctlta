#pragma once

#include "Algorithms/TTL/TopologyCentricTreeHierarchy.h"

class TruncatedTreeLabelling {

public:
    TruncatedTreeLabelling(const TopologyCentricTreeHierarchy &hierarchy)
            : hierarchy(hierarchy), upLabelData(), downLabelData() {}

    // Initializes tree labelling with underlying tree hierarchy. (Make sure to preprocess tree hierarchy before calling).
    void init() {
        uint64_t offset = 0;
        const auto &numHubs = hierarchy.getNumHubs();
        labelOffsets.resize(numHubs.size());
        for (auto i = 0; i < numHubs.size(); ++i) {
            labelOffsets[i] = offset;
            offset += 3 * numHubs[i]; // numHubs entries for distances and numHubs entries for path edge pointers and TODO: debug: rank in order of hub
        }
        upLabelData.resize(offset, INFTY);
        downLabelData.resize(offset, INFTY);
    }

    void reset() {
        std::fill(upLabelData.begin(), upLabelData.end(), INFTY);
        std::fill(downLabelData.begin(), downLabelData.end(), INFTY);
    }

    inline int32_t &upDist(const int32_t &v, const uint32_t &hubIdx) {
        return upLabelData[labelOffsets[v] + hubIdx];
    }

    const int32_t &upDist(const int32_t &v, const uint32_t &hubIdx) const {
        return upLabelData[labelOffsets[v] + hubIdx];
    }

    inline int32_t &upPathEdge(const int32_t &v, const uint32_t &hubIdx) {
        return upLabelData[labelOffsets[v] + hierarchy.getNumHubs()[v] + hubIdx];
    }

    const int32_t &upPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        return upLabelData[labelOffsets[v] + hierarchy.getNumHubs()[v] + hubIdx];
    }

    inline int32_t &upHub(const int32_t &v, const uint32_t &hubIdx) {
        return upLabelData[labelOffsets[v] + 2 * hierarchy.getNumHubs()[v] + hubIdx];
    }

    const  int32_t &upHub(const int32_t &v, const uint32_t &hubIdx) const {
        return upLabelData[labelOffsets[v] + 2 * hierarchy.getNumHubs()[v] + hubIdx];
    }

    inline int32_t &downDist(const int32_t &v, const uint32_t &hubIdx) {
        return downLabelData[labelOffsets[v] + hubIdx];
    }

    const int32_t &downDist(const int32_t &v, const uint32_t &hubIdx) const {
        return downLabelData[labelOffsets[v] + hubIdx];
    }

    inline int32_t &downPathEdge(const int32_t &v, const uint32_t &hubIdx) {
        return downLabelData[labelOffsets[v] + hierarchy.getNumHubs()[v] + hubIdx];
    }

    const int32_t &downPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        return downLabelData[labelOffsets[v] + hierarchy.getNumHubs()[v] + hubIdx];
    }

    inline int32_t &downHub(const int32_t &v, const uint32_t &hubIdx) {
        return downLabelData[labelOffsets[v] + 2 * hierarchy.getNumHubs()[v] + hubIdx];
    }

    const  int32_t &downHub(const int32_t &v, const uint32_t &hubIdx) const {
        return downLabelData[labelOffsets[v] + 2 * hierarchy.getNumHubs()[v] + hubIdx];
    }

    void assertFullyCustomized() const {
        for (int32_t v = 0; v < labelOffsets.size(); ++v) {
            for (auto i = 0; i < hierarchy.getNumHubs()[v]; ++i) {
                KASSERT(!(upDist(v, i) == INFTY || upPathEdge(v, i) == INFTY || downDist(v, i) == INFTY ||
                    downPathEdge(v, i) == INFTY));
            }
        }
    }

private:

    const TopologyCentricTreeHierarchy &hierarchy;

    std::vector<uint32_t> labelOffsets;
    std::vector<int32_t> upLabelData; // expects distances, and edge IDs to be int32_t.
    std::vector<int32_t> downLabelData; // expects distances, and edge IDs to be int32_t.

};