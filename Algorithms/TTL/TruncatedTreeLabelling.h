#pragma once

#include "Algorithms/TTL/BalancedTopologyCentricTreeHierarchy.h"

class TruncatedTreeLabelling {

    static constexpr uint32_t INVALID_OFFSET = static_cast<uint32_t>(-1);

public:

    struct ConstLabel {

        const int32_t &dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[hubIdx];
        }

        const int32_t &pathEdge(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[numHubs + hubIdx];
        }

        int32_t const *startOfLabel;
        uint32_t numHubs;
    };

    struct Label {
        const int32_t &dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[hubIdx];
        }

        inline int32_t &dist(const uint32_t &hubIdx) {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[hubIdx];
        }

        const int32_t &pathEdge(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[numHubs + hubIdx];
        }

        inline int32_t &pathEdge(const uint32_t &hubIdx) {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[numHubs + hubIdx];
        }

        int32_t *startOfLabel;
        uint32_t numHubs;
    };


    TruncatedTreeLabelling(const BalancedTopologyCentricTreeHierarchy &hierarchy)
            : hierarchy(hierarchy), upLabelData(), downLabelData() {}

    // Initializes tree labelling with underlying tree hierarchy. (Make sure to preprocess tree hierarchy before calling).
    void init() {
        uint64_t offset = 0;
        labelOffsets.resize(hierarchy.numVertices(), INVALID_OFFSET);
        for (auto v = 0; v < hierarchy.numVertices(); ++v) {
            // Do not initialize labels for truncated vertices.
            if (hierarchy.isVertexTruncated(v))
                continue;
            labelOffsets[v] = offset;
            offset += 2 * hierarchy.getNumHubs(v); // numHubs entries for distances and numHubs entries for path edge pointers
        }
        upLabelData.resize(offset, INFTY);
        downLabelData.resize(offset, INFTY);
    }

    void reset() {
        std::fill(upLabelData.begin(), upLabelData.end(), INFTY);
        std::fill(downLabelData.begin(), downLabelData.end(), INFTY);
    }

    ConstLabel upLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstLabel(upLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    ConstLabel cUpLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstLabel(upLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    Label upLabel(const int32_t &v) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return Label(upLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    ConstLabel downLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstLabel(downLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    ConstLabel cDownLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstLabel(downLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    Label downLabel(const int32_t &v) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return Label(downLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }


    inline int32_t &upDist(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return upLabelData[labelOffsets[v] + hubIdx];
    }

    const int32_t &upDist(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return upLabelData[labelOffsets[v] + hubIdx];
    }

    inline int32_t &upPathEdge(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return upLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    const int32_t &upPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return upLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    inline int32_t &downDist(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return downLabelData[labelOffsets[v] + hubIdx];
    }

    const int32_t &downDist(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return downLabelData[labelOffsets[v] + hubIdx];
    }

    inline int32_t &downPathEdge(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return downLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    const int32_t &downPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return downLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    void assertFullyCustomized() const {
        for (int32_t v = 0; v < labelOffsets.size(); ++v) {
            if (hierarchy.isVertexTruncated(v))
                continue;
            for (auto i = 0; i < hierarchy.getNumHubs(v); ++i) {
                KASSERT(!(upDist(v, i) == INFTY || upPathEdge(v, i) == INFTY || downDist(v, i) == INFTY ||
                          downPathEdge(v, i) == INFTY));
            }
        }
    }

    uint64_t sizeInBytes() const {
        return labelOffsets.size() * sizeof(decltype(labelOffsets)::value_type)
               + upLabelData.size() * sizeof(decltype(upLabelData)::value_type)
               + downLabelData.size() * sizeof(decltype(downLabelData)::value_type);
    }

private:

    const BalancedTopologyCentricTreeHierarchy &hierarchy;

    std::vector<uint32_t> labelOffsets;
    std::vector<int32_t> upLabelData; // expects distances, and edge IDs to be int32_t.
    std::vector<int32_t> downLabelData; // expects distances, and edge IDs to be int32_t.

};