#pragma once

#include "Algorithms/TTL/BalancedTopologyCentricTreeHierarchy.h"

#ifndef TTL_STORE_PATH_POINTERS
#define TTL_STORE_PATH_POINTERS true
#endif

template<bool StorePathPointers = TTL_STORE_PATH_POINTERS>
class TruncatedTreeLabelling {

    static constexpr uint64_t INVALID_OFFSET = static_cast<uint64_t>(-1);

public:

    static const bool StoresPathPointers = StorePathPointers;

    struct ConstLabel {

        const int32_t &dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[hubIdx];
        }

        template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
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

        template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
        const int32_t &pathEdge(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[numHubs + hubIdx];
        }

        template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
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
            offset += hierarchy.getNumHubs(v); // one distance entry per hub
            if constexpr (StorePathPointers)
                offset += hierarchy.getNumHubs(v); // one path edge per hub
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
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return upLabelData[labelOffsets[v] + hubIdx];
    }

    const int32_t &upDist(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return upLabelData[labelOffsets[v] + hubIdx];
    }

    template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
    inline int32_t &upPathEdge(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return upLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
    const int32_t &upPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return upLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    inline int32_t &downDist(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return downLabelData[labelOffsets[v] + hubIdx];
    }

    const int32_t &downDist(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return downLabelData[labelOffsets[v] + hubIdx];
    }

    template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
    inline int32_t &downPathEdge(const int32_t &v, const uint32_t &hubIdx) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        return downLabelData[labelOffsets[v] + hierarchy.getNumHubs(v) + hubIdx];
    }

    template<bool hasPathEdges = StorePathPointers, std::enable_if_t<hasPathEdges>...>
    const int32_t &downPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
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
        return labelOffsets.size() * sizeof(typename decltype(labelOffsets)::value_type)
               + upLabelData.size() * sizeof(typename decltype(upLabelData)::value_type)
               + downLabelData.size() * sizeof(typename decltype(downLabelData)::value_type);
    }

private:

    const BalancedTopologyCentricTreeHierarchy &hierarchy;

    std::vector<uint64_t> labelOffsets;
    AlignedVector<int32_t> upLabelData; // expects distances, and edge IDs to be int32_t.
    AlignedVector<int32_t> downLabelData; // expects distances, and edge IDs to be int32_t.

};