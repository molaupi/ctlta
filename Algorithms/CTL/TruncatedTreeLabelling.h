#pragma once

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/SimdLabelSet.h"

template<int K, bool KEEP_PARENT_EDGES>
class TruncatedTreeLabelling {

    static constexpr uint64_t INVALID_OFFSET = static_cast<uint64_t>(-1);

public:

    struct ConstBatchLabel {

        ConstBatchLabel(int const *const startOfDists, int const *const startOfEdges, const uint32_t numHubs)
                : startOfDists(startOfDists), startOfEdges(startOfEdges), numHubs(numHubs) {}

        int32_t const *startDists() const {
            return startOfDists;
        }

        int32_t const *startEdges() const {
            return startOfEdges;
        }

        const int &dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfDists[hubIdx];
        }

        const int &pathEdge(const uint32_t &hubIdx) const requires KEEP_PARENT_EDGES {
            KASSERT(hubIdx < numHubs);
            return startOfEdges[hubIdx];
        }

        int const *const startOfDists;
        int const *const startOfEdges;
        const uint32_t numHubs;
    };

    struct BatchLabel {

        BatchLabel(int *startOfDists, int *startOfEdges, uint32_t numHubs)
                : startOfDists(startOfDists), startOfEdges(startOfEdges), numHubs(numHubs) {}

        // Sets distance at last hub in label to 0 and all unused elements of last batch to INFTY.
        // Everything else is left untouched.
        void initializeLastHubDist() {
            startOfDists[numHubs - 1] = 0;
            for (int i = numHubs; i % K != 0; ++i) {
                startOfDists[i] = INFTY;
            }
        }

        // Sets path edge at last hub in label to INVALID_EDGE.
        void initializeLastHubPathEdge() requires KEEP_PARENT_EDGES {
            startOfEdges[numHubs - 1] = INVALID_EDGE;

        }

        int *startDists() {
            return startOfDists;
        }

        int const *startDists() const {
            return startOfDists;
        }

        int32_t *startEdges() {
            return startOfEdges;
        }

        int32_t const *startEdges() const {
            return startOfEdges;
        }

        const int &dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfDists[hubIdx];
        }

        inline int &dist(const uint32_t &hubIdx) {
            KASSERT(hubIdx < numHubs);
            return startOfDists[hubIdx];
        }

        const int &pathEdge(const uint32_t &hubIdx) const requires KEEP_PARENT_EDGES {
            KASSERT(hubIdx < numHubs);
            return startOfEdges[hubIdx];
        }

        inline int &pathEdge(const uint32_t &hubIdx) requires KEEP_PARENT_EDGES {
            KASSERT(hubIdx < numHubs);
            return startOfEdges[hubIdx];
        }

        int *startOfDists;
        int *startOfEdges;
        uint32_t numHubs;
    };

    static constexpr int padToNextMultipleOfK(const int &numHubs) {
        return (((numHubs - 1) / K) + 1) * K; // round up to next multiple of K
    }

public:

    explicit TruncatedTreeLabelling(const BalancedTopologyCentricTreeHierarchy &hierarchy)
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
            const auto numHubs = hierarchy.getNumHubs(v);
            // Pad each label to next multiple of K to simplify using SIMD operations.
            const auto paddedNumHubs = padToNextMultipleOfK(numHubs);
            offset += paddedNumHubs; // one distance entry per hub, K entries per vec
            if constexpr (KEEP_PARENT_EDGES)
                offset += paddedNumHubs; // one path edge per hub, K edges per vec
        }
        upLabelData.resize(offset, INFTY);
        downLabelData.resize(offset, INFTY);
    }

    void reset() {
        std::fill(upLabelData.begin(), upLabelData.end(), INFTY);
        std::fill(downLabelData.begin(), downLabelData.end(), INFTY);
    }

    ConstBatchLabel upLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        const int numHubs = hierarchy.getNumHubs(v);
        return ConstBatchLabel(upLabelData.data() + labelOffsets[v],
                               upLabelData.data() + labelOffsets[v] + padToNextMultipleOfK(numHubs),
                               numHubs);
    }

    ConstBatchLabel cUpLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        const int numHubs = hierarchy.getNumHubs(v);
        return ConstBatchLabel(upLabelData.data() + labelOffsets[v],
                               upLabelData.data() + labelOffsets[v] + padToNextMultipleOfK(numHubs),
                               numHubs);
    }

    BatchLabel upLabel(const int32_t &v) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        const int numHubs = hierarchy.getNumHubs(v);
        return BatchLabel(upLabelData.data() + labelOffsets[v],
                          upLabelData.data() + labelOffsets[v] + padToNextMultipleOfK(numHubs),
                          numHubs);
    }

    ConstBatchLabel downLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        const int numHubs = hierarchy.getNumHubs(v);
        return ConstBatchLabel(downLabelData.data() + labelOffsets[v],
                               downLabelData.data() + labelOffsets[v] + padToNextMultipleOfK(numHubs),
                               numHubs);
    }

    ConstBatchLabel cDownLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        const int numHubs = hierarchy.getNumHubs(v);
        return ConstBatchLabel(downLabelData.data() + labelOffsets[v],
                               downLabelData.data() + labelOffsets[v] + padToNextMultipleOfK(numHubs),
                               numHubs);
    }

    BatchLabel downLabel(const int32_t &v) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        const int numHubs = hierarchy.getNumHubs(v);
        return BatchLabel(downLabelData.data() + labelOffsets[v],
                               downLabelData.data() + labelOffsets[v] + padToNextMultipleOfK(numHubs),
                               numHubs);
    }

    // Convenience method to directly access an up path edge of a vertex without explicitly getting the vertex's label.
    int32_t upPathEdge(const int32_t &v, const uint32_t &hubIdx) const requires KEEP_PARENT_EDGES {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        const auto &numHubs = hierarchy.getNumHubs(v);
        return upLabelData[labelOffsets[v] + padToNextMultipleOfK(numHubs) + hubIdx];
    }

    // Convenience method to directly access a down path edge of a vertex without explicitly getting the vertex's label.
    int32_t downPathEdge(const int32_t &v, const uint32_t &hubIdx) const requires KEEP_PARENT_EDGES {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        const auto &numHubs = hierarchy.getNumHubs(v);
        return downLabelData[labelOffsets[v] + padToNextMultipleOfK(numHubs) + hubIdx];
    }

    uint64_t sizeInBytes() const {
        return sizeof(TruncatedTreeLabelling<K, KEEP_PARENT_EDGES>)
               + labelOffsets.size() * sizeof(typename decltype(labelOffsets)::value_type)
               + upLabelData.size() * sizeof(typename decltype(upLabelData)::value_type)
               + downLabelData.size() * sizeof(typename decltype(downLabelData)::value_type);
    }

private:

    const BalancedTopologyCentricTreeHierarchy &hierarchy;

    std::vector<uint64_t> labelOffsets;
    AlignedVector<int32_t> upLabelData; // expects distances, and edge IDs to be int32_t.
    AlignedVector<int32_t> downLabelData; // expects distances, and edge IDs to be int32_t.

};