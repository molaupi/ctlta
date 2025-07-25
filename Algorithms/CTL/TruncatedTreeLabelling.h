#pragma once

#include "Algorithms/CTL/BalancedTopologyCentricTreeHierarchy.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/SimdLabelSet.h"

template<typename LabelSetT>
class TruncatedTreeLabelling {

    static constexpr uint64_t INVALID_OFFSET = static_cast<uint64_t>(-1);

public:

    using LabelSet = LabelSetT;
    static constexpr uint32_t K = LabelSet::K;
    using Batch = typename LabelSet::DistanceLabel; // A batch of K distances or path edges
    using BatchMask = typename LabelSet::LabelMask; // A mask representing K boolean flags

    static_assert(LabelSet::logK != 1, "logK cannot be 1. Use 0 for no SIMD or >= 2 for SIMD.");

    struct ConstBatchLabel {

        ConstBatchLabel(Batch const *startOfLabel, uint32_t numHubs)
                : startOfLabel(startOfLabel), numHubs(numHubs), numBatches(numHubs / K + (numHubs % K != 0)) {}

        int32_t dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[hubIdx / K][hubIdx % K];
        }

        const Batch &distBatch(const uint32_t &batchIdx) const {
            KASSERT(batchIdx < numBatches);
            return startOfLabel[batchIdx];
        }

        template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
        int32_t pathEdge(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs);
            return startOfLabel[numBatches + hubIdx / K][hubIdx % K];
        }

        template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
        const Batch &pathEdgeBatch(const uint32_t &batchIdx) const {
            KASSERT(batchIdx < numBatches);
            return startOfLabel[numBatches + batchIdx];
        }

        Batch const *startOfLabel;
        uint32_t numHubs;
        uint32_t numBatches;
    };

    struct BatchLabel {

        BatchLabel(Batch *startOfLabel, uint32_t numHubs)
                : startOfLabel(startOfLabel), numHubs(numHubs), numBatches(numHubs / K + (numHubs % K != 0)) {}

        // Sets distance at last hub in label to 0 and all unused elements of last batch to INFTY.
        // Everything else is left untouched.
        void initializeLastHubDist() {
            startOfLabel[(numHubs - 1) / K].setUpperElements((numHubs - 1) % K, 0);
            if (numHubs % K != 0)
                startOfLabel[(numHubs - 1) / K].setUpperElements((numHubs) % K, INFTY);
        }

        // Sets path edge at last hub in label to INVALID_EDGE.
        template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
        void initializeLastHubPathEdge() {
            // The call to setUpperElements also sets any unused elements of the last batch to INVALID_EDGE. This
            // should never have an effect as the unused elements of the last batch should never be used for
            // any computations.
            startOfLabel[numBatches + (numHubs - 1) / K].setUpperElements((numHubs - 1) % K, INVALID_EDGE);

        }

        const Batch &distBatch(const uint32_t &batchIdx) const {
            KASSERT(batchIdx < numBatches);
            return startOfLabel[batchIdx];
        }

        inline Batch &distBatch(const uint32_t &batchIdx) {
            KASSERT(batchIdx < numBatches);
            return startOfLabel[batchIdx];
        }

        template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
        const Batch &pathEdgeBatch(const uint32_t &batchIdx) const {
            KASSERT(batchIdx < numBatches);
            return startOfLabel[numBatches + batchIdx];
        }

        template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES, std::enable_if_t<hasPathEdges, bool> = true>
        inline Batch &pathEdgeBatch(const uint32_t &batchIdx) {
            KASSERT(batchIdx < numBatches);
            return startOfLabel[numBatches + batchIdx];
        }

        Batch *startOfLabel;
        uint32_t numHubs;
        uint32_t numBatches;
    };

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
            const auto numBatches = numHubs / K + (numHubs % K != 0);
            offset += numBatches; // one distance entry per hub, K entries per vec
            if constexpr (LabelSet::KEEP_PARENT_EDGES)
                offset += numBatches; // one path edge per hub, K edges per vec
        }
        upLabelData.resize(offset, Batch(INFTY));
        downLabelData.resize(offset, Batch(INFTY));
    }

    void reset() {
        std::fill(upLabelData.begin(), upLabelData.end(), Batch(INFTY));
        std::fill(downLabelData.begin(), downLabelData.end(), Batch(INFTY));
    }

    ConstBatchLabel upLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstBatchLabel(upLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    ConstBatchLabel cUpLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstBatchLabel(upLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    BatchLabel upLabel(const int32_t &v) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return BatchLabel(upLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    ConstBatchLabel downLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstBatchLabel(downLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    ConstBatchLabel cDownLabel(const int32_t &v) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return ConstBatchLabel(downLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    BatchLabel downLabel(const int32_t &v) {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        return BatchLabel(downLabelData.data() + labelOffsets[v], hierarchy.getNumHubs(v));
    }

    // Convenience method to directly access an up path edge of a vertex without explicitly getting the vertex's label.
    template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES,
            std::enable_if_t<hasPathEdges, bool> = true>
    int32_t upPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        const auto &numHubs = hierarchy.getNumHubs(v);
        const auto numBatches = numHubs / K + (numHubs % K != 0);
        return upLabelData[labelOffsets[v] + numBatches + hubIdx / K][hubIdx % K];
    }

    // Convenience method to directly access a down path edge of a vertex without explicitly getting the vertex's label.
    template<bool hasPathEdges = LabelSet::KEEP_PARENT_EDGES,
            std::enable_if_t<hasPathEdges, bool> = true>
    int32_t downPathEdge(const int32_t &v, const uint32_t &hubIdx) const {
        KASSERT(labelOffsets[v] != INVALID_OFFSET);
        KASSERT(hubIdx < hierarchy.getNumHubs(v));
        const auto &numHubs = hierarchy.getNumHubs(v);
        const auto numBatches = numHubs / K + (numHubs % K != 0);
        return downLabelData[labelOffsets[v] + numBatches + hubIdx / K][hubIdx % K];
    }

    uint64_t sizeInBytes() const {
        return labelOffsets.size() * sizeof(typename decltype(labelOffsets)::value_type)
               + upLabelData.size() * sizeof(typename decltype(upLabelData)::value_type)
               + downLabelData.size() * sizeof(typename decltype(downLabelData)::value_type);
    }

private:

    const BalancedTopologyCentricTreeHierarchy &hierarchy;

    std::vector<uint64_t> labelOffsets;
    AlignedVector<Batch> upLabelData; // expects distances, and edge IDs to be int32_t.
    AlignedVector<Batch> downLabelData; // expects distances, and edge IDs to be int32_t.

};