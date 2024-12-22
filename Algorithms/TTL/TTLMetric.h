#pragma once

#include "Algorithms/TTL/TruncatedTreeLabelling.h"
#include "Algorithms/TTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCHMetric.h"
#include <vectorclass.h>
#include "Tools/Constants.h"

#ifndef TTL_USE_PERFECT_CUSTOMIZATION
#define TTL_USE_PERFECT_CUSTOMIZATION false
#endif

// An individual metric for truncated tree labelling. Uses a CCH to build customized hub labelling.
template<typename LabellingT, bool USE_PERFECT_CUSTOMIZATION = TTL_USE_PERFECT_CUSTOMIZATION>
class TTLMetric {

public:

    using SearchGraph = std::conditional_t<USE_PERFECT_CUSTOMIZATION, CH::SearchGraph, CCH::UpGraph>;

    // Constructs an individual metric incorporating the specified input weights in the specified
    // BalancedTopologyCentricTreeHierarchy on the basis of the specified CCH.
    TTLMetric(const BalancedTopologyCentricTreeHierarchy &hierarchy, const CCH &cch, const int32_t *const inputWeights)
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights), minimumWeightedCH() {}

    void buildCustomizedTTL(LabellingT &ttl) {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            minimumWeightedCH = cchMetric.buildMinimumWeightedCH();
        else
            cchMetric.customize();
        customizeLabelling(ttl);
    }

    // Returns the upward graph of the metric, which is either the CCH graph or the upward graph after perfect
    // customization.
    const SearchGraph &upwardGraph() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return minimumWeightedCH.upwardGraph();
        else
            return cch.getUpwardGraph();
    }

    // Returns the downward graph of the metric, which is either the CCH graph or the downward graph after perfect
    // customization.
    const SearchGraph &downwardGraph() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return minimumWeightedCH.downwardGraph();
        else
            return cch.getUpwardGraph();
    }

    // Returns current weights on edges of upward graph.
    int const *upwardWeights() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return &minimumWeightedCH.upwardGraph().template get<CH::Weight>(0);
        else
            return cchMetric.upWeights.data();
    }

    // Returns current weights on edges of downward graph.
    int const *downwardWeights() const {
        if constexpr (USE_PERFECT_CUSTOMIZATION)
            return &minimumWeightedCH.downwardGraph().template get<CH::Weight>(0);
        else
            return cchMetric.downWeights.data();
    }

private:

    using LabelSet = typename LabellingT::LabelSet;
    static constexpr int K = LabelSet::K;
    using Batch = typename LabelSet::DistanceLabel;
    using BatchMask = typename LabelSet::LabelMask;

    void customizeLabelling(LabellingT &ttl) {
        ttl.reset();

        const auto &upGraph = upwardGraph();
        const auto &downGraph = downwardGraph(); // reverse downward graph
        const auto upWeights = upwardWeights();
        const auto downWeights = downwardWeights();

#pragma omp parallel // parallelizes forEachVertexTopDown
#pragma omp single nowait
        cch.forEachVertexTopDown([&](const int u) {
            // Do not build labels for truncated vertices.
            if (hierarchy.isVertexTruncated(u))
                return;

            BatchMask improved;
            const auto numHubsU = hierarchy.getNumHubs(u);

            // Customize upward label of u using upper neighbors
            auto uUpLabel = ttl.upLabel(u);
            uUpLabel.initializeLastHubDist(); // distance to self
            if constexpr (LabelSet::KEEP_PARENT_EDGES)
                uUpLabel.initializeLastHubPathEdge(); // edge to self
            FORALL_INCIDENT_EDGES(upGraph, u, e) {
                const auto v = upGraph.edgeHead(e);
                const auto upWeight = upWeights[e];
                const auto numHubsV = hierarchy.getNumHubs(v);
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                const auto vUpLabel = ttl.cUpLabel(v);
                const uint32_t numBlocks = numHubsV / K + (numHubsV % K != 0);
                const Batch weightAsBatch(upWeight);
                const Batch eAsBatch(e);
                for (uint32_t b = 0; b < numBlocks; ++b) {
                    const auto dNew = weightAsBatch + vUpLabel.distBatch(b);
                    auto &dU = uUpLabel.distBatch(b);
                    improved = dNew < dU;
                    dU = select(improved, dNew, dU);
                    if constexpr (LabelSet::KEEP_PARENT_EDGES) {
                        auto &eU = uUpLabel.pathEdgeBatch(b);
                        eU = select(improved, eAsBatch, eU);
                    }
                }
            }

            // Customize (reverse) downward label of u using upper neighbors
            auto uDownLabel = ttl.downLabel(u);
            uDownLabel.initializeLastHubDist(); // distance to self
            if constexpr (LabelSet::KEEP_PARENT_EDGES)
                uDownLabel.initializeLastHubPathEdge(); // edge to self
            FORALL_INCIDENT_EDGES(downGraph, u, e) {
                const auto v = downGraph.edgeHead(e);
                const auto downWeight = downWeights[e];
                const auto numHubsV = hierarchy.getNumHubs(v);
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                const auto vDownLabel = ttl.cDownLabel(v);
                const uint32_t numBlocks = numHubsV / K + (numHubsV % K != 0);
                const Batch weightAsBatch(downWeight);
                const Batch eAsBatch(e);
                for (uint32_t b = 0; b < numBlocks; ++b) {
                    const auto dNew = weightAsBatch + vDownLabel.distBatch(b);
                    auto &dU = uDownLabel.distBatch(b);
                    improved = dNew < dU;
                    dU = select(improved, dNew, dU);
                    if constexpr (LabelSet::KEEP_PARENT_EDGES) {
                        auto &eU = uDownLabel.pathEdgeBatch(b);
                        eU = select(improved, eAsBatch, eU);
                    }
                }
            }
        });
    }

    const BalancedTopologyCentricTreeHierarchy &hierarchy;
    const CCH &cch;
    CCHMetric cchMetric;

    CH minimumWeightedCH;
};

