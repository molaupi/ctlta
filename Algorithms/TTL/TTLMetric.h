#pragma once

#include "Algorithms/TTL/TruncatedTreeLabelling.h"
#include "Algorithms/TTL/BalancedTopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCHMetric.h"

// An individual metric for truncated tree labelling. Uses a CCH to build customized hub labelling.
class TTLMetric {

public:

    // Constructs an individual metric incorporating the specified input weights in the specified
    // BalancedTopologyCentricTreeHierarchy on the basis of the specified CCH.
    TTLMetric(const BalancedTopologyCentricTreeHierarchy& hierarchy, const CCH& cch, const int32_t* const inputWeights)
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights), minimumWeightedCH() {}

    void buildCustomizedTTL(TruncatedTreeLabelling& ttl) {
        minimumWeightedCH = cchMetric.buildMinimumWeightedCH();
        customizeLabelling(ttl);
    }

//    // Returns current weights on up edges of CCH.
//    const std::vector<int32_t>& getUpWeights() const {
//        return cchMetric.upWeights;
//    }
//
//    // Returns current weights on up edges of CCH.
//    const std::vector<int32_t>& getDownWeights() const {
//        return cchMetric.downWeights;
//    }

    // Returns reference to minimum weighted CH. Only call after buildCustomizedTTL.
    const CH& getMinimumWeightedCH() const {
        return minimumWeightedCH;
    }

private:

    void customizeLabelling(TruncatedTreeLabelling& ttl) {
        ttl.reset();
        const auto& upGraph = minimumWeightedCH.upwardGraph();
        const auto& downGraph = minimumWeightedCH.downwardGraph(); // reverse downward graph

        const auto& numHubs = hierarchy.getNumHubs();
#pragma omp parallel // parallelizes callbacks within cch.forEachVertexTopDown.
#pragma omp single nowait
        cch.forEachVertexTopDown([&](const int32_t& u) {

            // Do not build labels for truncated vertices.
            if (hierarchy.isVertexTruncated(u))
                return;

            const auto numHubsU = numHubs[u];

            // Customize upward label of u using upper neighbors
            ttl.upDist(u, numHubsU - 1) = 0; // distance to self
            ttl.upPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            FORALL_INCIDENT_EDGES(upGraph, u, e) {
                const auto v = upGraph.edgeHead(e);
                const auto upWeight = upGraph.traversalCost(e);
                const auto numHubsV = numHubs[v];
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                // TODO: SIMD-ify
                auto uUpLabel = ttl.upLabel(u);
                const auto vUpLabel = ttl.cUpLabel(v);
                for (uint32_t i = 0; i < numHubsV; ++i) {
                    if (uUpLabel.dist(i) > upWeight + vUpLabel.dist(i)) {
                        uUpLabel.dist(i) = upWeight + vUpLabel.dist(i);
                        uUpLabel.pathEdge(i) = e;
                    }
                }
            }

            // Customize (reverse) downward label of u using upper neighbors
            ttl.downDist(u, numHubsU - 1) = 0; // distance to self
            ttl.downPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            FORALL_INCIDENT_EDGES(downGraph, u, e) {
                const auto v = downGraph.edgeHead(e);
                const auto downWeight = downGraph.traversalCost(e);
                const auto numHubsV = numHubs[v];
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                // TODO: SIMD-ify
                auto uDownLabel = ttl.downLabel(u);
                const auto vDownLabel = ttl.cDownLabel(v);
                for (uint32_t i = 0; i < numHubsV; ++i) {
                    if (uDownLabel.dist(i) > downWeight + vDownLabel.dist(i)) {
                        uDownLabel.dist(i) = downWeight + vDownLabel.dist(i);
                        uDownLabel.pathEdge(i) = e;
                    }
                }
            }
        });
//        ttl.assertFullyCustomized();
    }

    const BalancedTopologyCentricTreeHierarchy& hierarchy;
    const CCH& cch;
    CCHMetric cchMetric;

    CH minimumWeightedCH;
};

