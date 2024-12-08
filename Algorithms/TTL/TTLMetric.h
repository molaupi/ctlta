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
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights) {}

    void buildCustomizedTTL(TruncatedTreeLabelling& ttl) {
        cchMetric.customize();
        customizeLabelling(ttl);
    }

    // Returns current weights on up edges of CCH.
    const std::vector<int32_t>& getUpWeights() const {
        return cchMetric.upWeights;
    }

    // Returns current weights on up edges of CCH.
    const std::vector<int32_t>& getDownWeights() const {
        return cchMetric.downWeights;
    }

private:

    void customizeLabelling(TruncatedTreeLabelling& ttl) {
        ttl.reset();
        const auto& cchGraph = cch.getUpwardGraph();
        const auto& upWeights = cchMetric.upWeights;
        const auto& downWeights = cchMetric.downWeights;
        const auto& numHubs = hierarchy.getNumHubs();
#pragma omp parallel // parallelizes callbacks within cch.forEachVertexTopDown.
#pragma omp single nowait
        cch.forEachVertexTopDown([&](const int32_t& u) {

            // Do not build labels for truncated vertices.
            if (hierarchy.isVertexTruncated(u))
                return;

            const auto numHubsU = numHubs[u];
            ttl.upDist(u, numHubsU - 1) = 0; // distance to self
            ttl.upPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            ttl.downDist(u, numHubsU - 1) = 0; // distance to self
            ttl.downPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            FORALL_INCIDENT_EDGES(cchGraph, u, e) {
                const auto v = cchGraph.edgeHead(e);
                const auto upWeight = upWeights[e];
                const auto downWeight = downWeights[e];
                const auto numHubsV = numHubs[v];
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                // TODO: SIMD-ify
                auto uUpLabel = ttl.upLabel(u);
                auto uDownLabel = ttl.downLabel(u);
                const auto vUpLabel = ttl.cUpLabel(v);
                const auto vDownLabel = ttl.cDownLabel(v);
                for (uint32_t i = 0; i < numHubsV; ++i) {
                    if (uUpLabel.dist(i) > upWeight + vUpLabel.dist(i)) {
                        uUpLabel.dist(i) = upWeight + vUpLabel.dist(i);
                        uUpLabel.pathEdge(i) = e;
                    }
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
};

