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
        cch.forEachVertexTopDown([&](const int32_t& u) {

            // Do not build labels for truncated vertices.
            if (hierarchy.isVertexTruncated(u))
                return;

            const auto& numHubsU = hierarchy.getNumHubs()[u];
            ttl.upDist(u, numHubsU - 1) = 0; // distance to self
            ttl.upPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            ttl.downDist(u, numHubsU - 1) = 0; // distance to self
            ttl.downPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), u, e) {
                const auto v = cch.getUpwardGraph().edgeHead(e);
                const auto& numHubsV = hierarchy.getNumHubs()[v];
                const auto& upWeight = cchMetric.upWeights[e];
                const auto& downWeight = cchMetric.downWeights[e];
                KASSERT(numHubsV < numHubsU);
                const auto lch = hierarchy.getLowestCommonHub(u, v);
                // TODO: SIMD-ify
                for (uint32_t i = 0; i < lch; ++i) {
                    if (ttl.upDist(u, i) > upWeight + ttl.upDist(v, i)) {
                        ttl.upDist(u, i) = upWeight + ttl.upDist(v, i);
                        ttl.upPathEdge(u, i) = e;
                    }
                    if (ttl.downDist(u, i) > downWeight + ttl.downDist(v, i)) {
                        ttl.downDist(u, i) = downWeight + ttl.downDist(v, i);
                        ttl.downPathEdge(u, i) = e;
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

