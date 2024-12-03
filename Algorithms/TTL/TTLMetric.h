#pragma once

#include "Algorithms/TTL/TruncatedTreeLabelling.h"
#include "Algorithms/TTL/TopologyCentricTreeHierarchy.h"
#include "Algorithms/CCH/CCHMetric.h"

// An individual metric for truncated tree labelling. Uses a CCH to build customized hub labelling.
class TTLMetric {

public:

    // Constructs an individual metric incorporating the specified input weights in the specified
    // TopologyCentricTreeHierarchy on the basis of the specified CCH.
    TTLMetric(const TopologyCentricTreeHierarchy& hierarchy, const CCH& cch, const int32_t* const inputWeights)
            : hierarchy(hierarchy), cch(cch), cchMetric(cch, inputWeights) {}

    void buildCustomizedTTL(TruncatedTreeLabelling& ttl) {
        cchMetric.customize();
        customizeLabelling(ttl);
    }

private:

    void customizeLabelling(TruncatedTreeLabelling& ttl) {
        ttl.reset();
        cch.forEachVertexTopDown([&](const int32_t& u) {
            const auto& numHubsU = hierarchy.getNumHubs()[u];
            ttl.upDist(u, numHubsU - 1) = 0; // distance to self
            ttl.upPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            ttl.upHub(u, numHubsU - 1) = u;
            ttl.downDist(u, numHubsU - 1) = 0; // distance to self
            ttl.downPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            ttl.downHub(u, numHubsU - 1) = u;
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
                        ttl.upHub(u, i) = ttl.upHub(v, i);
                    }
                    if (ttl.downDist(u, i) > downWeight + ttl.downDist(v, i)) {
                        ttl.downDist(u, i) = downWeight + ttl.downDist(v, i);
                        ttl.downPathEdge(u, i) = e;
                        ttl.downHub(u, i) = ttl.downHub(v, i);
                    }
                }
            }
        });
//        ttl.assertFullyCustomized();
    }

    const TopologyCentricTreeHierarchy& hierarchy;
    const CCH& cch;
    CCHMetric cchMetric;
};

