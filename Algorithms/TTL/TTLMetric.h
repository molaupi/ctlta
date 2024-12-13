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
template<typename LabellingT, uint32_t logK = TTL_SIMD_LOGK, bool USE_PERFECT_CUSTOMIZATION = TTL_USE_PERFECT_CUSTOMIZATION>
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

    static_assert(logK != 1, "logK==1 is not supported for TTLMetric customization. "
                             "Use logK >= 2 for SIMD or logK = 0 for no SIMD.");
    static constexpr uint32_t K = 1 << logK;

    // Vectors of multiple data elements for use with SIMD instructions.
    using BooleanVector = std::conditional_t<logK == 2, Vec4ib, std::conditional_t<logK == 3, Vec8ib, Vec16ib>>;
    using IntegerVector = std::conditional_t<logK == 2, Vec4i, std::conditional_t<logK == 3, Vec8i, Vec16i>>;

    void customizeLabelling(LabellingT &ttl) {
        ttl.reset();
        const auto &upGraph = upwardGraph();
        const auto &downGraph = downwardGraph(); // reverse downward graph
        const auto upWeights = upwardWeights();
        const auto downWeights = downwardWeights();

#pragma omp parallel // parallelizes callbacks within cch.forEachVertexTopDown.
#pragma omp single nowait
        cch.forEachVertexTopDown([&](const int32_t &u) {

            // Do not build labels for truncated vertices.
            if (hierarchy.isVertexTruncated(u))
                return;

            const auto numHubsU = hierarchy.getNumHubs(u);

            // Customize upward label of u using upper neighbors
            ttl.upDist(u, numHubsU - 1) = 0; // distance to self
            if constexpr (LabellingT::StoresPathPointers)
                ttl.upPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            FORALL_INCIDENT_EDGES(upGraph, u, e) {
                const auto v = upGraph.edgeHead(e);
                const auto upWeight = upWeights[e];
                const auto numHubsV = hierarchy.getNumHubs(v);
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                auto uUpLabel = ttl.upLabel(u);
                const auto vUpLabel = ttl.cUpLabel(v);
                if constexpr (LabellingT::StoresPathPointers) {
                    if constexpr (logK >= 2) {
                        const uint32_t numFullBlocks = numHubsV / K;
                        const IntegerVector weightAsVec(upWeight);
                        const IntegerVector eAsVec(e);
                        IntegerVector dV;
                        IntegerVector dU;
                        IntegerVector eU;
                        BooleanVector improved;
                        uint32_t b = 0;
                        for (; b < numFullBlocks; ++b) {
                            dU.load(&uUpLabel.dist(b * K));
                            eU.load(&uUpLabel.pathEdge(b * K));
                            dV.load(&vUpLabel.dist(b * K));
                            dV += weightAsVec;
                            improved = dV < dU;
                            dU = select(improved, dV, dU);
                            eU = select(improved, eAsVec, eU);
                            dU.store(&uUpLabel.dist(b * K));
                            eU.store(&uUpLabel.pathEdge(b * K));
                        }
                        if (numHubsV != b * K) {
                            // Remaining partially filled block:
                            dU.load_partial(numHubsV - b * K, &uUpLabel.dist(b * K));
                            eU.load_partial(numHubsV - b * K, &uUpLabel.pathEdge(b * K));
                            dV.load_partial(numHubsV - b * K, &vUpLabel.dist(b * K));
                            dV += weightAsVec;
                            improved = dV < dU;
                            dU = select(improved, dV, dU);
                            eU = select(improved, eAsVec, eU);
                            dU.store_partial(numHubsV - b * K, &uUpLabel.dist(b * K));
                            eU.store_partial(numHubsV - b * K, &uUpLabel.pathEdge(b * K));
                        }
                    } else {
                        // For logK < 2, we do not use SIMD
                        for (uint32_t i = 0; i < numHubsV; ++i) {
                            if (uUpLabel.dist(i) > upWeight + vUpLabel.dist(i)) {
                                uUpLabel.dist(i) = upWeight + vUpLabel.dist(i);
                                uUpLabel.pathEdge(i) = e;
                            }
                        }
                    }
                } else {
                    // A good compiler will unroll and vectorize this loop better than we can manually.
                    // Do not manually use SIMD.
                    for (uint32_t i = 0; i < numHubsV; ++i)
                        uUpLabel.dist(i) = std::min(uUpLabel.dist(i), upWeight + vUpLabel.dist(i));
                }
            }

            // Customize (reverse) downward label of u using upper neighbors
            ttl.downDist(u, numHubsU - 1) = 0; // distance to self
            if constexpr (LabellingT::StoresPathPointers)
                ttl.downPathEdge(u, numHubsU - 1) = INVALID_EDGE; // edge to self
            FORALL_INCIDENT_EDGES(downGraph, u, e) {
                const auto v = downGraph.edgeHead(e);
                const auto downWeight = downWeights[e];
                const auto numHubsV = hierarchy.getNumHubs(v);
                KASSERT(numHubsV < numHubsU);
                KASSERT(numHubsV == hierarchy.getLowestCommonHub(u, v));
                auto uDownLabel = ttl.downLabel(u);
                const auto vDownLabel = ttl.cDownLabel(v);
                if constexpr (LabellingT::StoresPathPointers) {
                    if constexpr (logK >= 2) {

                        const uint32_t numFullBlocks = numHubsV / K;
                        const IntegerVector weightAsVec(downWeight);
                        const IntegerVector eAsVec(e);
                        IntegerVector dV;
                        IntegerVector dU;
                        IntegerVector eU;
                        BooleanVector improved;
                        uint32_t b = 0;
                        for (; b < numFullBlocks; ++b) {
                            dU.load(&uDownLabel.dist(b * K));
                            eU.load(&uDownLabel.pathEdge(b * K));
                            dV.load(&vDownLabel.dist(b * K));
                            dV += weightAsVec;
                            improved = dV < dU;
                            dU = select(improved, dV, dU);
                            eU = select(improved, eAsVec, eU);
                            dU.store(&uDownLabel.dist(b * K));
                            eU.store(&uDownLabel.pathEdge(b * K));
                        }
                        if (numHubsV != b * K) {
                            // Remaining partially filled block:
                            dU.load_partial(numHubsV - b * K, &uDownLabel.dist(b * K));
                            eU.load_partial(numHubsV - b * K, &uDownLabel.pathEdge(b * K));
                            dV.load_partial(numHubsV - b * K, &vDownLabel.dist(b * K));
                            dV += weightAsVec;
                            improved = dV < dU;
                            dU = select(improved, dV, dU);
                            eU = select(improved, eAsVec, eU);
                            dU.store_partial(numHubsV - b * K, &uDownLabel.dist(b * K));
                            eU.store_partial(numHubsV - b * K, &uDownLabel.pathEdge(b * K));
                        }
                    } else {
                        // For logK < 2, we do not use SIMD
                        for (uint32_t i = 0; i < numHubsV; ++i) {
                            if (uDownLabel.dist(i) > downWeight + vDownLabel.dist(i)) {
                                uDownLabel.dist(i) = downWeight + vDownLabel.dist(i);
                                uDownLabel.pathEdge(i) = e;
                            }
                        }
                    }
                } else {
                    // A good compiler will unroll and vectorize this loop better than we can manually.
                    // Do not manually use SIMD.
                    for (uint32_t i = 0; i < numHubsV; ++i)
                        uDownLabel.dist(i) = std::min(uDownLabel.dist(i), downWeight + vDownLabel.dist(i));
                }
            }
        });
//        ttl.assertFullyCustomized();
    }

    const BalancedTopologyCentricTreeHierarchy &hierarchy;
    const CCH &cch;
    CCHMetric cchMetric;

    CH minimumWeightedCH;
};

