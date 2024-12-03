#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>
#include <kassert/kassert.hpp>
#include <routingkit/nested_dissection.h>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Simd/AlignedVector.h"


#include "Algorithms/TTL/TopologyCentricTreeHierarchy.h"
#include "Algorithms/TTL/TTLMetric.h"
#include "Algorithms/TTL/TTLQuery.h"

namespace trafficassignment {

// An adapter that makes TTLs usable in the all-or-nothing assignment procedure.
    template<typename InputGraphT, typename WeightT>
    class TTLAdapter {
    public:

        // The number of source-target pairs per call to the query algorithm.
        // In the case of TTL, the K searches are computed sequentially.
        static constexpr int K = 1 << TA_LOG_K;

        using InputGraph = InputGraphT;

        // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
        // Multiple instances can work on the same data concurrently.
        class QueryAlgo {
        public:
            // Constructs a query algorithm instance working on the specified data.
            QueryAlgo(const InputGraphT &inputGraph, const CCH &cch,
                      const TopologyCentricTreeHierarchy &hierarchy, const TruncatedTreeLabelling &ttl,
                      AlignedVector<int> &flowsOnUpEdges, AlignedVector<int> &flowsOnDownEdges) :
                    inputGraph(inputGraph), cch(cch), ttlQuery(hierarchy, ttl),
                    flowsOnUpEdges(flowsOnUpEdges),
                    flowsOnDownEdges(flowsOnDownEdges),
                    localFlowsOnUpEdges(flowsOnUpEdges.size(), 0),
                    localFlowsOnDownEdges(flowsOnDownEdges.size(), 0) {
                assert(inputGraph.numEdges() == globalFlow.size());
                distances.fill(INFTY);
            }

            // Computes shortest paths from each source to its target simultaneously.
            // Parameter k is actual number of source-target pairs if batch is partially filled,
            // i.e., only pairs 0..k are relevant.
            void run(std::array<int, K> &sources, std::array<int, K> &targets, const int k) {

                // No facilities for centralized searches in TTL. Run each search individually:
                for (auto j = 0; j < k; ++j) {
                    ttlQuery.run(cch.getRanks()[sources[j]], cch.getRanks()[targets[j]]);
                    distances[j] = static_cast<int>(ttlQuery.getDistance());

                    // Assign flow to the edges on the computed paths.
                    const auto &upEdgePath = ttlQuery.getUpEdgePath();
                    const auto &downEdgePath = ttlQuery.getDownEdgePath();
                    KASSERT(cch.getUpwardGraph().edgeHead(upEdgePath.back()) ==
                            cch.getUpwardGraph().edgeHead(downEdgePath.front()));
                    for (const auto e: upEdgePath) {
                        KASSERT(e >= 0);
                        KASSERT(e < localFlowsOnUpEdges.size());
                        ++localFlowsOnUpEdges[e];
                    }
                    for (const auto e: downEdgePath) {
                        KASSERT(e >= 0);
                        KASSERT(e < localFlowsOnDownEdges.size());
                        ++localFlowsOnDownEdges[e];
                    }

                }
            }

            // Returns the length of the i-th shortest path.
            int getDistance(const int /*dst*/, const int i) {
                return distances[i];
            }

            // Adds the local flow counters to the global ones. Must be synchronized externally.
            void addLocalToGlobalFlows() {
                FORALL_EDGES(cch.getUpwardGraph(), e) {
                    flowsOnUpEdges[e] += localFlowsOnUpEdges[e];
                    flowsOnDownEdges[e] += localFlowsOnDownEdges[e];
                }
            }

        private:

            const InputGraphT &inputGraph; // The TA representation of the graph
            const CCH &cch;
            TTLQuery ttlQuery;

            std::array<int, K> distances; // distances computed in last call to run()

            AlignedVector<int> &flowsOnUpEdges;     // The flows in the upward graph.
            AlignedVector<int> &flowsOnDownEdges;   // The flows in the downward graph.
            std::vector<int> localFlowsOnUpEdges;   // The local flows in the upward graph.
            std::vector<int> localFlowsOnDownEdges; // The local flows in the downward graph.
        };

        // Constructs an adapter for CCHs.
        explicit TTLAdapter(const InputGraphT &inputGraph)
                : inputGraph(inputGraph) {
            assert(inputGraph.numEdges() > 0);
            assert(inputGraph.isDefrag());
        }

        // Invoked before the first iteration.
        void preprocess() {
            // Convert the input graph to RoutingKit's graph representation.
            std::vector<float> lats(inputGraph.numVertices());
            std::vector<float> lngs(inputGraph.numVertices());
            std::vector<unsigned int> tails(inputGraph.numEdges());
            std::vector<unsigned int> heads(inputGraph.numEdges());
            FORALL_VERTICES(inputGraph, u) {
                lats[u] = inputGraph.latLng(u).latInDeg();
                lngs[u] = inputGraph.latLng(u).lngInDeg();
                FORALL_INCIDENT_EDGES(inputGraph, u, e) {
                    tails[e] = u;
                    heads[e] = inputGraph.edgeHead(e);
                }
            }

            // Compute a separator decomposition for the input graph.
            const auto graph = RoutingKit::make_graph_fragment(inputGraph.numVertices(), tails, heads);
            auto computeSep = [&](const RoutingKit::GraphFragment &fragment) {
                const auto cut = inertial_flow(fragment, 30, lats, lngs);
                return derive_separator_from_cut(fragment, cut.is_node_on_side);
            };
            const auto decomp = compute_separator_decomposition(graph, computeSep);

            // Convert the separator decomposition to our representation.
            SeparatorDecomposition sepDecomp;
            for (const auto &n: decomp.tree) {
                SeparatorDecomposition::Node node;
                node.leftChild = n.left_child;
                node.rightSibling = n.right_sibling;
                node.firstSeparatorVertex = n.first_separator_vertex;
                node.lastSeparatorVertex = n.last_separator_vertex;
                sepDecomp.tree.push_back(node);
            }
            sepDecomp.order.assign(decomp.order.begin(), decomp.order.end());

            // Build the CCH.
            cch.preprocess(inputGraph, sepDecomp);

            // Build the tree hierarchy.
            treeHierarchy.preprocess(inputGraph, sepDecomp);

            // Allocate labels.
            ttl.init();
        }

        // Invoked before each iteration.
        void customize() {
            currentMetric.buildCustomizedTTL(ttl);
            flowsOnUpEdges.assign(cch.getUpwardGraph().numEdges(), 0);
            flowsOnDownEdges.assign(cch.getUpwardGraph().numEdges(), 0);

        }

        // Returns an instance of the query algorithm.
        QueryAlgo getQueryAlgoInstance() {
            return {inputGraph, cch, treeHierarchy, ttl, flowsOnUpEdges, flowsOnDownEdges};
        }

        // Propagates the flows on the edges in the search graphs to the edges in the input graph.
        void propagateFlowsToInputEdges(AlignedVector<int> &flowsOnInputEdges) {
            unused(flowsOnInputEdges);
            // TODO: how to unpack without actually having CH?
//            const auto &upGraph = cch.getUpwardGraph();
//            for (auto u = inputGraph.numVertices() - 1; u >= 0; --u) {
//                FORALL_INCIDENT_EDGES(upGraph, u, e) {
//                    if (upGraph.unpackingInfo(e).second == INVALID_EDGE) {
//                        flowsOnInputEdges[upGraph.unpackingInfo(e).first] = flowsOnUpEdges[e];
//                    } else {
//                        flowsOnDownEdges[upGraph.unpackingInfo(e).first] += flowsOnUpEdges[e];
//                        flowsOnUpEdges[upGraph.unpackingInfo(e).second] += flowsOnUpEdges[e];
//                    }
//                    if (upGraph.unpackingInfo(e).second == INVALID_EDGE) {
//                        flowsOnInputEdges[upGraph.unpackingInfo(e).first] = flowsOnDownEdges[e];
//                    } else {
//                        flowsOnDownEdges[upGraph.unpackingInfo(e).first] += flowsOnDownEdges[e];
//                        flowsOnUpEdges[upGraph.unpackingInfo(e).second] += flowsOnDownEdges[e];
//                    }
//                }
//            }
        }

    private:

        const InputGraphT &inputGraph; // The input graph in TA format.
        TopologyCentricTreeHierarchy treeHierarchy; // Tree hierarchy underlying TTL
        CCH cch;                      // The metric-independent CCH.
        TTLMetric currentMetric;      // The current metric for the CCH.
        TruncatedTreeLabelling ttl;   // The customized tree labelling.

        AlignedVector<int> flowsOnUpEdges;   // The flows on the edges in the upward graph.
        AlignedVector<int> flowsOnDownEdges; // The flows on the edges in the downward graph.
    };

}
