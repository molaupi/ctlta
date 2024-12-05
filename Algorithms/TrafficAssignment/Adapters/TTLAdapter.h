#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>
#include <kassert/kassert.hpp>
#include "DataStructures/Partitioning/nested_strict_dissection.h"

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Simd/AlignedVector.h"


#include "Algorithms/TTL/BalancedTopologyCentricTreeHierarchy.h"
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
            QueryAlgo(const CCH &cch, const BalancedTopologyCentricTreeHierarchy &hierarchy, const TruncatedTreeLabelling &ttl,
                      AlignedVector<int> &flowsOnUpEdges, AlignedVector<int> &flowsOnDownEdges) :
                    inputGraph(inputGraph), cch(cch), ttlQuery(hierarchy, cch.getUpwardGraph(), ttl),
                    flowsOnUpEdges(flowsOnUpEdges),
                    flowsOnDownEdges(flowsOnDownEdges),
                    localFlowsOnUpEdges(flowsOnUpEdges.size(), 0),
                    localFlowsOnDownEdges(flowsOnDownEdges.size(), 0) {
                assert(cch.getUpwardGraph().numEdges() == flowsOnUpEdges.size());
                assert(cch.getUpwardGraph().numEdges() == flowsOnDownEdges.size());
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
                    const auto &upEdgePath = ttlQuery.getUpEdgePath(); // In reverse order by convention of Dijkstra-based approaches
                    const auto &downEdgePath = ttlQuery.getDownEdgePath(); // In forward order by convention of Dijkstra-based approaches
                    KASSERT((upEdgePath.empty() && downEdgePath.empty())
                            || (upEdgePath.empty() &&
                                cch.getRanks()[sources[j]] == cch.getUpwardGraph().edgeHead(downEdgePath.front()))
                            || (downEdgePath.empty() &&
                                cch.getUpwardGraph().edgeHead(upEdgePath.front()) == cch.getRanks()[targets[j]])
                            || (!upEdgePath.empty() && !downEdgePath.empty() &&
                                cch.getUpwardGraph().edgeHead(upEdgePath.front()) ==
                                cch.getUpwardGraph().edgeHead(downEdgePath.front())));
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

        // Constructs an adapter for TTLs.
        explicit TTLAdapter(const InputGraphT &inputGraph)
                : inputGraph(inputGraph), treeHierarchy(), cch(),
                  metric(treeHierarchy, cch, &inputGraph.template get<WeightT>(0)), ttl(treeHierarchy) {
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
            auto graph = RoutingKit::make_graph_fragment(inputGraph.numVertices(), tails, heads);
            auto computeCut = [&](const RoutingKit::GraphFragment &fragment) {
                return inertial_flow(fragment, 30, lats, lngs);
            };
            RoutingKit::BitVector all(graph.node_count(), true);
            RoutingKit::BitVector none = ~all;
            const auto decomp = compute_separator_decomposition_with_strict_dissection(std::move(graph), computeCut, std::move(none), std::move(all));

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
            metric.buildCustomizedTTL(ttl);
            flowsOnUpEdges.assign(cch.getUpwardGraph().numEdges(), 0);
            flowsOnDownEdges.assign(cch.getUpwardGraph().numEdges(), 0);

        }

        // Returns an instance of the query algorithm.
        QueryAlgo getQueryAlgoInstance() {
            return {cch, treeHierarchy, ttl, flowsOnUpEdges, flowsOnDownEdges};
        }

        // Propagates the flows on the edges in the search graphs to the edges in the input graph.
        void propagateFlowsToInputEdges(AlignedVector<int> &flowsOnInputEdges) {
            unused(flowsOnInputEdges);

            // For each vertex v in top-down order and for each upward edge e out of v, we propagate the upward and
            // downward flow of e down into the lower triangles that e shortcuts (or add the flow to the corresponding
            // input edge if e is not a shortcut).
            const auto &cchGraph = cch.getUpwardGraph();
            cch.forEachVertexTopDown([&](const int &v) {
                FORALL_INCIDENT_EDGES(cchGraph, v, e) {
                    const auto tail = cchGraph.edgeTail(e);
                    const auto head = cchGraph.edgeHead(e);

                    // Propagate upward flow on e to input edge if e is not an up shortcut or down into the
                    // lower triangle it shortcuts.
                    if (flowsOnUpEdges[e] > 0) {
                        int upInputEdge = INVALID_EDGE;
                        if (!isUpShortCut(e, upInputEdge)) {
                            // If e is not a shortcut edge in upwards direction, add its upward flow to the according
                            // upwards input edge.
                            flowsOnInputEdges[upInputEdge] += flowsOnUpEdges[e];
                        } else {
                            // If e is a shortcut edge in upwards direction, find lower triangle that it shortcuts,
                            // and add upward flow of e to the downward flow of the lower edge in the triangle and to
                            // the upward flow of the middle edge in the triangle.
                            const auto noTriangleFound = cch.forEachLowerTriangle(
                                    tail, head, e, [&](int, const int lower, const int inter) {
                                        if (metric.getDownWeight(lower) + metric.getUpWeight(inter) ==
                                            metric.getUpWeight(e)) {
                                            flowsOnDownEdges[lower] += flowsOnUpEdges[e];
                                            flowsOnUpEdges[inter] += flowsOnUpEdges[e];
                                            return false;
                                        }
                                        return true;
                                    });
                            unused(noTriangleFound);
                            KASSERT(!noTriangleFound);
                        }
                    }


                    // Propagate downward flow on e to input edge if e is not a down shortcut or down into the
                    // lower triangle it shortcuts.
                    if (flowsOnDownEdges[e] > 0) {
                        int downInputEdge = INVALID_EDGE;
                        if (!isDownShortCut(e, downInputEdge)) {
                            // If e is not a shortcut edge in downwards direction, add its downward flow to the according
                            // downwards input edge.
                            flowsOnInputEdges[downInputEdge] += flowsOnDownEdges[e];
                        } else {
                            // If e is a shortcut edge in downwards direction, find lower triangle that it shortcuts,
                            // and add downward flow of e to the downward flow of the middle edge in the triangle and to
                            // the upward flow of the lower edge in the triangle.
                            const auto noTriangleFound = cch.forEachLowerTriangle(
                                    tail, head, e, [&](int, const int lower, const int inter) {
                                        if (metric.getDownWeight(inter) + metric.getUpWeight(lower) ==
                                            metric.getDownWeight(e)) {
                                            flowsOnDownEdges[inter] += flowsOnDownEdges[e];
                                            flowsOnUpEdges[lower] += flowsOnDownEdges[e];
                                            return false;
                                        }
                                        return true;
                                    });
                            unused(noTriangleFound);
                            KASSERT(!noTriangleFound);
                        }
                    }
                }
            });
        }

    private:

        // Returns true iff edge e in the CCH graph is an upward shortcut edge.
        // If e is not an upward shortcut edge, sets eInInputGraph to the according edge ID in the input graph.
        bool isUpShortCut(const int &e, int &eInInputGraph) {
            return cch.forEachUpwardInputEdge(e, [&](const int inputEdge) {
                if (inputGraph.traversalCost(inputEdge) == metric.getUpWeight(e)) {
                    eInInputGraph = inputEdge;
                    return false;
                }
                return true;
            });
        }

        // Returns true iff edge e in the CCH graph is a downward shortcut edge.
        // If e is not a downward shortcut edge, sets eInInputGraph to the according edge ID in the input graph.
        bool isDownShortCut(const int &e, int &eInInputGraph) {
            return cch.forEachDownwardInputEdge(e, [&](const int inputEdge) {
                if (inputGraph.traversalCost(inputEdge) == metric.getDownWeight(e)) {
                    eInInputGraph = inputEdge;
                    return false;
                }
                return true;
            });
        }

        const InputGraphT &inputGraph; // The input graph in TA format.
        BalancedTopologyCentricTreeHierarchy treeHierarchy; // Tree hierarchy underlying TTL
        CCH cch;                      // The metric-independent CCH.
        TTLMetric metric;      // The current metric for the CCH.
        TruncatedTreeLabelling ttl;   // The customized tree labelling.

        AlignedVector<int> flowsOnUpEdges;   // The flows on the edges in the upward graph.
        AlignedVector<int> flowsOnDownEdges; // The flows on the edges in the downward graph.
    };

}
