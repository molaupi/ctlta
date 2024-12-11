#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>
#include <kassert/kassert.hpp>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Simd/AlignedVector.h"

#include "Algorithms/TTLSA/road_network.h"

namespace trafficassignment {

// An adapter that makes CCHs usable in the all-or-nothing assignment procedure.
    template<typename InputGraphT, typename WeightT>
    class TTLSACCHAdapter {
    public:

        // The number of source-target pairs per call to the query algorithm.
        // In the case of TTLSA, the K searches are computed sequentially.
        static constexpr int K = 1 << TA_LOG_K;

        using InputGraph = InputGraphT;

        // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
        // Multiple instances can work on the same data concurrently.
        class QueryAlgo {
        public:
            // Constructs a query algorithm instance working on the specified data.
            QueryAlgo(
                    const InputGraphT &inputGraph,
                    const ttlsa::road_network::Graph &ttlsaGraph,
                    const ttlsa::road_network::ContractionHierarchy &ch,
                    AlignedVector<int> &globalFlow)
                    : inputGraph(inputGraph),
                      ttlsaGraph(ttlsaGraph),
                      ch(ch),
                      distances(),
                      globalFlow(globalFlow),
                      localFlow(inputGraph.numEdges()) {
                assert(inputGraph.numEdges() == globalFlow.size());
                distances.fill(INFTY);
            }

            // Computes shortest paths from each source to its target simultaneously.
            // Parameter k is actual number of source-target pairs if batch is partially filled,
            // i.e., only pairs 0..k are relevant.
            void run(std::array<int, K> &sources, std::array<int, K> &targets, const int k) {

                // No facilities for centralized searches in TTLSA. Run each search individually:
                for (auto i = 0; i < k; ++i) {

                    // TTLSA node IDs start at 1:
                    const auto ttlsaSource = sources[i] + 1;
                    const auto ttlsaTarget = targets[i] + 1;

                    ttlsa::road_network::distance_t dist;
                    const auto vertexPath = ttlsaGraph.query_contraction_hierarchy(ch, ttlsaSource, ttlsaTarget, dist);

                    // Find edges on path and add flow for every edge:
                    int recomputedDist = 0;
                    auto tail = vertexPath[0] - 1; // -1 because inputGraph vertex IDs start at 0
                    for (auto j = 1; j < vertexPath.size(); ++j) {
                        const auto head = vertexPath[j] - 1; // -1 because inputGraph vertex IDs start at 0
                        const auto e = inputGraph.uniqueEdgeBetween(tail, head);
                        KASSERT(e != -1);
//                        paths[i].push_back(e);
                        ++localFlow[e];
                        recomputedDist += inputGraph.template get<WeightT>(e);

                        tail = head;
                    }
                    KASSERT(dist == recomputedDist, "dist returned by TTLSACCH is not correct");

                    distances[i] = static_cast<int>(dist);
                }
            }

            // Returns the length of the i-th shortest path.
            int getDistance(const int /*dst*/, const int i) {
                return distances[i];
            }

            // Adds the local flow counters to the global ones. Must be synchronized externally.
            void addLocalToGlobalFlows() {
                FORALL_EDGES(inputGraph, e) globalFlow[e] += localFlow[e];
            }

        private:

            const InputGraphT &inputGraph; // The TA representation of the graph
            const ttlsa::road_network::Graph &ttlsaGraph; // The TTLSA representation of the graph
            const ttlsa::road_network::ContractionHierarchy &ch; // TTLSA underlying contraction hierarchy

            std::array<int, K> distances; // distances computed in last call to run()

            AlignedVector<int> &globalFlow; // The global flows on edges.
            std::vector<int> localFlow; // The local flows on the input graph (computed by this search, i.e., this thread).
        };

        // Constructs an adapter for CCHs.
        explicit TTLSACCHAdapter(const InputGraphT &inputGraph)
                : inputGraph(inputGraph) {
            assert(inputGraph.numEdges() > 0);
            assert(inputGraph.isDefrag());
            verifyUndirectedTopology(inputGraph);
        }

// Invoked before the first iteration.
        void preprocess() {

            // Convert the input graph to TTLSA representation.
            ttlsaGraph.resize(inputGraph.numVertices());
            FORALL_VALID_EDGES(inputGraph, u, e) {
                    // TTLSA graph node IDs start at 1
                    ttlsaGraph.add_edge(u + 1, inputGraph.edgeHead(e) + 1, inputGraph.template get<WeightT>(e), false);
                }

            // Contract degree 1 nodes
            std::vector<ttlsa::road_network::Neighbor> closest;
            ttlsaGraph.contract(closest, false);

            // Build balanced tree hierarchy
            std::vector<ttlsa::road_network::CutIndex> cutIndex;
            static constexpr double CUT_BALANCE = 0.2;
            static constexpr size_t LEAF_SIZE_THRESHOLD = 0; // CCH only works with THETA = 0.
            ttlsaGraph.create_cut_index(cutIndex, CUT_BALANCE, LEAF_SIZE_THRESHOLD);

            // Reset graph to original form before contractions
            ttlsaGraph.reset();

            // Initialize shortcut graph and labels
            ttlsaGraph.initialize(ch, cutIndex, closest);
        }

        // Invoked before each iteration.
        void customize() {

            // TTLSA can only deal with undirected graphs. Make sure traversal costs are the same in both directions.
            verifyUndirectedWeights(inputGraph);

            ttlsaGraph.resetSpDatastructures(ch);

            // Customize CCH and HL with new metric on edges:
            std::vector<ttlsa::road_network::Edge> edges;
            ttlsaGraph.get_edges(edges);
            for (auto &e: edges) {
                // TTLSA graph node IDs start at 1
                const auto tail = e.a - 1;
                const auto head = e.b - 1;
                const auto eInInputGraph = inputGraph.uniqueEdgeBetween(tail, head);
                KASSERT(eInInputGraph >= 0 && eInInputGraph < inputGraph.numEdges());
                e.d = inputGraph.traversalCost(eInInputGraph);

            }
            ttlsaGraph.customise_shortcut_graph(ch, edges);

            // Initialize flows vector
            flow.assign(inputGraph.numEdges(), 0);
        }

        // Returns an instance of the query algorithm.
        QueryAlgo getQueryAlgoInstance() {
            return {inputGraph, ttlsaGraph, ch, flow};
        }

        // Propagates the flows on the edges in the search graphs to the edges in the input graph.
        void propagateFlowsToInputEdges(AlignedVector<int> &flowsOnInputEdges) {
            KASSERT(flowsOnInputEdges.size() == inputGraph.numEdges());
            for (auto i = 0; i < inputGraph.numEdges(); ++i) {
                flowsOnInputEdges[i] += flow[i];
            }
        }

    private:

        static void verifyUndirectedTopology(const InputGraphT &g) {
            // Graph can be considered undirected if every edge exists both ways.
            FORALL_VALID_EDGES(g, u, e) {
                    const auto eBack = g.uniqueEdgeBetween(g.edgeHead(e), u);
                    KASSERT(eBack >= 0 && eBack < g.numEdges());
                }
        }

        static void verifyUndirectedWeights(const InputGraphT &g) {
            // Graph can be considered undirected if every edge exists both ways and the weight is the same both ways.
            FORALL_VALID_EDGES(g, u, e) {
                    KASSERT(g.template get<WeightT>(e) != WeightT::defaultValue());
                    const auto eBack = g.uniqueEdgeBetween(g.edgeHead(e), u);
                    KASSERT(eBack >= 0 && eBack < g.numEdges());
                    KASSERT(g.template get<WeightT>(e) == g.template get<WeightT>(eBack));
                }
        }

        const InputGraphT &inputGraph; // The input graph in TA format.
        ttlsa::road_network::Graph ttlsaGraph; // The graph in TTLSA format
        ttlsa::road_network::ContractionHierarchy ch; // TTLSA underlying contraction hierarchy

        AlignedVector<int> flow;   // The flows on the edges in the graph.
    };

}
