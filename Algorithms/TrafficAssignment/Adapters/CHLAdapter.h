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

#include "Algorithms/CHL/road_network.h"

namespace trafficassignment {

// An adapter that makes CCHs usable in the all-or-nothing assignment procedure.
    template<typename InputGraphT, typename WeightT>
    class CHLAdapter {
    public:

        // The number of source-target pairs per call to the query algorithm.
        // In the case of CHL, the K searches are computed sequentially.
        static constexpr int K = 1 << TA_LOG_K;

        using InputGraph = InputGraphT;

        // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
        // Multiple instances can work on the same data concurrently.
        class QueryAlgo {
        public:
            // Constructs a query algorithm instance working on the specified data.
            QueryAlgo(
                    const InputGraphT &inputGraph,
                    const chl::road_network::Graph &chlGraph,
                    const chl::road_network::ContractionHierarchy &ch,
                    const chl::road_network::ContractionIndex &ci,
                    AlignedVector<int> &globalFlow)
                    : inputGraph(inputGraph),
                      chlGraph(chlGraph),
                      ch(ch),
                      ci(ci),
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

                // No facilities for centralized searches in CHL. Run each search individually:
                for (auto i = 0; i < k; ++i) {

                    // CHL node IDs start at 1:
                    const auto chlSource = sources[i] + 1;
                    const auto chlTarget = targets[i] + 1;

                    chl::road_network::distance_t dist;
                    const auto vertexPath = chlGraph.query_path(ch, ci, chlSource, chlTarget, dist);

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
                    KASSERT(dist == recomputedDist, "dist returned by CHL is not correct", 20);

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
            const chl::road_network::Graph &chlGraph; // The CHL representation of the graph
            const chl::road_network::ContractionHierarchy &ch; // CHL underlying contraction hierarchy
            const chl::road_network::ContractionIndex &ci; // CHL underlying label structure

            std::array<int, K> distances; // distances computed in last call to run()

            AlignedVector<int> &globalFlow; // The global flows on edges.
            std::vector<int> localFlow; // The local flows on the input graph (computed by this search, i.e., this thread).
        };

        // Constructs an adapter for CCHs.
        explicit CHLAdapter(const InputGraphT &inputGraph)
                : inputGraph(inputGraph) {
            assert(inputGraph.numEdges() > 0);
            assert(inputGraph.isDefrag());
            verifyUndirectedTopology(inputGraph);
        }
// Invoked before the first iteration.
        void preprocess() {

            // Convert the input graph to CHL representation.
            chlGraph.resize(inputGraph.numVertices());
            FORALL_VALID_EDGES(inputGraph, u, e) {
                    // CHL graph node IDs start at 1
                    chlGraph.add_edge(u + 1, inputGraph.edgeHead(e) + 1, inputGraph.template get<WeightT>(e), false);
                }

            // Contract degree 1 nodes
            std::vector<chl::road_network::Neighbor> closest;
            chlGraph.contract(closest, false);

            // Build balanced tree hierarchy
            static constexpr double CUT_BALANCE = 0.2;
            static constexpr size_t LEAF_SIZE_THRESHOLD = 0;
            std::vector<chl::road_network::CutIndex> cutIndex;
            chlGraph.create_cut_index(cutIndex, CUT_BALANCE, LEAF_SIZE_THRESHOLD);

            // Reset graph to original form before contractions
            chlGraph.reset();

            // Initialize shortcut graph and labels
            chlGraph.initialize(ch, cutIndex, closest);
            ci = std::make_unique<chl::road_network::ContractionIndex>(cutIndex, closest);
        }

        // Invoked before each iteration.
        void customize() {

            // CHL can only deal with undirected graphs. Make sure traversal costs are the same in both directions.
            verifyUndirectedWeights(inputGraph);

            chlGraph.resetSpDatastructures(ch, *ci);

            // Customize CCH and HL with new metric on edges:
            std::vector<chl::road_network::Edge> edges;
            chlGraph.get_edges(edges);
            for (auto& e : edges) {
                // CHL graph node IDs start at 1
                const auto tail = e.a - 1;
                const auto head = e.b - 1;
                const auto eInInputGraph = inputGraph.uniqueEdgeBetween(tail, head);
                KASSERT(eInInputGraph >= 0 && eInInputGraph < inputGraph.numEdges());
                e.d = inputGraph.traversalCost(eInInputGraph);

            }
            chlGraph.customise_shortcut_graph(ch, *ci, edges);
            chlGraph.customise_hub_labelling(ch, *ci);

            // Initialize flows vector
            flow.assign(inputGraph.numEdges(), 0);
        }

        // Returns an instance of the query algorithm.
        QueryAlgo getQueryAlgoInstance() {
            KASSERT((bool) ci);
            return {inputGraph, chlGraph, ch, *ci, flow};
        }

        // Propagates the flows on the edges in the search graphs to the edges in the input graph.
        void propagateFlowsToInputEdges(AlignedVector<int> &flowsOnInputEdges) {
            KASSERT(flowsOnInputEdges.size() == inputGraph.numEdges());
            for (auto i = 0; i < inputGraph.numEdges(); ++i) {
                flowsOnInputEdges[i] += flow[i];
            }
        }

    private:

        static void verifyUndirectedTopology(const InputGraphT& g) {
            // Graph can be considered undirected if every edge exists both ways.
            FORALL_VALID_EDGES(g, u, e) {
                const auto eBack = g.uniqueEdgeBetween(g.edgeHead(e), u);
                KASSERT(eBack >= 0 && eBack < g.numEdges());
            }
        }

        static void verifyUndirectedWeights(const InputGraphT& g) {
            // Graph can be considered undirected if every edge exists both ways and the weight is the same both ways.
            FORALL_VALID_EDGES(g, u, e) {
                    KASSERT(g.template get<WeightT>(e) != WeightT::defaultValue());
                    const auto eBack = g.uniqueEdgeBetween(g.edgeHead(e), u);
                    KASSERT(eBack >= 0 && eBack < g.numEdges());
                    KASSERT(g.template get<WeightT>(e) == g.template get<WeightT>(eBack));
                }
        }

        const InputGraphT &inputGraph; // The input graph in TA format.
        chl::road_network::Graph chlGraph; // The graph in CHL format
        chl::road_network::ContractionHierarchy ch; // CHL underlying contraction hierarchy
        std::unique_ptr<chl::road_network::ContractionIndex> ci; // CHL underlying label structure

        AlignedVector<int> flow;   // The flows on the edges in the graph.
    };

}
