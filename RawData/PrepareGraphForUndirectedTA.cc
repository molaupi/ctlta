#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"

inline void printUsage() {
    std::cout <<
              "Usage: PrepareGraphForUndirectedTA -g <file> -o <file>\n"
              "This program makes sure that a graph is ready for undirected traffic assignment.\n"
              "It makes sure that every edge exists in both directions and has the same (non-zero!) travel time in both directions.\n"
              "  -g <file>          graph file in binary format\n"
              "  -o <file>         output file(s) without file extension\n"
              "  -help             display this help and exit\n";
}

// A graph data structure encompassing all vertex and edge attributes available for output.
using VertexAttributes = VertexAttrs<LatLngAttribute, SequentialVertexIdAttribute>;
using EdgeAttributes = EdgeAttrs<
        CapacityAttribute, EdgeIdAttribute, LengthAttribute, TravelTimeAttribute,
        TraversalCostAttribute>;
using graphT = DynamicGraph<VertexAttributes, EdgeAttributes>;

int main(int argc, char *argv[]) {
    try {
        CommandLineParser clp(argc, argv);
        if (clp.isSet("help")) {
            printUsage();
            return EXIT_SUCCESS;
        }

        std::cout << "Reading the input file(s)..." << std::flush;

        auto infile = clp.getValue<std::string>("g");
        if (!endsWith(infile, ".gr.bin"))
            infile += ".gr.bin";
        std::ifstream in(infile, std::ios::binary);
        if (!in.good())
            throw std::invalid_argument("file not found -- '" + infile + ".gr.bin'");
        auto graph = graphT(in);
        std::cout << " done." << std::endl;

        std::cout << "Unifying edges in both directions..." << std::flush;

        FORALL_VALID_EDGES(graph, u, e) {
                const auto v = graph.edgeHead(e);
                if (graph.travelTime(e) == 0) // eliminate travel time zero
                    graph.travelTime(e) = 1;
                if (graph.traversalCost(e) == 0)
                    graph.traversalCost(e) = 1;
                auto eBack = graph.uniqueEdgeBetween(v, u);
                if (eBack == -1) {
                    // If edge in back direction does not exist, add it
                    eBack = graph.insertEdge(v, u);
                    graph.capacity(eBack) = graph.capacity(e);
                    graph.length(eBack) = graph.length(e);
                    graph.travelTime(eBack) = graph.travelTime(e);
                    graph.traversalCost(eBack) = graph.traversalCost(e);
                } else {
                    // If edge does exist in back direction, unify values of attributes (by taking worse value).
                    graph.capacity(eBack) = std::min(graph.capacity(e), graph.capacity(eBack));
                    graph.length(eBack) = std::max(graph.length(e), graph.length(eBack));
                    graph.travelTime(eBack) = std::max(graph.travelTime(e), graph.travelTime(eBack));
                    graph.traversalCost(eBack) = std::max(graph.traversalCost(e), graph.traversalCost(eBack));
                }
            }
        graph.defrag();
        std::cout << " done." << std::endl;


        std::cout << "Writing the output file..." << std::flush;
        auto outfile = clp.getValue<std::string>("o");
        if (!endsWith(outfile, ".gr.bin"))
            outfile += ".gr.bin";
        std::ofstream out(outfile, std::ios::binary);
        if (!out.good())
            throw std::invalid_argument("file cannot be opened -- '" + outfile + ".gr.bin'");
        graph.writeTo(out);
        std::cout << " done." << std::endl;

    } catch (std::invalid_argument &e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        std::cerr << "Try '" << argv[0] << " -help' for more information." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
