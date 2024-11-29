#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"

inline void printUsage() {
    std::cout <<
              "Usage: GraphToSimpleDimacs -g <file> -o <file>\n"
              "This program takes a graph in binary format and outputs a DIMACS format graph with only the travel time attribute.\n"
              "  -g <file>         graph file in binary format\n"
              "  -o <file>         output file without file extension\n"
              "  -help             display this help and exit\n";
}

// A graph data structure encompassing all vertex and edge attributes available for output.
using VertexAttributes = VertexAttrs<>;
using EdgeAttributes = EdgeAttrs<TravelTimeAttribute>;
using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;

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
        auto graph = GraphT(in);
        std::cout << " done." << std::endl;

        std::cout << "Writing the output file..." << std::flush;
        auto outfile = clp.getValue<std::string>("o");
        if (!endsWith(outfile, ".gr.dimacs"))
            outfile += ".gr.dimacs";
        std::ofstream out(outfile, std::ios::binary);
        if (!out.good())
            throw std::invalid_argument("file cannot be opened -- '" + outfile + ".gr.bin'");

        out << "p sp " << graph.numVertices() << " " << graph.numEdges() << std::endl;
        FORALL_VALID_EDGES(graph, u, e) {
            out << "a " << u + 1 << " " << graph.edgeHead(e) + 1 << " " << graph.travelTime(e) << std::endl;
        }
        std::cout << " done." << std::endl;

    } catch (std::invalid_argument &e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        std::cerr << "Try '" << argv[0] << " -help' for more information." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
