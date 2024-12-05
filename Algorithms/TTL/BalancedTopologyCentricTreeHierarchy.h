#pragma once

#include <kassert/kassert.hpp>
#include <stack>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Algorithms/CCH/CCHMetric.h"

class BalancedTopologyCentricTreeHierarchy {

public:

    BalancedTopologyCentricTreeHierarchy() = default;

    // Builds the metric-independent CCH for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp) {

        KASSERT(hasStrictDissectionStructure(sepDecomp),
                "BalancedTopologyCentricTreeHierarchy requires strict dissection structure of separator decomposition.");
        const auto sdDepth = computeSepDecompDepth(sepDecomp);
        std::cout << "Depth of sepDecomp is " << sdDepth << std::endl;

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        denseSepSizeSum.clear();
        denseSepSizeSum.reserve((1 << (sdDepth + 1)));
        packedSideIds.resize(inputGraph.numVertices(), static_cast<uint64_t>(-1));
        std::vector<uint32_t> sepSizeSum;

        numHubs.resize(inputGraph.numVertices(), 0);
        initializeTreeHierarchy(sepDecomp, 0, 0, 0, 0);
        KASSERT(std::all_of(packedSideIds.begin(), packedSideIds.end(),
                            [](const uint64_t &id) { return id != static_cast<uint64_t>(-1); }));
    }

    const std::vector<uint32_t> &getNumHubs() const {
        return numHubs;
    }

    // Expects ranks in the CCH-order as inputs.
    uint32_t getLowestCommonHub(const int32_t &s, const int32_t &t) const {
        const auto minNumHubs = std::min(numHubs[s], numHubs[t]);

        // XOR packed side IDs to find out lowest common level in separator hierarchy.
        const auto l = lowestOneBit(packedSideIds[s] ^ packedSideIds[t]);

        // If packed side IDs of s and t are exactly the same, the branch of s subsumes the branch of t or vice
        // versa. In this case, the lowest common ancestor is the lower one of the two vertices.
        if (l == -1)
            return minNumHubs;

        const auto twoToTheL = 1 << l;

        // packedSideIds[s] and packedSideIds[t] are identical from bit 0 to bit l-1 (counting from least
        // significant).
        const auto commonLowerBits = (packedSideIds[s] & (twoToTheL - 1));
        KASSERT(commonLowerBits == (packedSideIds[t] & (twoToTheL - 1)));

        // Number of common hubs is defined by sum of size of separators up to and including separator that
        // separates s and t (LCA in separator hierarchy).
        // We find the separator hierarchy node representing this separator and read this value which is stored
        // for every separator node in sepSizeSum.

        // Find index of separator hierarchy node that contains the LCA:
        // Level of SH node is l. Packed side ID of LCA SH node is commonLowerBits.
        const auto &sepSizeSum = getSepSizeSum(l, commonLowerBits);
        return std::min(sepSizeSum, minNumHubs);
    }

private:

    // Returns true if every separator node in decomposition has at most two children, or false otherwise.
    static bool hasStrictDissectionStructure(const SeparatorDecomposition &sd) {
        for (const auto &n: sd.tree)
            if (n.rightSibling != 0 && sd.tree[n.rightSibling].rightSibling != 0)
                return false;
        return true;
    }

    static size_t computeSepDecompDepth(const SeparatorDecomposition &sd) {
        size_t maxDepth = 0;
        std::stack<uint32_t> sdNodesStack;
        sdNodesStack.push(0);
        bool returnedFromChildren = false;
        while (!sdNodesStack.empty()) {
            maxDepth = std::max(maxDepth, sdNodesStack.size());
            const auto node = sdNodesStack.top();
            if (!returnedFromChildren && sd.leftChild(node) != 0) {
                sdNodesStack.push(sd.leftChild(node));
                continue;
            }
            if (sd.rightSibling(node) != 0) {
                sdNodesStack.pop();
                sdNodesStack.push(sd.rightSibling(node));
                returnedFromChildren = false;
                continue;
            }
            sdNodesStack.pop(); // all children done, go back up in tree
            returnedFromChildren = true;
        }
        return maxDepth;
    }

    // Maps a depth and the packed side ID of a node to a unique node index.
    static uint64_t idxInSepSizeSum(const uint32_t depth, const uint64_t packedSideIdAtNode) {
        KASSERT((packedSideIdAtNode & (static_cast<uint64_t>(-1) << depth)) == 0);
        return (1 << depth) - 1 + packedSideIdAtNode;
    }

    // Set sum of sizes of separators on hierarchy branch up to separator hierarchy node (including size of separator
    // of node itself). Node is specified by depth and packedSideId.
    void setSepSizeSum(const uint32_t depth, const uint64_t packedSideId, const uint32_t sumSizes) {
        const auto idx = idxInSepSizeSum(depth, packedSideId);
        KASSERT(idx < denseSepSizeSum.capacity());
        if (idx >= denseSepSizeSum.size())
            denseSepSizeSum.resize(idx + 1, 0);
        denseSepSizeSum[idx] = sumSizes;
    }

    // Get sum of sizes of separators on hierarchy branch up to separator hierarchy node (including size of separator
    // of node itself). Node is specified by depth and packedSideId.
    const uint32_t &getSepSizeSum(const uint32_t depth, const uint64_t packedSideId) const {
        const auto idx = idxInSepSizeSum(depth, packedSideId);
        KASSERT(idx < denseSepSizeSum.size());
        return denseSepSizeSum[idx];
    }

    void initializeTreeHierarchy(const SeparatorDecomposition &sd,
                                 const int sdNode,
                                 const uint32_t depth,
                                 const uint64_t packedSideIdAtNode,
                                 const uint32_t sepSizeSumUntilHere) {

        KASSERT(depth == 0 || packedSideIdAtNode < (1 << depth));

        uint32_t inSepIdx = 0;
        for (auto v = sd.lastSeparatorVertex(sdNode) - 1; v >= sd.firstSeparatorVertex(sdNode); --v) {
            numHubs[v] = sepSizeSumUntilHere + inSepIdx + 1;
            packedSideIds[v] = packedSideIdAtNode;
            ++inSepIdx;
        }

        const uint32_t sepSizeSumIncludingThis = sepSizeSumUntilHere + inSepIdx;
        setSepSizeSum(depth, packedSideIdAtNode, sepSizeSumIncludingThis);

        // If this is a leaf node, return.
        if (sd.leftChild(sdNode) == 0) {
            return;
        }

        // Count children.
        uint32_t numChildren = 0;
        for (auto child = sd.leftChild(sdNode); child != 0; child = sd.rightSibling(child))
            ++numChildren;

        KASSERT(1 <= numChildren && numChildren <= 2);

        // For every child, call recursively with correct packed side ID
        uint32_t childIdx = 0;
        for (auto child = sd.leftChild(sdNode); child != 0; child = sd.rightSibling(child)) {
            const uint64_t childSideIdMask = childIdx == 0? 0 : 1 << depth;
            const uint64_t packedSideIdOfChild = packedSideIdAtNode | childSideIdMask;
            initializeTreeHierarchy(sd, child, depth + 1, packedSideIdOfChild, sepSizeSumIncludingThis);
            ++childIdx;
        }
        KASSERT(childIdx <= 2);
    }

    std::vector<uint32_t> numHubs; // For each vertex, stores number of hubs in label of vertex
    std::vector<uint64_t> packedSideIds; // store which side each vertex is on in each level of separator hierarchy

    // For each separator node, we have to store the sum of sizes of separators on the hierarchy branch up to the node.
    // These values need to be retrieved on the basis of a depth and the packedSideSum of the node.
    // We map depth and packedSideSum of a node to a unique index and store the sep size sum at that index.
    std::vector<uint32_t> denseSepSizeSum;
};

