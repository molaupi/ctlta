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
        if (!hasStrictDissectionStructure(sepDecomp))
            throw std::invalid_argument("BalancedTopologyCentricTreeHierarchy requires strict dissection "
                                        "structure of separator decomposition.");
        const auto sdDepth = computeSepDecompDepth(sepDecomp);
        std::cout << "Depth of sepDecomp is " << sdDepth << std::endl;

        const auto sdNodesToTruncateAtUsingGeneric = getSepDecompNodesToTruncateAt(sepDecomp, 5);
        unused(sdNodesToTruncateAtUsingGeneric);

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        denseSepSizeSum.clear();
        denseSepSizeSum.reserve((1 << (sdDepth + 1)));
        packedSideIds.resize(inputGraph.numVertices(), static_cast<uint64_t>(-1));
        std::vector<uint32_t> sepSizeSum;

        numHubs.resize(inputGraph.numVertices(), 0);
        initializeTreeHierarchy(sepDecomp);

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


    template<typename RecurseCallbacKT,
            typename BacktrackCallbackT>
    static void forEachSepDecompNodeInDfsOrder(const SeparatorDecomposition &sd,
                                               RecurseCallbacKT recurse,
                                               BacktrackCallbackT backtrack) {
        std::stack<uint32_t> sdNodesStack;
        sdNodesStack.push(0);
        bool returnedFromChildren = false;
        while (true) {
            const auto node = sdNodesStack.top();

            if (!returnedFromChildren && sd.leftChild(node) != 0) {
                recurse(node, sd.leftChild(node));
                sdNodesStack.push(sd.leftChild(node));
                continue;
            }

            // Done with this node. If there are siblings continue with siblings, otherwise return to parent.
            sdNodesStack.pop();
            if (sdNodesStack.empty())
                break; // Finished when stack becomes empty

            backtrack(node, sdNodesStack.top());
            if (sd.rightSibling(node) != 0) {
                recurse(sdNodesStack.top(), sd.rightSibling(node));
                sdNodesStack.push(sd.rightSibling(node));
                returnedFromChildren = false;
            } else {
                returnedFromChildren = true;
            }
        }
    }

    static size_t computeSepDecompDepth(const SeparatorDecomposition &sd) {
        uint32_t maxDepth = 0;
        uint32_t curDepth = 1;
        forEachSepDecompNodeInDfsOrder(sd,
                                       [&](const int, const int) {
                                           ++curDepth;
                                           maxDepth = std::max(maxDepth, curDepth);
                                       },
                                       [&](const int, const int) {
                                           --curDepth;
                                       });
        return maxDepth;
    }

    // Traverses separator decomposition and computes nodes to truncate at s.t. the number of vertices in the subtree
    // rooted at each marked node is at most theta.
    static std::vector<int>
    getSepDecompNodesToTruncateAt(const SeparatorDecomposition &sd, const uint32_t theta) {
        std::vector<int> nodesToTruncateAt;
        std::stack<uint32_t> numVerticesOnBranch;
        numVerticesOnBranch.push(0);

        forEachSepDecompNodeInDfsOrder(
                sd,
                [&](const int, const int) {
                    numVerticesOnBranch.push(0);
                },
                [&](const int child, const int parent) {
                    // Backtrack form child to parent when child is done.
                    // Add number of vertices in separator at child to number of vertices
                    auto numVerticesChild = numVerticesOnBranch.top();
                    numVerticesOnBranch.pop();
                    numVerticesChild += sd.lastSeparatorVertex(child) - sd.firstSeparatorVertex(child);

                    // If number of vertices at child does not exceed theta, remove all children of child from
                    // nodesToTruncateAt and add child instead. Children of child are at the back of nodesToTruncateAt.
                    if (numVerticesChild <= theta) {
                        for (int childChild = sd.leftChild(child);
                             childChild != 0; childChild = sd.rightSibling(childChild))
                            KASSERT(contains(nodesToTruncateAt.begin(), nodesToTruncateAt.end(), childChild));
                        for (int childChild = sd.leftChild(child);
                             childChild != 0; childChild = sd.rightSibling(childChild))
                            nodesToTruncateAt.pop_back();
                        nodesToTruncateAt.push_back(child);
                    }

                    // Add number of vertices at child to number of vertices at parent
                    numVerticesOnBranch.top() += numVerticesChild;
                });

        // If whole graph has fewer than theta nodes, return root node
        KASSERT(numVerticesOnBranch.size() == 1);
        if (numVerticesOnBranch.top() <= theta)
            return {0};
        return nodesToTruncateAt;
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


    void initializeTreeHierarchy(const SeparatorDecomposition &sd) {

        bool lastBacktrackedFromLeft = false;
        uint64_t packedSideId = 0;
        std::stack<uint32_t> sepSizeSumsOnBranch;
        sepSizeSumsOnBranch.push(sd.lastSeparatorVertex(0) - sd.firstSeparatorVertex(0));

        // Set num hubs for root node separator vertices
        uint32_t inSepIdx = 0;
        for (auto v = sd.lastSeparatorVertex(0) - 1; v >= sd.firstSeparatorVertex(0); --v) {
            numHubs[v] = inSepIdx + 1;
            packedSideIds[v] = packedSideId;
            ++inSepIdx;
        }

        const auto recurse = [&] (const int parent, const int child) {

            // Set next bit in packedSideId to 0 for recursion to left child and to 1 for recursion to right child.
            const auto depthOfChild = static_cast<int>(sepSizeSumsOnBranch.size());
            setBit(packedSideId, depthOfChild - 1, lastBacktrackedFromLeft);

            // Set number of hubs for separator vertices at child.
            const auto sepSumSizeOnBranchUntilParent = sepSizeSumsOnBranch.top();
            uint32_t inSepIdx = 0;
            for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                numHubs[v] = sepSumSizeOnBranchUntilParent + inSepIdx + 1;
                packedSideIds[v] = packedSideId;
                ++inSepIdx;
            }

            // Memorize separator size up to and including child
            sepSizeSumsOnBranch.push(sepSumSizeOnBranchUntilParent + inSepIdx);
            setSepSizeSum(depthOfChild, packedSideId, sepSumSizeOnBranchUntilParent + inSepIdx);
        };

        const auto backtrack = [&] (const int child, const int parent) {
            // Retreat one level in sepSizeSumsOnBranch and packedSideId
            sepSizeSumsOnBranch.pop();
            const auto depthOfChild = static_cast<int>(sepSizeSumsOnBranch.size());
            lastBacktrackedFromLeft = !getBit(packedSideId, depthOfChild - 1); // Remember whether to do right child next
            setBit(packedSideId, depthOfChild - 1, false);
        };

        forEachSepDecompNodeInDfsOrder(sd, recurse, backtrack);
    }

    std::vector<uint32_t> numHubs; // For each vertex, stores number of hubs in label of vertex
    std::vector<uint64_t> packedSideIds; // store which side each vertex is on in each level of separator hierarchy

    // For each separator node, we have to store the sum of sizes of separators on the hierarchy branch up to the node.
    // These values need to be retrieved on the basis of a depth and the packedSideSum of the node.
    // We map depth and packedSideSum of a node to a unique index and store the sep size sum at that index.
    std::vector<uint32_t> denseSepSizeSum;
};

