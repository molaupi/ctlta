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

        uint32_t maxUntruncatedDepth = 1;
        const auto sdNodesToTruncateAt = getSepDecompNodesToTruncateAt(sepDecomp, 0, maxUntruncatedDepth);
        std::cout << "Max depth of untruncated node in sepDecomp is " << maxUntruncatedDepth << std::endl;

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        sepSizeSum.clear();
        sepSizeSum.reserve((1 << (sdDepth + 1)));
        packedSideIds.resize(inputGraph.numVertices(), static_cast<uint64_t>(-1));

        numHubs.resize(inputGraph.numVertices(), 0);
        initializeTreeHierarchy(sepDecomp, sdNodesToTruncateAt);

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
    // rooted at each marked node is at most theta. Also sets maximum depth of an untruncated node.
    static BitVector
    getSepDecompNodesToTruncateAt(const SeparatorDecomposition &sd, const uint32_t theta, uint32_t& maxUntruncatedDepth) {
        BitVector doTruncate(sd.tree.size());
        std::stack<uint32_t> numVerticesOnBranch;
        numVerticesOnBranch.push(0);
        uint32_t curDepth = 1;
        maxUntruncatedDepth = 1;

        forEachSepDecompNodeInDfsOrder(
                sd,
                [&](const int, const int) {
                    numVerticesOnBranch.push(0);
                    ++curDepth;
                },
                [&](const int child, const int parent) {
                    // Backtrack from child to parent when child is done.
                    // Add number of vertices in separator at child to number of vertices
                    auto numVerticesChild = numVerticesOnBranch.top();
                    numVerticesOnBranch.pop();
                    numVerticesChild += sd.lastSeparatorVertex(child) - sd.firstSeparatorVertex(child);

                    // If number of vertices at child does not exceed theta, mark it to be treated as truncated.
                    if (numVerticesChild <= theta) {
                        for (int grandChild = sd.leftChild(child);
                             grandChild != 0; grandChild = sd.rightSibling(grandChild))
                            KASSERT(doTruncate[grandChild]);
                        doTruncate[child] = true;
                    } else {
                        maxUntruncatedDepth = std::max(maxUntruncatedDepth, curDepth);
                    }

                    --curDepth;

                    // Add number of vertices at child to number of vertices at parent
                    numVerticesOnBranch.top() += numVerticesChild;
                });

        // Potentially set doTruncate for root node
        KASSERT(numVerticesOnBranch.size() == 1);
        if (numVerticesOnBranch.top() <= theta) {
            KASSERT(maxUntruncatedDepth == 1);
            doTruncate[0] = true;
        }

        return doTruncate;
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
        KASSERT(idx < sepSizeSum.capacity());
        if (idx >= sepSizeSum.size())
            sepSizeSum.resize(idx + 1, 0);
        sepSizeSum[idx] = sumSizes;
    }

    // Get sum of sizes of separators on hierarchy branch up to separator hierarchy node (including size of separator
    // of node itself). Node is specified by depth and packedSideId.
    const uint32_t &getSepSizeSum(const uint32_t depth, const uint64_t packedSideId) const {
        const auto idx = idxInSepSizeSum(depth, packedSideId);
        KASSERT(idx < sepSizeSum.size());
        return sepSizeSum[idx];
    }


    void initializeTreeHierarchy(const SeparatorDecomposition &sd, const BitVector &truncateAtSdNode) {

        std::fill(numHubs.begin(), numHubs.end(), 0);

        // If root is truncated, give same packed side ID to all vertices and return. In this case, no tree hierarchy
        // is used at all.
        if (truncateAtSdNode[0]) {
            std::fill(packedSideIds.begin(), packedSideIds.end(), 0);
            return;
        }

        uint64_t packedSideId = 0;
        std::stack<uint32_t> sepSizeSumsOnBranch;
        std::stack<bool> doneWithLeftChild;
        sepSizeSumsOnBranch.push(sd.lastSeparatorVertex(0) - sd.firstSeparatorVertex(0));
        doneWithLeftChild.push(false);
        uint32_t depth = 1;

        // Set num hubs for root node separator vertices
        uint32_t inSepIdx = 0;
        for (auto v = sd.lastSeparatorVertex(0) - 1; v >= sd.firstSeparatorVertex(0); --v) {
            numHubs[v] = inSepIdx + 1;
            packedSideIds[v] = packedSideId;
            ++inSepIdx;
        }

        const auto recurse = [&](const int parent, const int child) {
            KASSERT(!truncateAtSdNode[parent] || truncateAtSdNode[child]);
            KASSERT(doneWithLeftChild.size() == depth);
            KASSERT(sepSizeSumsOnBranch.size() <= depth);

            // If the parent is not truncated, the child will get a new side ID with a new bit stating which
            // of the two children it is. If the parent is truncated, the child is also truncated and gets the same
            // side ID as the parent.
            if (!truncateAtSdNode[parent]) {
                // Set next bit in packedSideId to 0 for recursion to left child and to 1 for recursion to right child.
                setBit(packedSideId, depth - 1, doneWithLeftChild.top());
            }

            if (truncateAtSdNode[child]) {
                // If child is truncated, the vertices in the subtree rooted at child will have no labels, which means
                // they will have numHubs[v] = 0 and the sum of separator sizes on the branch is irrelevant. Set only the
                // side IDs of vertices.
                for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v)
                    packedSideIds[v] = packedSideId;
            } else {
                // If child is not truncated, set number of hubs and side IDs for separator vertices at child.

                const auto sepSumSizeOnBranchUntilParent = sepSizeSumsOnBranch.top();
                uint32_t inSepIdx = 0;
                for (auto v = sd.lastSeparatorVertex(child) - 1; v >= sd.firstSeparatorVertex(child); --v) {
                    numHubs[v] = sepSumSizeOnBranchUntilParent + inSepIdx + 1;
                    packedSideIds[v] = packedSideId;
                    ++inSepIdx;
                }

                // Memorize separator size up to and including child
                KASSERT(sepSizeSumsOnBranch.size() == doneWithLeftChild.size());
                sepSizeSumsOnBranch.push(sepSumSizeOnBranchUntilParent + inSepIdx);
                setSepSizeSum(depth, packedSideId, sepSumSizeOnBranchUntilParent + inSepIdx);
            }

            // Memorize that we are not done with left child of child
            doneWithLeftChild.push(false);
            ++depth;

        };

        const auto backtrack = [&](const int child, const int parent) {
            KASSERT(doneWithLeftChild.size() == depth);
            KASSERT(sepSizeSumsOnBranch.size() <= depth);

            // Remember that (at least) one child of parent is done. This way, after left child of parent is done,
            // we know that next recursion from parent is right child.
            doneWithLeftChild.pop();
            doneWithLeftChild.top() = true;
            --depth;

            // If parent was not truncated, packed ID of child has one bit more than that of parent. Reset this bit
            // since we are going back to parent.
            if (!truncateAtSdNode[parent]) {
                setBit(packedSideId, depth - 1, false);
            }

            // If child was not truncated, it has an entry on sepSizeSumsOnBranch. Remove this entry.
            if (!truncateAtSdNode[child])
                sepSizeSumsOnBranch.pop();
        };

        forEachSepDecompNodeInDfsOrder(sd, recurse, backtrack);
    }

    std::vector<uint32_t> numHubs; // For each vertex, stores number of hubs in label of vertex
    std::vector<uint64_t> packedSideIds; // store which side each vertex is on in each level of separator hierarchy

    // For each separator node, we have to store the sum of sizes of separators on the hierarchy branch up to the node.
    // These values need to be retrieved on the basis of a depth and the packedSideSum of the node.
    // We map depth and packedSideSum of a node to a unique index and store the sep size sum at that index.
    std::vector<uint32_t> sepSizeSum;
};

