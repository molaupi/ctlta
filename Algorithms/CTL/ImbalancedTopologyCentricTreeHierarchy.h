#pragma once

#include <kassert/kassert.hpp>
#include <stack>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Algorithms/CCH/CCHMetric.h"

// As opposed to BalancedTopologyCentricTreeHierarchy, this tree hierarchy can be based on a separator decomposition
// that is not strictly binary which may lead to an imbalanced tree hierarchy.
class ImbalancedTopologyCentricTreeHierarchy {

public:

    ImbalancedTopologyCentricTreeHierarchy() = default;

    // Builds the metric-independent CCH for the specified graph and separator decomposition.
    template<typename InputGraphT>
    void preprocess(const InputGraphT &inputGraph, const SeparatorDecomposition &sepDecomp) {
        const auto sdDepth = computeSepDecompDepth(sepDecomp);
        std::cout << "Depth of sepDecomp is " << sdDepth << std::endl;

        // Build labels and numCommonHubsComputer
        packedSideIds.clear();
        denseSepSizeSum.clear();
        denseSepSizeSum.reserve((1 << (MAX_DENSE_LEVEL + 1)));
        sparseSepSizeSum.clear();
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
        if (depth <= MAX_DENSE_LEVEL) {
            KASSERT(idx < denseSepSizeSum.capacity());
            if (idx >= denseSepSizeSum.size())
                denseSepSizeSum.resize(idx + 1, 0);
            denseSepSizeSum[idx] = sumSizes;
        } else {
            sparseSepSizeSum[idx] = sumSizes;
        }
    }

    // Get sum of sizes of separators on hierarchy branch up to separator hierarchy node (including size of separator
    // of node itself). Node is specified by depth and packedSideId.
    const uint32_t &getSepSizeSum(const uint32_t depth, const uint64_t packedSideId) const {
        const auto idx = idxInSepSizeSum(depth, packedSideId);
        if (depth <= MAX_DENSE_LEVEL) {
            KASSERT(idx < denseSepSizeSum.size());
            return denseSepSizeSum[idx];
        } else {
            KASSERT(sparseSepSizeSum.contains(idx));
            return sparseSepSizeSum.at(idx);
        }
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

        KASSERT(numChildren > 0);

        // sepDecomp can have more than two children per node if separator divides subgraph into more than two
        // parts. For the tree hierarchy, we want exactly two children per separator though, so LCA computation can
        // correctly rely on bits for left/right child.
        // Advance depth by l=ceil(log(x)) levels, giving enough (2^l >= x) combinations of bit strings to fit all
        // children. If l > 1, we pretend that there are separators of size 0 at levels in between. Those separators
        // on the non-existing levels get sepSizeSum equal to the parent of the x children so numCommonHubsComputer
        // correctly identifies the lowest common ancestor to be part of the parent.

        // depthDelta = ceil(log(numChildren))
        const uint32_t depthDelta = std::max(std::numeric_limits<uint32_t>::digits - numLeadingZeros(numChildren - 1),
                                             1);

        // Set sepSizeSum for empty levels (virtual empty separators):
        for (uint32_t l = 1; l < depthDelta; ++l) {
            const auto twoToTheL = (1 << l);
            // Virtual level has twoToTheL virtual separator nodes. Compute their packedSideIds and according
            // internal separator node indices, and set their sepSizeSum to that of the parent node (as their own
            // separators are empty).
            for (uint64_t idxInRelLevel = 0; idxInRelLevel < twoToTheL; ++idxInRelLevel) {
                const uint64_t packedSideIdOfVirtualSepNode = packedSideIdAtNode | (idxInRelLevel << depth);
                setSepSizeSum(depth + l, packedSideIdOfVirtualSepNode, sepSizeSumIncludingThis);
            }
        }

        // Get packedIds for level depthDelta and shift by depth to the left to get masks for packed side IDs of
        // children.
        std::vector<uint64_t> childSideIdMasks = computePackedIdsOnLevel(depthDelta);
        KASSERT(childSideIdMasks.size() >= numChildren);
        for (auto &m: childSideIdMasks)
            m <<= depth;

        // For every child, call recursively with correct packed side ID
        uint32_t childIdx = 0;
        for (auto child = sd.leftChild(sdNode); child != 0; child = sd.rightSibling(child)) {
            const uint64_t packedSideIdOfChild = packedSideIdAtNode | childSideIdMasks[childIdx];
            initializeTreeHierarchy(sd, child, depth + depthDelta, packedSideIdOfChild, sepSizeSumIncludingThis);
            ++childIdx;
        }
    }

    // Computes a vector of packed side IDs for all nodes on given lvl if higher levels have less significant bits.
    static std::vector<uint64_t> computePackedIdsOnLevel(const uint32_t lvl) {
        KASSERT(lvl >= 1);
        static std::vector<std::vector<uint64_t>> lookupTable = {
                {0, 1}, // 0, 1
                {0, 2, 1, 3}, // 00, 10, 01, 11
                {0, 4, 2, 6, 1, 5, 3, 7} // 000, 100, 010, 110, 001, 101, 011, 111
        };
        if (lvl - 1 < lookupTable.size())
            return lookupTable[lvl - 1];

        std::vector<uint64_t> res;
        computePackedIdsOnLevelRecursively(0, 0, lvl, res);
        return res;
    }

    static void
    computePackedIdsOnLevelRecursively(const uint32_t curLvl, const uint64_t curPackedId, const uint32_t tarLvl,
                                       std::vector<uint64_t> &res) {
        if (curLvl == tarLvl) {
            res.push_back(curPackedId);
            return;
        }
        computePackedIdsOnLevelRecursively(curLvl + 1, curPackedId, tarLvl, res); // Left child
        computePackedIdsOnLevelRecursively(curLvl + 1, curPackedId & (1 << curLvl), tarLvl, res); // Right child
    }

    std::vector<uint32_t> numHubs; // For each vertex, stores number of hubs in label of vertex
    std::vector<uint64_t> packedSideIds; // store which side each vertex is on in each level of separator hierarchy


    // For each separator node, we have to store the sum of sizes of separators on the hierarchy branch up to the node.
    // These values need to be retrieved on the basis of a depth and the packedSideSum of the node.
    // We can map depth and packedSideSum of a node to a unique value, however, this requires space in the order of
    // 2^d where d is the depth of the tree hierarchy. This is only worth it for dense parts of the hierarchy, i.e.,
    // levels at which most separator nodes exist. As the tree hierarchy may be highly skewed, this is not the case
    // everywhere and the lower levels may be very sparse.
    // Thus, we store values for a dense part of the hierarchy (up to level MAX_DENSE_LEVEL) in a full vector, and a
    // sparse part (everything at lower levels) using hashing.
    static constexpr size_t MAX_DENSE_LEVEL = 25;
    std::vector<uint32_t> denseSepSizeSum;
    std::unordered_map<uint64_t, uint32_t> sparseSepSizeSum;


};

