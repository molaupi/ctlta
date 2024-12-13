#pragma once

#include "DataStructures/Labels/BasicLabelSet.h"
#include "Algorithms/TTL/TruncatedTreeLabelling.h"
#include "Algorithms/Dijkstra/DagShortestPaths.h"

template<typename SearchGraphT, typename LabellingT, uint32_t logK = TTL_SIMD_LOGK>
class TTLQuery {


    // Temporary label being constructed for truncated vertices.
    struct TemporaryLabel {

        TemporaryLabel() = default;

        void init(const size_t numHubs) {
            dists.assign(numHubs, INFTY);
            if constexpr (LabellingT::StoresPathPointers)
                accessVertices.assign(numHubs, INVALID_VERTEX);
        }

        const int32_t &dist(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs());
            return dists[hubIdx];
        }

        int32_t &dist(const uint32_t &hubIdx) {
            KASSERT(hubIdx < numHubs());
            return dists[hubIdx];
        }

        template<bool hasPathEdges = LabellingT::StoresPathPointers, std::enable_if_t<hasPathEdges>...>
        const int32_t &accessVertex(const uint32_t &hubIdx) const {
            KASSERT(hubIdx < numHubs());
            return accessVertices[hubIdx];
        }

        template<bool hasPathEdges = LabellingT::StoresPathPointers, std::enable_if_t<hasPathEdges>...>
        int32_t &accessVertex(const uint32_t &hubIdx) {
            KASSERT(hubIdx < numHubs());
            return accessVertices[hubIdx];
        }

        uint32_t numHubs() const {
            return static_cast<uint32_t>(dists.size());
        }

    private:
        std::vector<int32_t> dists;
        std::vector<int32_t> accessVertices;
    };

    template<bool UP>
    struct PruneSearchAtUntruncatedVertices {

        PruneSearchAtUntruncatedVertices(TemporaryLabel &temporaryLabel,
                                         std::vector<int32_t> &truncatedSearchSpace,
                                         const BalancedTopologyCentricTreeHierarchy &hierarchy,
                                         const LabellingT &ttl)
                : temporaryLabel(temporaryLabel), truncatedSearchSpace(truncatedSearchSpace), hierarchy(hierarchy),
                  ttl(ttl) {}

        template<typename DistanceLabelT, typename DistanceLabelContainerT>
        bool operator()(const int v, DistanceLabelT &distToV, DistanceLabelContainerT &) {

            // If v is truncated, store it in truncated search space, and do not prune.
            if (hierarchy.isVertexTruncated(v)) {
                truncatedSearchSpace.push_back(v);
                return false;
            }

            // Update temporary label of source / target with label of v.
            const auto labelOfV = UP ? ttl.upLabel(v) : ttl.downLabel(v);
            const auto maxHub = std::min(temporaryLabel.numHubs(), labelOfV.numHubs);
            for (auto i = 0; i < maxHub; ++i) {
                const auto distViaV = distToV[0] + labelOfV.dist(i);
                if (distViaV < temporaryLabel.dist(i)) {
                    temporaryLabel.dist(i) = distViaV;
                    if constexpr (LabellingT::StoresPathPointers)
                        temporaryLabel.accessVertex(i) = v;
                }
            }

            // Prune at untruncated vertices.
            return true;
        }

        TemporaryLabel &temporaryLabel;
        std::vector<int32_t> &truncatedSearchSpace;
        const BalancedTopologyCentricTreeHierarchy &hierarchy;
        const LabellingT &ttl;
    };

    using TruncatedVertexUpwardSearch = DagShortestPaths<SearchGraphT, BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>, PruneSearchAtUntruncatedVertices<true>>;
    using TruncatedVertexDownwardSearch = DagShortestPaths<SearchGraphT, BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>, PruneSearchAtUntruncatedVertices<false>>;

public:
    TTLQuery(const BalancedTopologyCentricTreeHierarchy &hierarchy,
             const SearchGraphT &upGraph,
             const SearchGraphT &downGraph,
             int const * const upWeights,
             int const * const downWeights,
             const LabellingT &ttl)
            : hierarchy(hierarchy), upGraph(upGraph), downGraph(downGraph), ttl(ttl),
              buildUpLabelSearch(upGraph, upWeights,
                                 {tempUpLabel, upTruncatedSearchSpace, hierarchy, ttl}),
              buildDownLabelSearch(downGraph, downWeights,
                                   {tempDownLabel, downTruncatedSearchSpace, hierarchy, ttl}) {
    }

    // Expects ranks in the underlying separator decomposition order as inputs.
    void run(const int32_t s, const int32_t t) {
        lastDistance = INFTY;
        lastMeetingHubIdx = INVALID_INDEX;
        minCCHDistMeetingVertex = INVALID_VERTEX;
        lastS = s;
        lastT = t;
        const auto lch = hierarchy.getLowestCommonHub(s, t);

        if (!hierarchy.isVertexTruncated(s) && !hierarchy.isVertexTruncated(t)) {
            const auto sUpLabel = ttl.upLabel(s);
            const auto tDownLabel = ttl.downLabel(t);
            computeMinDistanceInLabels(sUpLabel, tDownLabel, lch);
        } else if (!hierarchy.isVertexTruncated(s)) {
            const auto sUpLabel = ttl.upLabel(s);
            buildTempDownLabel(t, lch);
            computeMinDistanceInLabels(sUpLabel, tempDownLabel, lch);
        } else if (!hierarchy.isVertexTruncated(t)) {
            buildTempUpLabel(s, lch);
            const auto tDownLabel = ttl.downLabel(t);
            computeMinDistanceInLabels(tempUpLabel, tDownLabel, lch);
        } else {
            buildTempUpLabel(s, lch);
            buildTempDownLabel(t, lch);
            computeMinDistanceInLabels(tempUpLabel, tempDownLabel, lch);

            // If both s and t are truncated, they may be within the same truncated separator subtree, which we can
            // check by comparing their number of hubs with the lch.
            // In this case, their shortest distance may not use a vertex that is high enough in the hierarchy
            // to be non-truncated. In this case, the shortest path is found using just the distances found during
            // the two elimination tree searches in the CCH. We only have to consider the truncated vertices in
            // the search space of the higher ranked vertex between s and t.
            if (hierarchy.getNumHubs(s) == lch && hierarchy.getNumHubs(t) == lch) {
                const auto& searchSpace = s > t? upTruncatedSearchSpace : downTruncatedSearchSpace;
                for (const auto& v : searchSpace) {
                    const auto cchDist = buildUpLabelSearch.getDistance(v) + buildDownLabelSearch.getDistance(v);
                    if (cchDist < lastDistance) {
                        lastDistance = cchDist;
                        minCCHDistMeetingVertex = v;
                    }
                }
            }
        }

    }

    // Returns the distance of the last query.
    int32_t getDistance() const {
        return lastDistance;
    }

    // Returns the CCH edges in the upward graph on the up segment of the up-down path (in reverse order to conform to
    // default orientation in graph-traversal-based searches).
    template<bool hasPathEdges = LabellingT::StoresPathPointers, std::enable_if_t<hasPathEdges>...>
    const std::vector<int32_t> &getUpEdgePath() {
        lastUpPath.clear();

        // If the best distance was found using not labels but the CCH directly, return the path to the meeting
        // vertex in the CCH:
        if (minCCHDistMeetingVertex != INVALID_VERTEX) {
            lastUpPath = buildUpLabelSearch.getReverseEdgePath(minCCHDistMeetingVertex);
            return lastUpPath;
        }

        int32_t v;
        if (hierarchy.isVertexTruncated(lastS)) {
            // If s was truncated, get path to access vertex, i.e., the first non-truncated vertex used on the up path,
            // using parent pointers in topo search.
            const auto accVertex = tempUpLabel.accessVertex(lastMeetingHubIdx);
            lastUpPath = buildUpLabelSearch.getReverseEdgePath(accVertex);
            std::reverse(lastUpPath.begin(), lastUpPath.end());
            v = accVertex;
        } else {
            // If s was not truncated, the access vertex is simply s itself.
            v = lastS;
        }

        // Enumerate path from access vertex to meeting hub using path edges stored in labels.
        int32_t e;
        while ((e = ttl.upPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            KASSERT(e >= 0 && e < upGraph.numEdges());
            lastUpPath.push_back(e);
            v = upGraph.edgeHead(e);
        }
        std::reverse(lastUpPath.begin(), lastUpPath.end());
        return lastUpPath;
    }

    // Returns the CCH edges in the upward graph on the down segment of the up-down path.
    template<bool hasPathEdges = LabellingT::StoresPathPointers, std::enable_if_t<hasPathEdges>...>
    const std::vector<int32_t> &getDownEdgePath() {
        lastDownPath.clear();

        // If the best distance was found using not labels but the CCH directly, return the path to the meeting
        // vertex in the CCH:
        if (minCCHDistMeetingVertex != INVALID_VERTEX) {
            lastDownPath = buildDownLabelSearch.getReverseEdgePath(minCCHDistMeetingVertex);
            return lastDownPath;
        }

        int32_t v;
        if (hierarchy.isVertexTruncated(lastT)) {
            // If t was truncated, get path to access vertex, i.e., the first non-truncated vertex used on the down path,
            // using parent pointers in topo search.
            const auto accVertex = tempDownLabel.accessVertex(lastMeetingHubIdx);
            lastDownPath = buildDownLabelSearch.getReverseEdgePath(accVertex);
            std::reverse(lastDownPath.begin(), lastDownPath.end());
            v = accVertex;
        } else {
            // If t was not truncated, the access vertex is simply t itself.
            v = lastT;
        }

        // Enumerate path from access vertex to meeting hub using path edges stored in labels.
        int32_t e;
        while ((e = ttl.downPathEdge(v, lastMeetingHubIdx)) != INVALID_EDGE) {
            KASSERT(e >= 0 && e < downGraph.numEdges());
            lastDownPath.push_back(e);
            v = downGraph.edgeHead(e);
        }
        std::reverse(lastDownPath.begin(), lastDownPath.end());
        return lastDownPath;
    }


private:

    // Populates tempUpLabel with up label for truncated vertex v by running topological upwards search from v.
    // Whenever search runs into a non-truncated vertex w, the label of w is used to update tempUpLabel and the
    // search is pruned.
    // Populates distances only for hubs 0..lch (inclusive).
    void buildTempUpLabel(const int32_t v, const uint32_t lch) {
        KASSERT(lch <= hierarchy.getNumHubs(v));
        tempUpLabel.init(lch);
        upTruncatedSearchSpace.clear();
        buildUpLabelSearch.run(v);
    }

    // Populates tempDownLabel with down label for truncated vertex v by running topological reverse-downwards search
    // from v. Whenever search runs into a non-truncated vertex w, the label of w is used to update tempDownLabel and
    // the search is pruned.
    // Populates distances only for hubs 0..lch (inclusive).
    void buildTempDownLabel(const int32_t v, const uint32_t lch) {
        KASSERT(lch <= hierarchy.getNumHubs(v));
        tempDownLabel.init(lch);
        downTruncatedSearchSpace.clear();
        buildDownLabelSearch.run(v);
    }

    static_assert(logK != 1, "logK==1 is not supported for TTLMetric customization. "
                             "Use logK >= 2 for SIMD or logK = 0 for no SIMD.");
    static constexpr uint32_t K = 1 << logK;

    // Vectors of multiple data elements for use with SIMD instructions.
    using BooleanVector = std::conditional_t<logK == 2, Vec4ib, std::conditional_t<logK == 3, Vec8ib, Vec16ib>>;
    using IntegerVector = std::conditional_t<logK == 2, Vec4i, std::conditional_t<logK == 3, Vec8i, Vec16i>>;

    template<typename UpLabel, typename DownLabel>
    inline void computeMinDistanceInLabels(const UpLabel& up, const DownLabel& down, const uint32_t lowestCommonHub) {


        int32_t b = 0;
        if constexpr (logK >= 2) {
            // Compute minimum values using SIMD instructions for full SIMD vectors:
            const uint32_t numFullBlocks = lowestCommonHub / K;
            IntegerVector dUp;
            IntegerVector dDown;
            IntegerVector dMin(lastDistance);
            IntegerVector minBlock;
            BooleanVector improved;
            for (; b < numFullBlocks; ++b) {
                dUp.load(&up.dist(b * K));
                dDown.load(&down.dist(b * K));
                dUp += dDown;
                improved = (dUp < dMin);
                dMin = select(improved, dUp, dMin);
                minBlock = select(improved, b, minBlock);
            }
            const auto minInBlocks = horizontal_min(dMin);
            if (minInBlocks < lastDistance) {
                lastDistance = minInBlocks;
                const auto idxWithinBlock = horizontal_find_first(dMin == minInBlocks);
                lastMeetingHubIdx = minBlock[idxWithinBlock] * K + idxWithinBlock;
            }
        }

        // For remaining partially filled block do not use SIMD:
        for (uint32_t i = b * K; i < lowestCommonHub; ++i) {
            const int32_t distForHubI = up.dist(i) + down.dist(i);
            if (distForHubI < lastDistance) {
                lastDistance = distForHubI;
                lastMeetingHubIdx = i;
            }
        }
//        for (uint32_t i = 0; i < lowestCommonHub; ++i) {
//            const int32_t distForHubI = up.dist(i) + down.dist(i);
//            if (distForHubI < lastDistance) {
//                lastDistance = distForHubI;
//                lastMeetingHubIdx = i;
//            }
//        }
    }

    const BalancedTopologyCentricTreeHierarchy &hierarchy;
    const SearchGraphT &upGraph;
    const SearchGraphT &downGraph;
    const LabellingT &ttl;

    int32_t lastS;
    int32_t lastT;
    int32_t lastDistance;
    uint32_t lastMeetingHubIdx;
    int32_t minCCHDistMeetingVertex; // if shortest path does not use labels but CCH distances, this is set to meeting vertex in CCH
    std::vector<int32_t> lastUpPath;
    std::vector<int32_t> lastDownPath;

    TemporaryLabel tempUpLabel;
    TemporaryLabel tempDownLabel;
    TruncatedVertexUpwardSearch buildUpLabelSearch;
    TruncatedVertexDownwardSearch buildDownLabelSearch;
    std::vector<int32_t> upTruncatedSearchSpace; // All truncated vertices that the last up search visited
    std::vector<int32_t> downTruncatedSearchSpace; // All truncated vertices that the last down search visited

};

