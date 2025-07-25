#include <ctlsa/road_network.h>
#include <ctlsa/util.h>

#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <bitset>
#include <unordered_set>
#include <atomic>
#include <cstring>
#include <random>

using namespace std;

#define DEBUG(X) //cerr << X << endl

// algorithm config
//#define CUT_REPEAT 3 // repeat whole cut computation multiple times (with random starting points for rough partition) and pick best result
#define MULTI_CUT // extract two different min-cuts from max-flow and pick more balanced result
static const bool weighted_furthest = false; // use edge weights for finding distant nodes during rough partitioning
static const bool weighted_diff = false; // use edge weights for computing rough partition

namespace ctlsa::road_network {

    static const NodeID NO_NODE = 0; // null value equivalent for integers identifying nodes
    static const SubgraphID NO_SUBGRAPH = 0; // used to indicate that node does not belong to any active subgraph
    static const uint16_t MAX_CUT_LEVEL = 58; // maximum height of decomposition tree; 58 bits to store binary path, plus 6 bits to store path length = 64 bit integer

// profiling
#ifndef NPROFILE
    static atomic<double> t_partition, t_label, t_shortcut;
#define START_TIMER util::start_timer()
#define STOP_TIMER(var) var += util::stop_timer()
#else
#define START_TIMER
#define STOP_TIMER(var)
#endif

// progress of 0 resets counter
    static bool log_progress_on = false;

    void log_progress(size_t p, ostream &os = cout) {
        static const size_t P_DIFF = 1000000L;
        static size_t progress = 0;
        if (log_progress_on) {
            size_t old_log = progress / P_DIFF;
            if (p == 0) {
                // terminate progress line & reset
                if (old_log > 0)
                    os << endl;
                progress = 0;
                return;
            }
            progress += p;
            size_t new_log = progress / P_DIFF;
            if (old_log < new_log) {
                for (size_t i = old_log; i < new_log; i++)
                    os << '.';
                os << flush;
            }
        } else
            progress += p;
    }

// half-matrix index for storing half-matrix in flat vector
    static size_t hmi(size_t a, size_t b) {
        assert(a != b);
        return a < b ? (b * (b - 1) >> 1) + a : (a * (a - 1) >> 1) + b;
    }

// offset by cut level
    uint16_t get_offset(const uint16_t *dist_index, size_t cut_level) {
        return cut_level ? dist_index[cut_level - 1] : 0;
    }

//--------------------------- CutIndex ------------------------------

    CutIndex::CutIndex() : truncated(false), partition(0), cut_level(0) {
#ifdef PRUNING
        pruning_2hop = pruning_3hop = pruning_tail = 0;
#endif
    }

#ifdef PRUNING
    void CutIndex::prune_tail()
    {
        assert(is_consistent(true));
        assert(dist_index.back() == distances.size());
        assert(distances.size() > 0);
        // cut_level may not be set yet
        size_t cl = dist_index.size() - 1;
        // only prune latest cut
        size_t last_unpruned = get_offset(&dist_index[0], cl);
        // nothing to prune for empty cuts
        if (last_unpruned == distances.size())
            return;
        // first node must never be pruned
        assert(distances[last_unpruned] & 1);
        // fix distances and recall last unprunded label
        for (size_t i = last_unpruned; i < distances.size(); i++)
        {
            if (distances[i] & 1)
                last_unpruned = i;
            distances[i] >>= 1;
        }
        size_t new_size = last_unpruned + 1;
        assert(new_size <= distances.size());
        if (new_size < distances.size())
        {
            pruning_tail += distances.size() - new_size;
            distances.resize(new_size);
            dist_index.back() = new_size;
            DEBUG("pruned tail: " << *this);
        }
    }
#endif

    bool CutIndex::is_consistent(bool partial) const {
        const uint64_t one = 1;
        if (cut_level > MAX_CUT_LEVEL) {
            cerr << "cut_level=" << (int) cut_level << endl;
            return false;
        }
        if (!partial && partition >= (one << cut_level)) {
            cerr << "partition=" << partition << " for cut_level=" << (int) cut_level << endl;
            return false;
        }
        if (!partial && dist_index.size() != cut_level + one) {
            cerr << "dist_index.size()=" << dist_index.size() << " for cut_level=" << (int) cut_level << endl;
            return false;
        }
        if (!is_sorted(dist_index.cbegin(), dist_index.cend())) {
            cerr << "unsorted dist_index: " << dist_index << endl;
            return false;
        }
        return true;
    }

    bool CutIndex::empty() const {
        return dist_index.empty();
    }

//--------------------------- PBV -----------------------------------

    namespace PBV {

        uint64_t from(uint64_t bits, uint16_t length) {
            if (length == 0)
                return 0;
            return (bits << (64 - length) >> (58 - length)) | length;
        }

        uint64_t partition(uint64_t bv) {
            // cutlevel is stored in lowest 6 bits
            return bv >> 6;
        }

        uint16_t cut_level(uint64_t bv) {
            // cutlevel is stored in lowest 6 bits
            return bv & 63ul;
        }

        uint16_t lca_level(uint64_t bv1, uint64_t bv2) {
            // find lowest level at which partitions differ
            uint16_t lca_level = min(cut_level(bv1), cut_level(bv2));
            uint64_t p1 = partition(bv1), p2 = partition(bv2);
            if (p1 != p2) {
                uint16_t diff_level = __builtin_ctzll(p1 ^ p2); // count trailing zeros
                if (diff_level < lca_level)
                    lca_level = diff_level;
            }
            return lca_level;
        }

        uint64_t lca(uint64_t bv1, uint64_t bv2) {
            uint64_t cut_level = lca_level(bv1, bv2);
            // shifting by 64 does not work
            if (cut_level == 0)
                return 0;
            return (bv1 >> 6) << (64 - cut_level) >> (58 - cut_level) | cut_level;
        }

        bool is_ancestor(uint64_t bv_ancestor, uint64_t bv_descendant) {
            uint16_t cla = cut_level(bv_ancestor), cld = cut_level(bv_descendant);
            // shifting by 64 does not work, so need to check for cla == 0
            return cla == 0 || (cla <= cld && (bv_ancestor ^ bv_descendant) >> 6 << (64 - cla) == 0);
        }

    }

//--------------------------- FlatCutIndex --------------------------

// helper function for memory alignment
    template<typename T>
    size_t aligned(size_t size);

    template<>
    size_t aligned<uint32_t>(size_t size) {
        size_t mod = size & 3ul;
        return mod ? size + (4 - mod) : size;
    }

    FlatCutIndex::FlatCutIndex() : data(nullptr) {
    }

    FlatCutIndex::FlatCutIndex(const CutIndex &ci) {
        assert(ci.is_consistent());
        // allocate memory for partition bitvector, distance_offset, label_count, dist_index, distances, and paths
        // distance_offset is redundant to speed up distance pointer calculation, label_count permits truncated labels to be stored
        size_t distance_offset =
                sizeof(uint64_t) + 2 * sizeof(uint16_t) + aligned<distance_t>(ci.dist_index.size() * sizeof(uint16_t));
        size_t data_size =
                distance_offset + ci.distances.size() * sizeof(distance_t) + ci.paths.size() * sizeof(CHNeighbor);
        data = (char *) calloc(data_size, 1);
        // copy partition bitvector, distance_offset, label_count, dist_index, distances, and paths into data
        *partition_bitvector() = PBV::from(ci.partition, ci.cut_level);
        *_distance_offset() = distance_offset;
        *_label_count() = ci.distances.size();
        memcpy(dist_index(), &ci.dist_index[0], ci.dist_index.size() * sizeof(uint16_t));
        memcpy(distances(), &ci.distances[0], ci.distances.size() * sizeof(distance_t));
        memcpy(paths(), &ci.paths[0], ci.paths.size() * sizeof(CHNeighbor));
    }

    bool FlatCutIndex::operator==(FlatCutIndex other) const {
        return data == other.data;
    }

    uint64_t *FlatCutIndex::partition_bitvector() {
        assert(!empty());
        return (uint64_t *) data;
    }

    const uint64_t *FlatCutIndex::partition_bitvector() const {
        assert(!empty());
        return (uint64_t *) data;
    }

    uint16_t *FlatCutIndex::_distance_offset() {
        assert(!empty());
        return (uint16_t *) (data + sizeof(uint64_t));
    }

    const uint16_t *FlatCutIndex::_distance_offset() const {
        assert(!empty());
        return (uint16_t *) (data + sizeof(uint64_t));
    }

    uint16_t *FlatCutIndex::_label_count() {
        assert(!empty());
        return (uint16_t *) (data + sizeof(uint64_t)) + 1;
    }

    const uint16_t *FlatCutIndex::_label_count() const {
        assert(!empty());
        return (uint16_t *) (data + sizeof(uint64_t)) + 1;
    }

    uint16_t *FlatCutIndex::dist_index() {
        assert(!empty());
        return (uint16_t *) (data + sizeof(uint64_t)) + 2;
    }

    const uint16_t *FlatCutIndex::dist_index() const {
        assert(!empty());
        return (uint16_t *) (data + sizeof(uint64_t)) + 2;
    }

    distance_t *FlatCutIndex::distances() {
        assert(!empty());
        return (distance_t *) (data + *_distance_offset());
    }

    const distance_t *FlatCutIndex::distances() const {
        assert(!empty());
        return (distance_t *) (data + *_distance_offset());
    }

    CHNeighbor *FlatCutIndex::paths() {
        assert(!empty());
        return (CHNeighbor *) (data + *_distance_offset() + *_label_count() * sizeof(distance_t));
    }

    const CHNeighbor *FlatCutIndex::paths() const {
        assert(!empty());
        return (CHNeighbor *) (data + *_distance_offset() + *_label_count() * sizeof(distance_t));
    }

    uint64_t FlatCutIndex::partition() const {
        return PBV::partition(*partition_bitvector());
    }

    uint16_t FlatCutIndex::cut_level() const {
        return PBV::cut_level(*partition_bitvector());
    }

    size_t FlatCutIndex::size() const {
        return *_distance_offset() + *_label_count() * sizeof(distance_t);
    }

#ifndef PRUNING

    size_t FlatCutIndex::ancestor_count() const {
        return dist_index()[cut_level()];
    }

#endif

    size_t FlatCutIndex::label_count() const {
        return *_label_count();
    }

    size_t FlatCutIndex::cut_size(size_t cl) const {
        return cl == 0 ? *dist_index() : dist_index()[cl] - dist_index()[cl - 1];
    }

    size_t FlatCutIndex::bottom_cut_size() const {
        return cut_size(cut_level());
    }

    bool FlatCutIndex::empty() const {
        return data == nullptr;
    }

    const distance_t *FlatCutIndex::cl_begin(size_t cl) const {
        uint16_t offset = get_offset(dist_index(), cl);
        return distances() + min(*_label_count(), offset);
    }

    const distance_t *FlatCutIndex::cl_end(size_t cl) const {
        uint16_t offset = dist_index()[cl];
        return distances() + min(*_label_count(), offset);
    }

    const CHNeighbor *FlatCutIndex::cl_begin_paths(size_t cl) const {
        uint16_t offset = get_offset(dist_index(), cl);
        return paths() + min(*_label_count(), offset);
    }

    const CHNeighbor *FlatCutIndex::cl_end_paths(size_t cl) const {
        uint16_t offset = dist_index()[cl];
        return paths() + min(*_label_count(), offset);
    }

    vector <vector<distance_t>> FlatCutIndex::unflatten() const {
        vector<vector<distance_t>> labels;
        if (!empty()) {
            for (size_t cl = 0; cl <= cut_level(); cl++) {
                vector<distance_t> cut_labels;
                for (const distance_t *l = cl_begin(cl); l != cl_end(cl); l++)
                    cut_labels.push_back(*l);
                labels.push_back(cut_labels);
            }
        }
        return labels;
    }

    std::vector<std::vector<CHNeighbor>> FlatCutIndex::unflatten_paths() const {
        vector<vector<CHNeighbor>> labels;
        if (!empty()) {
            for (size_t cl = 0; cl <= cut_level(); cl++) {
                vector<CHNeighbor> cut_labels;
                for (const CHNeighbor *l = cl_begin_paths(cl); l != cl_end_paths(cl); l++)
                    cut_labels.push_back(*l);
                labels.push_back(cut_labels);
            }
        }
        return labels;
    }

//--------------------------- ContractionLabel ----------------------

    ContractionLabel::ContractionLabel() : distance_offset(0), root(NO_NODE) {
    }

    size_t ContractionLabel::size() const {
        size_t total = sizeof(ContractionLabel);
        // only count index data if owned
        if (distance_offset == 0)
            total += cut_index.size();
        return total;
    }

//--------------------------- ContractionIndex ----------------------

    template<typename T>
    static void clear_and_shrink(vector <T> &v) {
        v.clear();
        v.shrink_to_fit();
    }

    ContractionIndex::ContractionIndex(vector <CutIndex> &ci, vector <Neighbor> &closest) {
        assert(ci.size() == closest.size());
        labels.resize(ci.size());
        // handle core nodes
        for (NodeID node = 1; node < closest.size(); node++) {
            if (closest[node].node == node) {
                assert(closest[node].distance == 0);
                labels[node].cut_index = FlatCutIndex(ci[node]);
            }
            // conserve memory
            clear_and_shrink(ci[node].dist_index);
            clear_and_shrink(ci[node].distances);
        }
        // handle periferal nodes
        for (NodeID node = 1; node < closest.size(); node++) {
            Neighbor n = closest[node];
            // isolated nodes got removed (n.node == NO_NODE)
            if (n.node != node && n.node != NO_NODE) {
                assert(n.distance > 0);
                // find root & distance
                NodeID root = n.node;
                distance_t root_dist = n.distance;
                while (closest[root].node != root) {
                    root_dist += closest[root].distance;
                    root = closest[root].node;
                }
                // copy index
                assert(!labels[root].cut_index.empty());
                labels[node].cut_index = labels[root].cut_index;
                labels[node].root = root;
                labels[node].distance_offset = root_dist;
            }
        }
        clear_and_shrink(ci);
        clear_and_shrink(closest);
    }

    ContractionIndex::ContractionIndex(std::vector<CutIndex> &ci) {
        labels.resize(ci.size());
        for (NodeID node = 1; node < ci.size(); node++)
            if (!ci[node].empty()) {
                labels[node].cut_index = FlatCutIndex(ci[node]);
                // conserve memory
                clear_and_shrink(ci[node].dist_index);
                clear_and_shrink(ci[node].distances);
            }
        clear_and_shrink(ci);
    }

    ContractionIndex::~ContractionIndex() {
        for (NodeID node = 1; node < labels.size(); node++)
            // not all labels own their cut index data
            if (!labels[node].cut_index.empty() && labels[node].distance_offset == 0)
                free(labels[node].cut_index.data);
    }

    size_t ContractionIndex::get_hoplinks(NodeID v, NodeID w) const {
        FlatCutIndex cv = labels[v].cut_index, cw = labels[w].cut_index;
        if (cv == cw)
            return 0;
        return get_hoplinks(cv, cw);
    }

    double ContractionIndex::avg_hoplinks(const std::vector<std::pair<NodeID, NodeID>> &queries) const {
        size_t hop_count = 0;
        for (pair<NodeID, NodeID> q: queries)
            hop_count += get_hoplinks(q.first, q.second);
        return hop_count / (double) queries.size();
    }

#ifndef PRUNING

    size_t ContractionIndex::common_ancestor_count(NodeID v, NodeID w) const {
        FlatCutIndex cv = labels[v].cut_index, cw = labels[w].cut_index;
        if (cv == cw)
            return 0;
        uint16_t lca_level = PBV::lca_level(*cv.partition_bitvector(), *cw.partition_bitvector());
        return min(cv.dist_index()[lca_level], cw.dist_index()[lca_level]);
    }

#endif

    size_t ContractionIndex::get_cut_level_hoplinks(FlatCutIndex a, FlatCutIndex b, size_t cut_level) {
        return min(a.cut_size(cut_level), b.cut_size(cut_level));
    }

    bool ContractionIndex::is_contracted(NodeID node) const {
        return labels[node].root != NO_NODE;
    }

    size_t ContractionIndex::uncontracted_count() const {
        size_t total = 0;
        for (NodeID node = 1; node < labels.size(); node++)
            if (!is_contracted(node))
                total++;
        return total;
    }

    bool ContractionIndex::in_partition_subgraph(NodeID node, uint64_t partition_bitvector) const {
        return !is_contracted(node) &&
               PBV::is_ancestor(partition_bitvector, *labels[node].cut_index.partition_bitvector());
    }

    uint16_t ContractionIndex::dist_index(NodeID node) const {
        FlatCutIndex const &ci = labels[node].cut_index;
        uint16_t index = get_offset(ci.dist_index(), ci.cut_level());
        while (ci.distances()[index] != 0)
            index++;
        return index;
    }

    ContractionLabel ContractionIndex::get_contraction_label(NodeID v) const {
        return labels[v];
    }

    void ContractionIndex::update_distance_offset(NodeID n, distance_t d) {
        labels[n].distance_offset = d;
    }

    size_t ContractionIndex::get_hoplinks(FlatCutIndex a, FlatCutIndex b) {
        // find lowest level at which partitions differ
        size_t cut_level = min(a.cut_level(), b.cut_level());
        uint64_t pa = a.partition(), pb = b.partition();
        if (pa != pb) {
            size_t diff_level = __builtin_ctzll(pa ^ pb); // count trailing zeros
            if (diff_level < cut_level)
                cut_level = diff_level;
        }
#ifdef NO_SHORTCUTS
        size_t hoplinks = 0;
        for (size_t cl = 0; cl <= cut_level; cl++)
            hoplinks += get_cut_level_hoplinks(a, b, cl);
        return hoplinks;
#else
        return get_cut_level_hoplinks(a, b, cut_level);
#endif
    }

    size_t ContractionIndex::size() const {
        size_t total = 0;
        for (NodeID node = 1; node < labels.size(); node++) {
            // skip isolated nodes (subgraph)
            if (!labels[node].cut_index.empty())
                total += labels[node].size();
        }
        return total;
    }

    double ContractionIndex::avg_cut_size() const {
        double cut_sum = 0, label_count = 0;
        for (NodeID node = 1; node < labels.size(); node++)
            if (!labels[node].cut_index.empty()) {
                cut_sum += labels[node].cut_index.cut_level() + 1;
                label_count += labels[node].cut_index.label_count();
                // adjust for label pruning
                label_count += labels[node].cut_index.bottom_cut_size() + 1;
            }
        return label_count / max(1.0, cut_sum);
    }

    size_t ContractionIndex::max_cut_size() const {
        size_t max_cut = 0;
        for (NodeID node = 1; node < labels.size(); node++)
            if (!labels[node].cut_index.empty())
                max_cut = max(max_cut, 1 + labels[node].cut_index.bottom_cut_size());
        return max_cut;
    }

    size_t ContractionIndex::height() const {
        uint16_t max_cut_level = 0;
        for (NodeID node = 1; node < labels.size(); node++)
            if (!labels[node].cut_index.empty())
                max_cut_level = max(max_cut_level, labels[node].cut_index.cut_level());
        return max_cut_level;
    }

    size_t ContractionIndex::max_label_count() const {
        size_t max_label_count = 0;
        for (NodeID node = 1; node < labels.size(); node++)
            if (!labels[node].cut_index.empty())
                max_label_count = max(max_label_count, labels[node].cut_index.label_count());
        return max_label_count;
    }

    size_t ContractionIndex::label_count() const {
        size_t total = 0;
        for (NodeID node = 1; node < labels.size(); node++)
            if (!labels[node].cut_index.empty() && labels[node].distance_offset == 0)
                total += labels[node].cut_index.label_count();
        return total;
    }

    size_t ContractionIndex::non_empty_cuts() const {
        size_t total = 0;
        for (NodeID node = 1; node < labels.size(); node++) {
            if (is_contracted(node))
                continue;
            // count nodes that come first within their cut
            FlatCutIndex const &ci = labels[node].cut_index;
            if (ci.distances()[get_offset(ci.dist_index(), ci.cut_level())] == 0)
                total++;
        }
        return total;
    }

    pair <NodeID, NodeID> ContractionIndex::random_query() const {
        assert(labels.size() > 1);
        NodeID node_count = labels.size() - 1;
        NodeID a = 1 + rand() % node_count;
        NodeID b = 1 + rand() % node_count;
        return make_pair(a, b);
    }

    void ContractionIndex::write(ostream &os) const {
        size_t node_count = labels.size() - 1;
        os.write((char *) &node_count, sizeof(size_t));
        for (NodeID node = 1; node < labels.size(); node++) {
            ContractionLabel cl = labels[node];
            os.write((char *) &cl.distance_offset, sizeof(distance_t));
            if (cl.distance_offset == 0) {
                size_t data_size = cl.cut_index.size();
                os.write((char *) &data_size, sizeof(size_t));
                os.write(cl.cut_index.data, data_size);
            } else
                os.write((char *) &cl.root, sizeof(NodeID));
        }
    }

    void ContractionIndex::write_json(std::ostream &os) const {
        ListFormat lf = get_list_format();
        set_list_format(ListFormat::plain);
        // print json
        os << '{' << endl;
        for (NodeID node = 1; node < labels.size(); node++) {
            os << node << ":";
            ContractionLabel cl = labels[node];
            if (cl.distance_offset == 0) {
                vector<vector<distance_t>> distances = cl.cut_index.unflatten();
                os << distances;
                if (!distances.empty())
                    os << "," << endl << "p" << node << ":" << cl.cut_index.unflatten_paths();
            } else
                os << "{\"r\":" << cl.root << ",\"d\":" << cl.distance_offset << "}";
            os << (node == labels.size() - 1 ? "" : ",") << endl;
        }
        os << '}' << endl;
        // reset formatting
        set_list_format(lf);
    }

    ContractionIndex::ContractionIndex(istream &is) {
        // read index data
        size_t node_count = 0;
        is.read((char *) &node_count, sizeof(size_t));
        labels.resize(node_count + 1);
        for (NodeID node = 1; node < labels.size(); node++) {
            ContractionLabel &cl = labels[node];
            is.read((char *) &cl.distance_offset, sizeof(distance_t));
            if (cl.distance_offset == 0) {
                size_t data_size = 0;
                is.read((char *) &data_size, sizeof(size_t));
                cl.cut_index.data = (char *) malloc(data_size);
                is.read(cl.cut_index.data, data_size);
            } else
                is.read((char *) &cl.root, sizeof(NodeID));
        }
        // fix data references
        for (NodeID node = 1; node < labels.size(); node++) {
            ContractionLabel &cl = labels[node];
            if (cl.distance_offset != 0) {
                NodeID root = cl.root;
                cl.cut_index = labels[root].cut_index;
            }
        }
    }

//--------------------------- ContractionHierarchy ------------------

    std::ostream &operator<<(std::ostream &os, const CHNode &ch) {
        return os << "CH(" << ch.dist_index << "," << ch.up_neighbors << "," << ch.down_neighbors << ")";
    }

    std::ostream &operator<<(std::ostream &os, const ContractionHierarchy &ch) {
        return os << ch.nodes;
    }

    size_t ContractionHierarchy::edge_count() const {
        size_t total = 0;
        for (CHNode const &node: nodes)
            total += node.up_neighbors.size();
        return total;
    }

//--------------------------- Graph ---------------------------------

    SubgraphID next_subgraph_id(bool reset) {
        static atomic<SubgraphID> next_id = 1;
        if (reset)
            next_id = 1;
        return next_id++;
    }

    Neighbor::Neighbor(NodeID node, distance_t distance) : node(node), distance(distance) {
    }

    CHNeighbor::CHNeighbor(NodeID neighbor, distance_t distance, path_data p) : neighbor(neighbor), distance(distance),
                                                                                p(p) {
    }

    composite_shortcut::composite_shortcut() {
        _no_node = NO_NODE;
        triangle_node = NO_NODE;
        lower = nullptr;
        upper = nullptr;
    }

    path_data::path_data() {
        intermediate.fill(NO_NODE);
        new(&cs) composite_shortcut();
    }

    uint32_t path_data::intermediate_count() {
        if (MAX_VALLEY_PATH_LENGTH > 0 && intermediate[0] == NO_NODE)
            return cs.triangle_node == NO_NODE ? 0 : infinity;
        uint32_t count = 1;
        while (count < MAX_VALLEY_PATH_LENGTH && intermediate[count] != NO_NODE)
            count++;
        return count;
    }

    QNode::QNode(CHNeighbor shortcut_used, uint16_t shortcut_start, NodeID node_id, distance_t distance)
            : shortcut_used(shortcut_used), shortcut_start(shortcut_start), node_id(node_id), distance(distance) {
    }

    bool Neighbor::operator<(const Neighbor &other) const {
        return node < other.node;
    }

    Node::Node(SubgraphID subgraph_id) : subgraph_id(subgraph_id) {
        distance = outcopy_distance = 0;
        inflow = outflow = NO_NODE;
        landmark_level = 0;
    }

    Node &MultiThreadNodeData::operator[](size_type pos) {
        if (pos == Graph::s)
            return s_data;
        if (pos == Graph::t)
            return t_data;
        return vector::operator[](pos);
    }

    const Node &MultiThreadNodeData::operator[](size_type pos) const {
        if (pos == Graph::s)
            return s_data;
        if (pos == Graph::t)
            return t_data;
        return vector::operator[](pos);
    }

    void MultiThreadNodeData::normalize() {
        vector::operator[](Graph::s) = s_data;
        vector::operator[](Graph::t) = t_data;
    }

    double Partition::rating() const {
        size_t l = left.size(), r = right.size(), c = cut.size();
        return min(l, r) / (c * c + 1.0);
    }

    Edge::Edge(NodeID a, NodeID b, distance_t d) : a(a), b(b), d(d) {
    }

    bool Edge::operator<(Edge other) const {
        return a < other.a
               || (a == other.a && b < other.b)
               || (a == other.a && b == other.b && d < other.d);
    }

    int32_t DiffData::diff() const {
        return static_cast<int32_t>(dist_a) - static_cast<int32_t>(dist_b);
    }

    distance_t DiffData::min() const {
        return std::min(dist_a, dist_b);
    }

    DiffData::DiffData(NodeID node, distance_t dist_a, distance_t dist_b) : node(node), dist_a(dist_a), dist_b(dist_b) {
    }

    bool DiffData::cmp_diff(DiffData x, DiffData y) {
        return x.diff() < y.diff();
    }

// definition of static members
    thread_local Node MultiThreadNodeData::s_data(NO_SUBGRAPH), MultiThreadNodeData::t_data(NO_SUBGRAPH);
#ifdef MULTI_THREAD
    MultiThreadNodeData Graph::node_data;
    size_t Graph::thread_threshold;
#else
    vector <Node> Graph::node_data;
#endif
    NodeID Graph::s, Graph::t;

    void Graph::show_progress(bool state) {
        log_progress_on = state;
    }

    bool Graph::contains(NodeID node) const {
        return node_data[node].subgraph_id == subgraph_id;
    }

    Graph::Graph(size_t node_count) {
        subgraph_id = next_subgraph_id(true);
        node_data.clear();
        resize(node_count);
        CHECK_CONSISTENT;
    }

    Graph::Graph(size_t node_count, const vector <Edge> &edges) : Graph(node_count) {
        for (Edge e: edges)
            add_edge(e.a, e.b, e.d, true);
    }

    void Graph::resize(size_t node_count) {
        assert(nodes.empty());
        // node numbering starts from 1, and we reserve two additional nodes for s & t
        node_data.clear();
        node_data.resize(node_count + 3, Node(subgraph_id));
        s = node_count + 1;
        t = node_count + 2;
        node_data[0].subgraph_id = node_data[s].subgraph_id = node_data[t].subgraph_id = NO_SUBGRAPH;
        nodes.reserve(node_count);
        for (NodeID node = 1; node <= node_count; node++)
            nodes.push_back(node);
#ifdef MULTI_THREAD
        thread_threshold = max(node_count / MULTI_THREAD, static_cast<size_t>(1000));
#endif
    }

    void Graph::add_edge(NodeID v, NodeID w, distance_t distance, bool add_reverse) {
        assert(v < node_data.size());
        assert(w < node_data.size());
        assert(distance > 0);
        // check for existing edge
        bool exists = false;
        for (Neighbor &n: node_data[v].neighbors)
            if (n.node == w) {
                exists = true;
                n.distance = min(n.distance, distance);
                break;
            }
        if (!exists)
            node_data[v].neighbors.push_back(Neighbor(w, distance));
        if (add_reverse)
            add_edge(w, v, distance, false);
    }

    void Graph::remove_edge(NodeID v, NodeID w) {
        std::erase_if(node_data[v].neighbors, [w](const Neighbor &n) { return n.node == w; });
        std::erase_if(node_data[w].neighbors, [v](const Neighbor &n) { return n.node == v; });
    }

    void Graph::remove_isolated() {
        unordered_set<NodeID> isolated;
        for (NodeID node: nodes)
            if (degree(node) == 0) {
                isolated.insert(node);
                node_data[node].subgraph_id = NO_SUBGRAPH;
            }
        std::erase_if(nodes, [&isolated](NodeID node) { return isolated.contains(node); });
    }

    void Graph::reset() {
        nodes.clear();
        for (NodeID node = 1; node < node_data.size() - 2; node++) {
            if (!node_data[node].neighbors.empty()) {
                nodes.push_back(node);
                node_data[node].subgraph_id = subgraph_id;
            }
        }
        node_data[s].subgraph_id = NO_SUBGRAPH;
        node_data[t].subgraph_id = NO_SUBGRAPH;
    }

    void Graph::add_node(NodeID v) {
        assert(v < node_data.size());
        nodes.push_back(v);
        node_data[v].subgraph_id = subgraph_id;
    }

    void Graph::remove_nodes(const vector <NodeID> &node_set) {
        util::remove_set(nodes, node_set);
        for (NodeID node: node_set)
            node_data[node].subgraph_id = NO_NODE;
    }

    size_t Graph::node_count() const {
        return nodes.size();
    }

    size_t Graph::edge_count() const {
        size_t ecount = 0;
        for (NodeID node: nodes)
            for (Neighbor n: node_data[node].neighbors)
                if (contains(n.node))
                    ecount++;
        return ecount / 2;
    }

    size_t Graph::degree(NodeID v) const {
        assert(contains(v));
        size_t deg = 0;
        for (Neighbor n: node_data[v].neighbors)
            if (contains(n.node))
                deg++;
        return deg;
    }

    Neighbor Graph::single_neighbor(NodeID v) const {
        assert(contains(v));
        Neighbor neighbor(NO_NODE, 0);
        for (Neighbor n: node_data[v].neighbors)
            if (contains(n.node)) {
                if (neighbor.node == NO_NODE)
                    neighbor = n;
                else
                    return Neighbor(NO_NODE, 0);
            }
        return neighbor;
    }

    size_t Graph::super_node_count() {
        return node_data.size() - 3;
    }

    const vector <NodeID> &Graph::get_nodes() const {
        return nodes;
    }

    void Graph::get_edges(vector <Edge> &edges) const {
        edges.clear();
        for (NodeID a: nodes)
            for (const Neighbor &n: node_data[a].neighbors)
                if (n.node > a && contains(n.node))
                    edges.push_back(Edge(a, n.node, n.distance));
    }

    void Graph::assign_nodes() {
        for (NodeID node: nodes)
            node_data[node].subgraph_id = subgraph_id;
    }

//--------------------------- Graph algorithms ----------------------

// helper struct to enque nodes by distance
    struct SearchNode {
        distance_t distance;
        NodeID node;

        // reversed for min-heap ordering
        bool operator<(const SearchNode &other) const { return distance > other.distance; }

        SearchNode(distance_t distance, NodeID node) : distance(distance), node(node) {}
    };

    void Graph::run_dijkstra(NodeID v) {
        CHECK_CONSISTENT;
        assert(contains(v));
        // init distances
        for (NodeID node: nodes)
            node_data[node].distance = infinity;
        node_data[v].distance = 0;
        // init queue
        priority_queue<SearchNode> q;
        q.push(SearchNode(0, v));
        // dijkstra
        while (!q.empty()) {
            SearchNode next = q.top();
            q.pop();

            for (Neighbor n: node_data[next.node].neighbors) {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // update distance and enque
                distance_t new_dist = next.distance + n.distance;
                if (new_dist < node_data[n.node].distance) {
                    node_data[n.node].distance = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                }
            }
        }
    }

    void Graph::run_bfs(NodeID v) {
        CHECK_CONSISTENT;
        assert(contains(v));
        // init distances
        for (NodeID node: nodes)
            node_data[node].distance = infinity;
        node_data[v].distance = 0;
        // init queue
        queue<NodeID> q;
        q.push(v);
        // BFS
        while (!q.empty()) {
            NodeID next = q.front();
            q.pop();

            distance_t new_dist = node_data[next].distance + 1;
            for (Neighbor n: node_data[next].neighbors) {
                // filter neighbors nodes not belonging to subgraph or already visited
                if (contains(n.node) && node_data[n.node].distance == infinity) {
                    // update distance and enque
                    node_data[n.node].distance = new_dist;
                    q.push(n.node);
                }
            }
        }
    }

// node in flow graph which splits nodes into incoming and outgoing copies
    struct FlowNode {
        NodeID node;
        bool outcopy; // outgoing copy of node?
        FlowNode(NodeID node, bool outcopy) : node(node), outcopy(outcopy) {}
    };

    ostream &operator<<(ostream &os, FlowNode fn) {
        return os << "(" << fn.node << "," << (fn.outcopy ? "T" : "F") << ")";
    }

// helper function
    bool update_distance(distance_t &d, distance_t d_new) {
        if (d > d_new) {
            d = d_new;
            return true;
        }
        return false;
    }

    void Graph::run_flow_bfs_from_s() {
        CHECK_CONSISTENT;
        assert(contains(s) && contains(t));
        // init distances
        for (NodeID node: nodes)
            node_data[node].distance = node_data[node].outcopy_distance = infinity;
        node_data[t].distance = node_data[t].outcopy_distance = 0;
        // init queue - start with neighbors of s as s requires special flow handling
        queue<FlowNode> q;
        for (Neighbor n: node_data[s].neighbors)
            if (contains(n.node) && node_data[n.node].inflow != s) {
                assert(node_data[n.node].inflow == NO_NODE);
                node_data[n.node].distance = 1;
                node_data[n.node].outcopy_distance = 1; // treat inner-node edges as length 0
                q.push(FlowNode(n.node, false));
            }
        // BFS
        while (!q.empty()) {
            FlowNode fn = q.front();
            q.pop();

            distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
            NodeID inflow = node_data[fn.node].inflow;
            // special treatment is needed for node with flow through it
            if (inflow != NO_NODE && !fn.outcopy) {
                // inflow is only valid neighbor
                if (update_distance(node_data[inflow].outcopy_distance, fn_dist + 1)) {
                    // need to set distance for 0-distance nodes immediately
                    // otherwise a longer path may set wrong distance value first
                    update_distance(node_data[inflow].distance, fn_dist + 1);
                    q.push(FlowNode(inflow, true));
                }
            } else {
                // when arriving at the outgoing copy of flow node, all neighbors except outflow are valid
                // outflow must have been already visited in this case, so checking all neighbors is fine
                for (Neighbor n: node_data[fn.node].neighbors) {
                    // filter neighbors nodes not belonging to subgraph
                    if (!contains(n.node))
                        continue;
                    // following inflow by inverting flow requires special handling
                    if (n.node == inflow) {
                        if (update_distance(node_data[n.node].outcopy_distance, fn_dist + 1)) {
                            // neighbor must be a flow node
                            update_distance(node_data[n.node].distance, fn_dist + 1);
                            q.push(FlowNode(n.node, true));
                        }
                    } else {
                        if (update_distance(node_data[n.node].distance, fn_dist + 1)) {
                            // neighbor may be a flow node
                            if (node_data[n.node].inflow == NO_NODE)
                                update_distance(node_data[n.node].outcopy_distance, fn_dist + 1);
                            q.push(FlowNode(n.node, false));
                        }
                    }
                }
            }
        }
    }

    void Graph::run_flow_bfs_from_t() {
        CHECK_CONSISTENT;
        assert(contains(s) && contains(t));
        // init distances
        for (NodeID node: nodes)
            node_data[node].distance = node_data[node].outcopy_distance = infinity;
        node_data[t].distance = node_data[t].outcopy_distance = 0;
        // init queue - start with neighbors of t as t requires special flow handling
        queue<FlowNode> q;
        for (Neighbor n: node_data[t].neighbors)
            if (contains(n.node) && node_data[n.node].outflow != t) {
                assert(node_data[n.node].outflow == NO_NODE);
                node_data[n.node].outcopy_distance = 1;
                node_data[n.node].distance = 1; // treat inner-node edges as length 0
                q.push(FlowNode(n.node, true));
            }
        // BFS
        while (!q.empty()) {
            FlowNode fn = q.front();
            q.pop();

            distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
            NodeID outflow = node_data[fn.node].outflow;
            // special treatment is needed for node with flow through it
            if (outflow != NO_NODE && fn.outcopy) {
                // outflow is only valid neighbor
                if (update_distance(node_data[outflow].distance, fn_dist + 1)) {
                    // need to set distance for 0-distance nodes immediately
                    // otherwise a longer path may set wrong distance value first
                    update_distance(node_data[outflow].outcopy_distance, fn_dist + 1);
                    q.push(FlowNode(outflow, false));
                }
            } else {
                // when arriving at the incoming copy of flow node, all neighbors except inflow are valid
                // inflow must have been already visited in this case, so checking all neighbors is fine
                for (Neighbor n: node_data[fn.node].neighbors) {
                    // filter neighbors nodes not belonging to subgraph
                    if (!contains(n.node))
                        continue;
                    // following outflow by inverting flow requires special handling
                    if (n.node == outflow) {
                        if (update_distance(node_data[n.node].distance, fn_dist + 1)) {
                            // neighbor must be a flow node
                            update_distance(node_data[n.node].outcopy_distance, fn_dist + 1);
                            q.push(FlowNode(n.node, false));
                        }
                    } else {
                        if (update_distance(node_data[n.node].outcopy_distance, fn_dist + 1)) {
                            // neighbor may be a flow node
                            if (node_data[n.node].outflow == NO_NODE)
                                update_distance(node_data[n.node].distance, fn_dist + 1);
                            q.push(FlowNode(n.node, true));
                        }
                    }
                }
            }
        }
    }

    distance_t Graph::get_distance(NodeID v, NodeID w, bool weighted) {
        assert(contains(v) && contains(w));
        weighted ? run_dijkstra(v) : run_bfs(v);
        return node_data[w].distance;
    }

    pair <NodeID, distance_t> Graph::get_furthest(NodeID v, bool weighted) {
        NodeID furthest = v;

        weighted ? run_dijkstra(v) : run_bfs(v);
        for (NodeID node: nodes)
            if (node_data[node].distance > node_data[furthest].distance)
                furthest = node;
        return make_pair(furthest, node_data[furthest].distance);
    }

    Edge Graph::get_furthest_pair(bool weighted) {
        assert(nodes.size() > 1);
        distance_t max_dist = 0;
        NodeID start = nodes[0];
        pair<NodeID, distance_t> furthest = get_furthest(start, weighted);
        while (furthest.second > max_dist) {
            max_dist = furthest.second;
            start = furthest.first;
            furthest = get_furthest(start, weighted);
        }
        return Edge(start, furthest.first, max_dist);
    }

    distance_t Graph::diameter(bool weighted) {
        if (nodes.size() < 2)
            return 0;
        return get_furthest_pair(weighted).d;
    }

    void Graph::get_diff_data(std::vector<DiffData> &diff, NodeID a, NodeID b, bool weighted, bool pre_computed) {
        CHECK_CONSISTENT;
        assert(diff.empty());
        assert(!pre_computed || node_data[a].distance == 0);
        diff.reserve(nodes.size());
        // init with distances to a
        if (!pre_computed)
            weighted ? run_dijkstra(a) : run_bfs(a);
        for (NodeID node: nodes)
            diff.push_back(DiffData(node, node_data[node].distance, 0));
        // add distances to b
        weighted ? run_dijkstra(b) : run_bfs(b);
        for (DiffData &dd: diff)
            dd.dist_b = node_data[dd.node].distance;
    }

// helper function for sorting connected components by size
    static bool cmp_size_desc(const vector <NodeID> &a, const vector <NodeID> &b) {
        return a.size() > b.size();
    };

// helper function for adding nodes to smaller of two sets
    static void add_to_smaller(vector <NodeID> &pa, vector <NodeID> &pb, const vector <NodeID> &cc) {
        vector<NodeID> &smaller = pa.size() <= pb.size() ? pa : pb;
        smaller.insert(smaller.begin(), cc.cbegin(), cc.cend());
    }

    bool Graph::get_rough_partition(Partition &p, double balance, bool disconnected) {
        DEBUG("get_rough_partition, p=" << p << ", disconnected=" << disconnected << " on " << *this);
        CHECK_CONSISTENT;
        assert(p.left.empty() && p.cut.empty() && p.right.empty());
        if (disconnected) {
            vector<vector<NodeID>> cc;
            get_connected_components(cc);
            if (cc.size() > 1) {
                DEBUG("found multiple connected components: " << cc);
                sort(cc.begin(), cc.end(), cmp_size_desc);
                // for size zero cuts we loosen the balance requirement
                if (cc[0].size() < nodes.size() * (1 - balance / 2)) {
                    for (vector<NodeID> &c: cc)
                        add_to_smaller(p.left, p.right, c);
                    return true;
                }
                // get rough partion over main component
                Graph main_cc(cc[0].begin(), cc[0].end());
                bool is_fine = main_cc.get_rough_partition(p, balance, false);
                // reset subgraph ids
                for (NodeID node: main_cc.nodes)
                    node_data[node].subgraph_id = subgraph_id;
                if (is_fine) {
                    // distribute remaining components
                    for (size_t i = 1; i < cc.size(); i++)
                        add_to_smaller(p.left, p.right, cc[i]);
                }
                return is_fine;
            }
        }
        // graph is connected - find two extreme points
#ifdef NDEBUG
        NodeID a = get_furthest(random_node(), weighted_furthest).first;
#else
        NodeID a = get_furthest(nodes[0], weighted_furthest).first;
#endif
        NodeID b = get_furthest(a, weighted_furthest).first;
        DEBUG("furthest nodes: a=" << a << ", b=" << b);
        // get distances from a and b and sort by difference
        vector<DiffData> diff;
        get_diff_data(diff, a, b, weighted_diff, weighted_furthest);
        sort(diff.begin(), diff.end(), DiffData::cmp_diff);
        DEBUG("diff=" << diff);
        // get parition bounds based on balance; round up if possible
        size_t max_left = min(nodes.size() / 2, static_cast<size_t>(ceil(nodes.size() * balance)));
        size_t min_right = nodes.size() - max_left;
        DEBUG("max_left=" << max_left << ", min_right=" << min_right);
        assert(max_left <= min_right);
        // check for corner case where most nodes have same distance difference
        if (diff[max_left - 1].diff() == diff[min_right].diff()) {
            // find bottleneck(s)
            const int32_t center_diff_value = diff[min_right].diff();
            distance_t min_dist = infinity;
            vector<NodeID> bottlenecks;
            for (DiffData dd: diff)
                if (dd.diff() == center_diff_value) {
                    if (dd.min() < min_dist) {
                        min_dist = dd.min();
                        bottlenecks.clear();
                    }
                    if (dd.min() == min_dist)
                        bottlenecks.push_back(dd.node);
                }
            sort(bottlenecks.begin(), bottlenecks.end());
            DEBUG("bottlenecks=" << bottlenecks);
            // try again with bottlenecks removed
            remove_nodes(bottlenecks);
            bool is_fine = get_rough_partition(p, balance, true);
            // add bottlenecks back to graph and to center partition
            for (NodeID bn: bottlenecks) {
                add_node(bn);
                p.cut.push_back(bn);
            }
            // if bottlenecks are the only cut vertices, they must form a minimal cut
            return is_fine && p.cut.size() == bottlenecks.size();
        }
        // ensure left and right pre-partitions are connected
        while (diff[max_left - 1].diff() == diff[max_left].diff())
            max_left++;
        while (diff[min_right - 1].diff() == diff[min_right].diff())
            min_right--;
        // assign nodes to left/cut/right
        for (size_t i = 0; i < diff.size(); i++) {
            if (i < max_left)
                p.left.push_back(diff[i].node);
            else if (i < min_right)
                p.cut.push_back(diff[i].node);
            else
                p.right.push_back(diff[i].node);
        }
        return false;
    }

    void Graph::min_vertex_cuts(vector <vector<NodeID>> &cuts) {
        DEBUG("min_vertex_cut over " << *this);
        CHECK_CONSISTENT;
        assert(contains(s) && contains(t));
        // set flow to empty
        for (NodeID node: nodes)
            node_data[node].inflow = node_data[node].outflow = NO_NODE;
#ifndef NDEBUG
        size_t last_s_distance = 1; // min s_distance is 2
#endif
        // find max s-t flow using Dinitz' algorithm
        while (true) {
            // construct BFS tree from t
            run_flow_bfs_from_t();
            DEBUG("BFS-tree: " << distances());
            const distance_t s_distance = node_data[s].outcopy_distance;
            if (s_distance == infinity)
                break;
            assert(s_distance > last_s_distance && (last_s_distance = s_distance));
            // run DFS from s along inverse BFS tree edges
            vector<NodeID> path;
            vector<FlowNode> stack;
            // iterating over neighbors of s directly simplifies stack cleanup after new s-t path is found
            for (Neighbor sn: node_data[s].neighbors) {
                if (!contains(sn.node) || node_data[sn.node].distance != s_distance - 1)
                    continue;
                // ensure edge from s to neighbor exists in residual graph
                if (node_data[sn.node].inflow != NO_NODE) {
                    assert(node_data[sn.node].inflow == s);
                    continue;
                }
                stack.push_back(FlowNode(sn.node, false));
                while (!stack.empty()) {
                    FlowNode fn = stack.back();
                    stack.pop_back();
                    DEBUG("fn=" << fn);
                    // clean up path (back tracking)
                    distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
                    // safeguard against re-visiting node during DFS (may have been enqueued before first visit)
                    if (fn_dist == infinity)
                        continue;
                    assert(fn_dist < s_distance && s_distance - fn_dist - 1 <= path.size());
                    path.resize(s_distance - fn_dist - 1);
                    // increase flow when s-t path is found
                    if (fn.node == t) {
                        DEBUG("flow path=" << path);
                        assert(node_data[path.front()].inflow == NO_NODE);
                        node_data[path.front()].inflow = s;
                        for (size_t path_pos = 1; path_pos < path.size(); path_pos++) {
                            NodeID from = path[path_pos - 1];
                            NodeID to = path[path_pos];
                            // we might be reverting existing flow
                            // from.inflow may have been changed already => check outflow
                            if (node_data[to].outflow == from) {
                                node_data[to].outflow = NO_NODE;
                                if (node_data[from].inflow == to)
                                    node_data[from].inflow = NO_NODE;
                            } else {
                                node_data[from].outflow = to;
                                node_data[to].inflow = from;
                            }
                        }
                        assert(node_data[path.back()].outflow == NO_NODE);
                        node_data[path.back()].outflow = t;
                        // skip to next neighbor of s
                        stack.clear();
                        path.clear();
                        DEBUG("new flow=" << flow());
                        break;
                    }
                    // ensure vertex is not re-visited during current DFS iteration
                    if (fn.outcopy)
                        node_data[fn.node].outcopy_distance = infinity;
                    else
                        node_data[fn.node].distance = infinity;
                    // continue DFS from node
                    path.push_back(fn.node);
                    distance_t next_distance = fn_dist - 1;
                    // when arriving at outgoing copy of a node with flow through it,
                    // we are inverting outflow, so all neighbors are valid (except outflow)
                    // otherwise inverting the inflow is the only possible option
                    NodeID inflow = node_data[fn.node].inflow;
                    if (inflow != NO_NODE && !fn.outcopy) {
                        if (node_data[inflow].outcopy_distance == next_distance)
                            stack.push_back(FlowNode(inflow, true));
                    } else {
                        for (Neighbor n: node_data[fn.node].neighbors) {
                            if (!contains(n.node))
                                continue;
                            // inflow inversion requires special handling
                            if (n.node == inflow) {
                                if (node_data[inflow].outcopy_distance == next_distance)
                                    stack.push_back(FlowNode(inflow, true));
                            } else {
                                if (node_data[n.node].distance == next_distance)
                                    stack.push_back(FlowNode(n.node, false));
                            }
                        }
                    }
                }
            }
        }
        // find min cut
        assert(cuts.empty());
        cuts.resize(1);
        // node-internal edge appears in cut iff outgoing copy is reachable from t in inverse residual graph and incoming copy is not
        // for node-external edges reachability of endpoint but unreachability of starting point is only possible if endpoint is t
        // in that case, starting point must become the cut vertex
        for (NodeID node: nodes) {
            NodeID outflow = node_data[node].outflow;
            // distance already stores distance from t in inverse residual graph
            if (outflow != NO_NODE) {
                assert(node_data[node].inflow != NO_NODE);
                if (node_data[node].outcopy_distance < infinity) {
                    // check inner edge
                    if (node_data[node].distance == infinity)
                        cuts[0].push_back(node);
                } else {
                    // check outer edge
                    if (outflow == t)
                        cuts[0].push_back(node);
                }
            }
        }
#ifdef MULTI_CUT
        // same thing but w.r.t. reachability from s in residual graph
        run_flow_bfs_from_s();
        cuts.resize(2);
        // distance now stores distance from s in residual graph
        for (NodeID node: nodes) {
            NodeID inflow = node_data[node].inflow;
            if (inflow != NO_NODE) {
                assert(node_data[node].outflow != NO_NODE);
                if (node_data[node].distance < infinity) {
                    // check inner edge
                    if (node_data[node].outcopy_distance == infinity)
                        cuts[1].push_back(node);
                } else {
                    // check outer edge
                    if (inflow == s)
                        cuts[1].push_back(node);
                }
            }
        }
        // eliminate potential duplicate
        if (cuts[0] == cuts[1])
            cuts.resize(1);
#endif
        DEBUG("cuts=" << cuts);
    }

    void Graph::get_connected_components(vector <vector<NodeID>> &components) {
        CHECK_CONSISTENT;
        components.clear();
        for (NodeID start_node: nodes) {
            // visited nodes are temporarily removed
            if (!contains(start_node))
                continue;
            node_data[start_node].subgraph_id = NO_SUBGRAPH;
            // create new connected component
            components.push_back(vector<NodeID>());
            vector<NodeID> &cc = components.back();
            vector<NodeID> stack;
            stack.push_back(start_node);
            while (!stack.empty()) {
                NodeID node = stack.back();
                stack.pop_back();
                cc.push_back(node);
                for (Neighbor n: node_data[node].neighbors)
                    if (contains(n.node)) {
                        node_data[n.node].subgraph_id = NO_SUBGRAPH;
                        stack.push_back(n.node);
                    }
            }
        }
        // reset subgraph IDs
        assign_nodes();
        DEBUG("components=" << components);
        assert(util::size_sum(components) == nodes.size());
    }

    void Graph::rough_partition_to_cuts(vector <vector<NodeID>> &cuts, const Partition &p) {
        // build subgraphs for rough partitions
        Graph left(p.left.cbegin(), p.left.cend());
        Graph center(p.cut.cbegin(), p.cut.cend());
        Graph right(p.right.cbegin(), p.right.cend());
        // construct s-t flow graph
        center.add_node(s);
        center.add_node(t);
        // handle corner case of edges between left and right partition
        // do this first as it can eliminate other s/t neighbors
        vector<NodeID> s_neighbors, t_neighbors;
        for (NodeID node: left.nodes)
            for (Neighbor n: node_data[node].neighbors)
                if (right.contains(n.node)) {
                    s_neighbors.push_back(node);
                    t_neighbors.push_back(n.node);
                }
        util::make_set(s_neighbors);
        util::make_set(t_neighbors);
        // update pre-partition
        DEBUG("moving " << s_neighbors << " and " << t_neighbors << " to center");
        left.remove_nodes(s_neighbors);
        for (NodeID node: s_neighbors)
            center.add_node(node);
        right.remove_nodes(t_neighbors);
        for (NodeID node: t_neighbors)
            center.add_node(node);
        DEBUG("pre-partition=" << left.nodes << "|" << center.nodes << "|" << right.nodes);
        // identify additional neighbors of s and t
        for (NodeID node: left.nodes)
            for (Neighbor n: node_data[node].neighbors)
                if (center.contains(n.node))
                    s_neighbors.push_back(n.node);
        for (NodeID node: right.nodes)
            for (Neighbor n: node_data[node].neighbors)
                if (center.contains(n.node))
                    t_neighbors.push_back(n.node);
        util::make_set(s_neighbors);
        util::make_set(t_neighbors);
        // add edges incident to s and t
        for (NodeID node: s_neighbors)
            center.add_edge(s, node, 1, true);
        for (NodeID node: t_neighbors)
            center.add_edge(t, node, 1, true);
        // find minimum cut
        center.min_vertex_cuts(cuts);
        // revert s-t addition
        for (NodeID node: t_neighbors) {
            assert(node_data[node].neighbors.back().node == t);
            node_data[node].neighbors.pop_back();
        }
        node_data[t].neighbors.clear();
        for (NodeID node: s_neighbors) {
            assert(node_data[node].neighbors.back().node == s);
            node_data[node].neighbors.pop_back();
        }
        node_data[s].neighbors.clear();
        // repair subgraph IDs
        assign_nodes();
    }

    void Graph::complete_partition(Partition &p) {
        CHECK_CONSISTENT;
        util::make_set(p.cut);
        remove_nodes(p.cut);
        // create left/right partitions
        p.left.clear();
        p.right.clear();
        vector<vector<NodeID>> components;
        get_connected_components(components);
        sort(components.begin(), components.end(), cmp_size_desc);
        for (const vector<NodeID> &cc: components)
            add_to_smaller(p.left, p.right, cc);
        // add cut vertices back to graph
        for (NodeID node: p.cut)
            add_node(node);
        assert(p.left.size() + p.right.size() + p.cut.size() == nodes.size());
    }

    void Graph::create_partition(Partition &p, double balance) {
        CHECK_CONSISTENT;
        assert(nodes.size() > 1);
        DEBUG("create_partition, p=" << p << " on " << *this);
        // find initial rough partition
#ifdef NO_SHORTCUTS
        bool is_fine = get_rough_partition(p, balance, true);
#else
        bool is_fine = get_rough_partition(p, balance, false);
#endif
        if (is_fine) {
            DEBUG("get_rough_partition found partition=" << p);
            return;
        }
        // find minimum cut
        vector<vector<NodeID>> cuts;
        rough_partition_to_cuts(cuts, p);
        assert(cuts.size() > 0);
        // create partition
        p.cut = cuts[0];
        complete_partition(p);
        for (size_t i = 1; i < cuts.size(); i++) {
            Partition p_alt;
            p_alt.cut = cuts[i];
            complete_partition(p_alt);
            if (p.rating() < p_alt.rating())
                p = p_alt;
        }
        DEBUG("partition=" << p);
    }

    void Graph::extend_on_partition(vector <CutIndex> &ci, double balance, uint8_t cut_level, const vector <NodeID> &p,
                                    [[maybe_unused]] const vector <NodeID> &cut, size_t leaf_size) {
        DEBUG("extend_on_partition, p=" << p << ", cut=" << cut);
        if (p.size() > 1) {
            Graph g(p.begin(), p.end());
#ifndef NO_SHORTCUTS
            START_TIMER;
            g.add_shortcuts(cut, ci);
            STOP_TIMER(t_shortcut);
#endif
            g.extend_cut_index(ci, balance, cut_level + 1, leaf_size);
        } else if (p.size() == 1) {
            ci[p[0]].cut_level = cut_level + 1;
            ci[p[0]].dist_index.push_back(ci[p[0]].dist_index[cut_level] + 1);
            assert(ci[p[0]].is_consistent());
        }
    }

    void Graph::extend_cut_index(vector <CutIndex> &ci, double balance, uint8_t cut_level, size_t leaf_size) {
        DEBUG("extend_cut_index at level " << (int) cut_level << " on " << *this);
        DEBUG("cut index=" << ci);
        CHECK_CONSISTENT;
        assert(cut_level <= MAX_CUT_LEVEL);
        if (node_count() < 2) {
            assert(cut_level == 0);
            for (NodeID node: nodes) {
                ci[node].cut_level = 0;
                ci[node].dist_index.push_back(0);
            }
            return;
        }
        // find balanced cut
        Partition p;
        if (cut_level < MAX_CUT_LEVEL) {
            START_TIMER;
            create_partition(p, balance);
#ifdef CUT_REPEAT
            for (size_t i = 1; i < CUT_REPEAT; i++)
            {
                Partition p_new;
                create_partition(p_new, balance);
                if (p_new.rating() > p.rating())
                    p = p_new;
            }
#endif
            STOP_TIMER(t_partition);
        } else
            p.cut = nodes;

        for (size_t c = 0; c < p.cut.size(); c++)
            node_data[p.cut[c]].landmark_level = p.cut.size() - c;

        // update dist_index
        for (NodeID node: nodes) {
            assert(ci[node].dist_index.size() == cut_level);
            const auto prev_lvl_val =
                    (cut_level == 0 || ci[node].dist_index.size() < cut_level) ? 0 : ci[node].dist_index[cut_level - 1];
            if (node_data[node].landmark_level == 0)
                ci[node].dist_index.push_back(prev_lvl_val + p.cut.size());
            else
                ci[node].dist_index.push_back(prev_lvl_val + (p.cut.size() - node_data[node].landmark_level + 1));
        }

        // track truncated labels
        if (!ci[nodes[0]].truncated && nodes.size() <= leaf_size) {
            for (NodeID node: nodes) {
                ci[node].truncated = true;
            }
        }

        // set cut_level
        for (NodeID c: p.cut) {
            ci[c].cut_level = cut_level;
            assert(ci[c].is_consistent());
        }

        // update partition bitstring
        for (NodeID node: p.right)
            ci[node].partition |= (static_cast<uint64_t>(1) << cut_level);
        DEBUG("cut index extended to " << ci);

        // reset landmark flags
        for (NodeID c: p.cut)
            node_data[c].landmark_level = 0;
        STOP_TIMER(t_label);

        // add shortcuts and recurse
#ifdef MULTI_THREAD
        if (nodes.size() > thread_threshold)
        {
            std::thread t_left(extend_on_partition, std::ref(ci), balance, cut_level, std::cref(p.left), std::cref(p.cut), leaf_size);
            extend_on_partition(ci, balance, cut_level, p.right, p.cut, leaf_size);
            t_left.join();
        }
        else
#endif
        {
            extend_on_partition(ci, balance, cut_level, p.left, p.cut, leaf_size);
            extend_on_partition(ci, balance, cut_level, p.right, p.cut, leaf_size);
        }
    }

    void Graph::create_cut_index(std::vector<CutIndex> &ci, double balance, size_t leaf_size) {
#ifndef NPROFILE
        t_partition = t_label = t_shortcut = 0;
#endif
        assert(is_undirected());
#ifndef NDEBUG
        // sort neighbors to make algorithms deterministic
        for (NodeID node: nodes)
            sort(node_data[node].neighbors.begin(), node_data[node].neighbors.end());
#endif
        // store original neighbor counts
        vector<NodeID> original_nodes = nodes;
        vector<size_t> original_neighbors(node_data.size());
        for (NodeID node: nodes)
            original_neighbors[node] = node_data[node].neighbors.size();
        // create index
        ci.clear();
        ci.resize(node_data.size() - 2);
        // reduce memory fragmentation by pre-allocating sensible values
        for (NodeID node: nodes) {
            ci[node].dist_index.reserve(32);
        }

        extend_cut_index(ci, balance, 0, leaf_size);
        log_progress(0);
        // reset nodes (top-level cut vertices got removed)
        nodes = original_nodes;

/*#ifndef NDEBUG
    for (NodeID node : nodes)
        if (!ci[node].is_consistent())
            cerr << "inconsistent cut index for node " << node << ": "<< ci[node] << endl;
#endif*/
#ifndef NPROFILE
        cerr << "partitioning took " << t_partition << "s" << endl;
#endif
    }

    void Graph::get_redundant_edges(std::vector<Edge> &edges) {
        CHECK_CONSISTENT;
        assert(edges.empty());
        // reset distances for all nodes
        for (NodeID node: nodes)
            node_data[node].distance = infinity;
        // run localized Dijkstra from each node
        vector<NodeID> visited;
        priority_queue<SearchNode> q;
        for (NodeID v: nodes) {
            node_data[v].distance = 0;
            visited.push_back(v);
            distance_t max_dist = 0;
            // init queue - starting from neighbors ensures that only paths of length 2+ are considered
            for (Neighbor n: node_data[v].neighbors)
                if (contains(n.node)) {
                    q.push(SearchNode(n.distance, n.node));
                    if (v < n.node)
                        max_dist = max(max_dist, n.distance);
                }
            // dijkstra
            while (!q.empty()) {
                SearchNode next = q.top();
                q.pop();

                for (Neighbor n: node_data[next.node].neighbors) {
                    // filter neighbors nodes not belonging to subgraph
                    if (!contains(n.node))
                        continue;
                    // update distance and enque
                    distance_t new_dist = next.distance + n.distance;
                    if (new_dist <= max_dist && new_dist < node_data[n.node].distance) {
                        node_data[n.node].distance = new_dist;
                        q.push(SearchNode(new_dist, n.node));
                        visited.push_back(n.node);
                    }
                }
            }
            // identify redundant edges
            for (Neighbor n: node_data[v].neighbors)
                // only add redundant edges once
                if (v < n.node && contains(n.node) && node_data[n.node].distance <= n.distance)
                    edges.push_back(Edge(v, n.node, n.distance));
            // cleanup
            for (NodeID w: visited)
                node_data[w].distance = infinity;
            visited.clear();
        }
    }

    void Graph::contract(vector <Neighbor> &closest, bool flag) {
        closest.resize(node_data.size() - 2, Neighbor(NO_NODE, 0));
        for (NodeID node: nodes)
            closest[node] = Neighbor(node, 0);

        if (flag) {
            // helper function to identify degree one nodes and associated neighbors
            auto find_degree_one = [this, &closest](const vector<NodeID> &nodes, vector<NodeID> &degree_one,
                                                    vector<NodeID> &neighbors) {
                degree_one.clear();
                neighbors.clear();
                for (NodeID node: nodes) {
                    Neighbor neighbor = single_neighbor(node);
                    if (neighbor.node != NO_NODE) {
                        // avoid complete contraction (screws with testing)
                        if (single_neighbor(neighbor.node).node == NO_NODE) {
                            closest[node] = neighbor;
                            degree_one.push_back(node);
                            neighbors.push_back(neighbor.node);
                        }
                    }
                }
            };
            // remove nodes
            vector<NodeID> degree_one, neighbors;
            find_degree_one(nodes, degree_one, neighbors);
            while (!degree_one.empty()) {
                sort(degree_one.begin(), degree_one.end());
                remove_nodes(degree_one);
                vector<NodeID> old_neighbors = neighbors;
                find_degree_one(old_neighbors, degree_one, neighbors);
            }
        }
    }

//--------------------------- ContractionHierarchy ------------------

    ContractionHierarchy::ContractionHierarchy() {

    }

/*ContractionHierarchy::ContractionHierarchy(istream &is) {

    size_t count, node_count, i, j;
    is.read((char*)&node_count, sizeof(size_t));
    nodes.resize(node_count);
    for(i = 1; i < node_count; i++) {
        is.read((char*)&nodes[i].dist_index, sizeof(uint16_t));

        is.read((char*)&count, sizeof(size_t));
        nodes[i].up_neighbors.reserve(count);
        for(j = 0; j < count; j++) {
            Neighbor node(NO_NODE, 0);
            is.read((char*)&node.node, sizeof(NodeID));
            is.read((char*)&node.distance, sizeof(distance_t));
            nodes[i].up_neighbors.push_back(node);
        }

        is.read((char*)&count, sizeof(size_t));
        nodes[i].down_neighbors.reserve(count);
        for(j = 0; j < count; j++) {
            NodeID node;
            is.read((char*)&node, sizeof(NodeID));
            nodes[i].down_neighbors.push_back(node);
        }
    }

    is.read((char*)&count, sizeof(size_t));
    bottom_up_nodes.reserve(count);
    for(i = 0; i < count; i++) {
        NodeID node;
        is.read((char*)&node, sizeof(NodeID));
        bottom_up_nodes.push_back(node);
    }
}

void ContractionHierarchy::write(ostream &os) {

    size_t count = nodes.size(), i;
    os.write((char*)&count, sizeof(size_t));
    for(i = 1; i < nodes.size() ; i++) {
        os.write((char*)&nodes[i].dist_index, sizeof(uint16_t));

        count = nodes[i].up_neighbors.size();
        os.write((char*)&count, sizeof(size_t));
        for(Neighbor n: nodes[i].up_neighbors) {
            os.write((char*)&n.node, sizeof(NodeID));
            os.write((char*)&n.distance, sizeof(distance_t));
        }

        count = nodes[i].down_neighbors.size();
        os.write((char*)&count, sizeof(size_t));
        for(NodeID n: nodes[i].down_neighbors)
            os.write((char*)&n, sizeof(NodeID));
    }

    // storing ordering information
    count = bottom_up_nodes.size();
    os.write((char*)&count, sizeof(size_t));
    for(NodeID node: bottom_up_nodes)
        os.write((char*)&node, sizeof(NodeID));
}*/

    size_t ContractionHierarchy::size() const {
        size_t total = 0;
        for (NodeID i = 1; i < nodes.size(); i++) {
            total += sizeof(uint64_t);
            total += nodes[i].up_neighbors.size() * sizeof(NodeID);
            total += nodes[i].up_neighbors.size() * sizeof(distance_t);
            total += nodes[i].down_neighbors.size() * sizeof(NodeID);
        }
        return total;
    }

    int Graph::UpNeighbor_position(const ContractionHierarchy &ch, NodeID v, NodeID w) const {
        for (size_t i = 0; i < ch.nodes[v].up_neighbors.size(); i++) {
            if (ch.nodes[v].up_neighbors[i].neighbor == w) return i;
        }
        assert(false);
        return -1;
    }

    void Graph::initialize(ContractionHierarchy &ch, std::vector<CutIndex> &ci, vector <Neighbor> &closest) const {
        ch.bottom_up_nodes.reserve(nodes.size() + 1);
        ch.contracted_nodes.reserve(1000);

        // setting up path_data
        const path_data blank_path;

        // initialize distance index to determine edge direction
        ch.nodes.resize(nodes.size() + 1);
        for (NodeID node: nodes) {
            if (closest[node].node == node) {
                ch.bottom_up_nodes.push_back(node);
                ch.nodes[node].dist_index = ci[node].dist_index[ci[node].cut_level] - 1;
                if (!ci[node].truncated) {
                    ci[node].distances.resize(ch.nodes[node].dist_index + 1, infinity);
                    ci[node].paths.resize(ch.nodes[node].dist_index + 1, nullptr);
                }
            } else {
                ch.contracted_nodes.push_back(node);
                ch.nodes[node].up_neighbors.push_back(CHNeighbor(closest[node].node, infinity, blank_path));
            }
        }

        // initialize with upwards graph edges
        for (NodeID node: ch.bottom_up_nodes) {
            for (Neighbor &n: node_data[node].neighbors)
                if (closest[n.node].node == n.node && ch.nodes[n.node].dist_index < ch.nodes[node].dist_index)
                    ch.nodes[node].up_neighbors.push_back(CHNeighbor(n.node, infinity, blank_path));
        }

        // add shortcuts bottom-up
        auto di_order = [&ch](NodeID a, NodeID b) -> bool {
            return ch.nodes[a].dist_index > ch.nodes[b].dist_index;
        };
        auto di_order1 = [&ch](CHNeighbor a, CHNeighbor b) -> bool {
            return ch.nodes[a.neighbor].dist_index > ch.nodes[b.neighbor].dist_index;
        };

        std::sort(ch.bottom_up_nodes.begin(), ch.bottom_up_nodes.end(), di_order);
        for (NodeID node: ch.bottom_up_nodes) {
            vector<CHNeighbor> &up = ch.nodes[node].up_neighbors;
            util::make_set(up, di_order1);

            for (size_t i = 0; i + 1 < up.size(); i++) {
                for (size_t j = i + 1; j < up.size(); j++)
                    ch.nodes[up[i].neighbor].up_neighbors.push_back(CHNeighbor(up[j].neighbor, infinity, blank_path));
            }

            // create downward neighbors from upward ones
            for (CHNeighbor upn: up)
                ch.nodes[upn.neighbor].down_neighbors.push_back(node);
        }
    }

    void Graph::customise_shortcut_graph(ContractionHierarchy &ch, ContractionIndex &tcl, vector <Edge> &edges) const {

        // update graph weights
        for (Edge &e: edges) {
            NodeID a = e.a, b = e.b;

            ContractionLabel ca = tcl.get_contraction_label(a);
            if (ca.root != NO_NODE && ch.nodes[a].up_neighbors[0].neighbor == b) {
                ch.nodes[a].up_neighbors[0].distance = e.d;
                continue;
            }
            ContractionLabel cb = tcl.get_contraction_label(b);
            if (cb.root != NO_NODE && ch.nodes[b].up_neighbors[0].neighbor == a) {
                ch.nodes[b].up_neighbors[0].distance = e.d;
                continue;
            }

            if (ch.nodes[a].dist_index < ch.nodes[b].dist_index) swap(a, b);
            for (CHNeighbor &n: ch.nodes[a].up_neighbors) {
                if (n.neighbor == b) {
                    n.distance = e.d;
                    ca = tcl.get_contraction_label(a);
                    if (ca.cut_index.label_count() != 0) {
                        ca.cut_index.distances()[ch.nodes[b].dist_index] = e.d;
                        ca.cut_index.paths()[ch.nodes[b].dist_index] = n;
                    }
                    break;
                }
            }
        }

#ifdef MULTI_THREAD_DISTANCES
        vector<thread> threads;
        auto chp = [this](ContractionHierarchy &ch, ContractionIndex &tcl, util::par_bucket_list<NodeID, MULTI_THREAD_DISTANCES> &que, size_t id) {

            NodeID x;
            while (que.next(x, id)) {
                FlatCutIndex cx = tcl.get_contraction_label(x).cut_index;

                for(NodeID z: ch.nodes[x].down_neighbors) {
                    int i = ch.nodes[x].up_neighbors.size() - 1, j = ch.nodes[z].up_neighbors.size() - 1, k = UpNeighbor_position(ch, z, x);

                    while (j > k) {
                        NodeID a = ch.nodes[x].up_neighbors[i].neighbor, b = ch.nodes[z].up_neighbors[j].neighbor;

                        if(a != b) i--;
                        else {
                            distance_t m = ch.nodes[z].up_neighbors[j].distance + ch.nodes[z].up_neighbors[k].distance;
                            if(ch.nodes[x].up_neighbors[i].distance > m) {
                                ch.nodes[x].up_neighbors[i].distance = m;

                                // memory reset
//                                memset(ch.nodes[x].up_neighbors[i].p.intermediate, NO_NODE, sizeof(ch.nodes[x].up_neighbors[i].p.intermediate));
                                ch.nodes[x].up_neighbors[i].p.intermediate.fill(NO_NODE);
                                ch.nodes[x].up_neighbors[i].p.cs.lower = nullptr;
                                ch.nodes[x].up_neighbors[i].p.cs.upper = nullptr;
                                ch.nodes[x].up_neighbors[i].p.cs.triangle_node = NO_NODE;

                                uint32_t lower_count = ch.nodes[z].up_neighbors[k].p.intermediate_count();
                                uint32_t upper_count = lower_count == infinity ? 0 : ch.nodes[z].up_neighbors[j].p.intermediate_count();
                                if (lower_count + upper_count + 1 > MAX_VALLEY_PATH_LENGTH)
                                {
                                    ch.nodes[x].up_neighbors[i].p.cs.lower = &ch.nodes[z].up_neighbors[k];
                                    ch.nodes[x].up_neighbors[i].p.cs.upper = &ch.nodes[z].up_neighbors[j];
                                    ch.nodes[x].up_neighbors[i].p.cs.triangle_node = z;
                                }
                                else
                                {
                                    // copying intermediate nodes for lower leg z-x
                                    size_t idx = 0, idx1;
                                    for(idx1 = lower_count; idx1 > 0; idx1--, idx++)
                                        ch.nodes[x].up_neighbors[i].p.intermediate[idx] = ch.nodes[z].up_neighbors[k].p.intermediate[idx1 - 1];

                                    // copying triangle node
                                    ch.nodes[x].up_neighbors[i].p.intermediate[idx] = z; idx++;

                                    // copying intermediate nodes for upper leg z-a
                                    for(idx1 = 0; idx1 < upper_count; idx1++, idx++)
                                        ch.nodes[x].up_neighbors[i].p.intermediate[idx] = ch.nodes[z].up_neighbors[j].p.intermediate[idx1];
                                }

                                if(cx.label_count() != 0)
                                {
                                    cx.paths()[ch.nodes[a].dist_index] = ch.nodes[x].up_neighbors[i];
                                    cx.distances()[ch.nodes[a].dist_index] = ch.nodes[x].up_neighbors[i].distance;
                                }
                            }
                            i--; j--;
                        }
                    }
                }
            }
        };

        util::par_bucket_list<NodeID, MULTI_THREAD_DISTANCES> que(ch.nodes[ch.bottom_up_nodes[0]].dist_index, false);
        for (NodeID node : ch.bottom_up_nodes)
            que.push(node, ch.nodes[node].dist_index);

        // add nodes to queue (pre-compute)
        for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
            threads.push_back(thread(chp, std::ref(ch), std::ref(tcl), std::ref(que), i));
        for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
            threads[i].join();
#else
        for (NodeID x: ch.bottom_up_nodes) {
            FlatCutIndex cx = tcl.get_contraction_label(x).cut_index;
            for (NodeID z: ch.nodes[x].down_neighbors) {
                int i = ch.nodes[x].up_neighbors.size() - 1, j =
                        ch.nodes[z].up_neighbors.size() - 1, k = UpNeighbor_position(ch, z, x);
                while (j > k) {
                    NodeID a = ch.nodes[x].up_neighbors[i].neighbor, b = ch.nodes[z].up_neighbors[j].neighbor;
                    if (a != b) {
                        i--;
                    } else {
                        distance_t m = ch.nodes[z].up_neighbors[j].distance + ch.nodes[z].up_neighbors[k].distance;
                        if (ch.nodes[x].up_neighbors[i].distance > m) {
                            ch.nodes[x].up_neighbors[i].distance = m;

                            // memory reset
//                            memset(ch.nodes[x].up_neighbors[i].p.intermediate, NO_NODE,
//                                   sizeof(ch.nodes[x].up_neighbors[i].p.intermediate));
                            ch.nodes[x].up_neighbors[i].p.intermediate.fill(NO_NODE);
                            ch.nodes[x].up_neighbors[i].p.cs.lower = nullptr;
                            ch.nodes[x].up_neighbors[i].p.cs.upper = nullptr;
                            ch.nodes[x].up_neighbors[i].p.cs.triangle_node = NO_NODE;

                            uint32_t lower_count = ch.nodes[z].up_neighbors[k].p.intermediate_count();
                            uint32_t upper_count =
                                    lower_count == infinity ? 0 : ch.nodes[z].up_neighbors[j].p.intermediate_count();
                            if (lower_count + upper_count + 1 > MAX_VALLEY_PATH_LENGTH) {
                                ch.nodes[x].up_neighbors[i].p.cs.lower = &ch.nodes[z].up_neighbors[k];
                                ch.nodes[x].up_neighbors[i].p.cs.upper = &ch.nodes[z].up_neighbors[j];
                                ch.nodes[x].up_neighbors[i].p.cs.triangle_node = z;
                            } else {
                                // copying intermediate nodes for lower leg z-x
                                size_t idx = 0, idx1;
                                for (idx1 = lower_count; idx1 > 0; idx1--, idx++)
                                    ch.nodes[x].up_neighbors[i].p.intermediate[idx] = ch.nodes[z].up_neighbors[k].p.intermediate[
                                            idx1 - 1];

                                // copying triangle node
                                ch.nodes[x].up_neighbors[i].p.intermediate[idx] = z;
                                idx++;

                                // copying intermediate nodes for upper leg z-a
                                for (idx1 = 0; idx1 < upper_count; idx1++, idx++)
                                    ch.nodes[x].up_neighbors[i].p.intermediate[idx] = ch.nodes[z].up_neighbors[j].p.intermediate[idx1];
                            }

                            if (cx.label_count() != 0) {
                                cx.paths()[ch.nodes[a].dist_index] = ch.nodes[x].up_neighbors[i];
                                cx.distances()[ch.nodes[a].dist_index] = ch.nodes[x].up_neighbors[i].distance;
                            }
                        }
                        i--;
                        j--;
                    }
                }
            }
        }
#endif
    }

    void Graph::customise_hub_labelling(ContractionHierarchy &ch, ContractionIndex &tcl) const {

        // setting up path_data
        const path_data blank_path;

#ifdef MULTI_THREAD_DISTANCES
        vector<thread> threads;
        auto tclp = [this](ContractionHierarchy &ch, ContractionIndex &tcl, util::par_bucket_list<NodeID, MULTI_THREAD_DISTANCES> &que, size_t id, path_data blank_path) {

            NodeID x;
            while (que.next(x, id)) {
                FlatCutIndex cx = tcl.get_contraction_label(x).cut_index;
                if(cx.label_count() != 0) {
                    for(CHNeighbor &n: ch.nodes[x].up_neighbors) {
                        FlatCutIndex cn = tcl.get_contraction_label(n.neighbor).cut_index;
                        for(size_t anc = 0; anc + 1 < cn.dist_index()[cn.cut_level()]; anc++) {
                            distance_t m = n.distance + cn.distances()[anc];
                            if(cx.distances()[anc] > m) {
                                cx.paths()[anc] = n;
                                cx.distances()[anc] = m;
                            }
                        }
                    }
                    cx.distances()[cx.dist_index()[cx.cut_level()] - 1] = 0;
                    cx.paths()[cx.dist_index()[cx.cut_level()] - 1] = CHNeighbor(x, 0, blank_path);
                }
            }
        };

        util::par_bucket_list<NodeID, MULTI_THREAD_DISTANCES> que(ch.nodes[ch.bottom_up_nodes[0]].dist_index, true);
        for (NodeID node : ch.bottom_up_nodes)
            que.push(node, ch.nodes[node].dist_index);

        for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
            threads.push_back(thread(tclp, std::ref(ch), std::ref(tcl), std::ref(que), i, blank_path));
        for (size_t i = 0; i < MULTI_THREAD_DISTANCES; i++)
            threads[i].join();
#else
        for (auto x = ch.bottom_up_nodes.rbegin(); x != ch.bottom_up_nodes.rend(); x++) {
            FlatCutIndex cx = tcl.get_contraction_label(*x).cut_index;

            if (cx.label_count() != 0) {
                for (CHNeighbor &n: ch.nodes[*x].up_neighbors) {
                    FlatCutIndex cn = tcl.get_contraction_label(n.neighbor).cut_index;
                    for (size_t anc = 0; anc < cn.dist_index()[cn.cut_level()] - 1; anc++) {
                        distance_t m = n.distance + cn.distances()[anc];
                        if (cx.distances()[anc] > m) {
                            cx.paths()[anc] = n;
                            cx.distances()[anc] = m;
                        }
                    }
                }

                cx.distances()[cx.dist_index()[cx.cut_level()] - 1] = 0;
                cx.paths()[cx.dist_index()[cx.cut_level()] - 1] = CHNeighbor(*x, 0, blank_path);
            }
        }
#endif

        // customising contracted part
        for (NodeID x: ch.contracted_nodes) {

            NodeID root = ch.nodes[x].up_neighbors[0].neighbor;
            distance_t root_dist = ch.nodes[x].up_neighbors[0].distance;
            while (tcl.get_contraction_label(root).root != NO_NODE) {
                root_dist += ch.nodes[root].up_neighbors[0].distance;
                root = ch.nodes[root].up_neighbors[0].neighbor;
            }

            tcl.update_distance_offset(x, root_dist);
        }
    }

    void Graph::reset(ContractionHierarchy &ch) const {
        static const path_data blank_path;

        // reseting ch
        for (NodeID x: ch.bottom_up_nodes) {
            for (CHNeighbor &n: ch.nodes[x].up_neighbors) {
                n.distance = infinity;
                n.p = blank_path;
            }
        }

        for (NodeID x: ch.contracted_nodes) {
            ch.nodes[x].up_neighbors[0].distance = infinity;
            ch.nodes[x].up_neighbors[0].p = blank_path;
        }
    }

    void Graph::reset(ContractionHierarchy &ch, ContractionIndex &tcl) const {
        // setting up path_data
        static const path_data blank_path;

        // reseting ch
        for (NodeID x: ch.bottom_up_nodes) {
            for (CHNeighbor &n: ch.nodes[x].up_neighbors) {
                n.distance = infinity;
                n.p = blank_path;
            }
        }

        for (NodeID x: ch.contracted_nodes) {
            ch.nodes[x].up_neighbors[0].distance = infinity;
            ch.nodes[x].up_neighbors[0].p = blank_path;
        }

        // reseting tcl
        for (NodeID x: ch.bottom_up_nodes) {
            FlatCutIndex cx = tcl.get_contraction_label(x).cut_index;
            for (size_t anc = 0; anc < cx.label_count(); anc++)
                cx.distances()[anc] = infinity;
        }

        for (NodeID x: ch.contracted_nodes)
            tcl.update_distance_offset(x, infinity);
    }

    void Graph::get_anc_dist(const ContractionHierarchy &ch, const ContractionIndex &tcl, NodeID v, uint16_t h_max,
                             vector <QNode> &result) const {
        // setting up path_data
        const path_data blank_path;
        const QNode blank_qn(CHNeighbor(NO_NODE, infinity, blank_path), 0, NO_NODE, infinity);

        uint16_t v_height = ch.nodes[v].dist_index;
        result.resize(v_height + 1, blank_qn);
        result[v_height] = QNode(CHNeighbor(v, 0, blank_path), v_height, v, 0);

        // no need to check ancestors of root, so h > 0 is ok here
        for (uint16_t h = v_height; h > 0; h--) {
            QNode current = result[h];
            if (current.node_id == NO_NODE)
                continue;
            FlatCutIndex ci = tcl.get_contraction_label(current.node_id).cut_index;
            if (ci.label_count() == 0) {
                for (const CHNeighbor &n: ch.nodes[current.node_id].up_neighbors) {
                    QNode &anc_result = result[ch.nodes[n.neighbor].dist_index];
                    distance_t new_dist = current.distance + n.distance;
                    if (anc_result.distance > new_dist) {
                        anc_result.shortcut_used = n; // for back-tracking how we found this distance value
                        anc_result.shortcut_start = h;
                        anc_result.node_id = n.neighbor;
                        anc_result.distance = new_dist;
                    }
                }
            } else {
                uint16_t anc_max = min(h, uint16_t(h_max + 1));
                for (uint16_t anc = 0; anc < anc_max; anc++) {
                    distance_t new_dist = current.distance + ci.distances()[anc];
                    if (new_dist <= result[anc].distance) {
                        result[anc].shortcut_used = ci.paths()[anc]; // first shortcut in chain to anc
                        result[anc].shortcut_start = h; // start of first shortcut in chain to anc
                        result[anc].node_id = NO_NODE; // stop searching from here
                        result[anc].distance = new_dist;
                    }
                }
            }
        }
        DEBUG("get_anc_dist(" << v << "," << h_max << ")=" << result);
    }

// appends path from w to v (both inclusive)
    void Graph::path_from_anc(const ContractionHierarchy &ch, const ContractionIndex &tcl, NodeID v, uint16_t w,
                              std::vector<QNode> &anc_dist, std::vector<NodeID> &result) const {
        const uint16_t v_height = ch.nodes[v].dist_index;
        uint16_t current = w;

        // unpack convex path segment between last recorded CH-shortcut and w
        if (anc_dist[w].node_id == NO_NODE) // reached via label
        {
            const size_t old_size = result.size();
            unpack_convex_path(ch, tcl, anc_dist[w].shortcut_used.neighbor, w, result);
            reverse(result.begin() + old_size, result.end());
            unpack(anc_dist[w].shortcut_used, result, true);
            current = anc_dist[w].shortcut_start;
            DEBUG("path_from_anc: label-based path unpacked as " << result);
        }

        // unpack recorded CH-shortcuts
        DEBUG("path_from_anc: starting CH-unpacking from height " << current);
        while (current < v_height) {
            result.push_back(anc_dist[current].node_id);
            unpack(anc_dist[current].shortcut_used, result, true);
            current = anc_dist[current].shortcut_start;
        }
        result.push_back(v);
        DEBUG("path_from_anc(" << v << "," << w << ")=" << result);
    }

// appends path from v to w (both inclusive)
    void Graph::path_to_anc(const ContractionHierarchy &ch, const ContractionIndex &tcl, NodeID v, uint16_t w,
                            std::vector<QNode> &anc_dist, std::vector<NodeID> &result) const {
        const size_t old_size = result.size();
        path_from_anc(ch, tcl, v, w, anc_dist, result);
        reverse(result.begin() + old_size, result.end());
    }

// appends intermediate nodes on path used to create shortcut
    void Graph::unpack(CHNeighbor const &n, vector <NodeID> &path, bool reverse) const {
        if (reverse) {
            if (MAX_VALLEY_PATH_LENGTH == 0 || n.p.intermediate[0] == NO_NODE) {
                if (n.p.cs.upper != nullptr)
                    unpack(*n.p.cs.upper, path, true);
                if (n.p.cs.triangle_node != NO_NODE)
                    path.push_back(n.p.cs.triangle_node);
                if (n.p.cs.lower != nullptr)
                    unpack(*n.p.cs.lower, path, false);
            } else {
                for (size_t i = MAX_VALLEY_PATH_LENGTH; i > 0; i--) {
                    if (n.p.intermediate[i - 1] != NO_NODE)
                        path.push_back(n.p.intermediate[i - 1]);
                }
            }
        } else {
            if (MAX_VALLEY_PATH_LENGTH == 0 || n.p.intermediate[0] == NO_NODE) {
                if (n.p.cs.lower != nullptr)
                    unpack(*n.p.cs.lower, path, true);
                if (n.p.cs.triangle_node != NO_NODE)
                    path.push_back(n.p.cs.triangle_node);
                if (n.p.cs.upper != nullptr)
                    unpack(*n.p.cs.upper, path, false);
            } else {
                for (size_t i = 0; i < MAX_VALLEY_PATH_LENGTH; i++) {
                    if (n.p.intermediate[i] == NO_NODE)
                        break;
                    path.push_back(n.p.intermediate[i]);
                }
            }
        }
    }

// appends path from v to w (both inclusive)
    void Graph::unpack_convex_path(const ContractionHierarchy &ch, const ContractionIndex &tcl, NodeID v, uint16_t w,
                                   vector <NodeID> &result) const {
        DEBUG("unpack_convex_path(" << v << "," << w << "):");
        FlatCutIndex ci = tcl.get_contraction_label(v).cut_index;
        DEBUG("\tci(v)=" << ci);
        result.push_back(v);
        while (ci.ancestor_count() - 1 != w) {
            DEBUG("\tunpacking " << ci.paths()[w]);
            unpack(ci.paths()[w], result, false);
            v = ci.paths()[w].neighbor;
            ci = tcl.get_contraction_label(v).cut_index;
            result.push_back(v);
        }
    }

    vector <NodeID> Graph::query_path(const ContractionHierarchy &ch, const ContractionIndex &tcl, NodeID v, NodeID w,
                                      distance_t &dist) const {
        ContractionLabel cv = tcl.get_contraction_label(v), cw = tcl.get_contraction_label(w);
        vector<NodeID> source, target;
        NodeID v_anc, w_anc;
        assert(!cv.cut_index.empty() && !cw.cut_index.empty());

        dist = 0;
        if (cv.cut_index == cw.cut_index) {
            if (v == w) {
                source.push_back(v);
                return source;
            }
            if (cv.distance_offset == 0) {
                w_anc = w;
                while (w_anc != v) {
                    source.push_back(w_anc);
                    dist += ch.nodes[w_anc].up_neighbors[0].distance;
                    w_anc = ch.nodes[w_anc].up_neighbors[0].neighbor;
                }
                source.push_back(w_anc);
                reverse(source.begin(), source.end());
                return source;
            }
            if (cw.distance_offset == 0) {
                v_anc = v;
                while (v_anc != w) {
                    source.push_back(v_anc);
                    dist += ch.nodes[v_anc].up_neighbors[0].distance;
                    v_anc = ch.nodes[v_anc].up_neighbors[0].neighbor;
                }
                source.push_back(v_anc);
                return source;
            }
            if (ch.nodes[v].up_neighbors[0].neighbor == w) {
                source.push_back(v);
                source.push_back(w);
                dist += ch.nodes[v].up_neighbors[0].distance;
                return source;
            }
            if (ch.nodes[w].up_neighbors[0].neighbor == v) {
                source.push_back(v);
                source.push_back(w);
                dist += ch.nodes[w].up_neighbors[0].distance;
                return source;
            }

            // find lowest common ancestor
            v_anc = v;
            w_anc = w;
            source.push_back(v_anc);
            target.push_back(w_anc);
            distance_t v_dist_offset = cv.distance_offset, w_dist_offset = cw.distance_offset;
            while (v_anc != w_anc) {
                if (v_dist_offset < w_dist_offset) {
                    w_dist_offset = w_dist_offset - ch.nodes[w_anc].up_neighbors[0].distance;
                    dist += ch.nodes[w_anc].up_neighbors[0].distance;
                    w_anc = ch.nodes[w_anc].up_neighbors[0].neighbor;
                    target.push_back(w_anc);
                } else if (v_dist_offset > w_dist_offset) {
                    v_dist_offset = v_dist_offset - ch.nodes[v_anc].up_neighbors[0].distance;
                    dist += ch.nodes[v_anc].up_neighbors[0].distance;
                    v_anc = ch.nodes[v_anc].up_neighbors[0].neighbor;
                    source.push_back(v_anc);
                } else {
                    w_dist_offset = w_dist_offset - ch.nodes[w_anc].up_neighbors[0].distance;
                    v_dist_offset = v_dist_offset - ch.nodes[v_anc].up_neighbors[0].distance;
                    dist += ch.nodes[w_anc].up_neighbors[0].distance;
                    dist += ch.nodes[v_anc].up_neighbors[0].distance;
                    w_anc = ch.nodes[w_anc].up_neighbors[0].neighbor;
                    v_anc = ch.nodes[v_anc].up_neighbors[0].neighbor;
                    target.push_back(w_anc);
                    source.push_back(v_anc);
                }
            }

            for (auto it = target.rbegin() + 1; it != target.rend(); it++) source.push_back(*it);
            return source;
        }

        distance_t offset = cv.distance_offset + cw.distance_offset;
        // contracted case
        if (cv.root != NO_NODE) {
            v_anc = v;
            while (v_anc != cv.root) {
                source.push_back(v_anc);
                v_anc = ch.nodes[v_anc].up_neighbors[0].neighbor;
            }
            v = cv.root;
        }

        if (cw.root != NO_NODE) {
            w_anc = w;
            while (w_anc != cw.root) {
                target.push_back(w_anc);
                w_anc = ch.nodes[w_anc].up_neighbors[0].neighbor;
            }
            w = cw.root;
        }

        vector<NodeID> temp;
        uint16_t lca_height = tcl.common_ancestor_count(v, w), anc, hub_height = 0;
        distance_t min_dist = infinity;

        vector<QNode> v_label, w_label;
        if (cv.cut_index.label_count() != 0 && cw.cut_index.label_count() != 0) {
            for (anc = 0; anc < lca_height; anc++) {
                distance_t temp_dist = cv.cut_index.distances()[anc] + cw.cut_index.distances()[anc];
                if (min_dist > temp_dist) {
                    min_dist = temp_dist;
                    hub_height = anc;
                }
            }

            unpack_convex_path(ch, tcl, v, hub_height, source);
            unpack_convex_path(ch, tcl, w, hub_height, target);

        } else if (cv.cut_index.label_count() != 0 && cw.cut_index.label_count() == 0) {
            get_anc_dist(ch, tcl, w, lca_height, w_label);
            for (anc = 0; anc < lca_height; anc++) {
                distance_t temp_dist = cv.cut_index.distances()[anc] + w_label[anc].distance;
                if (min_dist > temp_dist) {
                    min_dist = temp_dist;
                    hub_height = anc;
                }
            }

            unpack_convex_path(ch, tcl, v, hub_height, source);
            path_to_anc(ch, tcl, w, hub_height, w_label, target);

        } else if (cv.cut_index.label_count() == 0 && cw.cut_index.label_count() != 0) {
            get_anc_dist(ch, tcl, v, lca_height, v_label);
            for (anc = 0; anc < lca_height; anc++) {
                distance_t temp_dist = v_label[anc].distance + cw.cut_index.distances()[anc];
                if (min_dist > temp_dist) {
                    min_dist = temp_dist;
                    hub_height = anc;
                }
            }

            path_to_anc(ch, tcl, v, hub_height, v_label, source);
            unpack_convex_path(ch, tcl, w, hub_height, target);

        } else {
            get_anc_dist(ch, tcl, w, lca_height, w_label);
            get_anc_dist(ch, tcl, v, lca_height, v_label);
            for (anc = 0; anc < lca_height; anc++) {
                distance_t temp_dist = v_label[anc].distance + w_label[anc].distance;
                if (min_dist > temp_dist) {
                    min_dist = temp_dist;
                    hub_height = anc;
                }
            }

            path_to_anc(ch, tcl, v, hub_height, v_label, source);
            path_to_anc(ch, tcl, w, hub_height, w_label, target);
        }

        // combing paths s-lca-t
        for (auto it = target.rbegin() + 1; it != target.rend(); it++) source.push_back(*it);

//        // Make sure path contains only valid edges
//        size_t i = 0;
//        NodeID s = source[i], t;
//        distance_t weight = 0;
//        for (i = 1; i < source.size(); i++) {
//            t = source[i];
//            bool found = false;
//            for (Neighbor n: node_data[s].neighbors) {
//                if (n.node == t) {
//                    weight += n.distance;
//                    found = true;
//                    break;
//                }
//            }
//            if (!found)
//                throw std::runtime_error(
//                        "Edge not found between " + std::to_string(s) + " and " + std::to_string(t) + ", source " +
//                        std::to_string(v) + ", target " + std::to_string(w));
//            s = t;
//        }

        dist = min_dist + offset;
        return source;
    }

//--------------------------- Graph debug ---------------------------

    bool Graph::is_consistent() const {
        // all nodes in subgraph have correct subgraph ID assigned
        for (NodeID node: nodes)
            if (node_data[node].subgraph_id != subgraph_id) {
                DEBUG("wrong subgraph ID for " << node << " in " << *this);
                return false;
            }
        // number of nodes with subgraph_id of subgraph equals number of nodes in subgraph
        size_t count = 0;
        for (NodeID node = 0; node < node_data.size(); node++)
            if (contains(node))
                count++;
        if (count != nodes.size()) {
            DEBUG(count << "/" << nodes.size() << " nodes contained in " << *this);
            return false;
        }
        return true;
    }

    bool Graph::is_undirected() const {
        for (NodeID node: nodes)
            for (Neighbor n: node_data[node].neighbors) {
                bool found = false;
                for (Neighbor nn: node_data[n.node].neighbors)
                    if (nn.node == node && nn.distance == n.distance) {
                        found = true;
                        break;
                    }
                if (!found)
                    return false;
            }
        return true;
    }

    vector <pair<distance_t, distance_t>> Graph::distances() const {
        vector<pair<distance_t, distance_t>> d;
        for (const Node &n: node_data)
            d.push_back(pair(n.distance, n.outcopy_distance));
        return d;
    }

    vector <pair<NodeID, NodeID>> Graph::flow() const {
        vector<pair<NodeID, NodeID>> f;
        for (const Node &n: node_data)
            f.push_back(pair(n.inflow, n.outflow));
        return f;
    }

    NodeID Graph::random_node() const {
        return nodes[rand() % nodes.size()];
    }

    pair <NodeID, NodeID> Graph::random_pair(size_t steps) const {
        if (steps < 1)
            return make_pair(random_node(), random_node());
        NodeID start = random_node();
        NodeID stop = start;
        for (size_t i = 0; i < steps; i++) {
            NodeID n = NO_NODE;
            do {
                n = util::random(node_data[stop].neighbors).node;
            } while (!contains(n));
            stop = n;
        }
        return make_pair(start, stop);
    }

    void Graph::randomize() {
        shuffle(nodes.begin(), nodes.end(), default_random_engine());
        for (NodeID node: nodes)
            shuffle(node_data[node].neighbors.begin(), node_data[node].neighbors.end(), default_random_engine());
    }

    void print_graph(const Graph &g, ostream &os) {
        vector<Edge> edges;
        g.get_edges(edges);
        sort(edges.begin(), edges.end());
        os << "p sp " << Graph::super_node_count() << " " << edges.size() << endl;
        for (Edge e: edges)
            os << "a " << e.a << ' ' << e.b << ' ' << e.d << endl;
    }

    void read_graph(Graph &g, istream &in) {
        char line_id;
        uint32_t v, w, d;

        while (in >> line_id) {
            switch (line_id) {
                case 'p':
                    in.ignore(3);
                    in >> v;
                    in.ignore(1000, '\n');
                    g.resize(v);
                    break;
                case 'a':
                    in >> v >> w >> d;
                    g.add_edge(v, w, d, true);
                    break;
                default:
                    in.ignore(1000, '\n');
            }
        }
        g.remove_isolated();
    }

//--------------------------- ostream -------------------------------

// for easy distance printing
    struct Dist {
        distance_t d;

        Dist(distance_t d) : d(d) {}
    };

    static ostream &operator<<(ostream &os, Dist distance) {
        if (distance.d == infinity)
            return os << "inf";
        else
            return os << distance.d;
    }

// for easy bit string printing
    struct BitString {
        uint64_t bs;

        BitString(uint64_t bs) : bs(bs) {}
    };

    static ostream &operator<<(ostream &os, BitString bs) {
        size_t len = bs.bs & 63ul;
        uint64_t bits = bs.bs >> 6;
        for (size_t i = 0; i < len; i++) {
            os << (bits & 1);
            bits >>= 1;
        }
        return os;
    }

    ostream &operator<<(ostream &os, const CutIndex &ci) {
        return os << "CI(p=" << bitset<8>(ci.partition) << ",c=" << (int) ci.cut_level
                  << ",di=" << ci.dist_index << ",d=" << ci.distances << ")";
    }

    ostream &operator<<(ostream &os, const FlatCutIndex &ci) {
        uint64_t partition_bitvector = *ci.partition_bitvector();
        vector<uint16_t> dist_index(ci.dist_index(), ci.dist_index() + ci.cut_level() + 1);
        vector<distance_t> distances(ci.distances(), ci.distances() + ci.label_count());
        return os << "FCI(pb=" << BitString(partition_bitvector) << ",di=" << dist_index << ",d=" << distances << ")";
    }

    ostream &operator<<(ostream &os, const ContractionLabel &cl) {
        return os << "CL(" << cl.cut_index << ",d=" << cl.distance_offset << ",p=" << cl.root << ")";
    }

    ostream &operator<<(ostream &os, const Neighbor &n) {
        if (n.distance == 1)
            return os << n.node;
        else
            return os << n.node << "@" << Dist(n.distance);
    }

    std::ostream &operator<<(std::ostream &os, const composite_shortcut &s) {
        if (s.triangle_node != NO_NODE)
            return os << "CS(" << s.triangle_node << "," << *s.lower << "," << *s.upper << ")";
        else
            return os << "[]";
    }

    std::ostream &operator<<(std::ostream &os, const path_data &p) {
        if (p.cs._no_node == NO_NODE) // cs is in use
            return os << p.cs;
        vector<NodeID> nodes;
        for (size_t i = 0; i < MAX_VALLEY_PATH_LENGTH; i++)
            if (p.intermediate[i] != NO_NODE)
                nodes.push_back(p.intermediate[i]);
        return os << nodes;
    }

    std::ostream &operator<<(std::ostream &os, const CHNeighbor &n) {
        return os << "S(" << n.neighbor << "," << n.distance << "," << n.p << ")";
    }

    std::ostream &operator<<(std::ostream &os, const QNode &n) {
        return os << "Q(" << n.shortcut_used << "," << n.shortcut_start << "," << n.node_id << "," << n.distance << ")";
    }

    ostream &operator<<(ostream &os, const Node &n) {
        return os << "N(" << n.subgraph_id << "#" << n.neighbors << ")";
    }

    ostream &operator<<(ostream &os, const Partition &p) {
        return os << "P(" << p.left << "|" << p.cut << "|" << p.right << ")";
    }

    ostream &operator<<(ostream &os, const DiffData &dd) {
        return os << "D(" << dd.node << "@" << dd.dist_a << "-" << dd.dist_b << "=" << dd.diff() << ")";
    }

    ostream &operator<<(ostream &os, const Graph &g) {
#ifdef MULTI_THREAD
        g.node_data.normalize();
#endif
        return os << "G(" << g.subgraph_id << "#" << g.nodes << " over " << g.node_data << ")";
    }

} // road_network
