#pragma once

#include <routingkit/nested_dissection.h>
#include <routingkit/constants.h>
#include <routingkit/graph_util.h>
#include <routingkit/permutation.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/filter.h>
#include <routingkit/id_mapper.h>
#include <routingkit/bit_vector.h>
#include <routingkit/timer.h>

#include "Tools/Workarounds.h"

namespace RoutingKit {

    // Debug method copied from RoutinKit/include/routingkit/nested_dissection.h
    void assert_fragment_is_valid(const GraphFragment &fragment) {
        unused(fragment);
#ifndef NDEBUG
        unsigned node_count = fragment.node_count();
        unsigned arc_count = fragment.arc_count();

        assert(fragment.tail.size() == arc_count);
        assert(fragment.head.size() == arc_count);
        assert(fragment.back_arc.size() == arc_count);
        assert(fragment.global_node_id.size() == node_count);

        if (arc_count != 0) {
            assert(max_element_of(fragment.tail) < node_count);
            assert(max_element_of(fragment.head) < node_count);
        }

        assert(invert_inverse_vector(fragment.first_out) == fragment.tail);

        for (unsigned a = 0; a < arc_count; ++a) {
            assert(fragment.back_arc[a] < arc_count);
            assert(fragment.back_arc[fragment.back_arc[a]] == a);
            assert(fragment.tail[a] == fragment.head[fragment.back_arc[a]]);
            assert(fragment.head[a] == fragment.tail[fragment.back_arc[a]]);
        }
#endif
    }

    GraphFragment build_vertex_induced_subgraph_by_copy(const GraphFragment &fragment, const BitVector &keep_vertex) {

        KASSERT(fragment.node_count() == keep_vertex.size());
        auto subgraph = fragment;

        // Remove arcs that do not belong to subgraph
        BitVector is_subgraph_arc = make_bit_vector(
                fragment.arc_count(),
                [&](unsigned a) {
                    return keep_vertex.is_set(fragment.tail[a]) && keep_vertex.is_set(fragment.head[a]);
                }
        );

        inplace_keep_element_of_vector_if(is_subgraph_arc, subgraph.tail);
        inplace_keep_element_of_vector_if(is_subgraph_arc, subgraph.head);
        inplace_keep_element_of_vector_if(is_subgraph_arc, subgraph.back_arc);

        // Map back-arcs to new position of edge
        {
            LocalIDMapper map(is_subgraph_arc);
            for (auto &x: subgraph.back_arc)
                x = map.to_local(x);
        }

        // Remove global_node_id entries of vertices that are not in subgraph
        inplace_keep_element_of_vector_if(keep_vertex, subgraph.global_node_id);

        // Map vertex IDs in edges to new vertex IDs:
        {
            LocalIDMapper map(keep_vertex);
            for (auto &t: subgraph.tail)
                t = map.to_local(t);
            for (auto &h: subgraph.head)
                h = map.to_local(h);
        }

        // Restore first_out based on tail array
        subgraph.first_out = invert_vector(subgraph.tail, keep_vertex.population_count());
        KASSERT(subgraph.node_count() == keep_vertex.population_count());

        assert_fragment_is_valid(subgraph);
        return subgraph;
    }

    std::vector<GraphFragment> decompose_into_at_most_two_subgraphs_along_separator(GraphFragment &&fragment, BitVector &&is_on_left_side, BitVector&& is_on_right_side) {

        // Split graph into two sub-graphs along cut.
        std::vector<GraphFragment> res;
        res.push_back(build_vertex_induced_subgraph_by_copy(fragment, is_on_left_side));
        if (res.back().node_count() == 0)
            res.pop_back();
        res.push_back(build_vertex_induced_subgraph_by_copy(fragment, is_on_right_side));
        if (res.back().node_count() == 0)
            res.pop_back();

        return res;
    }

    // Modification of compute_separator_decomposition from RoutingKit/include/routingkit/nested_dissection.h.
    // Computes separator decomposition using strict dissection, i.e., each separator has exactly two children in the
    // nested decomposition. Sub-graphs may not be connected (unlike in general separator decomposition) but the
    // balance of the decomposition is improved.
    SeparatorDecomposition compute_separator_decomposition_with_strict_dissection(
            GraphFragment&& fragment,
            const std::function<CutSide(const GraphFragment &)> &compute_cut,
            BitVector&& prev_is_separator_node,
            BitVector&& prev_cut_side,
            const std::function<void(const std::string &)> &log_message = [](const std::string &) {}
    ) {
        assert_fragment_is_valid(fragment);

        long long timer = 0;

        SeparatorDecomposition decomp;
        decomp.order.resize(fragment.node_count());

        if (fragment.node_count() == 1) {
            decomp.tree.push_back({0, 0, 0, 1});
            decomp.order = std::move(fragment.global_node_id);
        } else {

            unsigned pred = 0;
            unsigned order_begin = 0, order_end = fragment.node_count();

            decomp.tree.push_back({0, 0, 0, order_end});

            // Enter all separator nodes into order:
            for (auto i = 0; i < prev_is_separator_node.size(); ++i) {
                if (prev_is_separator_node.is_set(i))
                    decomp.order[--order_end] = fragment.global_node_id[i];
            }

            // Decompose graph into two subgraphs
            auto parts = decompose_into_at_most_two_subgraphs_along_separator(std::move(fragment), ~prev_cut_side & ~prev_is_separator_node, prev_cut_side & ~prev_is_separator_node);
            for (auto &part: parts) {
                assert(part.node_count() != 0);
                if (part.node_count() == 1) {
                    decomp.order[--order_end] = part.global_node_id[0];
                } else {
                    if (log_message && part.node_count() > 1000) {
                        log_message("Computing decomposition for top level component with " +
                                    std::to_string(part.node_count()) + " nodes");
                        timer = -get_micro_time();
                        log_message("Start computing top level separator");
                    }
                    auto cut = compute_cut(part);
                    auto is_separator_node = derive_separator_from_cut(part, cut.is_node_on_side);
                    if (log_message && part.node_count() > 1000) {
                        timer += get_micro_time();
                        log_message("Finished computing top level separator, its size is " +
                                    std::to_string(is_separator_node.population_count()) + " nodes needed " +
                                    std::to_string(timer) + "musec");
                    }

                    BitVector keep_arc = make_bit_vector(
                            part.arc_count(),
                            [&](unsigned a) {
                                return !is_separator_node.is_set(part.tail[a]) &&
                                       !is_separator_node.is_set(part.head[a]);
                            }
                    );

                    inplace_keep_element_of_vector_if(keep_arc, part.tail);
                    inplace_keep_element_of_vector_if(keep_arc, part.head);
                    inplace_keep_element_of_vector_if(keep_arc, part.back_arc);

                    {
                        LocalIDMapper map(keep_arc);
                        for (auto &x: part.back_arc)
                            x = map.to_local(x);
                    }

                    part.first_out = invert_vector(part.tail, part.node_count());
                    assert_fragment_is_valid(part);


                    if (log_message && part.node_count() > 1000) {
                        timer = -get_micro_time();
                        log_message("Start computing remaining separator decomposition using recursion");
                    }

                    auto sub_decomp = compute_separator_decomposition_with_strict_dissection(std::move(part), compute_cut, std::move(is_separator_node), std::move(cut.is_node_on_side));
                    if (log_message && part.node_count() > 1000) {
                        timer += get_micro_time();
                        log_message("Finished recursion, needed " + std::to_string(timer) + "musec");
                    }

                    for (auto &node: sub_decomp.tree) {
                        if (node.left_child != 0)
                            node.left_child += decomp.tree.size();
                        if (node.right_sibling != 0)
                            node.right_sibling += decomp.tree.size();
                        node.first_separator_vertex += order_begin;
                        node.last_separator_vertex += order_begin;
                    }
                    if (pred == 0)
                        decomp.tree[pred].left_child = decomp.tree.size();
                    else
                        decomp.tree[pred].right_sibling = decomp.tree.size();
                    pred = decomp.tree.size();
                    decomp.tree.insert(decomp.tree.end(), sub_decomp.tree.begin(), sub_decomp.tree.end());
                    std::copy(sub_decomp.order.begin(), sub_decomp.order.end(), decomp.order.begin() + order_begin);
                    order_begin += sub_decomp.order.size();
                }
            }
            decomp.tree[0].first_separator_vertex = order_begin;
        }

        return decomp; // NVRO
    }

}