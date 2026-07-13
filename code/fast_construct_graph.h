#pragma once

#include "fast_construct_hashing.h"
#include <lemon/list_graph.h>
#include <map>
#include <tuple>
#include <unordered_map>
#include <utility>

// @brief Transform the given LableMap and the Graph based on the given bands and sequences.
// @param  bands_, seqs_ : The number of bands and sequences used.
// @param graph : The Graph to be construct_graphed
// @param LabelMap : A corresponding LableMap. It maps every Lable to a corresponding Node
//
// @note For the following, the Mapping Node -> Lable is bijective, so we can call Lable -> Node as well.
// @note Nodes (seqi,bandi), (seqj,bandj) only connected if seqi == seqj
// @note Number_Nodes = seqs_*bands_
void construct_graph_tower(const size_t seqs_, const size_t bands_,
                           lemon::ListGraph &graph,
                           std::map<std::tuple<size_t, size_t>, lemon::ListGraph::Node> &LabelMap);

// @brief Transform the given LabelMap and the Graph based on the given sequences.
// @param seqs_ : The number of sequences used.
// @param graph : The Graph to be construct_graphed
// @param LabelMap : A corresponding LableMap. It maps every Lable to a corresponding Node
//
// @note For the following, the Mapping Node -> Lable is bijective, so we can call Lable -> Node as well.
// @idea We want to represent every Sequence here as one Node. This might be computationally more efficient, since we just have seqs_ Nodes instead of seqs_*bands_. This speeds up connected components.
void construct_graph(const size_t seqs_,
                     lemon::ListGraph &graph,
                     std::unordered_map<size_t, lemon::ListGraph::Node> &LabelMap);

// @brief Transform the given LabelMap and Graph based on the given Labels.
// @param labs_ : The Labels that the Nodes will have.
// @param graph : The graph to transform.
// @param LabelMap : The LabelMap for the corresponding graph.
//
// @note use this overload for deeper levels in the Datastructure, because the Labels now might not be linear anymore (e.g. the level can contain sequences {1,3,5,6,7}).
void construct_graph(const std::vector<size_t> &labs_, lemon::ListGraph &graph, std::unordered_map<size_t, lemon::ListGraph::Node> &LabelMap);

// @brief Transform the given LabelMap and Graph based on the sets of Labels
// @param labs_ : The sets of Labels. Each Node is represented by one set.
// @param graph : The graph to transform.
// @param LabelMap : The LabelMap for the Graph
//
// @note use this overload for representing sets of labels (e.g. connected components) as one Node
template <typename Hasher> 
void construct_graph(const std::vector<std::vector<size_t>> &labs_, lemon::ListGraph &graph, std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher> &LabelMap);

#include "templates/fast_construct_graph_comp.tpp"

struct standardHasher {
    template <typename T>
    std::uint64_t operator()(const std::vector<T>& vec) const {return hash_vector(vec);}
};

struct debugHasher {
    template <typename T>
    std::uint64_t operator()(const std::vector<T>& vec) const {return hash_vector_debug(vec);}
};

// @brief Given a set of One-Permutation-Hashes and a corresponding graph, after setting the bandwidth, connect two Nodes when the bands are the same
// @param hashes : the set of OPHS computed before
// @param bwidth : The size of one band (e.g. bwidth = 5 would take the first 5 values and identify them as one band).
// @param bcount : The number of bands to be used in computation.
// @param graph : the graph to be modified
// @param LabelMap : Mapping Label -> Node
//
// @note hashes[i].size() > bwidth * bcount for any i
// @note not using the OneNode computation, so bcount = bands_
void connect_bands_tower(const std::vector<std::vector<std::uint64_t>> &hashes,
                         const size_t bcount, const size_t bwidth,
                         lemon::ListGraph &graph,
                         std::map<std::tuple<size_t, size_t>, lemon::ListGraph::Node> &LabelMap,
                         std::function<std::uint64_t(const std::vector<std::uint64_t> &)> hashFunc = standardHasher{});
                         
// @brief Given a set of One-permutation-Hashes and the information, which one are important, connect two bands when they are the same, considering bandwidth/ and count
// @param hashes : A set of OPHS (or comparable)
// @param sig_labs : A set describing, which Hashes should be considered
// @param bwidth : The size of one band (e.g. bwidth = 5 would take the first 5 values and identify them as one band).
// @param bcount : The number of bands to be used in computation.
// @param graph : the graph to be modified
// @param LabelMap : Mapping Label -> Node
// @param hasher : A constructor containing a Hash Functions for Vectors, as they will be the keys
//
// @note use this override for deeper levels, as relevant sequences might not be linear. On first level, use function above. It is, however just a wrapper for this function
// @note sig_labs.size() == LabelMap.count()
template <typename Hasher>
void connect_bands(const std::vector<std::vector<std::uint64_t>> &hashes,
                const std::vector<size_t> &sig_labs,
                const size_t bcount, const size_t bwidth,
                lemon::ListGraph &graph,
                std::unordered_map<size_t, lemon::ListGraph::Node> &LabelMap,
                Hasher hasher = standardHasher{});
#include "templates/fast_construct_connect_bands_given_seqs.tpp"

// @brief Given a set of One-Permutation-Hashes and a corresponding graph, after setting the bandwidth, connect two Nodes when the bands are the same
// @param hashes : the set of OPHS computed before
// @param bwidth : The size of one band (e.g. bwidth = 5 would take the first 5 values and identify them as one band).
// @param bcount : The number of bands to be used in computation.
// @param graph : the graph to be modified
// @param LabelMap : Mapping Label -> Node
// @param hasher : A constructor containing a Hash Functions for Vectors, as they will be the keys
//
// @note hashes[i].size() > bwidth * bcount for any i
// @note The implementation is nearly as same as connect_bands_tower with the only difference being the call of the singular Nodes in the LabelMap
template <typename Hasher>
void connect_bands_all_seqs(const std::vector<std::vector<std::uint64_t>> &hashes,
                   const size_t bcount, const size_t bwidth,
                   lemon::ListGraph &graph,
                   std::unordered_map<size_t, lemon::ListGraph::Node> &LabelMap,
                   Hasher hasher = standardHasher{});
#include "templates/fast_construct_connect_bands_all_seqs.tpp"


// @brief Gather every Disjoint Subgraph of the given Graph. Here, every vector in the vector represents one Subgraph by the given Nodes.
// @param graph : The graph the connected components should be gathered.
// @param labMap : The Label Map corresponding, in order to reconstruct the nodeMap.
//
// @note The given Nodes are still the same as the ones used the original graph. Lables given before that will yet apply.
// @note function might be minutely changed, st. the LableMap is adjusted (split into multiple maps with every map just working for every component) or and Index, Map returned (Mapping Lable -> Component Number)
std::vector<std::vector<size_t>> connected_components(lemon::ListGraph &graph, const std::unordered_map<size_t, lemon::ListGraph::Node> &labMap);

// @brief connect the parent node with all its children nodes gathered in the level below
// @param sigs : All Hash Signatures
// @param graph : The graph containing the parent node to be connected
// @param labMap : The graphs corresponding lableMap
// @param comp : The connected component to be worked on recursively
// @param lvls : a pair, containing the tuple (bcount, bwidth) for each lvl
// @param lvl_ : A global accumulator
template <typename Hasher>
void recursive_step(const std::vector<std::vector<std::uint64_t>> &sigs, lemon::ListGraph &graph, std::vector<std::unordered_map<std::vector<size_t>,lemon::ListGraph::Node, Hasher>>& labMaps,
                    const std::vector<size_t> &comp,
                    const std::vector<std::pair<size_t,size_t>>& lvls, size_t lvl_);
#include "templates/fast_construct_recursive_step.tpp"

// @brief Generate the Graph based on the Signatures and the levels wanted.
// @param sigs : The Signatures (OPH) used for the sequences
// @param lvls : A set of tuples (bcount, bwidth). Every tuple defines a level in the tree
// @param graph : The graph upon which is build.
template <typename Hasher>
std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>> generate_all(const std::vector<std::vector<std::uint64_t>>& sigs, const std::vector<std::pair<size_t,size_t>>& lvls,
                 lemon::ListGraph& graph);
#include "templates/fast_construct_generate_all.tpp"