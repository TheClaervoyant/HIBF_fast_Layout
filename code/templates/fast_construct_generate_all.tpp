#include "../fast_construct_graph.h"

template <typename Hasher>
std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>> generate_all(const std::vector<std::vector<std::uint64_t>>& sigs, const std::vector<std::pair<size_t,size_t>>& lvls,lemon::ListGraph& graph){
    std::unordered_map<size_t, lemon::ListGraph::Node> labMap0;
    lemon::ListGraph graph0;

    auto [bcount0, bwidth0] = lvls[0];
    construct_graph(sigs.size(), graph0, labMap0);
    connect_bands_all_seqs(sigs, bcount0, bwidth0, graph0, labMap0, standardHasher{});

    std::vector<std::vector<size_t>> con_comps = connected_components(graph0, labMap0);

    std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>> labMaps(lvls.size());
    construct_graph(con_comps, graph, labMaps[0]);

    for(const std::vector<size_t>& comp : con_comps)
        recursive_step(sigs, graph, labMaps, comp, lvls, 1);
    return labMaps;
}