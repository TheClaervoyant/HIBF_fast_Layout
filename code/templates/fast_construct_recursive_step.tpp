#include "../fast_construct_graph.h"

template <typename Hasher>
void recursive_step(const std::vector<std::vector<std::uint64_t>>& sigs, lemon::ListGraph& graph, 
    std::unordered_map<std::vector<size_t>,lemon::ListGraph::Node, Hasher>& labMap, 
    const std::vector<size_t>& comp, 
    const std::vector<std::pair<size_t,size_t>>& lvls, size_t lvl_){

    if(comp.size() == 1) return; // Can't work on a singleton.
    if(lvl_ >= lvls.size()) return; // Can't do more than that

    auto [bcount, bwidth] = lvls[lvl_];
    
    lemon::ListGraph tmp;
    std::unordered_map<size_t, lemon::ListGraph::Node> tmp_Map;
    construct_graph(comp, tmp, tmp_Map);
    connect_bands(sigs, comp, bcount, bwidth, tmp, tmp_Map, standardHasher{});

    std::vector<std::vector<size_t>> connect_comps_next_lvl = connected_components(tmp, tmp_Map);

    lemon::ListGraph::Node parent = (labMap.find(comp))->second;
    for(const std::vector<size_t>& sub_comp : connect_comps_next_lvl){
        lemon::ListGraph::Node child = graph.addNode();
        labMap[sub_comp] = child;
        graph.addEdge(parent, child);

        recursive_step(sigs, graph, labMap, sub_comp, lvls, lvl_ + 1);
    }
}