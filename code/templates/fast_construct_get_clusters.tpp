#include "../fast_construct_graph.h"

template <typename Hasher>
std::vector<std::unordered_map<size_t,const std::vector<size_t>*>> get_clusters(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps){
    std::vector<std::unordered_map<size_t,const std::vector<size_t>*>> seq_clusters(labMaps.size()); // Since we don't plan on changing the labMaps now, we can simply safe the pointer to the cluster itself.

    for(size_t lvl = 0; lvl < labMaps.size(); lvl++){
        for(auto& [comp, node] : labMaps[lvl]){
            for(size_t seq : comp) seq_clusters[lvl][seq] = &comp;
        }
    }
    return seq_clusters;
}