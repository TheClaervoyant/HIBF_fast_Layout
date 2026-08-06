#include "../fast_construct_graph.h"

template <typename Hasher>
LSH_Filtered<Hasher> filter_LSH(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, const std::vector<std::unordered_map<size_t, const std::vector<size_t>*>>& level_clusters, const std::vector<size_t>& relevant_seqs){
    const size_t deepest_lvl = level_clusters.size()-1;

    LSH_Filtered<Hasher> res;
    res.filtered_labMaps.resize(deepest_lvl + 1);
    res.filtered_level_clusters.resize(deepest_lvl + 1);
    res.filtered_cluster_storage.resize(deepest_lvl + 1);

    std::unordered_set<size_t> relevants(relevant_seqs.begin(), relevant_seqs.end()); // Enables efficient lookup

    for(size_t lvl = 0; lvl <= deepest_lvl; lvl++){
        auto& storage = res.filtered_cluster_storage[lvl];
        storage.reserve(level_clusters[lvl].size());

        std::unordered_set<const std::vector<size_t>*> seen; // We just need to know which cluster we've seen already.

        for(const auto& [seq, cluster] : level_clusters[lvl]){
            if(!seen.insert(cluster).second) continue;

            std::vector<size_t> filtered;
            filtered.reserve(cluster->size());

            for(size_t s : *cluster){
                if(relevants.count(s)) filtered.push_back(s);
            }

            if(filtered.empty()) continue;

            storage.push_back(std::move(filtered));

            const std::vector<size_t>* filtered_ptr = &storage.back();

            for(size_t s : *filtered_ptr) res.filtered_level_clusters[lvl][s] = filtered_ptr;
        }
    }

    for(size_t lvl = 0; lvl <= deepest_lvl; ++lvl){
        for(const auto& [cluster, node] : labMaps[lvl]){
            std::vector<size_t> filtered;
            filtered.reserve(cluster.size());

            for(size_t s : cluster){
                if(relevants.count(s)) filtered.push_back(s);
            }

            if(filtered.empty()) continue;

            res.filtered_labMaps[lvl].emplace(std::move(filtered), node);
        }
    }
    return res;
}