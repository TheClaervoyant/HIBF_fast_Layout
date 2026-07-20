#include "../fast_construct_graph.h"

template <typename Hasher>
std::vector<std::vector<size_t>> binning(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters, const size_t bins, const size_t t_max){
    std::vector<std::vector<size_t>> res(bins);
    std::unordered_set<size_t> binned; // Always check which sequences are already assigned to a bin.

    // For starters, we want to assign the smallest clusters into bins so we can use them as representatives.
    // Since we also want them to be as distant as possible, we go into the deepest level and from there, add singletons and whatnot.

    size_t deepest_lvl = labMaps.size() - 1;

    std::vector<const std::vector<size_t>*> candidates;
    for(auto& [component, node] : labMaps[deepest_lvl]) candidates.push_back(&component);

    // Sort them now in ascending order; this enables us to easily get the smallest clusters.
    std::sort(candidates.begin(),candidates.end(), [](const std::vector<size_t>* a, const std::vector<size_t>* b){
        if(a->size() != b->size()) return a->size() < b->size();
        return a->front() < b->front(); // Tie breaker to make binning deterministic.
    });

    std::vector<std::unordered_set<const std::vector<size_t>*>> used_clusters_per_level(deepest_lvl + 1); // To check whether the mother cluster was already used on each level. 

    size_t bin = 0;

    // TODO: This version only goes up one level to check whether a cluster was already used, this could be changed to make the cluster check top down, such that two sequences come from very different clusters.
    for(size_t lvl = 0; lvl <= deepest_lvl && bin < bins; lvl++){
        for(const std::vector<size_t>* cluster : candidates){
            if(bin >= bins) break;
    
            size_t representative = cluster->front();
            if(binned.count(representative)) continue; // Since we entered this representative already, we skip it, since we don't need it twice.
    
            auto it = level_clusters[lvl].find(representative);
            if(it == level_clusters[lvl].end()) continue;
            const std::vector<size_t>* super_cluster = it->second;

            if(used_clusters_per_level[lvl].count(super_cluster)) continue; // We used this supercluster already, but we want more diversity. So we skip this Cluster and thus, its representative.

            res[bin].push_back(representative);
            binned.insert(representative);

            for(size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++){
                auto iterator = level_clusters[lvl_].find(representative);
                if(iterator != level_clusters[lvl_].end()) used_clusters_per_level[lvl_].insert(iterator->second);
            }

            bin += 1;
        }
    }

    for(size_t b = 0; b < bins; b++){
        if(res[b].empty()) continue; // Can't add similar sequences to an empty bin, since there is no representative.
        if(res[b].size() >= t_max) continue; // We will not add elements to an already full bin.
        size_t repr = res[b][0];

        // We can now use this representative to go up one level and add every not already used sequence to the bin.
        auto it = level_clusters[deepest_lvl].find(repr);
        if(it == level_clusters[deepest_lvl].end()) continue; // If, for some reason, the representative should not be found, use this as a guard.
        const std::vector<size_t>& cluster = *(it->second);

        for(size_t seq : cluster){
            if(res[b].size() >= t_max) continue; // Will not add more than we already added.
            if(binned.count(seq)) continue;
            res[b].push_back(seq);
            binned.insert(seq);
        }
    }

    // TODO: This goes naively and randomly enters its stuff. Maybe think of something smart instead of being the stupidest form of greedy algorithm.
    for(auto& [seq, _] : level_clusters[0]){
        if(binned.count(seq)) continue;
        for(size_t b = 0; b < bins; b++){
            if(res[b].size() < t_max){
                res[b].push_back(seq);
                binned.insert(seq);
                break;
            }
        }
    }
    return res;
}