#include "../fast_construct_graph.h"

template <typename Hasher>
std::vector<std::vector<size_t>> binning(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters, const size_t bins, const size_t t_max){
    std::vector<std::vector<size_t>> res(bins);
    std::vector<size_t> track_fill(bins, 0); // Used later on to track the amount of elements in each bin.
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
            track_fill[bin] = 1;
            binned.insert(representative);

            for(size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++){
                auto iterator = level_clusters[lvl_].find(representative);
                if(iterator != level_clusters[lvl_].end()) used_clusters_per_level[lvl_].insert(iterator->second);
            }

            bin += 1;
        }
    }

    // TODO: the primary focus here lies on filling every bucket. Change it, such that we fill rather bottom up from the tree; group more similar sequences together rather than being set on filling a bin.
    std::vector<size_t> bin_lvls(bins, deepest_lvl); // We want to store the level each bin represents. This enables us to freely climb the tree upwards.
    std::queue<size_t> valid_bins; // Instead of iterating over every bin, we track which are valid.
    for(size_t b =  0; b < bins; b++){
        if(!res[b].empty() && track_fill[b] < t_max) valid_bins.push(b);
    }

    while(!valid_bins.empty()){
        size_t b = valid_bins.front();
        valid_bins.pop();

        size_t repr = res[b][0];
        size_t curr_lvl = bin_lvls[b];

        // We can now use this representative to go up one level and add every not already used sequence to the bin.
        auto it = level_clusters[curr_lvl].find(repr);
        if(it == level_clusters[curr_lvl].end()){ // Maybe the representative can be found in a higher level, just to be safe.
            if(curr_lvl > 0){
                bin_lvls[b] -= 1;
                valid_bins.push(b); // Enables a new itteration
            }
            continue;
        } 
        const std::vector<size_t>& cluster = *(it->second);

        for(size_t seq : cluster){
            if(track_fill[b] >= t_max) break; // Will not add more than we already added.
            if(binned.count(seq)) continue;
            res[b].push_back(seq);
            binned.insert(seq);
            track_fill[b] += 1;
        }

        if(track_fill[b] < t_max && curr_lvl > 0){
            bin_lvls[b] -= 1;
            valid_bins.push(b);
        }
    }
    // TODO: We can't requeue here yet, since this would result in an endless loop. Implement an increasing level, such that always something new happens.

    // TODO: This goes naively and randomly enters its stuff. Maybe think of something smart instead of being the stupidest form of greedy algorithm.
    std::queue<size_t> still_fillable;
    for(size_t b = 0; b < bins; b++){
        if(track_fill[b] < t_max) still_fillable.push(b);
    }

    for(auto& [seq, _] : level_clusters[0]){
        if(binned.count(seq)) continue;

        while(!still_fillable.empty() && track_fill[still_fillable.front()] >= t_max) still_fillable.pop(); // We filled these buckets already along the way. 
        if(still_fillable.empty()) break; // Every bucket is full.

        size_t b = still_fillable.front();
        res[b].push_back(seq);
        binned.insert(seq);
        track_fill[b] += 1;

        if(track_fill[b] >= t_max) still_fillable.pop(); // Remove it ASAP.
    }
    return res;
}