#include "../fast_construct_graph.h"

template <typename Hasher>
std::vector<std::vector<size_t>> binning(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters, const size_t bins, const size_t t_max){
    std::vector<std::vector<size_t>> res(bins); 
    std::vector<size_t> track_fill(bins, 0); // Used later on to track the amount of elements in each bin.
    std::unordered_set<size_t> binned; // Always check which sequences are already assigned to a bin.

    // For starters, we want to assign the smallest clusters into bins so we can use them as representatives.
    // Since we also want them to be as distant as possible, we go into the deepest level and from there, add singletons and whatnot.

    size_t deepest_lvl = labMaps.size() - 1;
    size_t par_lvl = (deepest_lvl > 0) ? deepest_lvl -1 : deepest_lvl; // Safeguard

    std::vector<const std::vector<size_t>*> candidates;
    for(auto& [component, node] : labMaps[deepest_lvl]) candidates.push_back(&component);

    std::unordered_map<const std::vector<size_t>*, size_t> parent_size;
    parent_size.reserve(candidates.size());
    for(const std::vector<size_t>* cand: candidates){
        auto it = level_clusters[par_lvl].find(cand->front());
        parent_size[cand] = (it != level_clusters[par_lvl].end()) ? it->second->size() : cand->size(); // Save how big the corresponding super cluster is
    }
    // Sort them now in ascending order; this enables us to easily get the smallest clusters.
    std::sort(candidates.begin(),candidates.end(), [&parent_size](const std::vector<size_t>* a, const std::vector<size_t>* b){
        if(a->size() != b->size()) return a->size() < b->size();
        size_t size_a = parent_size[a];
        size_t size_b = parent_size[b];
        if(size_a != size_b) return size_a > size_b; // We want to prefer small clusters with a greater super cluster.
        return a->front() < b->front(); // Tie breaker to make binning deterministic.
    });

    std::vector<std::unordered_set<const std::vector<size_t>*>> used_clusters_per_level(deepest_lvl + 1); // To check whether the mother cluster was already used on each level. 
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> bin_for_cluster(deepest_lvl + 1); // We want to save, for every cluster, is there already an element in a bin and if yes - which?
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> itteration_progress(deepest_lvl + 1); // We save for every cluster how far we've scanned it already.

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
                if(iterator != level_clusters[lvl_].end()){
                    used_clusters_per_level[lvl_].insert(iterator->second);
                    bin_for_cluster[lvl_].try_emplace(iterator->second, bin);
                } 
            }

            bin += 1;
        }
    }

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
        if(curr_lvl == 0) continue;

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
        bin_for_cluster[curr_lvl].try_emplace(it->second, b);

        size_t& start_it = itteration_progress[curr_lvl][it->second];

        for(; start_it < cluster.size(); start_it++){
            if(track_fill[b] >= t_max) break; // Will not add more than we already added.
            size_t seq = cluster[start_it];
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

    // TODO: This goes naively and randomly enters its stuff. Maybe think of something smart instead of being the stupidest form of greedy algorithm.
    std::priority_queue<std::pair<size_t,size_t>, std::vector<std::pair<size_t,size_t>>, std::greater<std::pair<size_t,size_t>>> still_fillable;
    for(size_t b = 0; b < bins; b++){
        if(track_fill[b] < t_max) still_fillable.push({track_fill[b],b});
    }

    for(auto& [seq, _] : level_clusters[0]){
        if(binned.count(seq)) continue;

        bool entered = false;
        for(size_t lvl = deepest_lvl; lvl > 0; --lvl){
            auto it = level_clusters[lvl].find(seq);
            if(it == level_clusters[lvl].end()) continue;
            auto iterator = bin_for_cluster[lvl].find(it->second);
            if(iterator != bin_for_cluster[lvl].end() && track_fill[iterator->second] < t_max){
                size_t b = iterator->second;
                res[b].push_back(seq);
                binned.insert(seq);
                track_fill[b] += 1;

                for(size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++){
                    auto it_ = level_clusters[lvl_].find(seq);
                    if(it_ != level_clusters[lvl_].end()) bin_for_cluster[lvl_].try_emplace(it_->second, b);
                }
                entered = true;
                break;
            }
        }
        if(!entered){
            while(!still_fillable.empty() && track_fill[still_fillable.top().second] >= t_max) still_fillable.pop(); // We filled these buckets already along the way. 
            if(still_fillable.empty()) break; // Every bucket is full.
            size_t b = still_fillable.top().second;
            still_fillable.pop();
            res[b].push_back(seq);
            binned.insert(seq);
            track_fill[b] += 1;

            for(size_t lvl_ = deepest_lvl; lvl_ > 0; lvl_--){
                auto it_ = level_clusters[lvl_].find(seq);
                auto cluster = it_->second;

                auto old_bin = bin_for_cluster[lvl_].find(cluster);

                if(old_bin == bin_for_cluster[lvl_].end() || track_fill[old_bin->second] >= t_max) bin_for_cluster[lvl_][cluster] = b;
            }

            if(track_fill[b] < t_max) still_fillable.push({track_fill[b],b}); // Remove it ASAP.
        }
    }

    return res;
}