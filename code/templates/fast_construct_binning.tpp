#include "../fast_construct_graph.h"

template <typename Hasher>
std::vector<std::vector<size_t>> binning(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters, const size_t bins, const size_t t_max){
    size_t deepest_lvl = labMaps.size() - 1;
    size_t par_lvl = (deepest_lvl > 0) ? deepest_lvl -1 : deepest_lvl; // Safeguard
    size_t bin = 0;
    size_t merge_bins = 0;
    
    std::vector<std::vector<size_t>> res(bins); 
    std::vector<size_t> track_fill(bins, 0); // Used later on to track the amount of elements in each bin.
    std::vector<const std::vector<size_t>*> candidates;
    std::vector<const std::vector<size_t>*> top_clusters;
    std::vector<size_t> bin_lvls(bins, deepest_lvl); // We want to store the level each bin represents. This enables us to freely climb the tree upwards.


    std::unordered_set<size_t> binned; // Always check which sequences are already assigned to a bin.
    std::unordered_map<const std::vector<size_t>*, size_t> parent_size;

    std::vector<std::unordered_set<const std::vector<size_t>*>> used_clusters_per_level(deepest_lvl + 1); // To check whether the mother cluster was already used on each level. 
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> bin_for_cluster(deepest_lvl + 1); // We want to save, for every cluster, is there already an element in a bin and if yes - which?
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> itteration_progress(deepest_lvl + 1); // We save for every cluster how far we've scanned it already.

    std::queue<size_t> valid_bins; // Instead of iterating over every bin, we track which are valid.
    std::priority_queue<std::pair<size_t,size_t>, std::vector<std::pair<size_t,size_t>>, std::greater<std::pair<size_t,size_t>>> still_fillable; // For the last loop, check which bins are still fillable.

    // @brief gieven update the bin for a given cluster, if the cluster either doesn't exist or it is full.
    // @param lvl : The level on which the update occurs
    // @param cluster : The cluster to be updated
    // @param size_t : The potential new cluster
    auto update_bin = [&](size_t lvl, const std::vector<size_t>* cluster, size_t b){
        auto old = bin_for_cluster[lvl].find(cluster);
        if(old == bin_for_cluster[lvl].end() || track_fill[old->second] >= t_max) bin_for_cluster[lvl][cluster] = b;

    };


    /// ============================================ ///


    /// @note this is all setup, such as sorting the highest and deepest lvl, saving parent sizes for the deepest lvl-
    for(auto& [component, node] : labMaps[deepest_lvl]) candidates.push_back(&component);
    for(auto& [component, node] : labMaps[0]) top_clusters.push_back(&component);

    std::sort(top_clusters.begin(), top_clusters.end(), [](const std::vector<size_t>* a, const std::vector<size_t>* b){
        if(a->size() != b->size()) return a->size() < b->size();
        return a->front() < b->front();
    });

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


    /// ============================================ ///


    /// @note This part is the merging
    /// We sort the Top level in ascending Order and search for the BIGGEST Cluster to merge and work from "front to end" by adding every cluster into bins.
    auto boundary = std::partition_point(top_clusters.begin(), top_clusters.end(), [t_max](const std::vector<size_t>* c){return c->size() < t_max;}); // Get the biggest cluster to be merged

    for(auto it = boundary; it != top_clusters.begin();){
        it = std::prev(it);

        while(merge_bins < bins && track_fill[merge_bins] >= t_max) merge_bins+= 1; // Find a merge-bin that can hold items.
        if(merge_bins >= bins) break;
        const std::vector<size_t>& cluster = **it;
        for(size_t merging : cluster){
            while(merge_bins < bins && track_fill[merge_bins] >= t_max) merge_bins+= 1; // Find a merge-bin that can hold items.
            if(merge_bins >= bins) break;
            if(binned.count(merging)) continue;

            res[merge_bins].push_back(merging);
            track_fill[merge_bins] += 1;
            binned.insert(merging);

            used_clusters_per_level[0].insert(&cluster);

            auto it_deep = level_clusters[deepest_lvl].find(merging);
            if(it_deep != level_clusters[deepest_lvl].end()){
                used_clusters_per_level[deepest_lvl].insert(it_deep->second);
                update_bin(deepest_lvl, it_deep->second, merge_bins);
            }
            update_bin(0, &cluster, merge_bins);
        }

    }

    
    /// ============================================ ///

    
    /// @note This is the original seeding Algorithm. We enter small clusters by first picking a representative (cluster.front). 
    /// The Seeding variety is top down - we first want representatives from every mother cluster, then every child cluster can contain one representative, etc.
    while(bin < bins && !res[bin].empty()) bin += 1;
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
            track_fill[bin] += 1;
            binned.insert(representative);

            for(size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++){
                auto iterator = level_clusters[lvl_].find(representative);
                if(iterator != level_clusters[lvl_].end()){
                    used_clusters_per_level[lvl_].insert(iterator->second);
                    update_bin(lvl_, iterator->second, bin);
                } 
            }

            bin += 1;
        }
    }


    /// ============================================ ///


    /// @note This part is the climbing mechanism. We take every representative and climb its cluster upwards by recursively going up one level after finishing. 
    /// This goes to lvl 0 (root), but elements are not taken from the mother clusters.
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
        update_bin(curr_lvl, it->second, b);

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


    /// ============================================ ///


    /// @note The Fallback mechanism.
    /// For every not entered sequence, we check whether there is already a cluster once checked with this element. If so, we can enter this element in the bin corresponding to the cluster.
    /// If there is no such cluster, take the emptiest bin and enter the sequence as a "new seed" and mark the cluster/s as seen for later.
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
                    if(it_ != level_clusters[lvl_].end()) update_bin(lvl_, it_->second, b);
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

                if(it_!= level_clusters[lvl_].end()) update_bin(lvl_, it_->second, b);
            }

            if(track_fill[b] < t_max) still_fillable.push({track_fill[b],b}); // Remove it ASAP.
        }
    }


    /// ============================================ ///


    return res;
}