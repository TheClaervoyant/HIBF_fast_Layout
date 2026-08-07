#include "../fast_construct_graph.h"

template <typename Hasher>
std::pair<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>> binning_core(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, 
                                            const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters,
                                            const std::vector<std::vector<std::uint64_t>>& fracmin_sketches,
                                            const double s, const size_t bins, const size_t t_max, const double f){

    size_t deepest_lvl = labMaps.size() - 1;
    size_t par_lvl = (deepest_lvl > 0) ? deepest_lvl -1 : deepest_lvl; // Safeguard
    size_t bin = 0;
    size_t merge_bins = 0;
    size_t split_bins = 0;

    std::vector<std::vector<size_t>> res(bins); 
    if(bins == 0) return {res, {0,0,0}};

    std::vector<size_t> track_fill(bins, 0); // Used later on to track the amount of elements in each bin.
    std::vector<std::unordered_set<std::uint64_t>> bin_sketches(bins); // We need to efficiently keep track of the growing fracmin sketch for every bin
    std::vector<const std::vector<size_t>*> candidates;
    std::vector<const std::vector<size_t>*> top_clusters;
    std::vector<const std::vector<std::uint64_t>*> supersketch; // We want to see the cardinality of all smaller clusters combined in order to determine how many merge_bins are truly needed. We don't want them to potentially spread across all bins.
    std::vector<size_t> bin_lvls(bins, deepest_lvl); // We want to store the level each bin represents. This enables us to freely climb the tree upwards.


    std::unordered_set<size_t> binned; // Always check which sequences are already assigned to a bin.
    std::unordered_map<const std::vector<size_t>*, size_t> parent_size;
    std::unordered_map<const std::vector<size_t>*, size_t> top_sizes;

    std::vector<std::unordered_set<const std::vector<size_t>*>> used_clusters_per_level(deepest_lvl + 1); // To check whether the mother cluster was already used on each level. 
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> bin_for_cluster(deepest_lvl + 1); // We want to save, for every cluster, is there already an element in a bin and if yes - which?
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> itteration_progress(deepest_lvl + 1); // We save for every cluster how far we've scanned it already.

    std::queue<size_t> valid_bins; // Instead of iterating over every bin, we track which are valid.
    std::priority_queue<std::pair<size_t,size_t>, std::vector<std::pair<size_t,size_t>>, std::greater<std::pair<size_t,size_t>>> still_fillable; // For the last loop, check which bins are still fillable.

    // @brief Given update the bin for a given cluster, if the cluster either doesn't exist or it is full.
    // @param lvl : The level on which the update occurs
    // @param cluster : The cluster to be updated
    // @param size_t : The potential new cluster
    auto update_bin = [&](size_t lvl, const std::vector<size_t>* cluster, size_t b){
        bin_for_cluster[lvl][cluster] = b;
    };

    // @brief Insert the sequence into the desired bin and update trackfill and the corresponding bin sketch.
    // @param seq : The Sequence ID (0,....n)
    // @param b : The bin where the sequence will be inserted.
    // @param force :If true, insert the sequence even if the bin was full.
    // @return true if the element was entered in the bin, false if not (overshoot).
    auto try_insert_sequence = [&](size_t seq, size_t b, bool force){
        std::unordered_set<std::uint64_t>& sketch = bin_sketches[b];
        std::vector<std::uint64_t> new_elems;

        for(std::uint64_t elem : fracmin_sketches[seq]) if(!sketch.count(elem)) new_elems.push_back(elem);

        size_t potential_size = static_cast<size_t>(static_cast<double>(sketch.size() + new_elems.size())/s);
        if(!res[b].empty() && potential_size > t_max && !force) return false; // We don't enter it here, since it would overshoot the bin.

        for(std::uint64_t new_elem : new_elems) sketch.insert(new_elem);
        res[b].push_back(seq);
        binned.insert(seq);

        track_fill[b] = potential_size;
        return true; // Managed to fill
    };

    // @brief Try inserting the sequence into the desired bin and update trackfill and the sketch. If the sequence does not fit in the bin, put it also in other bins until it was inserted.
    // @param seq : The Sequence ID (0,....n)
    // @param b : The bin where the sequence will be inserted.
    // @param naive : Question: are we looking for the nearest elements or not? If naive = true, we just take the nearest fitting bin
    // @note The sequence seq can be in multiple res[b], since it got split properly.
    // @note This function is redacted. Remove it.
    auto try_insert_sequence_split = [&](size_t seq, size_t start, bool naive){
        std::vector<std::uint64_t> desired(fracmin_sketches[seq].begin(), fracmin_sketches[seq].end());
        std::unordered_set<size_t> tried_bins;
        size_t b = start;
        bool placed = false;

        while(b < bins && !desired.empty()){
            tried_bins.insert(b);
            std::unordered_set<std::uint64_t>& sketch = bin_sketches[b];

            std::vector<std::uint64_t> new_elems;
            new_elems.reserve(desired.size());
            for(std::uint64_t elem : desired) if(!sketch.count(elem)) new_elems.push_back(elem);

            bool was_empty = res[b].empty(); // We want to insert either way at least one element, if the bin was empty to begin with.
            size_t fitted = 0; // How many elements were put into this bin without overshooting?
            for(; fitted < new_elems.size(); fitted++){
                size_t potential_size = static_cast<size_t>(static_cast<double>(sketch.size() + fitted + 1)/s);
                if(potential_size > t_max) break;
            }
            if(fitted == 0 && was_empty && !new_elems.empty()) fitted = 1;

            for(size_t i = 0; i < fitted; i++) sketch.insert(new_elems[i]);

            std::vector<std::uint64_t> still_desired;
            still_desired.reserve(desired.size());
            for(std::uint64_t elem : desired) if(!sketch.count(elem)) still_desired.push_back(elem);

            bool done = still_desired.empty();
            bool contributed = (fitted > 0) || done;
            desired.swap(still_desired); // For the next itteration

            
            if(contributed){
                if(std::find(res[b].begin(), res[b].end(), seq) == res[b].end()) res[b].push_back(seq);
                track_fill[b] = static_cast<size_t>(static_cast<double>(sketch.size())/s);
                placed = true;
            }
            
            if(desired.empty()) break;

            size_t next_b = bins;
            if(naive && b + 1 < bins){next_b = b+1;}
            // Just like in the Fallback mechanism, we want to first try and find the closest relative for the sequence.
            else{
                for(size_t lvl = deepest_lvl; lvl > 0 && next_b == bins; --lvl){
                    auto it = level_clusters[lvl].find(seq);
                    if(it == level_clusters[lvl].end()) continue;
                    auto iterator = bin_for_cluster[lvl].find(it->second);
                    if(iterator == bin_for_cluster[lvl].end()) continue;
                    size_t cand = iterator->second;
                    if(cand == b || tried_bins.count(cand) || track_fill[cand] >= t_max) continue; // Tried it already.
                    if(track_fill[cand] < t_max) next_b = cand;
                }

                // Also, just like in the Fallback mechanism; if we don't find a relative, we just enter it in the least empty one.
                if(next_b == bins){
                    size_t best_take = t_max + 1;
                    for(size_t cand = 0; cand < bins; cand++){
                        if(tried_bins.count(cand)) continue;
                        if(track_fill[cand] < t_max && track_fill[cand] < best_take){best_take = track_fill[cand]; next_b = cand;}
                    }
                }

                // There is no empty bin left, so we paste the last remaining part where we started
                if(next_b == bins){
                    size_t best_b = 0;
                    size_t best_take = track_fill[0];
                    for(size_t cand = 1; cand < bins; cand++){
                        if(track_fill[cand] < best_take){best_take = track_fill[cand]; best_b = cand;}
                    }

                    std::unordered_set<std::uint64_t>& sketch = bin_sketches[best_b];
                    for(std::uint64_t elem : desired) sketch.insert(elem);
                    if(std::find(res[best_b].begin(), res[best_b].end(), seq) == res[best_b].end()) res[best_b].push_back(seq);
                    track_fill[best_b] = static_cast<size_t>(static_cast<double>(sketch.size())/s);
                    placed = true;
                    desired.clear();
                    break;

                }
            }
            b = next_b;
        }
        if(placed) binned.insert(seq);
        return b;
    };


    /// ============================================ ///


    /// @note this is all setup, such as sorting the highest level and splitting sequences too large.

    for(auto& [component, node] : labMaps[0]){
        for(size_t seq : component){
            if(fracmin_sketches[seq].size() >= t_max){
                split_bins += 1;
                size_t split_bins_ = static_cast<size_t>(std::ceil(static_cast<double>(fracmin_sketches[seq].size())/(f * static_cast<double>(t_max))));
                size_t already_split = 0;
                split_bins_ = (split_bins_ > bins - bin) ? (bins - bin) : split_bins_;
                for(std::uint64_t elem : fracmin_sketches[seq]){
                    std::unordered_set<std::uint64_t>& sketch = bin_sketches[bin + already_split];
                    size_t sketch_size = track_fill[bin + already_split];
                    size_t potential_sketch_size = sketch_size + static_cast<size_t>(static_cast<double>(1)/s);
                    if(potential_sketch_size > t_max*f && already_split + 1 < split_bins_){res[bin + already_split].push_back(seq); already_split += 1; split_bins += 1;} 
                    bin_sketches[bin + already_split].insert(elem);
                    track_fill[bin + already_split] = track_fill[bin + already_split]  +  static_cast<size_t>(static_cast<double>(1)/s);
                }
                binned.insert(seq);
                bin += (already_split + 1);
                res[bin - 1].push_back(seq);
            }
        }
    }

    for(auto& [component, node] : labMaps[0]) top_clusters.push_back(&component);

    top_sizes.reserve(top_clusters.size());
    for(const std::vector<size_t>* clust : top_clusters){
        std::vector<const std::vector<std::uint64_t>*> subset;
        subset.reserve(clust->size());
        for(size_t seq : *clust){
            subset.push_back(&fracmin_sketches[seq]);
        } 

        size_t union_size = get_union_size_ptr(subset);
        size_t estimated_size = static_cast<std::uint64_t>(static_cast<double>(union_size)/s);
        top_sizes[clust] = estimated_size;

        if(estimated_size < t_max) for(size_t seq : *clust) supersketch.push_back(&fracmin_sketches[seq]);

    }

    std::sort(top_clusters.begin(), top_clusters.end(), [&top_sizes](const std::vector<size_t>* a, const std::vector<size_t>* b){
        size_t size_a = top_sizes[a];
        size_t size_b = top_sizes[b];

        if(size_a != size_b) return size_a < size_b;
        return a->front() < b->front(); // Deterministic Tie Breaker.
    });


    /// ============================================ ///


    /// @note This part is the merging
    /// We sort the Top level in ascending Order and search for the BIGGEST Cluster to merge and work from "front to end" by adding every cluster into bins.
    merge_bins = bin;
    size_t merge_start = merge_bins;
    auto boundary = std::partition_point(top_clusters.begin(), top_clusters.end(), [&top_sizes, t_max](const std::vector<size_t>* c){return top_sizes[c] < t_max;}); // Get the biggest cluster to be merged
    size_t allowed_merge = static_cast<size_t>(static_cast<double>(get_union_size_ptr(supersketch))/(s * static_cast<double>(t_max)));  // Take the supersketch size and divide it by s (Union size) and by t_max to know how much needs to be filled.
    size_t limit = (bins > merge_bins + allowed_merge) ? (merge_bins + allowed_merge) : bins;

    for(auto it = boundary; it != top_clusters.begin();){
        it = std::prev(it);
        size_t seen = 0;

        const std::vector<size_t>& cluster = **it;
        for(size_t merging : cluster){
            
            if(binned.count(merging)) continue;
            while(limit > merge_bins && !try_insert_sequence(merging, merge_bins, false)) merge_bins+= 1; // Find a merge-bin that can hold items.
            if(merge_bins >= limit) break;

            used_clusters_per_level[0].insert(&cluster);

            auto it_deep = level_clusters[deepest_lvl].find(merging);
            if(it_deep != level_clusters[deepest_lvl].end()){
                used_clusters_per_level[deepest_lvl].insert(it_deep->second);
                update_bin(deepest_lvl, it_deep->second, merge_bins);
            }
            update_bin(0, &cluster, merge_bins);
            seen += 1;
            
        }
        // Every allowed merge_bin is full. We will just enter every element now greedily in the least empty bin.
        for(; seen < cluster.size(); seen++){
            size_t merging = cluster[seen];
            size_t best_b = merge_bins;
            size_t best_fill = track_fill[bin];
            for(size_t b = bin + 1; b < limit; b++){if(track_fill[b] < best_fill){best_fill = track_fill[b]; best_b = b;}}
            try_insert_sequence(cluster[seen],best_b,true);
            bool overshoot = track_fill[best_b] > t_max;

            used_clusters_per_level[0].insert(&cluster);

            auto it_deep = level_clusters[deepest_lvl].find(merging);
            if(it_deep != level_clusters[deepest_lvl].end()){
                used_clusters_per_level[deepest_lvl].insert(it_deep->second);
                update_bin(deepest_lvl, it_deep->second, best_b);
            }
            if(!overshoot) update_bin(0, &cluster, best_b);
        }

    }


    /// ============================================ ///


    /// @note We are splitting Clusters that are too big here.
    while(bin < bins && !res[bin].empty()) bin += 1;
    size_t start = bin;
    size_t merge_end = start;

    for(auto& [component, node] : labMaps[deepest_lvl]){
        bool all_binned = std::all_of(component.begin(), component.end(), [&](size_t s){return binned.count(s) != 0;});
        if(all_binned) continue; 
        
        std::vector<const std::vector<std::uint64_t>*> subset;
        subset.reserve(component.size());
        for(size_t seq : component) if(!binned.count(seq)) subset.push_back(&fracmin_sketches[seq]);
        size_t union_size = get_union_size_ptr(subset);
        size_t estimator = static_cast<std::uint64_t>(static_cast<double>(union_size)/s);

        if(estimator <= t_max){ // Normal component, just insert into candidates.
            candidates.push_back(&component);
            continue;
        }

        // estimator > t_max -> Split
        for(size_t seq : component){
            while(bin < bins && !try_insert_sequence(seq,bin,false)) bin += 1;
            if(bin >= bins) break;

            for(size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++){
                auto it = level_clusters[lvl_].find(seq);
                if(it != level_clusters[lvl_].end()){
                    used_clusters_per_level[lvl_].insert(it->second);
                    update_bin(lvl_, it->second, bin);
                }
            }
        }
        candidates.push_back(&component);
    }

    /// ============================================ ///

    // @note this is setup for the seeding algorithm
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


    /// =========================================== ///

    /// @note This is the original seeding Algorithm. We enter small clusters by first picking a representative (cluster.front). 
    /// The Seeding variety is top down - we first want representatives from every mother cluster, then every child cluster can contain one representative, etc.

    std::vector<size_t> repr_progress(candidates.size(),0);
    for(size_t lvl = 0; lvl <= deepest_lvl && bin < bins; lvl++){
        for(size_t cl = 0; cl < candidates.size(); cl++){
            const std::vector<size_t>* cluster = candidates[cl];
            if(bin >= bins) break;

            size_t& idx = repr_progress[cl];
            while(idx < cluster->size() && binned.count((*cluster)[idx])) idx += 1;
            if(idx >= cluster->size()) continue; // Every Element already looked up.
            size_t representative = (*cluster)[idx];
            if(binned.count(representative)) continue; // Since we entered this representative already, we skip it, since we don't need it twice.

            auto it = level_clusters[lvl].find(representative);
            if(it == level_clusters[lvl].end()) continue;
            const std::vector<size_t>* super_cluster = it->second;

            if(used_clusters_per_level[lvl].count(super_cluster)) continue; // We used this supercluster already, but we want more diversity. So we skip this Cluster and thus, its representative.

            while(bin < bins && !res[bin].empty()) bin += 1;
            if(bin >= bins) break;

            try_insert_sequence(representative, bin, false);

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
    for(size_t b =  start; b < bins; b++){
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
            if(!try_insert_sequence(seq,b, false)) break; // That is too much
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

    for(size_t b = start; b < bins; b++){
        if(track_fill[b] < t_max) still_fillable.push({track_fill[b],b});
    }

    for(auto& [seq, _] : level_clusters[0]){
        if(binned.count(seq)) continue;

        bool entered = false;
        for(size_t lvl = deepest_lvl; lvl > 0; --lvl){
            auto it = level_clusters[lvl].find(seq);
            if(it == level_clusters[lvl].end()) continue;
            auto iterator = bin_for_cluster[lvl].find(it->second);
            if(iterator == bin_for_cluster[lvl].end()) continue;
            size_t b = iterator->second;
            if(b < start) continue; // We don't want to put stuff into merge / splitbins
            if(!try_insert_sequence(seq,b, false)) continue; // Maybe another level and its corresponding cluster.

            bool overshoot = track_fill[b] > t_max;

            if(!overshoot){
                for(size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++){
                    auto it_ = level_clusters[lvl_].find(seq);
                    if(it_ != level_clusters[lvl_].end()) update_bin(lvl_, it_->second, b);
                }
            }
            entered = true;
            break;
        }
        if(!entered){
            std::vector<std::pair<size_t, size_t>> skipped; // We want to know, which Buckets didn't fit the sequence but still have place for others potentially

            while(!still_fillable.empty()){
                auto [fill, b] = still_fillable.top();
                still_fillable.pop(); 
                if(track_fill[b] >= t_max) continue;

                if(try_insert_sequence(seq,b, false)){
                    entered = true;
                    bool overshoot = track_fill[b] > t_max;
                    if(!overshoot){
                        for(size_t lvl_ = deepest_lvl; lvl_ > 0; lvl_--){
                            auto it_ = level_clusters[lvl_].find(seq);
            
                            if(it_!= level_clusters[lvl_].end()) update_bin(lvl_, it_->second, b);
                    }
                }
                if(track_fill[b] < t_max) still_fillable.push({track_fill[b],b});
                break;
            }
            skipped.push_back({track_fill[b], b}); // Couldn't enter the sequence in that bin, but it HAS place
            } 

            for(auto& s : skipped) still_fillable.push(s);

            if(!entered){
                // FINAL FINAL WAY. We will enter every sequence and when there really was no bin able to hold it, we'll just overshoot one bin, but greedily the least full one.
                size_t best_b = start;
                size_t best_fill = track_fill[start];
                for(size_t b = start + 1; b < bins; b++){if(track_fill[b] < best_fill){best_fill = track_fill[b]; best_b = b;}}

                try_insert_sequence(seq, best_b, true);
                bool overshoot = track_fill[best_b] > t_max;
                if(!overshoot){
                    for(size_t lvl_ = deepest_lvl; lvl_ > 0; lvl_--){
                        auto it_ = level_clusters[lvl_].find(seq);
            
                        if(it_!= level_clusters[lvl_].end()) update_bin(lvl_, it_->second, best_b);
                    }
                }
                if(track_fill[best_b] < t_max) still_fillable.push({track_fill[best_b], best_b});
            }
        }
    }

    /// ============================================ ///

    std::sort(res.begin(), res.end(), [] (const std::vector<std::size_t>& a, const std::vector<std::uint64_t>& b){
        return a.size() < b.size();
    });

    auto merge_it = std::partition_point(res.begin(), res.end(), [](const std::vector<size_t>& bin){return bin.size() < 2;});
    auto split_it = std::partition_point(res.begin(), res.end(), [](const std::vector<size_t>& bin){return bin.size() == 0;}); // First Element with size 1; goes against empty bins

    merge_start = std::distance(res.begin(), merge_it);
    size_t split_start = std::distance(res.begin(), split_it);


    return {res, {split_start, split_bins, merge_start}};
}