#include "../fast_construct_bin.h"

template <typename Hasher>
std::tuple<std::vector<std::vector<IBF>>, std::vector<std::vector<std::tuple<size_t,size_t,size_t>>>, std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>>, std::vector<std::vector<size_t>>, std::vector<std::vector<std::pair<size_t,size_t>>>> generate_hibf(const std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>>& signatures,
                                            const std::vector<std::pair<size_t,size_t>>& levels,
                                            const double s, const size_t bins, const double f, const size_t p, const size_t max_level){

    size_t max = std::numeric_limits<size_t>::max();
    const std::vector<std::vector<std::uint64_t>>& oph_sigs = std::get<0>(signatures);
    const std::vector<std::vector<std::uint64_t>>& fracmin_sigs = std::get<1>(signatures);

    std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>> seq_layout; // We already want to store information according to Layout standards.
    std::vector<std::vector<size_t>> max_bin_ids;
    std::vector<std::vector<std::pair<size_t,size_t>>> parents;

    lemon::ListGraph graph;
    std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>> labMaps = generate_all<Hasher>(oph_sigs, levels, graph);
    std::vector<std::unordered_map<size_t, const std::vector<size_t>*>> clusts = get_clusters(labMaps);

    auto refine_and_bin = [&](const std::vector<size_t>& seqs, size_t sub_bins, size_t lower, size_t upper, bool all_seqs) {
        std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> res;
        size_t curr_lower = lower;
        size_t curr_upper = upper;
        size_t curr_t_max = (curr_lower + curr_upper)/(2*sub_bins*s);
        size_t old_t_max = 0;

        for(size_t it = 0; it <= p; it){
            if(curr_t_max == old_t_max) break; //  reached convergence.
            res = all_seqs ? binning(labMaps, clusts, fracmin_sigs, s, sub_bins, curr_t_max, f) : binning_given_seqs(labMaps, clusts, fracmin_sigs, seqs, s, sub_bins, curr_t_max, f);

            bool overflow = std::get<3>(res);
            if(overflow){
                curr_lower = curr_t_max;
                if(curr_t_max == 0) break;
                old_t_max = curr_t_max;
                curr_t_max = (curr_lower + curr_upper)/(2*sub_bins*s);
                continue;
            }
            if(it == p) break; // Last iteration done

            const auto& [split_start, split_bins, merge_start] = std::get<1>(res);
            const std::vector<std::vector<size_t>>& result = std::get<0>(res);

            size_t split_bin_amt = merge_start - split_start;
            size_t merge_bin_amt = result.size() - merge_start;

            if(split_bin_amt == 0 && merge_bin_amt == 0) break; // can't refine  if the IBF is empty

            size_t split_avg = split_bin_amt ? splitting_average(result, fracmin_sigs, split_start, merge_start, s, f) : 0;
            size_t merge_avg = merge_bin_amt ? merge_average(result, fracmin_sigs, merge_start, s) : 0;

            if(split_avg > merge_avg) curr_upper = curr_t_max;
            else curr_lower = curr_t_max;

            if(curr_t_max == 0) break; // If t_max should get really small
            old_t_max = curr_t_max;
            curr_t_max = (curr_lower + curr_upper)/(2*sub_bins*s);
            it += 1;
        }
        return res;
    };

    auto record_bins = [&](const IBF& result, size_t split_start, size_t merge_start, size_t index){
        for(size_t b = split_start; b < merge_start;){
            if(result[b].empty()){b += 1; continue;}
            size_t seq = result[b][0];
            size_t start = b;
            size_t count = 0;
            while(b < merge_start && !result[b].empty() && result[b][0] == seq){count += 1; b += 1;}
            seq_layout[seq].push_back({index,start, count});
        }

        for(size_t b = merge_start; b < result.size(); b++){
            for(size_t seq : result[b]) seq_layout[seq].push_back({index,b,1});
        }
    };

    auto get_upper_lower = [&](const std::vector<size_t>& seqs, size_t sub_bins){
        std::vector<const std::vector<std::uint64_t>*> ptrs;
        ptrs.reserve(seqs.size());
        for(size_t seq : seqs) ptrs.push_back(&fracmin_sigs[seq]);

        size_t union_size = get_union_size_ptr(ptrs);
        size_t sum = 0;
        for(size_t seq : seqs) sum += fracmin_sigs[seq].size();

        return std::pair{union_size, sum};
    };

    auto get_sub_bins = [&](size_t N){
        if(N == 0) return size_t{0};
        if(N <= 64) return N; // No need to make a fuss out of this minute stuff.
        double raw = std::sqrt(static_cast<double>(N));
        size_t rounded = static_cast<size_t>(std::ceil(raw/64.0))*64; // round to nearest multiple of 64
        if (rounded == 0) rounded = 64;
        return std::min(rounded, static_cast<size_t>(2000));
    };

    size_t union_size = get_union_size(fracmin_sigs);
    size_t sum_size = 0;
    for(const std::vector<std::uint64_t>& sketch : fracmin_sigs) sum_size += sketch.size();

    std::vector<std::vector<IBF>> hibf_levels; // We save every single level here.
    std::vector<std::vector<std::tuple<size_t,size_t,size_t>>> ranges; // In order to reconstruct, we need to know the Merge ranges for every IBF

    std::vector<size_t> dummy =  {0};

    auto root = refine_and_bin(dummy, bins, union_size, sum_size, true); // Since refine_and_bin doesnt need a spefific vector when calling binning, this dummy will do the trick.
    hibf_levels.push_back({std::get<0>(root)});
    auto [split_start, split_bins, merge_start] = std::get<1>(root);
    ranges.push_back({std::get<1>(root)});
    record_bins(std::get<0>(root), split_start, merge_start, 0);
    std::vector<size_t> trackfill = std::get<2>(root);
    size_t max_bin_id = trackfill.empty() ? 0 : std::distance(trackfill.begin(), std::max_element(trackfill.begin(), trackfill.end()));
    max_bin_ids.push_back({max_bin_id});
    parents.push_back({{max, max}});

    for(size_t lvl = 0; lvl + 1 < max_level; lvl++){
        std::vector<IBF> next_lvl;
        std::vector<std::tuple<size_t,size_t,size_t>> next_ranges;
        std::vector<size_t> next_max_bin_ids;
        std::vector<std::pair<size_t,size_t>> next_parents;

        for(size_t ibf_index = 0; ibf_index < hibf_levels[lvl].size(); ibf_index++){
            const IBF& ibf = hibf_levels[lvl][ibf_index];
            auto [split_start, split_bins, merge_start] = ranges[lvl][ibf_index];

            for(size_t b = merge_start; b < ibf.size(); b++){
                if(ibf[b].empty()) continue; // can't do stuff on an empty IBF.

                const std::vector<size_t>& sub_seqs = ibf[b];
                size_t sub_bins = get_sub_bins(sub_seqs.size());
                const auto& [sub_lower, sub_higher] = get_upper_lower(sub_seqs, sub_bins);

                auto child = refine_and_bin(sub_seqs, sub_bins, sub_lower, sub_higher, false);
                auto [split_start, split_bins, merge_start] = std::get<1>(child);
                std::vector<size_t> trackfill = std::get<2>(child);
                next_lvl.push_back(std::move(std::get<0>(child)));
                max_bin_id = trackfill.empty() ? 0 : std::distance(trackfill.begin(), std::max_element(trackfill.begin(), trackfill.end()));
                next_max_bin_ids.push_back(max_bin_id);
                next_ranges.push_back(std::get<1>(child));
                next_parents.push_back({ibf_index, b});
                record_bins(next_lvl.back(), split_start, merge_start, next_lvl.size() -1);

            }
        }
        if(next_lvl.empty()) break;
        hibf_levels.push_back(std::move(next_lvl));
        ranges.push_back(std::move(next_ranges));
        max_bin_ids.push_back(std::move(next_max_bin_ids));
        parents.push_back(std::move(next_parents));
    }

    return {hibf_levels, ranges, seq_layout, max_bin_ids, parents};
}