#include "fast_construct_bin.h"
#include "fast_construct_hashing.h"

size_t merge_average(const std::vector<std::vector<size_t>>& res, const std::vector<std::vector<std::uint64_t>>& fracmin_sketches, const size_t merge_start, const double s){
    size_t union_sum = 0;
    size_t amount = 0;
    for(size_t b = merge_start; b <res.size(); b++){
        if(res[b].empty()) continue; // Should by definition not happen but just in case.
        std::vector<const std::vector<std::uint64_t>*> ptrs;
        ptrs.reserve(res[b].size());
        for(size_t seq : res[b]) ptrs.push_back(&fracmin_sketches[seq]);
        union_sum += get_union_size_ptr(ptrs);
        amount += 1;
    }
    return amount ? static_cast<size_t>(static_cast<double>(union_sum)/(s*amount)) : 0;
}

size_t splitting_average(const std::vector<std::vector<size_t>>& res, const std::vector<std::vector<std::uint64_t>>& fracmin_sketches, const size_t split_start, const size_t split_end, const double s, const double f){
    std::unordered_map<size_t,size_t> occurences; // We need to know how many times each sequence occurs.
    size_t amount = 0;
    for(size_t b = split_start; b < split_end; b++){
        if(res[b].empty()) continue;
        amount += 1;
        occurences[res[b][0]] += 1; // We can assume that res[b] only contains one element by definition.
    }
    double split_sum = 0.0;
    for(auto& [seq, c] : occurences) split_sum += static_cast<double>(1.5*fracmin_sketches[seq].size()) / c;
    return amount ? static_cast<size_t>(split_sum/(s*amount)) : 0;
}

void write_header(std::ostream& out, const std::vector<std::vector<IBF>>& hibf_levels, const std::vector<std::vector<size_t>>& max_bin_ids, const std::vector<std::vector<std::pair<size_t,size_t>>>& parents){
    auto merge_bin_label = [](size_t level, size_t index, const std::vector<std::vector<std::pair<size_t,size_t>>>& parents){
        std::vector<size_t> chain;
        while(level > 0){
            auto [par_idx, merge_bin] = parents[level][index];
            chain.push_back(merge_bin);
            index = par_idx;
            level -= 1;
        }
        std::reverse(chain.begin(), chain.end());
        std::string lab;
        for(size_t i = 0; i < chain.size(); i++){
            if(i > 0) lab += ";";
            lab += std::to_string(chain[i]);
        }
        return lab;
    };

    for(size_t level = 0; level < hibf_levels.size(); level++){
        for(size_t ibf_idx = 0; ibf_idx < hibf_levels[level].size(); ibf_idx++){
            size_t const max_bin_id = max_bin_ids[level][ibf_idx];

            if(level == 0) out << "#HIGH_LEVEL_IBF max_bin_id:" << max_bin_id << "\n";
            else{
                std::string const label = merge_bin_label(level, ibf_idx, parents);
                out << "#MERGED_BIN_" << label << " max_bin_id:" << max_bin_id << "\n";
            }
        }
    }

    out << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n";
}

void write_content(std::ofstream& out, const std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>>& seq_layout, const std::unordered_map<size_t, std::string>& seq_to_path){
    for(auto const& [seq, entries] : seq_layout){
        auto path_it = seq_to_path.find(seq);
        if(path_it == seq_to_path.end()) continue; // Maybe throw? Shouldn't happen, though.

        std::string bin_indices;
        std::string number_of_bins;

        for(size_t i = 0; i < entries.size(); i++){
            auto const& [index, start, count] = entries[i];
            if(i > 0) {bin_indices += ";"; number_of_bins += ";";}
            bin_indices += std::to_string(start);
            number_of_bins += std::to_string(count);
        }
        out << path_it->second << "\t" << bin_indices << "\t" << number_of_bins << "\n";
    }
}