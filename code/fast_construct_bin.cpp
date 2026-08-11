#include "fast_construct_bin.h"
#include "fast_construct_hashing.h"

std::vector<double> compute_fcorrs(double fpr, size_t h){
    std::vector<double> res(101, 1.0);
    const double log_num = std::log(1.0 - std::pow(fpr, 1.0 / static_cast<double>(h)));

    for(size_t s = 1; s <= 100; s++){
        double pcorr = 1.0 - std::pow(1.0 - fpr, 1.0 / static_cast<double>(s));
        double log_denom = std::log(1.0 - std::pow(pcorr, 1.0 / static_cast<double>(h)));
        res[s] = log_num / log_denom;
    }
    return res;
}

size_t merge_average(const std::vector<size_t>& track_fill, const size_t merge_start){
    size_t sum = 0;
    size_t amount = 0;
    for(size_t b = merge_start; b < track_fill.size(); b++){
        if(track_fill[b] == 0) continue;
        sum += track_fill[b];
        amount += 1;  
    }
    return amount ? sum / amount : 0;
}

size_t splitting_average(const std::vector<size_t>& track_fill, const size_t split_start, const size_t split_end){
    size_t sum = 0;
    size_t amount = 0;
    for(size_t b = split_start; b < split_end; b++){
        if(track_fill[b] == 0) continue;
        sum += track_fill[b];
        amount += 1;
    }
    return amount ? sum / amount : 0;
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