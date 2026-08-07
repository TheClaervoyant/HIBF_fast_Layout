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