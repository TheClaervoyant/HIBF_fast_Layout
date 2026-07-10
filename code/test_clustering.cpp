#include "test_clustering.h"
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>

std::vector<std::uint64_t> gen_random(size_t size){
    std::vector<std::uint64_t> res(size);
    static std::mt19937_64 mt64(100000000000);
    for(std::uint64_t& elem: res) elem = mt64();
    return res;
}

std::vector<std::vector<std::uint64_t>> get_randoms(std::vector<std::uint64_t> base, size_t amt, size_t size, double j_sim, std::uint64_t seed){
    std::vector<std::vector<std::uint64_t>> res;
    size_t same_frac = static_cast<size_t>(std::round(j_sim * 2.0 * size)/(1.0 + j_sim));
    size_t replace_frac = size - same_frac;
    for(size_t i = 0; i < amt; i++){ 
        std::mt19937_64 mt64(seed + i);
        std::vector<std::uint64_t> tmp = base;

        std::vector<size_t> positions(size);
        std::iota(positions.begin(), positions.end(), 0);

        std::shuffle(positions.begin(), positions.end(), mt64);
        for(size_t j = 0; j < replace_frac; j++){
            tmp[positions[j]] = mt64(); // Change only the given fraction again, the other (1 - dist) integers should stay the same across every vector
        }
        res.push_back(tmp);
    }
    return res;
}

std::vector<std::vector<std::uint64_t>> get_clusters(size_t clusters, size_t clust_size,
                                                     size_t size, double j_sim, double clust_j_sim){

    std::vector<std::uint64_t> global = gen_random(size);
    std::vector<std::vector<std::uint64_t>> res;
    res.reserve(clusters * clust_size);

    // Generate "Centromers" for each cluster
    std::vector<std::vector<std::uint64_t>> bases = get_randoms(global, clusters, size, clust_j_sim, 1);

    // Generate Clusters around every centromer
    for(std::uint64_t clust_num = 0; clust_num < clusters; clust_num++){
        std::vector<std::vector<std::uint64_t>> cluster = get_randoms(bases[clust_num], clust_size, size, j_sim, clust_num*clust_size + clusters + 1);

        res.insert(res.end(), cluster.begin(), cluster.end());
    }
    return res;
}


std::vector<std::vector<std::uint64_t>> expand(const std::vector<std::uint64_t>& base, size_t size, const std::vector<double>& j_sims, const std::vector<size_t> counts, size_t seed){
    if(j_sims.empty()) return {base}; // basecase

    std::vector<std::vector<std::uint64_t>> res;

    std::vector<std::vector<std::uint64_t>> children = get_randoms(base, counts[0], size, j_sims[0], seed);

    // Shift the vector one further, since the j_sim and counts here are done
    std::vector<double> remaining_j_sims(j_sims.begin() + 1, j_sims.end());
    std::vector<size_t> reamaining_counts(counts.begin() + 1, counts.end());

    for(size_t i = 0; i < children.size(); i++){
        std::vector<std::vector<std::uint64_t>> leaves = expand(children[i], size, remaining_j_sims, reamaining_counts, seed + (i+1)*1000); // Generate a new seed for every new batch
        res.insert(res.end(), leaves.begin(), leaves.end());
    }
    return res;
}

std::vector<std::vector<std::uint64_t>> get_any_cluster(size_t size, const std::vector<double>& j_sims, const std::vector<size_t>& counts, size_t seed){
    std::vector<std::uint64_t> global = gen_random(size);
    return expand(global, size, j_sims, counts, seed);
}