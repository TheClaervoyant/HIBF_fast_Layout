#pragma once

#include <cstdint>
#include <vector>

// @brief generate a vector of a given size with random integers.
// @param size : The size of the vector
std::vector<std::uint64_t> gen_random(size_t size);

// @brief generate a set of vectors of integers with a given Jaccard Similarity.
// @param base : The vector to base of the calculation
// @param amt : The number of vectors to be generated
// @param site : The size of every of the returned vectors.
// @param j_sim : The Jaccard Similarity every vector might have to one another
// @param seed : The seed used in the random number generator
std::vector<std::vector<std::uint64_t>> get_randoms(std::vector<std::uint64_t> base, size_t amt, size_t size, double j_sim, std::uint64_t seed = 0);

// @brief generate random clusters of integers with preset distances
// @param clusters : The amount of Clusters to be generated.
// @param clust_size : The amount of entries in each cluster.
// @param size : The size of each vector.
// @param dist : The similarity of every element in each cluster
// @param clust_dist : The Similarity between clusters
std::vector<std::vector<std::uint64_t>> get_clusters(size_t clusters, size_t clust_size,
                                                     size_t size, double j_sim, double clust_j_sim);


// @brief expand a level of the cluster.
// @param base : The "centromer" of the cluster
// @param size : The size of the vector
// @param j_sims : the similarity on each lvl.
// @param counts : The amount of clusters on each lvl
// @param seed : The seed to be used here
std::vector<std::vector<std::uint64_t>> expand(const std::vector<std::uint64_t>& base, size_t size, const std::vector<double>& j_sims, const std::vector<size_t> counts, size_t seed);

// @brief Generate Clusters with any depth and any amount of sequences.
// @param size : The size of each vector
// @param j_sims : The similarity on each lvl 
// @param counts : The amount of clusters on each level
// @param seed : The seed to be used by default.
std::vector<std::vector<std::uint64_t>> get_any_cluster(size_t size, const std::vector<double>& j_sims, const std::vector<size_t>& counts, size_t seed = 0);