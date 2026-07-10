#pragma once

#include <vector>
#include <cstdint>
#include <functional>
#include <filesystem>
#define XXH_INLINE_ALL
#include "xxhash.h"

// @brief a wrapper around the uint64 hash of xxHash.
// @param input : the integer to hash
std::uint64_t xxhash_wrap(std::uint64_t input);

// @brief Given a set of hashes, compute the one permutation hash of the given collection
// @param hashes : the set of hashes to be hashed
// @param k : the intervall the function is cut in.
// @param hashFunc : the Hash Function to be used
//
// @note in the end, returns a vector of size 2^k.
// @note this function is rather for debugging purposes, when actually using stuff, refer to @func ophs
std::vector<std::vector<std::uint64_t>> one_permutation_hash(const std::vector<std::vector<std::uint64_t>>& hashes,
                                                const std::uint8_t k,
                                                std::function<std::uint64_t(std::uint64_t)> hashFunc = xxhash_wrap);

// @brief Given a set of sequences, k and a hashfunction, compute the one permutation hash of the given sequences.
// @param filepath : Path to the file containing the sequences to be hashed
// @param q : the size of the Q-Grams applied to the given #include <cstdint>sequences
// @param k : the Intervall the function is cut in.
// @param hashFunc : the Hash Function to be used.
//
// @note every vector inside of the returned one has a size of 2^k.
std::vector<std::vector<std::uint64_t>> ophs(const std::filesystem::path& filepath,
                                              const std::uint8_t q, const std::uint8_t k,
                                              std::function<std::uint64_t(std::uint64_t)> hashFunc = xxhash_wrap);


// @brief Given a set of One-Permutation Hashes, extract the wanted number of bands, given the bandwidth
// @param hashes : The set of One-Permutation Hashes.
// @param bwidth : The size of one band.
// @param bcount : The number of bands to be extracted.
//
// @note for any i : hashes[i].size() > bwidth * bcount.
// @note EveryHash(CollectionOfHashes(SingularBand)) is the nesting
std::vector<std::vector<std::vector<std::uint64_t>>> extract_bands(const std::vector<std::vector<std::uint64_t>>& hashes,
                                                                   const size_t bwidth, const size_t bcount);

// @brief Hash the given Vector. Using XXHash (https://github.com/cyan4973/xxhash)
// @param vec the vector to be hashed
template <typename T>
std::uint64_t hash_vector(const std::vector<T>& vec){
    return XXH64(vec.data(), vec.size() * sizeof(T), 0);
}

// @brief Hash the given vector by returning its size.
// @param vec the vector to be hashed.
// @note this function is debugging purposes and for little examples only!!!
template <typename T>
std::uint64_t hash_vector_debug(const std::vector<T>& vec){
  std::uint64_t res = 0;
  for(std::uint64_t elem: vec){
    res += elem;
  }
  return res;
}