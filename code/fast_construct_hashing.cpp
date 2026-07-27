#include "fast_construct_hashing.h"
#include <limits>
#include <cmath>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>


std::uint64_t xxhash_wrap(std::uint64_t input){ return XXH64(&input, sizeof(input), 0); }

std::pair<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>> one_permutation_fracmin_hash(const std::vector<std::vector<std::uint64_t>>& hashes, const std::uint8_t k, const double s, std::function<std::uint64_t(std::uint64_t)> hashFunc){
  std::vector<std::vector<std::uint64_t>> res_oph;
  res_oph.reserve(hashes.size());
  std::vector<std::vector<std::uint64_t>> res_fracmin;
  res_fracmin.reserve(hashes.size());

  for(const std::vector<std::uint64_t>& hashes_ : hashes){
    std::uint64_t max = std::numeric_limits<std::uint64_t>::max();
    std::vector<std::uint64_t> tmp_oph(std::pow(2,k), max);
    std::vector<std::uint64_t> tmp_fracmin;
    std::uint64_t thresh = static_cast<std::uint64_t>(static_cast<double>(max*s));
    for(std::uint64_t hash : hashes_){
      std::uint64_t hash_ = hashFunc(hash);

      // Fracmin Part
      if(hash_ < thresh) tmp_fracmin.push_back(hash_);

      // OPH Part
      std::uint64_t index = hash_ & ((1ULL << k) -1); //Extract the first k bits of the initial hash value.
      std::uint64_t val = hash_ >> k; // remove the k bits, they'd just be a bias since evry value in the bucket would ALWAYS have these bits same.
      if(val < tmp_oph[index]) tmp_oph[index] = val;
    }
    std::sort(tmp_fracmin.begin(), tmp_fracmin.end());

    res_oph.push_back(tmp_oph);
    res_fracmin.push_back(tmp_fracmin);
  }
  return {std::move(res_oph), std::move(res_fracmin)};
}

std::pair<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>> ophs_fmhs(const std::filesystem::path& filepath, const std::uint8_t q, const std::uint8_t k, const double s, std::function<std::uint64_t(std::uint64_t)> hashFunc){
  std::vector<std::vector<std::uint64_t>> res_oph;
  std::vector<std::vector<std::uint64_t>> res_fmh;

  auto fin = seqan3::sequence_file_input{filepath};
  auto qmer_view = seqan3::views::kmer_hash(seqan3::ungapped{q});
  
  for(auto & record : fin){
    std::uint64_t max = std::numeric_limits<std::uint64_t>::max();
    std::vector<std::uint64_t> tmp_oph(std::pow(2,k), max);
    std::vector<std::uint64_t> tmp_fracmin;
    std::uint64_t thresh = static_cast<std::uint64_t>(static_cast<double>(max*s));

    for(auto&& qgram : record.sequence() | qmer_view){
      std::uint64_t hash = hashFunc(qgram);

      if(hash < thresh) tmp_fracmin.push_back(hash);

      std::uint8_t index = hash & ((1ULL << k) -1);
      std::uint64_t val = hash >> k; // remove the k bits, they'd just be a bias since every value in the bucket would ALWAYS have these bits same.

      if(val < tmp_oph[index]) tmp_oph[index] = val;
    }
    std::sort(tmp_fracmin.begin(), tmp_fracmin.end());

    res_fmh.push_back(tmp_fracmin);
    res_oph.push_back(tmp_oph);  
  }
  return {std::move(res_oph), std::move(res_fmh)};
}

size_t get_union_size(const std::vector<std::vector<std::uint64_t>>& sketches){
  std::unordered_set<std::uint64_t> elems;
  size_t tots = 0; // Used for reserving space
  for(const std::vector<std::uint64_t>& sketch : sketches) tots += sketch.size();
  elems.reserve(tots); // "Worst Case" every Element is unique.

  for(const std::vector<std::uint64_t>& sketch : sketches) for(std::uint64_t elem : sketch) elems.insert(elem);

  return elems.size();
}

// @note this function is kind of redacted, sind its entire functionality is now in connect_bands. With the advantage, that we dont have to store every band respective. 
//  It is still used in the tests as the loop stays the same and is used in the connect_bands function.
std::vector<std::vector<std::vector<std::uint64_t>>> extract_bands(const std::vector<std::vector<std::uint64_t>>& hashes, const size_t bwidth, const size_t bcount){
  std::vector<std::vector<std::vector<std::uint64_t>>> res;
  for(std::vector<std::uint64_t> oph : hashes){ // itterate over every given hash
    std::vector<std::vector<std::uint64_t>> tmp;
    for(size_t count = 0; count < bcount; count++){ // We extract the Bands from the given signature
      auto start = oph.begin() + count * bwidth;
      auto end = oph.begin() + (count + 1) * bwidth;
      tmp.push_back(std::vector<std::uint64_t>(start,end)); // Work with itterators - efficiency
    }
    res.push_back(tmp);
  }
  return res;
}
