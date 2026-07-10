#include "../fast_construct_graph.h"
#include <numeric>

// We use this function as a special case wrapper.
template <typename Hasher>
void connect_bands_all_seqs(const std::vector<std::vector<std::uint64_t>>& hashes,
                   const size_t bcount, const size_t bwidth,
                   lemon::ListGraph& graph,
                   std::unordered_map<size_t,lemon::ListGraph::Node>& LabelMap,
                   Hasher hasher){
  std::vector<size_t> all_seqs(hashes.size());
  std::iota(all_seqs.begin(),all_seqs.end(),0); // This generates a vector of size_t with 0...n-1. Whichh is what we want, since EVERY sequence is important on the first level.

  connect_bands(hashes, all_seqs, bcount, bwidth, graph, LabelMap, hasher);
}