#include "../fast_construct_graph.h"

template <typename Hasher>
void connect_bands(const std::vector<std::vector<std::uint64_t>>& hashes,
                   const std::vector<size_t>& sig_labs,
                   const size_t bcount, const size_t bwidth,
                   lemon::ListGraph& graph,
                   std::unordered_map<size_t,lemon::ListGraph::Node>& LabelMap,
                   Hasher hasher){

  // @idea connect Node a with node b iff there Exist i,j st. sequence a band i == sequence b band j
  std::unordered_map<std::vector<std::uint64_t>,std::vector<lemon::ListGraph::Node>,Hasher> buckets;
  for(size_t lab : sig_labs){
        for(size_t band = 0; band < bcount; band++ ){
      auto start = hashes[lab].begin() + band * bwidth;
      auto end = hashes[lab].begin() + (band + 1) * bwidth;
      std::vector<std::uint64_t> curr_band(start,end);
      auto& bucket = buckets[curr_band];
      // Possible Test but probably not needed: if(std::find(bucket.begin(), bucket.end(), LabelMap[seq]) == bucket.end())
      bucket.push_back(LabelMap[lab]);
    }
  }
  for(const auto& [hashed,nodes] : buckets){ // Every bucket
    for(size_t i = 1; i < nodes.size() ; i++){ // Due to transitivity it is completely fine to connect every Node with Node 0; if Node 0 == Node i and Node 0 == Node j then it follows that Node i == Node j.
      if(nodes[0] != nodes[i]) graph.addEdge(nodes[0],nodes[i]);
    }
  }
}