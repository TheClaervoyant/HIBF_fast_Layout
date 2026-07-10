#include "fast_construct_graph.h"
#include <lemon/core.h>
#include <lemon/connectivity.h>

void construct_graph_tower(const size_t seqs_, const size_t bands_,
               lemon::ListGraph& graph, 
               std::map<std::tuple<size_t,size_t>,lemon::ListGraph::Node>& LabelMap){
  for(size_t seq = 0; seq < seqs_; seq++){
  lemon::ListGraph::Node precursor = lemon::INVALID; // We use this in order to simply be able to connect a node with its precursor.
  for(size_t band = 0; band < bands_; band++){
    lemon::ListGraph::Node nd = graph.addNode();
    LabelMap[{seq,band}] = nd;
    if (band != 0){ // This is not the first Node of a sequence, so it is entirely possible to connect it with its precursor.
      graph.addEdge(precursor,nd);
    }
    precursor = nd; // We itterate of the next node; our node now will be the new precursor.
    }
  }  
} 

// Just like construct_graph_tower, but no connecting needed. And since we dont save (Seq,band), the map only takes the size_t instead of the tuple.
void construct_graph(const size_t seqs_,
               lemon::ListGraph& graph,
               std::unordered_map<size_t,lemon::ListGraph::Node>& LabelMap){
  for(size_t seq = 0; seq < seqs_; seq++){
    lemon::ListGraph::Node sequence = graph.addNode();
    LabelMap[seq] = sequence;
  }
}

void construct_graph(const std::vector<size_t>& labs_, lemon::ListGraph& graph, std::unordered_map<size_t, lemon::ListGraph::Node>& LabelMap){
  for(size_t lab : labs_){
    lemon::ListGraph::Node node = graph.addNode();
    LabelMap[lab] = node;
  }
}

void connect_bands_tower(const std::vector<std::vector<std::uint64_t>>& hashes,
                   const size_t bcount, const size_t bwidth,
                   lemon::ListGraph& graph,
                   std::map<std::tuple<size_t,size_t>,lemon::ListGraph::Node>& LabelMap,
                   std::function<std::uint64_t(const std::vector<std::uint64_t>&)> hashFunc){
    
  // Create a Wrapper Struct ad hoc inside the function, so that it can be used in the HashMap
  struct FuncWrapper {
    std::function<std::uint64_t(const std::vector<std::uint64_t>&)> func;
    std::uint64_t operator()(const std::vector<std::uint64_t>& vec) const { return func(vec); }
  };

  // @idea Connect Nodes {a,b} with {x,y} if "b == y", as in they collide in a hashmap.
  std::unordered_map<std::vector<std::uint64_t>,std::vector<std::tuple<size_t,size_t>>,FuncWrapper> buckets(0, FuncWrapper{hashFunc});

  // Insert every Node into the vector. In the end, itterate over every vector and connect the Nodes.
  for(size_t seq = 0; seq < hashes.size() ; seq++){
    for(size_t band = 0; band < bcount; band++ ){
      // Extract the band needed
      auto start = hashes[seq].begin() + band * bwidth;
      auto end = hashes[seq].begin() + (band +1) * bwidth;
      std::vector<std::uint64_t> curr_band(start,end);

      auto& bucket = buckets[curr_band];
      // Possible Test but probably not needed if(std::find(bucket.begin(), bucket.end(), std::make_tuple(seq,band)) == bucket.end())
      bucket.push_back({seq,band});
    }
  }
  for(const auto& [hashed, labels] : buckets){ // Every bucket
    for(size_t i = 1; i < labels.size() ; i++){ // Due to transitivity it is completely fine to connect every Node with Node 0; if Node 0 == Node i and Node 0 == Node j then it follows that Node i == Node j.
      if(std::get<0>(labels[0]) != std::get<0>(labels[i])) graph.addEdge(LabelMap[labels[0]],LabelMap[labels[i]]); //Only finally add the edge if they're not part of the same sequence, since this still could cause multi edges.
    }
  }
}

std::vector<std::vector<size_t>> connected_components(lemon::ListGraph& graph, const std::unordered_map<size_t, lemon::ListGraph::Node>& labMap){
  lemon::ListGraph::NodeMap<int> compMap(graph);
  int numComp = lemon::connectedComponents(graph, compMap);

  // Invert the Map. We need to Map every node onto its label.
  std::map<lemon::ListGraph::Node, size_t> nodeMap;
  for(auto& [label, node] : labMap) nodeMap[node] = label;

  std::vector<std::vector<size_t>> components(numComp);
  for(lemon::ListGraph::NodeIt n(graph); n != lemon::INVALID; ++n) components[compMap[n]].push_back(nodeMap[n]);

  return components;
}
