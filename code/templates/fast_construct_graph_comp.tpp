#include "../fast_construct_graph.h"

template <typename Hasher>
void construct_graph(const std::vector<std::vector<size_t>>& labs_, lemon::ListGraph& graph, std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>& LabelMap){
  for(const std::vector<size_t>& labs : labs_){
    lemon::ListGraph::Node node = graph.addNode();
    LabelMap[labs] = node;
  }
}