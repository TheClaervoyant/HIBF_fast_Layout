#include "fast_construct_graph.h"
#include "test_clustering.h"
#include "fast_construct_hashing.h"
#include <lemon/connectivity.h>
#include <fstream>
#include <iostream>

// Helper Function to parse the set of doubles
std::vector<double> parse_doubles(const std::string& s){
    std::vector<double> res;
    std::stringstream ss(s);
    std::string tok;
    while(std::getline(ss,tok,',')) res.push_back(std::stod(tok));
    return res;
}

// Helper Function to parse the set of size_ts
std::vector<size_t> parse_sizes(const std::string& s){
    std::vector<size_t> res;
    std::stringstream ss(s);
    std::string tok;
    while(std::getline(ss,tok,',')) res.push_back(std::stoul(tok));
    return res;
}

std::vector<std::pair<size_t,size_t>> parse_lvls(const std::string& s){
    std::vector<std::pair<size_t,size_t>> res;
    std::stringstream ss(s);
    std::string tok;
    while(std::getline(ss,tok,',')){
        auto colon = tok.find(':');
        size_t bcount = std::stoul(tok.substr(0,colon));
        size_t bwidth = std::stoul(tok.substr(colon + 1));
        res.push_back({bcount,bwidth});
    }
    return res;
}

std::string printout(const size_t seq){
    return "Sequenz " + std::to_string(seq);
}

std::string printout(const std::vector<size_t>& labels){
    std::string res = "Sequenzen ";
    for(size_t lab : labels) res += (std::to_string(lab) += ", ");
    return res;
}

void printgraph(lemon::ListGraph& graph, std::unordered_map<size_t, lemon::ListGraph::Node>& labMap, std::string filepath){
    // Get the connected Components in order to colorize them later on.
    lemon::ListGraph::NodeMap<int> compMap(graph);
    int numComp = lemon::connectedComponents(graph, compMap);

    // In Order to print, we invert the labMap
    std::map<lemon::ListGraph::Node, size_t> nodeMap;
    for(auto& [label, node] : labMap) nodeMap[node] = label;

    // Potential colors, cycle around if there are more than 10 con.comps.
    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
    };

    std::ofstream dot(filepath);
    dot << "graph G {\n" << "  node [style = filled];\n";
    for(lemon::ListGraph::NodeIt n(graph); n!= lemon::INVALID; ++n){ // Print every node and its label
        std::string color = colors[compMap[n] % colors.size()]; // Get the color, wrap around if more than colors contained.
        dot << "  " << graph.id(n) << " [label=\"" << printout(nodeMap[n]) << "\", color=\"" << color << "\"];\n";
    }
    for(lemon::ListGraph::EdgeIt e(graph); e!= lemon::INVALID; ++e){ // Print every edge and its nodes
        dot << "  " << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e)) << ";\n";
    }
    dot << "}\n";
    dot.close();
}

template <typename Hasher>
void printgraph(lemon::ListGraph& graph, std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>& labMap, std::string filepath){
    // Get the connected Components in order to colorize them later on.
    lemon::ListGraph::NodeMap<int> compMap(graph);
    int numComp = lemon::connectedComponents(graph, compMap);

    // In Order to print, we invert the labMap
    std::map<lemon::ListGraph::Node, std::vector<size_t>> nodeMap;
    for(auto& [label, node] : labMap) nodeMap[node] = label;

    // Potential colors, cycle around if there are more than 10 con.comps.
    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
    };

    std::ofstream dot(filepath);
    dot << "graph G {\n" << "  node [style = filled];\n";
    for(lemon::ListGraph::NodeIt n(graph); n!= lemon::INVALID; ++n){ // Print every node and its label
        std::string color = colors[compMap[n] % colors.size()]; // Get the color, wrap around if more than colors contained.
        dot << "  " << graph.id(n) << " [label=\"" << printout(nodeMap[n]) << "\", color=\"" << color << "\"];\n";
    }
    for(lemon::ListGraph::EdgeIt e(graph); e!= lemon::INVALID; ++e){ // Print every edge and its nodes
        dot << "  " << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e)) << ";\n";
    }
    dot << "}\n";
    dot.close();
}

template <typename Hasher>
void printgraph(lemon::ListGraph& graph, std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, std::string filepath){
    // Get the connected Components in order to colorize them later on.
    lemon::ListGraph::NodeMap<int> compMap(graph);
    int numComp = lemon::connectedComponents(graph, compMap);

    // In Order to print, we invert the labMap
    std::map<lemon::ListGraph::Node, std::vector<size_t>> nodeMap;
    for(auto& labMap : labMaps) // labMap : unordered_map<vector, Node, Hasher>.
        for(auto& [label, node] : labMap) nodeMap[node] = label;

    // Potential colors, cycle around if there are more than 10 con.comps.
    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
    };

    std::ofstream dot(filepath);
    dot << "graph G {\n" << "  node [style = filled];\n";
    for(lemon::ListGraph::NodeIt n(graph); n!= lemon::INVALID; ++n){ // Print every node and its label
        std::string color = colors[compMap[n] % colors.size()]; // Get the color, wrap around if more than colors contained.
        dot << "  " << graph.id(n) << " [label=\"" << printout(nodeMap[n]) << "\", color=\"" << color << "\"];\n";
    }
    for(lemon::ListGraph::EdgeIt e(graph); e!= lemon::INVALID; ++e){ // Print every edge and its nodes
        dot << "  " << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e)) << ";\n";
    }
    dot << "}\n";
    dot.close();
}

int main(int argc, char* argv[]){
    if(argc != 5){
        std::cerr << "Missing parameters: " << "Amount cluster + cluster size + vector_size + Similarity in cluster + Similarity berween clusters. \n";
        return 1;
    }

    size_t vec_size = std::stoul(argv[1]);
    std::vector<double> j_sims = parse_doubles(argv[2]);
    std::vector<size_t> counts = parse_sizes(argv[3]);
    std::vector<std::pair<size_t, size_t>> lvls = parse_lvls(argv[4]);

    std::vector<std::vector<std::uint64_t>> rand_clusts = get_any_cluster(vec_size, j_sims, counts, 0);

    lemon::ListGraph graph;
    
    // Generate One Permutation Hashes for each "sequence":
    std::vector<std::vector<std::uint64_t>> oph_sigs = one_permutation_hash(rand_clusts, 8);
    
    std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, standardHasher>> labMaps = generate_all<standardHasher>(oph_sigs, lvls, graph);

    printgraph(graph, labMaps, "../results/All_lvls_analyzed.dot");

    std::cout <<"Each level in one go can be seen in ../results/All_lvls_analyzed.dot. \n";
    return 0;
}