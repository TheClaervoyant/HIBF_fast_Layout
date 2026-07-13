#include "fast_construct_graph.h"
#include "fast_construct_hashing.h"
#include <iostream>
#include <lemon/core.h>
#include <lemon/connectivity.h>

std::string colored_bool(bool val){
  return val ? "\033[32mTRUE\033[0m" : "\033[31mFALSE\033[0m";
}

void test_construct_graph_tower(){
  // Overall setup
  lemon::ListGraph graph;
  std::map<std::tuple<size_t,size_t>, lemon::ListGraph::Node> labMap;
  construct_graph_tower(2, 3, graph, labMap);

  // @test since we have 2 sequences with 3 bands, check whether the labMap has 6 entries only.
  bool correct_size = (labMap.size() == 6);

  // @test we want an edge between every band in the same sequence but just transitive (~ a line). Test whether this is correct.
  bool correct_edges = true;
  for(size_t seq = 0; seq < 2; seq++){ // Everything should hold for both lines
    for(size_t band = 0; band < 2; band++){
      lemon::ListGraph::Node curr = labMap[{seq,band}];
      lemon::ListGraph::Node suc = labMap[{seq,band + 1}];
      bool connected = false;
      for(lemon::ListGraph::IncEdgeIt e(graph,curr); e != lemon::INVALID; ++e){ // itterate over every given edge of a given node and check whether the successor is connected -> the line
        if(graph.u(e) == suc || graph.v(e) == suc) {
          connected = true;
          break;
        }
      }
      if(!connected) correct_edges = false; // if there exists a pair that is NOT connected, set correct_edges with false
    }
  }

  // @test Between two Lines (= two sequences) there should not be any form of connection (an edge). If there exists one, the test fails.
  bool seperate_seqs = true;
  for(size_t band = 0; band < 3; band++){
    lemon::ListGraph::Node n0 = labMap[{0,band}];
    lemon::ListGraph::Node n1_0 = labMap[{1,0}];
    lemon::ListGraph::Node n1_1 = labMap[{1,1}];
    lemon::ListGraph::Node n1_2 = labMap[{1,2}];
    bool connected = false;
    for(lemon::ListGraph::IncEdgeIt e(graph,n0); e != lemon::INVALID; ++e){
      if(graph.u(e) == n1_0 || graph.v(e) == n1_0 ||
         graph.u(e) == n1_1 || graph.v(e) == n1_1 ||
         graph.u(e) == n1_2 || graph.v(e) == n1_2){
          connected = true;
          break;
         }
    }
    if(connected) seperate_seqs = false;
  }
  
  std::cout << "\n=== Test the construct_graph_tower function ===\n";
  std::cout << "Right amount of vertices: " << colored_bool(correct_size) << "\n";
  std::cout << "Sequences are a path: " << colored_bool(correct_edges) << "\n";
  std::cout << "Sequences are isolated: " << colored_bool(seperate_seqs) << "\n";
}

void test_construct_graph(){
  // Overall setup
  lemon::ListGraph graph;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap;
  construct_graph(2, graph, labMap);

  // @test since we have 2 sequences, check whether the labMap has 2 entries only.
  bool correct_size = (labMap.size() == 2);

  // @test Between two Nodes (= two sequences) there should not be any form of connection (an edge). If there exists one, the test fails.
  bool seperate_seqs = true;
  lemon::ListGraph::Node n0 = labMap[0];
  lemon::ListGraph::Node n1 = labMap[1];
  bool connected = false;
  for(lemon::ListGraph::IncEdgeIt e(graph,n0); e != lemon::INVALID; ++e){
  if(graph.u(e) == n1 || graph.v(e) == n1){
      connected = true;
      break;
      }
  }
  if(connected) seperate_seqs = false;
  
  std::cout << "\n=== Test the construct_graph_OneNode function ===\n";
  std::cout << "Right amount of vertices: " << colored_bool(correct_size) << "\n";
  std::cout << "Sequences are isolated: " << colored_bool(seperate_seqs) << "\n";
}

void test_extract_bands(){
  // Assume we only have one signature
  std::vector<std::vector<std::uint64_t>> test_signature = {{1,2,2,4,5,3,7,8,9}};

  // @test assuming we take a bandwidth of 0 and an arbirtrary bandcount. We would assume for every band to be empty.
  std::vector<std::vector<std::vector<std::uint64_t>>> test_1 = extract_bands(test_signature, 0, 3);
  std::vector<std::vector<std::vector<std::uint64_t>>> expected_1 = {{{},{},{}}};
  bool empty_bands = (test_1 == expected_1);

  // @test assuming we take a bandwidth of 9 in this case; we would assume, that test_2[0][0] equals test_signature[0] since one band is the whole signature here.
  std::vector<std::vector<std::vector<std::uint64_t>>> test_2 = extract_bands(test_signature, 9, 1);
  bool whole_band = (test_2[0][0] == test_signature[0]);

  // @test When both side cases work, one might test also a simple example.
  std::vector<std::vector<std::vector<std::uint64_t>>> test_3 = extract_bands(test_signature, 2, 3);
  std::vector<std::vector<std::vector<std::uint64_t>>> expected_3 = {{{1,2},{2,4},{5,3}}};
  bool example = (test_3 == expected_3);

  std::cout << "\n=== Test the extract_bands function ===\n";
  std::cout << "Bands empty when bwidth =  0: " << colored_bool(empty_bands) << "\n";
  std::cout << "Band is whole sequence when width = hash.size(): " << colored_bool(whole_band) << "\n";
  std::cout << "Test example case with width = 2, count = 3: " << colored_bool(example) << "\n";
}

void test_connect_bands_tower(){
  // @test we use the debug version here. That means, every vector gets hashed onto its size. If we generate hashes with different sums, the number of edges in the graph should stay the same.
  lemon::ListGraph graph_1;
  std::map<std::tuple<size_t,size_t>, lemon::ListGraph::Node> labMap_1;
  construct_graph_tower(2, 2, graph_1, labMap_1); //correct, tested beforehand.
  int edges_before = lemon::countEdges(graph_1);

  std::vector<std::vector<std::uint64_t>> hashes_1 = { {1,1,  2,2}, {3,3,   4,4} }; // since every band has a different sum, no collision should occur.
  connect_bands_tower(hashes_1, 2, 2, graph_1, labMap_1, debugHasher{});
  int edges_after = lemon::countEdges(graph_1);

  bool no_new_edges = (edges_before == edges_after); // We can simply compare before and after, since connect bands does not remove any edges, so a same number must mean those are the same edges.

  // @test We now want to test, if two bands get connected, when they hash onto the same value (= the vectors have the same sum)
  lemon::ListGraph graph_2;
  std::map<std::tuple<size_t,size_t>, lemon::ListGraph::Node> labMap_2;
  construct_graph_tower(2, 2, graph_2, labMap_2);
  edges_before = lemon::countEdges(graph_2);

  std::vector<std::vector<std::uint64_t>> hashes_2 = { {1,1,  2,2}, {1,1,  3,3} }; // The first band of both sequences have the same hash value.
  connect_bands_tower(hashes_2, 2, 2, graph_2, labMap_2, debugHasher{});
  // We'd expect now one edge more. And this edge should be between (0,0) and (1,0). So we test for both respectively.
  edges_after = lemon::countEdges(graph_2);
  bool one_more = (edges_after == (edges_before +1));

  lemon::ListGraph::Node seq0_band0 = labMap_2[{0,0}];
  lemon::ListGraph::Node seq1_band0 = labMap_2[{1,0}]; // Check whether both of them are connected
  bool connected = false;
  for(lemon::ListGraph::IncEdgeIt e(graph_2,seq0_band0); e != lemon::INVALID; ++e){
    if(graph_2.u(e) == seq1_band0 || graph_2.v(e) == seq1_band0) {
      connected = true;
      break;
    }
  }
  bool test_pass = (connected && one_more);

  std::cout << "\n=== Test the connect_bands_tower function ===\n";
  std::cout << "No new Edge when no collision: " << colored_bool(no_new_edges) << "\n";
  std::cout << "Added one new Edge when there is exactly one Collision: " << colored_bool(test_pass) << "\n";
} 

void test_connect_bands(){
  // @test we use the debug version here. That means, every vector gets hashed onto its sum. If we generate bands with different sums, the number of edges in the graph should stay the same.
  lemon::ListGraph graph_1;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap_1;
  construct_graph(2, graph_1, labMap_1); //correct, tested beforehand.
  int edges_before = lemon::countEdges(graph_1);

  std::vector<std::vector<std::uint64_t>> hashes_1 = { {1,1,   2,2}, {3,3,   4,4} }; // since every vector has a different sum, no collision should occur.
  connect_bands_all_seqs(hashes_1, 2, 2, graph_1, labMap_1, debugHasher{});
  int edges_after = lemon::countEdges(graph_1);

  bool no_new_edges = (edges_before == edges_after); // We can simply compare before and after, since connect bands does not remove any edges, so a same number must mean those are the same edges.

  // @test We now want to test, if two bands get connected, when they hash onto the same value (= the vectors have the same sum)
  lemon::ListGraph graph_2;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap_2;
  construct_graph(2, graph_2, labMap_2);
  edges_before = lemon::countEdges(graph_2);

  std::vector<std::vector<std::uint64_t>> hashes_2 = { {1,1,   2,2}, {1,1,    3,3} }; // The first band of both sequences have the same hash value.
  connect_bands_all_seqs(hashes_2, 2, 2, graph_2, labMap_2, debugHasher{});
  // We'd expect now one edge more. And this edge should be between (0,0) and (1,0). So we test for both respectively.
  edges_after = lemon::countEdges(graph_2);
  bool one_more = (edges_after == (edges_before +1));

  lemon::ListGraph::Node seq0 = labMap_2[0];
  lemon::ListGraph::Node seq1 = labMap_2[1]; // Check whether both of them are connected
  bool connected = false;
  for(lemon::ListGraph::IncEdgeIt e(graph_2,seq0); e != lemon::INVALID; ++e){
    if(graph_2.u(e) == seq1 || graph_2.v(e) == seq1) {
      connected = true;
      break;
    }
  }
  bool test_pass = (connected && one_more);

  // @test when two bands of the same sequence hash together, we do not want an Edge on the sequence itself. So when this happens, the number of nodes should stay the same.
  lemon::ListGraph graph_3;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap_3;
  construct_graph(2, graph_3, labMap_3); //correct, tested beforehand.
  edges_before = lemon::countEdges(graph_3);

  std::vector<std::vector<std::uint64_t>> hashes_3 = { {1,1,  1,1}, {2,2,   3,3} }; // band 0 and 1 of sequence one collide, we do not want an Edge here, though. So the number of Edges should stay the same.
  connect_bands_all_seqs(hashes_3, 2, 2, graph_3, labMap_3, debugHasher{});
  edges_after = lemon::countEdges(graph_3);

  bool no_edge_with_collsion = (edges_before == edges_after); // We can simply compare before and after, since connect bands does not remove any edges, so a same number must mean those are the same edges.

  std::cout << "\n=== Test the connect_bands_all_seqs function ===\n" << "==== Since this uses a wrapper for given seqs, it tests both ==== \n";
  std::cout << "No new Edge when no collision: " << colored_bool(no_new_edges) << "\n";
  std::cout << "Added one new Edge when there is exactly one Collision: " << colored_bool(test_pass) << "\n";
  std::cout << "There is no new Edge on the sequence itself, when two bands of the same sequence collide: " << colored_bool(no_edge_with_collsion) << "\n";
}

void test_connected_components(){
  // @test when there is no collision between bands, every connected Component should contain exactly one item (thus having a size of one)
  // @test also, both entries should be different as one is the first Node and the other the second one
  lemon::ListGraph graph_1;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap_1;
  construct_graph(2, graph_1, labMap_1);

  std::vector<std::vector<std::uint64_t>> hashes_1 = { {1,1,   2,2},   {3,3,    4,4}  };  // No Bands collide, there should be no edge drawn, meaning every Node is a singleton
  connect_bands_all_seqs(hashes_1, 2, 2, graph_1, labMap_1, debugHasher{});
  std::vector<std::vector<size_t>> con_comps = connected_components(graph_1, labMap_1);

  bool different_entries = (con_comps[0] != con_comps[1]); 
  bool only_singletons = (con_comps[0].size() == 1 && con_comps[1].size() == 1);

  // @test the same results as above should apply when Nodes of the same sequence collide
  lemon::ListGraph graph_2;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap_2;
  construct_graph(2, graph_2, labMap_2);

  std::vector<std::vector<std::uint64_t>> hashes_2 = { {1,1,   1,1},   {3,3,    4,4}  };  // No Bands collide, there should be no edge drawn, meaning every Node is a singleton
  connect_bands_all_seqs(hashes_2, 2, 2, graph_2, labMap_2, debugHasher{});
  con_comps = connected_components(graph_2, labMap_2);

  bool different_entries_own_col = (con_comps[0] != con_comps[1]); 
  bool only_singletons_own_col = (con_comps[0].size() == 1 && con_comps[1].size() == 1);
  
  // @test when two bands of different sequences collide, in this case the con_comps should have a size of one (there is only one) and this connected component has a size of 2
  lemon::ListGraph graph_3;
  std::unordered_map<size_t, lemon::ListGraph::Node> labMap_3;
  construct_graph(2, graph_3, labMap_3);

  std::vector<std::vector<std::uint64_t>> hashes_3 = { {1,1,   2,2},   {3,3,    1,1}  };  // No Bands collide, there should be no edge drawn, meaning every Node is a singleton
  connect_bands_all_seqs(hashes_3, 2, 2, graph_3, labMap_3, debugHasher{});
  con_comps = connected_components(graph_3, labMap_3);

  bool only_one_comp = (con_comps.size() == 1);
  bool both_nodes = (con_comps[0].size() == 2);

  std::cout << "\n=== Test the connected_components function ===\n";
  std::cout << "The two connected components are different: " << colored_bool(different_entries) << "\n";
  std::cout << "The two connected components contain one of the Nodes: " << colored_bool(only_singletons) << "\n";
  std::cout << "The connected components are different with own band collision: " << colored_bool(different_entries_own_col) << "\n";
  std::cout << "The conneczed components contain one of the nodes with own band collision: " << colored_bool(only_singletons_own_col) << "\n";
  std::cout << "When coliding the two sequences, there is only one component: " << colored_bool(only_one_comp) << "\n";
  std::cout << "The one connected component contains the only two nodes present: " << colored_bool(both_nodes) << "\n";
}

void test_recursive_step(){
  lemon::ListGraph graph1;
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> labMaps1(1);
  std::vector<std::vector<size_t>> comp = {{1}};
  std::vector<size_t> singleton_comp = {1};
  construct_graph(comp, graph1, labMaps1[0]);

  // @test we would await, since this is a singleton, that there will be no new nodes or edges. That means, before and after should stay the same.

  int nodes_before = lemon::countNodes(graph1);
  int edges_before = lemon::countEdges(graph1);

  std::vector<std::vector<std::uint64_t>> two_sigs = {{0,0,0,0}, {1,1,1,1}};
  std::vector<std::pair<size_t,size_t>> lvls1 = {{2,2}};
  recursive_step(two_sigs, graph1, labMaps1, singleton_comp, lvls1, 1);

  bool stayed_same = (lemon::countNodes(graph1) == nodes_before) && (lemon::countEdges(graph1) == edges_before);

  lemon::ListGraph graph2;
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> labMaps2(2);
  std::vector<std::vector<size_t>> comp2 = {{0,1,2}};
  construct_graph(comp2, graph2, labMaps2[0]);

  std::vector<std::vector<std::uint64_t>> signatures = {
    {1,1,  2,2},
    {1,1,  4,4},
    {5,5,  6,6}
  }; // Collision between bands 0 and 0 of comp 1 and 0, two stays singleton.

  std::vector<std::pair<size_t,size_t>> lvls2 = {{2,2}, {2,2}};
  // @test this little example if it works in theory.

  recursive_step(signatures, graph2, labMaps2, comp2[0], lvls2, 1);

  bool correct_nodes = (lemon::countNodes(graph2) == 3) && (lemon::countEdges(graph2) == 2);

  bool no_unlabeled = (lemon::countNodes(graph2) == (int)(labMaps2[0].size() + labMaps2[1].size()));


  std::cout << "\n=== Test the recursive_step function ===\n";
  std::cout << "Instant stop after encountering a singleton: " << colored_bool(stayed_same) << "\n";
  std::cout << "Example fits the expectation: " << colored_bool(correct_nodes) << "\n";
  std::cout << "No node without label: " << colored_bool(no_unlabeled) << "\n";
}

int main(){
  test_construct_graph_tower();
  test_construct_graph();
  test_extract_bands();
  test_connect_bands_tower();
  test_connect_bands();
  test_connected_components();
  test_recursive_step();
}