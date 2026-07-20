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
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> labMaps1(2);
  std::vector<std::vector<size_t>> comp1 = {{0,1,2}};
  construct_graph(comp1, graph1, labMaps1[0]);

  std::vector<std::vector<std::uint64_t>> signatures = {
    {1,1,  2,2},
    {1,1,  4,4},
    {5,5,  6,6}
  }; // Collision between bands 0 and 0 of comp 1 and 0, two stays singleton.

  std::vector<std::pair<size_t,size_t>> lvls1 = {{2,2}, {2,2}};
  // @test this little example if it works in theory.

  recursive_step(signatures, graph1, labMaps1, comp1[0], lvls1, 1);

  bool correct_nodes = (lemon::countNodes(graph1) == 3) && (lemon::countEdges(graph1) == 2);

  bool no_unlabeled = (lemon::countNodes(graph1) == (int)(labMaps1[0].size() + labMaps1[1].size()));


  std::cout << "\n=== Test the recursive_step function ===\n";
  std::cout << "Example fits the expectation: " << colored_bool(correct_nodes) << "\n";
  std::cout << "No node without label: " << colored_bool(no_unlabeled) << "\n";
}

void test_get_clusters(){
  // @test when there is no concrete labMap given, assume there are no clusters Thus, the returned sequence-to-cluster Map should be empty as well.
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> emptyLabMaps;
  std::vector<std::unordered_map<size_t,const std::vector<size_t>*>> empty_clusters = get_clusters(emptyLabMaps);

  bool is_empty = empty_clusters.empty();

  // @test a given circumstance.
  
  lemon::ListGraph graph1;
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> labMaps1(2);

  // Construct the first level.
  construct_graph(std::vector<std::vector<size_t>>({{0,1,2}, {3}}), graph1, labMaps1[0]);

  // Construct the second level. In this instance, the levels will not be connected, but that will not be important, as the function does not care about connectivity.
  construct_graph(std::vector<std::vector<size_t>>({{0,1},{2}}), graph1, labMaps1[1]); //{3} will be excluded, since it was a singleton to begin with.

  std::vector<std::unordered_map<size_t,const std::vector<size_t>*>> clusters = get_clusters(labMaps1);

  bool fit_expectations = ( (clusters.size() == 2) && (clusters[0].size() == 4) && (clusters[1].size() == 3) ); // We expect two levels and every Sequence to be one element.

  // @test when the expectations fit, we can use this to test other things, such as the identity of the clusters.

  bool correct_pointer_lvl0 = ( (clusters[0][0] == clusters[0][1]) && (clusters[0][1] == clusters[0][2]) );  // Sequence 0, 1 and 2 were in the same cluster.
  bool different_pointer_lvl0 = (clusters[0][0] != clusters[0][3]); // Sequence 3 was the singleton.

  bool correct_cluster_seq0_lvl0 = (*(clusters[0][0]) == std::vector<size_t>({0,1,2}));
  bool correct_cluster_seq0_lvl1 = (*(clusters[1][0]) == std::vector<size_t>({0,1}));

  std::cout << "\n=== Test the get_clusters function ===\n";
  std::cout << "No Nodes result in an empty cluster gathering: " << colored_bool(is_empty) << "\n";
  std::cout << "The overall structure of the cluster gathering fits the expectations: " << colored_bool(fit_expectations) << "\n";
  std::cout << "Same clustered have the same pointers: " << colored_bool(correct_pointer_lvl0) << "\n";
  std::cout << "Different clustered sequences have different pointers: " << colored_bool(different_pointer_lvl0) << "\n";
  std::cout << "Sequence 0 references to the correct cluster on Level 0: " << colored_bool(correct_cluster_seq0_lvl0) << "\n";
  std::cout << "Sequence 0 references to the correct cluster on Level 1: " << colored_bool(correct_cluster_seq0_lvl1) << "\n";
}

void test_binning(){
  // @test it is not important how things look; when bins == 0, the resulting vector should be empty, as there can be no binning.
  lemon::ListGraph graph0;
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> labMaps0(2);
  construct_graph(std::vector<std::vector<size_t>>({{0,1,2}, {3}}), graph0, labMaps0[0]);
  construct_graph(std::vector<std::vector<size_t>>({{0,1}, {2}, {3}}), graph0, labMaps0[1]);
  std::vector<std::unordered_map<size_t, const std::vector<size_t>*>> clusts0 = get_clusters(labMaps0);

  std::vector<std::vector<size_t>> zero_bins = binning(labMaps0, clusts0, 0, 1);
  bool no_bins = zero_bins.empty();

  // @test Just a little example to see if it fits the expectations.
  lemon::ListGraph graph1;
  std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, debugHasher>> labMaps1(2);
  construct_graph(std::vector<std::vector<size_t>>({{0,1}, {2,3,4}}), graph1, labMaps1[0]);
  construct_graph(std::vector<std::vector<size_t>>({{0},{1}, {2,3}, {4}}), graph1, labMaps1[1]);
  std::vector<std::unordered_map<size_t, const std::vector<size_t>*>> clusts1 = get_clusters(labMaps1);

  std::vector<std::vector<size_t>> example_bin_1 = binning(labMaps1, clusts1, 3, 1); 
  std::vector<std::vector<size_t>> expected_bin_1 = {{0}, {4}, {1}}; // Since we have 3 bins and three singletons, we expect only them and in exact this order (0, 1 are from the same supercluster)
  bool match_expectation_1 = (example_bin_1 == expected_bin_1);

  std::vector<std::vector<size_t>> example_bin_2 = binning(labMaps1, clusts1, 4, 1); 
  std::vector<std::vector<size_t>> expected_bin_2 = {{0}, {4}, {1}, {2}}; // Adding one more allows us to look at {2,3}, where we take the first element, 2
  bool match_expectation_2 = (example_bin_2 == expected_bin_2);

  // @test we want to check that no sequence gets binned multiple times.
  std::vector<std::vector<size_t>> no_dup_bin = binning(labMaps1, clusts1, 2, 3);
  std::unordered_map<size_t,size_t> occurence_count;
  for(const std::vector<size_t>& bin : no_dup_bin)
    for(size_t seq : bin)
      occurence_count[seq] += 1;
  
  bool no_dups = true;
  for(auto& [seq,count] : occurence_count){
    if(count > 1) no_dups = false;
  }

  std::cout << "\n=== Test the binning function ===\n";
  std::cout << "bins = 0 results in an empty binning: " << colored_bool(no_bins) << "\n";
  std::cout << "First, Only Singletons will get clusters and diversity is top down: " << colored_bool(match_expectation_1) << "\n";
  std::cout << "When given, the first element of a greater cluster is taken: " << colored_bool(match_expectation_2) << "\n";
  std::cout << "No Sequence is binned multiple times: " << colored_bool(no_dups) << "\n";
}

int main(){
  test_construct_graph_tower();
  test_construct_graph();
  test_extract_bands();
  test_connect_bands_tower();
  test_connect_bands();
  test_connected_components();
  test_recursive_step();
  test_get_clusters();
  test_binning();
}