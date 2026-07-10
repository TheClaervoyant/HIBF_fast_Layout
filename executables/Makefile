CXX = g++
CXXFLAGS = -std=c++23 -I../seqan3/include -I../xxHash
LDFLAGS = -llemon -lz -lbz2
SRC = ../code

all: function_tests generate_clusters

function_tests: ${SRC}/fast_construct_test.cpp ${SRC}/fast_construct_hashing.cpp ${SRC}/fast_construct_graph.cpp
	${CXX} ${CXXFLAGS} $^ -o $@ $(LDFLAGS)

generate_clusters: ${SRC}/fast_construct_graph.cpp ${SRC}/test_clustering.cpp ${SRC}/fast_construct_hashing.cpp ${SRC}/generate_clusters.cpp
	${CXX} ${CXXFLAGS} $^ -o $@ ${LDFLAGS}

clean:
	rm -f function_tests generate_clusters

.PHONY: all clean
