# Goal
The Hierarchical Interleaved Bloom Filter (HIBF), introduced by Mehringer, Seiler et al., is a space- and time-efficient data structure for approximate sequence queries.
Its layout computation assigns sequences to bins, where the quality of the assignment affects the HIBF's performance and space consumption.

This project aims to develop a faster layout computation by pre-clustering sequences based on their similarity. Grouping similar sequences into the same bin improves the filter's overall performance.

Clustering is performed via Locality-Sensitive Hashing (LSH) using One Permutation Hashing (OPH). This yields an efficient approximation of the Jaccard similarity between two sequences.
The resulting similarity graph is then analyzed hierarchically to identify clusters at multiple levels of granularity, enabling a multi-resolution view before constructing the layout.

# Dependencies
In order to compile the executables (by using make function_tests or make generate_clusters), following dependencies are needed:
1. **The Lemon Graph Library** - must be installed locally (`sudo apt install liblemon-dev`)
2. **seqan3** - must be cloned into the parent directory of `executables` (`git clone https://github.com/seqan/seqan3.git`)
3. **zlib** and **libbz2** - required by SeqAn3 (`sudo apt install zlib1g-dev libbz2-dev`)
4. **xxHash** must be cloned into the parent directory of `executables` (`git clone https://github.com/Cyan4973/xxHash.git`)

A **C++23-compatible compiler** is required.

The expected directory structure would look like this: 

```text
HIBF_fast_Layout
├── code
├── executables
    └── Makefile
├── seqan3
└── xxHash
```
