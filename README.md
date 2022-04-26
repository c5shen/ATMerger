# AlignmentTransitivityMerger
A standalone alignment merging tool that utilizes transitivity that merges a set of overlapping alignments. The original implementation comes from PASTA and UPP (see [original codes availibility](#original-codes-availability)).

# Inputs
1. A list of overlapping alignments to merge (all alignments need to overlap with each other on some taxa)
2. The output file path
3. The number of cores (for the purpose of multi-processing)


# Usage
```bash
$ python3 merger.py [directory with only input alignments] [output path] [number of cores]
```

# I am working on ...
1. Allowing users to specify the merging order of alignments. This should eliminate the issue with not all alignments having overlapped compartments (e.g., a set of alignments obtained by neighbor pairs of nodes in a spanning tree, each node representing a cluster of taxa).
2. In the case of 1, the merger should operate in sequential order (no multi-processing).

# Original Codes Availability
1. [PASTA](https://github.com/smirarab/pasta) (see `pasta/alignment.py`)
2. [SEPP/UPP](https://github.com/smirarab/sepp) (see `sepp/alignment.py`)
