# Alignment Transitivity Merger
This is a standalone alignment merging tool that utilizes transitivity to merge a set of overlapping alignments. The original implementation comes from PASTA and UPP (see [original codes availibility](#original-codes-availability)). The core idea is that, if two alignments have some overlapping taxa, we could use them as "anchors" to link the two alignments together to form a larger alignment with the union of taxa from both.

### Limitation
Arbitrary merging of two alignments is not trivial, and this program is only designed to merge two or more alignments that share some sub-alignments. For example, if two alignments, each having 501 taxa, share exactly a sub-alignment of size 500, we can use this merger to merge them together to form a 502-taxa alignment.

# Inputs
1. A list of overlapping alignments to merge (all alignments need to overlap with each other on some taxa)
2. The output file path
3. The number of cores (for the purpose of multi-processing)


# Usage
```bash
$ python3 merger.py -d [directory with only input alignments] -o [output path] -t [number of cores]
```
Use `python3 merger.py -h` to see more details.

### Example
```bash
$ python3 merger.py -d examples/data -o merged.fasta -t 1
```

# I am working on ...
1. Allowing users to specify the merging order of alignments. This should eliminate the issue with not all alignments having overlapping compartments (e.g., a set of alignments obtained by neighbor pairs of nodes in a spanning tree, each node representing a cluster of taxa).
2. In the case of 1, the merger should operate in sequential order (no multi-processing).

# Original Codes Availability
1. [PASTA](https://github.com/smirarab/pasta) (see `pasta/alignment.py`)
2. [SEPP/UPP](https://github.com/smirarab/sepp) (see `sepp/alignment.py`)
