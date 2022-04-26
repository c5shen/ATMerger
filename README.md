# Alignment Transitivity Merger (ATMerger)
This is a standalone alignment merging tool that utilizes transitivity to merge a set of overlapping alignments. The original implementation comes from PASTA and UPP (see [original codes availibility](#original-codes-availability)). The core idea is that, if two alignments have some overlapping taxa, we could use them as "anchors" to link the two alignments together to form a larger alignment with the union of taxa from both.

### Limitations
* Arbitrary merging of two alignments is not trivial, and this program is only designed to merge two or more alignments that share some sub-alignments. For example, if two alignments, each having 501 taxa, share exactly a sub-alignment of size 500, we can use this merger to merge them together to form a 502-taxa alignment.
* Additionally, the order of merging is important, if not all alignments are overlapping. For example, if two non-overlapping alignments are attempted to be merged together, the outcome alignment will be entirely gapped between the two sets of taxa. Any additional alignment information that connects these taxa will not be able to override.
* To extend on the previous point, a simpler way to put is that: **ATMerger needs all overlapping alignments to ONLY come from a spanning tree**, where each node represents a set of taxa, and each edge represents an alignment in our input. Therefore, the input alignments are overlapping. Our goal is to find the alignment on _all_ taxa that appear in the spanning tree by merging these input alignments represented by edges.

# Inputs
1. A list of overlapping alignments to merge (all alignments need to overlap with each other on some taxa)
2. The output file path
3. The number of cores (for the purpose of multi-processing)
4. (Optional) a file stating the specific order that the merging should occur


# Usage
```bash
$ python3 merger.py -d [directory with only input alignments] -o [output path] -t [number of cores]
```
Use `python3 merger.py -h` to see more details.

#### Example 1
```bash
$ python3 merger.py -d examples/exp1 -o exp1.fasta -t 1
```

#### Example 2 (i.e., given specific merging order)
```bash
$ python3 merger.py -d examples/exp2 --order examples/exp2_order.txt -o exp2.fasta -t 1
```

#### Example 2 wrong order
```bash
$ python3 merger.py -d examples/exp2 --order examples/exp2_wrong_order.txt -o exp2_wrong.fasta -t 1
```


# I am working on ...
1. (DONE) Allowing users to specify the merging order of alignments. This should eliminate the issue with not all alignments having overlapping compartments (e.g., a set of alignments obtained by neighbor pairs of nodes in a spanning tree, each node representing a cluster of taxa).
2. (COROLLARILY DONE) In the case of 1, the merger should operate in sequential order (no multi-processing).

# Original Codes Availability
1. [PASTA](https://github.com/smirarab/pasta) (see `pasta/alignment.py`)
2. [SEPP/UPP](https://github.com/smirarab/sepp) (see `sepp/alignment.py`)
