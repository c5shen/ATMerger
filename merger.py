'''
Created on 10.28.2021 by Chengze Shen

Merger of all query alignments to form the final alignment. The merging step
is exactly the same one from PASTA and UPP (by transitivity).
'''

import os, sys, re
import time
import argparse
from utility import Alignment, CompactAlignment, compact 
from multiprocessing import cpu_count
from concurrent.futures.process import ProcessPoolExecutor
from math import ceil

'''
function to merge a set of input paths (to alignments) sequentially
'''
def sequential_merger(inpaths):
    init_index = 0
    init_aln = Alignment(); init_aln.read_file_object(inpaths[init_index])
    new_aln = compact(init_aln)
    for i in range(init_index + 1, len(inpaths)):
        inpath = inpaths[i]
        
        frag_aln = Alignment(); frag_aln.read_file_object(inpath)
        new_aln.merge_in(compact(frag_aln))
    return new_aln

'''
function to take in a set of result paths for merging, and write
the merged alignment to an output path
'''
def mergeAlignments(inpaths, outpath, t, pool):
    # split paths into NUM_CPUS chunks
    chunks = []
    chunk_size = ceil(len(inpaths) / t)
    for i in range(0, len(inpaths), chunk_size):
        chunks.append(inpaths[i:min(i+chunk_size, len(inpaths))])

    merged_alns = list(pool.map(sequential_merger, chunks))

    # for the merged chunks, merge them into one alignment 
    final_aln = merged_alns[0]
    for i in range(1, len(merged_alns)):
        final_aln.merge_in(merged_alns[i])
    final_aln.write(outpath, 'FASTA')

def dummy():
    pass

def main():
    parser = argparse.ArgumentParser(description='Merge a set of overlapping'
            ' alignments')
    parser.add_argument('-d', '--indir', required=True, type=str,
            help='the directory that contains only the overlapping alignments')
    parser.add_argument('-o', '--outpath', required=True, type=str,
            help='the merged alignment output destination')
    parser.add_argument('-t', '--cores', required=False, type=int, default=1,
            help='the number of cores to use for multi-processing (default: 1)')

    optional_group = parser.add_argument_group(
            "Optional".upper(), ' '.join(["These are optional fields."]))
    parser.groups = dict()
    parser.groups['optional'] = optional_group
    optional_group.add_argument('--order', metavar='textfile',
            required=False, type=str, default=None,
            help='the order of files that the merger needs to follow. '
                'Make sure that all lines in the order file are in the '
                'following format: '
                '1) alignment file name, or 2) absolute path to the alignment file')
    args = parser.parse_args()
    indir, outpath, t = args.indir, args.outpath, int(args.cores)
    if t < 1:
        raise ValueError('Number of cpus needs to be greater than 0')
    if t > cpu_count(): 
        print('Input number cores > max, readjusting to {}'.format(cpu_count()))
        t = cpu_count()
    order = args.order

    # get the list of alignments to merge
    inpaths = os.popen('ls {}'.format(indir)).read().split('\n')[:-1]
    inpaths = ['{}/{}'.format(indir, x) for x in inpaths]
    assert len(inpaths) > 0, 'input directory is empty'

    # reorder the inpaths if order is given
    if order != None:
        if not os.path.exists(order):
            print('The order file does not exist, ignoring...')
        else:
            with open(order, 'r') as f:
                lines = f.read().split('\n')[:-1]
                lengths = [len(x.split('/')) for x in lines]
                if sum(lengths) == len(lines):
                    inpaths = ['{}/{}'.format(indir, x) for x in lines]
                else:
                    inpaths = lines
    start = time.time()
    if t > 1:
        pool = ProcessPoolExecutor(t)
        _ = pool.submit(dummy)
        mergeAlignments(inpaths, outpath, t, pool)
    else:
        merged = sequential_merger(inpaths)
        merged.write(outpath, 'FASTA')
    end = time.time()
    print('Time used to merge: {0:.2f} s'.format(end - start))

if __name__ == "__main__":
    main()
