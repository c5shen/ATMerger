'''
Created on 10.28.2021 by Chengze Shen

Merger of all query alignments to form the final alignment. The merging step
is exactly the same one from PASTA and UPP (by transitivity).
'''

import os, sys, re
import time
import argparse
from utility import Alignment, CompactAlignment, compact 
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
def mergeAlignments(indir, outpath, t, pool):
    start = time.time()
    inpaths = os.popen('ls {}'.format(indir)).read().split('\n')[:-1]
    inpaths = ['{}/{}'.format(indir, x) for x in inpaths]
    assert len(inpaths) > 0

    # split paths into NUM_CPUS chunks
    chunks = []
    chunk_size = ceil(len(inpaths) / t)
    for i in range(0, len(inpaths), chunk_size):
        chunks.append(inpaths[i:min(i+chunk_size, len(inpaths))])

    # initialize Pool for multiprocessing
    #pool = Pool(Configs.num_cpus)
    merged_alns = list(pool.map(sequential_merger, chunks))
    #pool.close()
    #pool.join()

    # for the merged chunks, merge them into one alignment 
    final_aln = merged_alns[0]
    for i in range(1, len(merged_alns)):
        final_aln.merge_in(merged_alns[i])
    
    final_aln.write(outpath, 'FASTA')
    end = time.time()
    print('Time used to merge: {0:.2f} s'.format(end - start))

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
    args = parser.parse_args()
    indir, outpath, t = args.indir, args.outpath, int(args.cores)

    pool = ProcessPoolExecutor(t)
    _ = pool.submit(dummy)

    mergeAlignments(indir, outpath, t, pool)

if __name__ == "__main__":
    main()
