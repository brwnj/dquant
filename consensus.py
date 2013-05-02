#!/usr/bin/env python
# encoding: utf-8
"""
Build consensus of sequences across all samples.
"""
import sys
import tools as t
from collections import Counter
from toolshed import reader

def process_exact(files, cutoff):
    c = Counter()
    for f in args.bins:
        for l in reader(f, header=['seq','count']):
            if l['count'] < cutoff: continue
            c.update([l['seq']])
    return c

def main(args):
    print >>sys.stderr, ">> collapsing identical sequences"
    seqs = process_exact(args.bins, args.cutoff)
    print >>sys.stderr, ">> collapsing identical subsequences"
    seqs = t.process_exact_substring(seqs)
    print >>sys.stderr, ">> collapsing similar sequences"
    seqs = t.process_similar(seqs, args.mismatch)
    s = list(seqs)
    s.sort(key = len, reverse = True)
    print "\n".join(s)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bins', metavar='BINS', nargs="+",
            help="results of `quantify.py`")
    p.add_argument('-c', dest='cutoff', type=int, default=10,
            help="minimum bin count to be considered [%(default)s]")
    p.add_argument('-m', dest='mismatch', type=int, default=3,
            help="total mismatch tolerance when grouping bins [%(default)s]")
    args = p.parse_args()
    main(args)
