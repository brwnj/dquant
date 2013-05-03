#!/usr/bin/env python
# encoding: utf-8
"""
De novo identification and quantification of sequence data.
"""
import sys
import iqseq_utils as du

__version__ = "0.1"

def run_quantify(args):
    print >>sys.stderr, ">> collapsing identical sequences (1/5)"
    reads = du.process_exact_fastq(args.fastq)
    print >>sys.stderr, ">> constructing suffix tree (2/5)"
    t = du.construct_trie(reads)
    print >>sys.stderr, ">> collapsing identical subsequences (3/5)"
    reads = du.process_exact_substring(reads, t)
    print >>sys.stderr, ">> optimizing suffix tree (4/5)"
    t = du.construct_trie(reads)
    print >>sys.stderr, ">> collapsing similar sequences (5/5)"
    reads = du.process_similar(reads, t, args.mismatch)
    for seq, count in reads.iteritems():
        print "%s\t%d" % (seq, count)

def run_consensus(args):
    # work in progress...
    print >>sys.stderr, ">> collapsing identical sequences"
    seqs = process_exact(args.bins, args.cutoff)
    print >>sys.stderr, ">> collapsing identical subsequences"
    seqs = t.process_exact_substring(seqs)
    print >>sys.stderr, ">> collapsing similar sequences"
    seqs = t.process_similar(seqs, args.mismatch)
    s = list(seqs)
    s.sort(key = len, reverse = True)
    print "\n".join(s)

def main(args):
    args.func(args)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    subp = p.add_subparsers(help='commands')
    
    fquant = subp.add_parser('quantify', description="Find and quantify unique and similar sequences within a FASTQ.", help="")
    fquant.add_argument("fastq", metavar="FASTQ", help="reads to process")
    fquant.add_argument("-m", dest="mismatch", type=int, default=3,
            help="mismatch tolerance when grouping bins [%(default)s]")
    fquant.set_defaults(func=run_quantify)
    
    fcons = subp.add_parser('consensus', description="Build consensus of sequences across all samples.", help="")
    fcons.add_argument('bins', metavar='BINS', nargs="+",
            help="results of `quantify.py`")
    fcons.add_argument('-c', dest='cutoff', type=int, default=10,
            help="minimum bin count to be considered [%(default)s]")
    fcons.add_argument('-m', dest='mismatch', type=int, default=3,
            help="mismatch tolerance when grouping bins [%(default)s]")
    
    args = p.parse_args()
    main(args)