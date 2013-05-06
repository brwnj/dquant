#!/usr/bin/env python
# encoding: utf-8
"""
De novo identification and quantification of sequence data.
"""
import sys
import iqseq_utils as iu

__version__ = "0.1"

def run_quantify(args):
    print >>sys.stderr, ">> collapsing identical sequences (1/5)"
    reads = iu.process_exact_fastq(args.fastq, args.length)
    print >>sys.stderr, ">> constructing suffix tree (2/5)"
    t = iu.construct_simple_trie(reads)
    print >>sys.stderr, ">> collapsing identical subsequences (3/5)"
    reads = iu.process_exact_substring(reads, t)
    print >>sys.stderr, ">> optimizing suffix tree (4/5)"
    t = iu.construct_complex_trie(reads)
    print >>sys.stderr, ">> collapsing similar sequences (5/5)"
    reads = iu.process_similar(reads, t, args.mismatch)
    for seq, count in reads.iteritems():
        print "%s\t%d" % (seq, count)

def run_consensus(args):
    """does basically the same thing as quantify, except doesn't print out the
    count.
    """
    print >>sys.stderr, ">> collapsing identical sequences (1/5)"
    seqs = iu.process_exact_txt(args.bins, args.cutoff)
    print >>sys.stderr, ">> constructing suffix tree (2/5)"
    t = iu.construct_simple_trie(seqs)
    print >>sys.stderr, ">> collapsing identical subsequences (3/5)"
    seqs = iu.process_exact_substring(seqs, t)
    print >>sys.stderr, ">> optimizing suffix tree (4/5)"
    t = iu.construct_complex_trie(seqs)
    print >>sys.stderr, ">> collapsing similar sequences (5/5)"
    seqs = iu.process_similar(seqs, t, args.mismatch)
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
    fquant.add_argument("-l", dest="length", type=int, default=18,
            help="minimum allowable seq length [%(default)s]")
    fquant.add_argument("-m", dest="mismatch", type=int, default=3,
            help="mismatch tolerance when grouping bins [%(default)s]")
    fquant.set_defaults(func=run_quantify)
    
    fcons = subp.add_parser('consensus', description="Build consensus of sequences across all samples.", help="")
    fcons.add_argument('bins', metavar='BINS', nargs="+",
            help="results of `quantify`")
    fcons.add_argument('-c', dest='cutoff', type=int, default=10,
            help="minimum allowable count [%(default)s]")
    fcons.add_argument('-m', dest='mismatch', type=int, default=3,
            help="mismatch tolerance when grouping bins [%(default)s]")
    fcons.set_defaults(func=run_consensus)
    
    args = p.parse_args()
    main(args)