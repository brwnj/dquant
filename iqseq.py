#!/usr/bin/env python
# encoding: utf-8
"""
De novo identification and quantification of sequence data.
"""
import sys
import os.path as op
import iqseq_utils as iu

__version__ = "0.2"

def run_quantify(args):
    print >>sys.stderr, ">> collapsing identical sequences (1/5)"
    reads = iu.process_exact_fastq(args.fastq, args.cutoff)
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
    s.sort(key=len, reverse=True)
    print "\n".join(s)

def run_matrix(args):
    seqbins, lengths = iu.get_seq_bins(args.consensus)
    d = {}
    samples = set()
    to_process = len(args.counts)
    for i, f in enumerate(args.counts, start=1):
        sample = op.splitext(op.basename(f))[0]
        samples.add(sample)
        assert len(samples) == i
        print >>sys.stderr, ">> processing sample %s (%d/%d)" % \
                                                        (sample, i, to_process)
        d[sample] = {}
        # the sequence counts of each sample
        seqs = iu.process_counted(f)
        t = iu.construct_complex_trie(seqs, lengths)
        # process the sequences
        seqbins = iu.process_similar_matrix(seqbins, seqs, t, args.mismatch)
        for k, v in seqbins.iteritems():
            d[sample][k] = v
    iu.write_table(d)

def main(args):
    args.func(args)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    subp = p.add_subparsers(help='commands')
    
    fquant = subp.add_parser('quantify',
            description="Find and quantify unique and similar sequences \
            within a FASTQ.",
            help="quantify unique and similar sequences")
    fquant.add_argument("fastq", metavar="FASTQ", help="reads to process")
    fquant.add_argument("-c", dest="cutoff", type=int, default=18,
            help="minimum allowable seq length [%(default)s]")
    fquant.add_argument("-m", dest="mismatch", type=int, default=3,
            help="mismatch tolerance when grouping bins [%(default)s]")
    fquant.set_defaults(func=run_quantify)
    
    fcons = subp.add_parser('consensus',
            description="Build consensus of sequences across all samples.",
            help="build observed sequence library")
    fcons.add_argument('bins', metavar='BINS', nargs="+",
            help="results of `quantify`")
    fcons.add_argument('-c', dest='cutoff', type=int, default=100,
            help="minimum allowable count [%(default)s]")
    fcons.add_argument('-m', dest='mismatch', type=int, default=3,
            help="mismatch tolerance when grouping bins [%(default)s]")
    fcons.set_defaults(func=run_consensus)
    
    fmat = subp.add_parser('matrix', description="Generate counts matrix",
            help="generate counts matrix")
    fmat.add_argument("consensus", metavar="CONSENSUS", help="result of `consensus`")
    fmat.add_argument("counts", metavar="COUNTS", nargs="+", help="results of `quantify`")
    fmat.add_argument("-m", dest="mismatch", type=int, default=3, help="")
    fmat.set_defaults(func=run_matrix)
    
    args = p.parse_args()
    main(args)
