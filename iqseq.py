#!/usr/bin/env python
# encoding: utf-8
"""
De novo identification and quantification of sequence data.
"""
import sys
import os.path as op
import iqseq_utils as iu

__version__ = "0.3"

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
    d = {}
    samples = set()
    to_process = len(args.counts)
    for i, f in enumerate(args.counts, start=1):
        sample = op.splitext(op.basename(f))[0]
        samples.add(sample)
        assert len(samples) == i
        print >>sys.stderr, (">> processing sample {sample} "
                                "({i}/{to_process})").format(**locals())
        d[sample] = {}
        # the sequence counts of current sample
        seqs = iu.process_counted(f, args.cutoff)
        seq_lengths = sorted(set([len(k) for k in list(seqs)]))
        seq_bins = iu.get_seq_bins(args.consensus)
        # trie based on sequences of bins at lengths of query sequences
        t = iu.construct_complex_trie(seq_bins, seq_lengths)
        # process the sequences
        counts = iu.process_similar_matrix(seq_bins, seqs, t, args.mismatch)
        for k, v in counts.iteritems():
            d[sample][k] = v
    iu.write_table(d, args.norm)

def main(args):
    args.func(args)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__)
    subp = p.add_subparsers(help='commands')
    
    fquant = subp.add_parser('quantify',
            description=("Find and quantify unique and similar sequences "
                            "within a FASTQ."),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help="quantify unique and similar sequences")
    fquant.add_argument("fastq", metavar="FASTQ", help="reads to process")
    fquant.add_argument("-c", dest="cutoff", type=int, default=18,
            help="minimum allowable seq length")
    fquant.add_argument("-m", dest="mismatch", type=int, default=3,
            help="mismatch tolerance when grouping bins")
    fquant.set_defaults(func=run_quantify)
    
    fcons = subp.add_parser('consensus',
            description="Build consensus of sequences across all samples.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help="build observed sequence library")
    fcons.add_argument('bins', metavar='BINS', nargs="+",
            help="results of `quantify`")
    fcons.add_argument('-c', dest='cutoff', type=int, default=100,
            help="minimum allowable count")
    fcons.add_argument('-m', dest='mismatch', type=int, default=3,
            help="mismatch tolerance when grouping bins")
    fcons.set_defaults(func=run_consensus)
    
    fmat = subp.add_parser('matrix', description="Generate counts matrix",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            help="generate counts matrix")
    fmat.add_argument("consensus", metavar="CONSENSUS",
            help="result of `consensus`")
    fmat.add_argument("counts", metavar="COUNTS", nargs="+",
            help="results of `quantify`")
    fmat.add_argument('-c', dest='cutoff', type=int, default=100,
            help=("minimum allowable count for individual sample "
                    "sequences"))
    fmat.add_argument("-m", dest="mismatch", type=int, default=3,
            help="mismatch tolerance when grouping bins")
    fmat.add_argument("-n", dest="norm", action="store_true",
            help=("output scaling factor normalized table using method "
                    "developed by Anders and Huber for DESeq: "
                    "http://genomebiology.com/2010/11/10/R106"))
    fmat.set_defaults(func=run_matrix)
    
    args = p.parse_args()
    main(args)
