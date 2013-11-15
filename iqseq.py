#!/usr/bin/env python
# encoding: utf-8
"""
De novo identification and quantification of sequence data.
"""
import sys
import numpy as np
import pandas as pd
import os.path as op
from toolshed import nopen, reader
from Bio import trie, triefind
from collections import Counter
from itertools import islice, ifilterfalse

__version__ = "0.4"

def read_fastq(fh):
    """FASTQ parser that yields name, seq, and qual."""
    while True:
        values = list(islice(fh, 4))
        if len(values) == 4:
            id1, seq, id2, qual = values
        elif len(values) == 0:
            raise StopIteration
        else:
            raise EOFError("unexpected end of file")
        assert id1.startswith('@')
        assert id2.startswith('+')
        assert len(seq) == len(qual)
        yield id1[1:-1], seq[:-1], qual[:-1]

def trim_seq(seq, base=5):
    """round length of sequence to nearest `base`"""
    return seq[:int(base * round(len(seq)/base))]

def process_exact_fastq(fastq, n):
    """Group identical reads using a Counter. Returns Counter."""
    c = Counter()
    with nopen(fastq) as fh:
        for name, seq, qual in read_fastq(fh):
            seq = trim_seq(seq, 4)
            if len(seq) < n: continue
            c.update([seq])
    return c

def process_exact_txt(files, cutoff):
    """returns Counter from multiple quantify runs"""
    c = Counter()
    for f in files:
        for l in reader(f, header=['seq','count']):
            if int(l['count']) < cutoff: continue
            c.update([l['seq']])
    return c

def process_counted(fp, sample_id, cutoff):
    """method called to get sequence counts during `run_matrix`."""
    sequence_counts = Counter()
    library_size = 0
    for l in reader(fp, header=['seq','count']):
        count = int(l['count'])
        library_size += count
        if count < cutoff: continue
        sequence_counts[l['seq']] = count
    return sequence_counts, library_size

def get_seq_bins(fp):
    """fp to text. no ints on this input."""
    c = Counter()
    for l in nopen(fp):
        c[l.strip()] = 0
    return c

def chunker(it, n):
    # chunker('AAAABBBC', 4) --> AAAA AAAB AABB ABBB BBBC
    return [it[i:i+n] for i in xrange(0, len(it)+1-n, 1)]

def construct_simple_trie(counter):
    t = trie.trie()
    for seq, count in counter.iteritems():
        t[seq] = count
    return t

def construct_complex_trie(counter, lengths=None):
    t = trie.trie()
    seqs = list(counter)
    seqs.sort(key=len, reverse=True)
    if lengths is None:
        lengths = sorted(set([len(k) for k in seqs]))
    for seq in seqs:
        seq_len = len(seq)
        for l in lengths:
            if l > seq_len: continue
            for subseq in chunker(seq, l):
                if t.has_key(subseq): continue
                if subseq == seq:
                    t[seq] = counter[seq]
                else:
                    t[subseq] = seq
    return t

def process_exact_substring(counter, t):
    """use triefind.find to gather identical substring matches"""
    seqs = list(counter)
    seqs.sort(key=len, reverse=True)
    for seq in seqs:
        l = len(seq)
        for (match, start, end) in triefind.find(seq, t):
            if len(match) == l: continue
            counter[seq] += counter[match]
            counter[match] = 0
    counter += Counter()
    return counter

def unique_everseen(iterable, key=None):
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in ifilterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def process_similar(counter, t, n):
    """trie is composed of sequences being compared."""
    seqs = list(counter)
    seqs.sort(key=len, reverse=True)
    lengths = sorted(set([len(k) for k in seqs]))
    progress = 100
    to_process = len(seqs)
    for i, seq in enumerate(seqs, start=1):
        if i % progress == 0:
            progress = int(progress * 1.5)
            print >>sys.stderr, "processed %d of %d" % (i, to_process)
        if counter[seq] == 0: continue
        for (k, v, dist) in unique_everseen(t.get_approximate(seq, n), lambda (m,c,d): m):
            if dist == 0 or k == seq: continue
            if type(v) is int:
                counter[seq] += counter[k]
                counter[k] = 0
            else:
                # k is a subsequence; therefore add seq to v
                counter[v] += counter[seq]
                counter[seq] = 0
    counter += Counter()
    return counter

def process_similar_matrix(bins, seqs, t, n):
    """
    bins - sequence bins
    seqs - sequences to bin
    t    - trie
    n    - mismatches
    
    returns Counter
    """
    sample_seqs = list(seqs)
    sample_seqs.sort(key=len, reverse=True)
    to_process = len(sample_seqs)
    progress = 100
    for i, seq in enumerate(sample_seqs, start=1):
        if i % progress == 0:
            progress = int(progress * 1.5)
            print >>sys.stderr, "    >> processed {i} of {to_process}".format(**locals())
        # returning bins to which the sequence belongs
        for (k, v, dist) in unique_everseen(t.get_approximate(seq, n), lambda (m,c,d): m):
            if type(v) is int:
                bins[k] += seqs[seq]
                # set to zero? avoids adding counts to multiple bins
                seqs[seq] = 0
            else:
                bins[v] += seqs[seq]
                seqs[seq] = 0
    return bins

def scalefactor(counts):
    # mask inf and nan
    ma = np.ma.masked_invalid(counts)
    return np.exp(np.ma.median(ma))

def write_table(d, library_sizes, norm=None):
    if norm == "deseq":
        # details: http://genomebiology.com/2010/11/10/R106
        df = pd.DataFrame(d)
        # log of counts
        lg = df.apply(np.log)
        # per sample: exponential(median(log(counts) - geometric mean))
        sf = lg.sub(lg.mean(axis=1), axis=0).apply(scalefactor, axis=0)
        # apply scaling
        df = df.div(sf, axis=1)
    elif norm == "totalcount":
        df = pd.DataFrame(d)
        mean_total_count = float(sum(library_sizes.values())) / len(library_sizes)
        # apply total count scaling
        for col in df.columns:
            denominator = float(library_sizes[col])
            assert denominator > 0, \
                    "No counts found in sample {sampleid}".format(sampleid=col)
            df[col] = (df[col] / denominator) * mean_total_count
    else:
        df = pd.DataFrame(d)
    df.to_csv(sys.stdout, sep="\t")

def run_quantify(args):
    print >>sys.stderr, ">> collapsing identical sequences (1/5)"
    reads = process_exact_fastq(args.fastq, args.cutoff)
    print >>sys.stderr, ">> constructing suffix tree (2/5)"
    t = construct_simple_trie(reads)
    print >>sys.stderr, ">> collapsing identical subsequences (3/5)"
    reads = process_exact_substring(reads, t)
    print >>sys.stderr, ">> optimizing suffix tree (4/5)"
    t = construct_complex_trie(reads)
    print >>sys.stderr, ">> collapsing similar sequences (5/5)"
    reads = process_similar(reads, t, args.mismatch)
    for seq, count in reads.iteritems():
        print "%s\t%d" % (seq, count)

def run_consensus(args):
    """does basically the same thing as quantify, except doesn't print out the
    count.
    """
    print >>sys.stderr, ">> collapsing identical sequences (1/5)"
    seqs = process_exact_txt(args.bins, args.cutoff)
    print >>sys.stderr, ">> constructing suffix tree (2/5)"
    t = construct_simple_trie(seqs)
    print >>sys.stderr, ">> collapsing identical subsequences (3/5)"
    seqs = process_exact_substring(seqs, t)
    print >>sys.stderr, ">> optimizing suffix tree (4/5)"
    t = construct_complex_trie(seqs)
    print >>sys.stderr, ">> collapsing similar sequences (5/5)"
    seqs = process_similar(seqs, t, args.mismatch)
    s = list(seqs)
    s.sort(key=len, reverse=True)
    print "\n".join(s)

def run_matrix(args):
    d = {}
    samples = set()
    to_process = len(args.counts)
    library_sizes = {}
    for i, f in enumerate(args.counts, start=1):
        sample = op.splitext(op.basename(f))[0]
        samples.add(sample)
        assert len(samples) == i
        print >>sys.stderr, (">> processing sample {sample} "
                                "({i}/{to_process})").format(
                                                        sample=sample,
                                                        i=i,
                                                        to_process=to_process)
        d[sample] = {}
        # the sequence counts of current sample and total library size
        seqs, library_size = process_counted(f, sample, args.cutoff)
        library_sizes[sample] = library_size
        seq_lengths = sorted(set([len(k) for k in list(seqs)]))
        seq_bins = get_seq_bins(args.consensus)
        # trie based on sequences of bins at lengths of query sequences
        t = construct_complex_trie(seq_bins, seq_lengths)
        # process the sequences
        counts = process_similar_matrix(seq_bins, seqs, t, args.mismatch)
        for k, v in counts.iteritems():
            d[sample][k] = v
    write_table(d, library_sizes, args.norm)

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
            help=("minimum allowable count for individual sample sequences"))
    fmat.add_argument("-m", dest="mismatch", type=int, default=3,
            help="mismatch tolerance when grouping bins")
    fmat.add_argument("-n", dest="norm", default=None, choices=['deseq', 'totalcount'],
            help=("output normalized table using either DESeq or total count method"))
    fmat.set_defaults(func=run_matrix)
    
    args = p.parse_args()
    main(args)
