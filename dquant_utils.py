#!/usr/bin/env python
# encoding: utf-8

from toolshed import nopen
from Bio import trie, triefind
from collections import Counter
from itertools import islice, ifilterfalse

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

def process_exact_fastq(fastq):
    """Group identical reads using a Counter. Returns Counter."""
    c = Counter()
    with nopen(fastq) as fh:
        for name, seq, qual in read_fastq(fh):
            c.update([trim_seq(seq, 4)])
    return c

def process_exact_txt(files, cutoff):
    """returns Counter from multiple quantify runs"""
    c = Counter()
    for f in args.bins:
        for l in reader(f, header=['seq','count']):
            if l['count'] < cutoff: continue
            c.update([l['seq']])
    return c

def chunker(it, n):
    return [it[i:i+n] for i in xrange(0, len(it)+1-n, 1)]

def construct_trie(counter):
    """build suffix tree from Counter"""
    t = trie.trie()
    for seq, count in counter.iteritems():
        t[seq] = count
    return t

def process_exact_substring(counter, tree):
    """use triefind.find to gather identical substring matches"""
    seqs = list(counter)
    seqs.sort(key=len, reverse=True)
    for seq in seqs:
        l = len(seq)
        for (match, start, end) in triefind.find(seq, tree):
            if len(match) == l: continue
            counter[seq] += counter[match]
            counter[match] = 0
    counter += Counter()
    return counter

def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
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

def process_similar(counter, tree, n):
    seqs = list(counter)
    seqs.sort(key=len, reverse=True)
    lengths = sorted(set([len(k) for k in seqs]))
    for seq in seqs:
        if counter[seq] == 0: continue
        for l in lengths:
            if l > len(seq): continue
            for subseq in chunker(seq, l):
                for (match, count, dist) in unique_everseen(\
                            tree.get_approximate(subseq, n), lambda (m,c,d): m):
                    if dist == 0 or count == 0 or match == seq: continue
                    counter[seq] += counter[match]
                    counter[match] = 0
                    tree[match] = 0
    counter += Counter()
    return counter
