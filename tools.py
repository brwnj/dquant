#!/usr/bin/env python
# encoding: utf-8

import sys
import editdist as ed
from itertools import islice
from collections import Counter

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

def distance(a, b):
    """find best edit distance between two strings of potentially uneven length.
    """
    la, lb = len(a), len(b)
    if la < lb:
        return distance(b, a)
    if la == lb:
        return ed.distance(a, b)
    else:
        dists = []
        for i in xrange(0, la - lb):
            dists.append(ed.distance(a[i:i+lb], b))
        return min(dists)

def _initialize_matching(counter):
    """docstring for initialize_matching"""
    seqs = list(counter)
    seqs.sort(key = len, reverse = True)
    return set(), seqs, len(seqs)

def process_exact_substring(counter):
    """docstring for process_exact_substring"""
    seen, seqs, leng = _initialize_matching(counter)
    for i, target in enumerate(seqs, start=1):
        if i % 10000 == 0: print >>sys.stderr, ">> processed %d/%d" % (i, leng)
        lt = len(target)
        if target in seen: continue
        for query in seqs:
            if len(query) == lt: continue
            if query in target:
                counter[target] += counter[query]
                counter[query] = 0
                seen.add(query)
    counter += Counter()
    return counter
         
def process_similar(counter, n):
    """For each matching sequence, add its count to longest sequence."""
    seen, seqs, leng = _initialize_matching(counter)
    for i, target in enumerate(seqs, start=1):
        if i % 10000 == 0: print >>sys.stderr, ">> processed %d/%d" % (i, leng)
        if target in seen: continue
        seen.add(target)
        for query in seqs:
            if query in seen: continue
            if distance(target, query) <= n:
                counter[target] += counter[query]
                counter[query] = 0
                seen.add(query)
    counter += Counter()
    return counter
