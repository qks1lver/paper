#!/usr/bin/env python3

import os
import sys
import subprocess
from datetime import datetime
from time import time


def local_alignment(p_query, p_ref):

    print('Query: %s' % p_query)
    print('Database: %s' % p_ref)

    p_file, _ = os.path.splitext(p_query)
    p_out = p_file + '_localalignment-%s' % datetime.now().strftime('%Y%m%d%H%M%S') + '.txt'

    n_cpu = os.cpu_count()
    run_query = ['blastp', '-query', p_query, '-db', p_ref, '-out', p_out, '-num_threads', str(n_cpu), '-outfmt', '10 qseqid sseqid bitscore evalue qcovs']
    print('Running (%d CPUs): %s ...' % (n_cpu, ' '.join(run_query)))
    subprocess.run(run_query)

    print('Output: %s' % p_out)

    return(p_out)

def parse_local_alignment(p_alignment):

    print('Parsing: %s' % p_alignment)

    p_file, _ = os.path.splitext(p_alignment)
    p_out = p_file + '_parsed.txt'

    q2h = {}

    with open(p_alignment, 'r') as fi, open(p_out, 'w+') as fo:
        query = ''
        hits = []
        for l in fi:
            if l.startswith('Query='):
                query = l.split('=')[1].strip().split()[0]
                l = fi.readline()
                while not l.startswith('Length='):
                    l = fi.readline()
                seq_length = l.split('=')[1].strip()
                query += '/%s' % seq_length
            if query and l.startswith('Sequences producing'):
                _ = fi.readline()
                l = fi.readline().strip()
                while l:
                    tmp = l.split()
                    hits.append((tmp[0], tmp[-2], tmp[-1]))
                    l = fi.readline().strip()
            if query and hits:
                if query in q2h:
                    print('query ID conflict: %s' % query)
                q2h[query.split('/')[0]] = hits[0]
                str_hits = '|'.join(['%s/%s/%s' % (a,b,c) for a,b,c in hits])
                _ = fo.write('%s|%s\n' % (query, str_hits))
                query = ''
                hits = []
                fo.flush()

    print('Parsed %d entires to: %s' % (len(q2h), p_out))

    return(p_out, q2h)

def assess_reciprocals(q2h1, q2h2, p_out='', title1='', title2=''):

    if not title1:
        title1 = 'Query_1'

    if not title2:
        title2 = 'Query_2'

    if not p_out:
        p_out = 'reciprocal_output-%s.csv' % datetime.now().strftime('%Y%m%d%H%M%S')
    else:
        _, p_ext = os.path.splitext(p_out)
        if p_ext != '.csv':
            p_out += '.csv'

    overlaps = []
    with open(p_out, 'w+') as f:
        _ = f.write('%s,%s,E-value\n' % (title1, title2))
        for query in q2h1:
            hit = q2h1[query][0]
            e1 = float(q2h1[query][2])
            if hit in q2h2:
                hit_query = q2h2[hit][0]
                e2 = float(q2h2[hit][2])
                if hit_query == query:
                    ev = min(e1,e2)
                    overlaps.append((query, hit, ev))
                    _ = f.write('%s,%s,%.2E\n' % (query, hit, ev))

    print('Wrote %d matches to: %s' % (len(overlaps), p_out))

    return(p_out, overlaps)

if __name__ == '__main__':

    t0 = time()

    args = sys.argv
    opts = ''
    val = ''

    p_assembly = args[1]
    p_control = args[2]
    if len(args) > 3:
        opts = args[3]
    if len(args) > 4:
        val = args[4]

    if 'b' in opts:
        if val == '0':
            _ = local_alignment(p_assembly, p_control)
        else:
            _ = local_alignment(p_control, p_assembly)

    if 'p' in opts:
        _, q2h_assembly = parse_local_alignment(p_assembly)
        _, q2h_control = parse_local_alignment(p_control)

        if 'a' in opts:
            p_assessment = '../data/reciprocal_local_alignments.csv'
            _ = assess_reciprocals(q2h_assembly, q2h_control, p_assessment, 'assembly', 'control')

    t = (time() - t0) / 3600

    print('Completed in %.1f hours' % t)
