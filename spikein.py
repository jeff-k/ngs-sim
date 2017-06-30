#!/usr/bin/env python3
"""
In-silico NGS data synthesis
"""

import itertools
import argparse
import os.path
import sys
import numpy as np

def write_fastq(seq, rname):
    """Return a 4 line string representing a fastq record"""
    qual = []
    for _ in range(len(seq)):
        #qual.append(chr(int(np.floor(np.random.normal(70, 5)))))
        qual.append('H')
    return "{}\n{}\n{}\n{}".format(rname, seq, '+', ''.join(qual))

def read_fasta(f):
    """Parse a single reference fasta file into a string"""
    lines = []
    for line in open(f):
        if line[0] == '>':
            continue
        lines.append(line.strip())
    return ''.join(lines)

def rc(seq):
    # complementary ambiguous nucleotide codes
    d = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D'}
    return ''.join((map(lambda x: d[x], reversed(seq))))

def mutate(ops, seq):
    """Apply a list of operations to a reference sequence"""
    
    for p, d, i in reversed(sorted(ops, key=lambda x: x[2])):
        if i == '.':
            seq = seq[:p] + seq[p+d:]
        else:
            seq = seq[:p] + i + seq[p+d:]
    return seq

class Op:
    def __init__(self, pos, d, ins):
        """Initialize an operation
        
        pos: the reference coordinate of the mutation
        d: the number of deleted bases
        ins: the string of inserted bases
        """
        self.d = d
        self.ins = ins
        self.pos = pos

    def __call__(self, template):
        return template[:self.pos] + self.ins + template[self.pos + self.d:]

class ChemProfile:
    def __init__(self, frag_mean, read_size, adapters=None, sigma=100,
                 error_profile=False):
        self.frag_mean = frag_mean
        self.frag_std = sigma
        self.read_size = read_size
        if error_profile:
            self.err_mat = np.array([
                    [0.97, 0.01, 0.01, 0.01],
                    [0.01, 0.97, 0.01, 0.01],
                    [0.01, 0.01, 0.97, 0.01],
                    [0.01, 0.01, 0.01, 0.97]
                ])
        else:
            self.err_mat = None

    def fragment(self, seq):
        s = int(np.floor(np.random.normal(self.frag_mean, self.frag_std)))
        p = int(np.floor(np.random.randint(0, len(seq) - s)))
        frag = seq[p:p+s]
        return frag[:self.read_size], rc(frag[-self.read_size:])

    def amplify(self, seq):
        seq = list(seq)
        for i in range(len(seq)):

            # Mutate sequence with the error matrix
            if self.err_mat is not None:
                base = seq[i]
                ps = self.err_mat[['A','C','G','T'].index(base)]
                seq[i] = np.random.choice(['A','C','G','T'], p=ps)
        return ''.join(seq)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('ref')
    parser.add_argument('-p', default=1)
    parser.add_argument('-n', '--num', default=10000)
    parser.add_argument('-m', '--mutations', '-m', nargs='+', 
                        default=[])
    parser.add_argument('-t', '--template', '-t', default=None)
    parser.add_argument('-o', '--out', default="out")
    args = parser.parse_args()

    prof = ChemProfile(2000, 151)
    reads = []

    fq1_path = "{}_R1.fastq".format(args.out)
    fq2_path = "{}_R2.fastq".format(args.out)
    if os.path.exists(fq1_path) or os.path.exists(fq2_path):
        print("Output path already exists!", file=sys.stderr)
        exit(1)

 
    mutations = []
    for m in args.mutations:
        p = 1
        if ',' in m:
            vcf, p = m.split(',')
        else:
            vcf = m
        if not os.path.exists(vcf):
            print("Variant file {} does not exist!".format(vcf),
                  file=sys.stderr)
            exit(1)
        mutations.append((vcf, int(p)))

    y = 0
    total_p = sum(map(lambda x: x[1], mutations))
    conseq_p = max(map(lambda x: x[1], mutations))
    total_reads = args.num
    conseq = None
    
    # Randomize a run id to avoid read name collisions
    id_space = "01234567890ABCDEF"
    m_id = ""
    for _ in range(6):
        m_id += id_space[np.random.randint(len(id_space))]

    fq1 = open(fq1_path, 'w')
    fq2 = open(fq2_path, 'w')

    for m, f in mutations:
        y += 1
        ops = []
        seq = read_fasta(args.ref).upper()
        for line in open(m):
            if line[0] == '#':
                continue
            p, d, i = line.strip().split()
            ops.append((int(p),int(d),i.upper()))
        proportion = f / total_p

        sub_seq = mutate(ops, seq)
        if f == conseq_p:
            conseq = str(seq)

        for count in range(int(total_reads * proportion)):
            r1, r2 = prof.fragment(prof.amplify(sub_seq))
            print(write_fastq(r1,
                "@SYNTH-{}:0:0000000:0:0:{}:{} 1:N:0:0".format(m_id, y, count)),
                 file=fq1)
            print(write_fastq(r2,
                "@SYNTH-{}:0:0000000:0:0:{}:{} 2:N:0:0".format(m_id, y, count)),
                file=fq2)
    print(conseq, file=open("conseq.fasta", 'w'))
