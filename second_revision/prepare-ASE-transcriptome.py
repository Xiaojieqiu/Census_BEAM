#!/usr/bin/env python

from __future__ import print_function

import sys
import gzip
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

if len(sys.argv) != 4:
    print("Usage: ./prepare-ASE-transcriptome.py transcriptome.fa[.bgz] snps.txt[.gz] output.fa[.gz]", file=sys.stderr)
    print("transcriptome.fa must be either uncompressed or bgzip-ed", file=sys.stderr)
    sys.exit(1)

print()
print("Note: assumes paternal=REF, maternal=REF+SNPs.", file=sys.stderr)
print("TODO: make this more flexible", file=sys.stderr)
print("Also Note: Python is slow, maybe get coffee and come back later...")
print()

transcriptome_file = sys.argv[1]

if sys.argv[2][-3:] == ".gz":
    snps_file = gzip.open(sys.argv[2], "r")
else:
    snps_file = open(sys.argv[2], "r")

if sys.argv[3][-3:] == ".gz":
    output = gzip.open(sys.argv[3], "w")
else:
    output = open(sys.argv[3], "w")

print("Reading reference transcriptome", file=sys.stderr)
transcriptome = SeqIO.index(transcriptome_file, "fasta")


snps = {} # transcript id -> [(cDNA pos, ref, alt)]
          # ref/alt on appropriate strand
complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }

print("Reading SNPs", file=sys.stderr)
for line in snps_file:
    fields = line.strip().split("\t")
    ref  = fields[1]
    alt  = fields[2]
    t_id = fields[3]
    try:
        pos = int(fields[4]) - 1 # 1-index to 0-index
    except:
        continue
    if fields[6] == "-1":
        ref = complement[ref]
        alt = complement[alt]
    if t_id in snps:
        snps[t_id].append((pos, ref, alt))
    else:
        snps[t_id] = [(pos, ref, alt)]

snps_file.close()

print()
print("Creating paternal transcriptome", file=sys.stderr)
print("(just adding _P to transcript ids, but I/O still takes a while)", file=sys.stderr)
for t_id in transcriptome:
    transcript = transcriptome[t_id]
    transcript.id   += "_P"
    SeqIO.write(transcript, output, "fasta")

print()
print("Creating maternal transcriptome", file=sys.stderr)
for t_id in transcriptome:
    transcript = transcriptome[t_id]
    transcript.id   += "_M"
    mat_seq = MutableSeq(str(transcript.seq), generic_dna)
    if t_id in snps:
        for (pos, ref, alt) in snps[t_id]:
            if mat_seq[pos] == ref:
                mat_seq[pos] = alt
            else:
                print("WARNING: %s:%d ref allele = %s, but SNP says it's %s. Skipping SNP..." % (t_id, pos+1, mat_seq[pos], ref), file=sys.stderr)
    transcript.seq = mat_seq
    SeqIO.write(transcript, output, "fasta")
    
output.close()

