#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Generate core genome from mmseqs2 clustering', epilog= "This code was written by Bart Cuypers and Alessandro Brozzi. Please make a comment in the GitHub repository for support.")
parser.add_argument('Protein_FASTA', metavar='allproteins-fasta', type=str,
                    help='concatenated FASTA file with all proteins')
parser.add_argument('Clusters_tsv', metavar='clusters-tsv', type=str,
                    help='tsv file with mmseqs2 clusters')
parser.add_argument('-t', metavar='threshold', type=float, default = 100,
                    help='[0-100] Minimum percentage of strains in which a "core protein" should appear')

args = vars(parser.parse_args())

clusters = args['Clusters_tsv']
threshold = args['t']

# Make dict with protein names as keys, and strain names as values
protein_strain = {}
# Make set with all strain names
allstrains = set([])
fasta_sequences = SeqIO.parse(args['Protein_FASTA'],'fasta')
for protein in fasta_sequences:
    strain = protein.description.split('[')[-1].strip(']')
    protein_strain[protein.id] = strain
    allstrains.add(strain)

# make dict with the cluster representative as key, and as value the strains in which the cluster is found
cluster_strains = {}
with open(clusters) as myreadfile:
    for line in myreadfile:
        cluster_representative, member = line.strip('\n').split('\t')
        if cluster_representative in cluster_strains:
            cluster_strains[cluster_representative].add(protein_strain[member])
        else:
            cluster_strains[cluster_representative] = set([protein_strain[member]])

towrite = []

fasta_sequences = SeqIO.parse(args['Protein_FASTA'],'fasta')
for protein in fasta_sequences:
    if protein.id in cluster_strains:
        if len(cluster_strains[protein.id])/len(allstrains) * 100 >= threshold:
            towrite.append(protein)

SeqIO.write(towrite, "coregenome.fasta", "fasta")