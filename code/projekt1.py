from Bio import SeqIO, Phylo, AlignIO
from Bio.Seq import Seq
from pathlib import Path
import random
import math
import numpy as np
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import sys, getopt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings


def sequence(input_filename):
    """Load sequence from fasta file."""
    with open(input_filename) as inputfile:
        for record in SeqIO.parse(inputfile, "fasta"):
            S = record.seq
    return S

def aa2na(seq):
    """Convert protein sequence to DNA sequence."""
    AA2NA = {
    "A": list("GCT,GCC,GCA,GCG".split(",")),
    "R": list("CGT,CGC,CGA,CGG,AGA,AGG".split(",")),
    "N": list("AAT,AAC".split(",")),
    "D": list("GAT,GAC".split(",")),
    "C": list("TGT,TGC".split(",")),
    "Q": list("CAA,CAG".split(",")),
    "E": list("GAA,GAG".split(",")),
    "G": list("GGT,GGC,GGA,GGG".split(",")),
    "H": list("CAT,CAC".split(",")),
    "I": list("ATT,ATC,ATA".split(",")),
    "L": list("TTA,TTG,CTT,CTC,CTA,CTG".split(",")),
    "K": list("AAA,AAG".split(",")),
    "M": list("ATG".split(",")),
    "F": list("TTT,TTC".split(",")),
    "P": list("CCT,CCC,CCA,CCG".split(",")),
    "S": list("TCT,TCC,TCA,TCG,AGT,AGC".split(",")),
    "T": list("ACT,ACC,ACA,ACG".split(",")),
    "W": list("TGG".split(",")),
    "Y": list("TAT,TAC".split(",")),
    "V": list("GTT,GTC,GTA,GTG".split(",")),
    "*": list("TAA,TGA,TAG".split(","))
    }
    na_seq = [random.choice(AA2NA.get(c, ["---"])) for c in seq]
    return "".join(na_seq)

def P_Aj(a,j,t, alfa=1):
    """a, j - nucleotides,
    t - time,
    alfa - JC parameter; 
    Calculate the transformation of nucleotide "a" into nucleotide "j" 
    at time "t" using the Jukes-Cantor model"""
    if  a==j:
        p=1/4+3/4*math.exp(-4*alfa*t)
    else:
        p=1/4-1/4*math.exp(-4*alfa*t)
    return p

def choose_nuc(l):
    """We draw a mutation. l-list with the probabilities of converting a given nucleotide into another. 
    The list in the following positions has the probabilities with which the nucleotide transforms into A, C, T or G.
    The function returns the nucleotide that we get as a result of the mutation."""
    r = random.uniform(0, 1)
    if 0 <= r < l[0]:
        return 'A'
    elif l[0] <= r < l[0]+l[1]:
        return 'C'
    elif l[0]+l[1] <= r < l[0]+l[1]+l[2]:
        return 'T'
    else:
        return 'G'

def step(s, d):
    """s - sequence, d - time quantum;
    Simulation of the one step of the algorithm from the draw_new_seq function."""
    nuc = 'ACTG'
    new_seq = ''
    for i in s: # go through the nucleotides in the input sequence
        L = [] 
        for j in nuc: # go through the nucleotides A,C,T,G (in that order)
            L.append(P_Aj(i,j,d))
        new_seq += choose_nuc(L) # create new sequence (change the nucleotide i into A, C, T, or G)
    return new_seq

def draw_new_seq(S, alfa, d, K):
    """S - input sequence,
    alpha - the alpha parameter of the JC69 model,
    d - time quantum,
    K - number of simulation steps from the input tree branch;
    Return sequence after K mutations and after the time d."""

    for k in range(K): # steps
        new_seq = step(S,d)
        S = new_seq
    return S


def get_leaves_names(T):
    """T - tree;
    Return list of leaves names from the tree"""
    leaves = T.get_terminals()
    L= []
    for l in leaves:
        L.append(l.name)
    return L

def get_branches_len(l, T):
    """
    l - leaf name,
    T - tree;
    Create list of branch lengths from root to leaf."""
    X = []
    path_to_leaf = T.get_path(l)
    for clade in path_to_leaf:
        X.append(clade.branch_length)
    return X


def generate_new_sequence(l, S, T, alpha, d):
    """
    l - leaf name
    S - root sequence
    T - tree
    alpha - JC parameter
    d - time
    Create sequence from root by mutations."""
    branches_len_to_leaf = get_branches_len(l,T)
    for b in branches_len_to_leaf:
        S = draw_new_seq(S, alpha, d, int(b))
    return S

def all_seq(S, T, alpha, d):
    """
    Draw sequences A - E.
    """
    leaves = get_leaves_names(T)
    list_of_sequences = []
    for l in leaves:
        list_of_sequences.append(generate_new_sequence(l,S,T,alpha,d))
    return list_of_sequences

def translate_DNA_to_aa(list_of_sequences):
    """Translate DNA sequences into protein sequences."""
    proteins = []
    for i in range(len(list_of_sequences)):
        coding_dna = Seq(list_of_sequences[i])
        t = coding_dna.translate()
        proteins.append(t)
    return(proteins)

def seqences_to_fasta_file(P, leaves, T):
    """
    P - list of sequences
    leaves - names of leaves in the tree
    T - tree
    Create file bialka.fa with protein sequences."""
    fileout = Path('bialka.fa')
    fileout.touch(exist_ok=True)  # will create file, if it exists will do nothing
    ofile = open(str(fileout), "w")
    for i in range(len(P)):
        ofile.write(">" + leaves[i] + "\n" + str(P[i]) + "\n")
    ofile.close()


def align_clustalw(filein):
    fileout = Path('bialka_aligned.aln')
    fileout.touch(exist_ok=True)
    cline = ClustalwCommandline("clustalw", infile=str(filein), outfile=str(fileout))
    cline()


def main():
    warnings.simplefilter("ignore")

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--fasta_input_file", default="PAH.fa", help="Input file with protein sequence.")
    parser.add_argument("-t", "--tree", default="tree", help="Input tree file.")
    parser.add_argument("-d", "--time", default=100000, help="Time")

    args = vars(parser.parse_args())
     
    # 1. 
    input_file = args["fasta_input_file"]
    pah_aa_seq = sequence(input_file)
    pah_na_seq = aa2na(pah_aa_seq)

    #2.
    tree_file = args["tree"]
    tree = Phylo.read(tree_file, "newick")

    # 3.
    d = int(args["time"])
    S = all_seq(pah_na_seq, tree, 0.4, d)

    # 4.
    P = translate_DNA_to_aa(S)
    leaves = get_leaves_names(tree)
    seqences_to_fasta_file(P,leaves,tree)

    # 5.
    align_clustalw('bialka.fa') 

    # 6.
    aln = AlignIO.read('bialka_aligned.aln', 'clustal')
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(aln)

    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dm)
    nj_tree = constructor.nj(dm)

    Phylo.write(upgma_tree, "upgma_tree.xml", "phyloxml")
    Phylo.write(nj_tree, "nj_tree.xml", "newick")


if __name__ == "__main__":
    main()
