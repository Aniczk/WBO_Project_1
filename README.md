## Usage:
```
python ./projekt1.py -i <.fa input file> -t <input tree> -d <time>
```
## Default:
-i PAH.fa <br />
-t tree <br />
-d 100000 <br />

## Output:
bialka.fa <br />
nj_tree.xml <br />
upgma_tree.xml <br />
bialka_aligned.aln <br />
bialka.dnd <br />

## Main problem:
Simulate and reconstruct the evolution of the sequence according to the phylogenetic tree contained in the input file.
Assumption: we have a sequence coding a human protein (default: phenylalanine hydroxylase in the PAH.fa file). We imagine that in the future 5 different proteins will evolve from it (by speciation lasting millions of years), A, B, C, D and E according to the scheme from the tree file.
(We assume that the protein coding sequence is at the root of the tree and the edge lengths correspond to % point mutations).

## Details:
1. Read a protein sequence from a .fa file and generate one of the possible DNA sequences that could encode it (random choice).
2. Load a given tree from a file in the Newick format.
3. Randomly generate A-E sequences using the Markov model with assumptions compatible with the Jukes-Cantor model. Expected number of mutations is related to the edge lengths in the tree.
4. Translate the generated sequences into proteins and save the generated sequences in the file fasta (bialka.fa).
5. Basing on these protein sequences, create a multi alignment using the clustalw program.
6. Basing on multi alignment, generate phylogenetic trees using the upgma and nj methods.
