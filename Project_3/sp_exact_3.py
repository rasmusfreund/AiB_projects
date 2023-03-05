from typing import Sequence, TextIO
import time
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""Commands for adding command line arguments"""

parser = argparse.ArgumentParser(prog = "MSA_exact.py",
                                 usage = "python3 %(prog)s sequences.fasta",
                                 description = "Align 3 sequences with an exact match")
parser.add_argument("seqs", help = "Name of .FASTA file containing the three sequences to be aligned")
parser.add_argument("--out", help = "Produces a .FASTA file containing the aligned sequences as output",
                    action = "store_true")
parser.add_argument("--runtime", help = "For runtime analysis purposes",
                    action = "store_true")
parser.add_argument("--score", help = "Prints the optimal alignment cost of the alignments",
                    action = "store_true")
args = parser.parse_args()



##############################################################
### Values in scoreMatrix can be changed to any integer,   ###
### however, the alignment will always try to find the     ###
### minimum score, meaning that MATCHES should always      ###
### have the lowest score.                                 ###
##############################################################

scoreMatrix = { 'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
                'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
                'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
                'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}}

##############################################################
### The GAPCOST value will be used for both linear gap     ###
### cost alignment, and affine gap cost alignment.         ###
### In the case of affine gap cost, the GAPCOST will be    ###
### used as the gap opening cost, while the gap extension  ###
### cost must be defined through the command line.         ###
### See "global_alignment.py -h" for help.                 ###
##############################################################

GAPCOST = 5

##############################################################
############# Modify below this at your own risk #############
##############################################################


def fastaParse(fasta_file: Sequence) -> list[Sequence]:
    sequenceList = list(SeqIO.parse(fasta_file, "fasta"))
    container = []
    for i in range(len(sequenceList)):
        container.append(sequenceList[i].seq)
    return container


def empty_matrix(m: Sequence, n: Sequence, o: Sequence) -> list[list[list]]:
    """Creates a matrix of size len(m + 1) x len(n + 1) x len(o + 1) and fills with None"""

    outer_list = [[[None for _ in range(len(o) + 1)] for _ in range(len(n) + 1)] for _ in range(len(m) + 1)]
    return outer_list


def initiate_matrix(m: Sequence, n: Sequence, o: Sequence) -> list[list[list]]:
    """Fills the three outer sides of the cube with gapscores"""

    matrix = empty_matrix(m, n, o)

    for i in range(len(m) + 1):
        for j in range(len(n) + 1):
            matrix[i][j][0] = i * GAPCOST + j * GAPCOST

    for i in range(1, len(m) + 1):
        for k in range(1, len(o) + 1):
            matrix[i][0][k] = i * GAPCOST + k * GAPCOST

    for j in range(len(n) + 1):
        for k in range(1, len(o) + 1):
            matrix[0][j][k] = j * GAPCOST + k * GAPCOST

    return matrix


def fill_matrix(seq1: Sequence, seq2: Sequence, seq3: Sequence, score_matrix: dict) -> list[list[list[int]]]:
    """Fills the remaining nodes of the matrix """

    S_matrix = initiate_matrix(seq1, seq2, seq3)

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            for k in range(1, len(seq3) + 1):
                score_diagonal_3d = S_matrix[i-1][j-1][k-1] + score_matrix[seq1[i-1]][seq2[j-1]]\
                    + score_matrix[seq1[i-1]][seq3[k-1]] + score_matrix[seq2[j-1]][seq3[k-1]]
                score_diagonal_ij = S_matrix[i-1][j-1][k] + score_matrix[seq1[i-1]][seq2[j-1]] + GAPCOST * 2
                score_diagonal_ik = S_matrix[i-1][j][k-1] + score_matrix[seq1[i-1]][seq3[k-1]] + GAPCOST * 2
                score_diagonal_jk = S_matrix[i][j-1][k-1] + score_matrix[seq2[j-1]][seq3[k-1]] + GAPCOST * 2
                score_prev_i = S_matrix[i-1][j][k] + GAPCOST * 2
                score_prev_j = S_matrix[i][j-1][k] + GAPCOST * 2
                score_prev_k = S_matrix[i][j][k-1] + GAPCOST * 2

                S_matrix[i][j][k] = min(score_diagonal_3d, score_diagonal_ij, score_diagonal_ik, score_diagonal_jk,\
                    score_prev_i, score_prev_j, score_prev_k)

    if args.score: # Get optimal alignment score if requested
        print(S_matrix[-1][-1][-1])

    return S_matrix

#########################################################################################################
def traceback_direction(matrix: list[list], row: int , col: int, match_score: int) -> str:
    """Finds the node from which the current node's score comes from."""

    diagonal_score = matrix[row - 1][col - 1]
    up_score = matrix[row - 1][col]
    left_score = matrix[row][col - 1]
    node_score = matrix[row][col]

    if node_score == diagonal_score + match_score:
        return 'diagonal'
    elif node_score == left_score + GAPCOST:
        return 'left'
    elif node_score == up_score + GAPCOST:
        return 'up'


def get_base(sequence: Sequence, position: int) -> str:
    """Simply grabs a single nucleotide from a sequence based on a positional argument"""
    return sequence[position - 1]


def alignment(seq1_file: TextIO, seq2_file: TextIO, score_matrix: list[list]) -> str:
    """Creates a possible alignment from two fasta files."""

    # Initial string to save alignment
    align1 = ""
    align2 = ""

    # Load sequences and fill out alignment scores
    seq1, seq2 = fastaParse(seq1_file), fastaParse(seq2_file)
    seq1str, seq2str = seq1[0].seq, seq2[0].seq
    filled_matrix = fill_matrix(seq1str, seq2str, score_matrix)

    # Idx for bottom right node in the matrix
    row = len(seq1str)
    col = len(seq2str)

    # Backtrack and create alignment
    while row > 0 and col > 0:
        match_score = score_matrix[get_base(seq1str, row)][get_base(seq2str, col)]

        trace_direction = traceback_direction(filled_matrix, row, col, match_score)
        match trace_direction:
            case 'diagonal':
                align1 = get_base(seq1str, row) + align1
                align2 = get_base(seq2str, col) + align2
                row -= 1
                col -= 1
            case 'up':
                align1 = get_base(seq1str, row) + align1
                align2 = '-' + align2
                row -= 1
            case 'left':
                align1 = '-' + align1
                align2 = get_base(seq2str, col) + align2
                col -= 1

    return align1, align2



##############################################################
##################### Run the algorithm ######################
##############################################################



def main():
    if args.runtime:
        st = time.time()

    # seq_file = args.seqs,

    # aligned = alignment(seq_file, scoreMatrix)
    # print(aligned[0] + '\n' + aligned[1])

    m = "GTTCCGAAAGGCTAGCGCTAGGCGCC"
    n = "ATGGATTTATCTGCTCTTCG"
    o = "TGCATGCTGAAACTTCTCAACCA"
    # matrix = initiate_matrix(m, n, o)
    matrix = fill_matrix(m, n, o, scoreMatrix)


    if args.runtime:
        et = time.time()
        elapsed_time = et - st
        print("Execution time:", elapsed_time, "seconds")

    if args.out:
        """Produces an output-file in .fasta format using the original information
        regarding the sequences."""
        file1, file2 = fastaParse(seq1_file), fastaParse(seq2_file)

        # Convert alignments to sequence objects
        seqs = [Seq(x) for x in aligned]

        #Convert sequence objects to sequence records (yes, it's this tedious...)
        seq1_R = SeqRecord(seqs[0],
                           id = file1[0].id,
                           name = file1[0].name,
                           description = file1[0].description)
        seq2_R = SeqRecord(seqs[1],
                           id = file2[0].id,
                           name = file2[0].name,
                           description = file2[0].description)

        records_to_write = [seq1_R, seq2_R]

        SeqIO.write(records_to_write, file1[0].name+"_"+file2[0].name+".fasta", "fasta")


if __name__ == '__main__':
    main()