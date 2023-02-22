from typing import Sequence, TextIO
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""Commands for adding command line arguments"""

parser = argparse.ArgumentParser(prog = "global_alignment.py",
                                 usage = "python3 %(prog)s seq1.fasta seq2.fasta",
                                 description = "Align two sequences using either global linear gap score or affine gap score")
parser.add_argument("seq1", help = "Name of .FASTA file containing the first sequence to be aligned")
parser.add_argument("seq2", help = "Name of .FASTA file containing the second sequence to be aligned")
parser.add_argument("--affine", help = "changes gap score from linear to affine; \
    the gap extension score must be provided after this argument. E.g., '--affine 5'.",
                    action = "store")
parser.add_argument("--out", help = "Produces a .FASTA file containing the aligned sequences as output",
                    action = "store_true")
args = parser.parse_args()



##############################################################
### Values in scoreMatrix can be changed to any integer,   ###
### however, the alignment will always try to find the     ###
### minimum score, meaning that MATCHES should always      ###
### have the lowest score.                                 ###
##############################################################

scoreMatrix = {'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
               'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
               'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
               'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}}

##############################################################
### The GAPCOST value will be used for both linear gap     ###
### cost alignment, and affine gap cost alignment          ###
##############################################################

GAPCOST = 5

##############################################################
############# Modify below this at your own risk #############
##############################################################

if args.affine: # If affine is called, get gap extend value
    GAP_EXTEND = int(args.affine)


def fastaParse(fasta_file: Sequence) -> list[Sequence]:
    sequenceList = list(SeqIO.parse(fasta_file, "fasta"))
    sequenceList[0].seq = sequenceList[0].seq.upper()
    return sequenceList


def empty_matrix(m: Sequence, n: Sequence) -> list[list]:
    """Creates a matrix of size len(m) x len(n) and fills with None"""

    outer_list = []
    for _ in range(len(m) + 1):
        inner_list = [None for _ in range(len(n) + 1)]
        outer_list.append(inner_list)
    return outer_list


def initiate_matrix(m: Sequence, n: Sequence) -> list[list]:
    """Fills out first row and first column
    of the matrix using the gapcost."""

    matrix = empty_matrix(m, n)
    matrix[0][0] = 0
    for i in range(len(m) + 1):
        if args.affine:
            matrix[i][0] = GAPCOST + (i-1) * GAP_EXTEND
        else:
            matrix[i][0] = i * GAPCOST
    for j in range(len(n) + 1):
        if args.affine:
            matrix[0][j] = GAPCOST + (j-1) * GAP_EXTEND
        else:
            matrix[0][j] = j * GAPCOST
    return matrix


def fill_matrix(seq1: Sequence, seq2: Sequence, score_matrix: dict) -> list[list[int]]:
    """Fills the remaining nodes of the matrix """

    S_matrix = initiate_matrix(seq1, seq2)

    if args.affine: # Affine gap cost matrices
        D_matrix = initiate_matrix(seq1, seq2)
        I_matrix = initiate_matrix(seq1, seq2)
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                D_matrix[i][j] = min(S_matrix[i-1][j] + GAPCOST + GAP_EXTEND, D_matrix[i-1][j] + GAP_EXTEND)
                I_matrix[i][j] = min(S_matrix[i][j-1] + GAPCOST + GAP_EXTEND, I_matrix[i][j-1] + GAP_EXTEND)
                S_matrix[i][j] = min(S_matrix[i-1][j-1] + score_matrix[seq1[i-1]][seq2[j-1]],
                                     D_matrix[i][j],
                                     I_matrix[i][j])

    else:
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                score_diagonal = S_matrix[i-1][j-1] + score_matrix[seq1[i-1]][seq2[j-1]]
                score_up = S_matrix[i-1][j] + GAPCOST
                score_left = S_matrix[i][j-1] + GAPCOST
                S_matrix[i][j] = min(score_diagonal, score_left, score_up)

    return S_matrix


def get_optimal_score(seq1_file: TextIO, seq2_file: TextIO, score_matrix: dict) -> int:
    """Takes two seperate .fasta files and converts entries intro a list.
    Only the sequence of the first entry of each file will be used for the alignment."""

    seq1, seq2 = fastaParse(seq1_file), fastaParse(seq2_file)
    filled_matrix = fill_matrix(seq1[0].seq, seq2[0].seq, score_matrix)
    return filled_matrix[-1][-1]


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


def affine_traceback_direction(matrix: list[list], row: int, col: int, match_score: int, k = 1) -> str:
    """Finds the node from which the current node's score comes from -
    used for affine gap score, meaning a pervasive search of the current
    row + column will be done to look for longer gaps."""
    diagonal_score = matrix[row - 1][col - 1]
    up_score = matrix[row - 1][col]
    left_score = matrix[row][col - 1]
    node_score = matrix[row][col]

    while k <= row + 1 and k <= col + 1:
        if node_score == up_score + GAPCOST + k * GAP_EXTEND:
            return ('up', k)
        elif node_score == left_score + GAPCOST + k * GAP_EXTEND:
            return ('left', k)
        else:
            up_score = matrix[row - 1 - k][col]
            left_score = matrix[row][col - 1 - k]
            k += 1
    if node_score == diagonal_score + match_score:
        return ('diagonal', k)


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

        if not args.affine: # Linear gap cost backtrace
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

        if args.affine: # Affine gap cost backtrace
            trace_direction, k = affine_traceback_direction(filled_matrix, row, col, match_score)
            match trace_direction:
                case 'diagonal':
                    align1 = get_base(seq1str, row) + align1
                    align2 = get_base(seq2str, col) + align2
                    row -= 1
                    col -= 1
                case 'up':
                    for _ in range(k):
                        align1 = get_base(seq1str, row) + align1
                        align2 = '-' + align2
                        row -= 1
                case 'left':
                    for _ in range(k):
                        align1 = '-' + align1
                        align2 = get_base(seq2str, col) + align2
                        col -= 1

    return align1, align2


##############################################################
##################### Run the algorithm ######################
##############################################################

def main():
    seq1_file, seq2_file = args.seq1, args.seq2

    aligned = alignment(seq1_file, seq2_file, scoreMatrix)
    print(aligned[0] + '\n' + aligned[1])

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