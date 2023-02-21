from typing import Sequence, TextIO
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--affine", help = "changes gap score from linear to affine",
                    action = "store")
args = parser.parse_args()


scoreMatrix = {'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
               'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
               'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
               'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}}
GAPCOST = 5
if args.affine:
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
            matrix[i][0] = i * GAP_EXTEND
        else:
            matrix[i][0] = i * GAPCOST
    for j in range(len(n) + 1):
        if args.affine:
            matrix[0][j] = j * GAP_EXTEND
        else:
            matrix[0][j] = j * GAPCOST
    return matrix


def fill_matrix(seq1: Sequence, seq2: Sequence, score_matrix: dict) -> list[list[int]]:
    """Fills the remaining nodes of the matrix """

    S_matrix = initiate_matrix(seq1, seq2)

    if args.affine:
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


def affine_alignment(seq1_file: TextIO, seq2_file: TextIO, score_matrix: list[list]) -> str:

    # Initial string to save alignment
    align1 = ""
    align2 = ""

    seq1, seq2 = fastaParse(seq1_file), fastaParse(seq2_file)
    seq1str, seq2str = seq1[0].seq, seq2[0].seq
    S_matrix = fill_matrix(seq1str, seq2str, score_matrix)
    print(S_matrix)

    row, col = len(seq1str), len(seq2str)

    while row > 0 and col > 0:
        base1 = seq1str[row - 1]
        base2 = seq2str[col - 1]

        if S_matrix[row][col] == S_matrix[row - 1][col - 1] + score_matrix[base1][base2]:
            align1 = base1 + align1
            align2 = base2 + align2
            row -= 1
            col -= 1

        else:
            k = 1
            while True:
                if row >= k and S_matrix[row][col] == (S_matrix[row - k][col] + (GAPCOST + k * GAP_EXTEND)):
                    align1 = base1 + align1
                    align2 = '-' + align2
                    row -= k
                    break

                elif col >= k and S_matrix[row][col] == (S_matrix[row][col - k] + (GAPCOST + k * GAP_EXTEND)):
                    align1 = '-' + align1
                    align2 = base2 + align2
                    col -= k
                    break

                else:
                    k += 1

    return align1 + "\n" + align2

def alignment(seq1_file: TextIO, seq2_file: TextIO, score_matrix: list[list]) -> str:
    """Creates a possible alignment from two fasta files."""

    if args.affine:
        affine_alignment(seq1_file, seq2_file, score_matrix)

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
        base1 = seq1str[row - 1]
        base2 = seq2str[col - 1]

        match_score = score_matrix[base1][base2]

        trace_direction = traceback_direction(filled_matrix, row, col, match_score)

        match trace_direction:
            case 'up':
                align1 = base1 + align1
                align2 = '-' + align2
                row -= 1
            case 'left':
                align1 = '-' + align1
                align2 = base2 + align2
                col -= 1
            case 'diagonal':
                align1 = base1 + align1
                align2 = base2 + align2
                row -= 1
                col -= 1

    return align1 + "\n" + align2

print(alignment("seq1_test.fasta", "seq2_test.fasta", scoreMatrix))