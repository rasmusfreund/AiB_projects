from typing import Sequence, TextIO
from Bio import SeqIO

scoreMatrix = {'A': {'A': 10, 'C': 2, 'G': 5, 'T': 2},
               'C': {'A': 2, 'C': 10, 'G': 2, 'T': 5},
               'G': {'A': 5, 'C': 2, 'G': 10, 'T': 2},
               'T': {'A': 2, 'C': 5, 'G': 2, 'T': 10}}
GAPCOST = -5


def fastaParse(fasta_file: Sequence) -> list[Sequence]:
    sequenceList = list(SeqIO.parse(fasta_file, "fasta"))
    return sequenceList


def empty_matrix(m: Sequence, n: Sequence) -> list[list]:
    outer_list = []
    for i in range(len(m) + 1):
        inner_list = [None for j in range(len(n) + 1)]
        outer_list.append(inner_list)
    return outer_list


def initiate_matrix(m: Sequence, n: Sequence) -> list[list]:
    matrix = empty_matrix(m, n)
    matrix[0][0] = 0
    for i in range(len(m) + 1):
        matrix[i][0] = i * GAPCOST
    for j in range(len(n) + 1):
        matrix[0][j] = j * GAPCOST
    return matrix


def fill_matrix(seq1: Sequence, seq2: Sequence, score_matrix: dict[str, dict[str, int]]) -> list[list[int]]:
    matrix = initiate_matrix(seq1, seq2)
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_diagonal = matrix[i-1][j-1] + score_matrix[seq1[i-1]][seq2[j-1]]
            score_up = matrix[i-1][j] + GAPCOST
            score_left = matrix[i][j-1] + GAPCOST
            matrix[i][j] = max(score_diagonal, score_left, score_up)
    return matrix


def get_optimal_score(seq1_file: TextIO, seq2_file: TextIO, score_matrix: dict[str, dict[str, int]]) -> int:
    """Takes two seperate .fasta files and converts entries intro a list.
    Only the sequence of the first entry of each file will be used for the alignment."""
    seq1, seq2 = fastaParse(seq1_file), fastaParse(seq2_file)
    filled_matrix = fill_matrix(seq1[0].seq, seq2[0].seq, score_matrix)
    return filled_matrix[-1][-1]


print(get_optimal_score("seq1.fasta", "seq2.fasta", scoreMatrix))