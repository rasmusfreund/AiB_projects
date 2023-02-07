from typing import Sequence
from Bio import SeqIO

scoreMatrix = {'A': {'A': 10, 'C': 2, 'G': 5, 'T': 2},
               'C': {'A': 2, 'C': 10, 'G': 2, 'T': 5},
               'G': {'A': 5, 'C': 2, 'G': 10, 'T': 2},
               'T': {'A': 2, 'C': 5, 'G': 2, 'T': 10}}
GAPCOST = -5

def fastaParse(fasta_file: Sequence) -> list[Sequence]:
    sequenceList = list(SeqIO.parse(fasta_file, "fasta"))
    return sequenceList


def empty_matrix(m, n):
    outer_list = []
    for i in range(len(m)):
        inner_list = [None for j in range(len(n))]
        outer_list.append(inner_list)
    return outer_list


def initiate_matrix(m, n):
    matrix = empty_matrix(m, n)
    matrix[0][0] = 0
    for i in range(len(m)):
        matrix[i][0] = i * GAPCOST
    for j in range(len(n)):
        matrix[0][j] = j * GAPCOST
    return matrix





seq1, seq2 = fastaParse("seq1.fasta"), fastaParse("seq2.fasta")

print(initiate_matrix(seq1[0].seq, seq2[0].seq))