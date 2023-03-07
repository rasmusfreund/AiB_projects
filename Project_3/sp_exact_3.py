from typing import Sequence, TextIO
import time
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""Commands for adding command line arguments"""

parser = argparse.ArgumentParser(
    prog="MSA_exact.py",
    usage="python3 %(prog)s sequences.fasta",
    description="Align 3 sequences with an exact match",
)
parser.add_argument(
    "seqs", help="Name of .FASTA file containing the three sequences to be aligned"
)
parser.add_argument(
    "--out",
    help="Produces a .FASTA file containing the aligned sequences as output",
    action="store_true",
)
parser.add_argument(
    "--runtime", help="For runtime analysis purposes", action="store_true"
)
parser.add_argument(
    "--score",
    help="Prints the optimal alignment cost of the alignments",
    action="store_true",
)
args = parser.parse_args()


##############################################################
### Values in scoreMatrix can be changed to any integer,   ###
### however, the alignment will always try to find the     ###
### minimum score, meaning that MATCHES should always      ###
### have the lowest score.                                 ###
##############################################################

scoreMatrix = {
    "A": {"A": 0, "C": 5, "G": 2, "T": 5},
    "C": {"A": 5, "C": 0, "G": 5, "T": 2},
    "G": {"A": 2, "C": 5, "G": 0, "T": 5},
    "T": {"A": 5, "C": 2, "G": 5, "T": 0},
}

GAPCOST = 5

##############################################################
############# Modify below this at your own risk #############
##############################################################


def fastaParse(fasta_file: Sequence) -> list[Sequence]:
    sequenceList = list(SeqIO.parse(fasta_file, "fasta"))
    container = []
    for i in range(len(sequenceList)):
        container.append(sequenceList[i].seq.upper())
    if args.out:
        return container, sequenceList
    return container


def empty_matrix(m: Sequence, n: Sequence, o: Sequence) -> list[list[list]]:
    """Creates a matrix of size len(m + 1) x len(n + 1) x len(o + 1) and fills with None"""

    outer_list = [
        [[None for _ in range(len(o) + 1)] for _ in range(len(n) + 1)]
        for _ in range(len(m) + 1)
    ]
    return outer_list


def fill_matrix(
    m: Sequence, n: Sequence, o: Sequence, score_matrix: dict
) -> list[list[list[int]]]:
    """Fills the nodes of the matrix"""

    S_matrix = empty_matrix(m, n, o)

    for i in range(len(m) + 1):
        for j in range(len(n) + 1):
            for k in range(len(o) + 1):
                # Define all scores as positive infinity in order to exclude them during
                # minimization
                score_diagonal_3d = (
                    score_diagonal_ij
                ) = (
                    score_diagonal_ik
                ) = (
                    score_diagonal_jk
                ) = score_prev_i = score_prev_j = score_prev_k = float("inf")

                if (i, j, k) == (0, 0, 0):
                    S_matrix[0][0][0] = 0
                    continue

                if i > 0 and j > 0 and k > 0:
                    score_diagonal_3d = (
                        S_matrix[i - 1][j - 1][k - 1]
                        + score_matrix[m[i - 1]][n[j - 1]]
                        + score_matrix[m[i - 1]][o[k - 1]]
                        + score_matrix[n[j - 1]][o[k - 1]]
                    )

                if i > 0 and j > 0:
                    score_diagonal_ij = (
                        S_matrix[i - 1][j - 1][k]
                        + score_matrix[m[i - 1]][n[j - 1]]
                        + GAPCOST * 2
                    )

                if i > 0 and k > 0:
                    score_diagonal_ik = (
                        S_matrix[i - 1][j][k - 1]
                        + score_matrix[m[i - 1]][o[k - 1]]
                        + GAPCOST * 2
                    )

                if j > 0 and k > 0:
                    score_diagonal_jk = (
                        S_matrix[i][j - 1][k - 1]
                        + score_matrix[n[j - 1]][o[k - 1]]
                        + GAPCOST * 2
                    )

                if i > 0:
                    score_prev_i = S_matrix[i - 1][j][k] + GAPCOST * 2

                if j > 0:
                    score_prev_j = S_matrix[i][j - 1][k] + GAPCOST * 2

                if k > 0:
                    score_prev_k = S_matrix[i][j][k - 1] + GAPCOST * 2

                S_matrix[i][j][k] = min(
                    score_diagonal_3d,
                    score_diagonal_ij,
                    score_diagonal_ik,
                    score_diagonal_jk,
                    score_prev_i,
                    score_prev_j,
                    score_prev_k,
                )

    if args.score:  # Get optimal alignment score if requested
        print(S_matrix[-1][-1][-1])

    return S_matrix


def traceback_direction(
    matrix: list[list[list]],
    row: int,
    col: int,
    depth: int,
    rowBase: str,
    colBase: str,
    depthBase: str,
    score_matrix: dict,
) -> str:
    """Finds the node from which the current node's score comes from."""

    node_score = matrix[row][col][depth]

    up2D_score = matrix[row - 1][col][depth]
    left2D_score = matrix[row][col - 1][depth]
    depth2D_score = matrix[row][col][depth - 1]
    diagonal_isoDepth_score = matrix[row - 1][col - 1][depth]
    diagonal_isoCol_score = matrix[row - 1][col][depth - 1]
    diagonal_isoRow_score = matrix[row][col - 1][depth - 1]
    diagonal3D_score = matrix[row - 1][col - 1][depth - 1]

    match_score3D = (
        score_matrix[rowBase][colBase]
        + score_matrix[rowBase][depthBase]
        + score_matrix[colBase][depthBase]
    )
    match_score_isoDepth = score_matrix[rowBase][colBase]
    match_score_isoCol = score_matrix[rowBase][depthBase]
    match_score_isoRow = score_matrix[colBase][depthBase]

    if node_score == diagonal3D_score + match_score3D:
        return "diagonal3D"
    elif node_score == diagonal_isoDepth_score + match_score_isoDepth + GAPCOST * 2:
        return "diagonal same depth"
    elif node_score == diagonal_isoCol_score + match_score_isoCol + GAPCOST * 2:
        return "diagonal same column"
    elif node_score == diagonal_isoRow_score + match_score_isoRow + GAPCOST * 2:
        return "diagonal same row"
    elif node_score == left2D_score + GAPCOST * 2:
        return "left"
    elif node_score == up2D_score + GAPCOST * 2:
        return "up"
    elif node_score == depth2D_score + GAPCOST * 2:
        return "depth"


def get_base(sequence: Sequence, position: int) -> str:
    """Simply grabs a single nucleotide from a sequence based on a positional argument"""
    return sequence[position - 1]


def alignment(
    seq1: Sequence, seq2: Sequence, seq3: Sequence, score_matrix: dict
) -> str:
    """Creates a possible alignment from two fasta files."""

    # Initial string to save alignment
    align1, align2, align3 = "", "", ""

    # Idx for bottom right node in the matrix
    row = len(seq1)
    col = len(seq2)
    depth = len(seq3)

    # Matrix containing alignment scores
    S_matrix = fill_matrix(seq1, seq2, seq3, score_matrix)

    # Backtrack and create alignment
    while row > 0 and col > 0 and depth > 0:
        rowBase, colBase, depthBase = (
            get_base(seq1, row),
            get_base(seq2, col),
            get_base(seq3, depth),
        )
        trace_direction = traceback_direction(
            S_matrix, row, col, depth, rowBase, colBase, depthBase, scoreMatrix
        )

        match trace_direction:
            case "diagonal3D":
                align1 = get_base(seq1, row) + align1
                align2 = get_base(seq2, col) + align2
                align3 = get_base(seq3, depth) + align3
                row, col, depth = row - 1, col - 1, depth - 1
            case "diagonal same depth":
                align1 = get_base(seq1, row) + align1
                align2 = get_base(seq2, col) + align2
                align3 = "-" + align3
                row, col = row - 1, col - 1
            case "diagonal same column":
                align1 = get_base(seq1, row) + align1
                align2 = "-" + align2
                align3 = get_base(seq3, depth) + align3
                row, depth = row - 1, depth - 1
            case "diagonal same row":
                align1 = "-" + align1
                align2 = get_base(seq2, col) + align2
                align3 = get_base(seq3, depth) + align3
                col, depth = col - 1, depth - 1
            case "up":
                align1 = get_base(seq1, row) + align1
                align2 = "-" + align2
                align3 = "-" + align3
                row -= 1
            case "left":
                align1 = "-" + align1
                align2 = get_base(seq2, col) + align2
                align3 = "-" + align3
                col -= 1
            case "depth":
                align1 = "-" + align1
                align2 = "-" + align2
                align3 = get_base(seq3, depth) + align3
                depth -= 1

    return align1, align2, align3


def create_output(aligned_sequences, sequence_list):
    """Produces an output-file in .fasta format using the original information
    regarding the sequences."""
    # Convert alignments to sequence objects
    seqs = [Seq(x) for x in aligned_sequences]

    # Convert sequence objects to sequence records
    recordList = []
    for i in range(3):
        recordList.append(
            SeqRecord(
                seqs[i],
                id=sequence_list[i].id,
                name=sequence_list[i].name,
                description=sequence_list[i].description,
            )
        )

    SeqIO.write(recordList, "aligned_sequences.fasta", "fasta")


##############################################################
##################### Run the algorithm ######################
##############################################################


def main():
    if args.runtime:
        st = time.time()

    # Parse sequences
    if args.out:
        parsed_seqs, sequenceList = fastaParse(args.seqs)

    if not args.out:
        parsed_seqs = fastaParse(args.seqs)

    seq1, seq2, seq3 = parsed_seqs[0], parsed_seqs[1], parsed_seqs[2]

    # Create alignment
    aligned = alignment(seq1, seq2, seq3, scoreMatrix)
    print(aligned[0] + "\n" + aligned[1] + "\n" + aligned[2])

    # Output runtime if requested
    if args.runtime:
        et = time.time()
        elapsed_time = et - st
        print("Execution time:", elapsed_time, "seconds")

    # Output FASTA file with aligned sequences if requested
    if args.out:
        create_output(aligned, sequenceList)


if __name__ == "__main__":
    main()
