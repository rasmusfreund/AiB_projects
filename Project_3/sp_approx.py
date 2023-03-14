from typing import Sequence, TextIO
import time
import argparse
import os.path
import itertools
from extend_msa import msa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


"""Commands for adding command line arguments"""

parser = argparse.ArgumentParser(
    prog="sp_exact_3.py",
    usage="python3 %(prog)s sequences.fasta",
    description="Align 3 sequences with an exact match",
)
parser.add_argument(
    "seqs", help="Name of .FASTA file containing the sequences to be aligned"
)
parser.add_argument(
    "nr_seqs", help="""Indicator for how many sequences in the FASTA file should be aligned.
    Default is all sequences, which can also be chosen by entering 0 as an argument.""",
    type = int, action = "store", default = 0, nargs = '?'
)
parser.add_argument(
    "--out",
    help="Produces a .FASTA file containing the aligned sequences as output",
    action="store_true",
)
parser.add_argument(
    "--runtime",
    help="""For runtime analysis purposes; if not present in the working directory,
    this creates a .txt-file (runtime.txt) containing two tab-separated columns:
    'Sequence_length' and 'Runtime'. Once the file exists, additional runtime calls
    will simply append the values of the two columns to the existing file.""",
    action="store_true",
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
    if args.nr_seqs == 0:
        for i in range(len(sequenceList)):
            container.append(sequenceList[i].seq.upper())
        if args.out:
            return container, sequenceList
        return container
    else:
        nr_seqs = args.nr_seqs
        for i in range(nr_seqs):
            container.append(sequenceList[i].seq.upper())
        if args.out:
            return container, sequenceList
        return container


def empty_matrix(m: Sequence, n: Sequence) -> list[list]:
    """Creates a matrix of size len(m) x len(n) and fills with None"""

    outer_list = [[None for _ in range(len(m) + 1)] for _ in range(len(n) + 1)]
    return outer_list


def initiate_matrix(m: Sequence, n: Sequence) -> list[list]:
    """Fills out first row and first column
    of the matrix using the gapcost."""

    matrix = empty_matrix(m, n)
    matrix[0][0] = 0
    for i in range(1, len(m) + 1):
        matrix[i][0] = i * GAPCOST
    for j in range(1, len(n) + 1):
        matrix[0][j] = j * GAPCOST

    return matrix


def fill_matrix(seq1: Sequence, seq2: Sequence, score_matrix: dict) -> list[list[int]]:
    """Fills the remaining nodes of the matrix"""

    S_matrix = initiate_matrix(seq1, seq2)

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_diagonal = (
                S_matrix[i - 1][j - 1] + score_matrix[seq1[i - 1]][seq2[j - 1]]
            )
            score_up = S_matrix[i - 1][j] + GAPCOST
            score_left = S_matrix[i][j - 1] + GAPCOST
            S_matrix[i][j] = min(score_diagonal, score_left, score_up)

    return S_matrix


def traceback_direction(
    matrix: list[list], row: int, col: int, match_score: int
) -> str:
    """Finds the node from which the current node's score comes from."""

    diagonal_score = matrix[row - 1][col - 1]
    up_score = matrix[row - 1][col]
    left_score = matrix[row][col - 1]
    node_score = matrix[row][col]

    if node_score == diagonal_score + match_score:
        return "diagonal"
    elif node_score == left_score + GAPCOST:
        return "left"
    elif node_score == up_score + GAPCOST:
        return "up"


def get_base(sequence: Sequence, position: int) -> str:
    """Simply grabs a single nucleotide from a sequence based on a positional argument"""
    return sequence[position - 1]


def alignment(seq1: TextIO, seq2: TextIO, score_matrix: list[list]) -> str:
    """Creates a possible alignment from two fasta files."""

    # Initial string to save alignment
    align1 = ""
    align2 = ""

    # Load sequences and fill out alignment scores
    filled_matrix = fill_matrix(seq1, seq2, score_matrix)

    # Idx for bottom right node in the matrix
    row = len(seq1)
    col = len(seq2)

    # Backtrack and create alignment
    while row > 0 and col > 0:
        match_score = score_matrix[get_base(seq1, row)][get_base(seq2, col)]
        trace_direction = traceback_direction(filled_matrix, row, col, match_score)

        match trace_direction:
            case "diagonal":
                align1 = get_base(seq1, row) + align1
                align2 = get_base(seq2, col) + align2
                row -= 1
                col -= 1
            case "up":
                align1 = get_base(seq1, row) + align1
                align2 = "-" + align2
                row -= 1
            case "left":
                align1 = "-" + align1
                align2 = get_base(seq2, col) + align2
                col -= 1

    return align1, align2


def seq_combinations(nr_seqs: int) -> list:
    combos = itertools.combinations(range(nr_seqs), 2)
    return [x for x in combos]


def alignment_dict(nr_seqs: int) -> dict:
    alignment_scores = {x: [] for x in range(nr_seqs)}
    return alignment_scores


def center_seq(seq_combinations: list, seq_list: list, align_dict: dict, score_matrix: list[list]) -> dict:
    min_score = (float('inf'), float('inf'))
    for combo in seq_combinations:
        i, j = combo
        filled_matrix = fill_matrix(seq_list[i], seq_list[j], score_matrix)
        align_dict[i].append(filled_matrix[-1][-1])

    for m in range(1, len(align_dict)):
        for n in range(m):
            align_dict[m].insert(0, align_dict[n][m-1])

    for k in range(len(align_dict)):
        if sum(align_dict[k]) < min_score[1]:
            min_score = (k, sum(align_dict[k]))

    return min_score


def center_combinations(center_seq: int, nr_seqs: int) -> list:
    combinations = []
    for i in range(nr_seqs):
        if i == center_seq:
            continue
        combinations.append((center_seq, i))

    return combinations


def convert_to_columns(alignments: list) -> list:
    columns = []
    for i in range(len(alignments)):
        columns.append([])
        #for j in range(len(alignments[0])):
            #columns[i].append([alignments[0][j], alignments[1][j]])
    return columns




def create_output(aligned_sequences, sequence_list):
    """Produces an output-file in .fasta format using the original information
    regarding the sequences."""
    # Convert alignments to sequence objects
    seqs = [Seq(x) for x in aligned_sequences]

    # Convert sequence objects to sequence records
    recordList = []
    for i in range(2):
        recordList.append(
            SeqRecord(
                seqs[i],
                id=sequence_list[i].id,
                name=sequence_list[i].name,
                description=sequence_list[i].description,
            )
        )

    SeqIO.write(recordList, "aligned_sequences.fasta", "fasta")


def output_runtime(st, et, seq1, seq2):
    "Produces a .txt-file containing runtime and sequence length"
    elapsed_time = et - st
    print("Execution time:", elapsed_time, "seconds")

    if not os.path.exists("runtime.txt"):
        with open("runtime.txt", "w") as f:
            f.write("Sequence_length" + "\t" + "Runtime\n")
        f.close()

    with open("runtime.txt", "a") as f:
        if len(seq1) == len(seq2):
            f.write(str(len(seq1)) + "\t" + str(elapsed_time) + "\n")
        else:
            long_seq = max(len(seq1), len(seq2))
            f.write(str(long_seq) + "\t" + str(elapsed_time) + "\n")
    f.close()


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

    # Sequence combinations
    if args.nr_seqs == 0:
        combinations = seq_combinations(len(parsed_seqs))
        align_dict = alignment_dict(len(parsed_seqs))
        min_seq, min_score = center_seq(combinations, parsed_seqs, align_dict, scoreMatrix)
        center_combos = center_combinations(min_seq, len(parsed_seqs))
    if args.nr_seqs != 0:
        combinations = seq_combinations(args.nr_seqs)
        align_dict = alignment_dict(args.nr_seqs)
        min_seq, min_score = center_seq(combinations, parsed_seqs, align_dict, scoreMatrix)
        center_combos = center_combinations(min_seq, args.nr_seqs)

    if args.score:
        print("Alignment of sequence", min_seq + 1, "to all other sequences, resulted in a sum of", min_score)

    # Get all alignments
    alignments = []
    for i, j in center_combos:
        align1, align2 = alignment(parsed_seqs[i], parsed_seqs[j], scoreMatrix)
        alignments.append([align1, align2])

    column_pairs = convert_to_columns(alignments)
    print(column_pairs)

    # Output runtime if requested
    if args.runtime:
        et = time.time()
        output_runtime(st, et, seq1, seq2)

    # Output FASTA file with aligned sequences if requested
    if args.out:
        create_output(aligned, sequenceList)


if __name__ == "__main__":
    main()
