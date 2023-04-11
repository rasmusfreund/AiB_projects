from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
import argparse
from typing import TextIO

parser = argparse.ArgumentParser(
    description="Compute the Robinson-Fould (RF) distance between two unrooted evolutionary trees."
)

parser.add_argument(
    "tree1", help="First tree to compute RF distance from - must be newick format."
)

parser.add_argument(
    "tree2", help="Second tree to compute RF distance from - must be newick format."
)

parser.add_argument(
    "-d",
    "--draw",
    action="store_true",
    help="Outputs the two trees as ascii art before calculating the RF distance.",
)

args = parser.parse_args()


def read_newick(tree: TextIO):
    return Phylo.read(tree, "newick")


def print_newick(tree1: TextIO, tree2: TextIO):
    """
    Draw the trees as ascii art in the console if requested.
    """
    trees = [next(Phylo.parse(tree1, 'newick')), next(Phylo.parse(tree2, 'newick'))]
    for tree in trees:
        Phylo.draw_ascii(tree)
        print("_" * 70 + "\n")


def bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString("".join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs


def compare(tree1, tree2):
    return len(tree1.symmetric_difference(tree2))


def main():
    tree1, tree2 = args.tree1, args.tree2

    if args.draw:
        print_newick(tree1, tree2)

    distance = compare(bitstrs(read_newick(tree1)), bitstrs(read_newick(tree2)))
    print("Robinson-Foulds distance between the two trees:", distance)

if __name__ == "__main__":
    main()
