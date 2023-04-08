from Bio import Phylo
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


def parse_newick(tree: TextIO):
    return next(Phylo.parse(tree, "newick"))


def print_newick(tree1: TextIO, tree2: TextIO):
    """
    Draw the trees as ascii art in the console if requested.
    """
    trees = [parse_newick(tree1), parse_newick(tree2)]
    for tree in trees:
        Phylo.draw_ascii(tree)
        print("_" * 70 + "\n")


def robinson_foulds_distance(tree1, tree2):
    """
    First parses the two input trees using the parse_newick() function.
    Then, uses the get_subtree() function to obtain a list of all the subtrees in each tree.
    For each subtree in the first tree, the function checks whether it is present in the second tree
    by comparing it to all the subtrees in the second tree using the equal_subtrees() function.
    If a subtree is not present in the second tree, the function increments the Robinson-Foulds distance by one.
    This is then repeated for each subtree in the second tree that is not present in the first tree.
    """
    t1, t2 = parse_newick(tree1), parse_newick(tree2)
    nodes1, nodes2 = [], []
    rfd = 0
    get_subtree(t1.root, nodes1)
    get_subtree(t2.root, nodes2)
    for n1 in nodes1:
        found = False
        for n2 in nodes2:
            if equal_subtrees(n1, n2):
                found = True
                break
        if not found:
            rfd += 1
    for n2 in nodes2:  # Might not need this loop?
        found = False
        for n1 in nodes1:
            if equal_subtrees(n1, n2):
                found = True
                break
        if not found:
            rfd += 1
    return rfd // 2


def equal_subtrees(n1, n2):
    """
    Determine if two subtrees, rooted at nodes n1 and n2, are equal. Returns True if they are equal, and False otherwise.

    The function first checks if n1 and n2 are both strings. If they are, it compares them and returns the result.
    The leaves of the tree are represented by strings, and if n1 and n2 are both strings, they must both be leaves and can be directly compared.

    If n1 and n2 are not both strings, it checks if they are both lists. If they are not both lists, it returns False,
    because the structure of the subtrees does not match. If they are both lists, it checks if they have the same length.
    If they don't have the same length, it returns False, because the subtrees do not have the same structure.

    If n1 and n2 are both lists of the same length, the function recursively calls itself on each pair of corresponding elements in n1 and n2.
    If any of these calls return False, the function immediately returns False, because it means that the corresponding subtrees do not match.
    If all of the calls return True, the function returns True, because it means that the subtrees match at this level.
    """

    if isinstance(n1, str):
        return n1 == n2
    elif isinstance(n1, list):
        if not isinstance(n2, list):
            return False
        if len(n1) != len(n2):
            return False
        for i in range(len(n1)):
            if not equal_subtrees(n1[i], n2[i]):
                return False
        return True
    else:
        return False


def get_subtree(n, nodes):
    """
    Takes a node (n) of a tree and appends all the non-terminal descendants of the node to a nodes list.

    The function first checks if the current node n is a terminal node (i.e., a leaf node).
    If n is not a terminal node, the function appends it to the nodes list, as it is a non-terminal node
    that has descendants that need to be included in the tree comparison.

    Next, the function calls itself recursively for each child of the current node n.
    This will continue until a terminal node is reached, at which point the recursion will stop and the
    function will return to the parent node.

    Overall, the get_subtree() function is used to collect all of the non-terminal nodes of a tree so
    that they can be compared against another tree's non-terminal nodes.
    """

    if not n.is_terminal():
        nodes.append(n)
    for child in n.clades:
        get_subtree(child, nodes)


def main():
    if args.draw:
        print_newick(args.tree1, args.tree2)

    distance = robinson_foulds_distance(args.tree1, args.tree2)
    print("Robinson-Foulds distance between the two trees:", distance)


if __name__ == "__main__":
    main()
