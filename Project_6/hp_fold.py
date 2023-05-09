import argparse
import sys
import time
import os

###########################################################
###################     Arguments     #####################
###########################################################

parser = argparse.ArgumentParser(description = """Implementation of the 1/4-approximation algorithm for 2D protein folding.
                                 Accepts a hp-string as input and outputs the relative energy of the fold, measured as the
                                 amount of non-local hydrophobic bonds.""")

parser.add_argument("hp",
                    action = "store",
                    help = "Input string as hp-format")

parser.add_argument("-r", "--runtime",
                    action = "store_true",
                    help = "For runtime analysis")

parser.add_argument("-o", "--output",
                    action = "store_true",
                    help = "Outputs a text file containing hp-strings and their corresponding folded path-string")

args = parser.parse_args()


###########################################################
###################     Functions     #####################
###########################################################

def check_string(hp_string: str):
    legal = ["h", "p"]
    for i in hp_string:
        if i not in legal:
            raise ValueError(f"{i} is not a valid amino acid type")


def match_helper(hp_string: str):
    """Helper function to determine which matching is better
    between even-to-odd and odd-to-even"""
    n = len(hp_string)
    # If string is an even length, the final index has an uneven number
    # and vice versa
    is_even = (n % 2 == 0)

    even_idxs = match_even_hydrophobics(hp_string, n, is_even)
    odd_idxs = match_odd_hydrophobics(hp_string, n, is_even)

    if len(even_idxs) >= len(odd_idxs):
        return even_idxs
    else:
        return odd_idxs


def match_even_hydrophobics(hp_string: str, n: int, is_even: bool) -> list:
    """Matches left side even to right side odd hydrophobic residues"""
    # Every time two indices are matched, we need to decrement the right
    # side of the string.
    # E.g., hphpppphph -> match first and last index:
    # HphpppphpH; the string should now be shortened as to not look at
    # the last index, when matching the second 'h':
    # Hph| pppphp | H; only what is between the '|' should be searched
    # for a matching 'h' from the right side.
    right_limit = 0

    matched_idxs = []

    for i in range(0, n - 1 - right_limit, 2):
        if hp_string[i] != "h":
            continue
        if is_even:
            for j in range(n - 1 - right_limit, i + 2, -2):
                if hp_string[i] == hp_string[j] == "h":
                    matched_idxs.append((i, j))
                    right_limit = n - j + (j % 2)
                    break
        else:
            for j in range(n - 2 - right_limit, i + 2, -2):
                if hp_string[i] == hp_string[j] == "h":
                    matched_idxs.append((i, j))
                    right_limit = n - j
                    break
    return matched_idxs


def match_odd_hydrophobics(hp_string: str, n: int, is_even: bool) -> list:
    """Matches left side odd to right side even hydrophobic residues"""
    # Every time two indices are matched, we need to decrement the right
    # side of the string.
    # E.g., hphpppphph -> match first and last index:
    # HphpppphpH; the string should now be shortened as to not look at
    # the last index, when matching the second 'h':
    # Hph| pppphp | H; only what is between the '|' should be searched
    # for a matching 'h' from the right side.
    right_limit = 0

    matched_idxs = []

    for i in range(1, n - right_limit, 2):
        if hp_string[i] != "h":
            continue
        if is_even:
            for j in range(n - 2 - right_limit, i + 2, -2):
                if hp_string[i] == hp_string[j] == "h":
                    matched_idxs.append((i, j))
                    right_limit = n - j
                    break
        else:
            for j in range(n - 1 - right_limit, i + 2, -2):
                if hp_string[i] == hp_string[j] == "h":
                    matched_idxs.append((i, j))
                    right_limit = n - j + (1 - j % 2)
                    break
    return matched_idxs

def listify(matches: list) -> list:
    match_list = []
    for match in matches:
        for position in match:
            match_list.append(position)
    return sorted(match_list)


def path_fold(direction: str, path: list, idx: int, matches: list):
    """Helper function for the compute_path function"""
    match direction:
        case "n":
            path.append("n" * ((matches[idx + 1] - matches[idx] - 1) // 2))
            path.append("e")
            path.append("s" * ((matches[idx + 1] - matches[idx] - 1) // 2))
            return path
        case "e":
            path.append("e" * ((matches[idx + 1] - matches[idx] + 1) // 2))
            path.append("s")
            path.append("w" * ((matches[idx + 1] - matches[idx]) // 2))
            return path
        case "s":
            path.append("s" * ((matches[idx + 1] - matches[idx] - 1) // 2))
            path.append("w")
            path.append("n" * ((matches[idx + 1] - matches[idx] - 1) // 2))
            return path


def compute_path(hp_string: str, matches:list):
    """Computes a path with absolute directions based on matches found
    by the 1/4-approximation algorithm."""
    path = []
    matches = listify(matches)
    turn = len(matches) // 2

    idx = 0
    i = 0

    while idx < turn - 1:
        if i != matches[idx]:
            path.append("e")
            i += 1
        elif i == matches[idx]:
            fold_length = (matches[idx + 1] - matches[idx]) // 2
            path = path_fold("n", path, idx, matches)
            i += fold_length * 2 - 1
            idx += 1

    fold_length = (matches[idx + 1] - matches[idx]) // 2
    path = path_fold("e", path, idx, matches)
    i += (fold_length + 1) * 2
    idx += 1

    while i < len(hp_string) and idx < (len(matches) - 1):
        if i != matches[idx]:
            path.append("w")
            i += 1
        elif i == matches[idx]:
            fold_length = (matches[idx + 1] - matches[idx]) // 2
            path = path_fold("s", path, idx, matches)
            i += fold_length * 2 - 1
            idx += 1
    remaining = len(hp_string) - matches[-1]
    path.append("w" * remaining)

    return "".join(path)


def output(hp_string, path_string, score, *str_number):
    "Produces a .txt-file containing hp-string, path-string and relative score"

    if not os.path.exists("output.txt"):
        with open("output.txt", "w") as f:
            None
        f.close()

    with open("output.txt", "a") as f:
        if str_number:
            f.write(str(str_number[0][0]) + ":" + " " + hp_string + "\t" + path_string + "\t" + str(score) + "\n")
        else:
            f.write(hp_string + "\t" + path_string + "\t" + str(score) + "\n")
    f.close()


def hp_fold(hp_string, *str_number):
    "Runs the relevant functions"
    if args.runtime:
            st = time.time()

    match_idx = match_helper(hp_string)
    print(match_idx)
    score = len(match_idx)
    path = compute_path(hp_string, match_idx)
    print("Approximate score: " + str(score))
    print("Absolute folding path: " + path)

    if args.runtime:
        et = time.time()
        runtime = et - st
        print(runtime)

    if args.output:
        if "txt" in args.hp[-3:]:
            output(hp_string, path, score, str_number)
        else:
            output(hp_string, path, score)

###########################################################
######################     Main     #######################
###########################################################

def main():

    if "txt" in args.hp[-3:]:
        hp_string = []
        with open(args.hp, "r") as f:
            for line in f:
                try:
                    check_string(line.lower().strip("\n"))
                except:
                    sys.exit("Exiting due to error: " + e.args[0])
                else:
                    hp_string.append(line.lower().strip("\n"))

    else:
        hp_string = args.hp.lower()
        try:
            check_string(hp_string)
        except Exception as e:
            sys.exit("Exiting due to error: " + e.args[0])

    if type(hp_string) is str:
        hp_fold(hp_string)

    if type(hp_string) is list:
        for i in range(len(hp_string)):
            hp_fold(hp_string[i], i)

if __name__ == "__main__":
    main()