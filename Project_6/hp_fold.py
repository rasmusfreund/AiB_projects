import argparse
import sys

###########################################################
###################     Arguments     #####################
###########################################################

parser = argparse.ArgumentParser(description = """Implementation of the 1/4-approximation algorithm for 2D protein folding.
                                 Accepts a hp-string as input and outputs the relative energy of the fold, measured as the
                                 amount of non-local hydrophobic bonds.""")

parser.add_argument("hp", 
                    action = "store",
                    help = "Input string as hp-format")

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


###########################################################
######################     Main     #######################
###########################################################

def main():
    hp_string = args.hp.lower()

    try:
        check_string(hp_string)
    except Exception as e:
        sys.exit("Exiting due to error: " + e.args[0])

    print(hp_string)

    match_idx = match_helper(hp_string)
    print(match_idx)


if __name__ == "__main__":
    main()