import argparse

parser = argparse.ArgumentParser(description = """Implementation of the 1/4-approximation algorithm for 2D protein folding.
                                 Accepts a hp-string as input and outputs the relative energy of the fold.""")

parser.add_argument("hp", 
                    action = "store",
                    help = "Input string as hp-format")

args = parser.parse_args()


def string_parse_helper(hp_string: str):
    n = len(hp_string)
    is_even = 1 - (n % 2)

    idx_even = string_parse_even(hp_string, is_even)
    idx_odd = string_parse_odd(hp_string, is_even)

    print(idx_even)
    print(idx_odd)

    if len(idx_even) > len(idx_odd):
        return idx_even
    if len(idx_even) < len(idx_odd):
        return idx_odd


def string_parse_even(hp_string, is_even):
    idx_even = []
    decrement_length = 0

    for i in range(0, len(hp_string), 2):
        if hp_string[i] != "h":
            continue
        match is_even:
            case 0:
                for j in range(len(hp_string) - 2 - decrement_length, i, -2):
                    if hp_string[j] != "h":
                        continue
                    elif -2 < j - i < 2:
                        continue
                    else:
                        idx_even.append((i, j))
                        decrement_length = (len(hp_string) - j + 1)
                        break
            case 1:
                for j in range(len(hp_string) - 1 - decrement_length, i, -2):
                    if hp_string[j] != "h":
                        continue
                    elif -2 < j - i < 2:
                        continue
                    else:
                        idx_even.append((i, j))
                        decrement_length = (len(hp_string) - j + 1)
                    break
    return idx_even

def string_parse_odd(hp_string, is_even):
    idx_odd = []
    decrement_length = 0

    for i in range(1, len(hp_string), 2):
        if hp_string[i] != "h":
            continue
        match is_even:
            case 0:
                for j in range(len(hp_string) - 1 - decrement_length, i, -2):
                    if hp_string[j] != "h":
                        continue
                    elif -2 < j - i < 2:
                        continue
                    else:
                        idx_odd.append((i, j))
                        decrement_length = (len(hp_string) - j + 1)
                        break
            case 1:
                for j in range(len(hp_string) - 2 - decrement_length, i, -2):
                    if hp_string[j] != "h":
                        continue
                    elif -2 < j - i < 2:
                        continue
                    else:
                        idx_odd.append((i, j))
                        decrement_length = (len(hp_string) - j + 1)
                        break
    return idx_odd


def main():
    hp_string = args.hp.lower()
    print(hp_string)
    match_idx = string_parse_helper(hp_string)
    print(match_idx)


if __name__ == "__main__":
    main()