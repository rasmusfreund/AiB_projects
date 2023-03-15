
# Let M be an multiple alignment (represented as a python a list of lists/columns).
#
# Let A be a pairwise alignment (represented as a python a list of lists/columns).
#
# The first row of M and A is the same string if gaps are removed
#
# Example for slide 11 in SP-MSA-Approx.pdf:

# M = [['a','a','c'], ['-','t','t'], ['-','t','-'], ['c','c','c'], ['g','-','g'], ['t','t','a']]

# A = [['a','a'], ['c','c'], ['g','g'], ['-','g'], ['t','t']]

# Let MA be an MSA that is the extension of M with A cf. Gusfields 2-approximation algorithm.
def msa(pairwise_columns: list):

    # Check that all columns have the same length
    # column_length = len(pairwise_columns[0])
    # for column in pairwise_columns[1:]:
    #     if len(column) != column_length:
    #         raise ValueError("All columns must have the same length")


    MA = []
    M = pairwise_columns[0]

    for k in range(1, len(pairwise_columns)):

        i = 0
        j = 0

        A = pairwise_columns[k]

        while i < len(M) and j < len(A):

            # Invariant: (1) MA is a valid merge of all columns before column i in M
            # and all columns before column in A, and (2) the first row of M and A up
            # to (but not including) column i and j respectively is the same string
            # if gaps are removed.

            if M[i][0] == '-' and A[j][0] == '-':
                # Case 1: The next column in MA is column i in M extended with the second symbol
                # in column j in A.
                c = ['-'] * (len(M[i]) + len(A[j]) - 1)
                MA.append(M[i] + A[j][1:] + c)
                i = i + 1
                j = j + 1

            elif M[i][0] == '-' and A[j][0] != '-':
                # Case 2: A[j][0] is a character, so the second symbol in column j in A, A[j][1],
                # must be in the column of MA that is the column in M where the first symbol corresponds
                # to A[j][0]. By the invariant, this column in M is the next column in M, where the first
                # symbol is a character, so we just moved forward in M until we find this column.
                found = False
                for l in range(i + 1, len(M)):
                    if M[l][0] == A[j][0]:
                        M[l].append(A[j][1])
                        MA.append(M[l])
                        i = l
                        found = True
                        break
                if not found:
                    M[i].append('-')
                    MA.append(M[i])
                    i += 1

            elif M[i][0] != '-' and A[j][0] == '-':
                # Case 3: M[i][0] is a character, so column i in M must be in the column of MA that also
                # contains the second symbol from the column in A, where the first symbol is the character
                # corresponding to M[i][0]. By the invariant, this column in A is the next column in A,
                # where the first symbol is a character, so we just add columns from A to MA until we
                # find this column.
                found = False
                for l in range(j + 1, len(A)):
                    if A[l][0] == M[i][0]:
                        c = ['-']*len(M[i])
                        c.append(A[l][1])
                        MA.append(c)
                        j = l
                        found = True
                        break
                if not found:
                    c = ['-'] * len(M[i])
                    c.append('-')
                    MA.append(c)
                    j += 1

            elif M[i][0] != '-' and A[j][0] != '-':
                # Case 4: By the invariant the characters M[i][0] and A[j][0] are at the same position
                # in the string spelled by the row of M and A if gaps are removed. The next column in
                # MA is thus column i in M extended with the second symbol in column j in A.
                if M[i][0] == A[j][0]:
                    M[i].append(A[j][1])
                    MA.append(M[i])
                    i += 1
                    j += 1
                else:
                    found_i, found_j = False, False
                    for l in range(i + 1, len(M)):
                        if M[l][0] == A[j][0]:
                            c = M[l][:j] + ['-'] * (len(M[i]) - j) + [A[j][1]]
                            MA.append(c)
                            i = l
                            found_i = True
                            break
                    if not found_i:
                        c = ['-'] * j + M[i][1:] + [A[j][1]]
                        MA.append(c)
                        i += 1
                    for l in range(j + 1, len(A)):
                        if A[l][0] == M[i][0]:
                            c = A[l][:i] + ['-'] * (len(A[j]) - i) + [M[i][1]]
                            MA.append(c)
                            j = l
                            found_j = True
                            break
                    if not found_j:
                        c = ['-'] * i + A[j][1:]
                        MA.append(c)
                        j += 1

        # add the remaining columns of M to MA
        while i < len(M):
            c = M[i] + ['-'] * (len(A[0]) - 1)
            MA.append(c)
            i += 1

        # add the remaining columns of A to MA
        while j < len(A):
            c = ['-'] * len(M[0])
            c.append(A[j][1])
            MA.append(c)
            j += 1

        M = MA
        MA = []

    return M