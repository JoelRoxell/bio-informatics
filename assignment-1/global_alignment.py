import numpy as np
import os
import string
import copy

MATCH_SCORE = 2
GAP_PENALTY = 2
MISMATCH_SCORE = -1

STOP = 0
UP = 1
LEFT = 2
DIAG = 3


def main():
    seq_1 = 'ATCGAT'  # "ATTA"  #
    seq_2 = 'ATACGT'  # "ATTTTA"  #

    seq_1_length = len(seq_1) + 1
    seq_2_length = len(seq_2) + 1

    score_m = np.zeros((seq_1_length, seq_2_length))
    trace_m = np.empty((seq_1_length, seq_2_length), dtype=object)

    optimal_score = 0
    optimal_location = None

    for i in range(1, seq_1_length):
        for j in range(1, seq_2_length):
            score_m[i][0] = score_m[i-1][0] - GAP_PENALTY
            score_m[0][j] = score_m[0][j-1] - GAP_PENALTY

    for i in range(1, seq_1_length):
        for j in range(1, seq_2_length):
            diagonal_score = 0

            if seq_1[i - 1] == seq_2[j-1]:
                diagonal_score = score_m[i-1, j-1] + MATCH_SCORE
            else:
                diagonal_score = score_m[i-1, j-1] + MISMATCH_SCORE

            pre_up_score = score_m[i - 1][j] - GAP_PENALTY
            pre_left_score = score_m[i][j - 1] - GAP_PENALTY

            score = max(diagonal_score,
                        pre_up_score,
                        pre_left_score)

            # TODO: if score is equal more than one dir, store arr
            # print(score, pre_up_score, pre_left_score, diagonal_score)

            dirs = []

            if(score == diagonal_score):
                dirs.append(DIAG)

            if score == pre_up_score:
                dirs.append(UP)

            if score == pre_left_score:
                dirs.append(LEFT)

            trace_m[i][j] = dirs
            score_m[i][j] = score

    print("Score matrix")

    first_row = ' ' * 6

    for char in seq_2:
        first_row += "{:>5}".format(char)

    print(first_row)

    for i in range(0, seq_1_length):
        if (i == 0):
            print(' ', end='')
        else:
            print(seq_1[i - 1], end='')
            # print(trace_m)
        for j in range(0, seq_2_length):
            print('{:>5}'.format(int(score_m[i][j])), end='')

        print()

    print("\nTrace matrix")

    first_row = ' ' * 6

    for char in seq_2:
        first_row += "{:>5}".format(char)

    print(first_row)

    for i in range(0, seq_1_length):
        if (i == 0):
            print(' ', end='')
        else:
            print(seq_1[i - 1], end='')
            # print(trace_m)

        for j in range(0, seq_2_length):
            directions = trace_m[i][j]

            if directions is None:
                pos = ' '
            elif len(directions) == 1:
                pos = directions[0]
            else:
                pos = '|'.join(str(x) for x in directions)
                # print(pos)

            print('{:>5}'.format(pos), end='')

        print()
    print()

    solutions = []

    def trace_back(location, trace, a, b):
        i, j = location
        current = trace[i][j]

        if (current is None):
            print(a, b)
            solutions.append([a, b])

            return

        if len(current) == 1:
            if current[0] == DIAG:
                a += seq_1[i - 1]
                b += seq_2[j - 1]
                i -= 1
                j -= 1
            elif current[0] == LEFT:
                a += '-'
                b += seq_2[j - 1]
                j -= 1
            else:  # up
                a += seq_1[i - 1]
                b += '-'
                i -= 1

            trace_back((i, j), trace, a, b)
        else:
            trace_copy = trace.copy()

            for direction in current:
                trace_copy[i][j] = [direction]
                trace_back(location, trace_copy, a, b)

    trace_back((seq_1_length - 1, seq_2_length - 1), trace_m, "", "")

    return
    # TODO: pretty print solutions

    # Reverse the collected results
    align_1 = align_1[::-1]
    align_2 = align_2[::-1]

    # Display results and matches
    print(align_1)
    matches = 0
    hamming_distance = 0

    for i in range(0, len(align_1)):
        if align_1[i] == align_2[i]:
            print('|', end='')
            matches += 1
        else:
            hamming_distance += 1
            print(' ', end='')
    print()
    print(align_2)

    # Matches / possible matches aka string len
    print("\n\x1b[7mPercent matching for DNA sequences: {0:0.2f}%\x1b[0m".format(
        matches/seq_1_length * 100))
    print("\x1b[48;5;57mHamming distance: {}\x1b[0m\n".format(hamming_distance))


if __name__ == "__main__":
    main()
