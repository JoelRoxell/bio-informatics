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
    seq_1 = "ATTA"  # "ATCGAT"  #
    seq_2 = "ATTTTA"  # "ATACGT"  #

    # Add one index for staring cells.
    seq_1_length = len(seq_1) + 1
    seq_2_length = len(seq_2) + 1

    # Initialize score matrix and trace.
    score_m = np.zeros((seq_1_length, seq_2_length))
    trace_m = np.empty((seq_1_length, seq_2_length), dtype=object)

    for i in range(1, seq_1_length):
        for j in range(1, seq_2_length):
            score_m[i][0] = score_m[i-1][0] - GAP_PENALTY
            score_m[0][j] = score_m[0][j-1] - GAP_PENALTY

    # Calculate scores for each cell.
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
        for j in range(0, seq_2_length):
            print('{:>5}'.format(int(score_m[i][j])), end='')

        print()
    print()

    solutions = []

    def trace_back(location, trace, a, b):
        i, j = location
        current = trace[i][j]

        if (current is None):
            # End of search, add resulting paths to solutions & exit.
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

    print("Optimal solutions: {}".format(len(solutions)))

    for solution in solutions:
        align_1 = solution[0][::-1]
        align_2 = solution[1][::-1]
        matches = 0
        hamming_distance = 0

        print(align_1)

        for i in range(0, len(align_1)):
            if align_1[i] == align_2[i]:
                print('|', end='')

                matches += 1
            else:
                hamming_distance += 1

                print(' ', end='')

        print()
        print(align_2)

    # Percent identity is caluclated by using # of matches divided by the largest sequence length.
    # The argument beeing that the longest sequence defines the maximum number of possible matches that exists
    # in the given scenario.
    match_p = matches/max(len(seq_1), len(seq_2)) * 100

    print(
        "\n\x1b[7mPercent matching for DNA sequences: {0:0.2f}%\x1b[0m".format(match_p))
    print("\x1b[48;5;57mHamming distance: {}\x1b[0m\n".format(hamming_distance))


if __name__ == "__main__":
    main()
