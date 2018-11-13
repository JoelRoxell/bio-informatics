import numpy as np
import os

from utils.print_m import print_matrix

MATCH_SCORE = 2
GAP_PENALTY = 2
MISMATCH_SCORE = 1

STOP = 0
UP = 1
LEFT = 2
DIAG = 3


def main():
    seq_1 = 'PAWHEAE'
    seq_2 = 'HDAGAWGHEQ'

    seq_1_length = len(seq_1) + 1
    seq_2_length = len(seq_2) + 1
    largest = max(seq_1_length, seq_2_length)

    score_m = np.zeros((largest, largest))
    trace_m = np.zeros((largest, largest))

    optimal_score = 0
    optimal_location = None

    for i in range(1, seq_1_length):
        for j in range(1, seq_2_length):
            diagonal_score = 0
            direction = DIAG

            if seq_1[i - 1] == seq_2[j-1]:
                diagonal_score = score_m[i-1, j-1] + MATCH_SCORE
            else:
                diagonal_score = score_m[i-1, j-1] - MISMATCH_SCORE

            pre_up_score = score_m[i - 1][j] - GAP_PENALTY
            pre_left_score = score_m[i][j - 1] - GAP_PENALTY

            score = max(diagonal_score,
                        pre_up_score,
                        pre_left_score)

            if score == pre_up_score:
                direction = UP
            elif score == pre_left_score:
                direction = LEFT

            trace_m[i][j] = direction
            score_m[i][j] = max(score, 0)

            if score > optimal_score:
                optimal_score = score
                optimal_location = (i, j)

    print_matrix(seq_1, seq_2, score_m, "score")
    print_matrix(seq_1, seq_2, trace_m, "trace")

    align_1 = ''
    align_2 = ''
    i, j = optimal_location

    while(score_m[i][j] != 0):
        current = trace_m[i][j]

        if current == DIAG:
            align_1 += seq_1[i - 1]
            align_2 += seq_2[j - 1]
            i -= 1
            j -= 1
        elif current == LEFT:
            align_1 += '-'
            align_2 += seq_2[j - 1]
            j -= 1
        else:  # up
            align_1 += seq_1[i - 1]
            align_2 += '-'
            i -= 1

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

    print("\noptimal score location: {}".format(optimal_location))
    print("optimal score: {}".format(optimal_score))

    # Matches / possible matches aka string len
    print("\n\x1b[7mPercent matching for DNA sequences: {0:0.2f}%\x1b[0m".format(
        matches/max(seq_1_length, seq_2_length) * 100))
    print("\x1b[48;5;57mHamming distance: {}\x1b[0m\n".format(hamming_distance))


if __name__ == "__main__":
    main()
