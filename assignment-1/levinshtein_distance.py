import numpy as np
import os

MATCH_SCORE = 2
GAP_PENALTY = 2
MISMATCH_SCORE = -1
PENALTY = 1

STOP = 0
UP = 1
LEFT = 2
DIAG = 3


def main():
    seq_1 = 'azced'
    seq_2 = 'abcdef'

    seq_1_length = len(seq_1) + 1
    seq_2_length = len(seq_2) + 1
    largest = max(seq_1_length, seq_2_length)

    score_m = np.zeros((largest, largest))
    trace_m = np.zeros((largest, largest))

    for i in range(1, seq_1_length):
        for j in range(1, seq_2_length):
            score_m[i][0] = i
            score_m[0][j] = j

    for i in range(1, seq_1_length):
        for j in range(1, seq_2_length):
            diagonal_score = 0
            direction = DIAG

            if seq_1[i - 1] == seq_2[j-1]:
                diagonal_score = score_m[i-1, j-1]
            else:
                diagonal_score = score_m[i-1, j-1] + PENALTY

            pre_up_score = score_m[i - 1][j] + PENALTY
            pre_left_score = score_m[i][j - 1] + PENALTY

            score = min(diagonal_score,
                        pre_up_score,
                        pre_left_score)

            if score == pre_up_score:
                direction = UP
            elif score == pre_left_score:
                direction = LEFT

            trace_m[i][j] = direction
            score_m[i][j] = max(score, 0)

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

    print("\n\x1b[48;5;23mOperations required: {}\x1b[0m\n".format(
        int(score_m[seq_1_length-1][seq_2_length-1])))


if __name__ == "__main__":
    main()
