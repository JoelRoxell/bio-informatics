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
    """ global alignment alg. (needleman-wunsch)
    """

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
            pos = direction = int(trace_m[i][j])
            if direction == LEFT:
                pos = 'L'
            elif direction == UP:
                pos = 'U'
            elif direction == DIAG:
                pos = 'D'

            print('{:>5}'.format(pos), end='')

        print()
    print()

    align_1 = ''
    align_2 = ''
    i = seq_1_length - 1
    j = seq_2_length - 1

    while(trace_m[i][j] != STOP):
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

    # print(seq_1)
    # print(seq_2)

    # Display results and matches
    print(align_1)
    matches = 0
    hamming_distance = 0

    for i in range(0, seq_1_length):
        if i >= seq_1_length or i >= seq_2_length:
            break

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
    print("\x1b[48;5;57mHamming distance: {}\x1b[0m".format(hamming_distance))
    print("\x1b[48;5;23mOps. required: {}\x1b[0m\n".format(
        int(score_m[seq_1_length-1][seq_2_length-1])))


if __name__ == "__main__":
    main()
