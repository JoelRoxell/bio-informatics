def print_matrix(sequence_1, sequence_2, matrix, title):
    print("Matrix: {}".format(title))

    first_row = ' ' * 6

    for char in sequence_2:
        first_row += "{:>5}".format(char)

    print(first_row)

    for i in range(0, len(sequence_1) + 1):
        if (i == 0):
            print(' ', end='')
        else:
            print(sequence_1[i - 1], end='')
            # print(trace_m)

        for j in range(0, len(sequence_2) + 1):
            print('{:>5}'.format(int(matrix[i][j])), end='')

        print()
