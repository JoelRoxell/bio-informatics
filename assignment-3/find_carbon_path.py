import os
import string
import numpy as np


def distance_between(atom_1, atom_2):
    current_atom = np.array(atom_1, dtype='float')
    next_atom = np.array(atom_2, dtype='float')
    distance = np.linalg.norm(current_atom - next_atom)

    return distance


def parse_coordinates(line):
    ID, coordinates = line.split('\t')
    x, y, z = coordinates.split(' ')

    return ID.strip(), (float(x), float(y), float(z))


def read_file(file_path, data=[]):
    file = open(file_path)

    while True:
        line = file.readline()

        if not line:
            break

        data.append(parse_coordinates(line))

    return data


def determine_path(atom, pre, atom_list, path=''):
    ID, coordinates = atom

    best_path = ''

    for next_atom in atom_list:
        next_ID, next_coordinates = next_atom

        if ID == next_ID:
            continue

        if pre is not None:
            if pre[0] == next_ID:
                continue

        distance = distance_between(coordinates, next_coordinates)

        if distance <= 3.9 and distance > 0:
            sub_path = determine_path(next_atom, atom, atom_list, path)

            if len(best_path) <= len(sub_path):
                best_path = sub_path

    return "{},{}".format(ID, best_path)


def main():
    atom_list = read_file('data_q1.txt')
    best_path = ''

    for atom in atom_list:
        ID = atom[0]

        test_path = determine_path(atom, None, atom_list, path=ID)

        print(test_path)

        if len(test_path) >= len(best_path):
            best_path = test_path

    print(best_path)


if __name__ == '__main__':
    main()
