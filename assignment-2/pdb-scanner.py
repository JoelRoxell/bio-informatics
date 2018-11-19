import os
import string
import numpy as np
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt


def find_atoms(atom_type, file_dbp_path):
    atom_list = []

    with open(file_dbp_path, 'r') as file:
        while True:
            line = file.readline()

            if not line:
                break

            line = line.split()

            if (line[0] == 'ATOM' and line[2] == 'CA'):
                atom_list.append(line)

    return atom_list


def get_positions(atom_list):
    atoms = []

    for atom in atom_list:
        candidate = atom[5:9]
        atoms.append(candidate)
        # print(candidate)

    return atoms


def distance(a, b):
    a_coord = np.array(a[1:], dtype='float')
    b_coord = np.array(b[1:], dtype='float')
    dist = np.linalg.norm(a_coord - b_coord)

    return dist


def contacts(positions, threshold=12):
    data = []

    for i in range(0, len(positions) - 1):
        current = positions[i]

        for j in range(0, len(positions) - 1):
            b = positions[j]
            dist = distance(current, b)

            if dist <= threshold:
                data.append([int(current[0]), int(b[0])])

    return data


def paint(data, size, split):
    n1 = np.asarray(data)
    x = n1[:, 0]
    y = n1[:, 1]

    line = np.linspace(0, size, num=size)

    print(line)
    split_arr = [split] * size

    plt.plot(x, y, 'bo', markersize=1)

    plt.plot(split_arr, line, 'r')
    plt.plot(line, split_arr, 'r')

    plt.axis([0, size, 0, size])
    plt.grid()
    plt.show()


def main():
    atom_list = find_atoms('CA', os.path.abspath('2CSN.pdb'))
    positions = get_positions(atom_list)
    contact_points = contacts(positions)
    # print(positions)
    split_value = 0
    split_index = 0

    dist = np.asarray(contact_points)

    for CAk in atom_list:
        current = int(CAk[5])
        print(current)
        a = np.where((dist[:, 0] < current) & (dist[:, 1] < current))[0].size
        b = np.where((dist[:, 0] > current) & (dist[:, 1] > current))[0].size
        ext = np.where((dist[:, 0] <= current) &
                       (dist[:, 1] >= current))[0].size

        if not ext:
            score = 0
        else:
            score = (a / ext) * (b / ext)

        if score > split_value:
            split_value = score
            split_index = current

    paint(contact_points, len(positions), split_index)


if __name__ == "__main__":
    main()
