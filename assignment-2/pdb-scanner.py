import os
import string
import numpy as np
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt
from itertools import cycle


cycol = cycle('bgrcmk')


def find_atoms(atom_type, file_dbp_path):
    """Finds all CA atoms.

    Arguments:
        atom_type {[type]} -- [description]
        file_dbp_path {[type]} -- [description]

    Returns:
        [type] -- [description]
    """

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
    """Extracts the x,y,z position for each atom.

    Arguments:
        atom_list {[type]} -- [description]

    Returns:
        [type] -- [description]
    """

    atoms = []

    for atom in atom_list:
        candidate = atom[5:9]
        atoms.append(candidate)
        # print(candidate)

    return atoms


def distance(a, b):
    """Returns distance for two points.

    Arguments:
        a {[type]} -- [description]
        b {[type]} -- [description]

    Returns:
        [type] -- [description]
    """

    a_coord = np.array(a[1:], dtype='float')
    b_coord = np.array(b[1:], dtype='float')
    dist = np.linalg.norm(a_coord - b_coord)

    return dist


def contacts(positions, threshold=12):
    """Find contacts between points with a specified threshold.

    Arguments:
        positions {[type]} -- [description]

    Keyword Arguments:
        threshold {int} -- [description] (default: {12})

    Returns:
        [type] -- [description]
    """

    data = []

    for i in range(0, len(positions) - 1):
        current = positions[i]

        for j in range(0, len(positions) - 1):
            b = positions[j]
            dist = distance(current, b)

            if dist <= threshold:
                data.append([int(current[0]), int(b[0])])

    return data


def paint(data, size, splits):
    n1 = np.asarray(data)
    x = n1[:, 0]
    y = n1[:, 1]

    line = np.linspace(0, size, num=size)

    # print(line)
    plt.plot(x, y, 'bo', markersize=1)

    for split in splits:

        color = next(cycol)

        split_arr = [split] * size
        plt.plot(split_arr, line, c=color)
        plt.plot(line, split_arr, c=color)

    plt.axis([0, size, 0, size])
    plt.grid()
    plt.show()


def find_split(positions, atom_list):
    contact_points = contacts(positions)

    # print(positions)
    split_value = 0
    split_residue_index = 0

    dist = np.asarray(contact_points)

    i = 0

    for CAk in atom_list:
        current = int(CAk[5])

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
            split_residue_index = current
            index = i
            # print(score, split_value, current)

        i += 1

    return index, split_value, split_residue_index


def find_rec(positions, atom_list):

    positions_inner = get_positions(atom_list)

    i, split_value, split_residue_index_i = find_split(
        positions_inner, atom_list)

    domain_a = atom_list[:i]
    domain_b = atom_list[i:]

    if (len(domain_a) < MIN_RESIDUE or len(domain_b) < MIN_RESIDUE):
        return

    if split_value <= MIN_SPLIT_SCORE:
        return

    splits.append(split_residue_index_i)

    positions_a = get_positions(domain_a)
    positions_b = get_positions(domain_b)

    find_rec(positions_a, domain_a)
    find_rec(positions_b, domain_b)


splits = []
MIN_RESIDUE = 40
MIN_SPLIT_SCORE = 10


def main():
    atom_list = find_atoms('CA', os.path.abspath('4GAF_B.pdb'))  # 1HZH_H.pdb
    contact_points = contacts(get_positions(atom_list))
    positions = get_positions(atom_list)

    find_rec(positions, atom_list)

    paint(contact_points, len(positions), splits)


if __name__ == "__main__":
    main()
