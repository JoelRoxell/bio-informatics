import math

# Joel Roxell


def parse_file_to_atoms(file_path):
    """Parse db-file to an atom array (atom coordinates).
    """

    atoms = []

    with open(file_path, 'r') as file:
        while True:
            line = file.readline()

            if not line:
                break

            line = line.split()

            if (line[0] == 'ATOM' or line[0] == 'HETATM'):
                atoms.append(line)

    return atoms


def to_float(x):
    return float(x)


def distance_between(atom_1, atom_2):
    x1, y1, z1 = map(to_float, atom_1[6:9])
    x2, y2, z2 = map(to_float, atom_2[6:9])

    return math.sqrt(((x2-x1) ** 2) + ((y2-y1) ** 2) + ((z2-z1) ** 2))


def find_overlaps(atom_arr_1, atom_arr_2, atom_radius):
    """ Finds all overlapping atoms in both lists; in consideration to the specified atom radius.
    """

    comparison_count = 0
    collision_count = 0
    overlaps = {}

    for atom_1 in atom_arr_1:
        for atom_2 in atom_arr_2:
            if int(atom_2[1]) in overlaps:
                continue

            # All atoms are assumed to have the same radius and if the distance of between two points
            # is less than the comprised radius of the two, that is, 2 * atom_radius; such a case determines an
            # overlap.
            if distance_between(atom_1, atom_2) < 2 * atom_radius:
                overlaps.update({int(atom_2[1]): atom_2[1:6]})
                collision_count += 1

            comparison_count += 1

    return comparison_count, collision_count, overlaps


def main():
    dataset_1 = 'assets/1CDH.pdb'
    dataset_2 = 'assets/2CSN.pdb'
    ATOM_RADIUS = 2

    atom_arr_1 = parse_file_to_atoms(dataset_1)
    atom_arr_2 = parse_file_to_atoms(dataset_2)

    comparison_count, collision_count, overlaps = find_overlaps(
        atom_arr_1,
        atom_arr_2,
        ATOM_RADIUS
    )

    for overlap in sorted(overlaps):
        print(overlaps.get(overlap))

    print('Number of overlapping atoms: {}\nNumber of comparisons: {}'.format(
        collision_count, comparison_count))


if __name__ == "__main__":
    main()
