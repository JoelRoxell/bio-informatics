import os
import string
import numpy as np
import sys
import math
import sys

# Joel Roxell

BOUND_LIMIT = 1


def distance_between(atom_1, atom_2):
    x1, y1, z1 = atom_1
    x2, y2, z2 = atom_2

    return math.sqrt(((x2-x1) ** 2) + ((y2-y1) ** 2) + ((z2-z1) ** 2))


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


# Will contain the entire statespace.
relations = []
LOWER_BOUND = 3.6
UPPER_BOUND = 3.9


def collect_states(atom, pre, atom_list, path=''):
    ID, coordinates = atom

    for next_atom in atom_list:
        next_ID, next_coordinates = next_atom

        if ID == next_ID:
            continue

        if pre is not None:
            if pre[0] == next_ID:
                continue

        distance = distance_between(coordinates, next_coordinates)

        if distance <= UPPER_BOUND and distance > LOWER_BOUND:
            relations.append((int(ID), int(next_ID)))


def get_children_of(target):
    children = []

    for ID, child in relations:
        if ID == target:
            children.append(int(child))

    return children


class Node:
    def __init__(self, ID, children, count, parent):
        self.id = ID
        self.children = children
        self.count = count
        self.parent = parent


def get_coordinate_of(atom_id):
    for atom in atom_list:
        if int(atom[0]) == atom_id:
            return atom[1]


def calculateCentroid(children):
    X = []
    Y = []
    Z = []

    for child in children:
        for atom in atom_list:
            if int(atom[0]) == child:
                X.append(atom[1][0])
                Y.append(atom[1][1])
                Z.append(atom[1][2])

    sum_x = sum(X)
    sum_y = sum(Y)
    sum_z = sum(Z)
    n = len(X)

    return (sum_x/n, sum_y/n, sum_z/n)


def find_path(node, explored):
    explored.update({
        node.id: node
    })

    if not node.children:
        return explored

    # Problem needs to be relaxed, there are too many connection and cannot be exhausted, that is, event with a lower amount of relations(~700).
    centroid = calculateCentroid(node.children)
    distances = []
    children_closest_to_centroid = []

    for child in node.children:
        distance_to_centroid = distance_between(
            centroid,
            get_coordinate_of(child)
        )
        distances.append(distance_to_centroid)

    closest_atom = min(distances)
    closest_index = distances.index(closest_atom)

    t = node.children[closest_index]

    if node.parent and node.parent.id == t:
        # deadend or parent
        del distances[closest_index]

        if len(distances) == 0:
            return

        closest_atom = round(min(distances), 5)

    for i in range(0, len(distances)):
        bound = abs(closest_atom - round(distances[i], 5))

        if bound < BOUND_LIMIT:
            # More than one node could be equally close to the centroid.
            children_closest_to_centroid.append(node.children[i])

    for child in children_closest_to_centroid:
        nextNode = Node(
            child,
            get_children_of(child),
            node.count + 1,
            node
        )
        previous = explored.get(nextNode.id)

        if node.parent and nextNode.id == node.parent.id:
            # deadend or parent
            continue

        if not previous:
            # prevent loops
            find_path(nextNode, explored)

    return explored


def trace_back(node, list):
    if not node:
        return

    list.append(node.id)
    trace_back(node.parent, list)


atom_list = read_file('data_q2.txt')


def main():
    for atom in atom_list:
        ID = atom[0]
        collect_states(atom, None, atom_list, path=ID)

    # For each node bruteforce the longest path.
    longest_chain = []

    for atom in atom_list:
        ID = int(atom[0])
        start = Node(ID, get_children_of(ID), 0, None)

        current = start
        res = find_path(start, {})

        for _, node in res.items():
            if node.count >= current.count:
                current = node

        stack = []
        trace_back(current, stack)

        if len(stack) > len(longest_chain):
            longest_chain = stack

    for atom in longest_chain:
        print(atom)

    print('total number of alpha-carbons: {}'.format(len(longest_chain)))


if __name__ == '__main__':
    main()
