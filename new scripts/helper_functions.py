import math
import numpy as np
import bonding_distances as bond_dis


def distance_checker(xyz1, xyz2):
    """ Returns distance between 2 threedimensional points"""
    return math.sqrt((xyz1[0] - xyz2[0])**2 + (xyz1[1] - xyz2[1])**2 +
                     (xyz1[2] - xyz2[2])**2)


def normaliser(vec):
    """Normalises a vector"""
    norm = np.linalg.norm(vec)
    for i in range(len(vec)):
        vec[i] = vec[i] / norm
    return vec


def angle_checker(vec1, vec2):
    """ Calculates angle in radians between two vectors """
    vec1 = normaliser(vec1)
    vec2 = normaliser(vec2)
    angle = np.arccos(np.clip(np.dot(vec1, vec2), -1, 1))
    return angle


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. Taken from
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def bond_checker(atom, dict):
    "Check for all atoms in bonding range"
    bound = []
    for item, values in dict.items():
        bond_range = getattr(bond_dis, atom[0])().distances[values["element"]] + 0.1
        if (math.sqrt((atom[1] - values['x'])**2 + (atom[2] - values['y'])**2 +
                      (atom[3] - values['z'])**2) <= bond_range):
            bound.append(item)
    return bound


def y2true(text):
    while True:
        if text == 'y':
            return True
        elif text == 'n':
            return False
        else:
            text = input("Wrong input, try again: ")
