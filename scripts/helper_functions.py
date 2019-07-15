import math
import numpy as np
import sys
import csv
import os


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
    """
    Calculates angle in radians between two vectors. Note: this is an absolute
    angle.
    """
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


def check_bond_len(dict, el_a, el_b):
    """ Make sure all elements are in bond_len_dict, and return the value"""
    if el_a in dict:
        if el_b in dict[el_a]:
            return dict[el_a][el_b]
    print()
    print(el_a + " and " + el_b + " bond length currently unsupported. Add value to the csv file.")
    sys.exit()


def csv2dict(filename):
    """
    Transforms the bond_lengths.csv to a dict
    """
    dis_dict = {}
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            el_a = row["Element Name"]
            dis_dict[el_a] = {}
            for entry in row:
                if entry != "Element Name":
                    dis_dict[el_a][entry] = float(row[entry])
    csvfile.close()
    return dis_dict


def bond_checker(atom, dict, bond_dict):
    """Check for all atoms in bonding range"""
    bound = []
    for item, values in dict.items():
        bond_range = check_bond_len(bond_dict, atom[0], values["element"]) + 0.2
        if distance_checker(atom[1:], values["coor"]) <= bond_range:
            bound.append(item)
    return bound


def closest_atom(dict, coor):
    "Given a dict and coordinates returns the closest atom"
    min_dis = math.inf
    for atom, values in dict.items():
        dis = distance_checker(coor, values["coor"])
        if dis < min_dis and dis > 0.01:
            min_dis = dis
            min_id = atom
    return min_id


def print_lig():
    """ Prints available ligands """
    lig_list = os.listdir("../Ligands")
    print()
    for ligs in lig_list:
        # Skip folders
        if ligs[-4:] == ".xyz":
            print(ligs[:-4])
    print()


def file2dict(file, dict, start_id):
    """
    Builds simple dict out of .xyz file, containing just id, elements and
    coordinates
    """
    id = start_id
    line_number = 0
    file.seek(0)
    for line in file:
        if line_number == 0:
            n_atoms = int(float(line.strip()))
        if line_number >= 2 and line_number < n_atoms + 2:
            values_list = line.split()
            for i in range(1, 4):
                values_list[i] = float(values_list[i])
            dict[id] = {
                    "coor": values_list[1:],
                    "element": values_list[0]
                    }
            id += 1
        line_number += 1
    return dict


def dict2file(dict, filename, foldername):
    """
    Takes an atom dict and writes it to an .xyz file in foldername in
    /Created_QD with filename as name for the file
    """
    if foldername:
        if not os.path.exists("../Created_QD/" + foldername):
            os.makedirs("../Created_QD/" + foldername)
        file = open("../Created_QD/" + foldername + "/" + filename + ".xyz", "w")
    else:
        file = open("../Created_QD/" + filename + ".xyz", "w")
    file.write("        \n\n")
    for atom, values in dict.items():
        file.write(values['element'] + "\t" + str(values['coor'][0]) + "\t\t" +
                   str(values['coor'][1]) + "\t\t" + str(values['coor'][2]) + "\n")
    file.seek(0)
    file.write(str(len(dict)))
    file.close()
    print("\nQuantum Dot created :)")


def base_atom(dict):
    """ Finds atoms at the origin in a dict, returns its id"""
    for atom, values in dict.items():
        xyz = values["coor"]
        if xyz[0] == xyz[1] == xyz[2] == 0:
            return atom


def y2true(text):
    """Converts strings y and n to boolean"""
    while True:
        if text == 'y':
            return True
        elif text == 'n':
            return False
        else:
            text = input("Wrong input, try again: ")
