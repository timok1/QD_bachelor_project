# Timo Koster
import numpy as np
import math
import random
import copy
import sys

import helper_functions as hf
import unit_cells


def crystal_builder(a, atom_a, atom_b, diameter):
    """
    This function takes in the lattice constant, elements in the crystal and
    the diameter of the nanocrystal to be built. The unit cell is copied in
    all directions and cutoffs are made at the (111), (110) and (100) planes.
    ID numbers are assigned and neighbours are calculated.
    Then coordinates, elements, neighbours and atom id are stored in a dict.
    """
    n_unit = math.ceil(diameter)
    n_range = np.arange(-n_unit/2.0 + 0.5, n_unit/2.0)
    # Get unit cell
    uc = unit_cells.zns(a, atom_a, atom_b)
    all_atoms = []
    atom_dict = {}
    cut_off_distance = diameter * a / 2.0
    boundary_111 = 3.0 / math.sqrt(3) * cut_off_distance
    boundary_110 = 2.0 / math.sqrt(2) * cut_off_distance
    boundary_100 = cut_off_distance
    id = 0
    for x in n_range:
        for y in n_range:
            for z in n_range:
                # Place atoms one by one in correct spot
                for atom in uc.atom_xyz:
                    atom_element = atom[0]
                    atom_x = round(x*a + atom[1], 4)
                    atom_y = round(y*a + atom[2], 4)
                    atom_z = round(z*a + atom[3], 4)
                    atom_details = [atom_element, atom_x, atom_y, atom_z]

                    # Make sure not to put multiple atoms in same spot, and make cut-offs
                    if atom_details not in all_atoms and (abs(atom_x) + abs(atom_y) + abs(atom_z) <= boundary_111 and
                            (abs(atom_x) + abs(atom_y) <= boundary_110 and abs(atom_x) + abs(atom_z) <= boundary_110 and
                            abs(atom_y) + abs(atom_z) <= boundary_110) and (abs(atom_x) <= boundary_100
                            and abs(atom_y) <= boundary_100 and abs(atom_z) <= boundary_100)):
                        all_atoms.append(atom_details)
                        bound = hf.bond_checker(atom_details, atom_dict, bond_len_dict)
                        atom_dict[id] = {
                                        "coor": [atom_x, atom_y, atom_z],
                                        "element": atom_element,
                                        "bound": bound,
                                        "type": "crystal"
                                        }
                        # Update bonds for already placed atoms
                        for item in bound:
                            atom_dict[item]["bound"].append(id)
                        id += 1
    # Remove all singly bound atoms, if wanted
    include_singles = hf.y2true(input("Include singly bound atoms? y/n: "))
    if not include_singles:
        while True:
            for test_id, values in atom_dict.items():
                stopped = False
                # Delete singly bonded atoms, update neighbour
                if len(values['bound']) == 1:
                    neighbour = values['bound'][0]
                    atom_dict[neighbour]['bound'].remove(test_id)
                    del dict[test_id]
                    stopped = True
                    break
            # Break if no singly bonded atoms remain
            if not stopped:
                break
    return atom_dict


def crystal_reader(filename):
    """
    Takes an .xyz file and processes it in such a way that the output is
    identical to crystal_builder()
    """
    file = open("../input_crystals/" + filename, 'r')
    element_list = []
    # Register which elements belong to core, in attempt to allow for .xyz
    # files to contain ligands for a future feature
    while True:
        el = input("Type element in crystal, or type 'done': ")
        if el == 'done':
            break
        element_list.append(el)

    id = 0
    atom_dict = {}
    line_number = 0
    for line in file:
        if line_number == 0:
            n_atoms = line.strip()
            n_atoms = int(float(n_atoms))
        # Go to first element
        if line_number >= 2 and line_number < n_atoms + 2:
            values_list = line.split()
            for i in range(1, 4):
                values_list[i] = float(values_list[i])
            bound = hf.bond_checker(values_list, atom_dict, bond_len_dict)
            if values_list[0] in element_list:
                type = "crystal"
            else:
                type = "ligand"
            atom_dict[id] = {
                        "coor": values_list[1:],
                        "element": values_list[0],
                        "bound": bound,
                        "type": type
                        }
            # update bound for already placed atoms
            for item in bound:
                atom_dict[item]["bound"].append(id)
            id += 1
        line_number += 1
    return atom_dict


def series_lig(atom_dict, sites, foldername):
    """
    Obtains the lig_dict containing the values of the ligands to be added to
    the quantum dot. Here this is done by combining all 'bases' and extensions
    with eachother. e.g. when you put in C10H21 and C9H19 as bases, and COOH
    and NH2 as extensions you'll get C10H21 + COOH, C10H21 + NH2, C9H19 + COOH
    C9H19 + NH2
    """
    base_list = []
    ext_list = []
    lig_dict = {}
    n_sites = len(sites)
    sites_copy = copy.deepcopy(sites)
    # Obtain bases of ligands
    while True:
        base = input("Add ligand to base list, or type 'files' to show available ligands, or type 'done': ")
        if base == 'files':
            hf.print_lig()
        if base == 'done':
            break
        base_list.append(base + ".xyz")
    # Obtain extensions
    while True:
        ext = input("Add extension to extension list, or type 'none', or 'done': ")
        if ext == 'done':
            break
        if ext == 'none':
            ext_list.append([False])
        # Extend extensions
        else:
            ext_list.append([ext + ".xyz"])
            extra_ext = input("Further extension to extension, or type n: ")
            if extra_ext != 'n':
                ext_list[-1].append(extra_ext + ".xyz")

    coverage = float(input("Coverage (fraction): "))
    cap = input("Cap remaining sites with (single atom): ") + ".xyz"
    sec_buffer = float(input("Buffer size (Angstrom)? Distance between ligands will be at least bonding length + buffer: "))
    buffer = float(input("Initial buffer (Buffer for first QD. Helps create space for later QD's): "))
    atom_dict_copy = copy.deepcopy(atom_dict)
    i = 0

    # Create all quantum dots
    for base in base_list:
        print(base)
        for extension in ext_list:
            filename = base[:-4]
            for item in extension:
                if item is not False:
                    filename += "+" + item[:-4]
            # Create first quantum dot
            if i == 0:
                lig_dict[0] = {
                        "ligand_type": base,
                        "extension": extension,
                        "coverage": coverage,
                        "n_ligands": round(coverage * n_sites)
                }

                lig_dict[1] = {
                        "ligand_type": cap,
                        "extension": [False],
                        "coverage": 1,
                        "n_ligands": n_sites
                    }
                while True:
                    stopped = False
                    for ligand, values in lig_dict.items():
                        atom_dict = place_ligands(atom_dict, values, sites, buffer, False, cap)
                        # Detect failure to place all ligands
                        if atom_dict == 1:
                            atom_dict = copy.deepcopy(atom_dict_copy)
                            sites = copy.deepcopy(sites_copy)
                            stopped = True
                            break
                    if not stopped:
                        i += 1
                        buffer = sec_buffer
                        hf.dict2file(atom_dict, filename, foldername)
                        break
            # Replace ligands
            else:
                sites = sites_copy.copy()
                input_rep = [base, extension]
                replaced = replace_ligands(atom_dict, lig_dict, cap, sites, buffer, input_rep, filename, foldername)
                atom_dict = replaced[0]
                lig_dict = replaced[1]
    # Allow for more replacements
    while True:
        replacement_ligands = hf.y2true(input("Replace ligands to create new quantum dot? y/n: "))
        if replacement_ligands:
            sites = sites_copy.copy()
            filename = input("Write to file: ")
            buffer = float(input("Buffer: "))
            replaced = replace_ligands(atom_dict, lig_dict, cap, sites, buffer, False, filename, foldername)
            atom_dict = replaced[0]
            lig_dict = replaced[1]
        else:
            break


def single_qd(atom_dict, sites, foldername):
    """
    Create a single quantum dot and allow for replacement of ligands for a new
    one.
    """
    filename = input("Save QD as: ")
    i = 0
    lig_list = []
    extension_list = []
    coverage_list = []
    n_ligands_list = []
    lig_dict = {}
    n_sites = len(sites)
    # Get all ligands
    while True:
        ligand_file = input("Filename of ligand to be added, or type files to see available ligands: ") + ".xyz"
        # Print available ligands
        if ligand_file[:-4] == "files":
            hf.print_lig()
            ligand_file = input("Filename of ligand to be added: ") + ".xyz"
        lig_list.append(ligand_file)
        extend = hf.y2true(input("Extend ligand? Note: this replaces the last atom in the ligand file. y/n: "))
        if extend:
                extension_list.append([input("Extend with: ") + ".xyz"])
                while True:
                    extra_ext = input("Further extension to extension, or type n: ")
                    if extra_ext != 'n':
                        extension_list[-1].append(extra_ext + ".xyz")
                    else:
                        break
        else:
            extension_list.append([False])
        coverage = float(input("Coverage (fraction): "))
        coverage_list.append(coverage)
        n_ligands_list.append(round(coverage * n_sites))
        lig_dict[i] = {
                "ligand_type": ligand_file,
                "extension": extension_list[-1],
                "coverage": coverage,
                "n_ligands": round(coverage * n_sites)
        }
        i += 1
        more = hf.y2true(input("Add another type? y/n: "))
        if not more:
            cap = input("Cap remaining sites with (single atom): ") + ".xyz"
            lig_dict[i] = {
                    "ligand_type": cap,
                    "extension": [False],
                    "coverage": 1,
                    "n_ligands": n_sites
            }
            break

    buffer = float(input("Buffer (Angstrom) (distance between ligands will be at least buffer + bonding length): "))
    sites_copy = copy.deepcopy(sites)
    atom_dict_copy = copy.deepcopy(atom_dict)
    # Place ligands
    while True:
        stopped = False
        for ligand, values in lig_dict.items():
            atom_dict = place_ligands(atom_dict, values, sites, buffer, False, cap)
            if atom_dict == 1:
                atom_dict = copy.deepcopy(atom_dict_copy)
                sites = copy.deepcopy(sites_copy)
                stopped = True
                break
        if not stopped:
            break
    hf.dict2file(atom_dict, filename, foldername)

    # Replace ligands
    while True:
        replacement_ligands = hf.y2true(input("Replace ligands to create new quantum dot? y/n: "))
        if replacement_ligands:
            sites = sites_copy.copy()
            filename = input("Write to file: ")
            buffer = float(input("Buffer: "))
            replaced = replace_ligands(atom_dict, lig_dict, cap, sites, buffer, False, filename, foldername)
            atom_dict = replaced[0]
            lig_dict = replaced[1]
        else:
            break


def replace_ligands(atom_dict, lig_dict, cap, sites, buffer, input_rep, filename, foldername):
    """
    Takes in an atom dict with ligands and replaces the ligands by the
    requested replacements, saves it to a new quantum dot and returns the
    new atom dict.
    """
    replacement_list = []
    rep_ext_list = []
    loc_dict = {cap[:-4]: []}
    rep_dict = copy.deepcopy(lig_dict)
    lig_dict_inv = {}
    # Store the replacements and replacee's in the right way
    for ligand, lig_val in lig_dict.items():
        if lig_val["ligand_type"] != cap:
            # Use existing input
            if input_rep:
                replacement_list = [input_rep[0]]
                rep_ext_list = [input_rep[1]]
                rep_dict[ligand]["extension"] = rep_ext_list[-1]
            # Get new input
            else:
                replacement_list.append(input("Replace " + lig_val["ligand_type"] + " with: ") + ".xyz")
                extend = hf.y2true(input("Extend ligand? y/n: "))
                if extend:
                    rep_ext_list.append([input("Extend with: ") + ".xyz"])
                    rep_dict[ligand]["extension"] = rep_ext_list[-1]
                else:
                    rep_ext_list.append([False])
                    rep_dict[ligand]["extension"] = [False]
            # Store info in dict
            loc_dict[ligand] = {
                                "loc_id": [],
                                "replacement": replacement_list[-1]
                            }
            rep_dict[ligand]["ligand_type"] = replacement_list[-1]
            rep_dict[ligand]["loc_info"] = {}
            # Reverses key and value
            lig_dict_inv[lig_val["ligand_type"]] = ligand

    # Get info on all atoms that need to be deleted, and get the relevant info
    # of the replacee's
    atom_del_list = []
    for atom, at_values in atom_dict.items():
        if at_values["type"] == 'ligand':
            atom_del_list.append(atom)
            if at_values["ligand_type"] != cap:
                rep_dict[lig_dict_inv[at_values["ligand_type"]]]["loc_info"][at_values["loc_id"]] = {
                    "loc_vec": at_values["loc_vec"],
                    "rotation": at_values["rotation"]
                }

    rep_dict_copy = copy.deepcopy(rep_dict)
    atom_dict_copy = copy.deepcopy(atom_dict)
    sites_copy = copy.deepcopy(sites)
    done = False
    # Actually delete and replace
    while not done:
        for item in atom_del_list:
            del atom_dict[item]
        for lig, values in rep_dict.items():
            atom_dict = place_ligands(atom_dict, values, sites, buffer, True, cap)
            if atom_dict == 1:
                atom_dict = copy.deepcopy(atom_dict_copy)
                sites = copy.deepcopy(sites_copy)
                rep_dict = copy.deepcopy(rep_dict_copy)
                break
            else:
                done = True
    hf.dict2file(atom_dict, filename, foldername)

    return atom_dict, rep_dict


def tetra_sites(dict):
    """
    Takes in the atom_dict, and determines all possible connection sites for a
    tetrahedral structure. Returns a dict with the atom they connect to as a
    key, that atom's coordinates, and a unit vector from that atom in the
    direction of the site.
    """
    # Create a standard tetrahedron to determine sites
    x = 1 / math.sqrt(3)
    sites = {}
    # Determine possible sites atom by atom
    for id, values in dict.items():
        if len(values['bound']) < 4 and values["type"] == "crystal":
            sites[id] = {}
            tetrahedron = [[x, x, x], [x, -x, -x], [-x, x, -x], [-x, -x, x]]
            primary_xyz = values["coor"]
            connection_list = []
            # Create unit vectors for all connected atoms
            for secondary in values['bound']:
                sec_xyz_rel = [c1 - c2 for c1, c2 in zip(dict[secondary]["coor"], primary_xyz)]
                sec_xyz_rel = hf.normaliser(sec_xyz_rel)
                connection_list.append(sec_xyz_rel)

            # line up first vector
            axis = np.cross(tetrahedron[0], connection_list[0])
            # Check if already lined up
            if round(np.linalg.norm(axis), 2) == 0:
                axis = [0, x, -x]
                angle = 0
                # check if opposite
                if round(np.dot(tetrahedron[0], connection_list[0]), 2) == -1:
                    angle = math.pi
            else:
                angle = hf.angle_checker(tetrahedron[0], connection_list[0])
            rot_mat = hf.rotation_matrix(axis, angle)
            for i in range(4):
                tetrahedron[i] = list(np.dot(rot_mat, tetrahedron[i]))

            # line up second vector
            axis = tetrahedron[0]
            # Random rotation if single bond, and skip to next atom
            if len(connection_list) == 1:
                angle = random.random() * 2 * math.pi
                rot_mat = hf.rotation_matrix(axis, angle)
                for i in range(4):
                    tetrahedron[i] = list(np.around(np.array(np.dot(rot_mat, tetrahedron[i])), 4))
                sites[id] = {
                            'primary_xyz': primary_xyz,
                            'sites_xyz': [tetrahedron[1], tetrahedron[2], tetrahedron[3]]}
                continue
            # No rotation if already aligned
            elif round(np.dot(tetrahedron[1], connection_list[1]), 2) == 1:
                angle = 0
            # Get angle between planes spanned by first vector, second
            # connection vector and second vector in tetrahedron
            else:
                normal_1 = np.cross(tetrahedron[0], connection_list[1])
                normal_2 = np.cross(tetrahedron[0], tetrahedron[1])
                angle = hf.angle_checker(normal_1, normal_2)

            # Make sure rotation is in right direction, as angle checker gives
            # absolute angle
            test_rot_mat = hf.rotation_matrix(axis, angle)
            test_rahedron = tetrahedron.copy()
            for test in range(4):
                test_rahedron[test] = list(np.dot(test_rot_mat, tetrahedron[test]))
            if round(np.dot(connection_list[1], test_rahedron[1]), 1) != 1:
                angle = -angle

            # Rotate tetrahedron to correct orientation, add points without atom to sites
            rot_mat = hf.rotation_matrix(axis, angle)
            for i in range(4):
                tetrahedron[i] = list(np.around(np.array(np.dot(rot_mat, tetrahedron[i])), 4))
            if len(values['bound']) == 2:
                sites[id] = {
                            'primary_xyz': primary_xyz,
                            'sites_xyz': [tetrahedron[2], tetrahedron[3]]}
            elif len(values['bound']) == 3:
                if round(np.dot(tetrahedron[2], connection_list[2]), 1) == 1:
                    sites[id] = {
                                'primary_xyz': primary_xyz,
                                'sites_xyz': [tetrahedron[3]]}
                else:
                    sites[id] = {
                                'primary_xyz': primary_xyz,
                                'sites_xyz': [tetrahedron[2]]}
    return sites


def prep_ligand_file(atom_dict, ligand_type, extension, cap):
    """Read ligand file, possibly extend with another ligand and add to a dict"""

    path = "../Ligands/" + ligand_type
    file = open(path, 'r')
    id = 0
    lig = {}
    # Get element of atom connected to crystal
    lig = hf.file2dict(file, lig, id)
    file.close()

    # Extend ligand if needed
    if extension[0] is not False and ligand_type != cap:
        for item in extension:
            lig = extend_ligand(atom_dict, lig, item)
    return lig


def extend_ligand(atom_dict, lig, extension):
    """
    Takes a ligand dict, removes the last atom and adds all atoms in the
    requested extension in the righ orientation. Returns a new ligand dict.
    """
    path_ext = "../Ligands/" + extension
    file_ext = open(path_ext, 'r')
    id = max(lig) + 1
    rep = lig[id - 1]
    rep_coor = rep["coor"]
    ext = {}
    ext = hf.file2dict(file_ext, ext, id)
    file_ext.close()
    # Check for atom to connect extension to, will be closest atom
    closest_atom = hf.closest_atom(lig, rep_coor)
    base_atom_el = ext[hf.base_atom(ext)]["element"]
    init_length_ext = hf.check_bond_len(bond_len_dict, base_atom_el, lig[closest_atom]["element"])
    # Remove last atom
    del lig[max(lig)]

    # Get correct rotations for extension
    unit_v = hf.normaliser([c1 - c2 for c1, c2 in zip(rep_coor, lig[closest_atom]['coor'])])
    angle = hf.angle_checker(unit_v, [0, 0, 1])
    axis = np.cross(unit_v, [0, 0, 1])
    if round(np.linalg.norm(axis), 2) == 0:
        axis = [1, 0, 0]
    rot_mat = hf.rotation_matrix(axis, angle)
    # test for correct angle
    test_v = np.dot(rot_mat, unit_v)
    if round(np.dot(test_v, unit_v), 2) != 1:
        rot_mat = hf.rotation_matrix(axis, -angle)

    # Add all atoms in extension to ligand
    for atom, values in ext.items():
        values["coor"][2] += init_length_ext
        ext_coor = np.dot(rot_mat, values["coor"])
        ext[id] = {
                "coor": [c1 + c2 for c1, c2 in zip(ext_coor, lig[closest_atom]["coor"])],
                "element": values["element"]
                }
        id += 1
    lig = {**lig, **ext}
    return lig


def place_ligands(atom_dict, lig_info, sites, buffer, fixed_loc, cap):
    """
    Adds the requested amount of ligands at random sites, returns updated
    atom_dict.
    """
    # Sorry for the length of this function, should probably be broken up in
    # multiple functions.
    ligand_type = lig_info["ligand_type"]
    n_ligands = lig_info["n_ligands"]
    extension = lig_info["extension"]
    fail_bool = False
    id = max(atom_dict) + 1
    original_ligand = ligand_type
    j = 0
    tried = []
    # Print which ligands are being placed
    if extension and extension[0] is not False:
        ex_string = ''
        for x in extension:
            if x is not False:
                ex_string += x[:-4] + ' + '
        ex_string = ex_string[:-3]
        print("\nPlacing " + ligand_type[:-4] + " + " + ex_string)
    else:
        print("\nPlacing " + ligand_type[:-4])

    # Add ligands until done
    while j < n_ligands:
        # preventive in case of rounding errors
        if len(sites) == 0:
            break
        extra_rotation = 0
        ligand_type = original_ligand

        # Randomly choose site, get relevant info, don't pick site that has already been tried
        sites_list = list(sites).copy()
        remaining_sites = [x for x in sites_list if x not in tried]

        if len(remaining_sites) == 0:
            print("\nUnable to place more ligands of type " + str(original_ligand) + ". Placed "
                  + str(j) + " out of " + str(n_ligands) + " requested ligands.\n")
            retry = hf.y2true(input("Try again (y), or continue (n)?: "))
            if retry:
                print("Trying again...")
                return 1
            break
        # Use predetermined sites if required, "loc" variables refer to
        # connection site and its primary atom
        if fixed_loc and ligand_type != cap:
            if len(lig_info["loc_info"]) == 0:
                break
            loc_id = random.choice(list(lig_info["loc_info"]))
            loc_vec = lig_info["loc_info"][loc_id]["loc_vec"]
            loc_sites = copy.deepcopy(sites[loc_id]['sites_xyz'])
            random_rotation = lig_info["loc_info"][loc_id]["rotation"]
            del lig_info["loc_info"][loc_id]
        # Use random sites
        else:
            loc_id = random.choice(remaining_sites)
            loc_sites = sites[loc_id]['sites_xyz'].copy()
            loc_vec = random.choice(loc_sites)
            random_rotation = random.random() * 2 * math.pi
        loc_primary_xyz = sites[loc_id]['primary_xyz']
        tried_loc_vec = []
        lig = prep_ligand_file(atom_dict, ligand_type, extension, cap)
        stop = False

        # Loop over every site connected to chosen atom
        while True:
            # Make sure you don't try same loc_vec twice
            for vec in tried_loc_vec:
                if hf.distance_checker(vec, loc_vec) < 0.001:
                    stop = True
            if stop:
                break
            # Get correct rotation for ligand relative to site
            axis = np.cross(loc_vec, [0, 0, 1])
            # Prevent dividing by 0 when vectors are already lined up
            if np.linalg.norm(axis) == 0:
                axis = [1, 0, 0]
                angle = 0
            else:
                angle = -hf.angle_checker(loc_vec, [0, 0, 1])
            rot_mat = hf.rotation_matrix(axis, angle)

            # Place ligand if possible, otherwise rotate
            while True:
                xyz_list = []
                temp_atom_dict = {}
                broken = False
                min_dist_lig = math.inf
                initial_length = hf.check_bond_len(bond_len_dict, lig[hf.base_atom(lig)]["element"], atom_dict[loc_id]["element"])
                for atom, lig_val in lig.items():
                    lig_coor = copy.deepcopy(lig_val["coor"])
                    lig_coor[2] += initial_length
                    # Rotate every atom in ligand according to random_rotation
                    rotated_atom = np.dot(hf.rotation_matrix([0, 0, 1], random_rotation + extra_rotation), lig_coor)
                    # Rotate in right direction (with respect to crystal)
                    new_v = np.dot(rot_mat, rotated_atom)
                    atom_xyz = [c1 + c2 for c1, c2 in zip(new_v, loc_primary_xyz)]
                    atom_element = lig[atom]['element']
                    temp_atom_dict[id] = {
                                    "coor": atom_xyz,
                                    "element": atom_element,
                                    "type": "ligand",
                                    "ligand_type": ligand_type,
                                    "loc_id": loc_id,
                                    "loc_vec": loc_vec,
                                    "rotation": random_rotation + extra_rotation
                                    }
                    xyz_list.append(atom_xyz)
                    id += 1
                    # Check for overlap
                    if ligand_type != cap:
                        for test_atom, values in atom_dict.items():
                            space = hf.check_bond_len(bond_len_dict, atom_element, values["element"]) + buffer
                            if test_atom != loc_id:
                                dist = hf.distance_checker(values["coor"], xyz_list[-1])
                                if dist < space:
                                    broken = True
                                    break
                    # Check for overlap when capping, use different rotation
                    if ligand_type == cap:
                        new_pos = new_v.copy()
                        changed = False
                        while True:
                            min_dist_lig = math.inf
                            for test_lig, values in atom_dict.items():
                                if values["type"] == "ligand":
                                    dist = hf.distance_checker(values["coor"], [c1 + c2 for c1, c2 in zip(new_pos, loc_primary_xyz)])
                                    bond_len_loc = hf.check_bond_len(bond_len_dict, atom_element, values["element"])
                                    # Sorry for magic number
                                    if dist < bond_len_loc + 0.3:
                                        if dist < min_dist_lig:
                                            min_dist_lig = dist
                                            min_dist_lig_coor = values["coor"]
                                            changed = True
                            if min_dist_lig == math.inf:
                                if changed:
                                    xyz_list = [c1 + c2 for c1, c2 in zip(new_pos, loc_primary_xyz)]
                                    temp_atom_dict[id-1] = {
                                                    "coor": xyz_list,
                                                    "element": atom_element,
                                                    "type": "ligand",
                                                    "ligand_type": ligand_type,
                                                    "loc_id": loc_id,
                                                    "loc_vec": loc_vec
                                                    }
                                break
                            # Rotate away from closest atom
                            vec_closest_element = [c1 - c2 for c1, c2 in zip(min_dist_lig_coor, loc_primary_xyz)]
                            axis = np.cross(vec_closest_element, new_pos)
                            angle = 0.1
                            rot_mat = hf.rotation_matrix(axis, angle)
                            new_pos = np.dot(rot_mat, new_pos.copy())

                    if broken:
                        break
                # Check for ligand collisions except when last ligand is added
                if ligand_type != cap:
                    # Rotate if there is a collision
                    if broken and extra_rotation < math.pi * 2:
                        # Rotate less with fixed loc
                        if fixed_loc:
                            extra_rotation += 0.1
                        else:
                            extra_rotation += 0.2
                        continue
                    # break if rotated fully
                    elif extra_rotation > math.pi * 2:
                        extra_rotation = 0
                        tried.append(loc_id)
                        tried_loc_vec.append(loc_vec)
                        fail_bool = True
                        break
                    # Merge dicts if no collision
                    elif not broken:
                        atom_dict = {**atom_dict, **temp_atom_dict}
                        if j % 10 == 0:
                            print(str(j) + " ligands added")
                        if loc_id in sites:
                            del sites[loc_id]
                        tried_loc_vec = [loc_vec]
                        ligand_type = cap
                        lig = prep_ligand_file(atom_dict, cap, [False], cap)
                        j += 1
                        break
                else:
                    atom_dict = {**atom_dict, **temp_atom_dict}
                    if loc_id in sites:
                        del sites[loc_id]
                    tried_loc_vec.append(loc_vec)
                    break
            # Add final ligand type to remaining sites at chosen atom
            try:
                loc_vec = loc_sites[(loc_sites.index(loc_vec) + 1) % len(loc_sites)]
            except ValueError:
                min_dist = math.inf
                loc_index = 0
                for item in loc_sites:
                    distance = hf.distance_checker(loc_vec, item)
                    if distance < min_dist:
                        min_dist = distance
                        loc_index_min = loc_index
                    loc_index += 1
                loc_vec = loc_sites[(loc_index_min + 1) % len(loc_sites)]

    # Detect failure to place all requested ligands if replacing
    if fixed_loc and fail_bool:
        retry = hf.y2true(input("Not all ligands were able to be put in the same spot, try again? y/n: "))
        if retry:
            print("Trying again...")
            return 1
        else:
            cont = hf.y2true(input("Continue with fewer ligands? y/n: "))
        if not cont:
            sys.exit()

    return atom_dict


def bridges(atom_dict, sites):
    """
    Tried to build a function that could add bridge ligands, but gave up.
    Just left it in here to maybe help others start, or offer some
    schadenfreude
    """
    couples = []
    bridge_dict = {}
    bridge_id = 0
    tried = []
    for primary_site, values_prim in sites.items():
        for secondary_site, values_sec in sites.items():
            if secondary_site != primary_site:
                for xyz1 in values_prim["sites_xyz"]:
                    for xyz2 in values_sec["sites_xyz"]:
                        dist = hf.distance_checker([c1 + c2 for c1, c2 in zip(values_prim["primary_xyz"], xyz1)], [c1 + c2 for c1, c2 in zip(values_sec["primary_xyz"], xyz2)])
                        if dist < 2.5:
                            couples.append([primary_site, secondary_site])
                            coor1 = values_prim["primary_xyz"]
                            coor2 = values_sec["primary_xyz"]
                            new_loc = [(c1 + c2) / 2 for c1, c2 in zip(coor1, coor2)]
                            temp_xyz = [(c1 + c2) / 2 for c1, c2 in zip(xyz1, xyz2)]
                            new_site = [(c1 - c2) for c1, c2 in zip(temp_xyz, new_loc)]
                            temp_xyz = hf.normaliser(temp_xyz)
                            if new_site not in tried:
                                bridge_dict[bridge_id] = {
                                                        "primary_xyz": new_loc,
                                                        "sites_xyz":  temp_xyz,
                                                        "connected": couples[-1]
                                }
                                bridge_id += 1
                                tried.append(new_site)
    bridge_dict_length = len(bridge_dict)
    n_bridges = int(input(str(bridge_dict_length) + " possible bridge sites found. How many ligands should be placed at these sites?: "))
    bridge_ligand = input("type of ligand to be placed at bridge sites: ")
    for i in range(bridge_dict_length - n_bridges):
        cut = random.choice(list(bridge_dict))
        del bridge_dict[cut]
    atom_dict = place_bridge_ligands(atom_dict, sites, bridge_dict, bridge_ligand)
    return atom_dict


def place_bridge_ligands(atom_dict, sites, bridge_dict, bridge_ligand):
    """
    See bridges()
    """
    id = max(atom_dict) + 1
    for bridge, values in bridge_dict.items():
        random_rotation = math.pi/4
        # Determine height of ligand. If bonding distance is too short, set to 0.
        # This is only correct for crystal with one kind of atom
        lig = prep_ligand_file(atom_dict, "/DB_Ligands/" + bridge_ligand + ".xyz", False)
        try:
            initial_length = math.sqrt(bond_len_dict[lig.base_element]["Si"]  - (3.84/2)**2)
        except ValueError:
            initial_length = 0
        # Get correct rotation for ligand relative to site
        axis = np.cross(values["sites_xyz"], [0, 0, 1])
        # Prevent dividing by 0 when vectors are already lined up
        if np.linalg.norm(axis) == 0:
            axis = [1, 0, 0]
            if np.dot(values["sites_xyz"], [0, 0, 1]) > 0:
                angle = 0
            else:
                angle = math.pi
        else:
            angle = -hf.angle_checker(values["sites_xyz"], [0, 0, 1])
        rot_mat = hf.rotation_matrix(axis, angle)
        if values["connected"][0] in sites:
            del sites[values["connected"][0]]
        if values["connected"][1] in sites:
            del sites[values["connected"][1]]
        while True:
            xyz_list = []
            temp_atom_dict = {}
            test_vx = [1, 0, 0]
            test_vx_rot = np.dot(rot_mat, [1, 0, 0])
            test_vz_rot = np.dot(rot_mat, [0, 0, 1])
            uhh = atom_dict[values["connected"][0]]["coor"]
            new_ding = [c1 - c2 for c1, c2 in zip(uhh, values["primary_xyz"])]
            new_ding = hf.normaliser(new_ding)
            normal_1 = np.cross(test_vx_rot, test_vz_rot)
            normal_2 = np.cross(new_ding, test_vz_rot)
            angle2 = hf.angle_checker(normal_1, normal_2)
            random_rotation = math.pi / 2 - angle2
            vars2 = np.dot(hf.rotation_matrix([0, 0, 1], random_rotation), test_vx)
            vars3 = np.dot(rot_mat, vars2)
            angle3 = round(hf.angle_checker(vars3, normal_2), 3)
            if angle3 != 0 and angle3 != round(math.pi, 3):
                random_rotation = random_rotation + math.pi / 2
            random_rotation = random_rotation + random.choice((0, math.pi))

            # Rotate every atom in ligand according to random_rotation
            for atom in lig.atoms:
                new_v_temp = np.dot(rot_mat, [lig.atoms[atom]['x'], lig.atoms[atom]['y'], lig.atoms[atom]['z']])
                new_v = np.dot(hf.rotation_matrix(np.dot(rot_mat, [0, 0, 1]), random_rotation), new_v_temp)
                atom_xyz = [c1 + c2 for c1, c2 in zip(new_v, values["primary_xyz"])]
                atom_element = lig.atoms[atom]['element']
                temp_atom_dict[id] = {
                                "x": atom_xyz[0],
                                "y": atom_xyz[1],
                                "z": atom_xyz[2],
                                "element": atom_element,
                                "type": "ligand",
                                "ligand_type": bridge_ligand
                                }
                xyz_list.append(atom_xyz)
                id += 1
                for test_atom, values2 in atom_dict.items():
                    space = bond_len_dict[atom_element][values2["element"]]
                    # Stay further away from the crystal atoms
                    if values2['type'] == "crystal":
                        space = space + 0.25
                    if test_atom not in values["connected"]:
                        dist = hf.distance_checker(values2["coor"], xyz_list[-1])
                        if dist < space + 0.25:
                            break
            atom_dict = {**atom_dict, **temp_atom_dict}
            break
    return atom_dict


if __name__ == "__main__":
    # Having a global dict with bonding lengths improves speed a lot
    global bond_len_dict
    bond_len_dict = hf.csv2dict("bonding_distances.csv")
    build = hf.y2true(input("Create new crystal (y) or use existing file (n)?: "))
    if build:
        a = float(input("Specify lattice constant (in Ångström): "))
        atom_a = input("Element for first element type: ")
        atom_b = input("Element for second element type: ")
        diameter = float(input("Diameter of quantum dot (in unit cells): "))
        atom_dict = crystal_builder(a, atom_a, atom_b, diameter)
    else:
        crystal_file = input("Crystal file to use (don't write the file extension): ") + ".xyz"
        atom_dict = crystal_reader(crystal_file)

    foldername = input("Save in folder (or main): ")
    if foldername == "main":
        foldername = False
    sites = tetra_sites(atom_dict)
    # Get a series of ligands to add to crystal
    series_b = hf.y2true(input("Add a series of ligands (y), or add them one by one(n)?: "))
    if series_b:
        series_lig(atom_dict, sites, foldername)
    # Add ligands one by one
    else:
        single_qd(atom_dict, sites, foldername)
