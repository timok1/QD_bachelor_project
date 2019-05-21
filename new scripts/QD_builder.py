import unit_cells
import ligands
import numpy as np
import math
import random
import bonding_distances as bond_dis
import os
import helper_functions as hf


# Build crystal with input
def crystal_builder(structure, a, atom_a, atom_b, diameter, filename):
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
    bond_range = getattr(bond_dis, atom_a)().distances[atom_b] + 0.1
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
                        bound = hf.bond_checker(atom_details, atom_dict, bond_range)
                        atom_dict[id] = {
                                        "x": atom_x,
                                        "y": atom_y,
                                        "z": atom_z,
                                        "element": atom_element,
                                        "bound": bound,
                                        "type": "crystal"
                                        }
                        # Update bonds for already placed atoms
                        for item in bound:
                            atom_dict[item]["bound"].append(id)
                        id += 1
    # Remove all singly bound atoms, if wanted
    include_singles = input("Include singly bound atoms? y/n: ")
    if include_singles == 'n':
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


# Read xyz-file and place atoms in a dict.
def crystal_reader(filename):
    file = open("../SiQD/" + filename, 'r')
    id = 0
    atom_dict = {}
    # Calculate this later
    bond_range = 3
    line_number = 0
    for line in file:
        if line_number == 0:
            n_atoms = line.strip()
            n_atoms = int(float(n_atoms))
        if line_number >= 2 and line_number < n_atoms + 2:
            values_list = line.split()
            for i in range(1, 4):
                values_list[i] = float(values_list[i])
            bound = hf.bond_checker(values_list, atom_dict, bond_range)
            atom_dict[id] = {
                            "x": values_list[1],
                            "y": values_list[2],
                            "z": values_list[3],
                            "element": values_list[0],
                            "bound": bound,
                            "type": "crystal"
                            }
            for item in bound:
                atom_dict[item]["bound"].append(id)
            id += 1
        line_number += 1
    return atom_dict


def builder(atom_dict, filename, standard):
    sites = tetra_sites(atom_dict)
    n_sites = len(sites)

    # Get ligand parameters
    ligand_types = []
    coverage_list = []
    n_ligands_list = []
    extension_list = []
    if standard == "n":
        while True:
            ligand_file = input("Filename of ligand to be added, or type files to see available ligands: ") + ".xyz"
            if ligand_file[:-4] == "files":
                lig_list = os.listdir("../Ligands")
                print()
                for ligs in lig_list:
                    print(ligs[:-4])
                print()
                ligand_file = input("Filename of ligand to be added: ") + ".xyz"
            ligand_types.append(ligand_file)
            extend = input("Extend ligand? Note: this replaces the last atom in the ligand file. y/n: ")
            if extend == "y":
                extension_list.append(input("Extend with: ") + ".xyz")
            else:
                extension_list.append(False)
            coverage = float(input("Coverage (fraction): "))
            coverage_list.append(coverage)
            n_ligands_list.append(round(coverage * n_sites))
            more = input("Add another type? y/n: ")
            if more == 'n':
                final_ligand = "H.xyz"
                extension_list.append(False)
                break
        space = 2
    else:
        ligand_types = ["C2H5.xyz"]
        coverage = 0.5
        coverage_list = [0.5]
        n_ligands_list.append(round(coverage * n_sites))
        final_ligand = "H.xyz"
        space = 1.8
        extension_list = [False]

    # Copy sites to use again later
    sites_copy = sites.copy()
    bridge_bool = input("Allow bridges (work in progress)? y/n: ")
    if bridge_bool == "y":
        atom_dict = bridges(atom_dict, sites)
    for i in range(len(ligand_types)):
        atom_dict = place_ligands(atom_dict, sites, n_ligands_list[i], ligand_types[i], extension_list[i], final_ligand, space, False)
    # Place final ligands at remaining sites
    n_ligands = len(sites)
    atom_dict = place_ligands(atom_dict, sites, n_ligands, final_ligand, extension_list[-1], final_ligand, space, False)

    # Write atoms to file
    dict2file(atom_dict, filename)

    # Replace ligands to create new QD
    while True:
        replacement_ligands = input("Replace ligands to create new quantum dot? y/n: ")
        if replacement_ligands == 'y':
            replace_ligands(atom_dict, ligand_types, sites_copy, n_ligands_list)
        else:
            break


# Replace ligands in QD with new ligands
def replace_ligands(atom_dict, ligand_types, sites, n_ligands_list):
    filename = input("Write to file: ")
    replacement_list = []
    rep_ext_list = []
    loc_dict = {"H": []}
    rep_dict = {}
    for ligand in ligand_types:
        replacement_list.append(input("Replace " + ligand + " with: ") + ".xyz")
        loc_dict[ligand] = {
                            "loc_id": [],
                            "replacement": replacement_list[-1]
                        }
        rep_dict[replacement_list[-1]] = {}
        extend = input("Extend ligand? y/n: ")
        if extend == "y":
            rep_ext_list.append(input("Extend with: ") + ".xyz")
        else:
            rep_ext_list.append(False)
    atom_del_list = []
    for atom, values in atom_dict.items():
        if values["type"] == 'ligand':
            atom_del_list.append(atom)
            if values["ligand_type"] != "H.xyz":
                if values["loc_id"] not in loc_dict[values["ligand_type"]]["loc_id"]:
                    loc_dict[values['ligand_type']]["loc_id"].append(values["loc_id"])
                    rep_dict[loc_dict[values['ligand_type']]["replacement"]][values["loc_id"]] = values["loc"]
    for item in atom_del_list:
        del atom_dict[item]
    print(rep_dict)
    i = 0

    for i in range(len(replacement_list)):
        atom_dict = place_ligands(atom_dict, sites, n_ligands_list[i], replacement_list[i], rep_ext_list[i], "H", 2, rep_dict)
    n_sites = len(sites)
    atom_dict = place_ligands(atom_dict, sites, n_sites, "H.xyz", False, "H", 2, False)

    dict2file(atom_dict, filename)


def dict2file(dict, filename):
    file = open("../Created_QD/" + filename + ".xyz", "w")
    file.write("        \n\n")
    for atom, values in dict.items():
        file.write(values['element'] + "\t" + str(values['x']) + "\t\t" +
                   str(values['y']) + "\t\t" + str(values['z']) + "\n")
    file.seek(0)
    file.write(str(len(dict)))
    file.close()
    print("\nQuantum Dot created :)")


# Calculate all sites for ligands for tetrahedral structure
def tetra_sites(dict):
    # Create a standard tetrahedron to determine sites
    x = 1 / math.sqrt(3)
    sites = {}
    # Determine possible sites atom by atom
    for id, values in dict.items():
        if len(values['bound']) < 4:
            sites[id] = {}
            tetrahedron = [[x, x, x], [x, -x, -x], [-x, x, -x], [-x, -x, x]]
            primary_xyz = [values['x'], values['y'], values['z']]
            connection_list = []
            # Create unit vectors for all connected atoms
            for secondary in values['bound']:
                sec_xyz_rel = [dict[secondary]['x'] - primary_xyz[0],
                               dict[secondary]['y'] - primary_xyz[1],
                               dict[secondary]['z'] - primary_xyz[2]]
                sec_xyz_rel = hf.normaliser(sec_xyz_rel)
                connection_list.append(sec_xyz_rel)

            # line up first vector
            axis = np.cross(tetrahedron[0], connection_list[0])
            if round(np.linalg.norm(axis), 2) == 0:
                axis = [0, x, -x]
                angle = 0
                if round(np.dot(tetrahedron[0], connection_list[0]), 2) == -1:
                    angle = math.pi
            else:
                angle = hf.angle_checker(tetrahedron[0], connection_list[0])
            rot_mat = hf.rotation_matrix(axis, angle)
            for k in range(4):
                tetrahedron[k] = list(np.dot(rot_mat, tetrahedron[k]))

            # line up second vector
            axis = tetrahedron[0]
            # Random rotation if single bond, and skip to next atom
            if len(connection_list) == 1:
                angle = random.random() * 2 * math.pi
                rot_mat = hf.rotation_matrix(axis, angle)
                for k in range(4):
                    tetrahedron[k] = list(np.around(np.array(np.dot(rot_mat, tetrahedron[k])), 4))
                sites[id] = {
                            'primary_xyz': primary_xyz,
                            'sites_xyz': [tetrahedron[1], tetrahedron[2], tetrahedron[3]]}
                continue
            # No rotation if already aligned
            elif round(np.dot(tetrahedron[1], connection_list[1]), 2) == 1:
                angle = 0
            # Get angle between planes spanned by first vector, second connection vector and second vector in tetrahedron
            else:
                normal_1 = np.cross(tetrahedron[0], connection_list[1])
                normal_2 = np.cross(tetrahedron[0], tetrahedron[1])
                angle = hf.angle_checker(normal_1, normal_2)

            # Make sure rotation is in right direction
            test_rot_mat = hf.rotation_matrix(axis, angle)
            testrahedron = tetrahedron.copy()
            for test in range(4):
                testrahedron[test] = list(np.dot(test_rot_mat, tetrahedron[test]))
            if round(np.dot(connection_list[1], testrahedron[1]), 1) != 1:
                angle = -angle

            # Rotate tetrahedron to correct orientation, add points without atom to sites
            rot_mat = hf.rotation_matrix(axis, angle)
            for k in range(4):
                tetrahedron[k] = list(np.around(np.array(np.dot(rot_mat, tetrahedron[k])), 4))
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


def bridges(atom_dict, sites):
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


def prep_ligand_file(atom_dict, ligand_type, loc_id, extension):
    path = "../Ligands/" + ligand_type
    file = open(path, 'r')
    id = 0
    line_number = 0
    lig = ligands.empty_lig()

    for line in file:
        if line_number >= 2:
            values_list = line.split()
            if float(values_list[1]) == float(values_list[2]) == float(values_list[3]) == 0:
                base_element = values_list[0]
                break
        line_number += 1

    initial_length = getattr(bond_dis, base_element)().distances[atom_dict[loc_id]["element"]]
    line_number = 0
    file.seek(0)
    for line in file:
        if line_number == 0:
            n_atoms = int(float(line.strip()))
        if line_number >= 2 and line_number < n_atoms + 2:
            values_list = line.split()
            for i in range(1, 4):
                values_list[i] = float(values_list[i])
            lig.atoms[id] = {
                            "x": float(values_list[1]),
                            "y": float(values_list[2]),
                            "z": float(values_list[3]) + initial_length,
                            "element": values_list[0]
                        }
            id += 1
        line_number += 1

    if extension is not False and ligand_type != "H.xyz":
        path_ext = "../Ligands/" + extension
        file_ext = open(path_ext, 'r')
        line_number = 0
        random_rotation = random.random() * 2 * math.pi * 0
        rep = lig.atoms[max(lig.atoms)]
        rep_coor = [rep["x"], rep["y"], rep["z"]]
        min_dist = math.inf
        ext = {}
        for atom, values in lig.atoms.items():
            dist = hf.distance_checker(rep_coor, [values['x'], values['y'], values['z']])
            if dist < min_dist and dist != 0:
                min_dist = dist
                closest_atom = values
        del lig.atoms[max(lig.atoms)]

        unit_v = hf.normaliser([c1 - c2 for c1, c2 in zip(rep_coor, [closest_atom['x'], closest_atom['y'], closest_atom['z']])])
        angle = hf.angle_checker(unit_v, [0, 0, 1])
        axis = np.cross(unit_v, [0, 0, 1])
        if round(np.linalg.norm(axis), 2) == 0:
            axis = [1, 0, 0]
        rot_mat1 = hf.rotation_matrix([0, 0, 1], random_rotation)
        rot_mat2 = hf.rotation_matrix(axis, angle)

        # test for correct angle
        test_v = np.dot(rot_mat2, unit_v)
        if round(np.dot(test_v, unit_v), 2) != 1:
            rot_mat2 = hf.rotation_matrix(axis, -angle)

        for line in file_ext:
            if line_number >= 2:
                values_list = line.split()
                ext_coor = [float(values_list[1]), float(values_list[2]), float(values_list[3])]
                if ext_coor[0] == ext_coor[1] == ext_coor[2] == 0:
                    initial_length_ext = getattr(bond_dis, values_list[0])().distances[closest_atom["element"]]
                    line_number = 0
                    file_ext.seek(0)
                    break
            line_number += 1
        for line in file_ext:
            if line_number == 0:
                n_atoms_ext = int(float(line.strip()))
            if line_number >= 2 and line_number < n_atoms_ext + 2:
                values_list = line.split()
                ext_coor = [float(values_list[1]), float(values_list[2]), float(values_list[3]) + initial_length_ext]
                ext_coor = np.dot(rot_mat1, ext_coor)
                ext_coor = np.dot(rot_mat2, ext_coor)
                ext[id] = {
                                "x": ext_coor[0] + closest_atom['x'],
                                "y": ext_coor[1] + closest_atom['y'],
                                "z": ext_coor[2] + closest_atom['z'],
                                "element": values_list[0]
                            }
            line_number += 1
            id += 1
        lig.atoms = {**lig.atoms, **ext}
    return lig


# Randomly choose a site to place ligand
def place_ligands(atom_dict, sites, n_ligands, ligand_type, extension, final_ligand, space, loc_dict):
    buffer = 0.5
    id = max(atom_dict) + 1
    original_ligand = ligand_type
    j = 0
    tried = []
    print("\nPlacing " + ligand_type)
    while j < n_ligands:
        # preventive in case of rounding errors
        if len(sites) == 0:
            break
        random_rotation = random.random() * 2 * math.pi
        extra_rotation = 0
        ligand_type = original_ligand

        # Randomly choose site, get relevant info, don't pick site that has already been tried
        sites_list = list(sites).copy()
        remaining_sites = [x for x in sites_list if x not in tried]
        if len(remaining_sites) == 0:
            print("\nUnable to place more ligands of type " + str(original_ligand) + ". Placed " + str(j) + " out of " + str(n_ligands) + " requested ligands. Continuing with other types.\n")
            break
        if loc_dict:
            if len(loc_dict[ligand_type]) == 0:
                break
            loc_id = random.choice(list(loc_dict[ligand_type]))
            loc = loc_dict[ligand_type][loc_id]
            del loc_dict[ligand_type][loc_id]
            loc_sites = sites[loc_id]['sites_xyz'].copy()
        else:
            loc_id = random.choice(remaining_sites)
            loc_sites = sites[loc_id]['sites_xyz'].copy()
            loc = random.choice(loc_sites)
        loc_primary_xyz = sites[loc_id]['primary_xyz']
        tried_loc = []

        # Loop over every site connected to chosen atom
        while loc not in tried_loc:
            lig = prep_ligand_file(atom_dict, ligand_type, loc_id, extension)
            # Get correct rotation for ligand relative to site
            axis = np.cross(loc, [0, 0, 1])
            # Prevent dividing by 0 when vectors are already lined up
            if np.linalg.norm(axis) == 0:
                axis = [1, 0, 0]
                angle = 0
            else:
                angle = -hf.angle_checker(loc, [0, 0, 1])
            rot_mat = hf.rotation_matrix(axis, angle)

            # Place ligand if possible, otherwise rotate
            while True:
                xyz_list = []
                temp_atom_dict = {}
                broken = False
                min_dist_lig = math.inf

                for atom in lig.atoms:
                    # Rotate every atom in ligand according to random_rotation
                    rotated_atom = np.dot(hf.rotation_matrix([0, 0, 1], random_rotation + extra_rotation), [lig.atoms[atom]['x'], lig.atoms[atom]['y'], lig.atoms[atom]['z']])
                    # Rotate in right direction
                    new_v = np.dot(rot_mat, rotated_atom)
                    atom_xyz = [c1 + c2 for c1, c2 in zip(new_v, loc_primary_xyz)]
                    atom_element = lig.atoms[atom]['element']
                    temp_atom_dict[id] = {
                                    "x": atom_xyz[0],
                                    "y": atom_xyz[1],
                                    "z": atom_xyz[2],
                                    "element": atom_element,
                                    "type": "ligand",
                                    "ligand_type": ligand_type,
                                    "loc_id": loc_id,
                                    "loc": loc
                                    }
                    xyz_list.append(atom_xyz)
                    id += 1
                    # Check for overlap
                    if ligand_type != final_ligand:
                        for test_atom, values in atom_dict.items():
                            try:
                                space2 = getattr(bond_dis, atom_element)().distances[values["element"]] + buffer
                            except KeyError:
                                space2 = float(input("Bond length between " + str(atom_element) + " and " + str(values["element"]) + " not available. Enter manually: ")) + buffer
                            if test_atom != loc_id:
                                dist = hf.distance_checker([values['x'], values['y'], values['z']], xyz_list[-1])
                                if dist < space2:
                                    broken = True
                                    break
                            if broken:
                                break
                    if ligand_type == final_ligand:
                        new_pos = new_v.copy()
                        changed = False
                        while True:
                            min_dist_lig = math.inf
                            for test_lig, values in atom_dict.items():
                                if values["type"] == "ligand":
                                    dist = hf.distance_checker([values['x'], values['y'], values['z']], [c1 + c2 for c1, c2 in zip(new_pos, loc_primary_xyz)])
                                    bond_len_loc = getattr(bond_dis, atom_element)().distances[values["element"]]
                                    if dist < bond_len_loc + 0.3:
                                        if dist < min_dist_lig:
                                            min_dist_lig = dist
                                            min_dist_lig_coor = [values['x'], values['y'], values['z']]
                                            changed = True
                            if min_dist_lig == math.inf:
                                if changed:
                                    xyz_list = [c1 + c2 for c1, c2 in zip(new_pos, loc_primary_xyz)]
                                    temp_atom_dict[id-1] = {
                                                    "x": xyz_list[0],
                                                    "y": xyz_list[1],
                                                    "z": xyz_list[2],
                                                    "element": "H",
                                                    "type": "ligand",
                                                    "ligand_type": ligand_type,
                                                    "loc_id": loc_id,
                                                    "loc": loc
                                                    }
                                break
                            vec_closest_element = [c1 - c2 for c1, c2 in zip(min_dist_lig_coor, loc_primary_xyz)]
                            axis = np.cross(vec_closest_element, new_pos)
                            angle = 0.1
                            rot_mat = hf.rotation_matrix(axis, angle)
                            new_pos = np.dot(rot_mat, new_pos.copy())

                    if broken:
                        break
                # Check for ligand collisions except when last ligand is added
                if ligand_type != final_ligand:
                    # Rotate if there is a collision
                    if broken and extra_rotation < math.pi * 2:
                        extra_rotation += 0.2
                        continue
                    # break if rotated fully
                    elif extra_rotation > math.pi * 2:
                        extra_rotation = 0
                        tried.append(loc_id)
                        tried_loc.append(loc)
                        print("fail")
                        break
                    # Merge dicts if no collision
                    elif not broken:
                        atom_dict = {**atom_dict, **temp_atom_dict}
                        if j % 10 == 0:
                            print(str(j) + " ligands added")
                        if loc_id in sites:
                            del sites[loc_id]
                        tried_loc = [loc]
                        ligand_type = final_ligand
                        j += 1
                        break
                else:
                    atom_dict = {**atom_dict, **temp_atom_dict}
                    if loc_id in sites:
                        del sites[loc_id]
                    tried_loc.append(loc)
                    break
            # Add final ligand type to remaining sites at chosen atom
            loc = loc_sites[(loc_sites.index(loc) + 1) % len(loc_sites)]

    return atom_dict


def place_bridge_ligands(atom_dict, sites, bridge_dict, bridge_ligand):
    id = max(atom_dict) + 1
    for bridge, values in bridge_dict.items():
        crys_val = atom_dict[values["connected"][0]]
        crys_vec = [crys_val["x"], crys_val["y"], crys_val["z"]]
        random_rotation = math.pi/4
        # Determine height of ligand. If bonding distance is too short, set to 0.
        # This is only correct for crystal with one kind of atom
        try:
            initial_length = math.sqrt(getattr(bond_dis, "C")().distances["Si"]**2 - (3.84/2)**2)
        except ValueError:
            initial_length = 0
        lig_object = getattr(ligands, bridge_ligand)
        lig = lig_object(initial_length)
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
            uhh = [atom_dict[values["connected"][0]]["x"], atom_dict[values["connected"][0]]["y"], atom_dict[values["connected"][0]]["z"]]
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
                                "ligand_type": ligand_type
                                }
                xyz_list.append(atom_xyz)
                id += 1
                for test_atom, values2 in atom_dict.items():
                    space2 = getattr(bond_dis, atom_element)().distances[values2["element"]]
                    # Stay further away from the crystal atoms
                    if values2['type'] == "crystal":
                        space2 = space2 + 0.25
                    if test_atom not in values["connected"]:
                        dist = hf.distance_checker([values2['x'], values2['y'], values2['z']], xyz_list[-1])
                        if dist < space2 + 0.25:
                            break
            atom_dict = {**atom_dict, **temp_atom_dict}
            break
    return atom_dict


if __name__ == "__main__":
    standard = input("Standard? y/n: ")
    if standard == "n":
        build = input("Create new crystal (y) or use existing file (n)?: ")
        if build == "y":
            structure = input("Specify structure type: ")
            a = float(input("Specify lattice constant (in Ångström): "))
            atom_a = input("Element for first element type: ")
            atom_b = input("Element for second element type: ")
            diameter = float(input("Diameter of quantum dot (in unit cells): "))
            filename = input("Write to file: ")
            atom_dict = crystal_builder(structure, a, atom_a, atom_b, diameter, filename)
        else:
            crystal_file = input("Crystal file to use (don't write the file extension): ") + ".xyz"
            filename = input("Write to file: ")
            atom_dict = crystal_reader(crystal_file)
    else:
        structure = "zns"
        a = 5.431
        atom_a = "Si"
        atom_b = "Si"
        diameter = 4
        filename = "test"
        atom_dict = crystal_builder(structure, a, atom_a, atom_b, diameter, filename)
    builder(atom_dict, filename, standard)
