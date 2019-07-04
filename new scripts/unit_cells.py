import numpy as np

class pc(object):
    """Creates a primitive cubic unit cell with lattice constant a and atom
    atom_a"""

    def __init__(self, a, atom_a):
        self.a = a
        self.atom_a = atom_a
        self.atom_xyz = []
        for x in np.arange(-.5, 1, 1):
            for y in np.arange(-.5, 1, 1):
                for z in np.arange(-.5, 1, 1):
                    self.atom_xyz.append([atom_a, x*a, y*a, z*a])


class fcc(object):
    """Creates a face centered cubic unit cell with lattice constant a and atom
    atom_a"""

    def __init__(self, a, atom_a):
        self.a = a
        self.atom_a = atom_a
        self.atom_xyz = []
        # Atoms on corners of unit cell
        for x in np.arange(-.5, 1, 1):
            for y in np.arange(-.5, 1, 1):
                for z in np.arange(-.5, 1, 1):
                    self.atom_xyz.append([atom_a, x*a, y*a, z*a])

        # Atoms on faces
        for x in np.arange(-0.5, 1, 0.5):
            for y in np.arange(-0.5, 1, 0.5):
                for z in np.arange(-0.5, 1, 0.5):
                    if abs(x) + abs(y) + abs(z) == 0.5:
                        self.atom_xyz.append([atom_a, x*a, y*a, z*a])


class zns(object):
    """Creates a zincblende or diamond cubic unit cell with lattice constant a
    and elements atom_a and atom_b"""
    def __init__(self, a, atom_a, atom_b):
        self.a = a
        self.atom_a = atom_a
        self.atom_xyz = []
        xyz_pre = []
        # Find coordinates according to x = y = z (mod 2), and x + y + z = 0 or 1 (mod 4)
        for x in range(5):
            for y in range(5):
                for z in range(5):
                    if ((x + y + z) % 4 == 0 or (x + y + z) % 4 == 1) and (x % 2 == y % 2 == z % 2):
                        xyz_pre.append([x, y, z])
                        if x == y == z == 0:
                            for i in range(2):
                                for j in range(2):
                                    for k in range(2):
                                        if not x == y == z == 0:
                                            xyz_pre.append([x + 4 * i, y + 4 * j, z + 4 * k])
                        elif x == 0:
                            xyz_pre.append([x + 4, y, z])
                        elif y == 0:
                            xyz_pre.append([x, y + 4, z])
                        elif z == 0:
                            xyz_pre.append([x, y, z + 4])
        # Find coordinates of atom_b
        for atom in xyz_pre:
            if 3 in atom:
                xyz = [atom_b]
            else:
                xyz = [atom_a]
            for coordinate in atom:
                xyz.append((coordinate / 4 - 0.5) * a)
            self.atom_xyz.append(xyz)
