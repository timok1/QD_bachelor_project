import bonding_distances as bond_dis
import math


class empty_lig(object):
    def __init__(self):
        self.atoms = {}

class H(object):
    def __init__(self, bound_to):
        self.atoms = {
            1: {
                "element": "H",
                "x": 0,
                "y": 0,
                "z": bond_dis.H().distances[bound_to]
            }
        }


class CH2(object):
    def __init__(self, initial_length):
        CH_length = bond_dis.C().distances["H"]
        self.atoms = {
            1: {
                "element": "C",
                "x": 0,
                "y": 0,
                "z": initial_length
            },
            2: {
                "element": "H",
                "x": CH_length * 0.5 * math.sqrt(2),
                "y": 0,
                "z": initial_length + CH_length * 0.5 * math.sqrt(2),
            },
            3: {
                "element": "H",
                "x": -CH_length * 0.5 * math.sqrt(2),
                "y": 0,
                "z": initial_length + CH_length * 0.5 * math.sqrt(2)
            },
        }
        self.replacement_atom = 3


class CH3(object):
    # Follows tetrahedral shape
    def __init__(self, bound_to):
        CH_length = bond_dis.C().distances["H"]
        initial_length = bond_dis.C().distances[bound_to]
        self.base_atom = [0, 0, initial_length]
        self.base_element = "C"
        self.atoms = {
            1: {
                "element": "C",
                "x": 0,
                "y": 0,
                "z": initial_length
            },
            2: {
                "element": "H",
                "x": CH_length * math.sqrt(8/9),
                "y": 0,
                "z": initial_length + 1/3 * CH_length
            },
            3: {
                "element": "H",
                "x": -CH_length * math.sqrt(2/9),
                "y": CH_length * math.sqrt(2/3),
                "z": initial_length + 1/3 * CH_length
            },
            4: {
                "element": "H",
                "x": -CH_length * math.sqrt(2/9),
                "y": -CH_length * math.sqrt(2/3),
                "z": initial_length + 1/3 * CH_length
            }
        }
        self.replacement_atom = 4


class C2H4(object):
    def __init__(self, bound_to):
        CH_length = bond_dis.C().distances["H"]
        CC_length = bond_dis.C().distances["C"]
        initial_length = bond_dis.C().distances[bound_to]
        del self.atoms[1]
        self.atoms = CH3(bound_to).atoms
        self.base_atom = [-CC_length * math.sqrt(2/9), -CC_length * math.sqrt(2/3), initial_length + 1/3 * CC_length]
        self.base_element = "C"
        self.atoms[4] = {
                "element": "C",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2]
            }
        self.atoms[5] = {
                "element": "H",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2] + CH_length
            }
        self.atoms[6] = {
                "element": "H",
                "x": self.base_atom[0] - CH_length * math.sqrt(8/9),
                "y": self.base_atom[1],
                "z": self.base_atom[2] - (1/3 * CH_length)
            }
        self.atoms[7] = {
                "element": "H",
                "x": self.base_atom[0] + CH_length * math.sqrt(2/9),
                "y": self.base_atom[1] - CH_length * math.sqrt(2/3),
                "z": self.base_atom[2] - (1/3 * CH_length)
            }


class C2H5(object):
    def __init__(self, bound_to):
        CH_length = bond_dis.C().distances["H"]
        CC_length = bond_dis.C().distances["C"]
        initial_length = bond_dis.C().distances[bound_to]
        self.atoms = CH3(bound_to).atoms
        self.base_atom = [-CC_length * math.sqrt(2/9), -CC_length * math.sqrt(2/3), initial_length + 1/3 * CC_length]
        self.base_element = "C"
        self.atoms[4] = {
                "element": "C",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2]
            }
        self.atoms[5] = {
                "element": "H",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2] + CH_length
            }

        self.atoms[6] = {
                "element": "H",
                "x": self.base_atom[0] - CH_length * math.sqrt(8/9),
                "y": self.base_atom[1],
                "z": self.base_atom[2] - (1/3 * CH_length)
            }

        self.atoms[7] = {
                "element": "H",
                "x": self.base_atom[0] + CH_length * math.sqrt(2/9),
                "y": self.base_atom[1] - CH_length * math.sqrt(2/3),
                "z": self.base_atom[2] - (1/3 * CH_length)
            }
        self.replacement_atom = 5


class C3H7(object):
    def __init__(self, bound_to):
        CH_length = bond_dis.C().distances["H"]
        CC_length = bond_dis.C().distances["C"]
        C2H5_basis = C2H5(bound_to)

        self.atoms = C2H5_basis.atoms
        self.base_atom = C2H5_basis.base_atom
        self.base_atom[2] = C2H5_basis.base_atom[2] + CC_length
        self.base_element = "C"

        self.atoms[5] = {
                "element": "C",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2]
            }

        self.atoms[8] = {
                "element": "H",
                "x": self.base_atom[0] + CH_length * math.sqrt(8/9),
                "y": self.base_atom[1],
                "z": self.base_atom[2] + 1/3 * CH_length

            }

        self.atoms[9] = {
                "element": "H",
                "x": self.base_atom[0] - CH_length * math.sqrt(2/9),
                "y": self.base_atom[1] + CH_length * math.sqrt(2/3),
                "z": self.base_atom[2] + 1/3 * CH_length

            }

        self.atoms[10] = {
                "element": "H",
                "x": self.base_atom[0] - CH_length * math.sqrt(2/9),
                "y": self.base_atom[1] - CH_length * math.sqrt(2/3),
                "z": self.base_atom[2] + 1/3 * CH_length
            }

        self.replacement_atom = 10


class C4H9(object):
    def __init__(self, bound_to):
        CH_length = bond_dis.C().distances["H"]
        CC_length = bond_dis.C().distances["C"]
        C3H7_basis = C3H7(bound_to)
        self.atoms = C3H7_basis.atoms
        self.base_atom = [C3H7_basis.base_atom[0] - CC_length * math.sqrt(2/9),
                            C3H7_basis.base_atom[1] - CC_length * math.sqrt(2/3),
                            C3H7_basis.base_atom[2] + + 1/3 * CC_length]
        self.base_element = "C"

        self.atoms[10] = {
                "element": "C",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2]
            }

        self.atoms[11] = {
                "element": "H",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2] + CH_length
            }

        self.atoms[12] = {
                "element": "H",
                "x": self.base_atom[0] - CH_length * math.sqrt(8/9),
                "y": self.base_atom[1],
                "z": self.base_atom[2] - (1/3 * CH_length)
            }

        self.atoms[13] = {
                "element": "H",
                "x": self.base_atom[0] + CH_length * math.sqrt(2/9),
                "y": self.base_atom[1] - CH_length * math.sqrt(2/3),
                "z": self.base_atom[2] - (1/3 * CH_length)
            }

        self.replacement_atom = 11


class CF3(object):
    def __init__(self, bound_to):
        CF_length = bond_dis.C().distances["F"]
        initial_length = bond_dis.C().distances[bound_to]
        self.base_atom = [0, 0, initial_length]
        self.base_element = "C"
        self.atoms = {
            "1": {
                "element": "C",
                "x": 0,
                "y": 0,
                "z": initial_length
            },
            2: {
                "element": "F",
                "x": CF_length * math.sqrt(8/9),
                "y": 0,
                "z": initial_length + 1/3 * CF_length
            },
            3: {
                "element": "F",
                "x": -CF_length * math.sqrt(2/9),
                "y": CF_length * math.sqrt(2/3),
                "z": initial_length + 1/3 * CF_length
            },
            4: {
                "element": "F",
                "x": -CF_length * math.sqrt(2/9),
                "y": -CF_length * math.sqrt(2/3),
                "z": initial_length + 1/3 * CF_length
            }
        }
        self.replacement_atom = 4

class C2F5(object):
    def __init__(self, bound_to):
        CF_length = bond_dis.C().distances["F"]
        CC_length = bond_dis.C().distances["C"]
        initial_length = bond_dis.C().distances[bound_to]
        self.atoms = CF3(bound_to).atoms
        self.base_atom = [-CC_length * math.sqrt(2/9), -CC_length * math.sqrt(2/3), initial_length + 1/3 * CC_length]
        self.base_element = "C"
        self.atoms[4] = {
                "element": "C",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2]
            }
        self.atoms[5] = {
                "element": "F",
                "x": self.base_atom[0],
                "y": self.base_atom[1],
                "z": self.base_atom[2] + CF_length
            }

        self.atoms[6] = {
                "element": "F",
                "x": self.base_atom[0] - CF_length * math.sqrt(8/9),
                "y": self.base_atom[1],
                "z": self.base_atom[2] - (1/3 * CF_length)
            }

        self.atoms[7] = {
                "element": "F",
                "x": self.base_atom[0] + CF_length * math.sqrt(2/9),
                "y": self.base_atom[1] - CF_length * math.sqrt(2/3),
                "z": self.base_atom[2] - (1/3 * CF_length)
            }
        self.replacement_atom = 5
