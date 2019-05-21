import bonding_distances as bond_dis
import math
import numpy as np


class COOH(object):
    def __init__(self, bound_to):
        OH_length = bond_dis.H().distances["O"]
        init_length = bond_dis.C().distances[bound_to]
        self.atoms = {
            "1e": {
                "element": "C",
                "x": 0,
                "y": 0,
                "z": init_length
            },
            "2e": {
                "element": "O",
                "x": - 1.23 * math.cos(np.deg2rad(34)),
                "y": 0,
                "z": init_length + 1.23 * math.sin(np.deg2rad(34))
            },
            "3e": {
                "element": "O",
                "x": 1.32 * math.cos(np.deg2rad(21)),
                "y": 0,
                "z": init_length + 1.32 * math.sin(np.deg2rad(21))
            },
            "4e": {
                "element": "H",
                "x": 1.32 * math.cos(np.deg2rad(21)) + OH_length * math.sin(np.deg2rad(51)),
                "y": 0,
                "z": init_length + 1.32 * math.sin(np.deg2rad(21)) - OH_length * math.cos(np.deg2rad(51))
            }
        }


class Na(object):
    def __init__(self, bound_to):
        init_length = bond_dis.Na().distances[bound_to]
        self.atoms = {
            "1e": {
                "element": "Na",
                "x": 0,
                "y": 0,
                "z": init_length
            }
        }


class NH2(object):
    def __init__(self, bound_to):
        init_length = bond_dis.N().distances[bound_to]
        NH_length = bond_dis.N().distances["H"]
        self.atoms = {
            "1e": {
                "element": "N",
                "x": 0,
                "y": 0,
                "z": init_length
            },
            "2e": {
                "element": "H",
                "x": NH_length * math.sqrt(8/9),
                "y": 0,
                "z": init_length + 1/3 * NH_length
            },
            "3e": {
                "element": "H",
                "x": -NH_length * math.sqrt(2/9),
                "y": NH_length * math.sqrt(2/3),
                "z": init_length + 1/3 * NH_length
            }
        }
