class C(object):
    def __init__(self):
        self.distances = {
                    "H":	1.09,
                    "Be":	1.93,
                    "Mg":	2.07,
                    "B": 	1.56,
                    "Al":	2.24,
                    "In":	2.16,
                    "C":	1.53,
                    "Si":	1.86,
                    "Sn":	2.14,
                    "Pb":	2.29,
                    "N": 	1.47,
                    "P":	1.87,
                    "As":	1.98,
                    "Sb":	2.20,
                    "Bi":	2.30,
                    "O":	1.43,
                    "S":	2.18,
                    "Cr":	1.92,
                    "Se":	2.35,
                    "Te":	2.05,
                    "Mo":	2.08,
                    "W":	2.06,
                    "F":	1.34,
                    "Cl":	1.76,
                    "Br":	1.94,
                    "I":	2.13,
                    "Cd":   2.13,
                    "Na":   2.24,
                    "Ge":   1.94
                    }

class H(object):
    def __init__(self):
        self.distances = {
                    "H":	0.74,
                    "C":	1.09,
                    "Si":   1.48,
                    "F":    0.97,
                    "O":    0.96,
                    "N":    1.02,
                    "Na":   1.97,
                    "Br":   1.45,
                    "Ge":   1.55
                    }

class O(object):
    def __init__(self):
        self.distances = {
                    "O":    1.45,
                    "Si":   1.60,
                    "C":    1.43,
                    "H":    0.96,
                    "N":    1.45,
                    "F":    1.41,
                    "Na":   2.04,
                    "Br":   1.80,
                    "Ge":   1.79
                    }

class F(object):
    def __init__(self):
        self.distances = {
                    "H":    0.97,
                    "C":    1.34,
                    "F":    1.41,
                    "Si":   1.59,
                    "O":    1.41,
                    "N":    1.37,
                    "Na":   1.98,
                    "Br":   1.85,
                    "Ge":   1.76
                    }

class Si(object):
    def __init__(self):
        self.distances = {
                    "O":    1.60,
                    "H":    1.48,
                    "C":    1.86,
                    "F":    1.59,
                    "Si":   2.35,
                    "N":    1.70,
                    "Na":   2.63,
                    "Br":   2.21,
                    "Ge":   2.31
                    }

class N(object):
    def __init__(self):
        self.distances = {
                    "C":    1.47,
                    "H":    1.02,
                    "N":    1.44,
                    "Si":   1.70,
                    "O":    1.45,
                    "F":    1.37,
                    "Na":   2.13,
                    "Br":   1.85,
                    "Ge":   1.87
        }

class Na(object):
    def __init__(self):
        self.distances = {
                    "H":    1.97,
                    "C":    2.24,
                    "N":    2.13,
                    "O":    2.04,
                    "F":    1.98,
                    "Na":   3.08,
                    "Si":   2.63,
                    "Br":   2.57,
                    "Ge":   2.72
        }

class Br(object):
    def __init__(self):
        self.distances = {
                    "H":    1.45,
                    "C":    1.94,
                    "N":    1.85,
                    "O":    1.80,
                    "F":    1.85,
                    "Na":   2.57,
                    "Si":   2.21,
                    "Br":   2.30,
                    "Ge":   2.33
        }

class Ge(object):
    def __init__(self):
        self.distances = {
                    "H":    1.55,
                    "C":    1.94,
                    "N":    1.87,
                    "O":    1.79,
                    "F":    1.76,
                    "Na":   2.72,
                    "Si":   2.31,
                    "Br":   2.33,
                    "Ge":   2.39
        }
