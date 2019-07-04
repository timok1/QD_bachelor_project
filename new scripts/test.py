import math
import numpy as np

dict = {
    "H": 0,
    "C": 0,
    "Si": 0
}

file = open("relaxed.xyz", 'r')
i = 0
for line in file:
    if i > 1:
        element = line.split()[0]
        dict[element] += 1
    i += 1
print(dict)
