import math

file = open("rad_rate_sorted.txt", 'r')

i = 0
rad = []
energy = []
for line in file:
    for x in line.split(","):
        x = x.strip("\n")
        x = float(x)
        if i % 2 == 0:
            energy.append(x)
        else:
            rad.append(x)
        i += 1
print(i)

temp = 300
teller = 0
noemer = 0

for j, k in zip(rad, energy):
    boltz = math.exp(- k / (8.6173 * 10**-5 * temp))
    teller += (j * boltz)
    noemer += boltz

print(teller/noemer)
