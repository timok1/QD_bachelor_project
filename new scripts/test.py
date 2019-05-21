f = open("Si235.xyz")
lines = f.readlines()
f.close()
f = open("Si235.xyz", 'w')
for line in lines:
    f.write(line[1:])
f.close()
