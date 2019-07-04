#!/usr/bin/python
import sys

(len(sys.argv)>1) or sys.exit("Please specify input file! F.e. ./charge.py cp2k.out")

#print 'Loading file:', sys.argv[1]
with open(sys.argv[1], "r") as ins:
    readlines = 0
    data = {}
    for line in ins:
        if ( readlines == 0 and line.find('#  Atom') > 0 ):
            readlines = 1
        elif ( readlines > 0 and line.find('# Total') > 0 ):
            readlines = 0
        elif ( readlines > 0 ):
            variable = line[10:17].strip(' ')
            value = line[42:79].strip(' ')
            if not variable in data:
                data[variable] = 0
            data[variable] += float(value) 

for key in data:
    print (key, '=', data[key])

#print 'Finish!'
    
