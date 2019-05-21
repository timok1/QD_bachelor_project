# -*- coding: utf-8 -*-
"""
convert .xyz to .bas
"""
#! /usr/bin/python
import sys
import math
from elements import ELEMENTS

inname =sys.argv[1];    infile = open(inname , "r" )
outname=sys.argv[2];   outfile = open(outname, "w" )

ELDICT = dict( ( ELEMENTS[i][1], ELEMENTS[i] ) for i in range(len(ELEMENTS)) )
#print " ELDICT : ", ELDICT

n=int(infile.readline())
outfile.write (" %i \n" % n)
infile.readline()
for i in range(n):
	l=infile.readline().split()
	rec = ( l[0]  ,  float(l[1])  , float(l[2])  , float(l[3]) )
	#print i, rec[0], rec[1], rec[2], rec[3], ELDICT[ rec[0] ]
	outfile.write(  " %s %f %f %f \n"  %    (   ELDICT[ rec[0] ][0] , rec[1], rec[2], rec[3]    )    )
infile.close()
outfile.close()