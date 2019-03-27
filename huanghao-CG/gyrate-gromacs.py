#!/usr/bin/python

from sys import argv
import os,math,numpy

data=[]


f = open("gyrate.xvg", 'r')
w = open("gyrate",'w')

def read_frame(f):
    line = f.readline()    
    while line:
        if line.startswith("@ s3 legend"):
            while line:
                line = f.readline()
                strs = line.strip().split()
                if strs == []:
                	break
                else:
                	data.append(strs[1])
        line = f.readline()

def out(w):
    for i in data:
        w.write(i)
        w.write("\n")

read_frame(f)
out(w)

f.close()