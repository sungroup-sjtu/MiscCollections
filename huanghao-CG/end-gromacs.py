#!/usr/bin/python

from sys import argv
import os,math,numpy

data=[]


f = open("dist.xvg", 'r')
w = open("end",'w')

def read_frame(f):
    line = f.readline()    
    while line:
        if line.startswith("@TYPE xy"):
            while line:
                line = f.readline()
                strs = line.strip().split()
                if strs == []:
                	break
                else:
                	data.append(strs[1:len(strs)])
        line = f.readline()


def out(w):
    for i in data:
        w.write(" ".join(i))
        w.write("\n")

read_frame(f)
out(w)

f.close()