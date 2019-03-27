#!/bin/env python

from sys import argv
import os,math,numpy

linevalue=[]
ave=[]

f = open("end",'r')
w = open("end-ave",'w')

def deal(f):
	for line in f.readlines():
		strs = line.strip().split()
		if strs != []:
			linevalue=map(float,strs)
			avevalue=numpy.mean(linevalue)
			ave.append(float(avevalue))

def out(w):
	average=numpy.mean(ave)
	st=numpy.std(ave)
	w.write("%.2f, %.2f" %(average,st))

deal(f)
out(w)

f.close