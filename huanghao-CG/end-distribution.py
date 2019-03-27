#!/bin/env python
from sys import argv
import numpy as np
import math

min=0
interv=100
value=[]
min=0
lines=[]

def read_end(f):
    global data, count
    for line in f.readlines():
        lines.append(line)
        strs = line.strip().split()
        for i in range(len(strs)):
            value.append(float(strs[i]))
#   maxvalue=math.ceil(max(value))+1.0
    maxvalue=20.0
    global data, count 
    data  = np.arange(0.5,maxvalue+0.5,(maxvalue-min)/interv)
    count = np.zeros(len(data))
    
    for i in range(len(value)):
        for j in range(len(data)):
        	if value[i] >= j*((maxvalue-min)/interv) and value[i] < (j+1)*((maxvalue-min)/interv):
        	    count[j]=count[j]+1.0

def convert(out):
	for i in range(len(data)):
		out.write(str(data[i]) + " " + str(count[i]/len(value)) + "\n")

f = open("end", 'r') 
w = open("end-distribution",'w')
read_end(f)
convert(w)
f.close()
