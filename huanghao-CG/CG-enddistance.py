#!/bin/env python
from sys import argv
import os,math

ATOMSbegin=[]
ATOMSend=[]
enddistance=[]
CN=20

def read_frame(f):
    line = f.readline()
    count = 0
    for line in f.readlines():
        strs = line.strip().split()
        length = len(strs)
        if length == 5 and float(strs[1]) == 3:
            if count % 2 == 0:
                ATOMSbegin.append((float(strs[2]), float(strs[3]), float(strs[4])))
            elif count % 2 == 1:
                ATOMSend.append((float(strs[2]), float(strs[3]), float(strs[4])))
            count = count + 1

def result():
    average=0
    for i in range(len(ATOMSbegin)):
        result=((ATOMSend[i][0]-ATOMSbegin[i][0])**2+(ATOMSend[i][1]- \
        ATOMSbegin[i][1])**2+(ATOMSend[i][2]-ATOMSbegin[i][2])**2)**0.5
        results= "%.3f" %result
        enddistance.append(results)
        
def convert(out):
	for i in range(len(enddistance)/CN):
		c= ' '.join(map(str, enddistance[i*CN:min((i+1)*CN, len(enddistance))])) + '\n'
		out.write(c)

f = open(argv[1], 'r') 
w = open("end-end",'w')
read_frame(f)
result()
convert(w)
f.close()
