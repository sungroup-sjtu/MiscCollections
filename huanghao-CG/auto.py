#!/usr/bin/env python 
#coding=utf-8

import os, sys

e_list = ['B-5','B-10','B-20','B-50']
d_list = ['200K','300K','400K','500K']

for y1 in e_list:
    os.chdir(y2)
    for y2 in e_list:
    	os.chdir(y2)
    	os.system("qsub a.sh")
    	os.chdir("..")
    os.chdir("..")
