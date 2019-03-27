#!/usr/bin/env python
# coding=utf-8

import os,sys,time
import numpy as np

dl = 0.1
l_list = np.arange(0,1.01,dl)

for l in l_list:
    lam = '%.2f' %l
    dirname = 'bar-%s' %lam
    if not os.path.exists('%s' %dirname):
        os.mkdir('%s' %dirname)
    if l == 0:
        cmd='sed \'s/%%lambda%%/%s/;s/%%dl%%/%.2f/;/v_ndl/d;/c_cFepndl/d\' PDMS.in > %s/PDMS.in' %(lam, dl, dirname)
    elif l == 1:
        cmd='sed \'s/%%lambda%%/%s/;s/%%dl%%/%.2f/;/v_dl/d;/c_cFepdl/d\' PDMS.in > %s/PDMS.in' %(lam, dl, dirname)
    else:
        cmd = 'sed \'s/%%lambda%%/%s/;s/%%dl%%/%.2f/;\' PDMS.in > %s/PDMS.in' %(lam, dl, dirname)
    print cmd
    os.system(cmd)

	## write .sh
    os.chdir(dirname)

    fname = dirname + ".sh"
    fout = open(fname, 'w')
    fout.write("#!/bin/bash\n")
    #fout.write("#PBS -q fast\n") 
    fout.write("#PBS -l nodes=1:ppn=8,walltime=120:00:00\n\n")
    fout.write("cd " + os.getcwd() + "\n\n")
    fout.write("mpirun -np 8 lmp-soft < PDMS.in > output.log")
    fout.close()

    ## qsub
    os.system("qsub " + fname)
    time.sleep(1)

    os.chdir("..")
