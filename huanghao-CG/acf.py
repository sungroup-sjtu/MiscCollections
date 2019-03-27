import os
import sys
input = sys.argv[1]
timestep = 0.1
#pxx = []
#pyy = []
#pzz = []
pxy = []
pxz = []
pyz = []
acf = []
inf = open(input,'r')
point_data = 0
for line in inf:
    point_data += 1
    string = line.strip().split()
#    pxx.append(float(string[4]))
#    pyy.append(float(string[5]))
#    pzz.append(float(string[6]))
    pxy.append(0.5*(float(string[2])+float(string[4])))
    pxz.append(0.5*(float(string[3])+float(string[7])))
    pyz.append(0.5*(float(string[6])+float(string[8])))
#    pxz.append(float(string[8]))
#    pyz.append(float(string[9]))
inf.close()
normal=0
half_point = (point_data+1)/2
for i in range(0,point_data):
    normal += pxy[i]*pxy[i]
normal = normal/half_point    
outf = open ("pacf.log",'w')
for i in range(0,100000):
    a= 0.0
    b= 0.0
    c= 0.0
    d= 0.0
    e= 0.0
    f= 0.0
    k = 0.0
    for j in range(0,point_data):
        if i+j < point_data:
            k += 1.0
       #    a += 0.25*(pxx[j]-pyy[j])*(pxx[j+i]-pyy[j+i]) 
       #    b += 0.25*(pxx[j]-pzz[j])*(pxx[j+i]-pzz[j+i]) 
       #    c += 0.25*(pyy[j]-pzz[j])*(pyy[j+i]-pzz[j+i]) 
            d += pxy[j]*pxy[j+i] 
            e += pxz[j]*pxz[j+i] 
            f += pyz[j]*pyz[j+i] 
        
    
    #ave_a = a/k
    #ave_b = b/k
    #ave_c = c/k
    ave_d = d/k
    #ave = d/k
    ave_e = e/k
    ave_f = f/k
    ave = (ave_d+ave_e+ave_f)/3.0
    #acf.append(ave)
    #outf.write("%s %s %s %s %s %s %s\n" %(i*timestep,ave_a,ave_b,ave_c,ave_d,ave_e,ave_f))
    outf.write("%s %s\n" %(i*timestep,ave))
outf.close()
