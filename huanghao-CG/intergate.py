inf = open("pacf.log",'r')
outf = open("vis.log",'w')
a = 0
step = 0.1
index =0.0003227
for line in inf:

	string=line.strip().split()
	a += float(string[1])*step	
	v = a*index
	outf.write("%s %f\n" %(string[0],v))

inf.close()
outf.close()

