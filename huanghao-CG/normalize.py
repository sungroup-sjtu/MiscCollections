inf = open ("pacf.log",'r')
outf = open("nacf.log",'w')
c = 0
for line in inf:
	c += 1
	string = line.strip().split()
	if c == 1:
		n = float(string[1])
	a = float(string[1])/n
	outf.write("%s %f\n" %(string[0],a))

inf.close()
outf.close()		
