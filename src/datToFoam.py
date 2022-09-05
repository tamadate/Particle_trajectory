import math
import numpy as np
import itertools

N=3761572	
E=4046987
F=11794710
NEI=11691478

r=np.zeros((N,3))
param=np.zeros((E,6))
cell=np.zeros((E,8),dtype=int)

rOF=np.zeros((N,3))
ownerOF=np.zeros(F,dtype=int)
neighbourOF=np.zeros(NEI,dtype=int)

readFile="solution_03Torr.dat"
print ("start to get data from Tecplot file: " + str(readFile))
f=open(readFile,"r")
data=f.readlines()
f.close()


################Get point (x, y, z)##################
print "start to get point data, (x, y, z)"
start=4

for i in np.arange(3):
	for j in np.arange(N):
		r[j][i]=data[start+j]
	start+=N

################Get cell data (p, T, v, etc...)##################
print "start to get field data, (u, v, w, p, T, etc...)"
for i in np.arange(6):
	for j in np.arange(E):
		param[j][i]=data[start+j]
	start+=E

################Get cell points (n1, n2, n3, ... n8)##################
print "start to get cell points (n1, n2, n3,... n8)"
for i in np.arange(E):
	if(np.size(data[start+i].split( ))!=8):
		print("error0")
	arr=data[start+i].split( )
	for j in range(len(arr)):
		cell[i][j]=int(arr[j])-1


################Get cell points (n1, n2, n3, ... n8)##################
readFile="constant/polyMesh/points"
print ("start to get OF file: " + str(readFile))
f=open(readFile,"r")
data=f.readlines()
f.close()

start=21

MAXDIFF=0
for i in np.arange(N):
	arr=data[start+i].split( )
	rOF[i][0]=arr[0].split("(")[1]
	rOF[i][1]=arr[1].split("(")[0]
	rOF[i][2]=arr[2].split(")")[0]
	if(np.linalg.norm(r[i]-rOF[i])>MAXDIFF):
		MAXDIFF=np.linalg.norm(r[i]-rOF[i])


################Get cell points (n1, n2, n3, ... n8)##################
readFile="constant/polyMesh/owner"
print ("start to get OF file: " + str(readFile))
f=open(readFile,"r")
data=f.readlines()
f.close()

start=22

for i in np.arange(F):
	ownerOF[i]=data[start+i]

################Get cell points (n1, n2, n3, ... n8)##################
readFile="constant/polyMesh/neighbour"
print ("start to get OF file: " + str(readFile))
f=open(readFile,"r")
data=f.readlines()
f.close()

for i in np.arange(NEI):
	neighbourOF[i]=data[start+i]

################Get cell points (n1, n2, n3, ... n8)##################
readFile="constant/polyMesh/faces"
print ("start to get OF file: " + str(readFile))
f=open(readFile,"r")
data=f.readlines()
f.close()

start=21
cellOF=range(E)
for i in np.arange(E):
	cellOF[i]=[]

for i in np.arange(F):
	arr=data[start+i].split( )
	npoint=int(arr[0].split("(")[0])
	arr[0]=int(arr[0].split("(")[1])
	arr[npoint-1]=int(arr[npoint-1].split(")")[0])
	for j in arr:
		cellOF[ownerOF[i]].append(int(j))

for i in np.arange(NEI):
	arr=data[start+i].split( )
	npoint=int(arr[0].split("(")[0])
	arr[0]=int(arr[0].split("(")[1])
	arr[npoint-1]=int(arr[npoint-1].split(")")[0])
	for j in arr:
		cellOF[neighbourOF[i]].append(int(j))

for i in np.arange(E):
	cellOF[i]=np.unique(cellOF[i])
	if set(cellOF[i])!=set(cell[i]):
		print "error0"

################Get cell points (n1, n2, n3, ... n8)##################
print ("WRITING PRESSURE")
f=open("20000/p","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\nobject\tp;}\n\n\ndimensions [1 -1 -2 0 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
for i in param.T[5]:
	f.write(str(i)+"\n")

f.write(")\n;\n\n")
f.close()

print ("WRITING DENSITY")
f=open("20000/rho","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\nobject\trho;}\n\n\ndimensions [1 -3 0 0 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
for i in param.T[0]:
	f.write(str(i)+"\n")

f.write(")\n;\n\n")
f.close()

print ("WRITING TEMPERATURE")
f=open("20000/T","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\nobject\tT;}\n\n\ndimensions [0 0 0 1 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
for i in param.T[4]:
	f.write(str(i)+"\n")

f.write(")\n;\n\n")
f.close()

print ("WRITING VELOCITY")
f=open("20000/U","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolVectorField;\nlocation\t\"20000\";\nobject\tU;}\n\n\ndimensions [0 0 0 1 0 0 0];\n\ninternalField nonuniform List<vector>\n"+str(E)+"\n(\n")
for i in np.arange(E):
	f.write("("+str(param.T[1][i])+" "+str(param.T[2][i])+" "+str(param.T[3][i])+")\n")

f.write(")\n;\n\n")
f.close()

print ("WRITING OMEGA")
f=open("20000/omega","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\nobject\tomega;}\n\n\ndimensions [0 2 -3 0 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
for i in np.arange(E):
	f.write("1\n")

f.write(")\n;\n\n")
f.close()


print ("WRITING K")
f=open("20000/k","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\nobject\tk;}\n\n\ndimensions [0 2 -2 0 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
for i in np.arange(E):
	f.write("1\n")

f.write(")\n;\n\n")
f.close()



