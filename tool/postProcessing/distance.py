import matplotlib.pyplot as plt
import numpy as np
import glob
import os
##  To organize figure style
def pltNormal():
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['figure.subplot.bottom'] = 0.15
    plt.rcParams['figure.subplot.left'] = 0.2
    plt.rcParams["font.size"]=10

def axNormal(ax):
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')


pltNormal()
fig, axs = plt.subplots(1,1,figsize=(5.5,5))
axNormal(axs)

direction=0
rangeStart=200e-6
range=100e-6
rangeDirection=1
saveFigure=1
bid=7

#directory="/home/tama3rdgen/MURI/NewNozzle/results/"
directory="/home/tama3rdgen/MURI/OldNozzleHe/results/"
f=open(directory+"Velocities.dat","w")
for fid in np.arange(1):
    vs=[]
    dps=[]
    for idp in np.arange(-5.2,-6.2,-0.1):
        dp=10**idp*1e9
        dataV=np.loadtxt(directory+str(int(dp))+"/penetrateVelocity_"+str(fid)+".dat")
        dataR=np.loadtxt(directory+str(int(dp))+"/penetratePosition_"+str(fid)+".dat")
        N=1e-10
        v=0.0
        for i in np.arange(np.size(dataV.T[0])):
            if(dataV[i][direction]!=0):
                v+=dataV[i][direction]
                N+=1.0
        vs.append(v/N)
        dps.append(dp*1e-3)
    axs.scatter(dps,vs,c="blue")
    f.write(str(dps[0])+" "+str(vs[0])+"\n")

f.close()


directory="/home/tama3rdgen/MURI/OldNozzle/results/"
f=open(directory+"Velocities.dat","w")
for fid in np.arange(1):
    vs=[]
    dps=[]
    for idp in np.arange(-5.1,-5.6,-0.05):
        dp=10**idp*1e9
        dataV=np.loadtxt(directory+str(int(dp))+"/penetrateVelocity_"+str(fid)+".dat")
        dataR=np.loadtxt(directory+str(int(dp))+"/penetratePosition_"+str(fid)+".dat")
        N=1e-10
        v=0.0
        for i in np.arange(np.size(dataV.T[0])):
            if(dataV[i][direction]!=0):
                v+=dataV[i][direction]
                N+=1.0
        vs.append(v/N)
        dps.append(dp*1e-3)
    axs.scatter(dps,vs,c="red")
    f.write(str(dps[0])+" "+str(vs[0])+"\n")

f.close()
axs.set_xlabel(r'Particle diameter, $D_p$ [$\mu$m]',size=12)
axs.set_ylabel(r'$V_z$ [m/s]',size=12)
axs.set_xscale("log")
axs.set_xlim([0.5,10])
axs.set_ylim([200,1500])
if(saveFigure):
    plt.savefig(directory+"../penVelocity.png", dpi=1000)
plt.show()
