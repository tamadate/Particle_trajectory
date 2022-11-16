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
fig, axs = plt.subplots(1,1,figsize=(7,5))
axNormal(axs)

direction=0
saveFigure=1
bid=7
idp=-5
vs=[]
dps=[]

directory="/home/tama3rdgen/MURI/NewNozzle/results/"
directory="/home/tama3rdgen/MURI/OldNozzleHe/results/"
f=open(directory+"../depositionVelocity.dat","w")
for idp in np.arange(-5.1,-7.1,-0.1):
    dp=10**idp*1e9
    data=np.loadtxt(directory+str(int(dp))+"/finalVelocity.dat")
    #data=np.loadtxt(dir+"/penetrateVelocity_0.dat")

    N=1e-10
    v=0.0
    for i in data:
    	if(i[4]==bid):
        #if(i[0]!=0):
            v+=i[direction]
            N+=1.0
    vs.append(v/N)
    dps.append(dp*1e-3)
    f.write(str(dp*1e-3)+" "+str(v/N)+"\n")

f.close()
axs.set_xlabel(r'$D_p$ [$\mu$m]',size=12)
axs.set_ylabel(r'$V_x$ [m/s]',size=12)
axs.set_xscale("log")
axs.set_xlim([0.09,10])
axs.scatter(dps,vs)

if(saveFigure):
    plt.savefig(directory+"../velocity.png", dpi=1000)
plt.show()
