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
nozzle=0.1787-6e-3

faces=np.arange(0.173,0.1787,0.0002)
faces=np.append(faces,0.1787)
idps=np.arange(-5.1,-7.1,-0.2)
vs=np.zeros((np.size(idps),np.size(faces)))

#directory="/home/tama3rdgen/MURI/NewNozzle/results/"
directory="/home/tama3rdgen/MURI/NewNozzle/results/"


ii=0
tii=np.size(idps)
for idp in idps:
    dp=10**idp*1e9
    for fid in np.arange(np.size(faces)):
        dataV=np.loadtxt(directory+str(int(dp))+"/penetrateVelocity_"+str(fid)+".dat")
        dataR=np.loadtxt(directory+str(int(dp))+"/penetratePosition_"+str(fid)+".dat")
        N=1e-10
        v=0.0
        for i in np.arange(np.size(dataV.T[0])):
            if(dataV[i][direction]!=0 and rangeStart-range<dataR[i][rangeDirection]<rangeStart+range):
                v+=dataV[i][direction]
                N+=1.0
        vs[ii][fid]=v/N
    axs.plot(1000*(faces-nozzle),vs[ii],label=r"$D$_p = {:.2f} $\mu$m".format(dp*1e-3),c=(0.8*float(ii/tii)+0.1,0,0))
    ii+=1
    print(dp)
axs.set_xlabel(r'Distance from nozzle, $z$ [mm]',size=12)
axs.set_ylabel(r'$V_z$ [m/s]',size=12)
#axs.set_xscale("log")
axs.set_xlim([0,6.2])
#plt.legend()
if(saveFigure):
    plt.savefig(directory+"../distanceFromNozzle.png", dpi=1000)
plt.show()
