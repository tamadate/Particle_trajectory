import matplotlib.pyplot as plt
import numpy as np

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

center=np.array((0,0,0))    # center of calculation domain
directory="/home/tama3rdgen/jetaxis/NewVersion/L18/25Torr/0.1um/"
directory="/home/tama3rdgen/jetaxis/NewMesh/100mm/10um/"
data=np.loadtxt(directory+"finalPosition.dat")

mappingInt=100
rout=0.04
iX=1
iY=2
saveFigure=1
bid=5

delta=rout*2/float(mappingInt)
x=np.arange(0.5*delta,mappingInt*delta,delta)-rout
dist=np.zeros((mappingInt,mappingInt))
Iout=np.sum(data.T[4]==bid)
Itot=np.size(data.T[4])

P=Iout/float(Itot)

for i in data:
	if(i[4]==bid):
		ix=int((i[iX]-center[iX]+rout)/delta)
		iy=int((i[iY]-center[iY]+rout)/delta)
		dist[iy][ix]+=1

rout*=1e3
X,Y=np.meshgrid(x,x)
axs.set_xlabel(r'$y$ [mm]',size=12)
axs.set_ylabel(r'$z$ [mm]',size=12)
axs.set_xlim(-rout,rout)
axs.set_ylim(-rout,rout)
DIST=np.max(dist)
if(DIST==0):
    DIST=1
pl=axs.contourf(X*1e3,Y*1e3,dist/DIST,levels=np.linspace(0,1,100),cmap="gray")
axs.set_title("Penetration: {:.01f} %".format(P*100),size=15)
plt.colorbar(pl,ticks=np.array([0,0.2,0.4,0.6,0.8,1.0])).set_label(label="Normalized particle concentration [-]",size=12)

if(saveFigure):
    plt.savefig(directory+"location.png", dpi=1000)
plt.show()
