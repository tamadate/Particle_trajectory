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
fig, axs = plt.subplots(1,1,figsize=(5,10))
axNormal(axs)

center=np.array((0,0,0))    # center of calculation domain
directory="/home/tama3rdgen/jetaxis/NewMesh/100nm/"
Nparticle=int(np.loadtxt(directory+"result/nparticle"))

skipLevel=2
X=0
Y=1
saveFigure=0

for I in np.arange(int(Nparticle/skipLevel)):
    i=I*skipLevel
    data=np.loadtxt(directory+"result/position."+str(i))
    axs.set_xlabel(r'$x$ [mm]',size=12)
    axs.set_ylabel(r'$y$ [mm]',size=12)
    axs.plot(data.T[X]-center[X],data.T[Y]-center[Y],linewidth=0.5)

if(saveFigure):
    plt.savefig(directory+"trajectory.png", dpi=1000)
plt.show()
