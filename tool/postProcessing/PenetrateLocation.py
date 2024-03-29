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



center=np.array((0,0,0))    # center of calculation domain
directory="/home/tama3rdgen/jetaxis/NewVersion/L18/25Torr/0.1um/"
directory="/home/tama3rdgen/jetaxis/NewMesh/100mm/5um/"
for penID in np.arange(9):
    pltNormal()
    fig, axs = plt.subplots(1,1,figsize=(7,5))
    axNormal(axs)
    data=np.loadtxt(directory+"penetratePosition_"+str(penID)+".dat")

    mappingInt=100
    rout=0.04
    iX=1
    iY=2
    saveFigure=1

    delta=rout*2/float(mappingInt)
    x=np.arange(0.5*delta,mappingInt*delta,delta)-rout
    dist=np.zeros((mappingInt,mappingInt))
    Iout=0
    Itot=np.size(data.T[0])
    for i in data:
        if(i[0]==0 and i[1]==0 and i[2]==0):
            Iout+=1
        else:
            ix=int((i[iX]-center[iX]+rout)/delta)
            iy=int((i[iY]-center[iY]+rout)/delta)
            dist[iy][ix]+=1
            distv[iy][ix]+=datav[index][0]

    P=Iout/float(Itot)
    rout*=1e3
    X,Y=np.meshgrid(x,x)
    axs.set_xlabel(r'$y$ [mm]',size=12)
    axs.set_ylabel(r'$z$ [mm]',size=12)
    axs.set_xlim(-rout,rout)
    axs.set_ylim(-rout,rout)
    DIST=np.max(distv)
    if(DIST==0):
        DIST=1
    pl=axs.contourf(X*1e3,Y*1e3,distv/DIST,levels=np.linspace(0,500,100),cmap="gray")
    axs.set_title("Penetration: {:.01f} %".format(P*100),size=15)
    #plt.colorbar(pl,ticks=np.array([0,0.2,0.4,0.6,0.8,1.0])).set_label(label="Normalized particle concentration [-]",size=12)
    plt.colorbar(pl,ticks=np.array([0,100,200,300,400,500])).set_label(label="Normalized particle concentration [-]",size=12)

    if(saveFigure):
        plt.savefig(directory+"location"+str(penID)+".png", dpi=1000)
    plt.show(block=False)
    plt.close()
