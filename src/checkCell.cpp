#include "trajectory.hpp"

int
trajectory::checkCell(int pid){
	int escapeFlag=0;
	int loopFlag;
	while(escapeFlag<100){
		loopFlag=0;	//0:continue, -1:out from bound, 1:move to next cell
		int icell=vars->particles[pid].cell;
		int faceSize=cells[icell].iface.size();
		for (int i=0; i<faceSize; i++){
			int b=cells[icell].iface[i];
			int c0=faces[b].iface[0];
			double x0=points[c0].x[0];
			double y0=points[c0].x[1];
			double z0=points[c0].x[2];
			double x=vars->particles[pid].x.x[0]-x0;
			double y=vars->particles[pid].x.x[1]-y0;
			double z=vars->particles[pid].x.x[2]-z0;
			point c=cells[icell].norm[i];
			double inProduct=c.x[0]*x+c.x[1]*y+c.x[2]*z;
		    if(b>boundaryStartID-1)	{
				if(inProduct<vars->particles[pid].dp*0.5){
			        loopFlag=boundAction(b,pid,c);
				    if(loopFlag<0) break;
                }
		    }
			else{
				if(inProduct<0){
					loopFlag=1;
					escapeFlag++;
					if(owners[b]==vars->particles[pid].cell) vars->particles[pid].cell=neighbors[b];
					else vars->particles[pid].cell=owners[b];
					int newID=vars->particles[pid].cell;
					double Kn=vars->lamda[newID]/vars->particles[pid].dp;
					vars->particles[pid].Kn=Kn;
					vars->particles[pid].Cc=1+Kn*(A1+A2*exp(-A3/Kn));
					double threePiMuDp_Cc=3*M_PI*vars->myu[icell]*vars->particles[pid].dp/vars->particles[pid].Cc;
					vars->particles[pid].Zp=1.6e-19/threePiMuDp_Cc;
					vars->particles[pid].beta=threePiMuDp_Cc/vars->particles[pid].m;
					double tau100=1/(100*vars->particles[pid].beta);

					vars->particles[pid].dt=1/(100*vars->particles[pid].beta);

					break;
				}
			}
		}
		if(loopFlag<1) break;
	}
	if(loopFlag==1) loopFlag=0;
	return loopFlag;
}


int
trajectory::boundAction(int faceID, int pid, point norm){
	int boundSize=boundaries.size();
	int returnInt=0;
	for(int i=0; i<boundSize; i++){
		if(faceID+1>boundaries[i].startFace&&faceID<boundaries[i].startFace+boundaries[i].nFaces) {
			if(boundaries[i].type=="patch"){
				outParticles[pid].pid=vars->particles[pid].id;
				outParticles[pid].r=vars->particles[pid].x;
				outParticles[pid].v=vars->particles[pid].v;
				outParticles[pid].bid=i;
				//cout<<"out from "<<boundaries[i].name<<" "<<faceID<<endl;
				returnInt=-1;
			}
			if(boundaries[i].type=="symmetryPlane"){
				double b0=vars->particles[pid].x.x[0]-points[faces[faceID].iface[0]].x[0];
				double b1=vars->particles[pid].x.x[1]-points[faces[faceID].iface[0]].x[1];
				double b2=vars->particles[pid].x.x[2]-points[faces[faceID].iface[0]].x[2];
				double nb0=norm.x[0]*b0+norm.x[1]*b1+norm.x[2]*b2;
				vars->particles[pid].x.x[0]-=2*nb0*norm.x[0];
				vars->particles[pid].x.x[1]-=2*nb0*norm.x[1];
				vars->particles[pid].x.x[2]-=2*nb0*norm.x[2];

				double v0=vars->particles[pid].v.x[0];
				double v1=vars->particles[pid].v.x[1];
				double v2=vars->particles[pid].v.x[2];
				double nv0=norm.x[0]*v0+norm.x[1]*v1+norm.x[2]*v2;

				vars->particles[pid].v.x[0]-=2*nv0*norm.x[0];
				vars->particles[pid].v.x[1]-=2*nv0*norm.x[1];
				vars->particles[pid].v.x[2]-=2*nv0*norm.x[2];


				v0=vars->particles[pid].F.x[0];
				v1=vars->particles[pid].F.x[1];
				v2=vars->particles[pid].F.x[2];
				nv0=norm.x[0]*v0+norm.x[1]*v1+norm.x[2]*v2;

				vars->particles[pid].F.x[0]-=2*nv0*norm.x[0];
				vars->particles[pid].F.x[1]-=2*nv0*norm.x[1];
				vars->particles[pid].F.x[2]-=2*nv0*norm.x[2];
				returnInt=1;
			}
			if(boundaries[i].type=="wall"){
				outParticles[pid].pid=vars->particles[pid].id;
				outParticles[pid].r=vars->particles[pid].x;
				outParticles[pid].v=vars->particles[pid].v;
				outParticles[pid].bid=i;
				//cout<<"hit to "<<boundaries[i].name<<endl;
				returnInt=-1;
			}
		}
	}
	return returnInt;

}