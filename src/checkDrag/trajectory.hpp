#pragma onece
#include "functions.hpp"
#include "output.hpp"


void trajectory(Variables *vars){
	int itime=0;
    int loop=0;
    int Floop=timestep/Observe;
	for (auto &a:vars->particles) updateDisp(vars, a);
	while(loop<Floop){
		velocityUpdate(vars);
		positionUpdate(vars);
        if(dimensionFlag>0) symmetricModification(vars);
		checkCell(vars);
		if(vars->particles.size()==0) break;
		forceUpdate(vars);
		velocityUpdate(vars);

		if(dispFlag==1){
			for (auto &a:vars->particles) {
				a.randLoop++;
				if(a.randLoop==a.randUpdate) updateDisp(vars, a);
			}
		}
		if(itime%Observe==0){
			vars->time=dt*loop*Observe;
            loop++;
            itime=0;
			output(vars);
			for(auto &a: vars->particles){
				f=fopen("out.dat", "a");
				fprintf(f,"%f\t%f\t%f\t%d\n", a.x.x[0], a.x.x[1], a.x.x[2], loop);
				fclose(f);
/*				int icell=a.cell;
				int faceSize=cells[icell].iface.size();
				for (int i=0; i<faceSize; i++){
					int b=cells[icell].iface[i];
					for (auto &c0:faces[b].iface){
						f=fopen("cell.dat", "a");
						fprintf(f,"%f\t%f\t%f\t%d\n", points[c0].x[0], points[c0].x[1], points[c0].x[2], itime);
						fclose(f);
					}
				}*/
			}
		}
		itime++;
	}
    cout<<"Final time: "<<dt*itime<<" s"<<endl;
    cout<<vars->particles.size()<<" might be trapped circulation"<<endl;
    for(auto &a: vars->particles){
		outParticle op;
		op.pid=a.id;
		op.r=a.x;
		op.bid=-1;
		outParticles.push_back(op);
    }
    
}


void velocityUpdate(Variables *vars){
	for (auto &a: vars->particles){
		for(int i=0; i<3; i++){
			a.v.x[i]+=0.5*dt*a.F.x[i];
		}
	}
}

void positionUpdate(Variables *vars){
	for (auto &a: vars->particles){
		for(int i=0; i<3; i++){
			a.x.x[i]+=dt*a.v.x[i];
		}
	}
}


void forceUpdate(Variables *vars){
	for (auto &a: vars->particles){
		double FD=computeFD(vars, a);
		int icell=a.cell;
		for(int i=0; i<3; i++){
			a.F.x[i]=FD*(vars->U[icell].x[i]+a.Urand.x[i]-a.v.x[i])-vars->dp[icell].x[i]/rho_p;
		}
	}
}


void symmetricModification(Variables *vars){
	int particleSize=vars->particles.size();
	for(int pid=0; pid<particleSize; pid++){
		if(vars->particles[pid].x.x[plane2D]<Axis){
            vars->particles[pid].x.x[plane2D]=2*Axis-vars->particles[pid].x.x[plane2D];
            vars->particles[pid].v.x[plane2D]*=-1.0;
            vars->particles[pid].F.x[plane2D]*=-1.0;
            vars->particles[pid].reflect*=-1;
        }
    }
}

void checkCell(Variables *vars){
	int particleSize=vars->particles.size();
	for(int pid=0; pid<particleSize; pid++){
		int escapeFlag=0;
		while(escapeFlag<100){
			int loopFlag=0;
			
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
				        loopFlag=boundAction(vars,b,pid,c);
    				    break;
                    }
			    }
				else{
					if(inProduct<0){
						loopFlag=1;
						escapeFlag++;
						if(owners[b]==vars->particles[pid].cell) vars->particles[pid].cell=neighbors[b];
						else vars->particles[pid].cell=owners[b];
						int newID=vars->particles[pid].cell;
						double Kn=vars->ramda[newID]/vars->particles[pid].dp;
						double Cc=1+Kn*(A1+A2*exp(-A3/Kn));
						vars->particles[pid].Kn=Kn;
						vars->particles[pid].Cc=Cc;
						vars->particles[pid].Zp=Cc*1.6e-19/(3*M_PI*vars->myu[newID]*vars->particles[pid].dp);
						vars->particles[pid].beta=3.0*M_PI*vars->myu[newID]*vars->particles[pid].dp/Cc/vars->particles[pid].m;
						vars->particles[pid].le=pow(vars->k[newID],0.5)/vars->omega[newID]/pow(0.09,0.25);
						
						break;
					}
				}
			}
			if(loopFlag<1) {pid+=loopFlag; particleSize+=loopFlag; break;}
		}
	}

}


int boundAction(Variables *vars, int faceID, int particleID, point norm){
	int boundSize=boundaries.size();
	int returnInt=0;
	for(int i=0; i<boundSize; i++){
		if(faceID+1>boundaries[i].startFace&&faceID<boundaries[i].startFace+boundaries[i].nFaces) {
			if(boundaries[i].type=="patch"){
				outParticle op;
				op.pid=vars->particles[particleID].id;
				op.r=vars->particles[particleID].x;
				op.bid=i;
				outParticles.push_back(op);
				vars->particles.erase(vars->particles.begin()+particleID);
				cout<<"out from "<<boundaries[i].name<<" "<<faceID<<endl;
				returnInt=-1;
			}
			if(boundaries[i].type=="symmetryPlane"){
				double b0=vars->particles[particleID].x.x[0]-points[faces[faceID].iface[0]].x[0];
				double b1=vars->particles[particleID].x.x[1]-points[faces[faceID].iface[0]].x[1];
				double b2=vars->particles[particleID].x.x[2]-points[faces[faceID].iface[0]].x[2];
				double nb0=norm.x[0]*b0+norm.x[1]*b1+norm.x[2]*b2;
				vars->particles[particleID].x.x[0]-=2*nb0*norm.x[0];
				vars->particles[particleID].x.x[1]-=2*nb0*norm.x[1];
				vars->particles[particleID].x.x[2]-=2*nb0*norm.x[2];

				double v0=vars->particles[particleID].v.x[0];
				double v1=vars->particles[particleID].v.x[1];
				double v2=vars->particles[particleID].v.x[2];
				double nv0=norm.x[0]*v0+norm.x[1]*v1+norm.x[2]*v2;

				vars->particles[particleID].v.x[0]-=2*nv0*norm.x[0];
				vars->particles[particleID].v.x[1]-=2*nv0*norm.x[1];
				vars->particles[particleID].v.x[2]-=2*nv0*norm.x[2];

				v0=vars->particles[particleID].Urand.x[0];
				v1=vars->particles[particleID].Urand.x[1];
				v2=vars->particles[particleID].Urand.x[2];
				nv0=norm.x[0]*v0+norm.x[1]*v1+norm.x[2]*v2;

				vars->particles[particleID].Urand.x[0]-=2*nv0*norm.x[0];
				vars->particles[particleID].Urand.x[1]-=2*nv0*norm.x[1];
				vars->particles[particleID].Urand.x[2]-=2*nv0*norm.x[2];

				v0=vars->particles[particleID].F.x[0];
				v1=vars->particles[particleID].F.x[1];
				v2=vars->particles[particleID].F.x[2];
				nv0=norm.x[0]*v0+norm.x[1]*v1+norm.x[2]*v2;

				vars->particles[particleID].F.x[0]-=2*nv0*norm.x[0];
				vars->particles[particleID].F.x[1]-=2*nv0*norm.x[1];
				vars->particles[particleID].F.x[2]-=2*nv0*norm.x[2];
				returnInt=1;
			}
			if(boundaries[i].type=="wall"){
				outParticle op;
				op.pid=vars->particles[particleID].id;
				op.r=vars->particles[particleID].x;
				op.bid=i;
				outParticles.push_back(op);
				vars->particles.erase(vars->particles.begin()+particleID);
				cout<<"hit to "<<boundaries[i].name<<endl;
				returnInt=-1;
			}
		}
	}
	return returnInt;
	
}


void updateDisp(Variables *vars, particle &par){
	if(dispFlag==1){
		for(int i=0;i<3;i++) par.Urand.x[i]=get_vTD(vars->k[par.cell]);
		double tini=get_tini(vars, par);
		par.randUpdate=int(tini/dt)+1;
		par.randLoop=0;
	}
	if(dispFlag==0){
		for(int i=0;i<3;i++) par.Urand.x[i]=0;
		par.randUpdate=0;
		par.randLoop=0;
	}
}


double get_vTD(double dispEnergy){
	random_device seed;
	mt19937 mt(seed());
	normal_distribution<> dist_vTD(0.0, sqrt(dispEnergy/1.5));
	return dist_vTD(mt);
}


double get_tini(Variables *vars, particle par){
    random_device seed;
	mt19937 mt(seed());
    uniform_real_distribution<> randZeroOne(0, 1);
    double dv0=(vars->U[par.cell].x[0]+par.Urand.x[0])-par.v.x[0];
    double dv1=(vars->U[par.cell].x[1]+par.Urand.x[1])-par.v.x[1];
    double dv2=(vars->U[par.cell].x[2]+par.Urand.x[2])-par.v.x[2];
    double dv=sqrt(dv0*dv0+dv1*dv1+dv2*dv2);
    double vrand=sqrt(par.Urand.x[0]*par.Urand.x[0]+par.Urand.x[1]*par.Urand.x[1]+par.Urand.x[2]*par.Urand.x[2]);
	double Cd=computeCd(vars, par);
    double tau=4*rho_p*par.dp*par.dp/(3.0*vars->myu[par.cell]*Cd*par.Re);
    double TL=0.15/vars->omega[par.cell]/0.09;
    double te=-TL*log(randZeroOne(mt));
    double Le=vrand*te;
    double CtR=Le/tau/dv;
    double tini;
    if(CtR>1) tini=te;
    else {
        double tcross=-tau*log(1-CtR);
        if(te<tcross) tini=te;
        else tini=tcross;
    }
    return tini;
}

