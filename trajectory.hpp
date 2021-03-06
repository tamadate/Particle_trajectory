#pragma once
#include "functions.hpp"
#include "output.hpp"


void trajectory(Variables *vars){
	for (auto &a:vars->particles) updateDisp(vars, a);
	const int particleSize=vars->particles.size();
	int trapParticle=0;
//	# pragma omp parallel for
	for (int pid=0; pid<particleSize; pid++){
		particle &a=vars->particles[pid];
		int itime=0;
		int loop=0;
		int Floop=timestep/Observe;
		int breakFlag=0;
		int penetrationFlag=0;
		while(loop<Floop){
			for(int i=0; i<3; i++) a.v.x[i]+=0.5*dt*a.F.x[i];
			for(int i=0; i<3; i++) a.x.x[i]+=dt*a.v.x[i];
		    if(vars->particles[pid].x.x[plane2D]<Axis && dimensionFlag>0) {
		        a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
		        a.v.x[plane2D]*=-1.0;
		        a.F.x[plane2D]*=-1.0;
		        a.reflect*=-1;
			}
			breakFlag=checkCell(vars,pid);
			double FD=computeFD(vars, a);
			for(int i=0; i<3; i++) a.F.x[i]=FD*(vars->U[a.cell].x[i]+a.Urand.x[i]-a.v.x[i])-vars->dp[a.cell].x[i]/rho_p;

			for(int i=0; i<3; i++) a.v.x[i]+=0.5*dt*a.F.x[i];

			if(dispFlag==1){
				a.randLoop++;
				if(a.randLoop==a.randUpdate) updateDisp(vars, a);
			}
			if(a.x.x[0]>0.173&&penetrationFlag==0){
				outParticle op;
				op.pid=a.id;
				op.r=a.x;
				op.v=a.v;
				op.bid=0;
				outputPenetration(op);
				penetrationFlag=1;
			}
			if(itime%Observe==0){
				vars->time=dt*loop*Observe;
		        loop++;
		        itime=0;
				output(a, vars->time);
				f=fopen("out.dat", "a");
				fprintf(f,"%f\t%f\t%f\t%d\n", a.x.x[0], a.x.x[1], a.x.x[2], loop);
				fclose(f);
/*				int icell=a.cell;
				int faceSize=cells[icell].iface.size();
				for (int i=0; i<faceSize; i++){
					int b=cells[icell].iface[i];
					int c0=faces[b].iface[0];
					double x0=points[c0].x[0];
					double y0=points[c0].x[1];
					double z0=points[c0].x[2];
					double x=a.x.x[0]-x0;
					double y=a.x.x[1]-y0;
					double z=a.x.x[2]-z0;
					point c=cells[icell].norm[i];
					double inProduct=c.x[0]*x+c.x[1]*y+c.x[2]*z;
					f=fopen("inpro.dat", "a");
					fprintf(f,"%f\t%f\t%f\t%e\n", c.x[0], c.x[1], c.x[2], inProduct);
					fclose(f);
					f=fopen("cell.dat", "a");
					for (auto &c0:faces[b].iface){
						fprintf(f,"%f\t%f\t%f\t%d\n", points[c0].x[0], points[c0].x[1], points[c0].x[2], b);
					}
					fprintf(f,"\n");
					fclose(f);
				}*/
			}
			if(breakFlag==-1) break;
			itime++;
		}
		if(breakFlag==0){
			outParticle op;
			op.pid=a.id;
			op.r=a.x;
			op.v=a.v;
			op.bid=-1;
			outParticles.push_back(op);
			trapParticle++;
			cout<<op.pid<<" might be trapped circulation"<<endl;    
		}
		outputFinalPosition(outParticles[pid]);
	}
    cout<<trapParticle<<" might be trapped circulation"<<endl;    
}


int checkCell(Variables *vars, int pid){
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
			        loopFlag=boundAction(vars,b,pid,c);
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
		if(loopFlag<1) break;
	}
	if(loopFlag==1) loopFlag=0;
	return loopFlag;
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
				op.v=vars->particles[particleID].v;
				op.bid=i;
				outParticles.push_back(op);
				//cout<<"out from "<<boundaries[i].name<<" "<<faceID<<endl;
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
				op.v=vars->particles[particleID].v;
				op.bid=i;
				outParticles.push_back(op);
				//cout<<"hit to "<<boundaries[i].name<<endl;
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

