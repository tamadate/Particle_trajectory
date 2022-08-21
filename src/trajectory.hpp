#pragma once
#include "dragForce.hpp"
#include "output.hpp"


void euler(Variables *vars, Flags *flags, particle &a){
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		dt=1e-6/vmag;
		//dt=a.dt;
	}
	for(auto &force : forces) force->computeFD(vars,flags,a);
	for(int i=0; i<3; i++) a.v.x[i]+=dt*a.F.x[i];
	for(int i=0; i<3; i++) a.x.x[i]+=dt*a.v.x[i];

    if(a.x.x[plane2D]<Axis && flags->dimensionFlag>0) {
        a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
        a.v.x[plane2D]*=-1.0;
        a.F.x[plane2D]*=-1.0;
        a.reflect*=-1;
	}
	vars->time+=dt;

}


void analytical(Variables *vars, Flags *flags, particle &a){
	if(flags->autoStep) dt=a.dt;
	double Ux=vars->U[a.cell].x[0];
	double Uy=vars->U[a.cell].x[1];
	double Uz=vars->U[a.cell].x[2];
	double dvx=a.v.x[0]-Ux;
	double dvy=a.v.x[1]-Uy;
	double dvz=a.v.x[2]-Uz;
	double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
	double vmag=sqrt(v2);

	double dt_anal=1e-6/vmag;


	double EXP=exp(-dt_anal*a.beta);
	double dvxEXP=dvx*EXP;
	double dvyEXP=dvy*EXP;
	double dvzEXP=dvz*EXP;

	a.v.x[0]=Ux+dvxEXP;
	a.v.x[1]=Uy+dvyEXP;
	a.v.x[2]=Uz+dvzEXP;

	a.x.x[0]+=Ux*dt_anal-dvxEXP/a.beta;
	a.x.x[1]+=Uy*dt_anal-dvyEXP/a.beta;
	a.x.x[2]+=Uz*dt_anal-dvzEXP/a.beta;

    if(a.x.x[plane2D]<Axis && flags->dimensionFlag>0) {
        a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
        a.v.x[plane2D]*=-1.0;
        a.F.x[plane2D]*=-1.0;
        a.reflect*=-1;
	}
	vars->time+=dt_anal;

}

void trajectory(Variables *vars, Flags *flags){
	timer=clock();
	const int particleSize=vars->particles.size();
	int trapParticle=0;	// Number of cases when the calculaiton is not finished in required time steps
//	# pragma omp parallel for
	for (int pid=0; pid<particleSize; pid++){
		particle &a=vars->particles[pid];
		initialize(vars,flags);
		outputInitial(a);

		while(vars->time<totalTime){
			if(vars->time - vars->preOutTime > observeTime) {
				output(a, vars->time);
				vars->preOutTime=vars->time;
			}
			computeReMach(vars,a);
			if(a.Re<0.01 && a.Mach<0.1 && flags->analytical==1){analytical(vars,flags,a);}
			else{euler(vars,flags,a);}

			flags->breakFlag=checkCell(vars,pid);
			if(flags->breakFlag==-1) break;
		}
		if(flags->breakFlag==0){
			outParticle op;
			op.pid=a.id;
			op.r=a.x;
			op.v=a.v;
			op.bid=-1;
			outParticles.push_back(op);
			trapParticle++;  
		}
		outputFinalPosition(outParticles[pid]);
	}
    cout<<trapParticle<<" might be trapped circulation"<<endl;    
	cout<<"Trajectory calculation time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
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


