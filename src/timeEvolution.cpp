#include "trajectory.hpp"


void
trajectory::euler(particle &a){
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


void
trajectory::analytical(particle &a){
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
