#include "trajectory.hpp"

void
trajectory::calculateMyu(void){
	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double myu=myu0*pow(vars->T[i]/ST0,1.5)*(ST0+SC0)/(vars->T[i]+SC0);	// calculate viscosity via Sutherland's equation
		vars->myu.push_back(myu);
	}
}


void
trajectory::calculatelamda(void){
	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double lamda=lamda_coeff*vars->myu[i]/vars->rho[i]/sqrt(vars->T[i]); // calculate mean free path
		vars->lamda.push_back(lamda);
	}
}


void
trajectory::calculateNonDimension(particle &par){
	// calculate particle relative velocity
	double dUx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
	double dUy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
	double dUz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
	double U2=dUx*dUx+dUy*dUy+dUz*dUz;
	double Umag=sqrt(U2);

	double cg=sqrt(gamkb_m*vars->T[par.cell]);	// speed of sound

	par.Re=vars->rho[par.cell]*Umag*par.dp/vars->myu[par.cell]+1e-100;
	par.Mach=Umag/cg;
	par.Kn=par.Mach/par.Re*Kn_coeff;
	par.Cc=1+par.Kn*(A1+A2*exp(-A3/par.Kn));
	par.fric=3*M_PI*vars->myu[par.cell]*par.dp/par.Cc;
	par.beta=par.fric/par.m;
}
