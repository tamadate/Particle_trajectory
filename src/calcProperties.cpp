#include "trajectory.hpp"

void
trajectory::calculateMyu(void){
	string str;
	int iflag=0;

	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double myu=myu0*pow(vars->T[i]/ST0,1.5)*(ST0+SC0)/(vars->T[i]+SC0);
		vars->myu.push_back(myu);
	}
}

void
trajectory::calculatelamda(void){
	string str;
	int iflag=0;

	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double lamda=lamda_coeff*vars->myu[i]/vars->rho[i]/sqrt(vars->T[i]);
		vars->lamda.push_back(lamda);
	}
}


void
trajectory::computeReMach(particle &par){
	double dUx=vars->U[par.cell].x[0]-par.v.x[0];
	double dUy=vars->U[par.cell].x[1]-par.v.x[1];
	double dUz=vars->U[par.cell].x[2]-par.v.x[2];
	double U2=dUx*dUx+dUy*dUy+dUz*dUz;
	double Umag=sqrt(U2);
	double cg=sqrt(gamkb_m*vars->T[par.cell]);

	par.Re=vars->rho[par.cell]*Umag*par.dp/vars->myu[par.cell]+1e-20;
	par.Mach=Umag/cg;
}
