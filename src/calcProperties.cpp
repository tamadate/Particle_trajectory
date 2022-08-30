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


double
trajectory::get_tini(particle &par){
    random_device seed;
	  mt19937 mt(seed());
    uniform_real_distribution<> randZeroOne(0,1);
    double dvx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
    double dvy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
    double dvz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
    double dv=sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
    double vrand=sqrt(par.Urand.x[0]*par.Urand.x[0]+par.Urand.x[1]*par.Urand.x[1]+par.Urand.x[2]*par.Urand.x[2]);
    double tau=rho_p*par.dp*par.dp/(18.0*vars->myu[par.cell]);
    double le=0.866*pow(vars->k[par.cell],1.5)/vars->epsilon[par.cell];
		double te=le/vrand;
    double inLog=1-le/(tau*dv);
    double tini=te;
    if(inLog<0) tini=te;
    else {
        double tcross=-tau*log(inLog);
        if(tcross<tini) tini=tcross;
    }
    return tini;
}


void
trajectory::updateDisp(particle &par){
  random_device seed;
	mt19937 mt(seed());
	normal_distribution<> dist_vTD(0.0, sqrt(vars->k[par.cell]*0.5));
  par.Urand.x[0]=dist_vTD(mt);
  par.Urand.x[1]=dist_vTD(mt);
  par.Urand.x[2]=dist_vTD(mt);
  par.tini=get_tini(par);
}
