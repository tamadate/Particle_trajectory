#include "trajectory.hpp"


double
trajectory::get_tini(particle &par){
    random_device seed;
	  mt19937 mt(seed());
    uniform_real_distribution<> randZeroOne(0,1);
    double dvx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
    double dvy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
    double dvz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
    double dv=sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
    double vrand=sqrt(par.Urand.x[0]*par.Urand.x[0]+par.Urand.x[1]*par.Urand.x[1]+par.Urand.x[2]*par.Urand.x[2]);	// scale of random velocity
    double tau=rho_p*par.dp*par.dp/(18.0*vars->myu[par.cell]);	// This equation is not exact at high Reynolds number
    double le=0.3*pow(vars->k[par.cell],1.5)/vars->epsilon[par.cell];	// le = Cm^0.5 * k^1.5 / epsilon, where Cm = 0.09 is used
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
	normal_distribution<> dist_vTD(0.0, sqrt(vars->k[par.cell]*2));	// standard diviation, sigma = (2k)^0.5
  par.Urand.x[0]=dist_vTD(mt);
  par.Urand.x[1]=dist_vTD(mt);
  par.Urand.x[2]=dist_vTD(mt);
  par.tini=get_tini(par);
}
