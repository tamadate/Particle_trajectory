#pragma onece
#include "functions.hpp"


double computeFD(Variables *vars, particle &par) {
	double dUx=vars->U.x[0]-par.v.x[0];
	double dUy=vars->U.x[1]-par.v.x[1];
	double dUz=vars->U.x[2]-par.v.x[2];
	double U2=dUx*dUx+dUy*dUy+dUz*dUz;
	double Umag=sqrt(U2);
	par.Re=rho*Umag*par.dp/mu0+1e-20;
	double Cd=24.0/(par.Re*par.Cc);

	return M_PI*0.125*mu0*par.dp*par.Re*Cd/par.m;
}

	

	

