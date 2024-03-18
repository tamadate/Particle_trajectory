#pragma once
#include "../force.hpp"


//Stokes-Millikan
class dragForceSM : public force{
	public:
		void compute(particle &par){
			double Cd=24/par.Re/par.Cc;;	// Here, this code use another equations for each drag models

			double FD=M_PI*0.125*vars->myu[par.cell]*par.dp*par.Re*Cd/par.m;

			double dUx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
			double dUy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
			double dUz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
			par.F.x[0]+=FD*dUx;
			par.F.x[1]+=FD*dUy;
			par.F.x[2]+=FD*dUz;
		};
		dragForceSM(void){};
		~dragForceSM(void){};
	private:
};