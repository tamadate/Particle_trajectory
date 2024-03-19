#pragma once
#include "TD.hpp"

// with turbulent dispersion
class TD_on : public TD{
	public:
		double rho_p;
		string name="turbulent dispersion ON";

		void tini(particle &par){
			double dvx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
			double dvy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
			double dvz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
			double dv=sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
			double vrand=sqrt(par.Urand.x[0]*par.Urand.x[0]+par.Urand.x[1]*par.Urand.x[1]+par.Urand.x[2]*par.Urand.x[2]);	// scale of random velocity
			double tau=rho_p*par.dp*par.dp/(18.0*vars->myu[par.cell]);	// This equation is not exact at high Reynolds number
			double le=0.16432*pow(vars->k[par.cell],1.5)/vars->epsilon[par.cell];	// le = Cm^0.5 * k^1.5 / epsilon, where Cm = 0.09 is used
			double te=le/vrand;
			double inLog=1-le/(tau*dv);
			double tini=te;
			if(inLog<0) tini=te;
			else {
				double tcross=-tau*log(inLog);
				if(tcross<tini) tini=tcross;
			}
			par.tini=tini;
		};

		void update(particle &par){
			// update velocity if the tini < 0
			if(par.tini<=0){
				random_device seed;
				mt19937 mt(seed());
				normal_distribution<> dist_vTD(0.0, sqrt(vars->k[par.cell]*2*0.1));	// standard diviation, k=1.5*v*v

				par.Urand.x[0]=dist_vTD(mt);
				par.Urand.x[1]=dist_vTD(mt);
				par.Urand.x[2]=dist_vTD(mt);
				tini(par);
			}
			par.tini-=vars->dt;
		};

		void initial(Variables *Vars, double rho){
			vars=Vars;
			for(auto &a : vars->particles){
				a.Urand.x[0]=0;
				a.Urand.x[1]=0;
				a.Urand.x[2]=0;
				a.tini=-1;
				update(a);
			}
			cout<<name<<endl;
		};

		TD_on(void){};
		~TD_on(void){};

	private:

};
