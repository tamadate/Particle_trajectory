#pragma once
#include "solver.hpp"


// base function of force
class Hybrid : public Solver{
	public:
		double timeFactor;

		void solve(particle &a){
			if(a.Re<0.1 && a.Mach<0.1) {analytical(a);}
			else {Euler(a);}
		}
		void analytical(particle &a){
			double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
			double vmag=sqrt(v2);
			vars->dt=vars->meshScale*vars->analFactor/vmag;

			if(dtMax<vars->dt) vars->dt=dtMax ;

			if(a.tini<vars->dt) vars->dt=a.tini;

			// update dispersion
			td->update(a);

			for(int i=0; i<3; i++) a.F.x[i]=0;
			for(auto &force : forces) force->compute(a);

			// fluid velocity including eddy velocity
			double Ux = vars->U[a.cell].x[0] + a.Urand.x[0];
			double Uy = vars->U[a.cell].x[1] + a.Urand.x[1];
			double Uz = vars->U[a.cell].x[2] + a.Urand.x[2];

			// relative velocity
			double dvx = a.v.x[0] - Ux;
			double dvy = a.v.x[1] - Uy;
			double dvz = a.v.x[2] - Uz;

			double EXP = exp(-vars->dt * a.beta);
			double C1 = (1 - EXP) / a.beta;
			double C2 = (vars->dt - C1) / a.beta;

			// update velocity
			a.v.x[0] = Ux + dvx * EXP + a.F.x[0] * C1;
			a.v.x[1] = Uy + dvy * EXP + a.F.x[1] * C1;
			a.v.x[2] = Uz + dvz * EXP + a.F.x[2] * C1;

			// update position
			a.x.x[0] += Ux * vars->dt + dvx * C1 + a.F.x[0] * C2;
			a.x.x[1] += Uy * vars->dt + dvy * C1 + a.F.x[1] * C2;
			a.x.x[2] += Uz * vars->dt + dvz * C1 + a.F.x[2] * C2;
		};


		Hybrid(void){
			drag = new dragForceSM();
			td = new TD();
		}
		~Hybrid(void){}
	private:
};
