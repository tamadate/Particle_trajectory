#pragma once
#include "solver.hpp"


// base function of force
class EulerIL : public Solver{
	public:

		void solve(particle &a){
			// calculate time step
			double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
			double vmag=sqrt(v2);
			double v2rand=a.Urand.x[0]*a.Urand.x[0]+a.Urand.x[1]*a.Urand.x[1]+a.Urand.x[2]*a.Urand.x[2];
			if(v2rand>v2) vmag=sqrt(v2rand);
			vars->dt=vars->meshScale/vmag;		
			if(dtMax<vars->dt) vars->dt=dtMax;
			if(a.tini<vars->dt)	vars->dt=a.tini;

			// update dispersion
			td->update(a);

			// compute forces
			for(int i=0; i<3; i++) a.F.x[i]=0;
			for(auto &force : forces) force->compute(a);
			
			// velocity update
			a.v.x[0]=vars->U[a.cell].x[0] + a.Urand.x[0];//+a.F.x[i]/a.beta;	
			a.v.x[1]=vars->U[a.cell].x[1] + a.Urand.x[1];//+a.F.x[i]/a.beta;	
			a.v.x[2]=vars->U[a.cell].x[2] + a.Urand.x[2];//+a.F.x[i]/a.beta;	

			// position update
			a.x.x[0]+=vars->dt*a.v.x[0];
			a.x.x[1]+=vars->dt*a.v.x[1];
			a.x.x[2]+=vars->dt*a.v.x[2];
		};

		EulerIL(void){
			drag = new dragForceSM();
			td = new TD();
		}
		~EulerIL(void){}
	private:
};
