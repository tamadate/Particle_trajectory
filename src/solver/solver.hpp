#pragma once
#include "forces/force.hpp"
#include "forces/drag/drag.hpp"
#include "forces/gravity.hpp"
#include "forces/Langevin.hpp"
#include "forces/Brownian.hpp"
#include "forces/Coulomb.hpp"
#include "turbDisp/TD_on.hpp"



// base function of force
class Solver{
	public:
		Variables *vars;
		TD *td;
		std::vector<force*> forces;	// pointer array for force functions working on the particle
		force *drag;	// drag force

		// virtual functions
		virtual void solve(particle &a)	{Euler(a);};

		// shared functions
		void Euler(particle &a){
			// calculate time step
			double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
			double vmag=sqrt(v2);
			double v2rand=a.Urand.x[0]*a.Urand.x[0]+a.Urand.x[1]*a.Urand.x[1]+a.Urand.x[2]*a.Urand.x[2];
			if(v2rand>v2) vmag=sqrt(v2rand);
			vars->dt=vars->meshScale/vmag;	
			//vars->dt=vars->fixTimeStep;	

			if(dtMax<vars->dt) vars->dt=dtMax ;
			if(a.dt<vars->dt) vars->dt=a.dt;
			if(a.tini<vars->dt)	vars->dt=a.tini;

			// update dispersion
			td->update(a);
			
			// compute forces
			for(int i=0; i<3; i++) a.F.x[i]=0;
			drag->compute(a);
			for(auto &force : forces) force->compute(a);
			
			// velocity update
			for(int i=0; i<3; i++) a.v.x[i]+=vars->dt*a.F.x[i];	

			// position update
			for(int i=0; i<3; i++) a.x.x[i]+=vars->dt*a.v.x[i];
		};

		void initial(Variables *Vars, double rho_p){
			vars=Vars;
			drag->initial(vars);
			for(auto &force : forces) force->initial(vars);  
			td->initial(vars,rho_p);
		};

		Solver(void){
			drag = new dragForceSM();
			td = new TD();
		};
		~Solver(void){};
	private:
};
