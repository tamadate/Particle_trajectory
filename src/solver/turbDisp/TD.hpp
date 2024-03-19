#pragma once
#include "../../functions.hpp"
#include "../../variables.hpp"

// without turbulent dispersion
class TD{
	public:
		Variables *vars;

		virtual void tini(particle &par){};
		virtual void update(particle &par){};
		virtual void initial(Variables *Vars, double rho){
			vars=Vars;
			for(auto &a : vars->particles){
				a.Urand.x[0]=0;
				a.Urand.x[1]=0;
				a.Urand.x[2]=0;
				a.tini=10;
				update(a);
			}
		};
		TD(void){};
		~TD(void){};
	private:
};