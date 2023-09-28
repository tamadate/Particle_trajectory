#pragma once
#include "force.hpp"



// base function of force
class coul : public force{
	public:
		Variables *vars;
		Flags *flags;
		double q;
		virtual void compute(particle &par){
			par.F.x[0]+=vars->E[par.cell].x[0]*q;
			par.F.x[1]+=vars->E[par.cell].x[1]*q;
			par.F.x[2]+=vars->E[par.cell].x[2]*q;

		};
		coul(double z){
			q=z*e;
		}
		~coul(void){}
	private:
};
