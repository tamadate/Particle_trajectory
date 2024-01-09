#pragma once
#include "force.hpp"



// base function of force
class Coulomb : public force{
	public:
		const double coeff=e/4.0/M_PI/e0;
		void compute(particle &par){
			par.F.x[0] += vars->dV[par.cell].x[0] * coeff;
			par.F.x[1] += vars->dV[par.cell].x[1] * coeff;
			par.F.x[2] += vars->dV[par.cell].x[2] * coeff;

		};
		Coulomb(void){}
		~Coulomb(void){}
	private:
};
