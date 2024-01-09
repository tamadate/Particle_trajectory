#pragma once
#include "force.hpp"



// base function of force
class Coulomb : public force{
	public:
		void compute(particle &par){
			double Ex=vars->dV[par.cell].x[0];
			double Ey=vars->dV[par.cell].x[1];
			double Ez=vars->dV[par.cell].x[2];
			par.F.x[0] += vars->dV[par.cell].x[0] * e / par.m;
			par.F.x[1] += vars->dV[par.cell].x[1] * e / par.m;
			par.F.x[2] += vars->dV[par.cell].x[2] * e / par.m;

		};
		Coulomb(void){}
		~Coulomb(void){}
	private:
};
