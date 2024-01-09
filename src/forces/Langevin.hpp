#pragma once
#include "force.hpp"



// base function of force
class Langevin : public force{
	public:
		void compute(particle &par){
			// Random number generation setup
			std::random_device rd;  // Random device
			std::mt19937 gen(rd()); // Mersenne Twister generator
			std::normal_distribution<> d(0, 1); // Normal distribution with mean=0 and std=sigma
			double coeff = sqrt(2 * kb * vars->T[par.cell] * par.fric / vars->dt) / par.m;
			par.F.x[0] += d(gen) * coeff;
			par.F.x[1] += d(gen) * coeff;
			par.F.x[2] += d(gen) * coeff;

		};
		Langevin(void){}
		~Langevin(void){}
	private:
};
