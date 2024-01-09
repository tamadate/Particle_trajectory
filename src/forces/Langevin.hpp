#pragma once
#include "force.hpp"



// base function of force
class Langevin : public force{
	public:
		std::vector<double> randoms;
		const int N=10000000;
		int n;
		void generateRandom(void){
			randoms.clear();
			// Random number generation setup
			std::random_device rd;  // Random device
			std::mt19937 gen(rd()); // Mersenne Twister generator
			std::normal_distribution<> d(0, 1); // Normal distribution with mean=0 and std=sigma
			for(int i=0;i<N*3;i++) randoms.push_back(d(gen));
			n=0;
		}
		void compute(particle &par){
			double coeff = sqrt(2 * kb * vars->T[par.cell] * par.fric / vars->dt) / par.m;
			if(n>N-1) generateRandom();
			par.F.x[0] += randoms[3*n] * coeff;
			par.F.x[1] += randoms[3*n+1] * coeff;
			//par.F.x[2] += randoms[3*n+2] * coeff;
			n++;

		};
		Langevin(void){generateRandom();}
		~Langevin(void){}
	private:
};
