#pragma once
#include "force.hpp"



// base function of force
class Brownian : public force{
	public:
		std::vector<std::vector<double>> randoms;
		const int N=10000000;
		std::vector<int> n;
		double Sch;
		void generateRandom(void){
			int nth=omp_get_thread_num();
			randoms[nth].clear();
			// Random number generation setup
			std::random_device rd;  // Random device
			std::mt19937 gen(rd()); // Mersenne Twister generator
			std::normal_distribution<> d(0, 1); // Normal distribution with mean=0 and std=sigma
			for(int i=0;i<N*3;i++) randoms[nth].push_back(d(gen));
			n[nth]=0;
		};
		void compute(particle &par){
			int nth=omp_get_thread_num();
			double k=vars->k[par.cell];
			double mut=0.09*vars->rho[par.cell]*k*k/vars->epsilon[par.cell];
			double Dturb=mut/Sch;
			double coeff = sqrt(2 * Dturb / vars->dt) * par.beta;

			if(n[nth]>N-1) generateRandom();
			par.F.x[0] += randoms[nth][3*n[nth]] * coeff;
			par.F.x[1] += randoms[nth][3*n[nth]+1] * coeff;
			par.F.x[2] += randoms[nth][3*n[nth]+2] * coeff;
			n[nth]++;

		};
		Brownian(double sch){
			int N_threads=omp_get_max_threads();
			for (int i=0;i<N_threads;i++){
				std::vector<double> arrays;
				// Random number generation setup
				std::random_device rd;  // Random device
				std::mt19937 gen(rd()); // Mersenne Twister generator
				std::normal_distribution<> d(0, 1); // Normal distribution with mean=0 and std=sigma
				for(int i=0;i<N*3;i++) arrays.push_back(d(gen));
				randoms.push_back(arrays);
				n.push_back(1);	
			}
			Sch=sch;
		}
		~Brownian(void){}
	private:
};
