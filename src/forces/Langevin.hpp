#pragma once
#include "force.hpp"



// base function of force
class Langevin : public force{
	public:
		std::vector<std::vector<double>> randoms;
		const int N=10000000;
		std::vector<int> n;
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
			double coeff = sqrt(2 * kb * vars->T[par.cell] * par.fric / vars->dt) / par.m;
			if(n[nth]>N-1) generateRandom();
			par.F.x[0] += randoms[nth][3*n[nth]] * coeff;
			par.F.x[1] += randoms[nth][3*n[nth]+1] * coeff;
			//par.F.x[2] += randoms[3*n+2] * coeff;
			n[nth]++;

		};
		Langevin(void){
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
		}
		~Langevin(void){}
	private:
};
