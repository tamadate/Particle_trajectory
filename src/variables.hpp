#pragma once
#include "functions.hpp"

class Variables {
	private:
  	public:
		int Nth;
		std::vector<double> p;
		std::vector<point> U;
		std::vector<double> T;
		std::vector<double> rho;
		std::vector<double> myu;
		std::vector<double> lamda;
		std::vector<double> k;
		std::vector<double> epsilon;
		std::vector<point> dp;
		std::vector<point> dT;
		std::vector<point> dV;
		std::vector<particle> particles;
		double delta_r2;
		std::vector<double> time;
		double meshScale;
		double analFactor;
		double dt;
		double fixTimeStep;

		Variables(void){
			#pragma omp parallel
			{
				#pragma omp single
				{
					Nth=omp_get_num_threads();
					cout << "OMP threads = " << Nth << endl;
				}
			}
		};
		~Variables(void){};
};
