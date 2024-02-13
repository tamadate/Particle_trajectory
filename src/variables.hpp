#pragma once
#include "functions.hpp"

class Variables {
	private:
		int Nth;
  	public:
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
		std::vector<double> pvs;
		std::vector<double> pv_inf;


		double delta_r2;
		std::vector<double> time;
		double meshScale;
		double analFactor;
		double dt;
		double fixTimeStep;


		Variables(void){};
		~Variables(void){};
};
