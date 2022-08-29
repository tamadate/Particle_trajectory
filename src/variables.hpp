#pragma once
#include "functions.hpp"

class Variables {
	private:
    public:
		std::vector<double> p;
		std::vector<point> U;
		std::vector<double> T;
		std::vector<double> rho;
		std::vector<double> myu;
		std::vector<double> lamda;
		std::vector<point> dp;
		std::vector<particle> particles;
		double delta_r2;
		double time;
		double preOutTime;
		double analyticalStep;

		Variables(void){time=0;preOutTime=0;};
		~Variables(void){};
};
