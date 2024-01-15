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
		std::vector<double> time;

		
		double delta_r2;
		double dt;

		double rho_p;									// particle density
		double observeTime;						// output interval in time scale
		double totalTime;							// total calculation time
		int Observe;									// output interval in time step
		double delta_r2_trajectory;		// squre of migration distance for trajectory output
		int plane2D;									// face of 2D simulation (0:y-z, 1:x-z, 2:x-y)
		int noUpdateAxis;
		double Axis;									// axis of 2D axi-symmetric simulation
		char filepath[100];						// file path
		int boundaryStartID;					// boundary start face id
		double timer;

		Variables(void){};
		~Variables(void){};
};
