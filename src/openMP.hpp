#pragma once
#include "variables.hpp"

class VariablesMP {
	private:
		int Nth;
  	public:
		double delta_r2;
		std::vector<double> time;
		double meshScale;
		double analFactor;
		double dt;
		double fixTimeStep;

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

		VariablesMP(void){};
		~VariablesMP(void){};
};
