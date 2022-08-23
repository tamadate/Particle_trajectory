#include "functions.hpp"
#include "trajectory.hpp"
#include "dragForce.hpp"

int main(void){

	//Generate variable class
	Variables *vars;
	vars = new Variables();

	//set initial position of particles
	//3 particles(0.1, 1, 10 um) at (x,y,z)=(0,0,0)
	particle par;
	par.x.x[0]=par.x.x[1]=par.x.x[2]=0;
	par.id=0;
	par.dp=1e-6;
	vars->particles.push_back(par);

	par.x.x[0]=par.x.x[1]=par.x.x[2]=0;
	par.id=1;
	par.dp=10e-6;
	vars->particles.push_back(par);

	par.x.x[0]=par.x.x[1]=par.x.x[2]=0;
	par.id=2;
	par.dp=0.1e-6;
	vars->particles.push_back(par);

	initialParticle(vars);	//Initialization of particle properties (e.g., velocity)
	outputInitial(vars);	//Initialization of output files
	trajectoryEuler(vars);	//trajectory calculation
	
	return 0;
}





