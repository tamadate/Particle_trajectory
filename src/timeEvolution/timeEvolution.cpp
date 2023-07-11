#include "../trajectory.hpp"


double
trajectory::timeEvolution(particle &a){

	calculateNonDimension(a); // Reynolds, Mach, and Kndsen numbers.

	double timeStep;
	// use analytical or euler?
	if(a.Re<0.1 && a.Mach<0.1 && flags->analytical==1) timeStep=analytical(a);
	else timeStep=euler(a);

	// if the system is symmetric
	if(a.x.x[plane2D]<Axis && flags->dimensionFlag>0) {
		a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
		a.v.x[plane2D]*=-1.0;
		a.F.x[plane2D]*=-1.0;
		a.reflect*=-1;
	}
	a.tini-=timeStep;


	return timeStep;
}

double
trajectory::euler(particle &a){
	double timeStep;
	// calculate time step
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		double v2rand=a.Urand.x[0]*a.Urand.x[0]+a.Urand.x[1]*a.Urand.x[1]+a.Urand.x[2]*a.Urand.x[2];
		if(v2rand>v2) vmag=sqrt(v2rand);
		timeStep=vars->meshScale/vmag;		
	}

	if(a.dt<timeStep) timeStep=a.dt;

	// if the eddy reached to its life time (a.tini<timeStep),
	// the eddy is updated
	if(a.tini<timeStep){
		timeStep=a.tini;
		a.update=1;
	}

	// compute forces
	for(int i=0; i<3; i++) a.F.x[i]=0;
	for(auto &force : forces) force->compute(a);

	// velocity update
	for(int i=0; i<3; i++) a.v.x[i]+=timeStep*a.F.x[i];

	// position update
	for(int i=0; i<3; i++) a.x.x[i]+=timeStep*a.v.x[i];

	return timeStep;
}


double
trajectory::analytical(particle &a){
	double timeStep;
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		timeStep=vars->meshScale*vars->analFactor/vmag;
	}

	if(dtMax<timeStep) timeStep=dtMax ;

	if(a.tini<timeStep){
		timeStep=a.tini;
		a.update=1;
	}

	// fluid velocity including eddy velocity
	double Ux=vars->U[a.cell].x[0]+a.Urand.x[0]+gAnal[0]/a.beta;;
	double Uy=vars->U[a.cell].x[1]+a.Urand.x[1]+gAnal[1]/a.beta;;
	double Uz=vars->U[a.cell].x[2]+a.Urand.x[2]+gAnal[2]/a.beta;
	// relative velocity
	double dvx=a.v.x[0]-Ux;
	double dvy=a.v.x[1]-Uy;
	double dvz=a.v.x[2]-Uz;

	// displacement of the velocities
	double EXP=exp(-timeStep*a.beta);
	double dvxEXP=dvx*EXP;
	double dvyEXP=dvy*EXP;
	double dvzEXP=dvz*EXP;

	// update velocity
	a.v.x[0] = Ux + dvx * EXP;
	a.v.x[1] = Uy + dvy * EXP;
	a.v.x[2] = Uz + dvz * EXP;

	// update position
	double coeff= (1 - EXP) / a.beta;
	a.x.x[0] += Ux * timeStep - dvx * coeff;
	a.x.x[1] += Uy * timeStep - dvy * coeff;
	a.x.x[2] += Uz * timeStep - dvz * coeff;

	return timeStep;
}
