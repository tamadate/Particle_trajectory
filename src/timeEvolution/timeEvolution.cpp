#include "../trajectory.hpp"


void
trajectory::timeEvolution(particle &a){

	calculateNonDimension(a); // Reynolds, Mach, and Kndsen numbers.

	if(flags->inertia){
		// use analytical or euler?
		if(a.Re<0.1 && a.Mach<0.1 && flags->analytical==1) analytical(a);
		else euler(a);
	}
	else{eulerInertiaLess(a);}

	if(a.dp>1e-6){
		double ddp=4*18e-3*0.22e-4/8.314/vars->T[a.cell]/rho_p/a.dp*(vars->pv_inf[a.cell]-vars->pvs[a.cell]);
		a.dp+=ddp*vars->dt;
	}

	// if the system is symmetric
	if(flags->dimensionFlag>0){
		a.x.x[noUpdateAxis]=0;
		if(a.x.x[plane2D]<Axis) {
			a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
			a.v.x[plane2D]*=-1.0;
			a.F.x[plane2D]*=-1.0;
			a.reflect*=-1;
		}
	}
	a.tini-=vars->dt;
}

void
trajectory::euler(particle &a){
	// calculate time step
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		double v2rand=a.Urand.x[0]*a.Urand.x[0]+a.Urand.x[1]*a.Urand.x[1]+a.Urand.x[2]*a.Urand.x[2];
		if(v2rand>v2) vmag=sqrt(v2rand);
		vars->dt=vars->meshScale/vmag;		
	}
	else vars->dt=vars->fixTimeStep;

	if(dtMax<vars->dt) vars->dt=dtMax ;

	if(a.dt<vars->dt) vars->dt=a.dt;
	

	// if the eddy reached to its life time (a.tini<timeStep),
	// the eddy is updated
	if(a.tini<vars->dt){
		vars->dt=a.tini;
		a.update=1;
	}
	
	// compute forces
	for(int i=0; i<3; i++) a.F.x[i]=0;
	drag->compute(a);
	for(auto &force : forces) force->compute(a);
	
	// velocity update
	for(int i=0; i<3; i++) a.v.x[i]+=vars->dt*a.F.x[i];	

	// position update
	for(int i=0; i<3; i++) a.x.x[i]+=vars->dt*a.v.x[i];
}


void
trajectory::analytical(particle &a){
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		vars->dt=vars->meshScale*vars->analFactor/vmag;
	}
	else vars->dt=vars->fixTimeStep;

	if(dtMax<vars->dt) vars->dt=dtMax ;

	if(a.tini<vars->dt){
		vars->dt=a.tini;
		a.update=1;
	}

	for(int i=0; i<3; i++) a.F.x[i]=0;
	for(auto &force : forces) force->compute(a);

	// fluid velocity including eddy velocity
	double Ux = vars->U[a.cell].x[0] + a.Urand.x[0];
	double Uy = vars->U[a.cell].x[1] + a.Urand.x[1];
	double Uz = vars->U[a.cell].x[2] + a.Urand.x[2];

	// relative velocity
	double dvx = a.v.x[0] - Ux;
	double dvy = a.v.x[1] - Uy;
	double dvz = a.v.x[2] - Uz;

	double EXP = exp(-vars->dt * a.beta);
	double C1 = (1 - EXP) / a.beta;
	double C2 = (vars->dt - C1) / a.beta;

	// update velocity
	a.v.x[0] = Ux + dvx * EXP + a.F.x[0] * C1;
	a.v.x[1] = Uy + dvy * EXP + a.F.x[1] * C1;
	a.v.x[2] = Uz + dvz * EXP + a.F.x[2] * C1;

	// update position
	a.x.x[0] += Ux * vars->dt + dvx * C1 + a.F.x[0] * C2;
	a.x.x[1] += Uy * vars->dt + dvy * C1 + a.F.x[1] * C2;
	a.x.x[2] += Uz * vars->dt + dvz * C1 + a.F.x[2] * C2;


}



void
trajectory::eulerInertiaLess(particle &a){
	// calculate time step
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		double v2rand=a.Urand.x[0]*a.Urand.x[0]+a.Urand.x[1]*a.Urand.x[1]+a.Urand.x[2]*a.Urand.x[2];
		if(v2rand>v2) vmag=sqrt(v2rand);
		vars->dt=vars->meshScale/vmag;		
	}
	else vars->dt=vars->fixTimeStep;

	if(dtMax<vars->dt) vars->dt=dtMax ;

	// if the eddy reached to its life time (a.tini<timeStep),
	// the eddy is updated
	if(a.tini<vars->dt){
		vars->dt=a.tini;
		a.update=1;
	}

	// compute forces
	for(int i=0; i<3; i++) a.F.x[i]=0;
	for(auto &force : forces) force->compute(a);
	
	// velocity update
	double U[3];
	U[0] = vars->U[a.cell].x[0] + a.Urand.x[0];
	U[1] = vars->U[a.cell].x[1] + a.Urand.x[1];
	U[2] = vars->U[a.cell].x[2] + a.Urand.x[2];


	for(int i=0; i<3; i++) a.v.x[i]=U[i]+a.F.x[i]/a.beta;	

	// position update
	for(int i=0; i<3; i++) a.x.x[i]+=vars->dt*a.v.x[i];
}