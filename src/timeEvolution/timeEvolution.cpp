#include "../trajectory.hpp"


/*double
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
}*/

double
trajectory::timeEvolution(particle &a){

	calculateNonDimension(a); // Reynolds, Mach, and Kndsen numbers.

	// compute forces
	for(int i=0; i<3; i++) a.F.x[i]=0;
	for(auto &force : forces) force->compute(a);
	double absF = a.F.x[0]*a.F.x[0]+a.F.x[1]*a.F.x[1]+a.F.x[2]*a.F.x[2];

	// calculate time step
	double fric=3*M_PI*vars->myu[a.cell]*a.dp/a.Cc;
	double timeStep=a.dp*fric/absF;
	double t0=a.dp*a.dp*fric/(6*kb*vars->T[a.cell]);
	if(t0<timeStep) timeStep=t0;

	// check time step
	double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
	double vmag=sqrt(v2);
	if(vars->meshScale/vmag<timeStep) cout<<"time step not correct"<<endl;

	double vt[3];
	vt[0]=a.v.x[0];
	vt[1]=a.v.x[1];
	vt[2]=a.v.x[2];


    std::random_device rd;
    std::mt19937 generator(rd());
	double expTerm=exp(-a.beta*timeStep);
	double kT=kb*vars->T[a.cell];
	double RRvsq=sqrt(kT/a.m*(1-expTerm*expTerm));
	std::normal_distribution<double> dist(0.0, 1.0);

	// velocity update
	for(int i=0; i<3; i++) {
		a.v.x[i]+=vt[i]*expTerm;
		a.v.x[i]+=a.F.x[i]/fric*(1-expTerm);
		a.v.x[i]+=RRvsq*dist(generator);
	}

	double C1=(1-expTerm)/(1+expTerm);
	double RRxsq=sqrt(2*a.m*kT/fric/fric*(fric*timeStep/a.m-2*C1));
	// position update
	for(int i=0; i<3; i++) {
		a.x.x[i]+=a.m/fric*(a.v.x[i]+vt[i]-2*a.F.x[i]/fric)*C1;
		a.x.x[i]+=a.F.x[i]/fric*timeStep;
		a.x.x[i]+=RRxsq*dist(generator);
	}

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
