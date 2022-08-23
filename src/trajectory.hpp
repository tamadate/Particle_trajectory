#pragma onece
#include "functions.hpp"
#include "output.hpp"


void trajectoryEuler(Variables *vars){
    int loop=0;
	while(loop<timestep){
		forceUpdate(vars);
		velocityUpdateEuler(vars);
		positionUpdateEuler(vars);

		if(loop%Observe==0){
			vars->time=dt*loop;
			output(vars);
		}
        loop++;
	}
    cout<<"Final time: "<<dt*loop<<" s"<<endl;    
}


void velocityUpdateEuler(Variables *vars){
	for (auto &a: vars->particles){
		for(int i=0; i<3; i++){
			a.v.x[i]+=dt*a.F.x[i];
		}
	}
}

void positionUpdateEuler(Variables *vars){
	for (auto &a: vars->particles){
		for(int i=0; i<3; i++){
			a.x.x[i]+=dt*a.v.x[i];
		}
	}
}

void forceUpdate(Variables *vars){
	for (auto &a: vars->particles){
		double FD=computeFD(vars, a);
		for(int i=0; i<3; i++){
			a.F.x[i]=FD*(vars->U.x[i]-a.v.x[i]);
		}
	}
}

