#pragma once

#include "functions.hpp"


void output(particle a, double time){
	sprintf(filepath, "result/position.%d", int(a.id)); 
	f=fopen(filepath, "a");
	fprintf(f,"%e\t%e\t%e\t%d\n", a.x.x[0], a.x.x[1], a.x.x[2], a.cell);
	fclose(f);

	sprintf(filepath, "result/U.%d", int(a.id)); 
	f=fopen(filepath, "a");
	fprintf(f,"%e\t%e\t%e\n", a.v.x[0], a.v.x[1], a.v.x[2]);
	fclose(f);

	sprintf(filepath, "result/Time"); 
	f=fopen(filepath, "a");
	fprintf(f,"%e\n", time);
	fclose(f);

}


void outputInitial(Variables *vars){
	for (auto &a:vars->particles){
		sprintf(filepath, "result/position.%d", int(a.id)); 
		f=fopen(filepath, "w");
		fclose(f);

		sprintf(filepath, "result/U.%d", int(a.id)); 
		f=fopen(filepath, "w");
		fclose(f);
	}
	sprintf(filepath, "result/Time"); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "out.dat"); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "result/nparticle"); 
	f=fopen(filepath, "w");
	fprintf(f,"%d\n", int(vars->particles.size()));
	fclose(f);
	sprintf(filepath, "finalPosition.dat"); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "finalVelocity.dat"); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "penetratePosition.dat"); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "penetrateVelocity.dat"); 
	f=fopen(filepath, "w");
	fclose(f);
}

void outputFinalPosition(outParticle a){
	sprintf(filepath, "finalPosition.dat"); 
	f=fopen(filepath, "a");
   	fprintf(f,"%f\t%f\t%f\t%d\t%d\n", a.r.x[0], a.r.x[1], a.r.x[2], a.pid, a.bid);
	fclose(f);
	sprintf(filepath, "finalVelocity.dat"); 
	f=fopen(filepath, "a");
   	fprintf(f,"%f\t%f\t%f\t%d\t%d\n", a.v.x[0], a.v.x[1], a.v.x[2], a.pid, a.bid);
	fclose(f);
}

void outputFinalPositionInitial(){

}

void outputPenetration(outParticle a){
	sprintf(filepath, "penetratePosition.dat"); 
	f=fopen(filepath, "a");
   	fprintf(f,"%f\t%f\t%f\t%d\t%d\n", a.r.x[0], a.r.x[1], a.r.x[2], a.pid, a.bid);
	fclose(f);
	sprintf(filepath, "penetrateVelocity.dat"); 
	f=fopen(filepath, "a");
   	fprintf(f,"%f\t%f\t%f\t%d\t%d\n", a.v.x[0], a.v.x[1], a.v.x[2], a.pid, a.bid);
	fclose(f);
}

void outputPenetrationInitial(){

}

