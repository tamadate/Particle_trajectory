#include "trajectory.hpp"

void
trajectory::output(particle a, double time){
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

void
trajectory::outputTrajectory(particle a){
	sprintf(filepath, "result/trajectory.%d", int(a.id));
	f=fopen(filepath, "a");
	fprintf(f,"%e\t%e\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2]);
	fclose(f);

}


void
trajectory::outputInitial(void){
	for (auto &a:vars->particles){
		sprintf(filepath, "result/position.%d", int(a.id));
		f=fopen(filepath, "w");
		fprintf(f,"%e\t%e\t%e\t%d\n", a.x.x[0], a.x.x[1], a.x.x[2], a.cell);
		fclose(f);

		sprintf(filepath, "result/U.%d", int(a.id));
		f=fopen(filepath, "w");
		fprintf(f,"%e\t%e\t%e\n", a.v.x[0], a.v.x[1], a.v.x[2]);
		fclose(f);

		sprintf(filepath, "result/trajectory.%d", int(a.id));
		f=fopen(filepath, "w");
		fprintf(f,"%e\t%e\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2]);
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
	sprintf(filepath, "initialVelocity.dat");
	f=fopen(filepath, "w");
	fclose(f);
}

void
trajectory::outputFinalPosition(outParticle a){
	sprintf(filepath, "finalPosition.dat");
	f=fopen(filepath, "a");
   	fprintf(f,"%e\t%e\t%e\t%d\t%d\n", a.r.x[0], a.r.x[1], a.r.x[2], a.pid, a.bid);
	fclose(f);
	sprintf(filepath, "finalVelocity.dat");
	f=fopen(filepath, "a");
   	fprintf(f,"%e\t%e\t%e\t%d\t%d\n", a.v.x[0], a.v.x[1], a.v.x[2], a.pid, a.bid);
	fclose(f);
}

void
trajectory::outputInitial(particle p){
	sprintf(filepath, "initialVelocity.dat");
	f=fopen(filepath, "a");
   	fprintf(f,"%e\t%e\t%e\t%d\n", p.v.x[0], p.v.x[1], p.v.x[2], p.id);
	fclose(f);
}

void
trajectory::outputPenetration(particle p){
	sprintf(filepath, "penetratePosition.dat");
	f=fopen(filepath, "a");
   	fprintf(f,"%e\t%e\t%e\t%d\n", p.x.x[0], p.x.x[1], p.x.x[2], p.id);
	fclose(f);
	sprintf(filepath, "penetrateVelocity.dat");
	f=fopen(filepath, "a");
   	fprintf(f,"%e\t%e\t%e\t%d\n", p.v.x[0], p.v.x[1], p.v.x[2], p.id);
	fclose(f);
}
