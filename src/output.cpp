#include "trajectory.hpp"

void
trajectory::output(particle a, double time, int nth){
	FILE*fMP;

	// output particle position
	sprintf(filepathMP[nth], "result/position.%d", int(a.id));
	fMP=fopen(filepathMP[nth], "a");
	fprintf(fMP,"%e\t%e\t%e\t%d\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2], a.cell, time);
	fclose(fMP);

	// output particle velocity
	sprintf(filepathMP[nth], "result/U.%d", int(a.id));
	fMP=fopen(filepathMP[nth], "a");
	fprintf(fMP,"%e\t%e\t%e\t%e\n", a.v.x[0], a.v.x[1], a.v.x[2], time);
	fclose(fMP);

}

void
trajectory::outputTrajectory(particle a){
	// output particle position (output interval is controled by the migration distance but not time)
	sprintf(filepathMP[0], "result/trajectory.%d", int(a.id));
	FILE*fMP=fopen(filepathMP[0], "a");
	fprintf(fMP,"%e\t%e\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2]);
	fclose(fMP);
}

void
trajectory::outputInitial(void){
	for (auto &a:vars->particles){
		// particle position
		sprintf(filepath, "result/position.%d", int(a.id));
		f=fopen(filepath, "w");
		fprintf(f,"%e\t%e\t%e\t%d\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2], a.cell,0.0);
		fclose(f);

		// particle velocity
		sprintf(filepath, "result/U.%d", int(a.id));
		f=fopen(filepath, "w");
		fprintf(f,"%e\t%e\t%e\t%e\n", a.v.x[0], a.v.x[1], a.v.x[2],0.0);
		fclose(f);

		// particle position (output interval is controled by the migration distance but not time)
		sprintf(filepath, "result/trajectory.%d", int(a.id));
		f=fopen(filepath, "w");
		fprintf(f,"%e\t%e\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2]);
		fclose(f);
	}

	// number of particles
	sprintf(filepath, "result/nparticle");
	f=fopen(filepath, "w");
	fprintf(f,"%d\n", int(vars->particles.size()));
	fclose(f);
	// particle final position
	sprintf(filepath, "finalPosition.dat");
	f=fopen(filepath, "w");
	fclose(f);
	// particle final velocity
	sprintf(filepath, "finalVelocity.dat");
	f=fopen(filepath, "w");
	fclose(f);
	// particle position at arbitral point
	sprintf(filepath, "penetratePosition.dat");
	f=fopen(filepath, "w");
	fclose(f);
	// particle velocity at arbitral point
	sprintf(filepath, "penetrateVelocity.dat");
	f=fopen(filepath, "w");
	fclose(f);
	// particle initial velocity
	sprintf(filepath, "initialVelocity.dat");
	f=fopen(filepath, "w");
	fclose(f);
}

void
trajectory::outputFinalPosition(void){
	for (auto &a:outParticles){
		// particle final position
		sprintf(filepath, "finalPosition.dat");
		f=fopen(filepath, "a");
	 	fprintf(f,"%e\t%e\t%e\t%d\t%d\n", a.r.x[0], a.r.x[1], a.r.x[2], a.pid, a.bid);
		fclose(f);
		// particle final velocity
		sprintf(filepath, "finalVelocity.dat");
		f=fopen(filepath, "a");
	 	fprintf(f,"%e\t%e\t%e\t%d\t%d\n", a.v.x[0], a.v.x[1], a.v.x[2], a.pid, a.bid);
		fclose(f);
	}
}


void
trajectory::outputPenetration(particle p){
	char filepathMP[100];
	FILE*fMP;

	sprintf(filepathMP, "penetratePosition.dat");
	fMP=fopen(filepathMP, "a");
 	fprintf(fMP,"%e\t%e\t%e\t%d\n", p.x.x[0], p.x.x[1], p.x.x[2], p.id);
	fclose(fMP);

	sprintf(filepathMP, "penetrateVelocity.dat");
	fMP=fopen(filepathMP, "a");
 	fprintf(fMP,"%e\t%e\t%e\t%d\n", p.v.x[0], p.v.x[1], p.v.x[2], p.id);
	fclose(fMP);
}
