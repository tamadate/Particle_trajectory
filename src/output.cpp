#include "trajectory.hpp"

void
trajectory::output(particle a, double time, int nth){
	FILE*fMP;
	char filePathMP[100];

	// output particle position
	sprintf(filePathMP, "result/position.%d", int(a.id));
	fMP=fopen(filePathMP, "a");
	fprintf(fMP,"%e\t%e\t%e\t%d\t%e\n", a.x.x[0], a.x.x[1], a.x.x[2], a.cell, time);
	fclose(fMP);

	// output particle velocity
	sprintf(filePathMP, "result/U.%d", int(a.id));
	fMP=fopen(filePathMP, "a");
	fprintf(fMP,"%e\t%e\t%e\t%e\n", a.v.x[0], a.v.x[1], a.v.x[2], time);
	fclose(fMP);

}

void
trajectory::outputTrajectory(particle a){
	FILE*fMP;
	char filePathMP[100];
	// output particle position (output interval is controled by the migration distance but not time)
	sprintf(filePathMP, "result/trajectory.%d", int(a.id));
	fMP=fopen(filePathMP, "a");
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

	for(int i=0;i<penetrates.size();i++){
		// particle position at arbitral point
		sprintf(filepath, "penetratePosition_%d.dat",i);
		f=fopen(filepath, "w");
		fclose(f);
		// particle velocity at arbitral point
		sprintf(filepath, "penetrateVelocity_%d.dat",i);
		f=fopen(filepath, "w");
		fclose(f);
	}
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
trajectory::outputInitialVelocity(void){
	sprintf(filepath, "initialVelocity.dat");
	f=fopen(filepath, "w");
	for (auto &a:vars->particles){
		// particle final velocity
		for(int i=0; i<3; i++) a.v.x[i]=vars->U[a.cell].x[i]*0.99;
	 	fprintf(f,"%e\t%e\t%e\t%d\n", a.v.x[0], a.v.x[1], a.v.x[2], a.id);
	}
	fclose(f);
}


void
trajectory::outputPenetrate(void){
	char filePathMP[100];
	FILE*fMP;
	int i=0;
	for (auto &pen:penetrates){
		for (auto &p:pen.outPositions){
			sprintf(filePathMP, "penetratePosition_%d.dat",i);
			fMP=fopen(filePathMP, "a");
		 	fprintf(fMP,"%e\t%e\t%e\t%d\n", p.r.x[0], p.r.x[1], p.r.x[2], p.pid);
			fclose(fMP);

			sprintf(filePathMP, "penetrateVelocity_%d.dat",i);
			fMP=fopen(filePathMP, "a");
		 	fprintf(fMP,"%e\t%e\t%e\t%d\n", p.v.x[0], p.v.x[1], p.v.x[2], p.pid);
			fclose(fMP);
		}
		i++;
	}
}
