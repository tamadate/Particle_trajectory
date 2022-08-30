#include "trajectory.hpp"

void
trajectory::run(void){
	timer=omp_get_wtime();
	const int particleSize=vars->particles.size();
	int trapParticle=0;	// Number of cases when the calculaiton is not finished in required time steps

	//# pragma omp parallel for 
	for (int pid=0; pid<particleSize; pid++){
		particle &a=vars->particles[pid];
		initialize(a);

		while(vars->time<totalTime){
			if(vars->time - vars->preOutTime > observeTime) {
				output(a, vars->time);
				vars->preOutTime=vars->time;
			}

			timeEvolution(a);

			flags->breakFlag=checkCell(pid);
			if(flags->breakFlag==-1) break;

			if(a.update==1) updateDisp(a);	// Update dispersion?
		}

		// if particle did not hit on any "bounday" during wall time
		if(flags->breakFlag==0){
			outParticles[pid].pid=a.id;
			outParticles[pid].r=a.x;
			outParticles[pid].v=a.v;
			outParticles[pid].bid=-1;
			//trapParticle++;
		}
	}

	for (int pid=0; pid<particleSize; pid++) outputFinalPosition(outParticles[pid]);
  cout<<trapParticle<<" might be trapped circulation"<<endl;
	cout<<"Trajectory calculation time: "<<omp_get_wtime()-timer<<" sec"<<endl;
}
