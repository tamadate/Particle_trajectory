#include "trajectory.hpp"

void
trajectory::run(void){
	timer=omp_get_wtime();
	const int particleSize=vars->particles.size();
	int trapParticle=0;	// Number of cases when the calculaiton is not finished in required time steps

	# pragma omp parallel for
	for (int pid=0; pid<particleSize; pid++){
		int nth=omp_get_thread_num();
		particle &a=vars->particles[pid];
		initialize(a,nth);

		int breakFlag=0;
		double preOutTime=0;
		double time=0;
		while(time<totalTime){
			if(time - preOutTime > observeTime) {
				output(a, time);
				preOutTime=time;
			}

			time+=timeEvolution(a);

			breakFlag=checkCell(pid);
			if(breakFlag==-1) break;

			if(a.update==1) updateDisp(a);	// Update dispersion?
		}

		// if particle did not hit on any "bounday" during wall time
		if(breakFlag==0){
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
