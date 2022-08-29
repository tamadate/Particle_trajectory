#include "trajectory.hpp"

void
trajectory::run(void){
	timer=clock();
	const int particleSize=vars->particles.size();
	int trapParticle=0;	// Number of cases when the calculaiton is not finished in required time steps
//	# pragma omp parallel for
	for (int pid=0; pid<particleSize; pid++){
		particle &a=vars->particles[pid];
		initialize();
		outputInitial(a);

		while(vars->time<totalTime){
			if(vars->time - vars->preOutTime > observeTime) {
				output(a, vars->time);
				vars->preOutTime=vars->time;
			}
			computeReMach(a);
			if(a.Re<0.01 && a.Mach<0.1 && flags->analytical==1) analytical(a);
			else euler(a);

			flags->breakFlag=checkCell(pid);
			if(flags->breakFlag==-1) break;
		}
		if(flags->breakFlag==0){
			outParticle op;
			op.pid=a.id;
			op.r=a.x;
			op.v=a.v;
			op.bid=-1;
			outParticles.push_back(op);
			trapParticle++;
		}
		outputFinalPosition(outParticles[pid]);
	}
    cout<<trapParticle<<" might be trapped circulation"<<endl;
	cout<<"Trajectory calculation time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
}
