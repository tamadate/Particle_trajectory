#include "trajectory.hpp"


trajectory::trajectory(void){

  #pragma omp parallel
  {
    #pragma omp single
    {
      Nth=omp_get_num_threads();
    }
  }


  rho_p=1000;
  observeTime=1e-5;
  totalTime=1;
  Observe=10000;	// output time steps
  delta_r2_trajectory=1e-4*1e-4;	// output migration distance for the trajectory
  plane2D=-1;
  startDir="20000";
  Axis=0;
  flag=1;
  filepath[100];
  constTemp=300;
  constRho=1.2;

  vars = new Variables();
  flags = new Flags();

  readCondition();
  readGeometry();
  readCFDresults();

  calculateMyu();
  calculatelamda();

  readParticles();
  makeCells();
  initialParticle();
  outputInitial();
}

void
trajectory::initialParticle(void){
	int ps=vars->particles.size();
	findParticle();
	for(auto &a:vars->particles){
		int icell=a.cell;
		double Kn=vars->lamda[icell]/a.dp;
		double Cc=1+Kn*(A1+A2*exp(-A3/Kn));
		a.Kn=Kn;
		a.Cc=Cc;
		double threePiMuDp_Cc=3*M_PI*vars->myu[icell]*a.dp/Cc;
		a.Zp=1.6e-19/threePiMuDp_Cc;
		a.beta=threePiMuDp_Cc/a.m;
		a.dt=1/(40*a.beta);
		for(int i=0; i<3; i++) {
			a.v.x[i]=vars->U[icell].x[i]*0.99;
		}
		for(auto &force : forces) force->computeFD(vars,flags,a);
    outParticle op;
    op.pid=a.id;
    outParticles.push_back(op);
	}
}
