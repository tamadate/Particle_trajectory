#include "../trajectory.hpp"



trajectory::trajectory(void){
// set default values
    setDefault();

// generate class
    vars = new Variables(); // generate variables class
    flags = new Flags();    // initialize flags (see flags.hpp)
    solver = new Solver();
  
// initialization
    for(int i=0;i<vars->Nth;i++) trapParticle.push_back(0);
    readGeometry();   // read geometry ./constant/polyMesh/
    readCondition();  // read condition file ./particle/condition

    calculateMyu();   // viscosity calculation from field data through Sutherland's equation
    calculatelamda(); // mean free path calculation form field data

    readParticles();  // read particle file ./particle/particleSet
    makeCells();      // make cells data (e.g., norm vector) from geometry file

    initialParticle();// initialize particle (initial cell id, initial velocity)
    solver->initial(vars,rho_p);
    outputInitial();  // initialization of output files
}

void
trajectory::setDefault(void){
    rho_p=1000; // particle density (variable)
    observeTime=1e-5; // output time interval
    totalTime=1;    // total calculation time
    Observe=10000;	// output time steps
    delta_r2_trajectory=1e-4*1e-4;	// output migration distance for the trajectory
    startDir="20000"; // CFD simulation directory
}




