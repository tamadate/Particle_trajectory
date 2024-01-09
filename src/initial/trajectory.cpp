#include "../trajectory.hpp"



trajectory::trajectory(void){
// Get number of thread for OpenMP (parallel)
    #pragma omp parallel
        {
            #pragma omp single
            {
                Nth=omp_get_num_threads();
            }
    }

// set default values
    setDefault();
    for(int i=0;i<Nth;i++) trapParticle.push_back(0);

// generate class
    vars = new Variables(); // generate variables class
    flags = new Flags();    // initialize flags (see flags.hpp)
    //drag = new dragForceSM();
  
// initialization
    readGeometry();   // read geometry ./constant/polyMesh/
    readCondition();  // read condition file ./particle/condition

    calculateMyu();   // viscosity calculation from field data through Sutherland's equation
    calculatelamda(); // mean free path calculation form field data

    readParticles();  // read particle file ./particle/particleSet
    makeCells();      // make cells data (e.g., norm vector) from geometry file
    //drag->initial(vars,flags);
    for(auto &force : forces) force->initial(vars,flags);   
    initialParticle();// initialize particle (initial cell id, initial velocity)
    outputInitial();  // initialization of output files
}

void
trajectory::setDefault(void){
    rho_p=1000; // particle density (variable)
    observeTime=1e-5; // output time interval
    totalTime=1;    // total calculation time
    Observe=10000;	// output time steps
    delta_r2_trajectory=1e-4*1e-4;	// output migration distance for the trajectory
    plane2D=-1; //
    startDir="20000"; // CFD simulation directory
    Axis=0; // 0 for 3D simulaiton
}




