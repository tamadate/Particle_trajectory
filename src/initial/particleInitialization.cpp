#include "../trajectory.hpp"

// initialization of the particle cell id and setting initial velocity
void
trajectory::initialParticle(void){

    findParticle();  // find initial cell index of each particle

    int ps=vars->particles.size(); // total number of particles
    for(auto &a:vars->particles){
        int icell=a.cell;   // cell index
        a.dt=1/(2*a.beta);    // time step
        for(int i=0; i<3; i++) a.v.x[i]=vars->U[icell].x[i]*0.99;   // initial velocity is 99% of fluid velocity
        if(flags->v0) setInitialVelocity(); // if v0 flag is, set initial velocity through velocity file
        calculateNonDimension(a); // calculate non dimensional parameters (Re, Kn, and Mach numbers)
        for(auto &force : forces) force->compute(a);   // force calculation

        double threePiMuDp_Cc=3*M_PI*vars->myu[icell]*a.dp/a.Cc;    // friction factor
        a.Zp=1.6e-19/threePiMuDp_Cc;    // electrical mobility
        a.beta=threePiMuDp_Cc/a.m;  // beta relaxation time (Stokes)

        // initialize outParticles array (particle profiles at the end of the trajectory calculation)
        outParticle op;
        op.pid=a.id;
        op.bid=0;   // initialize boundary id 
        point dum; // dummy, dum=(0,0,0)
        dum.x[0]=dum.x[1]=dum.x[2]=0;
        op.r=dum;
        op.v=dum;
        outParticles.push_back(op);

        // initialize outParticles array for the intermediate profiles
        // (particle profiles when the particle penetrate a surface indicated in input file)
        for(auto &pen:penetrates){
            pen.outPositions.push_back(op);
            pen.dx0.push_back(pen.loc-a.x.x[pen.face]);
        }
    }
}


// set initial particle velocities
void
trajectory::setInitialVelocity(void){
    // if initial velocity file does exist
    if(is_file_exist("particle/particleVelocity")){
        string str;
        bool iflag=0;
        int particleID=0;
        ifstream stream("particle/particleVelocity");
        while(getline(stream,str)) {
            if (str=="vx\tvy\tvz") {iflag=1; continue;} // detect header
            if (iflag==1){
                int loop=0;
                string tmp;
                istringstream stream(str);
                while(getline(stream,tmp,'\t')) {
                    vars->particles[particleID].v.x[loop]=stod(tmp);
                    loop++;
                }
                particleID++;
            }
        }
        stream.close();

        // check if all particle velocities are set properly
        int N=vars->particles.size();
        if (N!=particleID) cout<<"lenght of particle/particleVelocity file does not matched"<<endl;
    }
    else{
        cout<<"particle/particleVelocity file does not exist"<<endl;
    }
}