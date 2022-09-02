# Particle_trajectory
This code calculate particle trajectory coupled with the CFD simulation result (only 1 way coupling is avairable in current version) from OpenFOAM (https://www.openfoam.com/).  Fluent results are also avairable by transforming to OpenFOAM fomat via "fluentToFoam" function in OpenFOAM.
## Usage
### 1. Building source code
Download this code and prepare the environment (gcc compiler and OpenMP are required). On the turminal, move to src directory from your downloaded directory and make code by just typing
~~~
make
~~~
Wait for several seconds and a execute file trajectory.out is created in the src directory.  
### 2. Prepare CFD simulation
Run OpenFOAM or Fluent.  Fluent format need to be transformed to the OpenFOAM format via fluentToFoam command in the OpenFOAM, which mean OpenFOAM installation is required both of the ways.
### 3. Set conditions and run
#### 3.1 Creat directories
In the CFD simulation directory, you must have four directories which are named 1.constant, 2.CFD result, 3.result, and 4.particle. First and second directories are from early performed CFD simulation, e.g., there is a directory named 20000, when the ending time of the simulation was 20000.  A system file may be also here but it does not matter if you leave or keep it since the trajectory simulation does not never use it.  On the third one, a blanck directory named "result" is required to store the trajectory simulation results (make sure creat it).  On the fourth directory, the name is needed to be "particle" which should have trajectory simulation conditions file (/particle/conditions) and particle information file (/particle/particleSet).
#### 3.2 conditions file in particle directory
This code read this file to set the trajectory simulaiton.  In this file, generally, setting parameter is mentioned at first and the setting values or method are followed it in a line (the tab is only allowed as a separation).  For example,
~~~
dragModel	Stokes-Millikan
~~~
set the drag model to Stokes-Millikan model.  All syntaxes are shown below.
1. Time step <br>
"dt" mention the time step and following syntax is auto or value.  Default is auto.  When auto is mentioned, this code automatically determine the time step from the particle veloctiy and relaxation time (or duration of the eddy when the turbulent dispersion is ON).  When value is mentioned, the time step is fixed as that value (unit is second).  Here is one example
~~~
dt  auto
~~~
2. Total calculation time <br>
"totalTime" mention the total simulation time (wall simulation time) and following syntax is value in the unit of second.  Ideally, this time is infinity since the particle exit from the calculation domain when the time is long enough.  However, the particles are potentially trapped by the curculation and it make simulation time longer.  This wall time is set to prevent the situation that the simulation is never done. Here is one example
~~~
totalTime  100
~~~
3. Drag model <br>
"dragModel" mention the drag model that is choosable from four drag models (Singh, Stokes-Millikan, Morsi, and Loth) and the syntax is same as the model name (default is Singh).  One of the feature of this code is Singh's drag model (https://arc.aiaa.org/doi/10.2514/1.J060648) and this is avairable for wide range of Reynodls number and Mach number even if it is supersonic flow, hence the original paper call this model as general particle drag model. Loth model is also avairable for the wide range (https://arc.aiaa.org/doi/10.2514/1.28943).  Morsi model may be avairable for incomplessible flow (https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/an-investigation-of-particle-trajectories-in-twophase-flow-systems/5B1FF01A248EBF4988C156246EFF844A) and Stokes-Millikan is limitted to only low Reynolds number flow. Here is one example when Singh drag model is applied
~~~
dragModel Singh
~~~
4. Turbulent dispersion <br>
"Dispersion" mention the turbulent dispersion and syntax is Yes or No (default is No).  On the turbulent flow, the eddy repeat generation and dissipation which mean the flow becomes random.  One of the approach to treat such turbulent flow as a steady state is RANS model approach, e.g., k - $\epsilon$ model.  This code also able to treat the effect of this random eddy from the parameters used in RANS model that are turbulent kinetic energy, k and the dissipation rate of the eddy, $\epsilon$ by following this method (https://arc.aiaa.org/doi/10.2514/3.62687).  This simulaiton require k and $\epsilon$ files under CFD simulation result directory.  Here is one example when it is ON
~~~
Dispersion  Yes
~~~
5. Compressibility <br>
"compressible" mention the compressibility of the fluid and it is Yes or No (default is No).  When it is Yes, temperature and density files (T and rho) are needed to be located in the CFD simulation result directory.  Here is one example when it is compressible flow
~~~
compressible  Yes
~~~
6. Froude Krylov force <br>
Under prepearation
~~~
FroudeKrylov  Yes
~~~
7. Observation interval <br>
"observeTime" mentiond the interval of the output of the particle location and the velocity.  Following syntax is the interval in the unit of second (default value is 1e-6).  When you export it every 1 ms, it is written as
~~~
observeTime 1e-3
~~~
8. Start directory
"startDir" mention a directory which store CFD simulation results and syntax is the name of the directory (default is 20000). When the ending time of the CFD simulation is 20000, the it is written as
~~~
startDir  20000
~~~
9. Dimension
"dimension" set the dimension as 3D, 2D, and 2D axi-symetric that syntaxes are 3D, 2D, and 2Daxi, respectively.  Two more syntax is required for 2Daxi case to indicate the axis.  100, 010, and 001 respectively indicate the direction of the axis is x, y, and z cordinate.  Third syntax is position of the axis. This example is when the axis is y=0 case.
~~~
dimension 010 0
~~~
#### 3.3 particleSet file in particle directory
This file stores the initial particle location and each particle size.  Initial line need to be "x y z dp" and from next lines, you can give the particle initial positions and diameters.  When the first particle is start from $(x, y, z) = (0, 0, 0)$ and second is $(x, y, z) = (1, 1, 1)$ and the particle diamters are 100 nm and 200 nm, it is expressed as
~~~
x y z dp
0 0 0 1e-7
1 1 1 2e-7
...
~~~
#### 3.4 Run the simulation
When you run the simulation, put trajectory.out file in the running directory and run execute file on the terminal by typing
~~~
./trajectory
~~~
### 4 Postprocessing

## Author
* Dr. Tomoya Tamadate
* [LinkedIn](https://www.linkedin.com/in/tomoya-tamadate-953673142/)/[ResearchGate](https://www.researchgate.net/profile/Tomoya-Tamadate)/[Google Scholar](https://scholar.google.com/citations?user=XXSOgXwAAAAJ&hl=ja)
* University of Minnesota
* tamalab0109[at]gmail.com
