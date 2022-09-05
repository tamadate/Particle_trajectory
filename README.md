# Particle_trajectory
This code calculate particle trajectory coupled with the CFD simulation result (only 1 way coupling is avairable in current version) from OpenFOAM (https://www.openfoam.com/).  Fluent results are also avairable by transforming to OpenFOAM fomat via "fluentToFoam" function in OpenFOAM.
## Usage
### 1. Building source code
Download this code and prepare compilation environment (c++ compiler, c++11 standard library and OpenMP are required). On the turminal, move to Particle_trajectory/src directory from your download directory and make code by just typing
~~~
make
~~~
After waiting for several seconds, a execute file trajectory.out is created in the src directory.  
### 2. Run (steady state) CFD simulation
This code is able to calculate particle trajectory coupled with steady state CFD simulation results, specifically OpenFOAM.  In terms of Fluent format, the mesh geometry is needed to be transformed to the OpenFOAM mesh format via fluentToFoam command from the OpenFOAM, which mean OpenFOAM installation is required anyways. The Fluent simulation result also need to be trasformed to OpenFOAM format using python script "datToFoam.py" located under the src directory.  The Fluent result format is also limitted and this code does not insure that your format is transformable.  In case of OpenFOAM format, the format changing is not required. You can find OpenFOAM installation guide from web search and the learn how to run it from manual and abandant tutorials.  This code developper recommend to do the CFD simulation with OpenFOAM due to the prospect of understanding of simulation as well as the affinity with this simulator.
### 3. Set conditions
This simulator require additional directories and files to store the simulation results and to set the calculation condition including particle initial conditions as following explanation.
#### 3.1 Creat directories
Under the CFD simulation directory, you must have four directories which are named 1.constant, 2.20000 (this is not needed to be 20000 but default setting is 20000 in this simulator), 3.result, and 4.particle. First directory is from early performed CFD simulation and mesh geometory information is stored in constant/polyMesh/ directory.  Second directory is also from CFD simulation and the name of directory is typically time of simulation if it is OpenFOAM, e.g., there is a directory named 20000, when the ending time of the simulation is 20000.  The CFD result directory name can be changed from the calculation condition file as written later.  A system file may be also here for setting CFD simulation condition if you run CFD with OpenFOAM but it does not matter if you leave or keep it since the trajectory simulation does not never access system directory.  On the third directory named "result" is a blank directory that is required to store the trajectory simulation results.  Please make sure to creat this directory; this has been often forgotten to creat and it dump error "segmentation fault". On the fourth directory, the name is needed to be "particle" where should have trajectory simulation conditions file (/particle/conditions) and particle information file (/particle/particleSet).  Detail of these two files are described following sections 3.2 and 3.3, respectively.
#### 3.2 conditions file in particle directory
This code read a file "conditions" to set the trajectory simulaiton conditions.  This file generally set a simulation condition in a line, e.g.,
~~~
dragModel	Stokes-Millikan
~~~
This example set the particle drag coefficient as Stokes-Millikan model.  As this example, first syntax mention the setting parameter and then setting values or method are written as 2nd (or more) syntax(es). The separation of syntaxes are only tab is avairable. All setting parameters and methods are shown below.
1. Time step <br>
"dt" mention the time step and following syntax is auto or value.  Default is auto.  When auto is mentioned, this code automatically determine the time step from the particle veloctiy and relaxation time (or duration of the eddy when the turbulent dispersion is ON).  When value is mentioned, the time step is fixed as that value (unit is second).  Here is one example when auto time step is utilized:
~~~
dt  auto
~~~
2. Total calculation time <br>
"totalTime" mention the total simulation time (wall simulation time) and following syntax is value in the unit of second.  Ideally, this time is infinity since the particle exit from the calculation domain when the time is long enough.  However, the particles are potentially trapped by the curculation and it make simulation time longer.  This wall time is set to prevent the situation that the simulation is never done. Here is one example when the wall time is 100 s:
~~~
totalTime  100
~~~
3. Drag model <br>
"dragModel" mention the drag model that is choosable from four drag models (Singh, Stokes-Millikan, Morsi, and Loth) and the syntax is same as the model name (default is Singh).  One of the feature of this code is Singh's drag model ([Singh et al., 2021](https://arc.aiaa.org/doi/10.2514/1.J060648)) and this is available for wide range of Reynodls number and Mach number even if it is supersonic flow, hence it is called general drag model as used in the title of original paper. Loth model is also available for the wide range ([Loth, 2008](https://arc.aiaa.org/doi/10.2514/1.28943)).  Morsi model may be avairable for incomplessible flow ([Morsi and Alexander, 2006](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/an-investigation-of-particle-trajectories-in-twophase-flow-systems/5B1FF01A248EBF4988C156246EFF844A)) and Stokes-Millikan is more limitted to only low Reynolds number flow. Here is one example when Singh drag model is applied:
~~~
dragModel Singh
~~~
4. Turbulent dispersion <br>
"Dispersion" mention the turbulent dispersion and syntax is Yes or No (default is No).  On the turbulent flow, the eddy repeat to generate and dissipate which mean the flow becomes random.  One of the approach to treat such turbulent flow as a steady state is RANS model approach in CFD, e.g., k - $\epsilon$ model.  This code also able to treat the effect of this random eddy from the parameters used in RANS model that are turbulent kinetic energy, k and the dissipation rate of the eddy, $\epsilon$ by following this method ([Gosmas and Loannides]https://arc.aiaa.org/doi/10.2514/3.62687).  This simulaiton require k and $\epsilon$ prifiles under CFD simulation result directory.  Here is one example when it is ON
~~~
Dispersion  Yes
~~~
5. Compressibility <br>
"compressible" mention the compressibility of the fluid and it is Yes or No (default is No).  When it is Yes, temperature and density profiles (T and rho) are needed to be located in the CFD simulation result directory, hence the compressible CFD simulation must to be performed.  Here is one example when it is compressible flow:
~~~
compressible  Yes
~~~
6. Froude Krylov force <br>
Under prepearation
~~~
FroudeKrylov  Yes
~~~
7. Observation interval <br>
"observeTime" mentiond the interval of the intermediate trajectory simulation results that are the particle location and the velocity.  Following syntax is the interval in the unit of second (default value is 1e-6).  When you export particle location and velocity every 1 ms, it becomes:
~~~
observeTime 1e-3
~~~
8. Start directory
"startDir" mention a directory which store CFD simulation results and syntax is the name of the directory (default is 20000). When the CFD simulation result is stored a directory named 100, it is:
~~~
startDir  100
~~~
9. Dimension
"dimension" set the dimension as 3D, 2D, and 2D axi-symetric that syntaxes are 3D, 2D, and 2Daxi, respectively.  Two more syntax is required for 2Daxi case to indicate the axis.  100, 010, and 001 respectively indicate the direction of the axis is x, y, and z cordinate.  Third syntax is position of the axis. This example is when the axis is y=0 case:
~~~
dimension 010 0
~~~
#### 3.3 particleSet file in particle directory
This file stores the initial particle location and each particle size.  Initial line need to be "x y z dp" (again, the separation is tab) and from next lines, you can give the particle initial positions and diameters.  When the first particle is start from $(x, y, z) = (0, 0, 0)$ and second is $(x, y, z) = (1, 1, 1)$ and the particle diameters are 100 nm and 200 nm, it is expressed as
~~~
x y z dp
0 0 0 1e-7
1 1 1 2e-7
...
~~~
### 4 Run the simulation
When you run the simulation, put trajectory.out file in the running directory and run execute file on the terminal by typing
~~~
./trajectory
~~~
### 5 Postprocessing

## Author
* Dr. Tomoya Tamadate
* [LinkedIn](https://www.linkedin.com/in/tomoya-tamadate-953673142/)/[ResearchGate](https://www.researchgate.net/profile/Tomoya-Tamadate)/[Google Scholar](https://scholar.google.com/citations?user=XXSOgXwAAAAJ&hl=ja)
* University of Minnesota
* tamalab0109[at]gmail.com
