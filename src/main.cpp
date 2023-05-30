#include "trajectory.hpp"




int main(int argc,char *argv[]){
	cout<<"**************************************"<<endl;
	cout<<"*** Particle Trajectory Simulation ***"<<endl;
	cout<<"*** Last editing 2023/05/29        ***"<<endl;
	cout<<"*** University of Minnesota        ***"<<endl;
	cout<<"*** Kanazawa University            ***"<<endl;
	cout<<"**************************************"<<endl;
	
	// trajectory class is initialized with constractor (see trajectory.cpp)
	trajectory *tra = new trajectory( );
	// This is main part of this simulaiton
	tra->run();  
	return 0;
}
