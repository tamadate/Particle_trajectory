#include "trajectory.hpp"




int main(int argc,char *argv[]){
	// trajectory class is initialized with constractor (see trajectory.cpp)
	trajectory *tra = new trajectory( );		
	// This is main part of this simulaiton
	tra->run();  
	return 0;
}
