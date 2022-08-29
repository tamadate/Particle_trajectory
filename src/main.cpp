#include "trajectory.hpp"


int main(int argc,char *argv[]){

	trajectory *tra = new trajectory();

	if (tra->flags->initialBreakFlag==0) tra->run();

	return 0;
}
