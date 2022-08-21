#include "read.hpp"
#include "trajectory.hpp"


int main(int argc,char *argv[]){

	Variables *vars = new Variables();
	Flags *flags = new Flags();
	read *reader = new read();

	reader->readCondition(vars,flags);
	reader->readGeometry();
	reader->readCFDresults(vars,flags);


	calculateMyu(vars);
	calculateRamda(vars);

	reader->readParticles(vars);
	makeCells();
	initialParticle(vars,flags);
	outputInitial(vars);
	if (flags->initialBreakFlag==0) trajectory(vars,flags);
	


	return 0;
}





