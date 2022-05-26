#include "functions.hpp"
#include "read.hpp"
#include "trajectory.hpp"
#include "dragForce.hpp"



int main(int argc,char *argv[]){
	
	Variables *vars;
	vars = new Variables();

	readCondition(vars);

	readNeighbors();
	readFaces();
	readPoints();
	readOwners();
	readBoundaries();

	sprintf ( filepath, "%s/p", startDir.c_str());
	readScalar(filepath,vars->p);
	sprintf ( filepath, "%s/U", startDir.c_str());
	readVector(filepath,vars->U);

    /*  Dispersion? */
    if(dispFlag==1){
	    sprintf ( filepath, "%s/k", startDir.c_str());
	    readScalar(filepath,vars->k);
	    sprintf ( filepath, "%s/omega", startDir.c_str());
	    readScalar(filepath,vars->omega);
    }
    if(dispFlag==0){
	    sprintf ( filepath, "%s/p", startDir.c_str());
	    readScalarDum(filepath,vars->k, 1);
	    readScalarDum(filepath,vars->omega, 1);
    }

    /*  Compressible? */
    if(compressFlag==1){
	    sprintf ( filepath, "%s/T", startDir.c_str());
	    readScalar(filepath,vars->T);
	    sprintf ( filepath, "%s/rho", startDir.c_str());
	    readScalar(filepath,vars->rho);
    }
    if(compressFlag==0){
	    sprintf (filepath, "%s/p", startDir.c_str());
	    readScalarDum(filepath,vars->T, constTemp);
	    readScalarDum(filepath,vars->rho, constRho);
    }

    /*  Dispersion? */
    if(KFFlag==1){
	    sprintf ( filepath, "%s/dp", startDir.c_str());
	    readVector(filepath,vars->dp);
    }
    if(KFFlag==0){
	    sprintf ( filepath, "%s/p", startDir.c_str());
	    readVectorDum(filepath,vars->dp, 0);
    }

	calculateMyu(vars);
	calculateRamda(vars);

	readParticles(vars);
	makeCells();
	initialParticle(vars);
	outputInitial(vars);
	if (initialBreakFlag==0) trajectory(vars);
	


	return 0;
}





