#include "functions.hpp"
#include "read.hpp"
#include "trajectory.hpp"
#include "dragForce.hpp"



int main(int argc,char *argv[]){
	f=fopen("dragCoeff.dat", "w");
	for(double iRe=-4; iRe<4;iRe+=0.01){
		double Re=pow(10,iRe);
		for(double Mach=0.32;Mach<5;Mach+=100){
			double Kn=Mach/Re*sqrt(0.5*M_PI*gam);
			double Cc=1+Kn*(A1+A2*exp(-A3/Kn));
			double Cd=computeCd_Loth(Re,Mach,Cc);
			fprintf(f,"%e\t%e\t%e\t%e\n", Re, Mach, Kn, Cd);
		}
	}
	fclose(f);
	return 0;
}





