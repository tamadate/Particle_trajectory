#pragma once

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <stdlib.h>

using namespace std;

//	Constant
const double mu0=1.81e-5;
const double rho_p=1000;
const double e=1.6e-19;
const double ramda0=67.0e-9;
const double mgas=(28*0.8+32*0.2)*0.001/6.02e23;
const double kb=1.38e-23;
const double T=300;
const double rho=1.2;

// Calculation conditions (dt, total time step etc...)
const double dt=1e-8;
const double timestep=1e7;
const int Observe=10; //file output time step
char filepath[100];
FILE*f;

// Parameters for slip correction
const double A1=2.514;
const double A2=0.8;
const double A3=0.55;

struct point {
	double x[3];
};

struct particle {
	point x;	//position
	point v;	//velocity 
	point F;	//force
	int id;	//particle id
	int cell;	//cell number
	double Re;	//Reynolds number
	double Kn;	//Knudsen number
	double dp;	//diameter
	double m;	//mass [kg]
	double Cc;	//Cc
	double Zp;	//electric mobility
	double beta;	//relaxation time
};

class Variables {
	private:

    public:
		double p;	//field pressure
		point U;	//field velocity
		std::vector<particle> particles;	//particle array
		double time;	//time

		Variables(void){
			U.x[0]=1.0;
			U.x[1]=0;
			U.x[2]=0;
			p=1e5;
			time=0;
		};
		~Variables(void){};
};

void trajectoryEuler(Variables *vars);
void velocityUpdateEuler(Variables *vars);
void positionUpdateEuler(Variables *vars);
void forceUpdate(Variables *vars);
double computeFD(Variables *vars, particle &par);

void initialParticle(Variables *vars){
	for(auto &a:vars->particles){
		a.m=M_PI*a.dp*a.dp*a.dp/6.0*rho_p;
		a.Kn=ramda0/a.dp;
		a.Cc=1+a.Kn*(A1+A2*exp(-A3/a.Kn));
		a.Zp=a.Cc*e/(3*M_PI*mu0*a.dp);
		a.beta=3.0*M_PI*a.dp/a.Cc/a.m;

		for(int i=0; i<3; i++) {
			a.v.x[i]=1e-10;
			a.F.x[i]=computeFD(vars, a)*(vars->U.x[i]-a.v.x[i]);
		}
	}

}



	

