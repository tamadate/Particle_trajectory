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

#include "structure.hpp"
#include "constants.hpp"


// Calculation conditions (flag, file path, etc...)
double dt=1e-7;
double rho_p=1000;
double observeTime=1e-5;
double totalTime=1;
int Observe=10000;	// output time steps
double delta_r2_trajectory=1e-4*1e-4;	// output migration distance for the trajectory
int plane2D=-1;
string startDir="20000";
double Axis=0;
double flag=1;
char filepath[100];
int boundaryStartID;
double constTemp=300;
double constRho=1.2;
double timer;
FILE*f;
FILE*file;

class Variables {
	private:
    public:
		std::vector<double> p;
		std::vector<point> U;
		std::vector<double> T;
		std::vector<double> rho;
		std::vector<double> myu;
		std::vector<double> ramda;
		std::vector<point> dp;
		std::vector<particle> particles;
		double delta_r2;
		double time;
		double preOutTime;
		double analyticalStep;

		Variables(void){time=0;preOutTime=0;};
		~Variables(void){};
};


class Flags {
	private:
    public:
		int dragFlag;
		int KFFlag;
		int compressFlag;
		int dimensionFlag;
		int initialBreakFlag;
		int autoStep;
		int breakFlag;
		int analytical;

		Flags(void){
			int dragFlag=1;
			int autoStep=0;
			int KFFlag=0;
			int compressFlag=0;
			int dimensionFlag=0;
			int initialBreakFlag=0;
			int analytical=0;
		};
		~Flags(void){};

};


//Stokes-Millikan
class dragForceSM{
	public:
		void computeFD(Variables *vars, Flags *flags, particle &par);
		virtual double computeCd(double Re, double Mach, double Cc);
	private:
};
std::vector<dragForceSM*> forces;




void trajectory(Variables *vars,  Flags *flags);
int checkCell(Variables *vars, int pid);
int boundAction(Variables *vars, int faceID, int particleID, point norm);
void updateDisp(Variables *vars, Flags *flags, particle &par);
void initialize(Variables *vars, Flags *flags){
	vars->time=0;
	vars->preOutTime=0;
	flags->breakFlag=0;
}


std::vector<face> faces;
std::vector<point> points;
std::vector<boundary> boundaries;
std::vector<cell> cells;
std::vector<int> neighbors;
std::vector<int> owners;
std::vector<outParticle> outParticles;
std::vector<int> erasePID;

void calculateMyu(Variables *vars){
	string str;
	int iflag=0;

	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double myu=myu0*pow(vars->T[i]/ST0,1.5)*(ST0+SC0)/(vars->T[i]+SC0);
		vars->myu.push_back(myu);
	}
}

void calculateRamda(Variables *vars){
	string str;
	int iflag=0;

	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double ramda=ramda_coeff*vars->myu[i]/vars->rho[i]/sqrt(vars->T[i]);
		vars->ramda.push_back(ramda);
	}
}

void  computeReMach(Variables *vars, particle &par){
	double dUx=vars->U[par.cell].x[0]-par.v.x[0];
	double dUy=vars->U[par.cell].x[1]-par.v.x[1];
	double dUz=vars->U[par.cell].x[2]-par.v.x[2];
	double U2=dUx*dUx+dUy*dUy+dUz*dUz;
	double Umag=sqrt(U2);
	double cg=sqrt(gamkb_m*vars->T[par.cell]);

	par.Re=vars->rho[par.cell]*Umag*par.dp/vars->myu[par.cell]+1e-20;
	par.Mach=Umag/cg;
}

void makeCells(){
	timer=clock();
	int maxOwner=0;
	for (auto &a:owners) {if(a>maxOwner) maxOwner=a;}
	cell c;
	for (int i=0; i<maxOwner+1; i++) {cells.push_back(c);}
	for (int i=0; i<owners.size(); i++) {
		cells[owners[i]].iface.push_back(i);
	}
	for (int i=0; i<neighbors.size(); i++) {
		cells[neighbors[i]].iface.push_back(i);
	}

	for (auto &a:cells) {
		int faceSize=a.iface.size();	// how many faces in the cell
		int sumN=0;
		double X=0;
		double Y=0;
		double Z=0;
		for(auto &b: a.iface){
			for(auto &c: faces[b].iface){
				X+=points[c].x[0];
				Y+=points[c].x[1];
				Z+=points[c].x[2];
				sumN+=1;
			}
		}
		X/=double(sumN);
		Y/=double(sumN);
		Z/=double(sumN);
		a.r.x[0]=X;
		a.r.x[1]=Y;
		a.r.x[2]=Z;

		for (int i=0; i<faceSize; i++){
			int b=a.iface[i];
			int c0=faces[b].iface[0];
			int c1=faces[b].iface[1];
			int c2=faces[b].iface[2];

			double x0=points[c0].x[0];
			double y0=points[c0].x[1];
			double z0=points[c0].x[2];

			double x1=points[c1].x[0];
			double y1=points[c1].x[1];
			double z1=points[c1].x[2];

			double x2=points[c2].x[0];
			double y2=points[c2].x[1];
			double z2=points[c2].x[2];

			double Ax=x1-x0;
			double Ay=y1-y0;
			double Az=z1-z0;

			double Bx=x2-x0;
			double By=y2-y0;
			double Bz=z2-z0;

			point norm;
			norm.x[0]=Ay*Bz-Az*By;
			norm.x[1]=Az*Bx-Ax*Bz;
			norm.x[2]=Ax*By-Ay*Bx;
			
			double x3=X-x0;
			double y3=Y-y0;
			double z3=Z-z0;

			double naiseki=norm.x[0]*x3+norm.x[1]*y3+norm.x[2]*z3;
			if(naiseki<0){
				norm.x[0]=By*Az-Bz*Ay;
				norm.x[1]=Bz*Ax-Bx*Az;
				norm.x[2]=Bx*Ay-By*Ax;
			}
			naiseki=norm.x[0]*x3+norm.x[1]*y3+norm.x[2]*z3;
			if(naiseki<0){cout<<"error1011"<<endl;}
			double r2=norm.x[0]*norm.x[0]+norm.x[1]*norm.x[1]+norm.x[2]*norm.x[2];
			double r=sqrt(r2);
			norm.x[0]/=r;
			norm.x[1]/=r;
			norm.x[2]/=r;

			a.norm.push_back(norm);

				
		}

	}
	cout<<"Making cell time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
}


void findParticle(Variables *vars){
	timer=clock();
	int ps=vars->particles.size();
	for(int i=0; i<ps; i++){
		int foundCells=0;
		double x=vars->particles[i].x.x[0];
		double y=vars->particles[i].x.x[1];
		double z=vars->particles[i].x.x[2];
		int cellSize=cells.size();
        double minDist=1e15;
		for(int j=0; j<cellSize; j++) {
			cell a=cells[j];
			int faceSize=a.iface.size();
			int flag=0;
			/*int btest=a.iface[0];
			int c0test=faces[btest].iface[0];
			double y0test=points[c0test].x[1];
			if(y0test<0.05) continue;*/
			for (int k=0; k<faceSize; k++){
				int b=a.iface[k];
				int c0=faces[b].iface[0];
				double x0=points[c0].x[0];
				double y0=points[c0].x[1];
				double z0=points[c0].x[2];
				point c=a.norm[k];
				double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
				if(naiseki<0) {flag=1;break;}

			}
			if(flag==0){
				double Dist=0;
				for (int k=0; k<faceSize; k++){
					int b=a.iface[k];
					int c0=faces[b].iface[0];
					double x0=points[c0].x[0];
					double y0=points[c0].x[1];
					double z0=points[c0].x[2];
					point c=a.norm[k];
					double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
                    Dist+=naiseki;
				}			
				foundCells++;
				if (Dist<minDist){
					minDist=Dist;
					vars->particles[i].cell=j;		
				}
			}
		}
//		if (foundCells>0){cout<<"I found "<<foundCells<<" initial cell(s) for pid "<<i<<endl;}
		if (foundCells==0){
			cout<<"**Error: I could not find an initial cell for pid "<<i<<endl;
			for(int j=0; j<cellSize; j++) {
				cell a=cells[j];
				int faceSize=a.iface.size();
				int flag=0;
				for (int k=0; k<faceSize; k++){
					int b=a.iface[k];
					int c0=faces[b].iface[0];
					double x0=points[c0].x[0];
					double y0=points[c0].x[1];
					double z0=points[c0].x[2];
					point c=a.norm[k];
					double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
					if(naiseki<0) {flag=1;break;}

				}
				double Dist=0;
				for (int k=0; k<faceSize; k++){
					int b=a.iface[k];
					int c0=faces[b].iface[0];
					double x0=points[c0].x[0];
					double y0=points[c0].x[1];
					double z0=points[c0].x[2];
					point c=a.norm[k];
					double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
	                Dist+=naiseki;
				}			
				if (Dist<minDist){
					minDist=Dist;
					vars->particles[i].cell=j;		
				}
			}	
		}
	}
	cout<<"Find particle time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
}



void initialParticle(Variables *vars, Flags *flags){
	int ps=vars->particles.size();
	findParticle(vars);
	for(auto &a:vars->particles){
		int icell=a.cell;
		double Kn=vars->ramda[icell]/a.dp;
		double Cc=1+Kn*(A1+A2*exp(-A3/Kn));
		a.Kn=Kn;
		a.Cc=Cc;
		double threePiMuDp_Cc=3*M_PI*vars->myu[icell]*a.dp/Cc;
		a.Zp=1.6e-19/threePiMuDp_Cc;
		a.beta=threePiMuDp_Cc/a.m;
		a.dt=1/(40*a.beta);
		for(int i=0; i<3; i++) {
			a.v.x[i]=vars->U[icell].x[i]*0.99;
		}
		for(auto &force : forces) force->computeFD(vars,flags,a);
	}
}



	

