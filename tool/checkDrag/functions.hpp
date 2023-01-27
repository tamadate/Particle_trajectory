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
const double myu0=1.81e-5;
const double rho_p=1000;
const double e=1.6e-19;
const double lamda0=67.0e-9;
const double mgas=(28*0.8+32*0.2)*0.001/6.02e23;
const double kb=1.38e-23;

// Calculation conditions (flag, file path, etc...)
double dt=1e-7;
double timestep=1e8;
int Observe=10000;
int dispFlag=0;
int dragFlag=1;
int KFFlag=0;
int compressFlag=0;
int dimensionFlag=0;
int plane2D=-1;
int startTime=20000;
double Axis=0;
int symmerticFlag=1;
double flag=1;
char filepath[100];
int boundaryStartID;
int outputLoop=0;
double constTemp=300;
double constRho=1.2;
FILE*f;
FILE*file;



// Parameters for slip correction
const double A1=2.514;
const double A2=0.8;
const double A3=0.55;

//	Parameters for drag force
const double delta0=9.4; // or 9.06?
const double alpha0=0.356;
const double ita=1.8;
const double alpha_hoc=1.27;
const double C0=24/delta0/delta0; // C0*delta0*delta0=24
const double CdMach=0.9;
const double Cp=1.00;
const double Cv=0.718;
const double gam=Cp/Cv;
const double omega=0.74;
const double gam_m=gam-1;
const double gam_mSQ=gam_m*gam_m;
const double gam_p=gam+1;
const double gam_pSQ=gam_p*gam_p;
const double gamkb_m=gam*kb/mgas;

//  For high speed calculation...
const double lamda_coeff=sqrt(M_PI*mgas*0.5/kb);
const double C1C1=CdMach-C0*pow(1+gam_m*gam_m*0.25/gam,gam/gam_m);
const double C1C2=gam_m/alpha0/gam_p;
const double Re_C=gam_p*0.5/gam-gam_m/gam*omega;




struct face {
	std::vector<int> iface;
};

struct point {
	double x[3];
};
struct boundary {
	string name;
	string type;
	int nFaces;
	int startFace;
};
struct cell {
	std::vector<int> iface;
	std::vector<point> norm;
	point r;
};

struct outParticle {
	int pid;
    int bid;
	point r;
};


struct particle {
	point x;
	point v;
	point F;
	int cell;
	int id;
	double Re;
	double Mach;
	double Kn;
	double dp;
	double m;
	double Cc;
	double Zp;
	double beta;
	double le;
	int randUpdate;
	point Urand;
	int randLoop;
    int reflect;
};

class Variables {
	private:


    public:
		std::vector<double> p;
		std::vector<point> U;
		std::vector<double> T;
		std::vector<double> rho;
		std::vector<double> k;
		std::vector<double> omega;
		std::vector<double> myu;
		std::vector<double> lamda;
		std::vector<point> dp;
		std::vector<particle> particles;
		double time;
		Variables(void){time=0;};
		~Variables(void){};

};


void trajectory(Variables *vars);
void velocityUpdate(Variables *vars);
void positionUpdate(Variables *vars);
void forceUpdate(Variables *vars);
void checkCell(Variables *vars);
void symmetricModification(Variables *vars);
int boundAction(Variables *vars, int faceID, int particleID, point norm);
void updateDisp(Variables *vars, particle &par);
double get_vTD(double dispEnergy);
double get_tini(Variables *vars, particle par);

double computeFD(Variables *vars, particle &par);
double computeCd_AIAA(double Re, double Mach, double Cc);
double computeCd_Stokes(double Re, double Cc);
double computeCd_Re(double Re, double Cc);
double computeCd_Loth(double Re, double Mach, double Cc);
double computeCd(Variables *vars, particle &par);

std::vector<face> faces;
std::vector<point> points;
std::vector<boundary> boundaries;
std::vector<cell> cells;
std::vector<int> neighbors;
std::vector<int> owners;
std::vector<outParticle> outParticles;



void calculateMyu(Variables *vars){
	string str;
	int iflag=0;

	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double myu=myu0*pow(vars->T[i]/288.15,1.5)*(288.15+110.4)/(vars->T[i]+110.4);
		vars->myu.push_back(myu);
	}
}

void calculatelamda(Variables *vars){
	string str;
	int iflag=0;

	int fieldSize=vars->T.size();
	for(int i=0; i<fieldSize; i++){
		double lamda=lamda_coeff*vars->myu[i]/vars->rho[i]/sqrt(vars->T[i]);
		vars->lamda.push_back(lamda);
	}
}



void makeCells(){
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
}


void findParticle(Variables *vars){
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
		if (foundCells==0){cout<<"**Error: I could not find an initial cell for pid "<<i<<endl;breakFlag=1;}
	}
}



void initialParticle(Variables *vars){
	int ps=vars->particles.size();
	findParticle(vars);
	for(auto &a:vars->particles){
		int icell=a.cell;
		double Kn=vars->lamda[icell]/a.dp;
		double Cc=1+Kn*(A1+A2*exp(-A3/Kn));
		a.Kn=Kn;
		a.Cc=Cc;
		a.Zp=Cc*1.6e-19/(3*M_PI*vars->myu[icell]*a.dp);
		a.beta=3.0*M_PI*vars->myu[icell]*a.dp/Cc/a.m;
		a.le=pow(vars->k[icell],0.5)/vars->omega[icell]/pow(0.09,0.25);
		for(int i=0; i<3; i++) {
			a.v.x[i]=vars->U[icell].x[i]*0.99;
			a.F.x[i]=computeFD(vars, a)*(vars->U[icell].x[i]-a.v.x[i]);
		}
	}

}
