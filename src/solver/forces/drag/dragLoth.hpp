#pragma once
#include "../force.hpp"



// Loth model
class dragForceLoth : public force{
	public:
		void compute(particle &par){
			double Cd=computeCd(par);	// Here, this code use another equations for each drag models

			double FD=M_PI*0.125*vars->myu[par.cell]*par.dp*par.Re*Cd/par.m;

			double dUx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
			double dUy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
			double dUz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
			par.F.x[0]+=FD*dUx;
			par.F.x[1]+=FD*dUy;
			par.F.x[2]+=FD*dUz;
		};
		double computeCd(particle &par){
			double M2=par.Mach*par.Mach;
			double M4=M2*M2;
			double Cd;
			if(par.Re<=45){
				double s=par.Mach*sqrt(gam*0.5);
				double ss=s*s;
				double ssss=ss*ss;
				double Cdfmd=(1+2*ss)*exp(-ss)/(ss*s*sqrt(M_PI))+(4*ssss+4*ss-1)*erf(s)/(2*ssss);
				double Cdfm=Cdfmd+2/(3*s)*sqrt(M_PI);
				double CdfmRe_inv=1+(Cdfmd/1.63-1)*sqrt(par.Re/45.0);
				double CdfmRe=Cdfm/CdfmRe_inv;
				double CdKnRe=24/par.Re*(1+0.15*pow(par.Re,0.687))/par.Cc;
				Cd=(CdKnRe+M4*CdfmRe)/(1+M4);
			}
			if(par.Re>45){
				double CM, GM;
				if(par.Mach>1.45) CM=2.044+0.2*exp(-1.8*log(par.Mach/1.5)*log(par.Mach/1.5));
				if(par.Mach<=1.45) CM=5/3.0+2/3.0*tanh(3*log(par.Mach+0.1));
				if(par.Mach>0.89) GM=0.0002+0.0008*tanh(12.77*(par.Mach-2.02));
				if(par.Mach<=0.89) GM=1-1.525*M4;
				double HM=1-0.258*CM/(1+514*GM);
				Cd=24/par.Re*(1+0.15*pow(par.Re,0.687))*HM;
				Cd+=0.42*CM/(1+42500*GM/pow(par.Re,1.16));
			}
			return Cd;
		};
		dragForceLoth(void){};
		~dragForceLoth(void){};
	private:
};
