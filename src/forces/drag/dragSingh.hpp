#pragma once
#include "../force.hpp"


// Singh et al. (2021)
class dragForceSingh : public force{
	public:
		void compute(particle &par){
			double Cd=computeCd(par);	// Here, this code use another equations for each drag models

			double FD=M_PI*0.125*vars->myu[par.cell]*par.dp*par.Re*Cd/par.m;

			double dUx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
			double dUy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
			double dUz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
			par.F.x[0] += FD * dUx;
			par.F.x[1] += FD * dUy;
			par.F.x[2] += FD * dUz;

		};

		double computeCd(particle &par){
			double Ts_T, Us_U, Ms, alpha, C1, Br, fKnWr;
			double Cdfm=0;
			if(par.Mach<1){
				Ts_T=1;
				Us_U=1;
				Ms=par.Mach;
				alpha=1;
				C1=1;
				Br=0;
				fKnWr=1/par.Cc;
			}
			else{
				double MachSQ=par.Mach*par.Mach;
				double s_C1=(gam_m*MachSQ+2);
				double s_C2=(2*gam*MachSQ-gam_m);
				Ts_T=s_C1*s_C2/(gam_pSQ*MachSQ);
				Us_U=s_C1/(gam_p*MachSQ);
				Ms=sqrt(s_C1/s_C2);
				double alpha_inv=alpha0*par.Mach+1-alpha0;
				alpha=1/alpha_inv;
				C1=C1C1/(1-C1C2/par.Mach);
			}

			double Wr=pow(par.Mach,2*omega)/par.Re;
			double WrT=Wr*pow(1+1/Ts_T,omega);
			double Mpow_TOMO=pow(par.Mach,2*omega-1);
			Br=WrT*(Mpow_TOMO+1)/Mpow_TOMO;

			double fKnWr_inv=par.Cc*(1+alpha_hoc*WrT);
			fKnWr=1/fKnWr_inv;

			double s=par.Mach*sqrt(gam*0.5);
			double ss=s*s;
			double ssss=ss*ss;
			Cdfm+=(1+2*ss)*exp(-ss)/(ss*s*sqrt(M_PI));
			Cdfm+=(4*ssss+4*ss-1)*erf(s)/(2*ssss);
			Cdfm+=2/(3*s)*sqrt(M_PI);

			double Ths=pow(1+gam_m*Ms*Ms*0.5,gam/gam_m);
			double Re_=par.Re*pow(1/(alpha*alpha*Ts_T),omega)*pow(Ths,Re_C);
			double C2=1+delta0/sqrt(Re_);
			double Cdc=C1*(1-alpha*Us_U)+C0*Ths*C2*C2;
			double Brita=pow(Br,ita);

			return (Cdc*fKnWr+Cdfm*Brita)/(1+Brita);

		};
		dragForceSingh(void){};
		~dragForceSingh(void){};
	private:
		//	Parameters for drag force
		const double delta0=9.4; // or 9.06?
		const double alpha0=0.356;
		const double ita=1.8;
		const double alpha_hoc=1.27;
		const double C0=24/delta0/delta0; // C0*delta0*delta0=24
		const double CdMach=0.9;
		const double omega=0.74;
		const double gam_m=gam-1;
		const double gam_mSQ=gam_m*gam_m;
		const double gam_p=gam+1;
		const double gam_pSQ=gam_p*gam_p;

		const double C1C1=CdMach-C0*pow(1+gam_m*gam_m*0.25/gam,gam/gam_m);
		const double C1C2=gam_m/alpha0/gam_p;
		const double Re_C=gam_p*0.5/gam-gam_m/gam*omega;
};
