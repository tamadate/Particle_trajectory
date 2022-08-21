#pragma once
#include "functions.hpp"



// Morsi & Alexander (1972)
class dragForceMA : public dragForceSM{
	public:
		double computeCd(double Re, double Mach, double Cc);
	private:
};

// Singh et al. (2022)
class dragForceSingh : public dragForceSM{
	public:
		double computeCd(double Re, double Mach, double Cc);
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

// Loth model
class dragForceLoth : public dragForceSM{
	public:
		double computeCd(double Re, double Mach, double Cc);
	private:
};


void 
dragForceSM::computeFD(Variables *vars, Flags *flags, particle &par){
	double dUx=vars->U[par.cell].x[0]-par.v.x[0];
	double dUy=vars->U[par.cell].x[1]-par.v.x[1];
	double dUz=vars->U[par.cell].x[2]-par.v.x[2];
	double U2=dUx*dUx+dUy*dUy+dUz*dUz;
	double Umag=sqrt(U2);
	double cg=sqrt(gamkb_m*vars->T[par.cell]);

	par.Re=vars->rho[par.cell]*Umag*par.dp/vars->myu[par.cell]+1e-20;
	par.Mach=Umag/cg;
	double Cd=computeCd(par.Re, par.Mach, par.Cc);

	double FD=M_PI*0.125*vars->myu[par.cell]*par.dp*par.Re*Cd/par.m;

	par.F.x[0]=FD*dUx;
	par.F.x[1]=FD*dUy;
	par.F.x[2]=FD*dUz;
}

double 
dragForceSM::computeCd(double Re, double Mach, double Cc){
	return 24/Re/Cc;
}

double 
dragForceMA::computeCd(double Re, double Mach, double Cc){
	double CdRE;

	if (Re < 0.1){CdRE = 24.0 / Re;}
	else if (Re < 1.0){	CdRE = 22.73 / Re + 0.0903 / Re / Re + 3.69;}
	else if (Re < 10.0){CdRE = 29.1667 / Re - 3.8889 / Re / Re + 1.222;}
	else if (Re < 100.0){CdRE = 46.5 / Re - 116.67 / Re / Re + 0.6167;}
	else if (Re < 1000.0){CdRE = 98.33 / Re - 2778.0 / Re / Re + 0.3644;}
	else if (Re < 5000.0){CdRE = 148.62 / Re - 47500.0 / Re / Re + 0.357;}
	else if (Re < 10000.0){CdRE = -490.546 / Re + 578700.0 / Re / Re + 0.46;}
	else if (Re < 50000.0){CdRE = -1662.5 / Re - 5416700.0 / Re / Re + 0.5191;}
	else{CdRE = 0.4;}
	return CdRE / Cc;
}

double 
dragForceSingh::computeCd(double Re, double Mach, double Cc){
	double Ts_T, Us_U, Ms, alpha, C1, Br, fKnWr;
	double Cdfm=0;
    if(Mach<1){
		Ts_T=1;
		Us_U=1;
		Ms=Mach;
		alpha=1;
		C1=1;
		Br=0;
		fKnWr=1/Cc;
	}
	else{
		double MachSQ=Mach*Mach;
		double s_C1=(gam_m*MachSQ+2);
		double s_C2=(2*gam*MachSQ-gam_m);
		Ts_T=s_C1*s_C2/(gam_pSQ*MachSQ);
		Us_U=s_C1/(gam_p*MachSQ);
		Ms=sqrt(s_C1/s_C2);
		double alpha_inv=alpha0*Mach+1-alpha0;
		alpha=1/alpha_inv;
		C1=C1C1/(1-C1C2/Mach);
	}

	double Wr=pow(Mach,2*omega)/Re;
	double WrT=Wr*pow(1+1/Ts_T,omega);
	double Mpow_TOMO=pow(Mach,2*omega-1);
	Br=WrT*(Mpow_TOMO+1)/Mpow_TOMO;

	double fKnWr_inv=Cc*(1+alpha_hoc*WrT);
	fKnWr=1/fKnWr_inv;

	double s=Mach*sqrt(gam*0.5);
	double ss=s*s;
	double ssss=ss*ss;
	Cdfm+=(1+2*ss)*exp(-ss)/(ss*s*sqrt(M_PI));
	Cdfm+=(4*ssss+4*ss-1)*erf(s)/(2*ssss);
	Cdfm+=2/(3*s)*sqrt(M_PI);

	double Ths=pow(1+gam_m*Ms*Ms*0.5,gam/gam_m);
	double Re_=Re*pow(1/(alpha*alpha*Ts_T),omega)*pow(Ths,Re_C);
	double C2=1+delta0/sqrt(Re_);
	double Cdc=C1*(1-alpha*Us_U)+C0*Ths*C2*C2;
	double Brita=pow(Br,ita);

	return (Cdc*fKnWr+Cdfm*Brita)/(1+Brita);
}

double 
dragForceLoth::computeCd(double Re, double Mach, double Cc){
	double M4=Mach*Mach*Mach*Mach;
	double Cd;
	if(Re<=45){
		double s=Mach*sqrt(gam*0.5);
		double ss=s*s;
		double ssss=ss*ss;
		double Cdfmd=(1+2*ss)*exp(-ss)/(ss*s*sqrt(M_PI))+(4*ssss+4*ss-1)*erf(s)/(2*ssss);
		double Cdfm=Cdfmd+2/(3*s)*sqrt(M_PI);
		double CdfmRe_inv=1+(Cdfmd/1.63-1)*sqrt(Re/45.0);
		double CdfmRe=Cdfm/CdfmRe_inv;
		double CdKnRe=24/Re*(1+0.15*pow(Re,0.687))/Cc;
		Cd=(CdKnRe+M4*CdfmRe)/(1+M4);
	}
	if(Re>45){
		double CM, GM;
		if(Mach>1.45) CM=2.044+0.2*exp(-1.8*log(Mach/1.5)*log(Mach/1.5));
		if(Mach<=1.45) CM=5/3.0+2/3.0*tanh(3*log(Mach+0.1));
		if(Mach>0.89) GM=0.0002+0.0008*tanh(12.77*(Mach-2.02));
		if(Mach<=0.89) GM=1-1.525*M4;
		double HM=1-0.258*CM/(1+514*GM);
		Cd=24/Re*(1+0.15*pow(Re,0.687))*HM;
		Cd+=0.42*CM/(1+42500*GM/pow(Re,1.16));
	}
	return Cd;
}


