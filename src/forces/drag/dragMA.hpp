#pragma once
#include "../force.hpp"


// Morsi & Alexander (1972)
class dragForceMA : public force{
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
			double CdRE;
			double Re_inv=1/par.Re;
			double Re_inv2=Re_inv*Re_inv;
			if (par.Re < 0.1){CdRE = 24.0 * Re_inv;}
			else if (par.Re < 1.0){	CdRE = 22.73 * Re_inv + 0.0903 * Re_inv2 + 3.69;}
			else if (par.Re < 10.0){CdRE = 29.1667 * Re_inv - 3.8889 * Re_inv2 + 1.222;}
			else if (par.Re < 100.0){CdRE = 46.5 * Re_inv - 116.67 * Re_inv2 + 0.6167;}
			else if (par.Re < 1000.0){CdRE = 98.33 * Re_inv - 2778.0 * Re_inv2 + 0.3644;}
			else if (par.Re < 5000.0){CdRE = 148.62 * Re_inv - 47500.0 * Re_inv2 + 0.357;}
			else if (par.Re < 10000.0){CdRE = -490.546 * Re_inv + 578700.0 * Re_inv2 + 0.46;}
			else if (par.Re < 50000.0){CdRE = -1662.5 * Re_inv - 5416700.0 * Re_inv2 + 0.5191;}
			else{CdRE = 0.4;}
			return CdRE / par.Cc;
		};
		dragForceMA(void){};
		~dragForceMA(void){};
	private:
};