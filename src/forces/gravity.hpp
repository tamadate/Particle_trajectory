#pragma once
#include "force.hpp"



// base function of force
class gravity : public force{
	public:
		Variables *vars;
		Flags *flags;
		double a[3];
		virtual void compute(particle &par){
			par.F.x[0]+=a[0];
			par.F.x[1]+=a[1];
			par.F.x[2]=+a[2];

		};
		gravity(double x, double y, double z){
			a[0]=x;
			a[1]=y;
			a[2]=z;
		}
		~gravity(void){}
	private:
};
