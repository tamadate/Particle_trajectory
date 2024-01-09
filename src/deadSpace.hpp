#pragma once
#include "functions.hpp"
#include "variables.hpp"
#include "flags.hpp"


class deadSpace{
	public:
		double xmin;
		double xmax;
		double ymin;
		double ymax;
		double zmin;
		double zmax;
		deadSpace(double x0, double x1, double y0, double y1, double z0, double z1){
			xmin=x0;
			xmax=x1;
			ymin=y0;
			xmax=y1;
			zmin=z0;
			zmax=z1;
		}
		~deadSpace(void){}
	private:
};
