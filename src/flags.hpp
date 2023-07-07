#pragma once
#include "functions.hpp"

class Flags {
	private:

  public:
		int dispersionFlag;		// Turbulent dispersion?
		int dimensionFlag;		// 2D or 3D?
		bool autoStep;					// Time step is auto? or fix?
		int analytical;				// Use analytical solution?
		int inletFace;				// Particle inlet face number (this is not avairable now)
		bool v0;		// if v0=0, particle initial velocity is 99% of fluid velocity
					// if v0=1, particle initial velocity is given in the velocity setting file


		Flags(void){
			autoStep=true;
			dispersionFlag=0;
			dimensionFlag=0;
			analytical=0;
			inletFace=-1;
			v0=0;
		};
		~Flags(void){};
};
