#pragma once
#include "functions.hpp"

class Flags {
	private:

  public:
		int dimensionFlag;		// 2D or 3D?
		int inletFace;				// Particle inlet face number (this is not avairable now)
		bool v0;		// if v0=0, particle initial velocity is 99% of fluid velocity


		Flags(void){
			dimensionFlag=0;
			inletFace=-1;
			v0=0;
		};
		~Flags(void){};
};
