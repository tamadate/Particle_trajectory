#pragma once
#include "functions.hpp"

class Flags {
	private:

  public:
		int Nth;
		int dragFlag;					// Type of drag force function
		int KFFlag;						// Froude-Krylov force?
		int compressFlag;			// Compressible flow?
		int dispersionFlag;		// Turbulent dispersion?
		int dimensionFlag;		// 2D or 3D?
		int initialBreakFlag; // break flag at initial (in initial statement fault case)
		int autoStep;					// Time step is auto? or fix?
		int analytical;				// Use analytical solution?
		int inletFace;				// Particle inlet face number (this is not avairable now)


		Flags(void){
			int dragFlag=1;
			int autoStep=0;
			int KFFlag=0;
			int dispersionFlag=0;
			int compressFlag=0;
			int dimensionFlag=0;
			int initialBreakFlag=0;
			int analytical=0;
			int inletFace=-1;
		};
		~Flags(void){};
};
