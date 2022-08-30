#pragma once
#include "functions.hpp"

class Flags {
	private:
    public:
		int dragFlag;
		int KFFlag;
		int compressFlag;
		int dispersionFlag;
		int dimensionFlag;
		int initialBreakFlag;
		int autoStep;
		int breakFlag;
		int analytical;
		int inletFace;

		Flags(void){
			int dragFlag=1;
			int autoStep=0;
			int KFFlag=0;
			int dispersionFla=0;
			int compressFlag=0;
			int dimensionFlag=0;
			int initialBreakFlag=0;
			int analytical=0;
			int inletFace=-1;
		};
		~Flags(void){};
};
