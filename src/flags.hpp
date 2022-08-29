#pragma once
#include "functions.hpp"

class Flags {
	private:
    public:
		int dragFlag;
		int KFFlag;
		int compressFlag;
		int dimensionFlag;
		int initialBreakFlag;
		int autoStep;
		int breakFlag;
		int analytical;

		Flags(void){
			int dragFlag=1;
			int autoStep=0;
			int KFFlag=0;
			int compressFlag=0;
			int dimensionFlag=0;
			int initialBreakFlag=0;
			int analytical=0;
		};
		~Flags(void){};
};
