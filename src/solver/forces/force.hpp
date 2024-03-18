#pragma once
#include "../../functions.hpp"
#include "../../variables.hpp"
#include "../../flags.hpp"



// base function of force
class force{
	public:
		Variables *vars;
		Flags *flags;
		virtual void compute(particle &par){};
		void initial(Variables *Vars){
			vars=Vars;
		};
		force(void){}
		~force(void){}
	private:
};
