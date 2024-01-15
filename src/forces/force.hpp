#pragma once
#include "../functions.hpp"
#include "../variables.hpp"
#include "../openMP.hpp"
#include "../flags.hpp"



// base function of force
class force{
	public:
		Variables *vars;
		Flags *flags;
		VariablesMP *varsMP;
		virtual void compute(particle &par){};
		void initial(Variables *Vars, Flags *Flags, VariablesMP *VarsMP){
			vars=Vars;
			flags=Flags;
			varsMP=VarsMP;
		};
		force(void){}
		~force(void){}
	private:
};
