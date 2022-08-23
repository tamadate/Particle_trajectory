#pragma onece
#include "functions.hpp"


void output(Variables *vars){
	for (auto &a:vars->particles){
		sprintf(filepath, "result/position.%d", int(a.id)); 
		f=fopen(filepath, "a");
		fprintf(f,"%e\t%e\t%e\t%d\n", a.x.x[0], a.x.x[1], a.x.x[2], a.cell);
		fclose(f);

		sprintf(filepath, "result/U.%d", int(a.id)); 
		f=fopen(filepath, "a");
		fprintf(f,"%e\t%e\t%e\n", a.v.x[0], a.v.x[1], a.v.x[2]);
		fclose(f);
	}
	sprintf(filepath, "result/Time"); 
	f=fopen(filepath, "a");
	fprintf(f,"%e\n", vars->time);
	fclose(f);

}


void outputInitial(Variables *vars){
	for (auto &a:vars->particles){
		sprintf(filepath, "result/position.%d", int(a.id)); 
		f=fopen(filepath, "w");
		fclose(f);

		sprintf(filepath, "result/U.%d", int(a.id)); 
		f=fopen(filepath, "w");
		fclose(f);
	}
	sprintf(filepath, "result/Time"); 
	f=fopen(filepath, "w");
	fclose(f);
}

