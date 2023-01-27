#include "trajectory.hpp"

// Reading funcitons for OpenFOAM style simulation results

void
trajectory::readScalar(char *readFile, std::vector<double> &variable){
	string str;
	int iflag=0;

	ifstream stream(readFile);
	while(getline(stream,str)) {
		if (str=="(") {iflag=1; continue;}
		if (str==")") {iflag=0; continue;}
		if (iflag==1){
			variable.push_back(stod(str));
		}
	}
	stream.close();
}

void
trajectory::readScalarDum(char *readFile, std::vector<double> &variable, double value){
	string str;
	ifstream stream(readFile);
	while(getline(stream,str)) variable.push_back(value);
	stream.close();
}

void
trajectory::readVector(char *readFile, std::vector<point> &variable){
	string str;
	int iflag=0;

	ifstream stream(readFile);
	while(getline(stream,str)) {
		if (str=="(") {iflag=1; continue;}
		if (str==")") {iflag=0; continue;}
		if (iflag==1){
			str=str.substr(str.find("(")+1);
			str=str.erase(str.find(")"));
			string tmp;
			istringstream stream(str);
			int loop=0;
			point POINT;
			while(getline(stream,tmp,' ')) {
				POINT.x[loop]=stod(tmp);
				loop++;
			}
			variable.push_back(POINT);
		}
	}
	stream.close();
}

void
trajectory::readVectorDum(char *readFile, std::vector<point> &variable, double value){
	string str;
	ifstream stream(readFile);
	point POINT;
	POINT.x[0]=value;
	POINT.x[1]=value;
	POINT.x[2]=value;
	while(getline(stream,str)) variable.push_back(POINT);
	stream.close();
}
