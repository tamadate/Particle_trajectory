#include "trajectory.hpp"

void
trajectory::readFaces(void){
	string str;
	string strPre;
	int iflag=0;

	ifstream stream("constant/polyMesh/faces");
	while(getline(stream,str)) {
		if (str=="(") {iflag=1; continue;}
		if (str==")") {iflag=0; continue;}
		if (iflag==1){
			str=str.substr(str.find("(")+1);
			str=str.erase(str.find(")"));
			string tmp;
			istringstream stream(str);
			int loop=0;
			face FACE;
			while(getline(stream,tmp,' ')) {
				FACE.iface.push_back(stoi(tmp));
			}
			faces.push_back(FACE);
		}
	}
	stream.close();
}

void
trajectory::readBoundaries(void){
	string str;
	string strPre;
	int iflag=0;
	int jflag=0;

	ifstream stream("constant/polyMesh/boundary");
	boundary bound;
	while(getline(stream,str)) {
		string tmp;
		istringstream stream(str);
		int outType=0;
		char del=' ';
		while(getline(stream,tmp,del)) {
			if (tmp.size()==0) {continue;}
			if (tmp=="(") {iflag=1; continue;}
			if (tmp==")") {iflag=0; continue;}
			if (tmp=="{"&&iflag==1) {jflag=1; bound.name=strPre; continue;}
			if (tmp=="}"&&iflag==1) {jflag=0; boundaries.push_back(bound);continue;}
			if (iflag==1&&jflag==1){
				if(outType==0&&tmp=="type") {outType=1;continue;}
				if(outType==0&&tmp=="nFaces") {outType=2;continue;}
				if(outType==0&&tmp=="startFace") {outType=3;continue;}
				if(outType==1) {bound.type=tmp.erase(tmp.find(";"));continue;}
				if(outType==2) {bound.nFaces=stoi(tmp);continue;}
				if(outType==3) {bound.startFace=stoi(tmp);continue;}
			}
		}
		strPre=tmp;
	}
	boundaryStartID=boundaries[0].startFace;
	for(auto &a:boundaries){
		if(a.startFace<boundaryStartID) boundaryStartID=a.startFace;
	}
	stream.close();
}

void
trajectory::readPoints(void){
	string str;
	string strPre;
	int iflag=0;

	ifstream stream("constant/polyMesh/points");
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
			points.push_back(POINT);
		}
	}
	stream.close();
}


void
trajectory::readNeighbors(void){
	string str;
	string strPre;
	int iflag=0;

	ifstream stream("constant/polyMesh/neighbour");
	while(getline(stream,str)) {
		if (str=="(") {iflag=1; continue;}
		if (str==")") {iflag=0; continue;}
		if (iflag==1){
			if(stoi(str)>-1) neighbors.push_back(stoi(str));
		}
	}
	stream.close();
}

void
trajectory::readOwners(void){
	string str;
	string strPre;
	int iflag=0;

	ifstream stream("constant/polyMesh/owner");
	while(getline(stream,str)) {
		if (str=="(") {iflag=1; continue;}
		if (str==")") {iflag=0; continue;}
		if (iflag==1){
			owners.push_back(stoi(str));
		}
	}
	stream.close();
}


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



void
trajectory::readParticles(void){
	string str;
	string strPre;
	int iflag=0;
	int particleID=0;
	ifstream stream("particle/particleSet");
	while(getline(stream,str)) {
		if (str=="x\ty\tz\tdp") {iflag=1; continue;}
		if (iflag==1){
			particle p;
			int loop=0;
			string tmp;
			istringstream stream(str);
			while(getline(stream,tmp,'\t')) {
				p.x.x[loop]=stod(tmp);
				p.v.x[loop]=stod(tmp);
				p.F.x[loop]=stod(tmp);
				if(loop==3) p.dp=stod(tmp);
				loop++;
			}
			p.id=particleID;
			p.m=rho_p*M_PI*p.dp*p.dp*p.dp/6.0;
			p.Zp=0;
			p.beta=0;
			p.Cc=0;
			p.Re=0;
			p.Mach=0;
			p.Kn=0;
			p.reflect=1;

			vars->particles.push_back(p);
			particleID++;
		}
	}
	stream.close();
}

void
trajectory::readCondition(void){
	string str;
	string strPre;
	ifstream stream("particle/condition");
	while(getline(stream,str)) {
		int readFlag=0;
		string tmp;
		istringstream stream(str);
		while(getline(stream,tmp,'\t')) {
			if (readFlag==0) {
				if(tmp=="dt") {readFlag=1; continue;}
				if(tmp=="totalStep") {cout<<"you can not use this argument"<<endl;}
				if(tmp=="totalTime") {readFlag=2; continue;}
				if(tmp=="dragModel") {readFlag=3; continue;}
				if(tmp=="Dispersion") {readFlag=4; continue;}
				if(tmp=="Observe") {cout<<"you can not use this argument"<<endl;}
				if(tmp=="observeTime") {readFlag=5; continue;}
				if(tmp=="FroudeKrylov") {readFlag=6; continue;}
				if(tmp=="dimension") {readFlag=7; continue;}
				if(tmp=="compressible") {readFlag=9; continue;}
				if(tmp=="constTemp") {readFlag=10; continue;}
				if(tmp=="constRho") {readFlag=11; continue;}
				if(tmp=="startDir") {readFlag=12; continue;}
				if(tmp=="particleDensity") {readFlag=13; continue;}
				if(tmp=="deltaTrajectory") {readFlag=14; continue;}
				if(tmp=="analytical") {readFlag=15; continue;}
			}
			if (readFlag==1) {
				if(tmp=="auto") flags->autoStep=1;
				else dt=stod(tmp);
			}
			if (readFlag==2) totalTime=stod(tmp);
			if (readFlag==3) {
				if(tmp=="Singh") forces.push_back(new dragForceSingh);
				if(tmp=="Stokes") forces.push_back(new dragForceSM);
				if(tmp=="Morsi") forces.push_back(new dragForceMA);
				if(tmp=="Loth") forces.push_back(new dragForceLoth);
			}
			if (readFlag==5) observeTime=stod(tmp);
			if (readFlag==6) {
				if(tmp=="Yes") flags->KFFlag=1;
				if(tmp=="No") flags->KFFlag=0;
			}
			if (readFlag==7) {
				if(tmp=="3D") flags->dimensionFlag=0;
				if(tmp=="2Dplane") {flags->dimensionFlag=1; readFlag=72;continue;}
				if(tmp=="2Daxi") {flags->dimensionFlag=2; readFlag=72;continue;}
			}
			if (readFlag==72) {
				if(tmp=="100") {plane2D=0; readFlag=73;continue;}
				if(tmp=="010") {plane2D=1; readFlag=73;continue;}
				if(tmp=="001") {plane2D=2; readFlag=73;continue;}
			}
			if (readFlag==73) {
				Axis=stod(tmp);
			}
			if (readFlag==9) {
				if(tmp=="Yes") flags->compressFlag=1;
				if(tmp=="No") flags->compressFlag=0;
			}
			if (readFlag==10) {
				constTemp=stod(tmp);
			}
			if (readFlag==11) {
				constRho=stod(tmp);
			}
			if (readFlag==12) {
				startDir=tmp;
			}
			if (readFlag==13) {
				rho_p=stod(tmp);
			}
			if (readFlag==15) {
				if(tmp=="Yes") {flags->analytical=1;readFlag=151;continue;}
				if(tmp=="No") {flags->analytical=0;}
			}
			if (readFlag==151) {
				vars->analyticalStep=stod(tmp);
			}
		}
	}
	stream.close();
}

void
trajectory::readGeometry(void){
	timer=clock();
	readNeighbors();
	readFaces();
	readPoints();
	readOwners();
	readBoundaries();
	cout<<"Read geometry time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
}

void
trajectory::readCFDresults(void){
	timer=clock();
	sprintf ( filepath, "%s/p", startDir.c_str());
	readScalar(filepath,vars->p);
	sprintf ( filepath, "%s/U", startDir.c_str());
	readVector(filepath,vars->U);

	/*  Compressible? */
	if(flags->compressFlag==1){
		sprintf ( filepath, "%s/T", startDir.c_str());
		readScalar(filepath,vars->T);
		sprintf ( filepath, "%s/rho", startDir.c_str());
		readScalar(filepath,vars->rho);
	}
	if(flags->compressFlag==0){
		sprintf (filepath, "%s/p", startDir.c_str());
		readScalarDum(filepath,vars->T, constTemp);
		readScalarDum(filepath,vars->rho, constRho);
	}

	/*  Dispersion? */
	if(flags->KFFlag==1){
		sprintf ( filepath, "%s/dp", startDir.c_str());
		readVector(filepath,vars->dp);
	}
	if(flags->KFFlag==0){
		sprintf ( filepath, "%s/p", startDir.c_str());
		readVectorDum(filepath,vars->dp, 0);
	}
	cout<<"Read CFD result time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
}
