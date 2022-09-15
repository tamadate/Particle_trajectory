#include "trajectory.hpp"


// Read OpenFOAM style simulation data via below functions

/******************************************************************************/
// conditions file reading
// ./particle/condition file
void
trajectory::readCondition(void){
	string str;
	string strPre;
	ifstream stream("particle/condition");
	while(getline(stream,str)) {
		int readFlag=0;
		string tmp;
		istringstream stream(str);
		std::vector<string> readings;
		string reading;
		while(getline(stream,reading,'\t')) {
			readings.push_back(reading);
		}
		while(getline(stream,tmp,'\t')) {
			if (readFlag==0) {
				if(readings[0]=="dt") {
					if(readings[1]=="auto") flags->autoStep=1;
					else dt=stod(readings[1]);
				}
				else if(readings[0]=="dragModel") {
					if(readings[1]=="Singh") forces.push_back(new dragForceSingh);
					if(readings[1]=="Stokes") forces.push_back(new dragForceSM);
					if(readings[1]=="Morsi") forces.push_back(new dragForceMA);
					if(readings[1]=="Loth") forces.push_back(new dragForceLoth);
				}
				else if(readings[0]=="Dispersion") {
					if(readings[1]=="Yes") flags->dispersionFlag=1;
					if(readings[1]=="No") flags->dispersionFlag=0;
				}
				else if(readings[0]=="FroudeKrylov") {
					if(readings[1]=="Yes") flags->KFFlag=1;
					if(readings[1]=="No") flags->KFFlag=0;
				}
				else if(readings[0]=="dimension") {
					if(readings[1]=="3D") flags->dimensionFlag=0;
					if(readings[1]=="2Dplane") flags->dimensionFlag=1;
					if(readings[1]=="2Daxi") flags->dimensionFlag=2;
					if(flags->dimensionFlag>0){
						if(readings[2]=="100") plane2D=0;
						if(readings[2]=="010") plane2D=1;
						if(readings[2]=="001") plane2D=2;
						Axis=stod(readings[3]);
					}
				}
				else if(readings[0]=="compressible") {
					if(readings[1]=="Yes") flags->compressFlag=1;
					if(readings[1]=="No") flags->compressFlag=0;
				}
				else if(readings[0]=="analytical") {
					if(readings[1]=="Yes") {
						flags->analytical=1;
						vars->analyticalStep=stod(readings[2]);
					}
					else if(readings[1]=="No") flags->analytical=0;
				}
				else if(readings[0]=="totalTime") totalTime=stod(readings[1]);
				else if(readings[0]=="observeTime") observeTime=stod(readings[1]);
				else if(readings[0]=="constTemp") constTemp=stod(readings[1]);
				else if(readings[0]=="constRho") constRho=stod(readings[1]);
				else if(readings[0]=="startDir") startDir=readings[1];
				else if(readings[0]=="particleDensity") rho_p=stod(readings[1]);
				else if(readings[0]=="deltaTrajectory") {readFlag=14; continue;}
				else if(readings[0]=="inletFace") flags->inletFace=stoi(readings[1]);
				else cout<<"Could not find syntax "<<readings[0]<<endl;
			}

		}
	}
	stream.close();
}

/******************************************************************************/
// CFD results reading
// p, U, T, etc... in startDir (default 20000)
void
trajectory::readCFDresults(void){
	timer=clock();

// read CFD simulation results using readScalar and read readVector functions (see readingFunctions.cpp)

	// read pressure and velocity field
	sprintf ( filepath, "%s/p", startDir.c_str());
	readScalar(filepath,vars->p);
	sprintf ( filepath, "%s/U", startDir.c_str());
	readVector(filepath,vars->U);

	// read temperature and density field
	// if it is incompressible, read pressure in stead
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

	// read k and epsilon field
	// if it is laminar, read pressure in stead
	if(flags->dispersionFlag){
		sprintf ( filepath, "%s/k", startDir.c_str());
		readScalar(filepath,vars->k);
		sprintf ( filepath, "%s/epsilon", startDir.c_str());
		readScalar(filepath,vars->epsilon);
	}
	if(flags->dispersionFlag==0){
		sprintf (filepath, "%s/p", startDir.c_str());
		readScalarDum(filepath,vars->k, 0);
		readScalarDum(filepath,vars->epsilon, 0);
	}

	// read delta-p for Froude Krylov force
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


/******************************************************************************/
// particle initial locations and diameters reading
// particle/particleSet file
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

/******************************************************************************/
// Geometry reading
// ./constant/polyMesh/* files
void
trajectory::readGeometry(void){
	timer=clock();

	readNeighbors(); // ./constant/polyMesh/neighbor
	readFaces();		// ./constant/polyMesh/faces
	readPoints();		// ./constant/polyMesh/points
	readOwners();		// ./constant/polyMesh/owner
	readBoundaries();	// ./constant/polyMesh/boundary

	cout<<"Read geometry time: "<<(clock()-timer)*1e-6<<" sec"<<endl;
}


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
