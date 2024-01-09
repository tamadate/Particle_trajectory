#include "../trajectory.hpp"


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
		while(getline(stream,reading,' ')) {
			readings.push_back(reading);
		}

		// time step
		if(readings[0]=="dt") {
			if(readings[1]=="auto") {
				flags->autoStep=true;
				vars->meshScale=stod(readings[2]);
				cout<<"Auto time step: "<<vars->meshScale<<endl;
			}
			else if(readings[1]=="fix"){
				flags->autoStep=false;
				vars->fixTimeStep=stod(readings[2]);
				cout<<"Fixed time step: "<<vars->dt<<endl;
			}
		}
		// gas type (currently you can select He or Air)
		/*else if(readings[0]=="gasType") {
			if(readings[1]=="He") forces.push_back(new dragForceSingh);
			if(readings[1]=="Air") forces.push_back(new dragForceSM);
			cout<<"Gas type: "<<readings[1]<<endl;
		}*/
		// model for the drag coefficient model
		else if(readings[0]=="dragModel") {
			/*delete drag;
			if(readings[1]=="Singh") drag = new dragForceSingh;
			if(readings[1]=="Stokes") drag = new dragForceSM;
			if(readings[1]=="Morsi") drag = new dragForceMA;
			if(readings[1]=="Loth") drag = new dragForceLoth;*/
			if(readings[1]=="Singh") forces.push_back(new dragForceSingh);
			if(readings[1]=="Stokes") forces.push_back(new dragForceSM);
			if(readings[1]=="Morsi") forces.push_back(new dragForceMA);
			if(readings[1]=="Loth") forces.push_back(new dragForceLoth);
			cout<<"Drag model: "<<readings[1]<<endl;
		}
		// select diffusion on/off
		else if(readings[0]=="diffusion") {
			if(readings[1]=="Yes") forces.push_back(new Langevin);
			cout<<"Diffusion: "<<readings[1]<<endl;
		}
		// turbulent dispersion on or off
		// you can select k-e or k-w as a RANS model
		else if(readings[0]=="Dispersion") {
			if(readings[1]=="k-e") {
				flags->dispersionFlag=1;
				sprintf ( filepath, "%s/k", startDir.c_str());
				readScalar(filepath,vars->k);
				sprintf ( filepath, "%s/epsilon", startDir.c_str());
				readScalar(filepath,vars->epsilon);
			}
			if(readings[1]=="k-w") {
				flags->dispersionFlag=2;
				sprintf ( filepath, "%s/k", startDir.c_str());
				readScalar(filepath,vars->k);
				sprintf ( filepath, "%s/omega", startDir.c_str());
				readScalar(filepath,vars->epsilon);
				int Ncell=vars->epsilon.size();
				for(int i=0; i<Ncell; i++){
					vars->epsilon[i]*=vars->k[i];
				}
			}
			if(readings[1]=="No" || readings[1]=="no") {
				flags->dispersionFlag=0;
				sprintf (filepath, "%s/p", startDir.c_str());
				readScalarDum(filepath,vars->k, 0);
				readScalarDum(filepath,vars->epsilon, 0);
			}
			cout<<"Dispersion: "<<readings[1]<<endl;
		}
		else if(readings[0]=="CoulombForce") {
			flags->dispersionFlag=1;
			sprintf ( filepath, "%s/E", startDir.c_str());
			readVector(filepath,vars->dV);
			forces.push_back(new Coulomb);
			cout<<"Coulomb force ON"<<endl;
		}
		// Froude-Krylov force (force caused by the pressure difference)
      	else if(readings[0]=="FroudeKrylov") {
        	if(readings[1]=="Yes" || readings[1]=="yes") {
				sprintf ( filepath, "%s/dp", startDir.c_str());
				readVector(filepath,vars->dp);
			}
	        if(readings[1]=="No" || readings[1]=="no") {
				sprintf ( filepath, "%s/p", startDir.c_str());
				readVectorDum(filepath,vars->dp, 0);
			}
		}
		// dimension of the CFD simulation
      	else if(readings[0]=="dimension") {
        	if(readings[1]=="3D") flags->dimensionFlag=0;
			if(readings[1]=="2Dplane") flags->dimensionFlag=1;
			if(readings[1]=="2Daxi") flags->dimensionFlag=2;
			cout<<"Dimension: "<<readings[1]<<endl;
			if(flags->dimensionFlag>0){
				if(readings[2]=="100") plane2D=0;
				if(readings[2]=="010") plane2D=1;
				if(readings[2]=="001") plane2D=2;
				Axis=stod(readings[3]);
        	}
		}
		// tracking intermediate properties
		else if(readings[0]=="penetrate") {
			penetrate p;
			if(readings[1]=="100") p.face=0;
			if(readings[1]=="010") p.face=1;
			if(readings[1]=="001") p.face=2;
			p.loc=stod(readings[2]);
			penetrates.push_back(p);
		}
		// fluid compressibility yes or no
      	else if(readings[0]=="compressible") {
			if(readings[1]=="Yes" || readings[1]=="yes") {
				sprintf ( filepath, "%s/T", startDir.c_str());
				readScalar(filepath,vars->T);
				sprintf ( filepath, "%s/rho", startDir.c_str());
				readScalar(filepath,vars->rho);
				cout<<"Using T and rho distribution files"<<endl;
			}
			if(readings[1]=="No" || readings[1]=="no") {
				sprintf (filepath, "%s/p", startDir.c_str());
				readScalarDum(filepath,vars->T, stod(readings[2]));
				readScalarDum(filepath,vars->rho, stod(readings[3]));
				cout<<"Fixed T and rho: "<<readings[2]<<" and "<<readings[3]<<endl;
			}
		}
		// set particle initial velocities (default is fluid)
		// fluid: 99% of the fluid velocity
		// file: velocities are given in the file
		else if(readings[0]=="initialVelocity") {
			flags->v0=1;
			cout<<"Reading initial velocity file"<<endl;
		}
		// not using now
      	else if(readings[0]=="analytical") {
			if(readings[1]=="Yes" || readings[1]=="yes") {
				flags->analytical=1;
				vars->analFactor=stod(readings[2]);
				cout<<"Use analytical mode: "<<vars->analFactor<<endl;
			}
			else if(readings[1]=="No" || readings[1]=="no") {
				flags->analytical=0;
				cout<<"Not use analytical mode"<<endl;
			}
		}
		// total simulation time step
		else if(readings[0]=="totalTime") totalTime=stod(readings[1]);
		// not using now
		else if(readings[0]=="observeTime") observeTime=stod(readings[1]);
		// set name of CFD simulation result (20000 is the default)
		else if(readings[0]=="startDir") {
			startDir=readings[1];
			sprintf ( filepath, "%s/p", startDir.c_str());
			readScalar(filepath,vars->p);
			sprintf ( filepath, "%s/U", startDir.c_str());
			readVector(filepath,vars->U);
			cout<<"Start directory: "<<readings[1]<<endl;
		}
		// gravity setting
		else if(readings[0]=="gravity") {
			gravity *grav = new gravity(stod(readings[1]),stod(readings[2]),stod(readings[3])); // generate variables class
			forces.push_back(grav);
			cout<<"Gravity: "<<readings[1]<<" m2/s, "<<readings[2]<<" m2/s, "<<readings[3]<<" m2/s"<<endl;
		}
		// Interception
		else if(readings[0]=="interception" || readings[0]=="Interception") {
			if(readings[1]=="distance") {
				deadSpaceDist dead;
				dead.x[0]=stod(readings[2]);
				dead.x[1]=stod(readings[3]);
				dead.x[2]=stod(readings[4]);
				dead.r=stod(readings[5]);
				dead.bid=stoi(readings[6]);
				deadDistance.push_back(dead);
				cout<<"Set dead volume x=("<<dead.x[0]<<" "<<dead.x[1]<<" "<<dead.x[2]<<") r="<<dead.r<<endl;
			}
			else if(readings[1]=="box") {
				std::vector<double> dead;
				dead.push_back(stod(readings[1]));
				dead.push_back(stod(readings[2]));
				dead.push_back(stod(readings[3]));
				dead.push_back(stod(readings[4]));
				dead.push_back(stod(readings[5]));
				dead.push_back(stod(readings[6]));
				dead.push_back(stod(readings[7]));
				deadBox.push_back(dead);
				cout<<"set dead volume (box)"<<endl;
			}
		}
		// particle density (default is 1000 kg/m3)
		else if(readings[0]=="particleDensity") rho_p=stod(readings[1]);
		
		// error check
		else cout<<"Could not find syntax "<<readings[0]<<endl;
	}
	stream.close();
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
