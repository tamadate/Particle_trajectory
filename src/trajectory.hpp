#pragma once
#include "dragForce.hpp"
#include "variables.hpp"
#include "flags.hpp"

class trajectory{
	public:
		Variables *vars;
		Flags *flags;
		std::vector<dragForceSM*> forces;

		// Calculation conditions (flag, file path, etc...)
		double dt;
		double rho_p;
		double observeTime;
		double totalTime;
		int Observe;
		double delta_r2_trajectory;
		int plane2D;
		string startDir;
		double Axis;
		double flag;
		char filepath[100];
		int boundaryStartID;
		double constTemp;
		double constRho;
		double timer;
		FILE*f;
		FILE*file;

		double get_tini(particle &par);
		void updateDisp(particle &par);

		void run(void);
		double timeEvolution(particle &a);
		double euler(particle &a);
		double analytical(particle &a);

		int checkCell(int pid);
		int boundAction(int faceID, int pid, point norm);
		void calculateMyu(void);
		void calculatelamda(void);
		void computeReMach(particle &par);

		std::vector<face> faces;
		std::vector<point> points;
		std::vector<boundary> boundaries;
		std::vector<cell> cells;
		std::vector<int> neighbors;
		std::vector<int> owners;
		std::vector<outParticle> outParticles;
		std::vector<int> erasePID;
		void makeCells(void);
		void readFaces(void);
		void readBoundaries(void);
		void readPoints(void);
		void readNeighbors(void);
		void readOwners(void);
		void readScalar(char *readFile, std::vector<double> &variable);
		void readScalarDum(char *readFile, std::vector<double> &variable, double value);
		void readVector(char *readFile, std::vector<point> &variable);
		void readVectorDum(char *readFile, std::vector<point> &variable, double value);
		void readParticles(void);
		void readCondition(void);
		void readGeometry(void);
		void readCFDresults(void);
		void initialParticle(void);
		void findParticle(void);
		void findParticleFace(std::vector<cell> cells);

		void output(particle a, double time);
		void outputTrajectory(particle a);
		void outputInitial(void);
		void outputFinalPosition(outParticle a);
		void outputInitialPosition(particle p);
		void outputPenetration(particle p);

		void initialize(particle &par, int nth){
			if(flags->dispersionFlag) par.update=1;
			else par.tini=1e10;
			outputInitialPosition(par);
		}

		trajectory(void);
		~trajectory(void);
	private:
};
