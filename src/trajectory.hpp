#pragma once
#include "solver/Hybrid.hpp"
#include "solver/EulerIL.hpp"

class trajectory{
	public:
	Variables *vars;	// pointer for variables
	Flags *flags;			// pointer for flags
	Solver *solver;

	double rho_p;									// particle density
	double observeTime;						// output interval in time scale
	double totalTime;							// total calculation time
	int Observe;									// output interval in time step
	double delta_r2_trajectory;		// squre of migration distance for trajectory output
	int noUpdateAxis;
	string startDir;							// start directory name (20000 is default)
	char filepath[100];						// file path
	int boundaryStartID;					// boundary start face id
	double timer;
	FILE*f;
	FILE*file;
	std::vector<outParticle> outParticles; // position & velocity of the particles at the end of the simulation
	std::vector<activeFace> penetrates;
	std::vector<activeFace> reflects;
	std::vector<int> trapParticle;				 // number of particles trapped in the curculation
	std::vector<std::vector<double>> deadBox;
	std::vector<deadSpaceDist> deadDistance;

	// Main simulation functions
	void run(void);
	int checkCell(int pid);
	int checkBoundCell(int pid);
	int boundAction(int faceID, int pid, point norm, double dot);

	// Physical properties calculation functions (see calcProperties.cpp)
	void calculateMyu(void);
	void calculatelamda(void);
	void calculateNonDimension(particle &par);

	// Simualtion geometrical information
	std::vector<face> faces;
	std::vector<point> points;
	std::vector<boundary> boundaries;
	std::vector<cell> cells;
	std::vector<int> neighbors;
	std::vector<int> owners;
	void makeCells(void); // make cells from abouve geometrical arrays

	// Reading functions
	void readGeometry(void);
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

	// Initialization
	void initialParticle(void);
	void findParticle(void);
	void findParticleFace(std::vector<cell> cells);
	void setInitialVelocity(void);
	void setDefault(void);

	// Export functions
	void output(particle a, double time);
	void outputTrajectory(particle a);
	void outputInitial(void);
	void outputFinalPosition(void);
	void outputPenetrate(void);
	void outputInitialVelocity(void);

	bool is_file_exist(const char *fileName)
	{
		std::ifstream infile(fileName);
		return infile.good();
	}

	void initialize(particle &par){}

	trajectory(void);
	~trajectory(void){};

	private:
};
