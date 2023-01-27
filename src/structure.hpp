#pragma once

struct face {
	std::vector<int> iface;
};

struct point {
	double x[3];
};

struct boundary {
	string name;
	string type;
	int nFaces;
	int startFace;
};

struct cell {
	std::vector<int> iface;
	std::vector<point> norm;
	point r;
};

struct outParticle {
	int pid;
  int bid;
	point r;
	point v;
};

struct penetrate {
	std::vector<outParticle> outPositions; // position & velocity of the particles at the end of the simulation
	std::vector<double> dx0;
	double loc;
	int face;
};

struct particle {
	point x;	// position
	point v;	// velocity
	point F;	// force
	int cell;	// cell id
	int id;		// particle id
	double Re;	// Reynolds number
	double Mach;	// Mach number
	double Kn;	//Kndsen number
	double dp;	// diameter
	double m;		// mass
	double Cc;	// slip coefficient
	double Zp;	// ion mobility
	double beta;	// relaxation time (SM)
	double dt;  // time step

	point Urand;	// dispersion velocity
	double tini;	// dispersion velocity life time
	int update;		// update flag

  int reflect;	// reflect flag for symmetric simulation
};
