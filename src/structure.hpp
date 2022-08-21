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


struct particle {
	point x;
	point v;
	point F;
	int cell;
	int id;
	double Re;
	double Mach;
	double Kn;
	double dp;
	double m;
	double Cc;
	double Zp;
	double beta;
	double dt;
    int reflect;
};
