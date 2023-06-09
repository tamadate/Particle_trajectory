#include "../trajectory.hpp"

void
trajectory::makeCells(void){
	timer=omp_get_wtime();

	// get maximum number on owners = number of cells
	int maxOwner=0;
	for (auto &a:owners) {if(a>maxOwner) maxOwner=a;}

	// make cell array
	// add face id (i) into a cell (cell id is owners[i])
	// add face id (i) into a cell (cell id is neighbors[i])
	cell c;
	for (int i=0; i<maxOwner+1; i++) {cells.push_back(c);}
	for (int i=0; i<owners.size(); i++) cells[owners[i]].iface.push_back(i);
	for (int i=0; i<neighbors.size(); i++) cells[neighbors[i]].iface.push_back(i);

	// calculate center of cells
	for (auto &a:cells) {
		int faceSize=a.iface.size();	// number of faces in this cell
		int sumN=0;
		double X=0;
		double Y=0;
		double Z=0;
		for(auto &b: a.iface){
			for(auto &c: faces[b].iface){
				X+=points[c].x[0];
				Y+=points[c].x[1];
				Z+=points[c].x[2];
				sumN+=1;
			}
		}
		X/=double(sumN);
		Y/=double(sumN);
		Z/=double(sumN);
		a.r.x[0]=X;
		a.r.x[1]=Y;
		a.r.x[2]=Z;

		// calculate norm vectors of each face
		for (int i=0; i<faceSize; i++){
			int b=a.iface[i];
			int c0=faces[b].iface[0];
			int c1=faces[b].iface[1];
			int c2=faces[b].iface[2];

			double x0=points[c0].x[0];
			double y0=points[c0].x[1];
			double z0=points[c0].x[2];

			double x1=points[c1].x[0];
			double y1=points[c1].x[1];
			double z1=points[c1].x[2];

			double x2=points[c2].x[0];
			double y2=points[c2].x[1];
			double z2=points[c2].x[2];

			double Ax=x1-x0;
			double Ay=y1-y0;
			double Az=z1-z0;

			double Bx=x2-x0;
			double By=y2-y0;
			double Bz=z2-z0;

			point norm;
			norm.x[0]=Ay*Bz-Az*By;
			norm.x[1]=Az*Bx-Ax*Bz;
			norm.x[2]=Ax*By-Ay*Bx;

			double x3=X-x0;
			double y3=Y-y0;
			double z3=Z-z0;

			double dot=norm.x[0]*x3+norm.x[1]*y3+norm.x[2]*z3;
			if(dot<0){
				norm.x[0]=By*Az-Bz*Ay;
				norm.x[1]=Bz*Ax-Bx*Az;
				norm.x[2]=Bx*Ay-By*Ax;
			}
			dot=norm.x[0]*x3+norm.x[1]*y3+norm.x[2]*z3;
			if(dot<0){cout<<"error1011"<<endl;}
			double r2=norm.x[0]*norm.x[0]+norm.x[1]*norm.x[1]+norm.x[2]*norm.x[2];
			double r=sqrt(r2);
			norm.x[0]/=r;
			norm.x[1]/=r;
			norm.x[2]/=r;

			a.norm.push_back(norm); // norm vector is stored in cell
		}

	}
	cout<<"Making cell time: "<<omp_get_wtime()-timer<<" sec"<<endl;
}
