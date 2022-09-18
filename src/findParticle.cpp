#include "trajectory.hpp"

// Find cell id where particle is existing

void
trajectory::findParticle(void){
	if(flags->inletFace==-1) findParticleFace(cells); // currently this is not working
	else {
		std::vector<cell> boundCells;
		int cellSize=cells.size();
		int bID=flags->inletFace;
		int bID0=boundaries[bID].startFace;
		int bIDend=boundaries[bID].startFace+boundaries[bID].nFaces;
		for(int j=0; j<cellSize; j++) {
			int faceSize=cells[j].iface.size();
			for (int k=0; k<faceSize; k++){
				int b=cells[j].iface[k];
				if(b+1 > bID0  && b < bIDend)	boundCells.push_back(cells[j]);
			}
		}
		findParticleFace(boundCells);
	}
	outputInitialVelocity();
}

void
trajectory::findParticleFace(std::vector<cell> targetCells){
	timer=omp_get_wtime();
	int ps=vars->particles.size();
	# pragma omp parallel
	{
		# pragma omp for
		for(int i=0; i<ps; i++){
			int foundCells=0; // number of found cells (sometimes this code find multiple cells for a particle)
			double x=vars->particles[i].x.x[0];
			double y=vars->particles[i].x.x[1];
			double z=vars->particles[i].x.x[2];
			int cellSize=cells.size();
	    double minDist=1e15;

			/*
				This loop calculate the inner product of particle position vector with norm vector of each faces
				regarding with every cells.
				If the every inner product is positive on a cell, it means the particle exists in that cell and
				it becomes initial cell id.
			*/
			for(int j=0; j<cellSize; j++) {
				cell a=cells[j];
				int faceSize=a.iface.size();
				int flag=0;
				int btest=a.iface[0];
				int c0test=faces[btest].iface[0];
				double y0test=points[c0test].x[1];
				if(y0test<0.088) continue;
				for (int k=0; k<faceSize; k++){
					int b=a.iface[k];
					int c0=faces[b].iface[0];
					double x0=points[c0].x[0];
					double y0=points[c0].x[1];
					double z0=points[c0].x[2];
					point c=a.norm[k];
					double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
					if(naiseki<0) {flag=1;break;}
				}
				if(flag==0){
					double Dist=0;
					for (int k=0; k<faceSize; k++){
						int b=a.iface[k];
						int c0=faces[b].iface[0];
						double x0=points[c0].x[0];
						double y0=points[c0].x[1];
						double z0=points[c0].x[2];
						point c=a.norm[k];
						double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
	                    Dist+=naiseki;
					}
					foundCells++;
					if (Dist<minDist){
						minDist=Dist;
						vars->particles[i].cell=j;
					}
				}
			}


			/*
				Sometimes this code can not find any cells (foundCells=1).
				In that case, the initial cell id is determined by the distance between particle and cell.
				(initial cell may give minimum distance with particle).
			*/
			if (foundCells==0){
				cout<<"**Error: I could not find an initial cell for pid "<<i<<endl; // Error message for just in case
				for(int j=0; j<cellSize; j++) {
					cell a=targetCells[j];
					int faceSize=a.iface.size();
					int flag=0;
					for (int k=0; k<faceSize; k++){
						int b=a.iface[k];
						int c0=faces[b].iface[0];
						double x0=points[c0].x[0];
						double y0=points[c0].x[1];
						double z0=points[c0].x[2];
						point c=a.norm[k];
						double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
						if(naiseki<0) {flag=1;break;}
					}
					double Dist=0;
					for (int k=0; k<faceSize; k++){
						int b=a.iface[k];
						int c0=faces[b].iface[0];
						double x0=points[c0].x[0];
						double y0=points[c0].x[1];
						double z0=points[c0].x[2];
						point c=a.norm[k];
						double naiseki=c.x[0]*(x-x0)+c.x[1]*(y-y0)+c.x[2]*(z-z0);
            Dist+=naiseki;
					}
					if (Dist<minDist){
						minDist=Dist;
						vars->particles[i].cell=j;
					}
				}
			}
		}
	}
	cout<<"Find particle time: "<<(omp_get_wtime()-timer)<<" sec"<<endl;
}
