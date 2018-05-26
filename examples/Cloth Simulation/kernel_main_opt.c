#include "cloth_code.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void initMatrix(int n, double mass, double fcon, 
		int delta, double grav, double sep, double rball,
		double offset,
		double dt, double **x, double **y, double **z,
		double** cpx, double** cpy, double** cpz,
		double **fx, double **fy, double **fz,
		double **vx, double **vy, double **vz,
		double **oldfx, double **oldfy, double **oldfz)
{
  int i,nx,ny;
  
  // Free any existing
  free(*x); free(*y); free(*z);
  free(*cpx); free(*cpy); free(*cpz);
  
  // allocate arrays to hold locations of nodes
  *x=(double*)malloc(n*n*sizeof(double));
  *y=(double*)malloc(n*n*sizeof(double));
  *z=(double*)malloc(n*n*sizeof(double));
  //This is for opengl stuff
  *cpx=(double*)malloc(n*n*sizeof(double));
  *cpy=(double*)malloc(n*n*sizeof(double));
  *cpz=(double*)malloc(n*n*sizeof(double));

  //initialize coordinates of cloth
  for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny++){
      (*x)[n*nx+ny] = nx*sep-(n-1)*sep*0.5+offset;
      (*z)[n*nx+ny] = rball+1;
      (*y)[n*nx+ny] = ny*sep-(n-1)*sep*0.5+offset;
      (*cpx)[n*nx+ny] = 0;
      (*cpz)[n*nx+ny] = 1;
      (*cpy)[n*nx+ny] = 0;
    }
  }

  // Throw away existing arrays
  free(*fx); free(*fy); free(*fz);
  free(*vx); free(*vy); free(*vz);
  free(*oldfx); free(*oldfy); free(*oldfz);
  // Alloc new
  *fx=(double*)malloc(n*n*sizeof(double));
  *fy=(double*)malloc(n*n*sizeof(double));
  *fz=(double*)malloc(n*n*sizeof(double));
  *vx=(double*)malloc(n*n*sizeof(double));
  *vy=(double*)malloc(n*n*sizeof(double));
  *vz=(double*)malloc(n*n*sizeof(double));
  *oldfx=(double*)malloc(n*n*sizeof(double));
  *oldfy=(double*)malloc(n*n*sizeof(double));
  *oldfz=(double*)malloc(n*n*sizeof(double));
  for (i=0;i<n*n;i++){
    (*vx)[i]=0.0;
    (*vy)[i]=0.0;
    (*vz)[i]=0.0;
    (*fx)[i]=0.0;
    (*fy)[i]=0.0;
    (*fz)[i]=0.0;
  }
}

void loopcode(int n, double mass, double fcon,
	      int delta, double grav, double sep, double rball,
	      double xball, double yball, double zball,
	      double dt, double *x, double *y, double *z,
	      double *fx, double *fy, double *fz,
	      double *vx, double *vy, double *vz,
	      double *oldfx, double *oldfy, double *oldfz, 
	      double *pe, double *ke,double *te, int no_threads)
{
  int i,j;
  double rlen,xdiff,ydiff,zdiff,vmag,damp,rx,ry,rz,rmag;
	
  //update position as per MD simulation
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      x[j*n+i] +=dt*(vx[j*n+i]+dt*fx[j*n+i]*0.5);
      oldfx[j*n+i]=fx[j*n+i];
    }
  }
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      y[j*n+i] +=dt*(vy[j*n+i]+dt*fy[j*n+i]*0.5);
      oldfy[j*n+i]=fy[j*n+i];
    }
  }
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      z[j*n+i] +=dt*(vz[j*n+i]+dt*fz[j*n+i]*0.5);
      oldfz[j*n+i]=fz[j*n+i];
    }
  }

  //	apply constraints - push cloth outside of ball, set to zero velocity
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      xdiff = x[j*n+i]-xball;
      ydiff = y[j*n+i]-yball;
      zdiff = z[j*n+i]-zball;
      vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);

      if (vmag < rball){

          x[j*n+i]=xball+xdiff*rball/vmag;
          y[j*n+i]=yball+ydiff*rball/vmag;
          z[j*n+i]=zball+zdiff*rball/vmag;
          
          xdiff = -1*xdiff;
          ydiff = -1*ydiff;
          zdiff = -1*zdiff;

          rx = xdiff/vmag;
          ry = ydiff/vmag;
          rz = zdiff/vmag;

          rmag = sqrt(rx*rx - ry*ry - rz*rz);

          rx = rx/rball;
          ry = ry/rball;
          rz = rz/rball;


          vx[j*n+i]= vx[j*n+i] -(rx * vx[j*n+i] * rx);
          vy[j*n+i]= vy[j*n+i] - (ry * vy[j*n+i] * ry);
          vz[j*n+i]= vz[j*n+i]- (rz * vz[j*n+i] * rz);
      }
    }
  }

  *pe=eval_pef(n,delta,grav,sep,fcon,x,y,z,fx,fy,fz,no_threads);

  //Add a damping factor to eventually set velocity to zero
  damp=0.995;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      vx[j*n+i]=(vx[j*n+i]+dt*(fx[j*n+i]+oldfx[j*n+i])*0.5)*damp;
      vy[j*n+i]=(vy[j*n+i]+dt*(fy[j*n+i]+oldfy[j*n+i])*0.5)*damp;
      vz[j*n+i]=(vz[j*n+i]+dt*(fz[j*n+i]+oldfz[j*n+i])*0.5)*damp;
    }
  }
  *ke=0.0;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      *ke += vx[j*n+i]*vx[j*n+i] + vy[j*n+i]*vy[j*n+i] + vz[j*n+i]*vz[j*n+i];
    }
  }
  *ke=*ke/2.0;
  *te=*pe+*ke;

}

double eval_pef(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,		
                double *fx, double *fy, double *fz, int no_threads)
{
  double pe,rlen,xdiff,ydiff,zdiff,vmag;
  int nx, ny, dx, dy;
  pe = 0.0;
  //loop over particles
  for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny++){
      fx[nx*n+ny]=0.0;
    }
  } 
  for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny++){
      fy[nx*n+ny]=0.0;
    }
  }
  for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny++){
      fz[nx*n+ny]=-grav;
    }
  }
  
      
    
  for (ny=0;ny<n;ny++){
    for (nx=0;nx<n;nx++){
      
      //loop over displacements

	    for (dx=MAX(nx-delta,0);dx<MIN(nx+delta+1,n);dx++){
        for (dy=MAX(ny-delta,0);dy<MIN(ny+delta+1,n);dy++){
          // exclude self interaction
          if (nx!=dx || ny!=dy){
            // compute reference distance
            rlen=sqrt((double)((nx-dx)*(nx-dx)+(ny-dy)*(ny-dy)))*sep;
            // compute actual distance
            xdiff = x[dx*n+dy]-x[nx*n+ny];
            ydiff = y[dx*n+dy]-y[nx*n+ny];
            zdiff = z[dx*n+dy]-z[nx*n+ny];
            vmag=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);
            //potential energy and force
            pe += fcon*(vmag-rlen);
            fx[nx*n+ny]+=fcon*xdiff*(vmag-rlen)/vmag;
            fy[nx*n+ny]+=fcon*ydiff*(vmag-rlen)/vmag;
            fz[nx*n+ny]+=fcon*zdiff*(vmag-rlen)/vmag;
          }
        }
      }
    }
  }
  return pe;
} 


