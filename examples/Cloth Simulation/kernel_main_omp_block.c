#include "cloth_code.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <emmintrin.h>
#include<omp.h>

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

void initializeX(int n,double *x,double dt, double *oldfx,double *vx, double *fx){
  int i,j;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      x[j*n+i] +=dt*(vx[j*n+i]+dt*fx[j*n+i]*0.5);
      oldfx[j*n+i]=fx[j*n+i];
    }
  }
}

void initializeY(int n,double *y,double dt, double *oldfy,double *vy, double *fy){
  int i,j;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      y[j*n+i] +=dt*(vy[j*n+i]+dt*fy[j*n+i]*0.5);
      oldfy[j*n+i]=fy[j*n+i];
    }
  }
}

void initializeZ(int n,double *z,double dt, double *oldfz,double *vz, double *fz){
  int i,j;
  for (j=0;j<n;j++){
    for (i=0;i<n;i++){
      z[j*n+i] +=dt*(vz[j*n+i]+dt*fz[j*n+i]*0.5);
      oldfz[j*n+i]=fz[j*n+i];
    }
  }
}

void loopcode_sub(int i, int n, double mass, double fcon,
	      int delta, double grav, double sep, double rball,
	      double xball, double yball, double zball,
	      double dt, double *x, double *y, double *z,
	      double *fx, double *fy, double *fz,
	      double *vx, double *vy, double *vz,
	      double *oldfx, double *oldfy, double *oldfz, 
	      double *pe, double *ke,double *te){

  float ball_v[3] __attribute__ ((aligned (64))) = {(float)xball, (float)yball,(float)zball} ;
    int j;
  double rlen,xdiff,ydiff,zdiff,vmag;
	
  for(j=0;j<n;j++){
      float diff_v[3] __attribute__ ((aligned (64))) = {(float)x[j*n+i],(float) y[j*n+i],(float)z[j*n+i]} ;
      __m128 sse_current = _mm_load_ps(diff_v);
      __m128 sse_ball = _mm_load_ps(ball_v);
      __m128 sse_point_ball_diff = _mm_sub_ps(sse_current,sse_ball);
      __m128 sse_diff_mag = _mm_mul_ps(sse_point_ball_diff,sse_point_ball_diff);
      float sum_v[3] __attribute__ ((aligned (64))) = {0.0, 0.0, 0.0} ;
      _mm_store_ps(sum_v, sse_diff_mag);
      vmag = sqrt(sum_v[0]+sum_v[1]+sum_v[2]);
     
      if (vmag < rball){
          float tmp[3] __attribute__ ((aligned (64))) = {(float)(rball/vmag),(float) (rball/vmag),(float)(rball/vmag)} ;
          __m128 tmp_sse = _mm_load_ps(tmp);
          tmp_sse = _mm_mul_ps(tmp_sse,sse_point_ball_diff);
          sse_current = _mm_add_ps(sse_ball,tmp_sse);
          _mm_store_ps(sum_v, sse_current);
          
          x[j*n+i]=sum_v[0];
          y[j*n+i]=sum_v[1];
          z[j*n+i]= sum_v[2];
          
          float one_v[3] __attribute__ ((aligned (64))) = {1.0, 1.0, 1.0} ;
          __m128 sse_ones = _mm_load_ps(one_v); 
          sse_point_ball_diff = _mm_mul_ps(sse_point_ball_diff, sse_ones);

          float vmag_v[3] __attribute__ ((aligned (64))) = {vmag, vmag, vmag} ;

          __m128 sse_vmag = _mm_load_ps(vmag_v);
           __m128 sse_r = _mm_div_ps(sse_point_ball_diff, sse_vmag);

          float ball_v[3] __attribute__ ((aligned (64))) = {rball, rball, rball} ;
          sse_r = _mm_div_ps(sse_r, _mm_load_ps(ball_v));

          float v_v[3] __attribute__ ((aligned (64))) = {vx[j*n+i],  vy[j*n+i], vz[j*n+i]} ;

        __m128 sse_v = _mm_load_ps(v_v);
         sse_v = _mm_sub_ps(sse_v,(_mm_mul_ps(sse_v,_mm_mul_ps(sse_r,sse_r))));
          _mm_store_ps(sum_v, sse_v);

          vx[j*n+i] = sum_v[0];
          vy[j*n+i] = sum_v[1];
          vz[j*n+i] = sum_v[2];
      }
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
  double damp;
	omp_set_num_threads(no_threads);
  #pragma omp parallel 
  {
  #pragma omp single
  initializeX(n,x,dt,oldfx,vx,fx);
  #pragma omp single
  initializeY(n,y,dt,oldfy,vy,fy);
  #pragma omp single
  initializeZ(n,z,dt,oldfz,vz,fz);
  }   
   
  #pragma omp parallel for 
  for (i=0;i<n;i++){
    loopcode_sub(i,n,mass,fcon,delta,grav,sep,rball,xball,yball,zball,dt,x,y,z,fx,fy,fz,vx,vy,vz,oldfx,oldfy,oldfz,pe,ke,te);
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

void updateX(double *fx, int n){
  int nx, ny;
  double zero_v[2] __attribute__ ((aligned (64))) = {0.0, 0.0} ;
  for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny+=2){
     __m128d x = _mm_load_pd(zero_v); 
     _mm_store_pd(&fx[nx*n+ny], x);
    }
  } 
}

void updateY(double *fy, int n){
  int nx, ny;
   double zero_v[2] __attribute__ ((aligned (64))) = {0.0, 0.0} ;
   for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny+=2){
      __m128d x = _mm_load_pd(zero_v); 
     _mm_store_pd(&fy[nx*n+ny], x);
    }
  }
}
void updateZ(double *fz, int n, double grav){
  int nx, ny;
  double grav_v[2] __attribute__ ((aligned (64))) = {-grav, -grav} ;
    for (nx=0;nx<n;nx++){
    for (ny=0;ny<n;ny+=2){
     __m128d x = _mm_load_pd(grav_v); 
     _mm_store_pd(&fz[nx*n+ny], x);
    }
  }
}
double eval_pef_sub_level_2(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,		
                double *fx, double *fy, double *fz,int nx, int ny, int block_size
                ){
  double pe,rlen,xdiff,ydiff,zdiff,vmag;
  int   dx, dy;
  pe = 0.0;
  int n_end = ny+block_size;
  for (ny=ny;ny<n_end;ny++){
    for (dx=MAX(nx-delta,0);dx<MIN(nx+delta+1,n);dx++){
      for (dy=MAX(ny-delta,0);dy<MIN(ny+delta+1,n);dy++){
          if (nx!=dx || ny!=dy){
            float n_d[2] __attribute__ ((aligned (16))) = {nx,ny} ;
            float d_d[2] __attribute__ ((aligned (16))) = {dx,dy} ;

            float temp_v[3] __attribute__ ((aligned (64))) = {0.0, 0.0, 0.0} ;
            __m128 s_tmp  = _mm_sub_ps(_mm_load_ps(n_d),_mm_load_ps(d_d));
            _mm_store_ps(temp_v,(_mm_mul_ps(s_tmp,s_tmp)));
            rlen=sqrt(temp_v[0]+temp_v[1]);
            rlen = rlen * sep;


            float d_diff[3] __attribute__ ((aligned (64))) = {(float)(x[dx*n+dy]),(float) (y[dx*n+dy]),(float)(z[dx*n+dy])} ;
            float n_diff[3] __attribute__ ((aligned (64))) = {(float)(x[nx*n+ny]),(float) (y[nx*n+ny]),(float)(z[nx*n+ny])} ;
            __m128 sse_p = _mm_sub_ps(_mm_load_ps(d_diff),_mm_load_ps(n_diff));

            __m128 s  = _mm_mul_ps(sse_p,sse_p);
            _mm_store_ps(temp_v, s);
            vmag=sqrt(temp_v[1]+temp_v[0]+temp_v[2]);
            double d_vmag_rlen = (vmag-rlen);
            pe += fcon*d_vmag_rlen;
            d_vmag_rlen = fcon* d_vmag_rlen /vmag;
            float d_fcon[3] __attribute__ ((aligned (64))) = {(float)(d_vmag_rlen),(float) (d_vmag_rlen),(float)(d_vmag_rlen)} ;
            __m128 sse_fx = _mm_mul_ps(_mm_load_ps(d_fcon),sse_p);
            _mm_store_ps(temp_v, sse_fx);
            fx[nx*n+ny]+=temp_v[0];
            fy[nx*n+ny]+=temp_v[1];
            fz[nx*n+ny]+=temp_v[2];
          }
      }
    }
  }
  return pe;
}
double eval_pef_sub(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,		
                double *fx, double *fy, double *fz,int nx, int block_size, int no_of_threads_for_columns
                ){
  double pe;
  int  ny;
  pe = 0.0;
  int n_end = nx+block_size;
  for (nx=nx;nx<n_end;nx++){
    omp_set_num_threads(no_of_threads_for_columns);
    #pragma omp parallel for
    for (ny=0;ny<n;ny+=block_size){
      pe+=eval_pef_sub_level_2(n,delta,grav,sep,fcon,x,y,z,fx,fy,fz,nx,ny, block_size);
    }
  }
  return pe;
}

double eval_pef(int n, int delta, double grav, double sep,
                double fcon, double *x, double *y, double *z,		
                double *fx, double *fy, double *fz, int no_threads)
{
  double pe;
  int nx, ny, dx, dy;
  pe = 0.0;
  omp_set_num_threads(no_threads);
  #pragma omp parallel 
  {
    updateX(fx,n);
    updateY(fy,n);
    updateZ(fz,n,grav);
  }  
  int no_of_thread_given = no_threads;
  int no_of_threads_for_rows = MAX(1,no_of_thread_given/2);
  int no_of_threads_for_columns = MAX(1,no_of_thread_given - no_of_threads_for_rows);
  int row_block_size = n/no_of_threads_for_rows;
  //row_bloack and coloumn block size is same as shape is n*n
  omp_set_num_threads(no_of_threads_for_rows);
  int extra_index = n-(row_block_size*no_of_threads_for_rows);
  #pragma omp parallel for 
  for(nx=0;nx<n;nx+=row_block_size){
    if(nx==row_block_size && extra_index>0){
	//To account for extra data, when n is cannot be divided without remainders
        pe+=eval_pef_sub(n,delta,grav,sep,fcon,x,y,z,fx,fy,fz,nx,row_block_size+extra_index, no_of_threads_for_columns);
    }else{
        pe+=eval_pef_sub(n,delta,grav,sep,fcon,x,y,z,fx,fy,fz,nx,row_block_size, no_of_threads_for_columns);
    }
  }


  return pe;
} 


