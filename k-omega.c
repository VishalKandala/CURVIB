#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "petscsnes.h"

extern PetscInt immersed, NumberOfBodies, ti, tistart, wallfunction;
extern double find_utau_Cabot(double nu,  double u, double y, double guess, double dpdn);
extern PetscErrorCode MyKSPMonitor1(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy);
extern PetscInt inletprofile, block_number;

//extern PetscErrorCode  DMCompositeGetAccessVecs(DMComposite packer,Vec gvec, Vec lvec[]);
//extern PetscErrorCode  DMCompositeRestoreAccessVecs(DMComposite packer,Vec gvec, Vec lvec[]);
extern PetscInt k_periodic,j_periodic,i_periodic,blkpbc;


double solid=0.1;
double innerblank=7.;

const double C1=10.0, C2=2.0;	// C1=1~10, C2=2~5

double a1=0.31, beta_star = 0.09, sigma = 0.5, sigma_star = 0.5;
double beta1 = 0.075, beta2 = 0.0828;
double alpha1 = 5./9., alpha2 = 0.44;
const double sigma_k1 = 0.85, sigma_k2 = 1.0;
const double sigma_o1 = 0.50, sigma_o2 = 0.856;

#define wall_omega(dist)  ( 6 / beta1 / user->ren / pow (dist, 2.0 ) )
// See Wilcox p.381

extern int rans;	/* 1: Wilcox Low Re; 2: Wilcox High Re, 3: SST Menter */
extern double kappa;
extern PetscInt moveframe,rotateframe;

void Get_alpha_beta_star(UserCtx *user, double K, double O, double *alpha, double *alpha_star, double *beta_star);

double Upwind(double W, double E, double a)
{
  if(a>0) return W;
  else return E;
};

void Get_alpha_beta_star(UserCtx *user, double K, double O, double *alpha, double *alpha_star, double *beta_star)
{
	
  if(rans==1) {
    const double Rb=8.0, Rk = 6.0, Rw=2.7, alpha0_star = beta1/3.0, alpha0 = 0.1;
    double ReT = K/O*user->ren;
    *alpha_star = (alpha0_star + ReT/Rk)/(1.0+ReT/Rk);
    *alpha = 5./9. * ( alpha0 + ReT/Rw) / ( 1.0 + ReT/Rw ) / (*alpha_star);
    *beta_star = 0.09 * ( 5./18. + pow( ReT/Rb, 4.0 ) ) / ( 1.0 + pow( ReT/Rb, 4.0 ) );
  }
  else *alpha = 5./9., *alpha_star=1., *beta_star = 0.09;
};

void Compute_du_i (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
		   double *dudc, double *dvdc, double *dwdc, 
		   double *dude, double *dvde, double *dwde,
		   double *dudz, double *dvdz, double *dwdz)
{        
  /*if( i_periodic && (i==0 || i==mx-2) ) {
   *dudc = ucat[k][j][1].x - ucat[k][j][mx-2].x;
   *dvdc = ucat[k][j][1].y - ucat[k][j][mx-2].y;
   *dwdc = ucat[k][j][1].z - ucat[k][j][mx-2].z;
   }
   else if( ii_periodic && i==0) {
   *dudc = ucat[k][j][i+1].x - ucat[k][j][-2].x;
   *dvdc = ucat[k][j][i+1].y - ucat[k][j][-2].y;
   *dwdc = ucat[k][j][i+1].z - ucat[k][j][-2].z;
   }
   else if( ii_periodic && i==mx-2) {
   *dudc = ucat[k][j][mx+1].x - ucat[k][j][i].x;
   *dvdc = ucat[k][j][mx+1].y - ucat[k][j][i].y;
   *dwdc = ucat[k][j][mx+1].z - ucat[k][j][i].z;
   }
   else */{
     *dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
     *dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
     *dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;
   }
	
   if ((nvert[k][j+1][i]> solid && nvert[k][j+1][i]<innerblank)  || 
       (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
     *dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
     *dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
     *dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
   }
   else if ((nvert[k][j-1][i]> solid && nvert[k][j-1][i]<innerblank) || 
	    (nvert[k][j-1][i+1]> solid && nvert[k][j-1][i+1]<innerblank)) {
     *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
     *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
     *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
   }
   else {
     *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
     *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
     *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
   }

   //        if((j==1 || j==my-2) && user->_this==blkpbc && j_periodic){ 
   if((j==1 || j==my-2) && j_periodic){ 
     *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
     *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
     *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
   }

   if ((nvert[k+1][j][i]> solid && nvert[k+1][j][i]<innerblank)|| 
       (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
     *dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
     *dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
     *dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
   }
   else if ((nvert[k-1][j][i]> solid && nvert[k-1][j][i] <innerblank) || 
	    (nvert[k-1][j][i+1] > solid && nvert[k-1][j][i+1] < innerblank)) {
     *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
     *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
     *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
   }
   else {
     *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
     *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
     *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
   }

   //        if((k==1 || k==mz-2) && user->_this==blkpbc && k_periodic){ 
   if((k==1 || k==mz-2) && k_periodic){ 
       
     *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
     *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
     *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
   }
}

void Compute_du_j (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
		   double *dudc, double *dvdc, double *dwdc, 
		   double *dude, double *dvde, double *dwde,
		   double *dudz, double *dvdz, double *dwdz)
				
{ 
  if ((nvert[k][j][i+1]> solid && nvert[k][j][i+1]<innerblank) || 
      (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
    *dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
    *dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
    *dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
  }
  else if ((nvert[k][j][i-1]> solid && nvert[k][j][i-1]<innerblank) || 
	   (nvert[k][j+1][i-1] > solid && nvert[k][j+1][i-1]<innerblank)) {
    *dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
    *dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
    *dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
  }
  else {
    *dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
    *dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
    *dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
  }
  /*
    if( j_periodic && (j==0 || j==my-2) ) {
    *dude = ucat[k][1][i].x - ucat[k][my-2][i].x;
    *dvde = ucat[k][1][i].y - ucat[k][my-2][i].y;
    *dwde = ucat[k][1][i].z - ucat[k][my-2][i].z;
    }
    else if( jj_periodic && j==0) {
    *dude = ucat[k][j+1][i].x - ucat[k][-2][i].x;
    *dvde = ucat[k][j+1][i].y - ucat[k][-2][i].y;
    *dwde = ucat[k][j+1][i].z - ucat[k][-2][i].z;
    }
    else if( jj_periodic && j==my-2) {
    *dude = ucat[k][my+1][i].x - ucat[k][j][i].x;
    *dvde = ucat[k][my+1][i].y - ucat[k][j][i].y;
    *dwde = ucat[k][my+1][i].z - ucat[k][j][i].z;
    }
    else */{
      *dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
      *dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
      *dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;
    }
	
  if ((nvert[k+1][j][i] > solid && nvert[k+1][j][i]<innerblank) || 
      (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
    *dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
    *dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
    *dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
  }
  else if ((nvert[k-1][j][i]> solid && nvert[k-1][j][i]<innerblank) || 
	   (nvert[k-1][j+1][i]> solid && nvert[k-1][j+1][i]<innerblank)) {
    *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
    *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
    *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
  }
  else {
    *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
    *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
    *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
  }
 
  //         if((k==1||k==mz-2)&& user->_this==blkpbc && k_periodic) {
  if((k==1||k==mz-2) && k_periodic) {
    *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
    *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
    *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
  }


};

void Compute_du_k (int i, int j, int k, int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
		   double *dudc, double *dvdc, double *dwdc, 
		   double *dude, double *dvde, double *dwde,
		   double *dudz, double *dvdz, double *dwdz)
				
{         
  if ((nvert[k][j][i+1]> solid && nvert[k][j][i+1]<innerblank) || 
      (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
    *dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
    *dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
    *dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
  }
  else if ((nvert[k][j][i-1]> solid && nvert[k][j][i-1]<innerblank) || 
	   (nvert[k+1][j][i-1]> solid && nvert[k+1][j][i-1]<innerblank)) {
    *dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
    *dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
    *dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
  }
  else {
    *dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
    *dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
    *dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
  }

  if ((nvert[k][j+1][i]> solid && nvert[k][j+1][i]<innerblank) || 
      (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
    *dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
    *dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
    *dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
  }
  else if ((nvert[k][j-1][i]> solid && nvert[k][j-1][i]<innerblank) || 
	   (nvert[k+1][j-1][i]> solid && nvert[k+1][j-1][i]<innerblank)) {
    *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
    *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
    *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
  }
  else {
    *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
    *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
    *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
  }

  //         if((j==1||j==my-2)&& user->_this==blkpbc && y_periodic) { 
  if((j==1||j==my-2)&& j_periodic) { 
    *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
    *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
    *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
  }

  /*
    if( k_periodic && (k==0 || k==mz-2) ) {
    *dudz = ucat[1][j][i].x - ucat[mz-2][j][i].x;
    *dvdz = ucat[1][j][i].y - ucat[mz-2][j][i].y;
    *dwdz = ucat[1][j][i].z - ucat[mz-2][j][i].z;
    }
    else if( kk_periodic && k==0) {
    *dudz = ucat[1][j][i].x - ucat[-2][j][i].x;
    *dvdz = ucat[1][j][i].y - ucat[-2][j][i].y;
    *dwdz = ucat[1][j][i].z - ucat[-2][j][i].z;
    }
    else if( kk_periodic && k==mz-2) {
    *dudz = ucat[mz+1][j][i].x - ucat[mz-2][j][i].x;
    *dvdz = ucat[mz+1][j][i].y - ucat[mz-2][j][i].y;
    *dwdz = ucat[mz+1][j][i].z - ucat[mz-2][j][i].z;
    }
    else*/ {
      *dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
      *dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
      *dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;
    }
}

void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
			double *dudc, double *dvdc, double *dwdc, 
			double *dude, double *dvde, double *dwde,
			double *dudz, double *dvdz, double *dwdz)

//issues with indices ????

{
  if ((nvert[k][j][i+1]> solid && nvert[k][j][i+1]<innerblank) || (!i_periodic && i==mx-2) ) {
    *dudc = ( ucat[k][j][i].x - ucat[k][j][i-1].x );
    *dvdc = ( ucat[k][j][i].y - ucat[k][j][i-1].y );
    *dwdc = ( ucat[k][j][i].z - ucat[k][j][i-1].z );
  }
      else if ((nvert[k][j][i-1]> solid && nvert[k][j][i-1]<innerblank) || (!i_periodic  && i==1) ) {
    *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i].x );
    *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i].y );
    *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i].z );
  }
  else {
    /*if(i_periodic && i==1) {
     *dudc = ( ucat[k][j][i+1].x - ucat[k][j][mx-2].x ) * 0.5;
     *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][mx-2].y ) * 0.5;
     *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][mx-2].z ) * 0.5;
     }
     else if(i_periodic && i==mx-2) {
     *dudc = ( ucat[k][j][1].x - ucat[k][j][i-1].x ) * 0.5;
     *dvdc = ( ucat[k][j][1].y - ucat[k][j][i-1].y ) * 0.5;
     *dwdc = ( ucat[k][j][1].z - ucat[k][j][i-1].z ) * 0.5;
     }
     else if(ii_periodic && i==1) {
     *dudc = ( ucat[k][j][i+1].x - ucat[k][j][-2].x ) * 0.5;
     *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][-2].y ) * 0.5;
     *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][-2].z ) * 0.5;
     }
     else if(ii_periodic && i==mx-2) {
     *dudc = ( ucat[k][j][mx+1].x - ucat[k][j][i-1].x ) * 0.5;
     *dvdc = ( ucat[k][j][mx+1].y - ucat[k][j][i-1].y ) * 0.5;
     *dwdc = ( ucat[k][j][mx+1].z - ucat[k][j][i-1].z ) * 0.5;
     }
     else*/ {
       *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
       *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
       *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
     }
  }

  if ((nvert[k][j+1][i]> solid && nvert[k][j+1][i]<innerblank) || (!j_periodic && j==my-2) ) {
    *dude = ( ucat[k][j][i].x - ucat[k][j-1][i].x );
    *dvde = ( ucat[k][j][i].y - ucat[k][j-1][i].y );
    *dwde = ( ucat[k][j][i].z - ucat[k][j-1][i].z );
  }
  else if ((nvert[k][j-1][i]> solid && nvert[k][j-1][i]<innerblank) || (!j_periodic &&  j==1) ) {
    *dude = ( ucat[k][j+1][i].x - ucat[k][j][i].x );
    *dvde = ( ucat[k][j+1][i].y - ucat[k][j][i].y );
    *dwde = ( ucat[k][j+1][i].z - ucat[k][j][i].z );
  }
  else {
    /*if(j_periodic && j==1) {
     *dude = ( ucat[k][j+1][i].x - ucat[k][my-2][i].x ) * 0.5;
     *dvde = ( ucat[k][j+1][i].y - ucat[k][my-2][i].y ) * 0.5;
     *dwde = ( ucat[k][j+1][i].z - ucat[k][my-2][i].z ) * 0.5;
     }
     else if(j_periodic && j==my-2) {
     *dude = ( ucat[k][1][i].x - ucat[k][j-1][i].x ) * 0.5;
     *dvde = ( ucat[k][1][i].y - ucat[k][j-1][i].y ) * 0.5;
     *dwde = ( ucat[k][1][i].z - ucat[k][j-1][i].z ) * 0.5;
     }
     else if(jj_periodic && j==1) {
     *dude = ( ucat[k][j+1][i].x - ucat[k][-2][i].x ) * 0.5;
     *dvde = ( ucat[k][j+1][i].y - ucat[k][-2][i].y ) * 0.5;
     *dwde = ( ucat[k][j+1][i].z - ucat[k][-2][i].z ) * 0.5;
     }
     else if(jj_periodic && j==my-2) {
     *dude = ( ucat[k][my+1][i].x - ucat[k][j-1][i].x ) * 0.5;
     *dvde = ( ucat[k][my+1][i].y - ucat[k][j-1][i].y ) * 0.5;
     *dwde = ( ucat[k][my+1][i].z - ucat[k][j-1][i].z ) * 0.5;
     }
     else */{
       *dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
       *dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
       *dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
     }
  }

  if ((nvert[k+1][j][i]> solid && nvert[k+1][j][i]<innerblank) || ( !k_periodic && k==mz-2) ) {
    *dudz = ( ucat[k][j][i].x - ucat[k-1][j][i].x );
    *dvdz = ( ucat[k][j][i].y - ucat[k-1][j][i].y );
    *dwdz = ( ucat[k][j][i].z - ucat[k-1][j][i].z );
  }
  else if ((nvert[k-1][j][i]> solid && nvert[k-1][j][i]<innerblank) || (!k_periodic && k==0) ) {
    *dudz = ( ucat[k+1][j][i].x - ucat[k][j][i].x );
    *dvdz = ( ucat[k+1][j][i].y - ucat[k][j][i].y );
    *dwdz = ( ucat[k+1][j][i].z - ucat[k][j][i].z );
  }
  else {
    /*if(k_periodic && k==1) {
     *dudz = ( ucat[k+1][j][i].x - ucat[mz-2][j][i].x ) * 0.5;
     *dvdz = ( ucat[k+1][j][i].y - ucat[mz-2][j][i].y ) * 0.5;
     *dwdz = ( ucat[k+1][j][i].z - ucat[mz-2][j][i].z ) * 0.5;
     }
     else if(k_periodic && k==mz-2) {
     *dudz = ( ucat[1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
     *dvdz = ( ucat[1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
     *dwdz = ( ucat[1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
     }
     else if(kk_periodic && k==1) {
     *dudz = ( ucat[k+1][j][i].x - ucat[-2][j][i].x ) * 0.5;
     *dvdz = ( ucat[k+1][j][i].y - ucat[-2][j][i].y ) * 0.5;
     *dwdz = ( ucat[k+1][j][i].z - ucat[-2][j][i].z ) * 0.5;
     }
     else if(kk_periodic && k==mz-2) {
     *dudz = ( ucat[mz+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
     *dvdz = ( ucat[mz+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
     *dwdz = ( ucat[mz+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
     }
     else*/ {
       *dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
       *dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
       *dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
     }
  }
}

void Compute_dkdo_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts2 ***komega, PetscReal ***nvert, 
			  double *dkdc, double *dodc,  double *dkde, double *dode,  double *dkdz, double *dodz )
{
  if (nvert[k][j][i+1]> solid && nvert[k][j][i+1]<innerblank) {
    *dkdc = ( komega[k][j][i].x - komega[k][j][i-1].x );
    *dodc = ( komega[k][j][i].y - komega[k][j][i-1].y );
  }
  else if (nvert[k][j][i-1]> solid && nvert[k][j][i-1]<innerblank) {
    *dkdc = ( komega[k][j][i+1].x - komega[k][j][i].x );
    *dodc = ( komega[k][j][i+1].y - komega[k][j][i].y );
  }
  else {
    if(i_periodic && i==1) {
      *dkdc = ( komega[k][j][i+1].x - komega[k][j][mx-2].x ) * 0.5;
      *dodc = ( komega[k][j][i+1].y - komega[k][j][mx-2].y ) * 0.5;
    }
    else if(i_periodic && i==mx-2) {
      *dkdc = ( komega[k][j][1].x - komega[k][j][i-1].x ) * 0.5;
      *dodc = ( komega[k][j][1].y - komega[k][j][i-1].y ) * 0.5;
    }
    /* else if(ii_periodic && i==1) { */
/*       *dkdc = ( komega[k][j][i+1].x - komega[k][j][-2].x ) * 0.5; */
/*       *dodc = ( komega[k][j][i+1].y - komega[k][j][-2].y ) * 0.5; */
/*     } */
   /*  else if(ii_periodic && i==mx-2) { */
/*       *dkdc = ( komega[k][j][mx+1].x - komega[k][j][i-1].x ) * 0.5; */
/*       *dodc = ( komega[k][j][mx+1].y - komega[k][j][i-1].y ) * 0.5; */
/*     } */
    else {
      *dkdc = ( komega[k][j][i+1].x - komega[k][j][i-1].x ) * 0.5;
      *dodc = ( komega[k][j][i+1].y - komega[k][j][i-1].y ) * 0.5;
    }
  }

  if (nvert[k][j+1][i]> solid && nvert[k][j+1][i]<innerblank) {
    *dkde = ( komega[k][j][i].x - komega[k][j-1][i].x );
    *dode = ( komega[k][j][i].y - komega[k][j-1][i].y );
  }
  else if (nvert[k][j-1][i]> solid && nvert[k][j-1][i]<innerblank) {
    *dkde = ( komega[k][j+1][i].x - komega[k][j][i].x );
    *dode = ( komega[k][j+1][i].y - komega[k][j][i].y );
  }
  else {
    if(j_periodic && j==1) {
      *dkde = ( komega[k][j+1][i].x - komega[k][my-2][i].x ) * 0.5;
      *dode = ( komega[k][j+1][i].y - komega[k][my-2][i].y ) * 0.5;
    }
    else if(j_periodic && j==my-2) {
      *dkde = ( komega[k][1][i].x - komega[k][j-1][i].x ) * 0.5;
      *dode = ( komega[k][1][i].y - komega[k][j-1][i].y ) * 0.5;
    }
    /* else if(jj_periodic && j==1) { */
/*       *dkde = ( komega[k][j+1][i].x - komega[k][-2][i].x ) * 0.5; */
/*       *dode = ( komega[k][j+1][i].y - komega[k][-2][i].y ) * 0.5; */
/*     } */
/*     else if(jj_periodic && j==my-2) { */
/*       *dkde = ( komega[k][my+1][i].x - komega[k][j-1][i].x ) * 0.5; */
/*       *dode = ( komega[k][my+1][i].y - komega[k][j-1][i].y ) * 0.5; */
/*     } */
    else {
      *dkde = ( komega[k][j+1][i].x - komega[k][j-1][i].x ) * 0.5;
      *dode = ( komega[k][j+1][i].y - komega[k][j-1][i].y ) * 0.5;
    }
  }

  if (nvert[k+1][j][i]> solid && nvert[k+1][j][i]<innerblank) {
    *dkdz = ( komega[k][j][i].x - komega[k-1][j][i].x );
    *dodz = ( komega[k][j][i].y - komega[k-1][j][i].y );
  }
  else if (nvert[k-1][j][i]> solid && nvert[k-1][j][i]<innerblank) {
    *dkdz = ( komega[k+1][j][i].x - komega[k][j][i].x );
    *dodz = ( komega[k+1][j][i].y - komega[k][j][i].y );
  }
  else {
    if(k_periodic && k==1) {
      *dkdz = ( komega[k+1][j][i].x - komega[mz-2][j][i].x ) * 0.5;
      *dodz = ( komega[k+1][j][i].y - komega[mz-2][j][i].y ) * 0.5;
    }
    else if(k_periodic && k==mz-2) {
      *dkdz = ( komega[1][j][i].x - komega[k-1][j][i].x ) * 0.5;
      *dodz = ( komega[1][j][i].y - komega[k-1][j][i].y ) * 0.5;
    }
   /*  else if(kk_periodic && k==1) { */
/*       *dkdz = ( komega[k+1][j][i].x - komega[-2][j][i].x ) * 0.5; */
/*       *dodz = ( komega[k+1][j][i].y - komega[-2][j][i].y ) * 0.5; */
/*     } */
/*     else if(kk_periodic && k==mz-2) { */
/*       *dkdz = ( komega[mz+1][j][i].x - komega[k-1][j][i].x ) * 0.5; */
/*       *dodz = ( komega[mz+1][j][i].y - komega[k-1][j][i].y ) * 0.5; */
/*     } */
    else {
      *dkdz = ( komega[k+1][j][i].x - komega[k-1][j][i].x ) * 0.5;
      *dodz = ( komega[k+1][j][i].y - komega[k-1][j][i].y ) * 0.5;
    }
  }
}

void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz)
{
  if (nvert[k][j][i+1]> solid && nvert[k][j][i+1]<innerblank) {
    *dkdc = ( K[k][j][i] - K[k][j][i-1] );
  }
  else if (nvert[k][j][i-1]> solid && nvert[k][j][i-1]<innerblank) {
    *dkdc = ( K[k][j][i+1] - K[k][j][i] );
  }
  else {
    /*if(i_periodic && i==1) {
     *dkdc = ( K[k][j][i+1] - K[k][j][mx-2] ) * 0.5;
     }
     else if(i_periodic && i==mx-2) {
     *dkdc = ( K[k][j][1] - K[k][j][i-1] ) * 0.5;
     }
     else if(ii_periodic && i==1) {
     *dkdc = ( K[k][j][i+1] - K[k][j][-2] ) * 0.5;
     }
     else if(ii_periodic && i==mx-2) {
     *dkdc = ( K[k][j][mx+1] - K[k][j][i-1] ) * 0.5;
     }
     else*/ {
       *dkdc = ( K[k][j][i+1] - K[k][j][i-1] ) * 0.5;
     }
  }

  if (nvert[k][j+1][i]> solid && nvert[k][j+1][i]<innerblank) {
    *dkde = ( K[k][j][i] - K[k][j-1][i] );
  }
  else if (nvert[k][j-1][i]> solid && nvert[k][j-1][i]<innerblank) {
    *dkde = ( K[k][j+1][i] - K[k][j][i] );
  }
  else {
    /*if(j_periodic && j==1) {
     *dkde = ( K[k][j+1][i] - K[k][my-2][i] ) * 0.5;
     }
     else if(j_periodic && j==my-2) {
     *dkde = ( K[k][1][i] - K[k][j-1][i] ) * 0.5;
     }
     else if(jj_periodic && j==1) {
     *dkde = ( K[k][j+1][i] - K[k][-2][i] ) * 0.5;
     }
     else if(jj_periodic && j==my-2) {
     *dkde = ( K[k][my+1][i] - K[k][j-1][i] ) * 0.5;
     }
     else*/ {
       *dkde = ( K[k][j+1][i] - K[k][j-1][i] ) * 0.5;
     }
  }

  if (nvert[k+1][j][i]> solid && nvert[k+1][j][i]<innerblank) {
    *dkdz = ( K[k][j][i] - K[k-1][j][i] );
  }
  else if (nvert[k-1][j][i]> solid && nvert[k-1][j][i]<innerblank) {
    *dkdz = ( K[k+1][j][i] - K[k][j][i] );
  }
  else {
    /*if(k_periodic && k==1) {
     *dkdz = ( K[k+1][j][i] - K[mz-2][j][i] ) * 0.5;
     }
     else if(k_periodic && k==mz-2) {
     *dkdz = ( K[1][j][i] - K[k-1][j][i] ) * 0.5;
     }
     else if(kk_periodic && k==1) {
     *dkdz = ( K[k+1][j][i] - K[-2][j][i] ) * 0.5;
     }
     else if(kk_periodic && k==mz-2) {
     *dkdz = ( K[mz+1][j][i] - K[k-1][j][i] ) * 0.5;
     }
     else */{
       *dkdz = ( K[k+1][j][i] - K[k-1][j][i] ) * 0.5;
     }
  }	
}

void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
			double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
			double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz )
{
  *du_dx = (dudc * csi0 + dude * eta0 + dudz * zet0) * ajc;
  *du_dy = (dudc * csi1 + dude * eta1 + dudz * zet1) * ajc;
  *du_dz = (dudc * csi2 + dude * eta2 + dudz * zet2) * ajc;
  *dv_dx = (dvdc * csi0 + dvde * eta0 + dvdz * zet0) * ajc;
  *dv_dy = (dvdc * csi1 + dvde * eta1 + dvdz * zet1) * ajc;
  *dv_dz = (dvdc * csi2 + dvde * eta2 + dvdz * zet2) * ajc;
  *dw_dx = (dwdc * csi0 + dwde * eta0 + dwdz * zet0) * ajc;
  *dw_dy = (dwdc * csi1 + dwde * eta1 + dwdz * zet1) * ajc;	
  *dw_dz = (dwdc * csi2 + dwde * eta2 + dwdz * zet2) * ajc;
};

void Compute_dkdo_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
				double dkdc, double dodc, double dkde, double dode, double dkdz, double dodz,
				double *dk_dx, double *do_dx, double *dk_dy, double *do_dy, double *dk_dz, double *do_dz)
{
  *dk_dx = (dkdc * csi0 + dkde * eta0 + dkdz * zet0) * ajc;
  *dk_dy = (dkdc * csi1 + dkde * eta1 + dkdz * zet1) * ajc;
  *dk_dz = (dkdc * csi2 + dkde * eta2 + dkdz * zet2) * ajc;
  *do_dx = (dodc * csi0 + dode * eta0 + dodz * zet0) * ajc;
  *do_dy = (dodc * csi1 + dode * eta1 + dodz * zet1) * ajc;
  *do_dz = (dodc * csi2 + dode * eta2 + dodz * zet2) * ajc;
};

void Compute_dscalar_dxyz ( double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
			    double dkdc, double dkde, double dkdz, double *dk_dx, double *dk_dy, double *dk_dz)
{
  *dk_dx = (dkdc * csi0 + dkde * eta0 + dkdz * zet0) * ajc;
  *dk_dy = (dkdc * csi1 + dkde * eta1 + dkdz * zet1) * ajc;
  *dk_dz = (dkdc * csi2 + dkde * eta2 + dkdz * zet2) * ajc;
};

void RHS_K_Omega(UserCtx *user, Vec KOmega_RHS)
{
  DM		da = user->da, fda = user->fda, fda2 = user->fda2;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  double nu = 1./user->ren;
  PetscReal	***aj;
	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
  Cmpnts2	***K_Omega, ***K_Omega_o, ***komega_rhs;
  Cmpnts	***ucont, ***ucat, ***vcont;
  Cmpnts	***csi, ***eta, ***zet;
  PetscReal	***nvert, ***distance, ***lf1, ***lnu_t;
 
  Vec Fpk1, Fpk2, Fpk3;
  Vec Fpo1, Fpo2, Fpo3;
	
  PetscReal	***fp_K1, ***fp_Omega1;
  PetscReal	***fp_K2, ***fp_Omega2;
  PetscReal	***fp_K3, ***fp_Omega3;
	
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  PetscReal	***iaj, ***jaj, ***kaj;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
	
  VecDuplicate(user->lP, &Fpk1);
  VecDuplicate(user->lP, &Fpk2);
  VecDuplicate(user->lP, &Fpk3);
  VecDuplicate(user->lP, &Fpo1);
  VecDuplicate(user->lP, &Fpo2);
  VecDuplicate(user->lP, &Fpo3);
	
  VecSet(Fpk1,0);
  VecSet(Fpk2,0);
  VecSet(Fpk3,0);
  VecSet(Fpo1,0);
  VecSet(Fpo2,0);
  VecSet(Fpo3,0);
	
  if(rans==3) DMDAVecGetArray(da, user->lF1, &lf1);
  DMDAVecGetArray(da, user->lNu_t, &lnu_t);
	
  DMDAVecGetArray(da, Fpk1, &fp_K1);
  DMDAVecGetArray(da, Fpk2, &fp_K2);
  DMDAVecGetArray(da, Fpk3, &fp_K3);
  DMDAVecGetArray(da, Fpo1, &fp_Omega1);
  DMDAVecGetArray(da, Fpo2, &fp_Omega2);
  DMDAVecGetArray(da, Fpo3, &fp_Omega3);
	
  if (moveframe || rotateframe) 
    DMDAVecGetArray(fda, user->lVcont, &vcont);

  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
	
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
	
  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lIZet, &izet);
  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lJZet, &jzet);
  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lKZet, &kzet);
	
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->Distance, &distance);
	
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);
	
  DMDAVecGetArray(fda2, user->lK_Omega, &K_Omega);
  DMDAVecGetArray(fda2, user->lK_Omega_o, &K_Omega_o);
  DMDAVecGetArray(fda2, KOmega_RHS, &komega_rhs);
		
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs-1; i<lxe; i++) {

	double csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
	double eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
	double zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
	double ajc = iaj[k][j][i];
		
	double g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	double g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	double g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
		
	double dkdc, dkde, dkdz;
	double dodc, dode, dodz;
	double nu_t, ucon, vcon;
		
/* 	if (moveframe || rotateframe)  */
/* 	  ucon=ucont[k][j][i].x-vcont[k][j][i].x; */
/* 	else */
	  ucon=ucont[k][j][i].x;

	if (moveframe || rotateframe)
	  vcon=-vcont[k][j][i].x;
	else
	  vcon=0.;

	dkdc = K_Omega[k][j][i+1].x - K_Omega[k][j][i].x;
	dodc = K_Omega[k][j][i+1].y - K_Omega[k][j][i].y;

	if ((nvert[k][j+1][i]> 1.1 && nvert[k][j+1][i] < innerblank) || 
	    (nvert[k][j+1][i+1]> 1.1 && nvert[k][j+1][i+1] < innerblank)) {
	  dkde = (K_Omega[k][j  ][i+1].x + K_Omega[k][j  ][i].x - K_Omega[k][j-1][i+1].x - K_Omega[k][j-1][i].x) * 0.5;
	  dode = (K_Omega[k][j  ][i+1].y + K_Omega[k][j  ][i].y - K_Omega[k][j-1][i+1].y - K_Omega[k][j-1][i].y) * 0.5;
	}
	else if  ((nvert[k][j-1][i]> 1.1 && nvert[k][j-1][i] <innerblank) || 
		  (nvert[k][j-1][i+1]> 1.1 && nvert[k][j-1][i+1] < innerblank)) {
	  dkde = (K_Omega[k][j+1][i+1].x + K_Omega[k][j+1][i].x - K_Omega[k][j  ][i+1].x - K_Omega[k][j  ][i].x) * 0.5;
	  dode = (K_Omega[k][j+1][i+1].y + K_Omega[k][j+1][i].y - K_Omega[k][j  ][i+1].y - K_Omega[k][j  ][i].y) * 0.5;
	}
	else {
	  dkde = (K_Omega[k][j+1][i+1].x + K_Omega[k][j+1][i].x - K_Omega[k][j-1][i+1].x - K_Omega[k][j-1][i].x) * 0.25;
	  dode = (K_Omega[k][j+1][i+1].y + K_Omega[k][j+1][i].y - K_Omega[k][j-1][i+1].y - K_Omega[k][j-1][i].y) * 0.25;
	}	  

	if ((nvert[k+1][j][i]> 1.1 && nvert[k+1][j][i]<innerblank) || 
	    (nvert[k+1][j][i+1]> 1.1 && nvert[k+1][j][i+1]<innerblank)) {
	  dkdz = (K_Omega[k  ][j][i+1].x + K_Omega[k  ][j][i].x - K_Omega[k-1][j][i+1].x - K_Omega[k-1][j][i].x) * 0.5;
	  dodz = (K_Omega[k  ][j][i+1].y + K_Omega[k  ][j][i].y - K_Omega[k-1][j][i+1].y - K_Omega[k-1][j][i].y) * 0.5;
	}
	else if ((nvert[k-1][j][i] > 1.1 && nvert[k-1][j][i] < innerblank) || 
		 (nvert[k-1][j][i+1]> 1.1 && nvert[k-1][j][i+1] < innerblank)) {
	  dkdz = (K_Omega[k+1][j][i+1].x + K_Omega[k+1][j][i].x - K_Omega[k  ][j][i+1].x - K_Omega[k  ][j][i].x) * 0.5;
	  dodz = (K_Omega[k+1][j][i+1].y + K_Omega[k+1][j][i].y - K_Omega[k  ][j][i+1].y - K_Omega[k  ][j][i].y) * 0.5;
	}
	else {
	  dkdz = (K_Omega[k+1][j][i+1].x + K_Omega[k+1][j][i].x - K_Omega[k-1][j][i+1].x - K_Omega[k-1][j][i].x) * 0.25;
	  dodz = (K_Omega[k+1][j][i+1].y + K_Omega[k+1][j][i].y - K_Omega[k-1][j][i+1].y - K_Omega[k-1][j][i].y) * 0.25;
	}
		
/* 	double dk_dx = (dkdc * csi0  + dkde * eta0 + dkdz * zet0) * ajc; */
/* 	double dk_dy = (dkdc * csi1  + dkde * eta1 + dkdz * zet1) * ajc; */
/* 	double dk_dz = (dkdc * csi2  + dkde * eta2 + dkdz * zet2) * ajc; */
/* 	double do_dx = (dodc * csi0  + dode * eta0 + dodz * zet0) * ajc; */
/* 	double do_dy = (dodc * csi1  + dode * eta1 + dodz * zet1) * ajc; */
/* 	double do_dz = (dodc * csi2  + dode * eta2 + dodz * zet2) * ajc; */
		
	double wL=1./aj[k][j][i], wR=1./aj[k][j][i+1];
	double Km_o = ( K_Omega_o[k][j][i].x*wL + K_Omega_o[k][j][i+1].x*wR ) / (wL+wR);		Km_o = PetscMax ( Km_o, 0 );
	double Om_o = ( K_Omega_o[k][j][i].y*wL + K_Omega_o[k][j][i+1].y*wR ) / (wL+wR);	Om_o = PetscMax ( Om_o, 1.e-5 );
	double sigma_k=0, sigma_o=0;
		
	if(rans==3) {
	  /*
	    double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	    double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			
	    Compute_du_i (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	    Compute_du_dxyz ( csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
	    &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
	    double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	    double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	    double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	    double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
	    double vorticity = sqrt( pow( dw_dy - dv_dz, 2. ) + pow ( du_dz - dw_dx, 2. ) +  pow( dv_dx - du_dy, 2. ) );
	    double arg2 = PetscMax ( 2*sqrt(Km_o)/beta_star/Om_o/d, 500/user->ren/Om_o/(d*d) );
	    double f2 = tanh ( arg2*arg2 );
	    nu_t = a1* Km_o / PetscMax ( a1* Om_o, f2 * S );

	    double CD = PetscMax ( 2*sigma_o2/Om_o * ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ), 1.e-20 );
	    double arg1 = PetscMin ( PetscMax ( sqrt(Km_o)/beta_star/Om_o/d, 500./user->ren/Om_o/(d*d) ),  4*sigma_o2*Km_o/CD/(d*d) );
	    double f1 = tanh ( pow(arg1, 4 ) );
	  */
	  nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j][i+1]);
	  double f1 = 0.5 * ( lf1[k][j][i] + lf1[k][j][i+1] );
	  sigma_k = f1 * sigma_k1 + (1.0 - f1) * sigma_k2;
	  sigma_o = f1 * sigma_o1 + (1.0 - f1) * sigma_o2;
	}
	else nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j][i+1]);
	//else nu_t = Km_o/Om_o;
		
	nu_t = PetscMax(nu_t, 0);
		
	fp_K1[k][j][i] = -ucon * Upwind ( K_Omega[k][j][i].x, K_Omega[k][j][i+1].x, ucon) 
	  - vcon * Upwind ( K_Omega[k][j][i].x, K_Omega[k][j][i+1].x, vcon);
	fp_Omega1[k][j][i] = -ucon * Upwind ( K_Omega[k][j][i].y, K_Omega[k][j][i+1].y, ucon)
	  - vcon * Upwind ( K_Omega[k][j][i].y, K_Omega[k][j][i+1].y, vcon);
		
		
		
	if(rans==3) {
	  fp_K1[k][j][i] += (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_k * nu_t);
	  fp_Omega1[k][j][i] += (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_o * nu_t);
	  //fp_K1[k][j][i] += ( nu + sigma_k * nu_t ) * (dk_dx*csi0 + dk_dy*csi1 + dk_dz*csi2);
	  //fp_Omega1[k][j][i] += ( nu + sigma_o * nu_t ) * (do_dx*csi0 + do_dy*csi1 + do_dz*csi2);
	}
	else  {
	  fp_K1[k][j][i] += (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_star * nu_t);
	  fp_Omega1[k][j][i] += (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_star * nu_t);
	  //fp_K1[k][j][i] += ( nu + sigma_star * nu_t ) * (dk_dx*csi0 + dk_dy*csi1 + dk_dz*csi2);
	  //fp_Omega1[k][j][i] += ( nu + sigma_star * nu_t ) * (do_dx*csi0 + do_dy*csi1 + do_dz*csi2);
	}
		
		
      }
	
  for (k=lzs; k<lze; k++)
    for (j=lys-1; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
		
	double csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
	double eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
	double zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
	double ajc = jaj[k][j][i];
		
	double g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
	double g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	double g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
		
	double dkdc, dkde, dkdz;
	double dodc, dode, dodz;
	double nu_t, ucon, vcon;		
		
/* 	if (moveframe || rotateframe)  */
/* 	  ucon=ucont[k][j][i].y-vcont[k][j][i].y; */
/* 	else */
	  ucon=ucont[k][j][i].y;

	if (moveframe || rotateframe) 
	  vcon=-vcont[k][j][i].y;
	else
	  vcon=0.;

	//                if(user->_this==blkpbc && j_periodic && j==0) ucon = ucont[k][my-2][i].y;
	if(j_periodic && j==0) ucon = ucont[k][my-2][i].y;

	if ((nvert[k][j][i+1]> 1.1 && nvert[k][j][i+1]<innerblank) || 
	    (nvert[k][j+1][i+1]> 1.1 && nvert[k][j+1][i+1]<innerblank)) {
	  dkdc = (K_Omega[k][j+1][i  ].x + K_Omega[k][j][i  ].x - K_Omega[k][j+1][i-1].x - K_Omega[k][j][i-1].x) * 0.5;
	  dodc = (K_Omega[k][j+1][i  ].y + K_Omega[k][j][i  ].y - K_Omega[k][j+1][i-1].y - K_Omega[k][j][i-1].y) * 0.5;
	}
	else if ((nvert[k][j][i-1]> 1.1 && nvert[k][j][i-1]<innerblank) || 
		 (nvert[k][j+1][i-1]> 1.1 && nvert[k][j+1][i-1]<innerblank)) {
	  dkdc = (K_Omega[k][j+1][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k][j+1][i  ].x - K_Omega[k][j][i  ].x) * 0.5;
	  dodc = (K_Omega[k][j+1][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k][j+1][i  ].y - K_Omega[k][j][i  ].y) * 0.5;
	}
	else {
	  dkdc = (K_Omega[k][j+1][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k][j+1][i-1].x - K_Omega[k][j][i-1].x) * 0.25;
	  dodc = (K_Omega[k][j+1][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k][j+1][i-1].y - K_Omega[k][j][i-1].y) * 0.25;
	}

	dkde = K_Omega[k][j+1][i].x - K_Omega[k][j][i].x;
	dode = K_Omega[k][j+1][i].y - K_Omega[k][j][i].y;

	if ((nvert[k+1][j][i]> 1.1 && nvert[k+1][j][i]<innerblank) || 
	    (nvert[k+1][j+1][i]> 1.1 && nvert[k+1][j+1][i]<innerblank)) {
	  dkdz = (K_Omega[k  ][j+1][i].x + K_Omega[k  ][j][i].x - K_Omega[k-1][j+1][i].x - K_Omega[k-1][j][i].x) * 0.5;
	  dodz = (K_Omega[k  ][j+1][i].y + K_Omega[k  ][j][i].y - K_Omega[k-1][j+1][i].y - K_Omega[k-1][j][i].y) * 0.5;
	}
	else if ((nvert[k-1][j][i]> 1.1 && nvert[k-1][j][i]<innerblank) || 
		 (nvert[k-1][j+1][i] > 1.1 && nvert[k-1][j+1][i] < innerblank)) {
	  dkdz = (K_Omega[k+1][j+1][i].x + K_Omega[k+1][j][i].x - K_Omega[k  ][j+1][i].x - K_Omega[k  ][j][i].x) * 0.5;
	  dodz = (K_Omega[k+1][j+1][i].y + K_Omega[k+1][j][i].y - K_Omega[k  ][j+1][i].y - K_Omega[k  ][j][i].y) * 0.5;
	}
	else {
	  dkdz = (K_Omega[k+1][j+1][i].x + K_Omega[k+1][j][i].x - K_Omega[k-1][j+1][i].x - K_Omega[k-1][j][i].x) * 0.25;
	  dodz = (K_Omega[k+1][j+1][i].y + K_Omega[k+1][j][i].y - K_Omega[k-1][j+1][i].y - K_Omega[k-1][j][i].y) * 0.25;
	}
		
/* 	double dk_dx = (dkdc * csi0  + dkde * eta0 + dkdz * zet0) * ajc; */
/* 	double dk_dy = (dkdc * csi1  + dkde * eta1 + dkdz * zet1) * ajc; */
/* 	double dk_dz = (dkdc * csi2  + dkde * eta2 + dkdz * zet2) * ajc; */
/* 	double do_dx = (dodc * csi0  + dode * eta0 + dodz * zet0) * ajc; */
/* 	double do_dy = (dodc * csi1  + dode * eta1 + dodz * zet1) * ajc; */
/* 	double do_dz = (dodc * csi2  + dode * eta2 + dodz * zet2) * ajc; */
		
	double wL=1./aj[k][j][i], wR=1./aj[k][j+1][i];
	double Km_o = ( K_Omega_o[k][j][i].x*wL + K_Omega_o[k][j+1][i].x*wR ) / (wL+wR);		Km_o = PetscMax ( Km_o, 0 );
	double Om_o = ( K_Omega_o[k][j][i].y*wL + K_Omega_o[k][j+1][i].y*wR ) / (wL+wR);	Om_o = PetscMax ( Om_o, 1.e-5 );
	double sigma_k=0, sigma_o=0;
		
	if(rans==3) {
	  /*
	    double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	    double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			
	    Compute_du_j (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	    Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
	    &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
	    double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	    double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	    double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	    double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
	    double vorticity = sqrt( pow( dw_dy - dv_dz, 2. ) + pow ( du_dz - dw_dx, 2. ) +  pow( dv_dx - du_dy, 2. ) );
	    double arg2 = PetscMax ( 2*sqrt(Km_o)/beta_star/Om_o/d, 500/user->ren/Om_o/(d*d) );
	    double f2 = tanh ( arg2*arg2 );
	    nu_t = a1* Km_o / PetscMax ( a1* Om_o, f2 * S );
			
	    double CD = PetscMax ( 2*sigma_o2/Om_o * ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ), 1.e-20 );
	    double arg1 = PetscMin ( PetscMax ( sqrt(Km_o)/beta_star/Om_o/d, 500./user->ren/Om_o/(d*d) ),  4*sigma_o2*Km_o/CD/(d*d) );
	    double f1 = tanh ( pow(arg1, 4 ) );
	  */
	  nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i]);
	  double f1 = 0.5 * ( lf1[k][j][i] + lf1[k][j+1][i] );
			
	  sigma_k = f1 * sigma_k1 + (1.0 - f1) * sigma_k2;
	  sigma_o = f1 * sigma_o1 + (1.0 - f1) * sigma_o2;
	}
	else nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k][j+1][i]);
	//else nu_t = Km_o/Om_o;
	nu_t = PetscMax(nu_t, 0);
		
	fp_K2[k][j][i] = -ucon * Upwind ( K_Omega[k][j][i].x, K_Omega[k][j+1][i].x, ucon)
	  -vcon * Upwind ( K_Omega[k][j][i].x, K_Omega[k][j+1][i].x, vcon);
	fp_Omega2[k][j][i] = -ucon * Upwind ( K_Omega[k][j][i].y, K_Omega[k][j+1][i].y, ucon)
	  -vcon * Upwind ( K_Omega[k][j][i].y, K_Omega[k][j+1][i].y, vcon);
		
	if(rans==3) {
	  fp_K2[k][j][i] += (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_k * nu_t);
	  fp_Omega2[k][j][i] += (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_o * nu_t);
	  //fp_K2[k][j][i] += ( nu + sigma_k * nu_t ) * (dk_dx*eta0 + dk_dy*eta1 + dk_dz*eta2);
	  //fp_Omega2[k][j][i] += ( nu + sigma_o * nu_t ) * (do_dx*eta0 + do_dy*eta1 + do_dz*eta2);
	}
	else  {
	  fp_K2[k][j][i] += (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_star * nu_t);
	  fp_Omega2[k][j][i] += (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_star * nu_t);
	  //fp_K2[k][j][i] += ( nu + sigma_star * nu_t ) * (dk_dx*eta0 + dk_dy*eta1 + dk_dz*eta2);
	  //fp_Omega2[k][j][i] += ( nu + sigma_star * nu_t ) * (do_dx*eta0 + do_dy*eta1 + do_dz*eta2);
	}
      }
	
  for (k=lzs-1; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
		
	double csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
	double eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
	double zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
	double ajc = kaj[k][j][i];
		
	double g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
	double g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
	double g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
		
	double dkdc, dkde, dkdz;
	double dodc, dode, dodz;
	double nu_t, ucon, vcon;

/* 	if (moveframe || rotateframe)  */
/* 	  ucon=ucont[k][j][i].z-vcont[k][j][i].z; */
/* 	else */
	  ucon=ucont[k][j][i].z;

	if (moveframe || rotateframe) 
	  vcon=-vcont[k][j][i].z;
	else
	  vcon=0.;

	//                if(user->_this==blkpbc && k_periodic && k==0) ucon = ucont[mz-2][j][i].z;
	if(k_periodic && k==0) ucon = ucont[mz-2][j][i].z;
		
	if ((nvert[k][j][i+1]> 1.1 && nvert[k][j][i+1]<innerblank) || 
	    (nvert[k+1][j][i+1]> 1.1 && nvert[k+1][j][i+1]<innerblank)) {
	  dkdc = (K_Omega[k+1][j][i  ].x + K_Omega[k][j][i  ].x - K_Omega[k+1][j][i-1].x - K_Omega[k][j][i-1].x) * 0.5;
	  dodc = (K_Omega[k+1][j][i  ].y + K_Omega[k][j][i  ].y - K_Omega[k+1][j][i-1].y - K_Omega[k][j][i-1].y) * 0.5;
	}
	else if ((nvert[k][j][i-1]> 1.1 && nvert[k][j][i-1] < innerblank) || 
		 (nvert[k+1][j][i-1]> 1.1 && nvert[k+1][j][i-1]<innerblank)) {
	  dkdc = (K_Omega[k+1][j][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k+1][j][i  ].x - K_Omega[k][j][i  ].x) * 0.5;
	  dodc = (K_Omega[k+1][j][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k+1][j][i  ].y - K_Omega[k][j][i  ].y) * 0.5;
	}
	else {
	  dkdc = (K_Omega[k+1][j][i+1].x + K_Omega[k][j][i+1].x - K_Omega[k+1][j][i-1].x - K_Omega[k][j][i-1].x) * 0.25;
	  dodc = (K_Omega[k+1][j][i+1].y + K_Omega[k][j][i+1].y - K_Omega[k+1][j][i-1].y - K_Omega[k][j][i-1].y) * 0.25;
	}

	if ((nvert[k][j+1][i]> 1.1 && nvert[k][j+1][i]<innerblank) || 
	    (nvert[k+1][j+1][i]> 1.1 && nvert[k+1][j+1][i]<innerblank)) {
	  dkde = (K_Omega[k+1][j  ][i].x + K_Omega[k][j  ][i].x - K_Omega[k+1][j-1][i].x - K_Omega[k][j-1][i].x) * 0.5;
	  dode = (K_Omega[k+1][j  ][i].y + K_Omega[k][j  ][i].y - K_Omega[k+1][j-1][i].y - K_Omega[k][j-1][i].y) * 0.5;
	}
	else if ((nvert[k][j-1][i]> 1.1 && nvert[k][j-1][i]<innerblank) || 
		 (nvert[k+1][j-1][i]> 1.1 && nvert[k+1][j-1][i]<innerblank)) {
	  dkde = (K_Omega[k+1][j+1][i].x + K_Omega[k][j+1][i].x - K_Omega[k+1][j  ][i].x - K_Omega[k][j  ][i].x) * 0.5;
	  dode = (K_Omega[k+1][j+1][i].y + K_Omega[k][j+1][i].y - K_Omega[k+1][j  ][i].y - K_Omega[k][j  ][i].y) * 0.5;
	}
	else {
	  dkde = (K_Omega[k+1][j+1][i].x + K_Omega[k][j+1][i].x - K_Omega[k+1][j-1][i].x - K_Omega[k][j-1][i].x) * 0.25;
	  dode = (K_Omega[k+1][j+1][i].y + K_Omega[k][j+1][i].y - K_Omega[k+1][j-1][i].y - K_Omega[k][j-1][i].y) * 0.25;
	}

	dkdz = K_Omega[k+1][j][i].x - K_Omega[k][j][i].x;
	dodz = K_Omega[k+1][j][i].y - K_Omega[k][j][i].y;
		
/* 	double dk_dx = (dkdc * csi0  + dkde * eta0 + dkdz * zet0) * ajc; */
/* 	double dk_dy = (dkdc * csi1  + dkde * eta1 + dkdz * zet1) * ajc; */
/* 	double dk_dz = (dkdc * csi2  + dkde * eta2 + dkdz * zet2) * ajc; */
/* 	double do_dx = (dodc * csi0  + dode * eta0 + dodz * zet0) * ajc; */
/* 	double do_dy = (dodc * csi1  + dode * eta1 + dodz * zet1) * ajc; */
/* 	double do_dz = (dodc * csi2  + dode * eta2 + dodz * zet2) * ajc; */
		
	double wL=1./aj[k][j][i], wR=1./aj[k+1][j][i];
	double Km_o = ( K_Omega_o[k][j][i].x*wL + K_Omega_o[k+1][j][i].x*wR ) / (wL+wR);		Km_o = PetscMax ( Km_o, 0 );
	double Om_o = ( K_Omega_o[k][j][i].y*wL + K_Omega_o[k+1][j][i].y*wR ) / (wL+wR);	Om_o = PetscMax ( Om_o, 1.e-5 );
	double sigma_k=0, sigma_o=0;
		
	if(rans==3) {
	  /*
	    double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	    double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			
	    Compute_du_k (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	    Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
	    &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
			
	    double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	    double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	    double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	    double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
	    double vorticity = sqrt( pow( dw_dy - dv_dz, 2. ) + pow ( du_dz - dw_dx, 2. ) +  pow( dv_dx - du_dy, 2. ) );
	    double arg2 = PetscMax ( 2*sqrt(Km_o)/beta_star/Om_o/d, 500/user->ren/Om_o/(d*d) );
	    double f2 = tanh ( arg2*arg2 );
	    nu_t = a1* Km_o / PetscMax ( a1* Om_o, f2 * S );
			
	    double CD = PetscMax ( 2*sigma_o2/Om_o * ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ), 1.e-20 );
	    double arg1 = PetscMin ( PetscMax ( sqrt(Km_o)/beta_star/Om_o/d, 500./user->ren/Om_o/(d*d) ),  4*sigma_o2*Km_o/CD/(d*d) );
	    double f1 = tanh ( pow(arg1, 4 ) );
	  */
			
	  double f1 = 0.5 * ( lf1[k][j][i] + lf1[k+1][j][i] );
	  nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k+1][j][i]);
			
	  sigma_k = f1 * sigma_k1 + (1.0 - f1) * sigma_k2;
	  sigma_o = f1 * sigma_o1 + (1.0 - f1) * sigma_o2;
	}
	else nu_t = 0.5*(lnu_t[k][j][i]+lnu_t[k+1][j][i]);
	//else nu_t = Km_o/Om_o;
	nu_t = PetscMax(nu_t, 0);
		
	fp_K3[k][j][i] = -ucon * Upwind ( K_Omega[k][j][i].x, K_Omega[k+1][j][i].x, ucon)
	  -vcon * Upwind ( K_Omega[k][j][i].x, K_Omega[k+1][j][i].x, vcon);
	fp_Omega3[k][j][i] = -ucon * Upwind ( K_Omega[k][j][i].y, K_Omega[k+1][j][i].y, ucon)
	  -vcon * Upwind ( K_Omega[k][j][i].y, K_Omega[k+1][j][i].y, vcon);
		
	if(rans==3) {
	  fp_K3[k][j][i] += (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_k * nu_t);
	  fp_Omega3[k][j][i] += (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_o * nu_t);
	  //fp_K3[k][j][i] += ( nu + sigma_k * nu_t ) * (dk_dx*zet0 + dk_dy*zet1 + dk_dz*zet2);
	  //fp_Omega3[k][j][i] += ( nu + sigma_o * nu_t ) * (do_dx*zet0 + do_dy*zet1 + do_dz*zet2);
	}
	else  {
	  fp_K3[k][j][i] += (g11 * dkdc + g21 * dkde + g31 * dkdz) * ajc * (nu + sigma_star * nu_t);
	  fp_Omega3[k][j][i] += (g11 * dodc + g21 * dode + g31 * dodz) * ajc * (nu + sigma_star * nu_t);
	  //fp_K3[k][j][i] += ( nu + sigma_star * nu_t ) * (dk_dx*zet0 + dk_dy*zet1 + dk_dz*zet2);
	  //fp_Omega3[k][j][i] += ( nu + sigma_star * nu_t ) * (do_dx*zet0 + do_dy*zet1 + do_dz*zet2);
	}
      }

	
	
  DMDAVecRestoreArray(da, Fpk1, &fp_K1);
  DMDAVecRestoreArray(da, Fpk2, &fp_K2);
  DMDAVecRestoreArray(da, Fpk3, &fp_K3);
  DMDAVecRestoreArray(da, Fpo1, &fp_Omega1);
  DMDAVecRestoreArray(da, Fpo2, &fp_Omega2);
  DMDAVecRestoreArray(da, Fpo3, &fp_Omega3);
	
  DMLocalToLocalBegin(da, Fpk1, INSERT_VALUES, Fpk1);
  DMLocalToLocalEnd(da, Fpk1, INSERT_VALUES, Fpk1);
	
  DMLocalToLocalBegin(da, Fpk2, INSERT_VALUES, Fpk2);
  DMLocalToLocalEnd(da, Fpk2, INSERT_VALUES, Fpk2);
	
  DMLocalToLocalBegin(da, Fpk3, INSERT_VALUES, Fpk3);
  DMLocalToLocalEnd(da, Fpk3, INSERT_VALUES, Fpk3);
	
  DMLocalToLocalBegin(da, Fpo1, INSERT_VALUES, Fpo1);
  DMLocalToLocalEnd(da, Fpo1, INSERT_VALUES, Fpo1);
	
  DMLocalToLocalBegin(da, Fpo2, INSERT_VALUES, Fpo2);
  DMLocalToLocalEnd(da, Fpo2, INSERT_VALUES, Fpo2);
	
  DMLocalToLocalBegin(da, Fpo3, INSERT_VALUES, Fpo3);
  DMLocalToLocalEnd(da, Fpo3, INSERT_VALUES, Fpo3);
	
  DMDAVecGetArray(da, Fpk1, &fp_K1);
  DMDAVecGetArray(da, Fpk2, &fp_K2);
  DMDAVecGetArray(da, Fpk3, &fp_K3);
  DMDAVecGetArray(da, Fpo1, &fp_Omega1);
  DMDAVecGetArray(da, Fpo2, &fp_Omega2);
  DMDAVecGetArray(da, Fpo3, &fp_Omega3);
	
	
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if ( i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1 ) {
	  komega_rhs[k][j][i].x = komega_rhs[k][j][i].y = 0;
	  continue;
	}
			
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
	double Km = K_Omega[k][j][i].x, Om = K_Omega[k][j][i].y;
	double Km_o = K_Omega_o[k][j][i].x, Om_o = K_Omega_o[k][j][i].y;
		
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double ajc = aj[k][j][i];
		
	double dkdc, dodc, dkde, dode, dkdz, dodz;
	double dk_dx, dk_dy, dk_dz, do_dx,  do_dy, do_dz;
		
	Compute_dkdo_center (i, j, k, mx, my, mz, K_Omega, nvert,  &dkdc,  &dodc,  &dkde,  &dode,  &dkdz,  &dodz );
	Compute_dkdo_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,dkdc, dodc, dkde, dode, dkdz, dodz, &dk_dx, &do_dx, &dk_dy, &do_dy, &dk_dz, &do_dz);
		
	Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
				&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
							
	/*double dk_dx_o, dk_dy_o, dk_dz_o, do_dx_o,  do_dy_o, do_dz_o;
	  Compute_dkdo_center (i, j, k, mx, my, mz, K_Omega_o, nvert,  &dkdc,  &dodc,  &dkde,  &dode,  &dkdz,  &dodz );
	  Compute_dkdo_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,dkdc, dodc, dkde, dode, dkdz, dodz, &dk_dx_o, &do_dx_o, &dk_dy_o, &do_dy_o, &dk_dz_o, &do_dz_o);
	*/
		

	double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
		
	double nu_t, inv_nu_t, f1=1;
		
	if(rans==3) {
	  /*
	    double arg2 = PetscMax ( 2*sqrt(Km_o)/beta_star/Om_o/d, 500/user->ren/Om_o/(d*d) );
	    double f2 = tanh ( arg2*arg2 );
	    nu_t = a1* Km_o / PetscMax ( a1* Om_o, f2 * S );
			
	    double CD = PetscMax ( 2*sigma_o2/Om * ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ), 1.e-10 );	// implicit treatment is more stable !
	    double arg1 = PetscMin ( PetscMax ( sqrt(Km)/beta_star/Om/d, 500./user->ren/Om/(d*d) ),  4*sigma_o2*Km/CD/(d*d) );
	    f1 = tanh ( pow(arg1, 4 ) );
	  */
			
	  nu_t = lnu_t[k][j][i];
	  f1 = lf1[k][j][i];
			
	  if( !i_periodic && (i==0 || i==mx-2 ) ) f1=1;
	  if( !j_periodic && (j==0 || j==my-2 ) ) f1=1;
	  if( !k_periodic && (k==0 || k==mz-2 ) ) f1=1;
	  if( nvert[k][j][i-1]+nvert[k][j][i+1]>0.1 ) f1=1;
	  if( nvert[k][j-1][i]+nvert[k][j+1][i]>0.1 ) f1=1;
	  if( nvert[k-1][j][i]+nvert[k+1][j][i]>0.1 ) f1=1;
	}
	else nu_t = lnu_t[k][j][i];
	//else nu_t = Km_o/Om_o;
		
	if( nu_t<1.e-20 ) inv_nu_t = 0;
	else inv_nu_t=1./nu_t;
		
	nu_t = PetscMax(nu_t, 0);


		
	//lS[k][j][i] = S;
							
	if ( nvert[k][j][i] < 0.1 || nvert[k][j][i]>=innerblank) {
	  double alpha0 = alpha1, beta0 = beta1, alpha_star=1.0;
	  double Production = 0;
			
	  if(rans==1) Get_alpha_beta_star(user, Km_o, Om_o, &alpha0, &alpha_star, &beta_star);
	  double Dissipation = beta_star * Km * Om;
	  Production = nu_t * S * S;
			
			
	  /* K RHS */
	  // f~ucont*K
	  // rhs~f*ajc~U*K*ajc~u*L^2*K/L^3~uK/L
	  komega_rhs[k][j][i].x = ( fp_K1[k][j][i] - fp_K1[k][j][i-1] + fp_K2[k][j][i] - fp_K2[k][j-1][i] + fp_K3[k][j][i] - fp_K3[k-1][j][i] ) * ajc;	// advection, diffusion
			
	  if(rans==3) komega_rhs[k][j][i].x += PetscMax( PetscMin ( Production, 10*beta_star*Km_o*Om_o ), 0 );
	  else komega_rhs[k][j][i].x += Production;
			
	  komega_rhs[k][j][i].x -= Dissipation;
			
	  if(rans==3) {
	    //				assert( f1>=0 && f1<=1.0 );
	    alpha1 = beta1/beta_star - sigma_o1*kappa*kappa/sqrt(beta_star);
	    alpha2 = beta2/beta_star - sigma_o2*kappa*kappa/sqrt(beta_star);
				
	    alpha0 = alpha1*f1 + alpha2*(1.0-f1);
	    beta0 = beta1*f1 + beta2*(1.0-f1);	
	  }

	  /* Omega RHS */
	  komega_rhs[k][j][i].y = ( fp_Omega1[k][j][i] - fp_Omega1[k][j][i-1] + fp_Omega2[k][j][i] - fp_Omega2[k][j-1][i] + fp_Omega3[k][j][i] - fp_Omega3[k-1][j][i] ) * ajc;	// advection, diffusion
			
	  if(rans==3) komega_rhs[k][j][i].y += alpha0 * S * S;
	  else komega_rhs[k][j][i].y += alpha0 * alpha_star * S * S;
			
	  if(rans==3) komega_rhs[k][j][i].y += 2.0*(1.0-f1) * sigma_o2 / Om_o *  ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz );
	  komega_rhs[k][j][i].y -= beta0 * Om * Om;
	}
      }
	
  if(rans==3) DMDAVecRestoreArray(da, user->lF1, &lf1);
  DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
  DMDAVecRestoreArray(da, Fpk1, &fp_K1);
  DMDAVecRestoreArray(da, Fpk2, &fp_K2);
  DMDAVecRestoreArray(da, Fpk3, &fp_K3);
  DMDAVecRestoreArray(da, Fpo1, &fp_Omega1);
  DMDAVecRestoreArray(da, Fpo2, &fp_Omega2);
  DMDAVecRestoreArray(da, Fpo3, &fp_Omega3);
	
  if (moveframe || rotateframe) 
    DMDAVecRestoreArray(fda, user->lVcont, &vcont);	
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(fda, user->lUcat,  &ucat);
	
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
	
  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);
  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);
  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);
	
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->Distance, &distance);
	
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);
	
  DMDAVecRestoreArray(fda2, user->lK_Omega, &K_Omega);
  DMDAVecRestoreArray(fda2, user->lK_Omega_o, &K_Omega_o);
  DMDAVecRestoreArray(fda2, KOmega_RHS, &komega_rhs);
	
  VecDestroy(&Fpk1);
  VecDestroy(&Fpk2);
  VecDestroy(&Fpk3);
  VecDestroy(&Fpo1);
  VecDestroy(&Fpo2);
  VecDestroy(&Fpo3);
};


void K_Omega_IC(UserCtx *user)
{
  DM		da = user->da;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  double nu = 1./user->ren;
	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
  Cmpnts2 ***K_Omega;
  PetscReal	***nvert;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
	
  // BC for K, Omega
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {	// pressure node; cell center
		
	K_Omega[k][j][i].y = C1;
	K_Omega[k][j][i].x = pow ( 10.0, -C2 ) * nu * C1;
	//K_Omega[k][j][i].x = 1;
		
	if( (nvert[k][j][i]>solid && nvert[k][j][i]<innerblank) || i==0 || i==mx-1 
	    || j==0 || j==my-1 || k==0 || k==mz-1 )
	  K_Omega[k][j][i].x = K_Omega[k][j][i].y = 0;
      }
	
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);
	
  DMLocalToGlobalBegin(user->fda2, user->lK_Omega, INSERT_VALUES, user->K_Omega);
  DMLocalToGlobalEnd(user->fda2, user->lK_Omega, INSERT_VALUES, user->K_Omega);
	
};


void K_Omega_BC(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k, ibi;
  double nu = 1./user->ren;
  PetscReal	***aj, ***distance;
	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

	
  Cmpnts2 ***K_Omega;
  Cmpnts	***csi, ***eta, ***zet, ***ucat;
  PetscReal	***nvert, ***ustar;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(fda, user->lUcat, &ucat);
	
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lUstar, &ustar);

  /* 	PetscReal    norm; */
  /* 	VecNorm(user->lUstar, NORM_INFINITY, &norm); */
  /* 	PetscPrintf(PETSC_COMM_WORLD, "Ustar norm %le\n", norm); */

  DMDAVecGetArray(da, user->Distance, &distance);
	
  DMDAVecGetArray(user->fda2, user->lK_Omega, &K_Omega);
	
  // BC for K, Omega
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {	// pressure node; cell center
	// from saved inflow file
	if(inletprofile==100 && k==1) {
	  K_Omega[k-1][j][i] = user->komega_plane[j][i];
	}
	else if ( user->bctype[4] == 5 && k==1 ) {	// inflow
	  K_Omega[k-1][j][i].y = 2*C1 - K_Omega[k][j][i].y;
	  K_Omega[k-1][j][i].x = 2*pow ( 10.0, -C2 ) * nu * C1 - K_Omega[k][j][i].x;
	}
		
	// slip
	if ( user->bctype[0] == 3 && i==1 ) K_Omega[k][j][i-1] = K_Omega[k][j][i];
	if ( user->bctype[1] == 3 && i==mx-2 ) K_Omega[k][j][i+1] = K_Omega[k][j][i];
	if ( user->bctype[2] == 3 && j==1 ) K_Omega[k][j-1][i] = K_Omega[k][j][i];
	if ( user->bctype[3] == 3 && j==my-2 ) K_Omega[k][j+1][i] = K_Omega[k][j][i];
	if ( user->bctype[4] == 3 && k==1 ) K_Omega[k-1][j][i] = K_Omega[k][j][i];
	if ( user->bctype[5] == 3 && k==mz-2 ) K_Omega[k+1][j][i] = K_Omega[k][j][i];
			
	// outflow
	if ( user->bctype[5] == 4 && k==mz-2 ) {
	  K_Omega[k+1][j][i] = K_Omega[k][j][i];
	}
		
	// couette
	if ( user->bctype[3] == 12 && j==my-2 ) {
	  double dist = distance[k][j][i];
	  K_Omega[k][j][i].y = wall_omega(dist);
	  if(j==my-2) K_Omega[k][j+1][i].x = - K_Omega[k][j][i].x;
	}
		
		
	// wall
	if ( user->bctype[0] == 1 && i<=1 ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double dist = distance[k][j][i];
			
	  dist = 0.5/aj[k][j][i]/area;
			
	  K_Omega[k][j][i].y = wall_omega(dist);
	  if(i==1) K_Omega[k][j][i-1].x = - K_Omega[k][j][i].x + 1.e-5;
	  //K_Omega[k][j][i-1].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
	}
	if ( user->bctype[1] == 1 && i>=mx-2 ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			
	  dist = 0.5/aj[k][j][i]/area;
			
	  K_Omega[k][j][i].y = wall_omega(dist);
	  if(i==mx-2) K_Omega[k][j][i+1].x = - K_Omega[k][j][i].x + 1.e-5;
	  //K_Omega[k][j][i+1].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
	}
	if ( user->bctype[2] == 1 && j<=1 ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			
	  dist = 0.5/aj[k][j][i]/area;
			
	  K_Omega[k][j][i].y = wall_omega(dist);
	  if(j==1) K_Omega[k][j-1][i].x = - K_Omega[k][j][i].x + 1.e-5;
	  //K_Omega[k][j-1][i].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
	}
	if ( user->bctype[3] == 1 && j>=my-2 ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double dist = distance[k][j][i];// 0.5/aj[k][j][i]/area;
			
	  dist = 0.5/aj[k][j][i]/area;
			
	  K_Omega[k][j][i].y = wall_omega(dist);
	  if(j==my-2) K_Omega[k][j+1][i].x = - K_Omega[k][j][i].x + 1.e-5;
	  //K_Omega[k][j+1][i].y = 2*10*wall_omega(dist) - K_Omega[k][j][i].y;
	}
				
	// wall function && wall
	if( nvert[k][j][i]<1.1 && ( (user->bctype[0]==-1 && i==1) || 
				    (user->bctype[1]==-1 && i==mx-2) ) ) {
			
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double sb, sc; 
	  Cmpnts Uc;
			
	  sb = 0.5/aj[k][j][i]/area;
			
	  if(i==1) {
	    sc = 2*sb + 0.5/aj[k][j][i+1]/area;
	    Uc = ucat[k][j][i+1];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j][i-1]/area;
	    Uc = ucat[k][j][i-1];
	  }
			
	  double utau = ustar[k][j][i];
			
	  double Kc = utau*utau/sqrt(0.09);
	  double Ob = utau/sqrt(0.09)/(kappa*sb);
			
	  K_Omega[k][j][i].x = Kc;
	  K_Omega[k][j][i].y = Ob;
	}
		
	if( nvert[k][j][i]<1.1 && ( (user->bctype[2]==-1 && j==1) || 
				    (user->bctype[3]==-1 && j==my-2) ) ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc; 
	  Cmpnts Uc;
			
	  sb = 0.5/aj[k][j][i]/area;
			
	  if(j==1) {
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j-1][i]/area;
	    Uc = ucat[k][j-1][i];
	  }
			
	  double utau = ustar[k][j][i];
			
	  double Kc = utau*utau/sqrt(0.09);
	  double Ob = utau/sqrt(0.09)/(kappa*sb);
			
	  K_Omega[k][j][i].x = Kc;
	  K_Omega[k][j][i].y = Ob;
	}
	       
	if ( nvert[k][j][i] > 1.1 && nvert[k][j][i]<innerblank) K_Omega[k][j][i].x = K_Omega[k][j][i].y = 0.;
				
      }

	
  for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
	int flag=0, a=i, b=j, c=k;
		
	if(i_periodic && i==0) a=mx-2, flag=1;
	else if(i_periodic && i==mx-1) a=1, flag=1;
		
	if (j_periodic && j==0) b=my-2, flag=1;
	else if(j_periodic && j==my-1) b=1, flag=1;
		
	if(k_periodic && k==0) c=mz-2, flag=1;
	else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
/* 	if(ii_periodic && i==0) a=-2, flag=1; */
/* 	else if(ii_periodic && i==mx-1) a=mx+1, flag=1; */
		
/* 	if(jj_periodic && j==0) b=-2, flag=1; */
/* 	else if(jj_periodic && j==my-1) b=my+1, flag=1; */
		
/* 	if(kk_periodic && k==0) c=-2, flag=1; */
/* 	else if(kk_periodic && k==mz-1) c=mz+1, flag=1; */
	
	if(flag) K_Omega[k][j][i] = K_Omega[c][b][a];
      }
	
  if(immersed)
    for(ibi=0; ibi<NumberOfBodies; ibi++)
      {
	//		extern IBMNodes *ibm_ptr;
	IBMNodes *ibm = &(user->ibm[ibi]);
		
	IBMListNode *current;
	current = user->ibmlist[ibi].head;
	while (current) {
	  IBMInfo *ibminfo = &current->ibm_intp;
	  int ni = ibminfo->cell;
	  current = current->next;
	  double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
	  int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
	  int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
	  int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
	  i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
	  double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
	  double nfx=ibm->nf_x[ibminfo->cell];
	  double nfy=ibm->nf_y[ibminfo->cell];
	  double nfz=ibm->nf_z[ibminfo->cell];
			
	  double Kc = (K_Omega[kp1][jp1][ip1].x * sk1 + K_Omega[kp2][jp2][ip2].x * sk2 + K_Omega[kp3][jp3][ip3].x * sk3);
			
	  if(wallfunction && ti>0) {
	    double nx = ibm->nf_x[ni], ny = ibm->nf_y[ni], nz = ibm->nf_z[ni];
				
	    double u_c = (ucat[kp1][jp1][ip1].x * sk1 + ucat[kp2][jp2][ip2].x * sk2 + ucat[kp3][jp3][ip3].x * sk3);
	    double v_c = (ucat[kp1][jp1][ip1].y * sk1 + ucat[kp2][jp2][ip2].y * sk2 + ucat[kp3][jp3][ip3].y * sk3);
	    double w_c = (ucat[kp1][jp1][ip1].z * sk1 + ucat[kp2][jp2][ip2].z * sk2 + ucat[kp3][jp3][ip3].z * sk3);
				
	    double un = u_c * nx + v_c * ny + w_c * nz;
	 
	    double utau = ustar[k][j][i];
				
	    const double yplus_min = 5;
	   
	    sb = PetscMax( sb, 1./user->ren * yplus_min / utau );	// prevent the case sb=0
	    double K, Ob;
				
	    // Wilcox pp.108-109
	    K = utau*utau/sqrt(0.09);
	    Ob = utau/sqrt(0.09)/(kappa*sb);
				
	    K_Omega[k][j][i].x = K;
	    K_Omega[k][j][i].y = Ob;
			
	  }
	  else {
	    double u_b = ucat[k][j][i].x, v_b = ucat[k][j][i].y, w_b = ucat[k][j][i].z;
	    double un = u_b * nfx + v_b * nfy + w_b * nfz;
	    double ut = u_b - un * nfx, vt = v_b - un * nfy, wt = w_b - un * nfz;
	    double ut_mag = sqrt ( ut*ut + vt*vt + wt*wt );
				
	    const double yplus_min = 0.3;
	    double utau = sqrt ( 1./user->ren * ( ut_mag / sb ) );
	    sb = PetscMax( sb, 1./user->ren * yplus_min / utau );	// prevent the case sb=0
				
	    K_Omega[k][j][i].y = wall_omega(sb);	
	    K_Omega[k][j][i].x = PetscMax (Kc * sb / sc, 0);
	  }
	};
      }
	
		
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
	
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);
  DMDAVecRestoreArray(da, user->Distance, &distance);
	
  DMDAVecRestoreArray(user->fda2, user->lK_Omega, &K_Omega);
	
  DMLocalToLocalBegin(user->fda2, user->lK_Omega, INSERT_VALUES, user->lK_Omega);
  DMLocalToLocalEnd(user->fda2, user->lK_Omega, INSERT_VALUES, user->lK_Omega);
};

PetscErrorCode FormFunction_K_omega(SNES snes, Vec K_Omega, Vec Rhs, void *ptr)
{
  UserCtx *user = (UserCtx*)ptr;

  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
	
  PetscReal	***nvert;
  Cmpnts2 ***rhs;

  DMDAGetLocalInfo(user->da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
	

  DMGlobalToLocalBegin(user->fda2, K_Omega, INSERT_VALUES, user->lK_Omega);
  DMGlobalToLocalEnd(user->fda2, K_Omega, INSERT_VALUES, user->lK_Omega);
	
  K_Omega_BC(user);
  RHS_K_Omega(user, Rhs);
		
	
  VecAXPY(Rhs, -1/user->dt, K_Omega);
  VecAXPY(Rhs, 1/user->dt, user->K_Omega_o);
	
  DMDAVecGetArray(user->da, user->lNvert, &nvert);
  DMDAVecGetArray(user->fda2, Rhs, &rhs);
	
  for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
	if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1) {
	  rhs[k][j][i].x = rhs[k][j][i].y = 0;
	}
		
	// couette
	if ( user->bctype[3] == 12 && j==my-2 ) rhs[k][j][i].y = 0;
	
	//wall_omega
	if( i<=1 && user->bctype[0]==1 ) rhs[k][j][i].y = 0;
	if( i>=mx-2 && user->bctype[1]==1 ) rhs[k][j][i].y = 0;
		
	if( j<=1 && user->bctype[2]==1 ) rhs[k][j][i].y = 0;
	if( j>=my-2 && user->bctype[3]==1 ) rhs[k][j][i].y = 0;
		
	if( k==1 && user->bctype[4]==1 ) rhs[k][j][i].y = 0;
	if( k==mz-2 && user->bctype[5]==1 ) rhs[k][j][i].y = 0;
		
	// wall function k, omega
	if ( 	( i==1 && user->bctype[0] == -1 ) || ( i==mx-2 && user->bctype[1] == -1 ) ||
		( j==1 && user->bctype[2] == -1 ) || ( j==my-2 && user->bctype[3] == -1 ) ||
		( k==1 && user->bctype[4] == -1 ) || ( k==mz-2 && user->bctype[5] == -1 )  ) {
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	}

      }
  DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
  DMDAVecRestoreArray(user->fda2, Rhs, &rhs);
	
  return(0);
}

int snes_ko_created=0;
Vec r_ko;
Mat J_ko;
SNES snes_two_eq;

void Solve_K_Omega(UserCtx *user)
{
	
  KSP ksp;
  PC pc;
	
	
  double norm;
  Vec Norm;
	
  int bi=0;
  double tol=1.e-6;
	
  if(!snes_ko_created) {
    snes_ko_created=1;
		
    VecDuplicate(user[bi].K_Omega, &r_ko);
    SNESCreate(PETSC_COMM_WORLD,&snes_two_eq);
    SNESSetFunction(snes_two_eq,r_ko,FormFunction_K_omega,(void *)&user[bi]);
    MatCreateSNESMF(snes_two_eq, &J_ko);
    SNESSetJacobian(snes_two_eq,J_ko,J_ko,MatMFFDComputeJacobian,(void *)&user[bi]);
    SNESSetType(snes_two_eq, SNESNEWTONTR);			//SNESTR,SNESLS	
    SNESSetMaxLinearSolveFailures(snes_two_eq,10000);
    SNESSetMaxNonlinearStepFailures(snes_two_eq,10000);		
    SNESKSPSetUseEW(snes_two_eq, PETSC_TRUE);
    SNESKSPSetParametersEW(snes_two_eq,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    SNESSetTolerances(snes_two_eq,PETSC_DEFAULT,tol,PETSC_DEFAULT,50,50000);
			
    SNESGetKSP(snes_two_eq, &ksp);
    KSPSetType(ksp, KSPGMRES);
    //KSPGMRESSetPreAllocateVectors(ksp);
		
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCNONE);
		
    int maxits=10000;	double rtol=tol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;
    KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
  }
	
  extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy);
  SNESMonitorSet(snes_two_eq,MySNESMonitor,PETSC_NULL,PETSC_NULL);
	
  PetscPrintf(PETSC_COMM_WORLD, "\nSolving K-Omega 0...\n");
	
  SNESSolve(snes_two_eq, PETSC_NULL, user[bi].K_Omega);
	
  // SNESGetFunctionNorm(snes_two_eq, &norm);
  SNESGetFunction(user[bi].snes, &Norm,0,0);
  VecNorm(Norm, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "\nK_Omega SNES residual norm=%.5e\n\n", norm);

  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  DMDAGetLocalInfo(user->da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
	
  Cmpnts2 ***k_omega;
	
  DMDAVecGetArray(user->fda2, user->K_Omega, &k_omega);
  for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
	int flag=0, a=i, b=j, c=k;

	if(i_periodic && i==0) a=mx-2, flag=1;
	else if(i_periodic && i==mx-1) a=1, flag=1;
		
	if(j_periodic && j==0) b=my-2, flag=1;
	else if(j_periodic && j==my-1) b=1, flag=1;
		
	if(k_periodic && k==0) c=mz-2, flag=1;
	else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
	/* if(ii_periodic && i==0) a=-2, flag=1; */
/* 	else if(ii_periodic && i==mx-1) a=mx+1, flag=1; */
		
/* 	if(jj_periodic && j==0) b=-2, flag=1; */
/* 	else if(jj_periodic && j==my-1) b=my+1, flag=1; */
		
/* 	if(kk_periodic && k==0) c=-2, flag=1; */
/* 	else if(kk_periodic && k==mz-1) c=mz+1, flag=1; */
		
	/*
	  if(k_periodic && k==0) c=mz-2, flag=1;
	  else if(k_periodic && k==mz-1) c=1, flag=1;
	  if(j_periodic && j==0) b=my-2, flag=1;
	  else if(j_periodic && j==my-1) b=1, flag=1;
	  if(i_periodic && i==0) a=mx-2, flag=1;
	  else if(i_periodic && i==mx-1) a=1, flag=1;
		
	  if(kk_periodic && k==0) c=-2, flag=1;
	  else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
	*/
	if(flag) k_omega[k][j][i] = k_omega[c][b][a];
      }
  DMDAVecRestoreArray(user->fda2, user->K_Omega, &k_omega);
	
  DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
  DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
  K_Omega_BC(user);
  DMLocalToGlobalBegin(user->fda2, user->lK_Omega, INSERT_VALUES, user->K_Omega);
  DMLocalToGlobalEnd(user->fda2, user->lK_Omega, INSERT_VALUES, user->K_Omega);
};

void K_Omega_Set_Constant(UserCtx *user)
{
  DM		da = user->da;
  DM		fda = user->fda;
  DM		fda2 = user->fda2;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  Cmpnts2 ***K_Omega;
  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal	***nvert, ***nu_t, ***f1, ***aj, ***distance;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
	
  DMGlobalToLocalBegin(fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
  DMGlobalToLocalEnd(fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
	
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lNu_t, &nu_t);
  if(rans==3) DMDAVecGetArray(da, user->lF1, &f1);
  DMDAVecGetArray(da, user->Distance, &distance);
  DMDAVecGetArray(fda2, user->lK_Omega, &K_Omega);
  DMDAVecGetArray(fda, user->lUcat, &ucat);
	
  // BC for K, Omega
  for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {	// pressure node; cell center
	double Km = K_Omega[k][j][i].x, Om = K_Omega[k][j][i].y;
		
	if( (nvert[k][j][i]>0.1  && nvert[k][j][i]<innerblank) || i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 ) {
	  nu_t[k][j][i] = 0;
	  if(rans==3) f1[k][j][i] = 1;
	  continue;
	}
	
	if (Km<0.) {
	  PetscPrintf(PETSC_COMM_SELF, "K < 0 at %d %d %d \n",i,j,k);
	  K_Omega[k][j][i].x=0.;
	}
	
	if (Om<0.){
	  PetscPrintf(PETSC_COMM_SELF, "Omega < 0 at %d %d %d \n",i,j,k);
	  K_Omega[k][j][i].y=1.e-10;
	}

	if(rans==3) {
	  const double a1=0.31, beta_star=0.09;
	  double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	  double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	  double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	  double ajc = aj[k][j][i], d = distance[k][j][i];

	  double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	  double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
			
	  double dkdc, dodc, dkde, dode, dkdz, dodz;
	  double dk_dx, dk_dy, dk_dz, do_dx,  do_dy, do_dz;
			
	  Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	  Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
				&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
	  Compute_dkdo_center (i, j, k, mx, my, mz, K_Omega, nvert,  &dkdc,  &dodc,  &dkde,  &dode,  &dkdz,  &dodz );
	  Compute_dkdo_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,dkdc, dodc, dkde, dode, dkdz, dodz, &dk_dx, &do_dx, &dk_dy, &do_dy, &dk_dz, &do_dz);
			
	  double CD = PetscMax ( 2*sigma_o2/Om * ( dk_dx*do_dx + dk_dy*do_dy + dk_dz*do_dz ), 1.e-10 );
	  double arg1 = PetscMin ( PetscMax ( sqrt(Km)/beta_star/Om/d, 500./user->ren/Om/(d*d) ),  4*sigma_o2*Km/CD/(d*d) );
	  f1[k][j][i] = tanh ( pow(arg1, 4 ) );
			
	  //if(j==my/2 && i==mx/2 && k==mz/2) printf("(%d,%d,%d) d=%f, f1=%f, %f, %f, %f\n", i,j,k, d, f1[k][j][i], sqrt(Km)/beta_star/Om/d, 500./user->ren/Om/(d*d), 4*sigma_o2*Km/CD/(d*d) );
			
	  double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	  double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
	  double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
	  double S = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz) );
	  double vorticity = sqrt( pow( dw_dy - dv_dz, 2. ) + pow ( du_dz - dw_dx, 2. ) +  pow( dv_dx - du_dy, 2. ) );
	  double arg2 = PetscMax ( 2*sqrt(Km)/beta_star/Om/d, 500/user->ren/Om/(d*d) );
	  double f2 = tanh ( arg2*arg2 );
			
	  nu_t[k][j][i] = a1* Km / PetscMax ( a1* Om, f2 * S );
	}
	else {
	  double alpha0, alpha_star=1, beta_star;
	  if(rans==1) Get_alpha_beta_star(user, Km, Om, &alpha0, &alpha_star, &beta_star);
	  nu_t[k][j][i] = Km/Om * alpha_star;
	}
		
	nu_t[k][j][i] = PetscMax(nu_t[k][j][i], 0);
      }
	
  for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
	int flag=0, a=i, b=j, c=k;
		
	if(i_periodic && i==0) a=mx-2, flag=1;
	else if(i_periodic && i==mx-1) a=1, flag=1;
		
	if(j_periodic && j==0) b=my-2, flag=1;
	else if(j_periodic && j==my-1) b=1, flag=1;
		
	if(k_periodic && k==0) c=mz-2, flag=1;
	else if(k_periodic && k==mz-1) c=1, flag=1;
		
		
	/*		if(ii_periodic && i==0) a=-2, flag=1;
	  else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
	  if(jj_periodic && j==0) b=-2, flag=1;
	  else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
	  if(kk_periodic && k==0) c=-2, flag=1;
	  else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
	*/
	/*
	  if(k_periodic && k==0) c=mz-2, flag=1;
	  else if(k_periodic && k==mz-1) c=1, flag=1;
	  if(j_periodic && j==0) b=my-2, flag=1;
	  else if(j_periodic && j==my-1) b=1, flag=1;
	  if(i_periodic && i==0) a=mx-2, flag=1;
	  else if(i_periodic && i==mx-1) a=1, flag=1;
		
	  if(kk_periodic && k==0) c=-2, flag=1;
	  else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
	*/
	if(flag) {
	  nu_t[k][j][i] = nu_t[c][b][a];
	  if(rans==3) f1[k][j][i] = f1[c][b][a];
	}
      }
	
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lNu_t, &nu_t);
  if(rans==3) DMDAVecRestoreArray(da, user->lF1, &f1);
  DMDAVecRestoreArray(da, user->Distance, &distance);
  DMDAVecRestoreArray(fda2, user->lK_Omega, &K_Omega);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
	
  DMLocalToLocalBegin(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
  DMLocalToLocalEnd(da, user->lNu_t, INSERT_VALUES, user->lNu_t);
  if(rans==3)  {
    DMLocalToLocalBegin(da, user->lF1, INSERT_VALUES, user->lF1);
    DMLocalToLocalEnd(da, user->lF1, INSERT_VALUES, user->lF1);
  }
};

PetscErrorCode Block_Interface_K_Omega(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	hostU;
  Cmpnts2 ***itfc, kocent;
  PetscReal *hostu,***nvert;

  VecScatter tolocalall;
  PetscErrorCode ierr;
    
  /* Put All K-O values in a sequential vector available on all cpus */
  for (bi=0; bi<block_number; bi++) {
    /* hostU is a parallel PETSc vector that will hold vector values 
       in the natural numbering, rather than in the PETSc parallel 
       numbering associated with the DA */
    ierr = DMDACreateNaturalVector(user[bi].fda2, &hostU);
    
    // put the center node velocities in the natural ordering in nhostU 
    DMDAGlobalToNaturalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, hostU);
    DMDAGlobalToNaturalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, hostU);

    /* the velocities at cell centers are saved in user[bi].nhostU
       which is a sequential vector*/   
    VecScatterCreateToAll(hostU, &tolocalall, &(user[bi].nhostU));

    VecScatterBegin(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		    SCATTER_FORWARD);
    VecScatterEnd(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		  SCATTER_FORWARD);
   
    VecScatterDestroy(&tolocalall);
    VecDestroy(&hostU);
  }
  
  /**************************************************************************************************************************/
  /* Interpolate the K-omega on the interface nodes
     from the host nodes.
     itfc is the K-O at the cell corners
     hostU is the K-O at the cell centers of the host block */
  /**************************************************************************************************************************/
  
  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    //    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt    lmx, lmy, lmz;

    VecDuplicate(user[bi].lK_Omega, &(user[bi].lItfcKO));		
    VecSet(user[bi].lItfcKO,0.);

    DMDAVecGetArray(user[bi].fda2, user[bi].lItfcKO, &itfc);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];

      if ((x+y+z<-1.) &&
	  ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {
	itfc[ck][cj][ci].x = 0.;
	itfc[ck][cj][ci].y = 0.;
	//	PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %le DDDD\n", ci,cj ,ck,x+y+z);
      } else
	if (ci>=xs && ci<xe &&
	    cj>=ys && cj<ye &&
	    ck>=zs && ck<ze   ) {

	  VecGetArray(user[hb].nhostU, &hostu);
	  lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	  itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 2]
				* (1-x) * (1-y) * (1-z) + //i,j,k
				hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 2]
				* x     * (1-y) * (1-z) + //i+1,j,k
				hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 2]
				* (1-x) * y     * (1-z) + //i,j+1,k
				hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 2]
				* x     * y     * (1-z) + //i+1,j+1,k
				hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 2]
				* (1-x) * (1-y) * z     + //i,j,k+1
				hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 2]
				* x     * (1-y) * z     + //i+1,j,k+1
				hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 2]
				* (1-x) * y     * z     + //i,j+1,k+1
				hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 2]
				* x     * y     * z); //i+1,j+1,k+1

	  itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*2+1]
				* (1-x) * (1-y) * (1-z) + //i,j,k
				hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*2+1]
				* x     * (1-y) * (1-z) + //i+1,j,k
				hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*2+1]
				* (1-x) * y     * (1-z) + //i,j+1,k
				hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*2+1]
				* x     * y     * (1-z) + //i+1,j+1,k
				hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*2+1]
				* (1-x) * (1-y) * z     + //i,j,k+1
				hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*2+1]
				* x     * (1-y) * z     + //i+1,j,k+1
				hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*2+1]
				* (1-x) * y     * z     + //i,j+1,k+1
				hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*2+1]
				* x     * y     * z); //i+1,j+1,k+1

	  VecRestoreArray(user[hb].nhostU, &hostu);
	}
    } // for itfcnumber

    DMDAVecRestoreArray(user[bi].fda2, user[bi].lItfcKO, &itfc);
    DMLocalToLocalBegin(user[bi].fda2, user[bi].lItfcKO, INSERT_VALUES,user[bi].lItfcKO);
    DMLocalToLocalEnd(user[bi].fda2, user[bi].lItfcKO, INSERT_VALUES,user[bi].lItfcKO);
  } // for bi
  
  /**************************************************************************************************************************/
  /* 
     Create BC on cell centers based on the interpolated values of itfc
     at cell corners
  */
  /**************************************************************************************************************************/

  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    Cmpnts2 ***K_Omega;

    DMDAVecGetArray(user[bi].fda2, user[bi].lItfcKO, &itfc);
    DMDAVecGetArray(user[bi].fda2, user[bi].lK_Omega, &K_Omega);

    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  kocent.x = 0.25 * (itfc[k-1][j  ][i  ].x +
			     itfc[k  ][j-1][i  ].x +
			     itfc[k  ][j  ][i  ].x +
			     itfc[k-1][j-1][i  ].x) ;
	  kocent.y = 0.25 * (itfc[k-1][j  ][i  ].y +
			     itfc[k  ][j-1][i  ].y +
			     itfc[k  ][j  ][i  ].y +
			     itfc[k-1][j-1][i  ].y);
	  
	  K_Omega[k][j][i+1].x = (kocent.x *2./3. + K_Omega[k][j][i+2].x/3.);
	  K_Omega[k][j][i+1].y = (kocent.y *2./3. + K_Omega[k][j][i+2].y/3.);

	  K_Omega[k][j][i].x =  2.*kocent.x - K_Omega[k][j][i+1].x;
	  K_Omega[k][j][i].y =  2.*kocent.y - K_Omega[k][j][i+1].y;

	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  kocent.x = 0.25 * (itfc[k-1][j  ][i  ].x +
			     itfc[k  ][j-1][i  ].x +
			     itfc[k  ][j  ][i  ].x +
			     itfc[k-1][j-1][i  ].x);
	  kocent.y = 0.25 * (itfc[k-1][j  ][i  ].y +
			     itfc[k  ][j-1][i  ].y +
			     itfc[k  ][j  ][i  ].y +
			     itfc[k-1][j-1][i  ].y);

	  K_Omega[k][j][i].x = (kocent.x *2./3. + K_Omega[k][j][i-1].x/3.);
	  K_Omega[k][j][i].y = (kocent.y *2./3. + K_Omega[k][j][i-1].y/3.);

	  K_Omega[k][j][i+1].x =  2.*kocent.x - K_Omega[k][j][i].x;
	  K_Omega[k][j][i+1].y =  2.*kocent.y - K_Omega[k][j][i].y;

	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  kocent.x = 0.25 * (itfc[k-1][j  ][i  ].x +
			     itfc[k  ][j  ][i-1].x +
			     itfc[k  ][j  ][i  ].x +
			     itfc[k-1][j  ][i-1].x);
	  kocent.y = 0.25 * (itfc[k-1][j  ][i  ].y +
			     itfc[k  ][j  ][i-1].y +
			     itfc[k  ][j  ][i  ].y +
			     itfc[k-1][j  ][i-1].y);

	  K_Omega[k][j+1][i].x = (kocent.x *2./3. + K_Omega[k][j+2][i].x/3.);
	  K_Omega[k][j+1][i].y = (kocent.y *2./3. + K_Omega[k][j+2][i].y/3.);

	  K_Omega[k][j  ][i].x =  2.*kocent.x - K_Omega[k][j+1][i].x;
	  K_Omega[k][j  ][i].y =  2.*kocent.y - K_Omega[k][j+1][i].y;
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  kocent.x = 0.25 * (itfc[k-1][j  ][i  ].x +
			     itfc[k  ][j  ][i-1].x +
			     itfc[k  ][j  ][i  ].x +
			     itfc[k-1][j  ][i-1].x);
	  kocent.y = 0.25 * (itfc[k-1][j  ][i  ].y +
			     itfc[k  ][j  ][i-1].y +
			     itfc[k  ][j  ][i  ].y +
			     itfc[k-1][j  ][i-1].y);
	  K_Omega[k][j][i].x = (kocent.x *2./3. + K_Omega[k][j-1][i].x/3.);
	  K_Omega[k][j][i].y = (kocent.y *2./3. + K_Omega[k][j-1][i].y/3.);

	  K_Omega[k][j+1][i].x =  2.*kocent.x - K_Omega[k][j][i].x;
	  K_Omega[k][j+1][i].y =  2.*kocent.y - K_Omega[k][j][i].y;
	}
      }
    }


    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;//1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  kocent.x = 0.25 * (itfc[k  ][j  ][i  ].x +
			     itfc[k  ][j  ][i-1].x +
			     itfc[k  ][j-1][i  ].x +
			     itfc[k  ][j-1][i-1].x);
	  kocent.y = 0.25 * (itfc[k  ][j  ][i  ].y +
			     itfc[k  ][j  ][i-1].y +
			     itfc[k  ][j-1][i  ].y +
			     itfc[k  ][j-1][i-1].y);
	  
	  K_Omega[k+1][j][i].x = (kocent.x *2./3.+ K_Omega[k+2][j][i].x /3.);
	  K_Omega[k+1][j][i].y = (kocent.y *2./3.+ K_Omega[k+2][j][i].y /3.);

	  K_Omega[k][j][i].x =  2.*kocent.x - K_Omega[k+1][j][i].x;
	  K_Omega[k][j][i].y =  2.*kocent.y - K_Omega[k+1][j][i].y;
	}
      }
    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  kocent.x = 0.25 * (itfc[k  ][j  ][i  ].x +
			     itfc[k  ][j  ][i-1].x +
			     itfc[k  ][j-1][i  ].x +
			     itfc[k  ][j-1][i-1].x);
	  kocent.y = 0.25 * (itfc[k  ][j  ][i  ].y +
			     itfc[k  ][j  ][i-1].y +
			     itfc[k  ][j-1][i  ].y +
			     itfc[k  ][j-1][i-1].y);

	  K_Omega[k][j][i].x = (kocent.x *2./3.+ K_Omega[k-1][j][i].x /3.);
	  K_Omega[k][j][i].y = (kocent.y *2./3.+ K_Omega[k-1][j][i].y /3.);

	  K_Omega[k+1][j][i].x =  2.*kocent.x - K_Omega[k][j][i].x;
	  K_Omega[k+1][j][i].y =  2.*kocent.y - K_Omega[k][j][i].y;
	}
      }
    }


    // This part is for blanking
    PetscInt blank=0;
    PetscOptionsGetInt(PETSC_NULL, "-blk", &blank, PETSC_NULL);
    if(blank && bi==0){
      DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert);

      for (sb=1; sb<block_number; sb++) {
	  
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    for (i=lxs; i<lxe; i++) {
	      
	      if (((int)(nvert[k][j][i]+0.1) < sb*10) &&
		  ((int)(nvert[k][j][i]+0.1) > sb*10-3) ) {

		K_Omega[k][j][i].x = 0.125 * (itfc[k  ][j  ][i  ].x+
					      itfc[k  ][j  ][i-1].x+
					      itfc[k  ][j-1][i  ].x+
					      itfc[k  ][j-1][i-1].x+
					      itfc[k-1][j  ][i  ].x+
					      itfc[k-1][j  ][i-1].x+
					      itfc[k-1][j-1][i  ].x+
					      itfc[k-1][j-1][i-1].x);

		K_Omega[k][j][i].y = 0.125 * (itfc[k  ][j  ][i  ].y+
					      itfc[k  ][j  ][i-1].y+
					      itfc[k  ][j-1][i  ].y+
					      itfc[k  ][j-1][i-1].y+
					      itfc[k-1][j  ][i  ].y+
					      itfc[k-1][j  ][i-1].y+
					      itfc[k-1][j-1][i  ].y+
					      itfc[k-1][j-1][i-1].y);

	      }// if (nvert)
	    }
	  }
	}
      } //for sb
      DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert);

    } // if blank

    DMDAVecRestoreArray(user[bi].fda2, user[bi].lItfcKO, &itfc);
    DMDAVecRestoreArray(user[bi].fda2, user[bi].lK_Omega, &K_Omega);

    DMLocalToGlobalBegin(user[bi].fda2, user[bi].lK_Omega,INSERT_VALUES,user[bi].K_Omega);
    DMLocalToGlobalEnd(user[bi].fda2, user[bi].lK_Omega,INSERT_VALUES,user[bi].K_Omega);
    
    DMGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);
    DMGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega);

  } // bi

  PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation K-Omega  DDDD\n");

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostU);
    VecDestroy(&user[bi].lItfcKO);
  }

  return(0);
}

/* PetscErrorCode FormFunction_K_Omega_DMMG(SNES snes, Vec KOmega, Vec Rhs, void *ptr) */
/* { */
/*   PetscInt     bi; */
/*   Vec          U[block_number], R[block_number]; */
/*   DMMG         dmmg = (DMMG)ptr; */
/*   DMComposite  dm = (DMComposite)dmmg->dm; */
 
/*   UserCtx      *user = (UserCtx*)dmmg->user; */

/*   // Copy the individual global vecs from the composite vec */
/*   DMCompositeGetAccessVecs(dm, KOmega, U); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(U[bi], user[bi].K_Omega); */
    
/*     DMGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega); */
/*     DMGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega); */
/*   } */
/*   DMCompositeRestoreAccessVecs(dm,KOmega, U); */

/*   if (block_number>1)  */
/*     Block_Interface_K_Omega(user); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     K_Omega_BC(&(user[bi])); */
    
/*     VecDuplicate(user[bi].K_Omega, &(user[bi].Rhs)); */
/*     VecSet(user[bi].Rhs,0.); */

/*     RHS_K_Omega(&user[bi],user[bi].Rhs); */

/*     VecAXPY(user[bi].Rhs, -1/user[bi].dt, user[bi].K_Omega); */
/*     VecAXPY(user[bi].Rhs, 1/user[bi].dt, user[bi].K_Omega_o); */

/*     //Set RHS BC */
/*     DMDALocalInfo	info; */
/*     PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information */
/*     PetscInt	mx, my, mz; // Dimensions in three directions */
/*     PetscInt	i, j, k; */
/*     PetscInt	lxs, lxe, lys, lye, lzs, lze; */
    
/*     PetscReal	***nvert; */
/*     Cmpnts2     ***rhs; */
    
/*     DAGetLocalInfo(user[bi].da, &info); */
/*     mx = info.mx; my = info.my; mz = info.mz; */
/*     xs = info.xs; xe = xs + info.xm; */
/*     ys = info.ys; ye = ys + info.ym; */

/*     zs = info.zs; ze = zs + info.zm; */
    
/*     lxs = xs; lxe = xe; */
/*     lys = ys; lye = ye; */
/*     lzs = zs; lze = ze; */
    
/*     if (xs==0) lxs = xs+1; */
/*     if (ys==0) lys = ys+1; */
/*     if (zs==0) lzs = zs+1; */
    
/*     if (xe==mx) lxe = xe-1; */
/*     if (ye==my) lye = ye-1; */
/*     if (ze==mz) lze = ze-1;     */

/*     DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert); */
/*     DMDAVecGetArray(user[bi].fda2, user[bi].Rhs, &rhs); */
    
/*     for (k=zs; k<ze; k++) */
/*       for (j=ys; j<ye; j++) */
/* 	for (i=xs; i<xe; i++) { */
/* 	  if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i]>0.1) { */
/* 	    rhs[k][j][i].x = rhs[k][j][i].y = 0; */
/* 	  } */
	  
/* 	  // couette */
/* 	  if ( user[bi].bctype[3] == 12 && j==my-2 ) rhs[k][j][i].y = 0; */
	  
/* 	  //wall_omega */
/* 	  if( i<=1 && user[bi].bctype[0]==1 ) rhs[k][j][i].y = 0; */
/* 	  if( i>=mx-2 && user[bi].bctype[1]==1 ) rhs[k][j][i].y = 0; */
	  
/* 	  if( j<=1 && user[bi].bctype[2]==1 ) rhs[k][j][i].y = 0; */
/* 	  if( j>=my-2 && user[bi].bctype[3]==1 ) rhs[k][j][i].y = 0; */
	  
/* 	  if( k==1 && user[bi].bctype[4]==1 ) rhs[k][j][i].y = 0; */
/* 	  if( k==mz-2 && user[bi].bctype[5]==1 ) rhs[k][j][i].y = 0; */
	  
/* 	  // wall function k, omega */
/* 	  if ( 	( i==1 && user[bi].bctype[0] == -1 ) || ( i==mx-2 && user[bi].bctype[1] == -1 ) || */
/* 		( j==1 && user[bi].bctype[2] == -1 ) || ( j==my-2 && user[bi].bctype[3] == -1 ) || */
/* 		( k==1 && user[bi].bctype[4] == -1 ) || ( k==mz-2 && user[bi].bctype[5] == -1 )  ) { */
/* 	    rhs[k][j][i].x = 0; */
/* 	    rhs[k][j][i].y = 0; */
/* 	  } */

/* 	  // interface */
/* 	  if ( 	( i==1 && user[bi].bctype[0] == 0 ) || ( i==mx-2 && user[bi].bctype[1] == 0 ) || */
/* 		( j==1 && user[bi].bctype[2] == 0 ) || ( j==my-2 && user[bi].bctype[3] == 0 ) || */
/* 		( k==1 && user[bi].bctype[4] == 0 ) || ( k==mz-2 && user[bi].bctype[5] == 0 )  ) { */
/* 	    rhs[k][j][i].x = 0.; */
/* 	    rhs[k][j][i].y = 0.; */
/* 	  } */

/* 	} */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert); */
/*     DMDAVecRestoreArray(user[bi].fda2, user[bi].Rhs, &rhs); */

/*   } // for bi */
  
/*   DMCompositeGetAccessVecs(dm, Rhs, R); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(user[bi].Rhs, R[bi]); */
/*   } */
/*   DMCompositeRestoreAccessVecs(dm, Rhs, R); */

/*   // Destroy    */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecDestroy(&user[bi].Rhs); */
/*   }  */
/*   return 0; */
/* } */

/* PetscErrorCode FormInitialGuess_K_Omega_DMMG(DMMG dmmg,Vec X) */
/* { */
/*   DMComposite    dm = (DMComposite)dmmg->dm;   */
/*   UserCtx        *user = (UserCtx*)dmmg->user;  */
/*   PetscInt       bi; */
/*   Vec            U[block_number]; */

/*   DMCompositeGetAccessVecs(dm, X, U); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(user[bi].K_Omega,U[bi]); */
/*     /\*     PetscReal    norm; *\/ */
/*     /\*     VecNorm(user[bi].K_Omega, NORM_INFINITY, &norm); *\/ */
/*     /\*     PetscPrintf(PETSC_COMM_WORLD, "K-Omega norm %d %le\n",bi, norm); *\/ */
/*   } */
/*   DMCompositeRestoreAccessVecs(dm, X, U); */

/*   return(0); */
/* } */

/* PetscErrorCode FormFinalSolution_K_Omega_DMMG(DMMG dmmg,Vec X) */
/* { // copy the solution of DMMG on user->K_Omega */
/*   DMComposite    dm = (DMComposite)dmmg->dm; */
/*   Vec            U[block_number]; */
/*   UserCtx        *user = (UserCtx*)dmmg->user;  */
/*   PetscInt       bi; */

/*   DMCompositeGetAccessVecs(dm, X, U); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(U[bi],user[bi].K_Omega); */
/*     DMGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega); */
/*     DMGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega, INSERT_VALUES, user[bi].lK_Omega); */
/*     /\*     PetscReal    norm; *\/ */
/*     /\*     VecNorm(user[bi].K_Omega, NORM_INFINITY, &norm); *\/ */
/*     /\*     PetscPrintf(PETSC_COMM_WORLD, "K-Omega norm %d %le\n",bi, norm); *\/ */
/*   } */
/*   DMCompositeRestoreAccessVecs(dm, X, U); */

/*   // Apply BC */
/*   if (block_number>1)  */
/*     Block_Interface_K_Omega(user); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     K_Omega_BC(&user[bi]); */
/*     DALocalToGlobal(user[bi].fda2, user[bi].lK_Omega, INSERT_VALUES, user[bi].K_Omega); */
/*   } */

/*   return(0); */
/* } */

/* PetscErrorCode Solve_K_Omega_DMMG(UserCtx *user) */
/* { */
/*   PetscInt     bi; */
/*   DMMG         *dmmg;               /\* multilevel grid structure *\/ */
/*   DMComposite  packer; */
/*   //  UserCtx      *user; */

/*   DMCompositeCreate(PETSC_COMM_WORLD,&packer);     */

/*   for (bi=0; bi<block_number; bi++)  */
/*     DMCompositeAddDM(packer,(DM)user[bi].fda2); */

/*   PetscPrintf(PETSC_COMM_WORLD, "k-o packer created!\n");  */
  
/*   DMMGCreate(PETSC_COMM_WORLD,1,(void *)user,&dmmg); */
    
/*   PetscPrintf(PETSC_COMM_WORLD, "k-o dmmg created\n");  */
    
/*   DMMGSetDM(dmmg,(DM)packer); */
/*   DMMGSetISColoringType(dmmg,IS_COLORING_GLOBAL);//dmmg only works with this colring! */

/*   // set up the snes and the nonlinear function */
/*   DMMGSetSNES(dmmg,FormFunction_K_Omega_DMMG,0); */
/*   // transfer the user context */
/*   DMMGSetUser(dmmg,0,user); */
  
/*   // Set the initial solution from the previous time step */
/*   DMMGSetInitialGuess(dmmg,FormInitialGuess_K_Omega_DMMG); */

/*   /\* Set from options */
/*      IMPORTANT: Always use this option -dmmg_jacobian_mf_fd or it will be */
/*      very slow!!! */
/*      Other options: */
/*      -snes_mf */
/*      -snes_type tr */
/*      -snes_max_it 10 */
/*      -snes_monitor */
/*      -snes_ksp_ew */
/*      -snes_ksp_ew_version 3 */
/*   *\/ */
/*   DMMGSetFromOptions(dmmg); */

/*   //  int rank;   */
/*   //  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   PetscReal ts,te,cput; */
/*   PetscGetTime(&ts); */

/*   DMMGSolve(dmmg); */

/*   PetscGetTime(&te); */
/*   cput=te-ts; */
/*   double norm; */
/*   SNESGetFunctionNorm(DMMGGetSNES(dmmg), &norm); */
/*   PetscPrintf(PETSC_COMM_WORLD, "\nK_Omega SNES residual norm=%.5e time %.3e\n\n", norm, cput); */

/*   // Copy the solution of the DMMG onto user */
/*   FormFinalSolution_K_Omega_DMMG(*dmmg,DMMGGetx(dmmg)); */

/*   // Destroy */
/*   DMCompositeDestroy(packer); */
/*   DMMGDestroy(dmmg); */

/*   return(0); */
/* } */
