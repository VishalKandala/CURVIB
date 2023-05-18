#include "variables.h"

extern PetscInt ti;
extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed, invicid;
extern PetscInt TwoD;
extern PetscInt moveframe,rotateframe;
extern PetscInt  les, clark, rans;

Cmpnts MINUS(Cmpnts v1,Cmpnts v2);
Cmpnts AVERAGE4(Cmpnts v1,Cmpnts v2, Cmpnts v3, Cmpnts v5);
Cmpnts Cross(Cmpnts v1,Cmpnts v2);

void Calculate_dxdydz(PetscReal ajc, Cmpnts csi, Cmpnts eta, 
		      Cmpnts zet, double *dx, double *dy, double *dz);

/* Reconstruct Cartesian velocity components at cell centers from
   volume fluxes defined at cell surface centers
   Input: 
   *user
   Ucont --------- Local Vec of surface volume fluxes
   Ucat  --------- Global Vec of Cartesian velocity
*/

PetscErrorCode Contra2Cart(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      gxs, gxe, gys, gye, gzs, gze;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	mat[3][3], det, det0, det1, det2;

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  Vec		Aj = user->lAj;
  Vec           Ucat= user->Ucat, Ucont=user->lUcont;

  Cmpnts	***csi, ***eta, ***zet;
  PetscReal	***aj;
  Cmpnts	***ucont, ***ucat;

  PetscReal	***nvert;
  PetscReal     innerblank=7.;

  PetscReal	q[3]; //local working array
  PetscInt	i, j, k;
 
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  if (lxs==0) lxs++;
  if (lys==0) lys++;
  if (lzs==0) lzs++;
  
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

 
  DMDAVecGetArray(fda, Csi, &csi);
  DMDAVecGetArray(fda, Eta, &eta);
  DMDAVecGetArray(fda, Zet, &zet);
  DMDAVecGetArray(da,  Aj,  &aj);

  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ( nvert[k][j][i] < 0.1 || nvert[k][j][i] > innerblank) {
	  mat[0][0] = 0.5 * (csi[k][j][i-1].x + csi[k][j][i].x);
	  mat[0][1] = 0.5 * (csi[k][j][i-1].y + csi[k][j][i].y);
	  mat[0][2] = 0.5 * (csi[k][j][i-1].z + csi[k][j][i].z);
	         	      
	  mat[1][0] = 0.5 * (eta[k][j-1][i].x + eta[k][j][i].x);
	  mat[1][1] = 0.5 * (eta[k][j-1][i].y + eta[k][j][i].y);
	  mat[1][2] = 0.5 * (eta[k][j-1][i].z + eta[k][j][i].z);
	         	      
	  mat[2][0] = 0.5 * (zet[k-1][j][i].x + zet[k][j][i].x);
	  mat[2][1] = 0.5 * (zet[k-1][j][i].y + zet[k][j][i].y);
	  mat[2][2] = 0.5 * (zet[k-1][j][i].z + zet[k][j][i].z);


	  q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	    
	  
	  q[1] = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	    
	  
	  q[2] = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	    

	  det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

	  det0 = q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
	    q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

	  det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
	    q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

	  det2 = q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
	    q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
	    q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);


	  ucat[k][j][i].x = det0 / det;
	  ucat[k][j][i].y = det1 / det;
	  ucat[k][j][i].z = det2 / det;
	}
      }
    }
  }
   
  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMDAVecRestoreArray(fda, Csi, &csi);
  DMDAVecRestoreArray(fda, Eta, &eta);
  DMDAVecRestoreArray(fda, Zet, &zet);
  DMDAVecRestoreArray(da,  Aj,  &aj);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  return(0);
}

PetscErrorCode CalcFrameVelocity(UserCtx *user,  Cmpnts u_c, 
				  Cmpnts omega_c, Cmpnts a_c)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt	i, j, k;

  Vec           Coor;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***mareacsi, ***mareaeta, ***mareazet;
  Cmpnts        ***vcont,***coor,***cent, ***wcat,r;

  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);

  if (rotateframe) {
  DMDAVecGetArray(fda, user->lMAreaCsi, &mareacsi);
  DMDAVecGetArray(fda, user->lMAreaEta, &mareaeta);
  DMDAVecGetArray(fda, user->lMAreaZet, &mareazet);
  }
  DMDAVecGetArray(fda, user->Vcont, &vcont);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	  vcont[k][j][i].x= 
	     u_c.x*csi[k][j][i].x +
	     u_c.y*csi[k][j][i].y +
	     u_c.z*csi[k][j][i].z ;

	  if (rotateframe) {
	    vcont[k][j][i].x  = (mareacsi[k][j][i].x * omega_c.x +
				 mareacsi[k][j][i].y * omega_c.y +
				 mareacsi[k][j][i].z * omega_c.z );
	  }	  

      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	  vcont[k][j][i].y =
	    u_c.x*eta[k][j][i].x +
	    u_c.y*eta[k][j][i].y +
	    u_c.z*eta[k][j][i].z ;

	  if (rotateframe) {
	    vcont[k][j][i].y  = (mareaeta[k][j][i].x * omega_c.x +
				 mareaeta[k][j][i].y * omega_c.y +
				 mareaeta[k][j][i].z * omega_c.z );
	  }	  


      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	  
	  vcont[k][j][i].z =
	    u_c.x*zet[k][j][i].x +
	    u_c.y*zet[k][j][i].y +
	    u_c.z*zet[k][j][i].z ;

	  if (rotateframe) {
	    vcont[k][j][i].z  = (mareazet[k][j][i].x * omega_c.x +
				 mareazet[k][j][i].y * omega_c.y +
				 mareazet[k][j][i].z * omega_c.z );
	  }	  

      }
    }
  }

  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);

  if (rotateframe) {
  DMDAVecRestoreArray(fda, user->lMAreaCsi, &mareacsi);
  DMDAVecRestoreArray(fda, user->lMAreaEta, &mareaeta);
  DMDAVecRestoreArray(fda, user->lMAreaZet, &mareazet);
  }

  DMDAVecRestoreArray(fda, user->Vcont, &vcont);

  DMGlobalToLocalBegin(fda, user->Vcont, INSERT_VALUES, user->lVcont);
  DMGlobalToLocalEnd(fda, user->Vcont, INSERT_VALUES, user->lVcont);

  if (rotateframe) {
    DMDAVecGetArray(fda, user->Wcat, &wcat);
    DMDAVecGetArray(fda, user->Cent, &cent);

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 
   
      // BC
    /* frame velocity calculation for periodic boundary condition April 2014*/
  
    if (user->bctype[0]==7 && xs==0){
      i=xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 
    if (user->bctype[0]==7 && xe==mx){
      i=xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 

    if (user->bctype[2]==7 && ys==0){
      j = ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 
    if (user->bctype[2]==7 && ye==my){
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 
    if (user->bctype[4]==7 && zs==0){
      k = zs;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 
    if (user->bctype[4]==7 && ze==mz){
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  r = MINUS(cent[k][j][i],a_c);
	  wcat[k][j][i] = Cross(omega_c, r); 		  
	}
      }
    } 

    Cmpnts    p1,p2,p3,p4;
    if (xs==0 && user->bctype[0]!=7 ) {
      i = xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  p1 = coor[k  ][j  ][i];
	  p2 = coor[k-1][j  ][i];
	  p3 = coor[k  ][j-1][i];
	  p4 = coor[k-1][j-1][i];  

	  r=MINUS(AVERAGE4(p1,p2,p3,p4),a_c);	  
	  wcat[k][j][i] = Cross(omega_c, r);

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j][i+1].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j][i+1].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j][i+1].z;
	}
      }
    }
    
    if (xe==mx && user->bctype[0]!=7 ) {
      i = xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  p1 = coor[k  ][j  ][i-1];
	  p2 = coor[k-1][j  ][i-1];
	  p3 = coor[k  ][j-1][i-1];
	  p4 = coor[k-1][j-1][i-1];  

	  r=MINUS(AVERAGE4(p1,p2,p3,p4),a_c);	  
	  wcat[k][j][i] = Cross(omega_c, r);

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j][i-1].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j][i-1].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j][i-1].z;
	}
      }
    }
    
    
    if (ys==0 && user->bctype[2]!=7) {
      j = ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  p1 = coor[k  ][j][i  ];
	  p2 = coor[k-1][j][i  ];
	  p3 = coor[k  ][j][i-1];
	  p4 = coor[k-1][j][i-1];

	  r=MINUS(AVERAGE4(p1,p2,p3,p4),a_c);	  
	  wcat[k][j][i] = Cross(omega_c, r);

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j+1][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j+1][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j+1][i].z;
	}
      }
    }
    
    if (ye==my && user->bctype[2]!=7) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  p1 = coor[k  ][j-1][i  ];
	  p2 = coor[k-1][j-1][i  ];
	  p3 = coor[k  ][j-1][i-1];
	  p4 = coor[k-1][j-1][i-1];

	  r=MINUS(AVERAGE4(p1,p2,p3,p4),a_c);	  
	  wcat[k][j][i] = Cross(omega_c, r);

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j-1][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j-1][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j-1][i].z;
	}
      }
    }
    
    if (zs==0 &&  user->bctype[4]!=7) {
      k = zs;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  p1 = coor[k][j  ][i  ];
	  p2 = coor[k][j  ][i-1];
	  p3 = coor[k][j-1][i  ];
	  p4 = coor[k][j-1][i-1];

	  r=MINUS(AVERAGE4(p1,p2,p3,p4),a_c);	  
	  wcat[k][j][i] = Cross(omega_c, r);

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k+1][j][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k+1][j][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k+1][j][i].z;
	}
      }
    }
    
    if (ze==mz && user->bctype[4]!=7) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  p1 = coor[k-1][j  ][i  ];
	  p2 = coor[k-1][j  ][i-1];
	  p3 = coor[k-1][j-1][i  ];
	  p4 = coor[k-1][j-1][i-1];

	  r=MINUS(AVERAGE4(p1,p2,p3,p4),a_c);	  
	  wcat[k][j][i] = Cross(omega_c, r);

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k-1][j][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k-1][j][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k-1][j][i].z;
	}
      }
    }

    DMDAVecRestoreArray(fda, user->Cent, &cent);
    DMDAVecRestoreArray(fda, user->Wcat, &wcat);        
  } //if rotateframe

  DMDAVecRestoreArray(fda, Coor, &coor);
  
  
  if (rotateframe) {
  DMGlobalToLocalBegin(fda, user->Wcat, INSERT_VALUES, user->lWcat);
  DMGlobalToLocalEnd(fda, user->Wcat, INSERT_VALUES, user->lWcat);
  
  PetscReal wmax, wmin;
  VecMax(user->Wcat, &i, &wmax);
  VecMin(user->Wcat, &i, &wmin);
  PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!! Vec min %le max %le %le %le %le\n",wmin,wmax,omega_c.x, omega_c.y, omega_c.z);
  }
  return(0);
}

PetscErrorCode CalcFrameVelocity_old(UserCtx *user,  Cmpnts u_c, 
				     Cmpnts omega_c, Cmpnts a_c)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt	i, j, k;

  Vec           Coor;
  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***vcont,***coor,w_c,***wcat,***cent;
  PetscReal     rx,ry,rz;

  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);

  DMDAVecGetArray(fda, user->Vcont, &vcont);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {

	  if (rotateframe) {
	    rx = 0.25*(coor[k][j  ][i].x+coor[k-1][j  ][i].x+
		       coor[k][j-1][i].x+coor[k-1][j-1][i].x)-a_c.x;
	    ry = 0.25*(coor[k][j  ][i].y+coor[k-1][j  ][i].y+
		       coor[k][j-1][i].y+coor[k-1][j-1][i].y)-a_c.y;
	    rz = 0.25*(coor[k][j  ][i].z+coor[k-1][j  ][i].z+
		       coor[k][j-1][i].z+coor[k-1][j-1][i].z)-a_c.z;      
	    
	    w_c.x =-( ry*omega_c.z-omega_c.y*rz ) ;
	    w_c.y = ( rx*omega_c.z-omega_c.x*rz ) ;
	    w_c.z =-( rx*omega_c.y-omega_c.x*ry ) ; 	
	  } else {
	    w_c.x =0. ;
	    w_c.y =0. ;
	    w_c.z =0. ; 	
	  }	  

	  vcont[k][j][i].x= 
	     (u_c.x+w_c.x)*csi[k][j][i].x +
	     (u_c.y+w_c.y)*csi[k][j][i].y +
	     (u_c.z+w_c.z)*csi[k][j][i].z ;
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	  if (rotateframe) {	    
	    rx = 0.25*(coor[k][j][i  ].x+coor[k-1][j][i  ].x+
		       coor[k][j][i-1].x+coor[k-1][j][i-1].x)-a_c.x;
	    ry = 0.25*(coor[k][j][i  ].y+coor[k-1][j][i  ].y+
		       coor[k][j][i-1].y+coor[k-1][j][i-1].y)-a_c.y;
	    rz = 0.25*(coor[k][j][i  ].z+coor[k-1][j][i  ].z+
		       coor[k][j][i-1].z+coor[k-1][j][i-1].z)-a_c.z;      
	    
	    w_c.x =-( ry*omega_c.z-omega_c.y*rz ) ;
	    w_c.y = ( rx*omega_c.z-omega_c.x*rz ) ;
	    w_c.z =-( rx*omega_c.y-omega_c.x*ry ) ; 	
	  } else {
	    w_c.x =0. ;
	    w_c.y =0. ;
	    w_c.z =0. ; 	
	  }	  
	  
	  vcont[k][j][i].y =
	     (u_c.x+w_c.x)*eta[k][j][i].x +
	     (u_c.y+w_c.y)*eta[k][j][i].y +
	     (u_c.z+w_c.z)*eta[k][j][i].z ;

      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	  if (rotateframe) {
	    rx = 0.25*(coor[k][j][i  ].x+coor[k][j-1][i  ].x+
		       coor[k][j][i-1].x+coor[k][j-1][i-1].x)-a_c.x;
	    ry = 0.25*(coor[k][j][i  ].y+coor[k][j-1][i  ].y+
		       coor[k][j][i-1].y+coor[k][j-1][i-1].y)-a_c.y;
	    rz = 0.25*(coor[k][j][i  ].z+coor[k][j-1][i  ].z+
		       coor[k][j][i-1].z+coor[k][j-1][i-1].z)-a_c.z;      
	    
	    w_c.x =-( ry*omega_c.z-omega_c.y*rz ) ;
	    w_c.y = ( rx*omega_c.z-omega_c.x*rz ) ;
	    w_c.z =-( rx*omega_c.y-omega_c.x*ry ) ; 	
	  } else {
	    w_c.x =0. ;
	    w_c.y =0. ;
	    w_c.z =0. ; 	
	  }	  
	  
	  vcont[k][j][i].z =
	     (u_c.x+w_c.x)*zet[k][j][i].x +
	     (u_c.y+w_c.y)*zet[k][j][i].y +
	     (u_c.z+w_c.z)*zet[k][j][i].z ;

      }
    }
  }

  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);

  DMDAVecRestoreArray(fda, user->Vcont, &vcont);

  DMGlobalToLocalBegin(fda, user->Vcont, INSERT_VALUES, user->lVcont);
  DMGlobalToLocalEnd(fda, user->Vcont, INSERT_VALUES, user->lVcont);

  if (rotateframe) {
    //    VecSet(user->Wcat,0.);
    DMDAVecGetArray(fda, user->Wcat, &wcat);
    DMDAVecGetArray(fda, user->Cent, &cent);
    //PetscBarrier(PETSC_NULL);
    //    PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!!\n");

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	  rx = cent[k][j][i].x - a_c.x;
	  ry = cent[k][j][i].y - a_c.y;
	  rz = cent[k][j][i].z - a_c.z;      	

	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	
	  
	}
      }
    } 
    
    //    PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!!\n");

    // BC
    if (xs==0) {
      i = xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  rx = 0.25*(coor[k][j  ][i].x+coor[k-1][j  ][i].x+
		     coor[k][j-1][i].x+coor[k-1][j-1][i].x)-a_c.x;
	  ry = 0.25*(coor[k][j  ][i].y+coor[k-1][j  ][i].y+
		     coor[k][j-1][i].y+coor[k-1][j-1][i].y)-a_c.y;
	  rz = 0.25*(coor[k][j  ][i].z+coor[k-1][j  ][i].z+
		     coor[k][j-1][i].z+coor[k-1][j-1][i].z)-a_c.z;      

	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j][i+1].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j][i+1].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j][i+1].z;
	}
      }
    }
    
    if (xe==mx) {
      i = xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  rx = 0.25*(coor[k][j  ][i-1].x+coor[k-1][j  ][i-1].x+
		     coor[k][j-1][i-1].x+coor[k-1][j-1][i-1].x)-a_c.x;
	  ry = 0.25*(coor[k][j  ][i-1].y+coor[k-1][j  ][i-1].y+
		     coor[k][j-1][i-1].y+coor[k-1][j-1][i-1].y)-a_c.y;
	  rz = 0.25*(coor[k][j  ][i-1].z+coor[k-1][j  ][i-1].z+
		     coor[k][j-1][i-1].z+coor[k-1][j-1][i-1].z)-a_c.z;      

	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	


	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j][i-1].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j][i-1].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j][i-1].z;
	}
      }
    }
    
    
    if (ys==0) {
      j = ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  rx = 0.25*(coor[k][j][i  ].x+coor[k-1][j][i  ].x+
		     coor[k][j][i-1].x+coor[k-1][j][i-1].x)-a_c.x;
	  ry = 0.25*(coor[k][j][i  ].y+coor[k-1][j][i  ].y+
		     coor[k][j][i-1].y+coor[k-1][j][i-1].y)-a_c.y;
	  rz = 0.25*(coor[k][j][i  ].z+coor[k-1][j][i  ].z+
		     coor[k][j][i-1].z+coor[k-1][j][i-1].z)-a_c.z;      
	  
	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j+1][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j+1][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j+1][i].z;
	}
      }
    }
    
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  rx = 0.25*(coor[k][j-1][i  ].x+coor[k-1][j-1][i  ].x+
		     coor[k][j-1][i-1].x+coor[k-1][j-1][i-1].x)-a_c.x;
	  ry = 0.25*(coor[k][j-1][i  ].y+coor[k-1][j-1][i  ].y+
		     coor[k][j-1][i-1].y+coor[k-1][j-1][i-1].y)-a_c.y;
	  rz = 0.25*(coor[k][j-1][i  ].z+coor[k-1][j-1][i  ].z+
		     coor[k][j-1][i-1].z+coor[k-1][j-1][i-1].z)-a_c.z;      
	  
	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k][j-1][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k][j-1][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k][j-1][i].z;
	}
      }
    }
    
    if (zs==0) {
      k = zs;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  rx = 0.25*(coor[k][j][i  ].x+coor[k][j-1][i  ].x+
		     coor[k][j][i-1].x+coor[k][j-1][i-1].x)-a_c.x;
	  ry = 0.25*(coor[k][j][i  ].y+coor[k][j-1][i  ].y+
		     coor[k][j][i-1].y+coor[k][j-1][i-1].y)-a_c.y;
	  rz = 0.25*(coor[k][j][i  ].z+coor[k][j-1][i  ].z+
		     coor[k][j][i-1].z+coor[k][j-1][i-1].z)-a_c.z;      
	  
	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k+1][j][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k+1][j][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k+1][j][i].z;
	}
      }
    }
    
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  rx = 0.25*(coor[k-1][j][i  ].x+coor[k-1][j-1][i  ].x+
		     coor[k-1][j][i-1].x+coor[k-1][j-1][i-1].x)-a_c.x;
	  ry = 0.25*(coor[k-1][j][i  ].y+coor[k-1][j-1][i  ].y+
		     coor[k-1][j][i-1].y+coor[k-1][j-1][i-1].y)-a_c.y;
	  rz = 0.25*(coor[k-1][j][i  ].z+coor[k-1][j-1][i  ].z+
		     coor[k-1][j][i-1].z+coor[k-1][j-1][i-1].z)-a_c.z;      
	  //	  PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!!\n");
	  wcat[k][j][i].x =-( ry*omega_c.z-omega_c.y*rz ) ;
	  wcat[k][j][i].y = ( rx*omega_c.z-omega_c.x*rz ) ;
	  wcat[k][j][i].z =-( rx*omega_c.y-omega_c.x*ry ) ; 	

	  wcat[k][j][i].x = 2 * wcat[k][j][i].x - wcat[k-1][j][i].x;
	  wcat[k][j][i].y = 2 * wcat[k][j][i].y - wcat[k-1][j][i].y;
	  wcat[k][j][i].z = 2 * wcat[k][j][i].z - wcat[k-1][j][i].z;
	}
      }
    }

    DMDAVecRestoreArray(fda, user->Cent, &cent);
    DMDAVecRestoreArray(fda, user->Wcat, &wcat);        
  } //if rotateframe

  //  PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!!\n");

  DMDAVecRestoreArray(fda, Coor, &coor);


  DMGlobalToLocalBegin(fda, user->Wcat, INSERT_VALUES, user->lWcat);
  DMGlobalToLocalEnd(fda, user->Wcat, INSERT_VALUES, user->lWcat);

  VecMax(user->Wcat, &i, &rz);
  VecMin(user->Wcat, &i, &rx);
  PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!! Vec min %le max %le %le %le %le\n",rx,rz,omega_c.x, omega_c.y, omega_c.z);

  return(0);
}
PetscErrorCode Convection_MV(UserCtx *user, Vec Ucont, Vec Ucat, 
			      Vec Conv)
{
  
  
  Cmpnts	***ucont, ***ucat;
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  Cmpnts	***fp1, ***fp2, ***fp3;
  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***conv,***vcont,***wcat;
  PetscReal     ***aj;
 
  PetscReal	ucon, up, um;
  PetscReal	vcon, vp, vm;
  PetscReal	coef = 0.125, innerblank=7.;
 
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  PetscReal	***nvert;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);

  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, user->lVcont, &vcont);
  if (rotateframe)
    DMDAVecGetArray(fda, user->lWcat, &wcat);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(fda, Conv,  &conv);

  VecDuplicate(user->lUcont, &Fp1);
  VecDuplicate(user->lUcont, &Fp2);
  VecDuplicate(user->lUcont, &Fp3);

  DMDAVecGetArray(fda, Fp1, &fp1);
  DMDAVecGetArray(fda, Fp2, &fp2);
  DMDAVecGetArray(fda, Fp3, &fp3);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lAj, &aj);

  /* We have two different sets of node: 1. grid node, the physical points
     where grid lines intercross; 2. storage node, where we store variables.
     All node without explicitly specified as "grid node" refers to
     storage node.

     The integer node is defined at cell center while half node refers to
     the actual grid node. (The reason to choose this arrangement is we need
     ghost node, which is half node away from boundaries, to specify boundary
     conditions. By using this storage arrangement, the actual storage need
     is (IM+1) * (JM + 1) * (KM+1) where IM, JM, & KM refer to the number of
     grid nodes along i, j, k directions.)

     DA, the data structure used to define the storage of 3D arrays, is defined
     as mx * my * mz. mx = IM+1, my = JM+1, mz = KM+1.

     Staggered grid arrangement is used in this solver.
     Pressure is stored at interger node (hence the cell center) and volume
     fluxes defined on the center of each surface of a given control volume
     is stored on the cloest upper integer node. */

  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  //Mohsen Feb 12//
  /* First update the computational ghost points velocity for periodic boundary conditions
     just for this subroutine because of Quick scheme for velocity deravatives */
 if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i]=ucat[k][j][i-2];
	  nvert[k][j][i]=nvert[k][j][i-2];
	}
      }
    }
    if (xe==mx){
      i=mx;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i]=ucat[k][j][i+2];
	  nvert[k][j][i]=nvert[k][j][i+2];
	}
      }
    }
 }
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i]=ucat[k][j-2][i];
	  nvert[k][j][i]=nvert[k][j-2][i];
	}
      }
    }
    if (ye==my){
      j=my;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i]=ucat[k][j+2][i];
	  nvert[k][j][i]=nvert[k][j+2][i];
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i]=ucat[k-2][j][i];
	  nvert[k][j][i]=nvert[k-2][j][i];
	}
      }
    }
    if (ze==mz){
      k=mz;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i]=ucat[k+2][j][i];
	  nvert[k][j][i]=nvert[k+2][j][i];
	}
      }
    }
  }
  
  VecSet(Conv, 0.0);
  
  /* Calculating the convective terms on cell centers.
     First calcualte the contribution from i direction
     The flux is evaluated by QUICK scheme */

  

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs-1; i<lxe; i++){
	
	ucon = ( ucont[k][j][i].x )*0.5;
	
	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);
	
	vcon = (vcont[k][j][i].x) * 0.5;
	
	vm = vcon + fabs(vcon);
	vp = vcon - fabs(vcon);
	
	vm = um - vm;
	vp = up - vp;
	
	
	if (TwoD==1) {
	  fp1[k][j][i].x = 0.;
	  fp1[k][j][i].y = 0.;
	  fp1[k][j][i].z = 0.;
	} else if (i>0 && i<mx-2 &&
		   (nvert[k][j][i+1] < 0.1 || nvert[k][j][i+1]>innerblank) &&
		   (nvert[k][j][i-1] < 0.1 || nvert[k][j][i-1]>innerblank)) { // interial nodes
	  if (les) {
	    fp1[k][j][i].x =(ucon - vcon) * (ucat[k][j][i].x +ucat[k][j][i+1].x);
	    fp1[k][j][i].y =(ucon - vcon) * (ucat[k][j][i].y +ucat[k][j][i+1].y);
	    fp1[k][j][i].z =(ucon - vcon) * (ucat[k][j][i].z +ucat[k][j][i+1].z);
	    
	    if (rotateframe) {
	      
	      fp1[k][j][i].x +=ucon * (wcat[k][j][i+1].x +wcat[k][j][i].x);
	      fp1[k][j][i].y +=ucon * (wcat[k][j][i+1].y +wcat[k][j][i].y);
	      fp1[k][j][i].z +=ucon * (wcat[k][j][i+1].z +wcat[k][j][i].z);
	    }
	  } else {
	    fp1[k][j][i].x = (vm)*(coef*(-ucat[k][j][i+2].x-2.*ucat[k][j][i+1].x +3.*ucat[k][j][i  ].x)+ucat[k][j][i+1].x)
	                   + (vp)*(coef*(-ucat[k][j][i-1].x-2.*ucat[k][j][i  ].x +3.*ucat[k][j][i+1].x)+ucat[k][j][i  ].x);
	    fp1[k][j][i].y = (vm)*(coef*(-ucat[k][j][i+2].y-2.*ucat[k][j][i+1].y +3.*ucat[k][j][i  ].y)+ucat[k][j][i+1].y)
	                   + (vp)*(coef*(-ucat[k][j][i-1].y-2.*ucat[k][j][i  ].y +3.*ucat[k][j][i+1].y)+ucat[k][j][i  ].y);
	    fp1[k][j][i].z = (vm)*(coef*(-ucat[k][j][i+2].z-2.*ucat[k][j][i+1].z +3.*ucat[k][j][i  ].z)+ucat[k][j][i+1].z)
	                   + (vp)*(coef*(-ucat[k][j][i-1].z-2.*ucat[k][j][i  ].z +3.*ucat[k][j][i+1].z)+ucat[k][j][i  ].z);
	    if (rotateframe) {
	      
	      fp1[k][j][i].x += (um)*(coef*(-wcat[k][j][i+2].x-2.* wcat[k][j][i+1].x +3.*wcat[k][j][i  ].x)+wcat[k][j][i+1].x) +
	         	        (up)*(coef*(-wcat[k][j][i-1].x-2.* wcat[k][j][i  ].x +3.*wcat[k][j][i+1].x)+wcat[k][j][i  ].x);
	      fp1[k][j][i].y += (um)*(coef*(-wcat[k][j][i+2].y-2.* wcat[k][j][i+1].y +3.*wcat[k][j][i  ].y)+wcat[k][j][i+1].y) +
		                (up)*(coef*(-wcat[k][j][i-1].y-2.* wcat[k][j][i  ].y +3.*wcat[k][j][i+1].y)+wcat[k][j][i  ].y);
	      fp1[k][j][i].z += (um)*(coef*(-wcat[k][j][i+2].z-2.* wcat[k][j][i+1].z +3.*wcat[k][j][i  ].z)+wcat[k][j][i+1].z) +
		                (up)*(coef*(-wcat[k][j][i-1].z-2.* wcat[k][j][i  ].z +3.*wcat[k][j][i+1].z)+wcat[k][j][i  ].z);
	    }
	  }
	}
	else if (i==0 || nvert[k][j][i-1] > 0.1) {
	  if (les && i==0 && user->bctype[0]==7) {
	    fp1[k][j][i].x =(ucon - vcon) * (ucat[k][j][i].x +ucat[k][j][i+1].x);
	    fp1[k][j][i].y =(ucon - vcon) * (ucat[k][j][i].y +ucat[k][j][i+1].y);
	    fp1[k][j][i].z =(ucon - vcon) * (ucat[k][j][i].z +ucat[k][j][i+1].z);
	    
	    if (rotateframe) {
	      
	      fp1[k][j][i].x +=ucon * (wcat[k][j][i+1].x +wcat[k][j][i].x);
	      fp1[k][j][i].y +=ucon * (wcat[k][j][i+1].y +wcat[k][j][i].y);
	      fp1[k][j][i].z +=ucon * (wcat[k][j][i+1].z +wcat[k][j][i].z);
	    }
	  }else{
	    if (user->bctype[0]==7 && i==0){
	      if (nvert[k][j][i-1]<0.1 && nvert[k][j][i+1]<0.1) {//Mohsen Feb 12//
		fp1[k][j][i].x = (vm)*(coef*(-ucat[k][j][i+2].x-2.*ucat[k][j][i+1].x +3.*ucat[k][j][i  ].x)+ucat[k][j][i+1].x)
	                        +(vp)*(coef*(-ucat[k][j][i-1].x-2.*ucat[k][j][i  ].x +3.*ucat[k][j][i+1].x)+ucat[k][j][i  ].x);
		fp1[k][j][i].y = (vm)*(coef*(-ucat[k][j][i+2].y-2.*ucat[k][j][i+1].y +3.*ucat[k][j][i  ].y)+ucat[k][j][i+1].y)
	                        +(vp)*(coef*(-ucat[k][j][i-1].y-2.*ucat[k][j][i  ].y +3.*ucat[k][j][i+1].y)+ucat[k][j][i  ].y);
		fp1[k][j][i].z = (vm)*(coef*(-ucat[k][j][i+2].z-2.*ucat[k][j][i+1].z +3.*ucat[k][j][i  ].z)+ucat[k][j][i+1].z)
		                +(vp)*(coef*(-ucat[k][j][i-1].z-2.*ucat[k][j][i  ].z +3.*ucat[k][j][i+1].z)+ucat[k][j][i  ].z);
	      }
	    }else{
	      fp1[k][j][i].x =(vm)*(coef*(-ucat[k][j][i+2].x-2.*ucat[k][j][i+1].x+3.*ucat[k][j][i  ].x)+ucat[k][j][i+1].x) +
	                      (vp)*(coef*(-ucat[k][j][i  ].x-2.*ucat[k][j][i  ].x+3.*ucat[k][j][i+1].x)+ucat[k][j][i  ].x);
	      fp1[k][j][i].y =(vm)*(coef*(-ucat[k][j][i+2].y-2.*ucat[k][j][i+1].y+3.*ucat[k][j][i  ].y)+ucat[k][j][i+1].y) +
	                      (vp)*(coef*(-ucat[k][j][i  ].y-2.*ucat[k][j][i  ].y+3.*ucat[k][j][i+1].y)+ucat[k][j][i  ].y);
	      fp1[k][j][i].z =(vm)*(coef*(-ucat[k][j][i+2].z-2.*ucat[k][j][i+1].z+3.*ucat[k][j][i  ].z)+ucat[k][j][i+1].z) +
	                      (vp)*(coef*(-ucat[k][j][i  ].z-2.*ucat[k][j][i  ].z+3.*ucat[k][j][i+1].z)+ucat[k][j][i  ].z);
	    }	 
	    if (rotateframe){
		
	      fp1[k][j][i].x +=um*(coef*(-wcat[k][j][i+2].x-2.*wcat[k][j][i+1].x+3.*wcat[k][j][i  ].x)+wcat[k][j][i+1].x) +
	        	       up*(coef*(-wcat[k][j][i  ].x-2.*wcat[k][j][i  ].x+3.*wcat[k][j][i+1].x)+wcat[k][j][i  ].x);
	      fp1[k][j][i].y +=um*(coef*(-wcat[k][j][i+2].y-2.*wcat[k][j][i+1].y+3.*wcat[k][j][i  ].y)+wcat[k][j][i+1].y) +
		               up*(coef*(-wcat[k][j][i  ].y-2.*wcat[k][j][i  ].y+3.*wcat[k][j][i+1].y)+wcat[k][j][i  ].y);
	      fp1[k][j][i].z +=um*(coef*(-wcat[k][j][i+2].z-2.*wcat[k][j][i+1].z+3.*wcat[k][j][i  ].z)+wcat[k][j][i+1].z) +
		               up*(coef*(-wcat[k][j][i  ].z-2.*wcat[k][j][i  ].z+3.*wcat[k][j][i+1].z)+wcat[k][j][i  ].z);
	      
	    }
	  }
	}	
	else if ( i==mx-2 || nvert[k][j][i+1] > 0.1 ) {
	  if (les && user->bctype[0]==7 && i==mx-2 ) {
	    fp1[k][j][i].x =(ucon - vcon) * (ucat[k][j][i].x +ucat[k][j][i+1].x);
	    fp1[k][j][i].y =(ucon - vcon) * (ucat[k][j][i].y +ucat[k][j][i+1].y);
	    fp1[k][j][i].z =(ucon - vcon) * (ucat[k][j][i].z +ucat[k][j][i+1].z);
	    
	    if (rotateframe) {
	      
	      fp1[k][j][i].x +=ucon * (wcat[k][j][i+1].x +wcat[k][j][i].x);
	      fp1[k][j][i].y +=ucon * (wcat[k][j][i+1].y +wcat[k][j][i].y);
	      fp1[k][j][i].z +=ucon * (wcat[k][j][i+1].z +wcat[k][j][i].z);
	    }
	  } else{
	    if (user->bctype[0]==7 && i==mx-2) {
	      if(nvert[k][j][i-1]<0.1 && nvert[k][j][i+1]<0.1){//Mohsen Feb 2012
		
		fp1[k][j][i].x = (vm)*(coef*(-ucat[k][j][i+2].x-2.*ucat[k][j][i+1].x +3.*ucat[k][j][i  ].x)+ucat[k][j][i+1].x)
		               + (vp)*(coef*(-ucat[k][j][i-1].x-2.*ucat[k][j][i  ].x +3.*ucat[k][j][i+1].x)+ucat[k][j][i  ].x);
		fp1[k][j][i].y = (vm)*(coef*(-ucat[k][j][i+2].y-2.*ucat[k][j][i+1].y +3.*ucat[k][j][i  ].y)+ucat[k][j][i+1].y)
	                       + (vp)*(coef*(-ucat[k][j][i-1].y-2.*ucat[k][j][i  ].y +3.*ucat[k][j][i+1].y)+ucat[k][j][i  ].y);
		fp1[k][j][i].z = (vm)*(coef*(-ucat[k][j][i+2].z-2.*ucat[k][j][i+1].z +3.*ucat[k][j][i  ].z)+ucat[k][j][i+1].z)
	                       + (vp)*(coef*(-ucat[k][j][i-1].z-2.*ucat[k][j][i  ].z +3.*ucat[k][j][i+1].z)+ucat[k][j][i  ].z);
	      }
	    }else{
	      fp1[k][j][i].x =(vm)*(coef*(-ucat[k][j][i+1].x-2.*ucat[k][j][i+1].x+3.*ucat[k][j][i  ].x)+ucat[k][j][i+1].x) +
	                      (vp)*(coef*(-ucat[k][j][i-1].x-2.*ucat[k][j][i  ].x+3.*ucat[k][j][i+1].x)+ucat[k][j][i  ].x);
	      fp1[k][j][i].y =(vm)*(coef*(-ucat[k][j][i+1].y-2.*ucat[k][j][i+1].y+3.*ucat[k][j][i  ].y)+ucat[k][j][i+1].y) +
	                      (vp)*(coef*(-ucat[k][j][i-1].y-2.*ucat[k][j][i  ].y+3.*ucat[k][j][i+1].y)+ucat[k][j][i  ].y);
	      fp1[k][j][i].z =(vm)*(coef*(-ucat[k][j][i+1].z-2.*ucat[k][j][i+1].z+3.*ucat[k][j][i  ].z)+ucat[k][j][i+1].z) +
	                      (vp)*(coef*(-ucat[k][j][i-1].z-2.*ucat[k][j][i  ].z+3.*ucat[k][j][i+1].z)+ucat[k][j][i  ].z);
	    }
	    if (rotateframe) {
	      
	      fp1[k][j][i].x +=um*(coef*(-wcat[k][j][i+1].x-2.*wcat[k][j][i+1].x+3.*wcat[k][j][i  ].x)+wcat[k][j][i+1].x) +
		               up*(coef*(-wcat[k][j][i-1].x-2.*wcat[k][j][i  ].x+3.*wcat[k][j][i+1].x)+wcat[k][j][i  ].x);
	      fp1[k][j][i].y +=um*(coef*(-wcat[k][j][i+1].y-2.*wcat[k][j][i+1].y+3.*wcat[k][j][i  ].y)+wcat[k][j][i+1].y) +
	                       up*(coef*(-wcat[k][j][i-1].y-2.*wcat[k][j][i  ].y+3.*wcat[k][j][i+1].y)+wcat[k][j][i  ].y);
	      fp1[k][j][i].z +=um*(coef*(-wcat[k][j][i+1].z-2.*wcat[k][j][i+1].z+3.*wcat[k][j][i  ].z)+wcat[k][j][i+1].z) +
		               up*(coef*(-wcat[k][j][i-1].z-2.*wcat[k][j][i  ].z+3.*wcat[k][j][i+1].z)+wcat[k][j][i  ].z);
	      
	    }
	  }
	}
      }
    }
  }


  /* j direction */
  for (k=lzs; k<lze; k++) {
    for(j=lys-1; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {

	ucon = ( ucont[k][j][i].y )*0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	vcon = (vcont[k][j][i].y) * 0.5;  

	vm = vcon + fabs(vcon);
	vp = vcon - fabs(vcon);
	vm = um - vm;
	vp = up - vp;

	if (TwoD==2) {
	  fp2[k][j][i].x = 0.;
	  fp2[k][j][i].y = 0.;
	  fp2[k][j][i].z = 0.;
	}
	else if (j>0 && j<my-2 &&
		 (nvert[k][j+1][i] < 0.1 || nvert[k][j+1][i] > innerblank) &&
		 (nvert[k][j-1][i] < 0.1 || nvert[k][j-1][i] > innerblank)) {
	  if (les) {
	    fp2[k][j][i].x =(ucon - vcon) * (ucat[k][j][i].x +ucat[k][j+1][i].x);
	    fp2[k][j][i].y =(ucon - vcon) * (ucat[k][j][i].y +ucat[k][j+1][i].y);
	    fp2[k][j][i].z =(ucon - vcon) * (ucat[k][j][i].z +ucat[k][j+1][i].z);
	    
	    if (rotateframe) {
	      
	      fp2[k][j][i].x +=ucon * (wcat[k][j][i].x +wcat[k][j+1][i].x);
	      fp2[k][j][i].y +=ucon * (wcat[k][j][i].y +wcat[k][j+1][i].y);
	      fp2[k][j][i].z +=ucon * (wcat[k][j][i].z +wcat[k][j+1][i].z);
	    }
	  } else {

	    fp2[k][j][i].x = 
	    (vm) * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    (vp) * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j][i  ].x);
	    fp2[k][j][i].y = 
	    (vm) * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    (vp) * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);
	    fp2[k][j][i].z = 				     
	    (vm) * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    (vp) * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);

	    if (rotateframe){
	      fp2[k][j][i].x += 
		um * (coef * (-wcat[k][j+2][i].x -2. * wcat[k][j+1][i].x +3. * wcat[k][j  ][i].x) + wcat[k][j+1][i].x) +		     
		up * (coef * (-wcat[k][j-1][i].x -2. * wcat[k][j  ][i].x +3. * wcat[k][j+1][i].x) + wcat[k][j][i  ].x);
	      fp2[k][j][i].y += 				     
		um * (coef * (-wcat[k][j+2][i].y -2. * wcat[k][j+1][i].y +3. * wcat[k][j  ][i].y) + wcat[k][j+1][i].y) +		     
		up * (coef * (-wcat[k][j-1][i].y -2. * wcat[k][j  ][i].y +3. * wcat[k][j+1][i].y) + wcat[k][j][i  ].y);
	      fp2[k][j][i].z += 				     
		um * (coef * (-wcat[k][j+2][i].z -2. * wcat[k][j+1][i].z +3. * wcat[k][j  ][i].z) + wcat[k][j+1][i].z) +		     
		up * (coef * (-wcat[k][j-1][i].z -2. * wcat[k][j  ][i].z +3. * wcat[k][j+1][i].z) + wcat[k][j][i  ].z);
	    }
	  }
	}
	else if  ( j==0 || (nvert[k][j-1][i]) > 0.1) {
	  if (les && user->bctype[2]==7 && j==0) {
	    fp2[k][j][i].x = (ucon - vcon) * (ucat[k][j][i].x + ucat[k][j+1][i].x) ;
	    fp2[k][j][i].y = (ucon - vcon) * (ucat[k][j][i].y + ucat[k][j+1][i].y);
	    fp2[k][j][i].z = (ucon - vcon) * (ucat[k][j][i].z + ucat[k][j+1][i].z) ;
	    if (rotateframe) {
	      fp2[k][j][i].x +=ucon * (wcat[k][j  ][i].x +wcat[k][j+1][i].x);
	      fp2[k][j][i].y +=ucon * (wcat[k][j  ][i].y +wcat[k][j+1][i].y);
	      fp2[k][j][i].z +=ucon * (wcat[k][j  ][i].z +wcat[k][j+1][i].z);
	      
	    }
	  }else{
	    if (user->bctype[2]==7 && j==0){
	      if ( nvert[k][j-1][i]<0.1 && nvert[k][j+1][i]<0.1){//Mohsen Feb 2012
		
		fp2[k][j][i].x = 
		  (vm) * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
		  (vp) * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j][i  ].x);
		fp2[k][j][i].y = 				     
		  (vm) * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
		  (vp) * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);
		fp2[k][j][i].z = 				     
		  (vm) * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
		  (vp) * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);
	      }
	    }else{
	      fp2[k][j][i].x = 
		(vm) * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
		(vp) * (coef * (-ucat[k][j  ][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j][i  ].x);		     
	      fp2[k][j][i].y = 				     
		(vm) * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
		(vp) * (coef * (-ucat[k][j  ][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);		     
	      fp2[k][j][i].z = 				     
		(vm) * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
		(vp) * (coef * (-ucat[k][j  ][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);
	    }
	    if (rotateframe){
	      fp2[k][j][i].x += 
		um * (coef * (-wcat[k][j+2][i].x -2. * wcat[k][j+1][i].x +3. * wcat[k][j  ][i].x) + wcat[k][j+1][i].x) +		     
		up * (coef * (-wcat[k][j  ][i].x -2. * wcat[k][j  ][i].x +3. * wcat[k][j+1][i].x) +wcat[k][j][i  ].x);		     
	      fp2[k][j][i].y += 				     
		um * (coef * (-wcat[k][j+2][i].y -2. * wcat[k][j+1][i].y +3. * wcat[k][j  ][i].y) +wcat[k][j+1][i].y) +		     
		up * (coef * (-wcat[k][j  ][i].y -2. * wcat[k][j  ][i].y +3. * wcat[k][j+1][i].y) +wcat[k][j][i  ].y);		     
	      fp2[k][j][i].z += 				     
		um * (coef * (-wcat[k][j+2][i].z - 2. * wcat[k][j+1][i].z + 3. * wcat[k][j  ][i].z) + wcat[k][j+1][i].z) +		     
		up * (coef * (-wcat[k][j  ][i].z - 2. * wcat[k][j  ][i].z + 3. * wcat[k][j+1][i].z) + wcat[k][j][i  ].z);
	      
	    }
	  }
	}
	else if ( j==my-2 ||(nvert[k][j+1][i]) > 0.1) {
	  if (les && user->bctype[2]==7 && j==my-2) {
	    fp2[k][j][i].x = (ucon - vcon) * (ucat[k][j][i].x + ucat[k][j+1][i].x) ;
	    fp2[k][j][i].y = (ucon - vcon) * (ucat[k][j][i].y + ucat[k][j+1][i].y);
	    fp2[k][j][i].z = (ucon - vcon) * (ucat[k][j][i].z + ucat[k][j+1][i].z) ;
	    if (rotateframe) {
	      fp2[k][j][i].x +=ucon * (wcat[k][j  ][i].x +wcat[k][j+1][i].x);
	      fp2[k][j][i].y +=ucon * (wcat[k][j  ][i].y +wcat[k][j+1][i].y);
	      fp2[k][j][i].z +=ucon * (wcat[k][j  ][i].z +wcat[k][j+1][i].z);
	      
	    }
	  }else{
	    if (user->bctype[2]==7 && j==my-2){
	      if (nvert[k][j+1][i]<0.1 && nvert[k][j-1][i]>0.1){
	     
		fp2[k][j][i].x = 
		  (vm) * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
		  (vp) * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j][i  ].x);
		fp2[k][j][i].y = 				     
		  (vm) * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
		  (vp) * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);
		fp2[k][j][i].z = 				     
		  (vm) * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
		  (vp) * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);
	      }
	    }else{
	      fp2[k][j][i].x = 
		(vm) * (coef * (-ucat[k][j+1][i].x - 2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
		(vp) * (coef * (-ucat[k][j-1][i].x - 2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j][i  ].x);		     
	      fp2[k][j][i].y = 				     
		(vm) * (coef * (-ucat[k][j+1][i].y - 2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
		(vp) * (coef * (-ucat[k][j-1][i].y - 2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);		     
	      fp2[k][j][i].z = 				     
		(vm) * (coef * (-ucat[k][j+1][i].z - 2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
		(vp) * (coef * (-ucat[k][j-1][i].z - 2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);
	    }
	    if (rotateframe){
	      fp2[k][j][i].x += 
		um * (coef * (-wcat[k][j+1][i].x -2. * wcat[k][j+1][i].x +3. * wcat[k][j  ][i].x) +wcat[k][j+1][i].x) +		     
		up * (coef * (-wcat[k][j-1][i].x -2. * wcat[k][j  ][i].x +3. * wcat[k][j+1][i].x) +wcat[k][j][i  ].x);		     
	      fp2[k][j][i].y += 				     
		um * (coef * (-wcat[k][j+1][i].y -2. * wcat[k][j+1][i].y +3. * wcat[k][j  ][i].y) +wcat[k][j+1][i].y) +		     
		up * (coef * (-wcat[k][j-1][i].y -2. * wcat[k][j  ][i].y +3. * wcat[k][j+1][i].y) +wcat[k][j][i  ].y);		     
	      fp2[k][j][i].z += 				     
		um * (coef * (-wcat[k][j+1][i].z -2. * wcat[k][j+1][i].z +3. * wcat[k][j  ][i].z) +wcat[k][j+1][i].z) +		     
		up * (coef * (-wcat[k][j-1][i].z -2. * wcat[k][j  ][i].z +3. * wcat[k][j+1][i].z) +wcat[k][j][i  ].z);
	    }
	  }
	}
      }
    }
  }
  /* k direction */
  for (k=lzs-1; k<lze; k++) {
    for(j=lys; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	
	ucon = ( ucont[k][j][i].z )*0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);
	
	vcon = (vcont[k][j][i].z) * 0.5;  
	
	vm = vcon + fabs(vcon);
	vp = vcon - fabs(vcon);
	vm = um - vm;
	vp = up - vp;

	if (TwoD==3) {
	  fp3[k][j][i].x = 0.;
	  fp3[k][j][i].y = 0.;
	  fp3[k][j][i].z = 0.;
	} else if (k>0 && k<mz-2 &&
		   (nvert[k+1][j][i] < 0.1 || nvert[k+1][j][i] > innerblank) &&
		   (nvert[k-1][j][i] < 0.1 || nvert[k-1][j][i] > innerblank)) {
	  if (les) {
	    fp3[k][j][i].x =(ucon - vcon) * (ucat[k][j][i].x +ucat[k+1][j][i].x) ;
	    fp3[k][j][i].y =(ucon - vcon) * (ucat[k][j][i].y +ucat[k+1][j][i].y) ;
	    fp3[k][j][i].z =(ucon - vcon) * (ucat[k][j][i].z +ucat[k+1][j][i].z) ;
	    
	    if (rotateframe){
	      fp3[k][j][i].x +=ucon * (wcat[k][j][i].x +wcat[k+1][j][i].x) ;
	      fp3[k][j][i].y +=ucon * (wcat[k][j][i].y +wcat[k+1][j][i].y) ;
	      fp3[k][j][i].z +=ucon * (wcat[k][j][i].z +wcat[k+1][j][i].z) ;
	    }
	  } else {
	    fp3[k][j][i].x = 
	      (vm) * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	      (vp) * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	    fp3[k][j][i].y = 	       			     
	      (vm) * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	      (vp) * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	    fp3[k][j][i].z = 	       			     
	      (vm) * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	      (vp) * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	    if (rotateframe){
	      fp3[k][j][i].x += 
		um * (coef * (-wcat[k+2][j][i].x -2. * wcat[k+1][j][i].x +3. * wcat[k  ][j][i].x) +wcat[k+1][j][i].x)   +		     
		up * (coef * (-wcat[k-1][j][i].x -2. * wcat[k  ][j][i].x +3. * wcat[k+1][j][i].x) +wcat[k][j][i  ].x);  		     
	      fp3[k][j][i].y += 	       			     
		um * (coef * (-wcat[k+2][j][i].y -2. * wcat[k+1][j][i].y +3. * wcat[k  ][j][i].y) +wcat[k+1][j][i].y)   +		     
		up * (coef * (-wcat[k-1][j][i].y -2. * wcat[k  ][j][i].y +3. * wcat[k+1][j][i].y) +wcat[k][j][i  ].y);  		     
	      fp3[k][j][i].z += 	       			     
		um * (coef * (-wcat[k+2][j][i].z -2. * wcat[k+1][j][i].z +3. * wcat[k  ][j][i].z) +wcat[k+1][j][i].z)   +		     
		up * (coef * (-wcat[k-1][j][i].z -2. * wcat[k  ][j][i].z +3. * wcat[k+1][j][i].z) +wcat[k][j][i  ].z);
	    }
	  }	  
	}else if (k<mz-2 && ( k==0 ||(nvert[k-1][j][i]) > 0.1)) {
	  if (les && user->bctype[4]==7 && k==0) {
	    fp3[k][j][i].x =(ucon - vcon) * (ucat[k  ][j][i].x +ucat[k+1][j][i].x) ;
	    fp3[k][j][i].y =(ucon - vcon) * (ucat[k  ][j][i].y +ucat[k+1][j][i].y) ;
	    fp3[k][j][i].z =(ucon - vcon) * (ucat[k  ][j][i].z +ucat[k+1][j][i].z) ;
	    
	    if (rotateframe){
	      fp3[k][j][i].x +=ucon * (wcat[k  ][j][i].x +wcat[k+1][j][i].x) ;
	      fp3[k][j][i].y +=ucon * (wcat[k  ][j][i].y +wcat[k+1][j][i].y) ;
	      fp3[k][j][i].z +=ucon * (wcat[k  ][j][i].z +wcat[k+1][j][i].z) ;
	    }
	  } else {
	    if (user->bctype[4]==7 && k==0){
	      if ( nvert[k-1][j][i]<0.1 && nvert[k+1][j][i]<0.1){
		
		fp3[k][j][i].x = 
		  (vm) * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
		  (vp) * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
		fp3[k][j][i].y = 	       			     
		  (vm) * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
		  (vp) * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
		fp3[k][j][i].z = 	       			     
		  (vm) * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
		  (vp) * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
		
	      }
	    }else{
	      fp3[k][j][i].x = 
		(vm) * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
		(vp) * (coef * (-ucat[k  ][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	      fp3[k][j][i].y = 	       			     
		(vm) * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
		(vp) * (coef * (-ucat[k  ][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	      fp3[k][j][i].z = 	       			     
		(vm) * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
		(vp) * (coef * (-ucat[k  ][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	    }
	    if (rotateframe){
	      fp3[k][j][i].x += 
		um * (coef * (-wcat[k+2][j][i].x -2. * wcat[k+1][j][i].x +3. * wcat[k  ][j][i].x) +wcat[k+1][j][i].x)   +		     
		up * (coef * (-wcat[k  ][j][i].x -2. * wcat[k  ][j][i].x +3. * wcat[k+1][j][i].x) +wcat[k][j][i  ].x);  		     
	      fp3[k][j][i].y += 	       			     
		um * (coef * (-wcat[k+2][j][i].y -2. * wcat[k+1][j][i].y +3. * wcat[k  ][j][i].y) +wcat[k+1][j][i].y)   +		     
		up * (coef * (-wcat[k  ][j][i].y -2. * wcat[k  ][j][i].y +3. * wcat[k+1][j][i].y) +wcat[k][j][i  ].y);  		     
	      fp3[k][j][i].z += 	       			     
		um * (coef * (-wcat[k+2][j][i].z -2. * wcat[k+1][j][i].z +3. * wcat[k  ][j][i].z) +wcat[k+1][j][i].z)   +		     
		up * (coef * (-wcat[k  ][j][i].z -2. * wcat[k  ][j][i].z +3. * wcat[k+1][j][i].z) +wcat[k][j][i  ].z);
	    }
	  }
	}
	else if (k>0 && ( k==mz-2 ||(nvert[k+1][j][i]) > 0.1)) {
	 if (les && user->bctype[4]==7 && k==mz-2) {
	    fp3[k][j][i].x =(ucon - vcon) * (ucat[k  ][j][i].x +ucat[k+1][j][i].x) ;
	    fp3[k][j][i].y =(ucon - vcon) * (ucat[k  ][j][i].y +ucat[k+1][j][i].y) ;
	    fp3[k][j][i].z =(ucon - vcon) * (ucat[k  ][j][i].z +ucat[k+1][j][i].z) ;
	    
	    if (rotateframe){
	      fp3[k][j][i].x +=ucon * (wcat[k  ][j][i].x +wcat[k+1][j][i].x) ;
	      fp3[k][j][i].y +=ucon * (wcat[k  ][j][i].y +wcat[k+1][j][i].y) ;
	      fp3[k][j][i].z +=ucon * (wcat[k  ][j][i].z +wcat[k+1][j][i].z) ;
	    }
	  } else {
	    if(user->bctype[4]==7 && k==mz-2){
	      if(nvert[k-1][j][i]<0.1 && nvert[k+1][j][i]<0.1){
		
		fp3[k][j][i].x = 
		  (vm) * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
		  (vp) * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
		fp3[k][j][i].y = 	       			     
		  (vm) * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
		  (vp) * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
		fp3[k][j][i].z = 	       			     
		  (vm) * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
		  (vp) * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	      }
	    }else{
	      fp3[k][j][i].x = 
		(vm) * (coef * (-ucat[k+1][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
		(vp) * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	      fp3[k][j][i].y = 	       			     
		(vm) * (coef * (-ucat[k+1][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
		(vp) * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	      fp3[k][j][i].z = 	       			     
		(vm) * (coef * (-ucat[k+1][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	      (vp) * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	    }
	    if (rotateframe){
	      fp3[k][j][i].x += 
		um * (coef * (-wcat[k+1][j][i].x -2. * wcat[k+1][j][i].x +3. * wcat[k  ][j][i].x) +wcat[k+1][j][i].x)   +		     
		up * (coef * (-wcat[k-1][j][i].x -2. * wcat[k  ][j][i].x +3. * wcat[k+1][j][i].x) +wcat[k][j][i  ].x);  		     
	      fp3[k][j][i].y += 	       			     
		um * (coef * (-wcat[k+1][j][i].y -2. * wcat[k+1][j][i].y +3. * wcat[k  ][j][i].y) +wcat[k+1][j][i].y)   +		     
		up * (coef * (-wcat[k-1][j][i].y -2. * wcat[k  ][j][i].y +3. * wcat[k+1][j][i].y) +wcat[k][j][i  ].y);  		     
	      fp3[k][j][i].z += 	       		 
		um * (coef * (-wcat[k+1][j][i].z -2. * wcat[k+1][j][i].z +3. * wcat[k  ][j][i].z) +wcat[k+1][j][i].z)   +		     
		up * (coef * (-wcat[k-1][j][i].z -2. * wcat[k  ][j][i].z +3. * wcat[k+1][j][i].z) +wcat[k][j][i  ].z);
	    }
	 }
	}
      }
    }
  }

  /* Calculate the convective terms under cartesian coordinates */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	conv[k][j][i].x = 
	  fp1[k][j][i].x - fp1[k][j][i-1].x +
	  fp2[k][j][i].x - fp2[k][j-1][i].x +
	  fp3[k][j][i].x - fp3[k-1][j][i].x;
	conv[k][j][i].y = 
	  fp1[k][j][i].y - fp1[k][j][i-1].y +
	  fp2[k][j][i].y - fp2[k][j-1][i].y +
	  fp3[k][j][i].y - fp3[k-1][j][i].y;
	conv[k][j][i].z = 
	  fp1[k][j][i].z - fp1[k][j][i-1].z +
	  fp2[k][j][i].z - fp2[k][j-1][i].z +
	  fp3[k][j][i].z - fp3[k-1][j][i].z;

      }
    }
  }



  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);

  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->lVcont, &vcont);
  if (rotateframe)  DMDAVecRestoreArray(fda, user->lWcat, &wcat);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMDAVecRestoreArray(fda, Conv,  &conv);

  DMDAVecRestoreArray(fda, Fp1, &fp1);
  DMDAVecRestoreArray(fda, Fp2, &fp2);
  DMDAVecRestoreArray(fda, Fp3, &fp3);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lAj, &aj);

  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);

  return (0);
}


PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv)
{
   
  Cmpnts	***ucont, ***ucat;
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  Cmpnts	***fp1, ***fp2, ***fp3;
  Cmpnts	***conv;
 
  PetscReal	ucon, up, um;
  PetscReal	coef = 0.125, innerblank=7.;

  PetscInt	lxs, lxe, lys, lye, lzs, lze, gxs, gxe, gys, gye, gzs,gze;

  PetscReal	***nvert,***aj;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  PetscInt      central=0;
  PetscOptionsGetInt(PETSC_NULL, "-central", &central, PETSC_NULL);

  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(fda, Conv,  &conv);
  DMDAVecGetArray(da, user->lAj, &aj);

  VecDuplicate(Ucont, &Fp1);
  VecDuplicate(Ucont, &Fp2);
  VecDuplicate(Ucont, &Fp3);

  DMDAVecGetArray(fda, Fp1, &fp1);
  DMDAVecGetArray(fda, Fp2, &fp2);
  DMDAVecGetArray(fda, Fp3, &fp3);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  

  /* We have two different sets of node: 1. grid node, the physical points
     where grid lines intercross; 2. storage node, where we store variables.
     All node without explicitly specified as "grid node" refers to
     storage node.

     The integer node is defined at cell center while half node refers to
     the actual grid node. (The reason to choose this arrangement is we need
     ghost node, which is half node away from boundaries, to specify boundary
     conditions. By using this storage arrangement, the actual storage need
     is (IM+1) * (JM + 1) * (KM+1) where IM, JM, & KM refer to the number of
     grid nodes along i, j, k directions.)

     DA, the data structure used to define the storage of 3D arrays, is defined
     as mx * my * mz. mx = IM+1, my = JM+1, mz = KM+1.

     Staggered grid arrangement is used in this solver.
     Pressure is stored at interger node (hence the cell center) and volume
     fluxes defined on the center of each surface of a given control volume
     is stored on the cloest upper integer node. */

  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;
  
  //Mohsen Sep 2012//
/* First update the computational ghost points velocity for periodic boundary conditions
 just for this subroutine because of Quick scheme for velocity deravatives */
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==1) && (k==0 || k==1)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d u is %.15le v is %.15le  w is %.15le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z ); */
/*       } */
/*     } */
/*   } */
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i-2].x;
	  ucat[k][j][i].y=ucat[k][j][i-2].y;
	  ucat[k][j][i].z=ucat[k][j][i-2].z;
	  nvert[k][j][i]=nvert[k][j][i-2];
	}
      }
    }
    if (xe==mx){
      i=mx;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i+2].x;
	  ucat[k][j][i].y=ucat[k][j][i+2].y;
	  ucat[k][j][i].z=ucat[k][j][i+2].z;
	  nvert[k][j][i]=nvert[k][j][i+2];
	}
      }
    }
  }  
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j-2][i].x;
	  ucat[k][j][i].y=ucat[k][j-2][i].y;
	  ucat[k][j][i].z=ucat[k][j-2][i].z;
	  nvert[k][j][i]=nvert[k][j-2][i];
	}
      }
    }
    if (ye==my){
      j=my;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j+2][i].x;
	  ucat[k][j][i].y=ucat[k][j+2][i].y;
	  ucat[k][j][i].z=ucat[k][j+2][i].z;
	  nvert[k][j][i]=nvert[k][j+2][i];
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k-2][j][i].x;
	  ucat[k][j][i].y=ucat[k-2][j][i].y;
	  ucat[k][j][i].z=ucat[k-2][j][i].z;
	  nvert[k][j][i]=nvert[k-2][j][i];
	}
      }
    }
    if (ze==mz){
      k=mz;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k+2][j][i].x;
	  ucat[k][j][i].y=ucat[k+2][j][i].y;
	  ucat[k][j][i].z=ucat[k+2][j][i].z;
	  nvert[k][j][i]=nvert[k+2][j][i];
	}
      }
    }
  }

  VecSet(Conv, 0.0);
 
 /* Calculating the convective terms on cell centers.
    First calcualte the contribution from i direction
    The flux is evaluated by QUICK scheme */

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs-1; i<lxe; i++){
       
       
	ucon = ucont[k][j][i].x * 0.5;
       
	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);
	
	if (i>0  && i<mx-2 &&
		 (nvert[k][j][i+1] < 0.1 || nvert[k][j][i+1]>innerblank) &&
		 (nvert[k][j][i-1] < 0.1 || nvert[k][j][i-1]>innerblank)) { // interial nodes
	  if ((les || central)) {
	    fp1[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j][i+1].x );
	    fp1[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j][i+1].y );
	    fp1[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j][i+1].z );

	  } else {
	  fp1[k][j][i].x =  
	    um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	    up * (coef * (-ucat[k][j][i-1].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	  fp1[k][j][i].y = 
	    um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	    up * (coef * (-ucat[k][j][i-1].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	  fp1[k][j][i].z = 
	    um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	    up * (coef * (-ucat[k][j][i-1].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }
	}
	else  if ((les || central) && (i==0 || i==mx-2) &&
		  (nvert[k][j][i+1] < 0.1 || nvert[k][j][i+1]>innerblank) &&
		  (nvert[k][j][i  ] < 0.1 || nvert[k][j][i  ]>innerblank)) 
	  {
	    fp1[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j][i+1].x );
	    fp1[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j][i+1].y );
	    fp1[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j][i+1].z );
	  }
	else if (i==0 ||(nvert[k][j][i-1] > 0.1) ) {
	  if (user->bctype[0]==7 && i==0 && (nvert[k][j][i-1]<0.1 && nvert[k][j][i+1]<0.1)){//Mohsen Feb 12 
	    fp1[k][j][i].x =  
	      um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	      up * (coef * (-ucat[k][j][i-1].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	    fp1[k][j][i].y = 
	      um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	      up * (coef * (-ucat[k][j][i-1].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	    fp1[k][j][i].z = 
	      um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	      up * (coef * (-ucat[k][j][i-1].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }else{
	    fp1[k][j][i].x = 
	      um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	      up * (coef * (-ucat[k][j][i  ].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	    fp1[k][j][i].y = 
	      um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	      up * (coef * (-ucat[k][j][i  ].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	    fp1[k][j][i].z = 
	      um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	      up * (coef * (-ucat[k][j][i  ].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);	  
	  }
	}
	else if (i==mx-2 ||(nvert[k][j][i+1]) > 0.1) {
	  if (user->bctype[0]==7 && i==mx-2 &&(nvert[k][j][i-1]<0.1 && nvert[k][j][i+1]<0.1)){//Mohsen Feb 12 
	    fp1[k][j][i].x =  
	    um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	    up * (coef * (-ucat[k][j][i-1].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	  fp1[k][j][i].y = 
	    um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	    up * (coef * (-ucat[k][j][i-1].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	  fp1[k][j][i].z = 
	    um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	    up * (coef * (-ucat[k][j][i-1].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }else{
	  fp1[k][j][i].x = 
	    um * (coef * (-ucat[k][j][i+1].x -2. * ucat[k][j][i+1].x +3. * ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	    up * (coef * (-ucat[k][j][i-1].x -2. * ucat[k][j][i  ].x +3. * ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	  fp1[k][j][i].y = 
	    um * (coef * (-ucat[k][j][i+1].y -2. * ucat[k][j][i+1].y +3. * ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	    up * (coef * (-ucat[k][j][i-1].y -2. * ucat[k][j][i  ].y +3. * ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	  fp1[k][j][i].z = 
	    um * (coef * (-ucat[k][j][i+1].z -2. * ucat[k][j][i+1].z +3. * ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	    up * (coef * (-ucat[k][j][i-1].z -2. * ucat[k][j][i  ].z +3. * ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }
	}
      }
    }
  }

  /* j direction */
  for (k=lzs; k<lze; k++) {
    for(j=lys-1; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	ucon = ucont[k][j][i].y * 0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (j>0 && j<my-2 &&
	    (nvert[k][j+1][i] < 0.1 || nvert[k][j+1][i] > innerblank) &&
	    (nvert[k][j-1][i] < 0.1 || nvert[k][j-1][i] > innerblank))  {
	  if ((les || central)) {
	    fp2[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j+1][i].x );
	    fp2[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j+1][i].y );
	    fp2[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j+1][i].z );

	  } else {
	  fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	  }
	}
	else  if ((les || central) && (j==0 || i==my-2) &&
		  (nvert[k][j+1][i] < 0.1 || nvert[k][j+1][i]>innerblank) &&
		  (nvert[k][j  ][i] < 0.1 || nvert[k][j  ][i]>innerblank)) 
	  {
	    fp2[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j+1][i].x );
	    fp2[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j+1][i].y );
	    fp2[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j+1][i].z );
	  }
	else if (j==0 || (nvert[k][j-1][i]) > 0.1) {
	  if (user->bctype[2]==7 && j==0 && (nvert[k][j-1][i]<0.1 && nvert[k][j+1][i]<0.1 )){//Mohsen Feb 12 //
	    fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	  }else{
	  fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j  ][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j  ][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j  ][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	}
	}
	else if (j==my-2 ||(nvert[k][j+1][i]) > 0.1) {
	  if (user->bctype[2]==7 && j==my-2 && (nvert[k][j-1][i]<0.1 && nvert[k][j+1][i]<0.1 )){//Mohsen Feb 12// 
	    fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	  }else{
	  fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+1][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+1][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+1][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);
	  }
	}
      }
    }
  }
 

  /* k direction */
  for (k=lzs-1; k<lze; k++) {
    for(j=lys; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	ucon = ucont[k][j][i].z * 0.5;
       
	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (k>0 && k<mz-2 &&
	    (nvert[k+1][j][i] < 0.1 || nvert[k+1][j][i] > innerblank) &&
	    (nvert[k-1][j][i] < 0.1 || nvert[k-1][j][i] > innerblank)) {
	  if ((les || central)) {
	    fp3[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k+1][j][i].x );
	    fp3[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k+1][j][i].y );
	    fp3[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k+1][j][i].z );

	  } else {
	  fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k  ][j][i].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k  ][j][i].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k  ][j][i].z);
	  }
	}
	else if ((les || central) && (k==0 || k==mz-2) &&
		 (nvert[k+1][j][i] < 0.1 || nvert[k+1][j][i]>innerblank) &&
		 (nvert[k  ][j][i] < 0.1 || nvert[k  ][j][i]>innerblank)) 
	  {
	    fp3[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k+1][j][i].x );
	    fp3[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k+1][j][i].y );
	    fp3[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k+1][j][i].z );
	  }
	else if (k<mz-2 && (k==0 ||(nvert[k-1][j][i]) > 0.1)) {
	  if(user->bctype[4]==7 && k==0 && (nvert[k-1][j][i]<0.1 && nvert[k+1][j][i]<0.1)){//Mohsen Feb 12//
	    fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k  ][j][i].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k  ][j][i].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k  ][j][i].z);
	  }else{
	  fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k  ][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k  ][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k  ][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	  }
	}
	else if (k>0 && (k==mz-2 ||(nvert[k+1][j][i]) > 0.1)) {
	  if (user->bctype[4]==7 && k==mz-2 && (nvert[k-1][j][i]<0.1 && nvert[k+1][j][i]<0.1)){//Mohsen Feb 12//
	    fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k  ][j][i].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k  ][j][i].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k  ][j][i].z);
	  }else{
	  fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+1][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+1][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+1][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	  }
	}
      }
    }
  }
 
  /* Calculate the convective terms under cartesian coordinates */

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	conv[k][j][i].x = 
	  fp1[k][j][i].x - fp1[k][j][i-1].x +
	  fp2[k][j][i].x - fp2[k][j-1][i].x +
	  fp3[k][j][i].x - fp3[k-1][j][i].x;

	conv[k][j][i].y = 
	  fp1[k][j][i].y - fp1[k][j][i-1].y +
	  fp2[k][j][i].y - fp2[k][j-1][i].y +
	  fp3[k][j][i].y - fp3[k-1][j][i].y;

	conv[k][j][i].z =
	  fp1[k][j][i].z - fp1[k][j][i-1].z +
	  fp2[k][j][i].z - fp2[k][j-1][i].z +
	  fp3[k][j][i].z - fp3[k-1][j][i].z;
      }
    }
  } 
  /* for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==1) && (k==1 || k==21 || k==22|| k==200)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d conv.y is %.15le  conv.z is %.15le \n",i,j,k,conv[k][j][i].y,conv[k][j][i].z); */
/*       } */
/*     } */
/*   }  */

  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMDAVecRestoreArray(fda, Conv,  &conv);
  DMDAVecRestoreArray(da, user->lAj, &aj);

  DMDAVecRestoreArray(fda, Fp1, &fp1);
  DMDAVecRestoreArray(fda, Fp2, &fp2);
  DMDAVecRestoreArray(fda, Fp3, &fp3);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);
   
  return (0);
}

PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc)
{
  
  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;

  Cmpnts	***ucont, ***ucat;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***nvert;

  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  Cmpnts	***fp1, ***fp2, ***fp3;
  Cmpnts	***visc;
  PetscReal	***aj, ***iaj, ***jaj, ***kaj;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     ajc;

  PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g21, g31;
  PetscReal	r11, r21, r31, r12, r22, r32, r13, r23, r33;

  PetscScalar	solid,innerblank;

  solid = 0.5;
  innerblank = 7.;
  
  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(fda, Visc,  &visc);

  DMDAVecGetArray(fda, Csi, &csi);
  DMDAVecGetArray(fda, Eta, &eta);
  DMDAVecGetArray(fda, Zet, &zet);

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

  VecDuplicate(Ucont, &Fp1);
  VecDuplicate(Ucont, &Fp2);
  VecDuplicate(Ucont, &Fp3);

  DMDAVecGetArray(fda, Fp1, &fp1);
  DMDAVecGetArray(fda, Fp2, &fp2);
  DMDAVecGetArray(fda, Fp3, &fp3);

  DMDAVecGetArray(da, user->lAj, &aj);

  DMDAGetLocalInfo(da, &info);
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;


  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;
 
  VecSet(Visc,0.0);
  
  PetscReal ***lnu_t;
 
  if(les) {
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);
  } else if (rans) {
   
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);
  }

  /* The visc flux on each surface center is stored at previous integer node */

  DMDAVecGetArray(da, user->lIAj, &iaj);
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==0 ||j==1 || j==2) && (k==21 || k==22|| k==20)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d u is %.15le v is %.15le  w is %.15le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z ); */
/*       } */
/*     } */
/*   } */
  // i direction
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs-1; i<lxe; i++) {
	
	dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
	dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
	dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;
	
	if ((nvert[k][j+1][i  ]> solid && nvert[k][j+1][i  ]<innerblank)  ||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
	  dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x -
		  ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
	  dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y -
		  ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
	  dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z -
		  ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
	}
	else if  ((nvert[k][j-1][i  ]> solid && nvert[k][j-1][i  ]<innerblank)  ||
		  (nvert[k][j-1][i+1]> solid && nvert[k][j-1][i+1]<innerblank)) {
	  dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x -
		  ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
	  dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y -
		  ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
	  dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z -
		  ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
	}
	else {
	  dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x -
		  ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
	  dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y -
		  ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
	  dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z -
		  ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
	}
	
	if ((nvert[k+1][j][i  ]> solid && nvert[k+1][j][i  ]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
	  dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x -
		  ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
	  dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y -
		  ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
	  dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
	}
	else if ((nvert[k-1][j][i  ]> solid && nvert[k-1][j][i  ]<innerblank) ||
		 (nvert[k-1][j][i+1]> solid && nvert[k-1][j][i+1]<innerblank)) {

	  dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x -
		  ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
	  dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y -
		  ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
	  dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z -
		  ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
	}
	else {
	  dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x -
		  ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
	  dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y -
		  ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
	  dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
	}

 	csi0 = icsi[k][j][i].x;
	csi1 = icsi[k][j][i].y;
	csi2 = icsi[k][j][i].z;

	eta0 = ieta[k][j][i].x;
	eta1 = ieta[k][j][i].y;
	eta2 = ieta[k][j][i].z;
	
	zet0 = izet[k][j][i].x;
	zet1 = izet[k][j][i].y;
	zet2 = izet[k][j][i].z;

	g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

	ajc = iaj[k][j][i]; 

	double nu = 1./user->ren, nu_t=0;
	
	if( les || (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j][i+1]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);
	  if ( (user->bctype[0]==1 && i==0) || (user->bctype[1]==1 && i==mx-2) ) nu_t=0;    
	  fp1[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nu_t); 
	  fp1[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nu_t);
	  fp1[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nu_t);
	}
	else {
	  fp1[k][j][i].x = 0;
	  fp1[k][j][i].y = 0;
	  fp1[k][j][i].z = 0;
	}

	fp1[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz+ r11 * csi0 + r21 * csi1 + r31 * csi2 ) * ajc * (nu);
	fp1[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz+ r12 * csi0 + r22 * csi1 + r32 * csi2 ) * ajc * (nu);
	fp1[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz+ r13 * csi0 + r23 * csi1 + r33 * csi2 ) * ajc * (nu);


	if(clark) {
	  double dc, de, dz;
	  Calculate_dxdydz (ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dc, &de, &dz);
	  double dc2=dc*dc, de2=de*de, dz2=dz*dz;
	  
	  double t11 = ( dudc * dudc * dc2 + dude * dude * de2 + dudz * dudz * dz2 );
	  double t12 = ( dudc * dvdc * dc2 + dude * dvde * de2 + dudz * dvdz * dz2 );
	  double t13 = ( dudc * dwdc * dc2 + dude * dwde * de2 + dudz * dwdz * dz2 );
	  double t21 = t12;
	  double t22 = ( dvdc * dvdc * dc2 + dvde * dvde * de2 + dvdz * dvdz * dz2 );
	  double t23 = ( dvdc * dwdc * dc2 + dvde * dwde * de2 + dvdz * dwdz * dz2 );
	  double t31 = t13;
	  double t32 = t23;
	  double t33 = ( dwdc * dwdc * dc2 + dwde * dwde * de2 + dwdz * dwdz * dz2 );
	  
	  fp1[k][j][i].x -= ( t11 * csi0 + t12 * csi1 + t13 * csi2 ) / 12.;
	  fp1[k][j][i].y -= ( t21 * csi0 + t22 * csi1 + t23 * csi2 ) / 12.;
	  fp1[k][j][i].z -= ( t31 * csi0 + t32 * csi1 + t33 * csi2 ) / 12.;
	}

      }
    }
  }
  DMDAVecRestoreArray(da, user->lIAj, &iaj);  
 

  // j direction
  DMDAVecGetArray(da, user->lJAj, &jaj);
  for (k=lzs; k<lze; k++) {
    for (j=lys-1; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
 
	if ((nvert[k][j  ][i+1]> solid && nvert[k][j  ][i+1]<innerblank)||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
	  dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x -
		  ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
	  dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y -
		  ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
	  dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z -
		  ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
	}
	else if ((nvert[k][j  ][i-1]> solid && nvert[k][j  ][i-1]<innerblank) ||
		 (nvert[k][j+1][i-1]> solid && nvert[k][j+1][i-1]<innerblank)) {
	  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x -
		  ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
	  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y -
		  ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
	  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z -
		  ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
	}
	else {
	  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x -
		  ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
	  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y -
		  ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
	  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z -
		  ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
	}

	dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
	dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
	dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;

	if ((nvert[k+1][j  ][i]> solid && nvert[k+1][j  ][i]<innerblank)||
	    (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
	  dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x -
		  ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
	  dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y -
		  ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
	  dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z -
		  ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
	}
	else if ((nvert[k-1][j  ][i]> solid && nvert[k-1][j  ][i]<innerblank)||
		 (nvert[k-1][j+1][i]> solid && nvert[k-1][j+1][i]<innerblank)) {
	  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x -
		  ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
	  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y -
		  ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
	  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z -
		  ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
	}
	else {
	  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x -
		  ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
	  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y -
		  ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
	  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z -
		  ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
	}

 	csi0 = jcsi[k][j][i].x;
	csi1 = jcsi[k][j][i].y;
	csi2 = jcsi[k][j][i].z;

	eta0 = jeta[k][j][i].x;
	eta1 = jeta[k][j][i].y;
	eta2 = jeta[k][j][i].z;
	
	zet0 = jzet[k][j][i].x;
	zet1 = jzet[k][j][i].y;
	zet2 = jzet[k][j][i].z;


	g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
	g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d dvdc is %.15le dvde is %.15le dvdz is %.15le \n",i,j,k,dvdc,dvde,dvdz);
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d dwdc is %.15le dwde is %.15le dwdz is %.15le \n",i,j,k,dwdc,dwde,dwdz);
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d jcsi is %.15le jeta is %.15le jzet is %.15le \n",i,j,k,jcsi[k][j][i].z,jeta[k][j][i].z,jzet[k][j][i].z);
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d r13 is %.15le r23 is %.15le r33 is %.15le \n",i,j,k,r13,r23,r33);



	ajc = jaj[k][j][i];

	double nu = 1./user->ren, nu_t = 0;
		
	if( les || (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j+1][i]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
	  if ( (user->bctype[2]==1 && j==0) || (user->bctype[3]==1 && j==my-2) ) nu_t=0;
		
	  fp2[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nu_t);
	  fp2[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nu_t);
	  fp2[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nu_t);
	}
	else {
	  fp2[k][j][i].x = 0;
	  fp2[k][j][i].y = 0;
	  fp2[k][j][i].z = 0;
	}
		
	fp2[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz+ r11 * eta0 + r21 * eta1 + r31 * eta2 ) * ajc * (nu);
	fp2[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz+ r12 * eta0 + r22 * eta1 + r32 * eta2 ) * ajc * (nu);
	fp2[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz+ r13 * eta0 + r23 * eta1 + r33 * eta2 ) * ajc * (nu);
		
	if(clark) {
	  double dc, de, dz;
	  Calculate_dxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dc, &de, &dz);
	  double dc2=dc*dc, de2=de*de, dz2=dz*dz;
			
	  double t11 = ( dudc * dudc * dc2 + dude * dude * de2 + dudz * dudz * dz2 );
	  double t12 = ( dudc * dvdc * dc2 + dude * dvde * de2 + dudz * dvdz * dz2 );
	  double t13 = ( dudc * dwdc * dc2 + dude * dwde * de2 + dudz * dwdz * dz2 );
	  double t21 = t12;
	  double t22 = ( dvdc * dvdc * dc2 + dvde * dvde * de2 + dvdz * dvdz * dz2 );
	  double t23 = ( dvdc * dwdc * dc2 + dvde * dwde * de2 + dvdz * dwdz * dz2 );
	  double t31 = t13;
	  double t32 = t23;
	  double t33 = ( dwdc * dwdc * dc2 + dwde * dwde * de2 + dwdz * dwdz * dz2 );
			
	  fp2[k][j][i].x -= ( t11 * eta0 + t12 * eta1 + t13 * eta2 ) / 12.;
	  fp2[k][j][i].y -= ( t21 * eta0 + t22 * eta1 + t23 * eta2 ) / 12.;
	  fp2[k][j][i].z -= ( t31 * eta0 + t32 * eta1 + t33 * eta2 ) / 12.;
	}
      }
    }
  }

  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  // k direction

  DMDAVecGetArray(da, user->lKAj, &kaj);
  for (k=lzs-1; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((nvert[k  ][j][i+1]> solid && nvert[k  ][j][i+1]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
	  dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x -
		  ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
	  dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y -
		  ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
	  dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z -
		  ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
	}
	else if ((nvert[k  ][j][i-1]> solid && nvert[k  ][j][i-1]<innerblank) ||
		 (nvert[k+1][j][i-1]> solid && nvert[k+1][j][i-1]<innerblank)) {
	  dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x -
		  ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
	  dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y -
		  ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
	  dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z -
		  ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
	}
	else {
	  dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x -
		  ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
	  dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y -
		  ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
	  dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z -
		  ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
	}
	  
	if ((nvert[k  ][j+1][i]> solid && nvert[k  ][j+1][i]<innerblank)||
	    (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
	  dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x -
		  ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
	  dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y -
		  ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
	  dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z -
		  ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
	}
	else if ((nvert[k  ][j-1][i]> solid && nvert[k  ][j-1][i]<innerblank) ||
		 (nvert[k+1][j-1][i]> solid && nvert[k+1][j-1][i]<innerblank)){
	  dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x -
		  ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
	  dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y -
		  ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
	  dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z -
		  ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
	}
	else {
	  dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x -
		  ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
	  dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y -
		  ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
	  dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z -
		  ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
	}

	dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
	dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
	dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;


 	csi0 = kcsi[k][j][i].x;
	csi1 = kcsi[k][j][i].y;
	csi2 = kcsi[k][j][i].z;

	eta0 = keta[k][j][i].x;
	eta1 = keta[k][j][i].y;
	eta2 = keta[k][j][i].z;
	
	zet0 = kzet[k][j][i].x;
	zet1 = kzet[k][j][i].y;
	zet2 = kzet[k][j][i].z;


	g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
	g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
	g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

	ajc = kaj[k][j][i];

	double nu = 1./user->ren, nu_t =0;
		
	if( les || (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k+1][j][i]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
	  if ( (user->bctype[4]==1 && k==0) || (user->bctype[5]==1 && k==mz-2) ) nu_t=0;
		
	  fp3[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu_t);
	  fp3[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu_t);
	  fp3[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu_t);
	}
	else {
	  fp3[k][j][i].x = 0;
	  fp3[k][j][i].y = 0;
	  fp3[k][j][i].z = 0;
	}
	fp3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu);//
	fp3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu);//
	fp3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu);//

	if(clark) {
	  double dc, de, dz;
	  Calculate_dxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dc, &de, &dz);
	  double dc2=dc*dc, de2=de*de, dz2=dz*dz;
			
	  double t11 = ( dudc * dudc * dc2 + dude * dude * de2 + dudz * dudz * dz2 );
	  double t12 = ( dudc * dvdc * dc2 + dude * dvde * de2 + dudz * dvdz * dz2 );
	  double t13 = ( dudc * dwdc * dc2 + dude * dwde * de2 + dudz * dwdz * dz2 );
	  double t21 = t12;
	  double t22 = ( dvdc * dvdc * dc2 + dvde * dvde * de2 + dvdz * dvdz * dz2 );
	  double t23 = ( dvdc * dwdc * dc2 + dvde * dwde * de2 + dvdz * dwdz * dz2 );
	  double t31 = t13;
	  double t32 = t23;
	  double t33 = ( dwdc * dwdc * dc2 + dwde * dwde * de2 + dwdz * dwdz * dz2 );
			
	  fp3[k][j][i].x -= ( t11 * zet0 + t12 * zet1 + t13 * zet2 ) / 12.;
	  fp3[k][j][i].y -= ( t21 * zet0 + t22 * zet1 + t23 * zet2 ) / 12.;
	  fp3[k][j][i].z -= ( t31 * zet0 + t32 * zet1 + t33 * zet2 ) / 12.;
	}
      }
    }
  }

  DMDAVecRestoreArray(da, user->lKAj, &kaj);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	visc[k][j][i].x =
	  (fp1[k][j][i].x - fp1[k][j][i-1].x +
	   fp2[k][j][i].x - fp2[k][j-1][i].x +
	   fp3[k][j][i].x - fp3[k-1][j][i].x);
	
	visc[k][j][i].y =
	  (fp1[k][j][i].y - fp1[k][j][i-1].y +
	   fp2[k][j][i].y - fp2[k][j-1][i].y +
	   fp3[k][j][i].y - fp3[k-1][j][i].y);

	visc[k][j][i].z =
	  (fp1[k][j][i].z - fp1[k][j][i-1].z +
	   fp2[k][j][i].z - fp2[k][j-1][i].z +
	   fp3[k][j][i].z - fp3[k-1][j][i].z);
			
      }
    }
  }
/*   for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp1.z is %.15le \n",i,j,k,fp1[k][j][i].z); */
/* 	if (i==0 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp1.z is %.15le \n",i,j,k,fp1[k][j][i].z); */
/* 	if (i==1 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp2.z is %.15le \n",i,j,k,fp2[k][j][i].z); */
/* 	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp2.z is %.15le \n",i,j,k,fp2[k][j][i].z); */
/* 	if (i==1 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp3.z is %.15le \n",i,j,k,fp3[k][j][i].z); */
/* 	if (i==1 && j==1 && k==20) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp3.z is %.15le \n",i,j,k,fp3[k][j][i].z); */
	 
/*       } */
/*     } */
/*   }  */
  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMDAVecRestoreArray(fda, Visc,  &visc);

  DMDAVecRestoreArray(fda, Csi, &csi);
  DMDAVecRestoreArray(fda, Eta, &eta);
  DMDAVecRestoreArray(fda, Zet, &zet);

  DMDAVecRestoreArray(fda, Fp1, &fp1);
  DMDAVecRestoreArray(fda, Fp2, &fp2);
  DMDAVecRestoreArray(fda, Fp3, &fp3);

  DMDAVecRestoreArray(da, user->lAj, &aj);

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

  if(les) {
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  } else if (rans) {
  
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  }
 

  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);
  
  return(0);
}

PetscErrorCode FormFunction1(UserCtx *user, Vec Rhs)
{

  DM 		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt	mx, my, mz;
  PetscInt	i, j, k;
  PetscReal	dpdc, dpde, dpdz;

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  Vec		Aj  = user->lAj;
  Vec		ICsi = user->lICsi, IEta = user->lIEta, IZet = user->lIZet;
  Vec		JCsi = user->lJCsi, JEta = user->lJEta, JZet = user->lJZet;
  Vec		KCsi = user->lKCsi, KEta = user->lKEta, KZet = user->lKZet;
  Vec		IAj = user->lIAj, JAj = user->lJAj, KAj = user->lKAj;
  Vec		P = user->lP;
  Vec		Conv, Visc;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  Cmpnts	***ucont;
  PetscReal	***p, ***iaj, ***jaj, ***kaj, ***aj;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  Vec		Rc;
  Cmpnts	***rc;

  Vec		Rct;
  Cmpnts	***rct;
 
 
  Cmpnts	***rhs;

  PetscReal	***nvert, ***nvert_o;


  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  mx = info.mx; my = info.my; mz = info.mz;

  DMDAVecGetArray(fda, Csi, &csi);
  DMDAVecGetArray(fda, Eta, &eta);
  DMDAVecGetArray(fda, Zet, &zet);
  DMDAVecGetArray(da,  Aj,  &aj);
           
  DMDAVecGetArray(fda, ICsi, &icsi);
  DMDAVecGetArray(fda, IEta, &ieta);
  DMDAVecGetArray(fda, IZet, &izet);
           
  DMDAVecGetArray(fda, JCsi, &jcsi);
  DMDAVecGetArray(fda, JEta, &jeta);
  DMDAVecGetArray(fda, JZet, &jzet);
           
  DMDAVecGetArray(fda, KCsi, &kcsi);
  DMDAVecGetArray(fda, KEta, &keta);
  DMDAVecGetArray(fda, KZet, &kzet);

  DMDAVecGetArray(fda, Rhs, &rhs);

  DMDAVecGetArray(da, P, &p);
  DMDAVecGetArray(da, IAj, &iaj);
  DMDAVecGetArray(da, JAj, &jaj);
  DMDAVecGetArray(da, KAj, &kaj);
  // Get a working array rc to store the contravariant rhs on cell center
  VecDuplicate(user->lUcont, &Rc);

  VecDuplicate(Rc, &Rct);
  DMDAVecGetArray(fda, Rct, &rct);
  VecDuplicate(Rct, &Conv);
  VecDuplicate(Rct, &Visc);

  // Obtain the Cartesian velocity components from Contravariant velocity
  Contra2Cart(user);
  // Convective term
  if (moveframe || rotateframe) {
   
    Convection_MV(user, user->lUcont, user->lUcat, Conv);
  } else
    Convection(user, user->lUcont, user->lUcat, Conv);

  // Viscous term
  if (invicid)
    VecSet(Visc,0.);
  else
   
    Viscous(user, user->lUcont, user->lUcat, Visc);
  

  // Right hand side terms from convective and viscous terms
 
  VecWAXPY(Rc,-1.0, Conv, Visc);

  DMDAVecGetArray(fda, Rc, &rc);


  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

 
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;
  

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  /* Calculate the contravariant rhs from the cartesian rhs */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	rct[k][j][i].x = aj[k][j][i] * 
	  (0.5 * (csi[k][j][i].x + csi[k][j][i-1].x) * rc[k][j][i].x +
	   0.5 * (csi[k][j][i].y + csi[k][j][i-1].y) * rc[k][j][i].y +
	   0.5 * (csi[k][j][i].z + csi[k][j][i-1].z) * rc[k][j][i].z);
	rct[k][j][i].y = aj[k][j][i] *
	  (0.5 * (eta[k][j][i].x + eta[k][j-1][i].x) * rc[k][j][i].x +
	   0.5 * (eta[k][j][i].y + eta[k][j-1][i].y) * rc[k][j][i].y +
	   0.5 * (eta[k][j][i].z + eta[k][j-1][i].z) * rc[k][j][i].z);
	rct[k][j][i].z = aj[k][j][i] *
	  (0.5 * (zet[k][j][i].x + zet[k-1][j][i].x) * rc[k][j][i].x +
	   0.5 * (zet[k][j][i].y + zet[k-1][j][i].y) * rc[k][j][i].y +
	   0.5 * (zet[k][j][i].z + zet[k-1][j][i].z) * rc[k][j][i].z);
      }
    }
  }

 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==1) && (k==0 || k==1 || k==21 || k==22|| k==200 || k==201)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  rc.y is %.15le rc.z is %.15le \n",i,j,k,rc[k][j][i].y,rc[k][j][i].z ); */
/*       } */
/*     } */
/*   }  */

  DMDAVecRestoreArray(fda, Rct, &rct);
 
  DMDAVecRestoreArray(fda, Rc, &rc);

  PetscBarrier(PETSC_NULL);
 
  DMLocalToLocalBegin(fda, Rct, INSERT_VALUES, Rct);
  DMLocalToLocalEnd(fda, Rct, INSERT_VALUES, Rct);

  DMDAVecGetArray(fda, Rct, &rct);

  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  rct[k][j][i].x=rct[k][j][i-2].x;
	}
      }
    }

    if (xe==mx){
      i=mx-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  rct[k][j][i].x=rct[k][j][i+2].x;
	}
      }
    }
  }

  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  rct[k][j][i].y=rct[k][j-2][i].y;
	}
      }
    }

    if (ye==my){
      j=my-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  rct[k][j][i].y=rct[k][j+2][i].y;
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
     k=zs;
      for (j=lys; j<lye; j++) {
       for (i=lxs; i<lxe; i++) {
	 rct[k][j][i].z=rct[k-2][j][i].z;
       }
     }
    }
   if (ze==mz){
     k=mz-1;
     for (j=lys; j<lye; j++) {
       for (i=lxs; i<lxe; i++) {
	 rct[k][j][i].z=rct[k+2][j][i].z;
       }
     }
   }
  }
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && j==1 && k==1) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  rct.x is %.15le rct.y is %.15le rct.z is %.15le \n",i,j,k,rct[k][j][i].x,rct[k][j][i].y,rct[k][j][i].z ); */
/*       } */
/*     } */
/*   } */
  DMDAVecRestoreArray(fda, Rct, &rct);
  PetscBarrier(PETSC_NULL);
  DMDAVecGetArray(fda, Rct, &rct);
 
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lNvert_o, &nvert_o);
  
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	dpdc = p[k][j][i+1] - p[k][j][i];

	if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		  p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	}

	if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else  if ((k == mz-2 || k==1) &&  user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		  p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	}
	    
	rhs[k][j][i].x =0.5 * (rct[k][j][i].x + rct[k][j][i+1].x);
	
	
	rhs[k][j][i].x -=
	  (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		   icsi[k][j][i].y * icsi[k][j][i].y +
		   icsi[k][j][i].z * icsi[k][j][i].z)+
	   dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		   ieta[k][j][i].y * icsi[k][j][i].y +
		   ieta[k][j][i].z * icsi[k][j][i].z)+
	   dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		   izet[k][j][i].y * icsi[k][j][i].y +
		   izet[k][j][i].z * icsi[k][j][i].z)) * iaj[k][j][i];

	if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		  p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	}
	
	dpde = p[k][j+1][i] - p[k][j][i];
	
	if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		  p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	}

	rhs[k][j][i].y =0.5 * (rct[k][j][i].y + rct[k][j+1][i].y);
	  
	  
	rhs[k][j][i].y -=
	  (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		   jcsi[k][j][i].y * jeta[k][j][i].y +
		   jcsi[k][j][i].z * jeta[k][j][i].z) +
	   dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		   jeta[k][j][i].y * jeta[k][j][i].y +
		   jeta[k][j][i].z * jeta[k][j][i].z) +
	   dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		   jzet[k][j][i].y * jeta[k][j][i].y +
		   jzet[k][j][i].z * jeta[k][j][i].z)) * jaj[k][j][i];

	if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		  p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	}

	if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		  p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	}

	dpdz = (p[k+1][j][i] - p[k][j][i]);
	
	rhs[k][j][i].z =0.5 * (rct[k][j][i].z + rct[k+1][j][i].z);
	
	rhs[k][j][i].z -=
	  (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		   kcsi[k][j][i].y * kzet[k][j][i].y +
		   kcsi[k][j][i].z * kzet[k][j][i].z) +
	   dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		   keta[k][j][i].y * kzet[k][j][i].y +
		   keta[k][j][i].z * kzet[k][j][i].z) +
	   dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		   kzet[k][j][i].y * kzet[k][j][i].y +
		   kzet[k][j][i].z * kzet[k][j][i].z)) * kaj[k][j][i];
	
      }
    }
  }
  
  
  //Mohsen March 2012//
  
  // rhs.x at boundaries for periodic bc at i direction//
  if (user->bctype[0]==7 && xs==0){
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {	  
	i=xs;

	dpdc = p[k][j][i+1] - p[k][j][i];
	
	if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		  p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	}

	if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else  if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		  p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	}

	rhs[k][j][i].x =0.5 * (rct[k][j][i].x + rct[k][j][i+1].x);
	rhs[k][j][i].x -=
	  (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		   icsi[k][j][i].y * icsi[k][j][i].y +
		   icsi[k][j][i].z * icsi[k][j][i].z)+
	   dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		   ieta[k][j][i].y * icsi[k][j][i].y +
		   ieta[k][j][i].z * icsi[k][j][i].z)+
	   dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		   izet[k][j][i].y * icsi[k][j][i].y +
		   izet[k][j][i].z * icsi[k][j][i].z)) * iaj[k][j][i];
      }
    }
  }
  
// rhs.y at boundaries for periodic bc at j direction//
  if (user->bctype[2]==7 && ys==0){
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {	  

	j=ys;

	if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		  p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	}

	dpde = p[k][j+1][i] - p[k][j][i];

	if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		  p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	}

	rhs[k][j][i].y =0.5 * (rct[k][j][i].y + rct[k][j+1][i].y);
	 
       	rhs[k][j][i].y -=
	  (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		   jcsi[k][j][i].y * jeta[k][j][i].y +
		   jcsi[k][j][i].z * jeta[k][j][i].z)+
	   dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		   jeta[k][j][i].y * jeta[k][j][i].y +
		   jeta[k][j][i].z * jeta[k][j][i].z)+
	   dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		   jzet[k][j][i].y * jeta[k][j][i].y +
		   jzet[k][j][i].z * jeta[k][j][i].z)) * jaj[k][j][i];

      }
    }
  }
  
  // rhs.z at boundaries for periodic bc at k direction//
  if (user->bctype[4]==7&& zs==0){
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	k=zs; 

	if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		  p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	}

	if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		  p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	}

	dpdz = (p[k+1][j][i] - p[k][j][i]);
	
	rhs[k][j][i].z =0.5 * (rct[k][j][i].z + rct[k+1][j][i].z);
	
	rhs[k][j][i].z -=
	  (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		   kcsi[k][j][i].y * kzet[k][j][i].y +
		   kcsi[k][j][i].z * kzet[k][j][i].z)+
	   dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		   keta[k][j][i].y * kzet[k][j][i].y +
		   keta[k][j][i].z * kzet[k][j][i].z)+
	   dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		   kzet[k][j][i].y * kzet[k][j][i].y +
		   kzet[k][j][i].z * kzet[k][j][i].z)) * kaj[k][j][i];
	
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (TwoD==1)
	  rhs[k][j][i].x =0.;
	else if (TwoD==2)
	  rhs[k][j][i].y =0.;
	else if (TwoD==3)
	  rhs[k][j][i].z =0.;
	
	if (nvert[k][j][i]>0.1) {
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	  rhs[k][j][i].z = 0;
	}
	if (nvert[k][j][i+1]>0.1) {
	  rhs[k][j][i].x=0;
	}
	if (nvert[k][j+1][i]>0.1) {
	  rhs[k][j][i].y=0;
	}
	if (nvert[k+1][j][i]>0.1) {
	  rhs[k][j][i].z=0;
	}
      }
    }
  }

 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && j==1 && k==1) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d rhs.x is %.15le rhs.y is %.15le rhs.z is %.15le \n",i,j,k,rhs[k][j][i].x,rhs[k][j][i].y,rhs[k][j][i].z ); */
/*       } */
/*     } */
/*   } */
 
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);   
  DMDAVecRestoreArray(fda, Rhs, &rhs); 
  DMDAVecRestoreArray(fda, Rct, &rct);

  DMDAVecRestoreArray(fda, Csi, &csi);
  DMDAVecRestoreArray(fda, Eta, &eta);
  DMDAVecRestoreArray(fda, Zet, &zet);
  DMDAVecRestoreArray(da,  Aj,  &aj);

  DMDAVecRestoreArray(fda, ICsi, &icsi);
  DMDAVecRestoreArray(fda, IEta, &ieta);
  DMDAVecRestoreArray(fda, IZet, &izet);

  DMDAVecRestoreArray(fda, JCsi, &jcsi);
  DMDAVecRestoreArray(fda, JEta, &jeta);
  DMDAVecRestoreArray(fda, JZet, &jzet);
         
  DMDAVecRestoreArray(fda, KCsi, &kcsi);
  DMDAVecRestoreArray(fda, KEta, &keta);
  DMDAVecRestoreArray(fda, KZet, &kzet);

  DMDAVecRestoreArray(da, P, &p);
  DMDAVecRestoreArray(da, IAj, &iaj);
  DMDAVecRestoreArray(da, JAj, &jaj);
  DMDAVecRestoreArray(da, KAj, &kaj);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);
  
  VecDestroy(&Conv);
  VecDestroy(&Visc);
  VecDestroy(&Rc);
  VecDestroy(&Rct);
  
  return(0);
}

PetscErrorCode Spectral(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscReal cfl = 0.5;
  PetscReal vnn = 0.5;

  PetscOptionsGetReal(PETSC_NULL, "-cfl", &cfl, PETSC_NULL);
  PetscInt i, j, k;
  PetscReal abu, abv, abw, ren = user->ren;
  Cmpnts ***ucont, ***csi, ***eta, ***zet;
  PetscReal ***dt, ***aj;
  PetscReal g1, g2, g3, eigen1, eigen2, eigen3, temp, dtii, dtij, dtik;
  PetscReal dtvi, dtvj, dtvk;
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->Dt, &dt);
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	g1 = (csi[k][j][i].x * csi[k][j][i].x +
	      csi[k][j][i].y * csi[k][j][i].y +
	      csi[k][j][i].z * csi[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	g2 = (eta[k][j][i].x * eta[k][j][i].x +
	      eta[k][j][i].y * eta[k][j][i].y +
	      eta[k][j][i].z * eta[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	g3 = (zet[k][j][i].x * zet[k][j][i].x +
	      zet[k][j][i].y * zet[k][j][i].y +
	      zet[k][j][i].z * zet[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	abu = fabs(ucont[k][j][i].x + ucont[k][j][i-1].x) * 0.5 * aj[k][j][i];
	abv = fabs(ucont[k][j][i].y + ucont[k][j][i-1].y) * 0.5 * aj[k][j][i];
	abw = fabs(ucont[k][j][i].z + ucont[k][j][i-1].z) * 0.5 * aj[k][j][i];

	if (g1<1.e-10) g1 = 1.;
	if (g2<1.e-10) g2 = 1.;
	if (g3<1.e-10) g3 = 1.;
	eigen1 = abu + sqrt(abu*abu + g1);
	eigen2 = abv + sqrt(abv*abv + g2);
	eigen3 = abw + sqrt(abw*abw + g3);

	dtii = cfl / eigen1;
	dtij = cfl / eigen2;
	dtik = cfl / eigen3;

	temp = vnn * ren;
	dtvi = temp / g1;
	dtvj = temp / g2;
	dtvk = temp / g3;

	temp = PetscMin(dtii, dtij);
	temp = PetscMin(temp, dtik);
	temp = PetscMin(temp, dtvi);
	temp = PetscMin(temp, dtvj);
	dt[k][j][i] = PetscMin(temp, dtvk);
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->Dt, &dt);
  return 0;
}
