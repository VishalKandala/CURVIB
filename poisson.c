#include "variables.h"
#include "petscksp.h"
#include "petscpc.h"
//#include "petscpcmg.h"

extern PetscInt block_number;
extern PetscInt immersed,MHV,LV,aneurysm,visflg,LVAD;
extern PetscInt ti;

PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
//PetscInt lidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
PetscErrorCode Contra2Cart(UserCtx *user);
PetscErrorCode GhostNodeVelocity(UserCtx *user);
PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, 
			  PetscReal *ibm_Area, PetscInt flg);
PetscErrorCode VolumeFlux_rev(UserCtx *user, PetscReal *ibm_Flux, 
			      PetscReal *ibm_Area, PetscInt flg);

PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnom, void *dummy)
{
  UserCtx *user = (UserCtx*)dummy;
  Vec x;
  PetscInt tttt, mi, j, k;
  PetscReal norm;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  VecDuplicate(user->P, &x);
  KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &x);
  VecMax(x, &tttt, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "KSP Max %d %le\n", tttt, norm);
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==tttt) {
	  PetscPrintf(PETSC_COMM_SELF, "KspMMax %d %d %d %d %le\n", Gidx(mi,j,k,user),mi,j, k, norm);
	}
      }
    }
  }

  VecMin(x, &tttt, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "KSP Max %d %le\n", tttt, norm);

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==tttt) {
	  PetscPrintf(PETSC_COMM_SELF, "KspMMin %d %d %d %d %le\n", Gidx(mi,j,k,user),mi,j, k, norm);
	}
      }
    }
  }

  VecDestroy(&x);

  return 0;
}

/* PetscInt lidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user) */
/* { */
/*   DMDALocalInfo	info = user->info; */


/*   PetscInt	gxs, gxe, gys, gye, gzs, gze; */

/*   gxs = info.gxs; gxe = gxs + info.gxm; */
/*   gys = info.gys; gye = gys + info.gym; */
/*   gzs = info.gzs; gze = gzs + info.gzm; */


/*   if (!(user->aotopetsc)) { */

/*     user->aotopetsc = PETSC_TRUE; */

/*     DMDAGetGlobalIndices(user->da, PETSC_NULL, &user->idx_from); */
/*   } */
/*   return (user->idx_from[(k-gzs) * (info.gxm*info.gym) + (j-gys)*(info.gxm) + (i-gxs)]); */
  
/* } */
//Mohsen May 2012
PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)

{
  PetscInt nidx;
  DMDALocalInfo	info = user->info;

  PetscInt	mx = info.mx, my = info.my;
  
  AO ao;
  DMDAGetAO(user->da, &ao);
  nidx=i+j*mx+k*mx*my;
  
  AOApplicationToPetsc(ao,1,&nidx);
  
  return (nidx);
}
PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user)
{
  if ((user->isc)) {
    *ih = i;
  }
  else {
    *ih = 2 * i;
  }

  if ((user->jsc)) {
    *jh = j;
  }
  else {
    *jh = 2 * j;
  }

  if ((user->ksc)) {
    *kh = k;
  }
  else {
    *kh = 2 * k;
  }

  return 0;
}

PetscErrorCode MyRestriction(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, (void**)&user);

  
  DM	da = user->da;

  DM	da_f = *user->da_f;

  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  //  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert;
  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  DMDAVecGetArray(da,   F, &f);

  Vec lX;
 
  DMCreateLocalVector(da_f, &lX);
  DMGlobalToLocalBegin(da_f, X, INSERT_VALUES, lX);
  DMGlobalToLocalEnd(da_f, X, INSERT_VALUES, lX);  
  DMDAVecGetArray(da_f, lX, &x);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal ***nvert_f;
  DMDAVecGetArray(da_f, user->user_f->lNvert, &nvert_f);

  if ((user->isc)) ia = 0;
  else ia = 1;

  if ((user->jsc)) ja = 0;
  else ja = 1;

  if ((user->ksc)) ka = 0;
  else ka = 1;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (k==0) {
	  f[k][j][i] = 0.;
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;
	}
	else if (j==0) {
	  f[k][j][i] = 0.;
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;
	}
	else if (i==0) {
	  f[k][j][i] = 0.;
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;
	}
	else {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user);
	  f[k][j][i] = 0.125 *
	    (x[kh   ][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih   ]) +
	     x[kh   ][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih-ia]) +
	     x[kh   ][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih   ]) +
	     x[kh-ka][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih   ]) +
	     x[kh   ][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih-ia]) +
	     x[kh-ka][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih   ]) +
	     x[kh-ka][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih-ia]) +
	     x[kh-ka][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih-ia]));



	  if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;
	}
      }
    }
  }


  DMDAVecRestoreArray(da_f, user->user_f->lNvert, &nvert_f);

  DMDAVecRestoreArray(da_f, lX, &x);
  VecDestroy(&lX);
 
  DMDAVecRestoreArray(da,   F,  &f);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);


  return 0;
}

#define GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user) \
  if ((user->isc)) { \
    ic = i; \
    ia = 0; \
  } \
  else { \
    ic = (i+1) / 2; \
    ia = (i - 2 * (ic)) == 0 ? 1 : -1; \
    if (i==1 || i==mx-2) ia = 0; \
  }\
  if ((user->jsc)) { \
    jc = j; \
    ja = 0; \
  } \
  else { \
    jc = (j+1) / 2; \
    ja = (j - 2 * (jc)) == 0 ? 1 : -1; \
    if (j==1 || j==my-2) ja = 0; \
  } \
  if ((user->ksc)) { \
    kc = k; \
    ka = 0; \
  } \
  else { \
    kc = (k+1) / 2; \
    ka = (k - 2 * (kc)) == 0 ? 1 : -1; \
    if (k==1 || k==mz-2) ka = 0; \
  } \
  if (ka==-1 && nvert_c[kc-1][jc][ic] > 0.1) ka = 0; \
  else if (ka==1 && nvert_c[kc+1][jc][ic] > 0.1) ka = 0; \
  if (ja==-1 && nvert_c[kc][jc-1][ic] > 0.1) ja = 0; \
  else if (ja==1 && nvert_c[kc][jc+1][ic] > 0.1) ja = 0; \
  if (ia==-1 && nvert_c[kc][jc][ic-1] > 0.1) ia = 0; \
  else if (ia==1 && nvert_c[kc][jc][ic+1] > 0.1) ia = 0;



PetscErrorCode MyInterpolation(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, (void**)&user);

  
  
  DM	da = user->da;

  DM	da_c = *user->da_c;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert, ***nvert_c;
  PetscInt i, j, k, ic, jc, kc, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


  DMDAVecGetArray(da,   F, &f);


  Vec lX;
  DMCreateLocalVector(da_c, &lX);
 
  DMGlobalToLocalBegin(da_c, X, INSERT_VALUES, lX);
  DMGlobalToLocalEnd(da_c, X, INSERT_VALUES, lX);  
  DMDAVecGetArray(da_c, lX, &x);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da_c, *(user->lNvert_c), &nvert_c);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  f[k][j][i] = (x[kc   ][jc   ][ic   ] * 9 +
			x[kc   ][jc+ja][ic   ] * 3 +
			x[kc   ][jc   ][ic+ia] * 3 +
			x[kc   ][jc+ja][ic+ia]) * 3./64. +
	    (x[kc+ka][jc   ][ic   ] * 9 +
	     x[kc+ka][jc+ja][ic   ] * 3 +
	     x[kc+ka][jc   ][ic+ia] * 3 +
	     x[kc+ka][jc+ja][ic+ia]) /64.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {

	if (i==0) {
	  f[k][j][i] = 0.;//-f[k][j][i+1];
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;//-f[k][j][i-1];
	}
	else if (j==0) {
	  f[k][j][i] = 0.;//-f[k][j+1][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;//-f[k][j-1][i];
	}
	else if (k==0) {
	  f[k][j][i] = 0.;//-f[k+1][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;//-f[k-1][j][i];
	}
	if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;

      }
    }
  }

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da_c, *(user->lNvert_c), &nvert_c);

  DMDAVecRestoreArray(da_c, lX, &x);
 
  VecDestroy(&lX);
  DMDAVecRestoreArray(da,   F,  &f);



  return 0;
 
}

PetscErrorCode PoissonNullSpaceFunction(MatNullSpace nullsp,Vec X, void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;

  DM da = user->da;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	***x, ***nvert;
  PetscInt	i, j, k;

/*   /\* First remove a constant from the Vec field X *\/ */


  /* Then apply boundary conditions */
  DMDAVecGetArray(da, X, &x);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal lsum, sum;
  PetscReal  lnum, num;

  // if (user->multinullspace) PetscPrintf(PETSC_COMM_WORLD, "MultiNullSpace!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  if (!user->multinullspace) {
    lsum = 0;
    lnum = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }

    MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lnum,&num,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   /*  PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD); */
    sum = sum / (-1.0 * num);

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    x[k][j][i] +=sum;
	  }
	}
      }
    }
  }
  else {
    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lnum,&num,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   /*  PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD); */
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }

    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lnum,&num,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   /*  PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD); */
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }
			   
  } //if multinullspace
  if (zs == 0) {
    k = 0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k = mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ys == 0) {
    j = 0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ye == my) {
    j = my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xs == 0) {
    i = 0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xe == mx) {
    i = mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert[k][j][i] > 0.1)
	  x[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, X, &x);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  return 0;
}

PetscErrorCode mymatmultadd(Mat mat, Vec v1, Vec v2)
{

  Vec vt;
  VecDuplicate(v2, &vt);
  MatMult(mat, v1, vt);
  VecAYPX(v2, 1., vt);
  VecDestroy(&vt);
  return(0);
}


#define CP  0
#define EP  1
#define WP  2
#define NP  3
#define SP  4
#define TP  5
#define BP  6
#define NE  7
#define SE  8
#define NW  9
#define SW  10
#define TN  11
#define BN  12
#define TS  13
#define BS  14
#define TE  15
#define BE  16
#define TW  17
#define BW  18




PetscErrorCode PoissonLHSNew(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt      IM=user->IM, JM=user->JM, KM=user->KM;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  PetscInt      gxs, gxe, gys, gye, gzs, gze;

  Vec		G11, G12, G13, G21, G22, G23, G31, G32, G33;
  PetscReal	***g11, ***g12, ***g13, ***g21, ***g22, ***g23;
  PetscReal	***g31, ***g32, ***g33;

  PetscReal	***nvert;
  PetscScalar	vol[19];
  PetscInt	idx[19], row;

  PetscInt	i, j, k, N;
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (!user->assignedA) {
    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(user->Phi, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 19, PETSC_NULL, 19, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  MatZeroEntries(user->A);

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

  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(da, user->lNvert, &nvert);
 

  VecDuplicate(user->lAj, &G11);
  VecDuplicate(user->lAj, &G12);
  VecDuplicate(user->lAj, &G13);
  VecDuplicate(user->lAj, &G21);
  VecDuplicate(user->lAj, &G22);
  VecDuplicate(user->lAj, &G23);
  VecDuplicate(user->lAj, &G31);
  VecDuplicate(user->lAj, &G32);
  VecDuplicate(user->lAj, &G33);

  DMDAVecGetArray(da, G11, &g11);
  DMDAVecGetArray(da, G12, &g12);
  DMDAVecGetArray(da, G13, &g13);
  DMDAVecGetArray(da, G21, &g21);
  DMDAVecGetArray(da, G22, &g22);
  DMDAVecGetArray(da, G23, &g23);
  DMDAVecGetArray(da, G31, &g31);
  DMDAVecGetArray(da, G32, &g32);
  DMDAVecGetArray(da, G33, &g33);

  for (k=gzs; k<gze; k++) {
    for (j=gys; j<gye; j++) {
      for (i=gxs; i<gxe; i++) {
	if(i>-1 && j>-1 && k>-1 && i<IM+1 && j<JM+1 && k<KM+1){ //Mohsen April 2012
	g11[k][j][i] = (icsi[k][j][i].x * icsi[k][j][i].x +
			icsi[k][j][i].y * icsi[k][j][i].y +
			icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g12[k][j][i] = (ieta[k][j][i].x * icsi[k][j][i].x +
			ieta[k][j][i].y * icsi[k][j][i].y +
			ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g13[k][j][i] = (izet[k][j][i].x * icsi[k][j][i].x +
			izet[k][j][i].y * icsi[k][j][i].y +
			izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
					  			 
	g21[k][j][i] = (jcsi[k][j][i].x * jeta[k][j][i].x +
			jcsi[k][j][i].y * jeta[k][j][i].y +
			jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g22[k][j][i] = (jeta[k][j][i].x * jeta[k][j][i].x +
			jeta[k][j][i].y * jeta[k][j][i].y +
			jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g23[k][j][i] = (jzet[k][j][i].x * jeta[k][j][i].x +
			jzet[k][j][i].y * jeta[k][j][i].y +
			jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
					  			 
	g31[k][j][i] = (kcsi[k][j][i].x * kzet[k][j][i].x +
			kcsi[k][j][i].y * kzet[k][j][i].y +
			kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g32[k][j][i] = (keta[k][j][i].x * kzet[k][j][i].x +
			keta[k][j][i].y * kzet[k][j][i].y +
			keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g33[k][j][i] = (kzet[k][j][i].x * kzet[k][j][i].x +
			kzet[k][j][i].y * kzet[k][j][i].y +
			kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	}	
      }
    }
  }

  PetscInt m;
  PetscReal threshold;
  if (user->thislevel == 2) {
    threshold = 0.1;
  }
  else {
    threshold = 0.1;
  }
  //Mohsen April 2012 

  PetscInt x_str,x_end,y_str,y_end,z_str,z_end;

  if (user->bctype[0]==7) {
    x_end=mx-1;
    x_str=0;
  }
  if (user->bctype[2]==7) {
    y_end=my-1;
    y_str=0;
  }
  if (user->bctype[4]==7) {
    z_end=mz-1;
    z_str=0;
  }

  if (user->bctype[0]!=7){
    x_end=mx-2;
    x_str=1;
  }
  if (user->bctype[2]!=7){
    y_end=my-2;
    y_str=1;
  }
  if (user->bctype[4]!=7){
    z_end=mz-2;
    z_str=1;
  }
 
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  vol[CP] = 1.;	idx[CP] = Gidx(i, j, k, user);
	  MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	}
	else {
	  if (nvert[k][j][i] > 0.1) { // i, j, k is not a fluid point
	    vol[CP] = 1.; idx[CP] = Gidx(i, j, k, user);
	    MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	  }
	  else { // i, j, k is a fluid point
	    for (m=0; m<19; m++) {
	      vol[m] = 0.;
	    }
	    /* Contribution from i+1 - i */
	    if (nvert[k][j][i+1] < threshold && i != x_end) { // i+1, j, k is a fluid point
	      /* dpdc{i} = (p_{i+1} - p_{i}) * g11_{i} */
	      vol[CP] -= g11[k][j][i]; //i, j, k
	      vol[EP] += g11[k][j][i]; // i+1, j, k

	      /* dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1})
	       * 0.25 * g12[k][j][i] */
	      if ((j == my-2 && user->bctype[2]!=7) || nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] += g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] += g12[k][j][i] * 0.5; //i+1, j, k
		  vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
		  vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
		  vol[CP] += g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] += g12[k][j][i] * 0.5; //i+1, j, k
		  vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
		  vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5; //i, j+1, k
		  vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
		  vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5; //i, j+1, k
		  vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
		  vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
		}
	      }
	      else {
		vol[NP] += g12[k][j][i] * 0.25; // i, j+1, k
		vol[NE] += g12[k][j][i] * 0.25; // i+1, j+1, k
		vol[SP] -= g12[k][j][i] * 0.25; // i, j-1, k
		vol[SE] -= g12[k][j][i] * 0.25; // i+1, j-1, k
	      }
	      
	      /* dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] += g13[k][j][i] * 0.5; // i, j, k
		  vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
		  vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
		  vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
		  vol[CP] += g13[k][j][i] * 0.5; // i, j, k
		  vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
		  vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
		  vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7) || nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
		  vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
		  vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
		  vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
		  vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
		  vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
		  vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
		}
	      }
	      else {
		vol[TP] += g13[k][j][i] * 0.25; //i, j, k+1
		vol[TE] += g13[k][j][i] * 0.25; //i+1, j, k+1
		vol[BP] -= g13[k][j][i] * 0.25; //i, j, k-1
		vol[BE] -= g13[k][j][i] * 0.25; //i+1, j, k-1
	      }
	    }  // end of i+1 - i

	    /* Contribution from i - i-1 */
	    if (nvert[k][j][i-1] < threshold && i != x_str) { // i-1, j, k is a fluid point
	      /* -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i} */
	      vol[CP] -= g11[k][j][i-1];  //i, j, k
	      vol[WP] += g11[k][j][i-1];  //i-1, j, k

	      /* -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1})
	       * 0.25 * g12[k][j][i-1] */
	      if ((j == my-2 && user->bctype[2]!=7) || nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
		  vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
		  vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
		  vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
		  vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; // i, j+1, k
		  vol[NW] -= g12[k][j][i-1] * 0.5; // i-1, j+1, k
		  vol[CP] += g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g12[k][j][i-1] * 0.5; // i-1, j, k
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; // i, j+1, k
		  vol[NW] -= g12[k][j][i-1] * 0.5; // i-1, j+1, k
		  vol[CP] += g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g12[k][j][i-1] * 0.5; // i-1, j, k
		}
	      }
	      else {
		vol[NP] -= g12[k][j][i-1] * 0.25; // i, j+1, k
		vol[NW] -= g12[k][j][i-1] * 0.25; //i-1, j+1, k
		vol[SP] += g12[k][j][i-1] * 0.25; // i, j-1, k
		vol[SW] += g12[k][j][i-1] * 0.25; // i-1, j-1, k
	      }

	      /* -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
		  vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
		  vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
		  vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
		  vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7) || nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; // i, j, k+1
		  vol[TW] -= g13[k][j][i-1] * 0.5; //i-1, j, k+1
		  vol[CP] += g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g13[k][j][i-1] * 0.5; //i-1, j, k
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; // i, j, k+1
		  vol[TW] -= g13[k][j][i-1] * 0.5; //i-1, j, k+1
		  vol[CP] += g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g13[k][j][i-1] * 0.5; //i-1, j, k
		}
	      }
	      else {
		vol[TP] -= g13[k][j][i-1] * 0.25;  // i, j, k+1
		vol[TW] -= g13[k][j][i-1] * 0.25; // i-1, j, k+1
		vol[BP] += g13[k][j][i-1] * 0.25;  // i, j, k-1
		vol[BW] += g13[k][j][i-1] * 0.25; // i-1, j, k-1
	      }
	    } // end of i - i-1

	    /* Contribution from j+1 - j */
	    if (nvert[k][j+1][i] < threshold && j != y_end) {
	      /* dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) *
		 0.25 */
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] += g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
		  vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
		  vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
		  vol[CP] += g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
		  vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
		  vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7) || nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; // i+1, j, k
		  vol[NE] += g21[k][j][i] * 0.5; // i+1, j+1, k
		  vol[CP] -= g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] -= g21[k][j][i] * 0.5; // i, j+1, k
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 &&  nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; // i+1, j, k
		  vol[NE] += g21[k][j][i] * 0.5; // i+1, j+1, k
		  vol[CP] -= g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] -= g21[k][j][i] * 0.5; // i, j+1, k
		}
	      }
	      else {
		vol[EP] += g21[k][j][i] * 0.25; //i+1, j, k
		vol[NE] += g21[k][j][i] * 0.25; //i+1, j+1, k
		vol[WP] -= g21[k][j][i] * 0.25; //i-1, j, k
		vol[NW] -= g21[k][j][i] * 0.25; //i-1, j+1, k
	      }
	     
	      /* dpde{j} = (p{j+1} - p{j}) * g22[k][j][i] */
	      vol[CP] -= g22[k][j][i];
	      vol[NP] += g22[k][j][i];

	      /* dpdz{j} = (p{j, k+1}+p{j+1,k+1} - p{j,k-1}-p{j+1,k-1}) *0.25*/
	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] += g23[k][j][i] * 0.5; //i,j,k
		  vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
		  vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
		  vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
		}
	      }
	      else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[CP] += g23[k][j][i] * 0.5; //i,j,k
		  vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
		  vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
		  vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; //i, j, k+1
		  vol[TN] += g23[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g23[k][j][i] * 0.5;//i, j, k
		  vol[NP] -= g23[k][j][i] * 0.5;//i, j+1, k
		}
	      }
	      else if ((k == 1 || k==mz-2 ) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; //i, j, k+1
		  vol[TN] += g23[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g23[k][j][i] * 0.5;//i, j, k
		  vol[NP] -= g23[k][j][i] * 0.5;//i, j+1, k
		}
	      }
	      else {
		vol[TP] += g23[k][j][i] * 0.25; // i, j, k+1
		vol[TN] += g23[k][j][i] * 0.25; // i, j+1, k+1
		vol[BP] -= g23[k][j][i] * 0.25; // i, j, k-1
		vol[BN] -= g23[k][j][i] * 0.25; // i, j+1, k-1
	      }
	    } // End of j+1 - j

	    /* Contribution j - j-1 */
	    if (nvert[k][j-1][i] < threshold && j != y_str) {
	      /* -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) *
		 0.25 * g21[k][j-1][i] */
	      if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
		  vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
		  vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
		  vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
		}
	      }
	      else  if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1) {
		  vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
		  vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
		  vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
		  vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5;//i+1, j, k
		  vol[SE] -= g21[k][j-1][i] * 0.5;//i+1, j-1, k
		  vol[CP] += g21[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g21[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5;//i+1, j, k
		  vol[SE] -= g21[k][j-1][i] * 0.5;//i+1, j-1, k
		  vol[CP] += g21[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g21[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[EP] -= g21[k][j-1][i] * 0.25;// i+1, j, k
		vol[SE] -= g21[k][j-1][i] * 0.25;// i+1, j-1, k
		vol[WP] += g21[k][j-1][i] * 0.25;// i-1, j, k
		vol[SW] += g21[k][j-1][i] * 0.25;// i-1, j-1, k
	      }
	      
	      /* -dpde{j-1} = -(p{j} - p{j-1}) * g22[k][j-1][i] */
	      vol[CP] -= g22[k][j-1][i];
	      vol[SP] += g22[k][j-1][i];

	      /* -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) *
		 0.25 * g23[k][j-1][i] */
	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
		  vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
		  vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 ) {
		  vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
		  vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
		  vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
		}
	      }
	     
	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5;//i, j, k+1
		  vol[TS] -= g23[k][j-1][i] * 0.5;//i, j-1, k+1
		  vol[CP] += g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g23[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5;//i, j, k+1
		  vol[TS] -= g23[k][j-1][i] * 0.5;//i, j-1, k+1
		  vol[CP] += g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g23[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[TP] -= g23[k][j-1][i] * 0.25;//i, j, k+1
		vol[TS] -= g23[k][j-1][i] * 0.25;//i, j-1, k+1
		vol[BP] += g23[k][j-1][i] * 0.25;//i, j, k-1
		vol[BS] += g23[k][j-1][i] * 0.25;//i, j-1, k-1
	      }
	    } // End of j - j-1

	    /* contribution from k+1 - k */
	    if (nvert[k+1][j][i] < threshold && k != z_end) {
	      /* dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) *
		 0.25 * g31[k][j][i] */
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] += g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
		  vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
		  vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
		  vol[CP] += g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
		  vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
		  vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5;//i+1, j, k
		  vol[TE] += g31[k][j][i] * 0.5;//i+1, j, k+1
		  vol[CP] -= g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g31[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5;//i+1, j, k
		  vol[TE] += g31[k][j][i] * 0.5;//i+1, j, k+1
		  vol[CP] -= g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g31[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[EP] += g31[k][j][i] * 0.25;//i+1, j, k
		vol[TE] += g31[k][j][i] * 0.25;//i+1, j, k+1
		vol[WP] -= g31[k][j][i] * 0.25;//i-1, j, k
		vol[TW] -= g31[k][j][i] * 0.25;//i-1, j, k+1
	      }

	      /* dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 
		 0.25 * g32[k][j][i] */
	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] += g32[k][j][i] * 0.5;//i, j,k
		  vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
		  vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
		  vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[CP] += g32[k][j][i] * 0.5;//i, j,k
		  vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
		  vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
		  vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5;//i, j+1, k
		  vol[TN] += g32[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g32[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g32[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5;//i, j+1, k
		  vol[TN] += g32[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g32[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g32[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[NP] += g32[k][j][i] * 0.25;//i, j+1, k
		vol[TN] += g32[k][j][i] * 0.25;//i, j+1, k+1
		vol[SP] -= g32[k][j][i] * 0.25;//i, j-1, k
		vol[TS] -= g32[k][j][i] * 0.25;//i, j-1, k+1
	      }

	      /* dpdz{k} = p{k+1} - p{k} */
	      vol[CP] -= g33[k][j][i]; //i, j, k
	      vol[TP] += g33[k][j][i]; //i, j, k+1
	    } // End of k+1 - k

	    /* Contribution from k - k-1 */
	    if (nvert[k-1][j][i] < threshold && k != z_str) {
	      /* -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) *
		 0.25 * g31[k-1][j][i] */
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
		  vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
		  vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1) {
		  vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
		  vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
		  vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5;//i+1, j, k
		  vol[BE] -= g31[k-1][j][i] * 0.5;//i+1, j, k-1
		  vol[CP] += g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g31[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5;//i+1, j, k
		  vol[BE] -= g31[k-1][j][i] * 0.5;//i+1, j, k-1
		  vol[CP] += g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g31[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[EP] -= g31[k-1][j][i] * 0.25;//i+1, j, k
		vol[BE] -= g31[k-1][j][i] * 0.25;//i+1, j, k-1
		vol[WP] += g31[k-1][j][i] * 0.25;//i-1, j, k
		vol[BW] += g31[k-1][j][i] * 0.25;//i-1, j, k-1
	      }
	      
	      /* -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) * 
		 0.25 * g32[k-1][j][i] */
	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
		  vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
		  vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
		  vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1) {
		  vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
		  vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
		  vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
		  vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5;//i, j+1, k
		  vol[BN] -= g32[k-1][j][i] * 0.5;//i, j+1, k-1
		  vol[CP] += g32[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g32[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5;//i, j+1, k
		  vol[BN] -= g32[k-1][j][i] * 0.5;//i, j+1, k-1
		  vol[CP] += g32[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g32[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[NP] -= g32[k-1][j][i] * 0.25;//i, j+1, k
		vol[BN] -= g32[k-1][j][i] * 0.25;//i, j+1, k-1
		vol[SP] += g32[k-1][j][i] * 0.25;//i, j-1, k
		vol[BS] += g32[k-1][j][i] * 0.25;//i, j-1, k-1
	      }
	      
	      /* -dpdz{k-1} = -(p{k} - p{k-1}) * g33[k-1][j][i] */
	      vol[CP] -= g33[k-1][j][i]; // i, j, k
	      vol[BP] += g33[k-1][j][i]; //i, j, k-1
	    } // End of k - k-1
	    for (m=0; m<19; m++) {
	      vol[m] *= -aj[k][j][i];
	    }
	    //Mohsen April 2012
	    idx[CP] = Gidx(i  , j  , k  , user);
	    if (user->bctype[0]==7 && i==mx-2)
	      idx[EP] = Gidx(1, j  , k  , user);
	    else
	      idx[EP] = Gidx(i+1, j  , k  , user);
	    if (user->bctype[0]==7 && i==1)
	      idx[WP] = Gidx(mx-2, j  , k  , user);
	    else
	      idx[WP] = Gidx(i-1, j  , k  , user);
	    if (user->bctype[2]==7 && j==my-2)
	      idx[NP] = Gidx(i, 1  , k  , user);
	    else
	      idx[NP] = Gidx(i  , j+1, k  , user);
	    if (user->bctype[2]==7 && j==1)
	      idx[SP] = Gidx(i, my-2  , k  , user);
	    else
	      idx[SP] = Gidx(i  , j-1, k  , user);
	    
	    if (user->bctype[4]==7 && k==mz-2){
	      idx[TP] = Gidx(i  , j  , 1, user);
	    } else
	      idx[TP] = Gidx(i  , j  , k+1, user);
   	    if (user->bctype[4]==7 && k==1)
	      idx[BP] = Gidx(i  , j  , mz-2, user);
	    else
	      idx[BP] = Gidx(i  , j  , k-1, user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==mx-2 && j==my-2)
	      idx[NE] = Gidx(1, 1, k  , user);
	    else if (user->bctype[0]==7 && i==mx-2)
	      idx[NE] = Gidx(1, j+1, k  , user);
	    else if (user->bctype[2]==7 && j==my-2)
	      idx[NE] = Gidx(i+1, 1, k  , user);
	    else
	      idx[NE] = Gidx(i+1, j+1, k  , user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==mx-2 && j==1)
	      idx[SE] = Gidx(1, my-2, k  , user);
	    else if (user->bctype[0]==7 && i==mx-2)
	      idx[SE] = Gidx(1, j-1, k  , user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[SE] = Gidx(i+1, my-2, k  , user);
	    else
	      idx[SE] = Gidx(i+1, j-1, k  , user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==1 && j==my-2)
	      idx[NW] = Gidx(mx-2, 1, k  , user);
	    else if (user->bctype[0]==7 && i==1)
	      idx[NW] = Gidx(mx-2, j+1, k  , user);
	    else if (user->bctype[2]==7 && j==my-2)
	      idx[NW] = Gidx(i-1, 1, k  , user);
	    else
	      idx[NW] = Gidx(i-1, j+1, k  , user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==1 && j==1)
	      idx[SW] = Gidx(mx-2, my-2, k  , user);
	    else if (user->bctype[0]==7 && i==1)
	      idx[SW] = Gidx(mx-2, j-1, k  , user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[SW] = Gidx(i-1, my-2, k  , user);
	    else
	      idx[SW] = Gidx(i-1, j-1, k  , user);
	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==my-2 && k==mz-2)
	      idx[TN] = Gidx(i  , 1, 1, user);
	    else if (user->bctype[2]==7 && j==my-2)
	      idx[TN] = Gidx(i  , 1, k+1, user);
	    else if (user->bctype[4]==7 && k==mz-2)
	      idx[TN] = Gidx(i  , j+1, 1, user);
	    else
	      idx[TN] = Gidx(i  , j+1, k+1, user);
	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==my-2 && k==1)
	      idx[BN] = Gidx(i  , 1, mz-2 , user);
	    else if(user->bctype[2]==7 && j==my-2)
	      idx[BN] = Gidx(i  , 1, k-1, user);
	    else if (user->bctype[4]==7 && k==1)
	      idx[BN] = Gidx(i  , j+1, mz-2, user);
	    else
	      idx[BN] = Gidx(i  , j+1, k-1, user);
	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==1 && k==mz-2)
	      idx[TS] = Gidx(i  , my-2, 1, user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[TS] = Gidx(i  , my-2, k+1, user);
	    else if (user->bctype[4]==7 && k==mz-2)
	      idx[TS] = Gidx(i  , j-1, 1, user);
	    else
	      idx[TS] = Gidx(i  , j-1, k+1, user);
       	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==1 && k==1)
	      idx[BS] = Gidx(i  , my-2, mz-2, user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[BS] = Gidx(i  , my-2, k-1, user);
	    else if (user->bctype[4]==7 && k==1)
	      idx[BS] = Gidx(i  , j-1, mz-2, user);
	    else
	      idx[BS] = Gidx(i  , j-1, k-1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==mx-2 && k==mz-2)
	      idx[TE] = Gidx(1, j  , 1, user);
	    else if(user->bctype[0]==7 && i==mx-2)
	      idx[TE] = Gidx(1, j  , k+1, user);
	    else if(user->bctype[4]==7 && k==mz-2)
	      idx[TE] = Gidx(i+1, j  , 1, user);
	    else
	      idx[TE] = Gidx(i+1, j  , k+1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==mx-2 && k==1)
	      idx[BE] = Gidx(1, j  , mz-2, user);
	    else if(user->bctype[0]==7 && i==mx-2)
	      idx[BE] = Gidx(1, j  , k-1, user);
	    else if(user->bctype[4]==7 && k==1)
	      idx[BE] = Gidx(i+1, j  , mz-2, user);
	    else
	      idx[BE] = Gidx(i+1, j  , k-1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==1 && k==mz-2)
	      idx[TW] = Gidx(mx-2, j  , 1, user);
	    else if(user->bctype[0]==7 && i==1)
	      idx[TW] = Gidx(mx-2, j  , k+1, user);
	    else if (user->bctype[4]==7 && k==mz-2)
	      idx[TW] = Gidx(i-1, j  , 1, user);
	    else
	      idx[TW] = Gidx(i-1, j  , k+1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==1 && k==1)
	      idx[BW] = Gidx(mx-2, j  , mz-2, user);
	    else if (user->bctype[0]==7 && i==1)
	      idx[BW] = Gidx(mx-2, j  , k-1, user);
	    else if (user->bctype[4]==7 && k==1)
	      idx[BW] = Gidx(i-1, j  , mz-2, user);
	    else
	      idx[BW] = Gidx(i-1, j  , k-1, user);
	    MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);

	  } // End of fluid point
	} // End of interial points
      }
    }
  }
  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  DMDAVecRestoreArray(da, G11, &g11);
  DMDAVecRestoreArray(da, G12, &g12);
  DMDAVecRestoreArray(da, G13, &g13);
  DMDAVecRestoreArray(da, G21, &g21);
  DMDAVecRestoreArray(da, G22, &g22);
  DMDAVecRestoreArray(da, G23, &g23);
  DMDAVecRestoreArray(da, G31, &g31);
  DMDAVecRestoreArray(da, G32, &g32);
  DMDAVecRestoreArray(da, G33, &g33);
  
  VecDestroy(&G11);
  VecDestroy(&G12);
  VecDestroy(&G13);
  VecDestroy(&G21);
  VecDestroy(&G22);
  VecDestroy(&G23);
  VecDestroy(&G31);
  VecDestroy(&G32);
  VecDestroy(&G33);


  //  VecCopy(user->Phi, user->P);
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

  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  return 0;
}

PetscErrorCode PoissonRHS(UserCtx *user, Vec B)
{
  DMDALocalInfo info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
 
  PetscInt      i, j, k;
  PetscReal	***nvert, ***aj, ***rb, dt = user->dt;
  struct Components{
    PetscReal x;
    PetscReal y;
    PetscReal z;
  } *** ucont;
       

  DMDAVecGetArray(user->da, B, &rb);
  DMDAVecGetArray(user->fda, user->lUcont, &ucont);
  DMDAVecGetArray(user->da, user->lNvert, &nvert);
  DMDAVecGetArray(user->da, user->lAj, &aj);

  for (k=zs; k<ze; k++) { 
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {

	if (i==0 || i==mx-1 || j==0 || j==my-1 ||  k==0 || k==mz-1) {
	  rb[k][j][i] = 0.;
	}
	else if (nvert[k][j][i] > 0.1) {
	  rb[k][j][i] = 0;
	}
	else {
	  rb[k][j][i] = -(ucont[k][j][i].x - ucont[k][j][i-1].x +
			  ucont[k][j][i].y - ucont[k][j-1][i].y +
			  ucont[k][j][i].z - ucont[k-1][j][i].z) / dt
 	    * aj[k][j][i] / user->st * COEF_TIME_ACCURACY;
	 
	}
      }
    }
  }



  PetscReal lsum=0., sum=0.;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	
	lsum += rb[k][j][i] / aj[k][j][i]* dt/COEF_TIME_ACCURACY;

      }
    }
  }
  
  MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  // PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "Summation RHS %le\n", sum);
	
  DMDAVecRestoreArray(user->fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
  DMDAVecRestoreArray(user->da, user->lAj, &aj);
  DMDAVecRestoreArray(user->da, B, &rb);
 
  return 0;
}

PetscErrorCode FullyBlocked(UserCtx *user)
{
  DM da = user->da;
  Vec nNvert;
  DMDALocalInfo info = user->info;

  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;

  PetscInt *KSKE = user->KSKE;
  PetscReal ***nvert;
  PetscBool *Blocked;

  DMDACreateNaturalVector(da, &nNvert);
  DMDAGlobalToNaturalBegin(da, user->Nvert, INSERT_VALUES, nNvert);
  DMDAGlobalToNaturalEnd(da, user->Nvert, INSERT_VALUES, nNvert);

  VecScatter ctx;
  Vec Zvert;
  VecScatterCreateToZero(nNvert, &ctx, &Zvert);

  VecScatterBegin(ctx, nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD);

  VecScatterDestroy(&ctx);
  VecDestroy(&nNvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {

    VecGetArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    PetscMalloc(mx*my*sizeof(PetscBool), &Blocked);
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	Blocked[j*mx+i] = PETSC_FALSE;
	for (k=0; k<mz; k++) {
	  if (nvert[k][j][i] > 0.1) {
	    if (!Blocked[j*mx+i]) {
	      KSKE[2*(j*mx+i)] = k;
	      Blocked[j*mx+i] = PETSC_TRUE;
	    }
	    else {
	      KSKE[2*(j*mx+i)] = PetscMin(KSKE[2*(j*mx+i)], k);
	    }
	  }
	}
      }
    }


    user->multinullspace = PETSC_TRUE;
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	if (!Blocked[j*mx+i]) {
	  user->multinullspace = PETSC_FALSE;
	  break;
	}
      }
    }
    PetscFree(Blocked);
    VecRestoreArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);

    }
  }
  else {
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }



  VecDestroy(&Zvert);
  return 0;
}

PetscErrorCode MyNvertRestriction(UserCtx *user_h, UserCtx *user_c)
{
  //  DA		da = user_c->da, fda = user_c->fda;



  DMDALocalInfo	info = user_c->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i,j,k;
  PetscInt ih, jh, kh, ia, ja, ka;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert, ***nvert_h;

  DMDAVecGetArray(user_h->da, user_h->lNvert, &nvert_h);
  DMDAVecGetArray(user_c->da, user_c->Nvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if ((user_c->isc)) ia = 0;
  else ia = 1;

  if ((user_c->jsc)) ja = 0;
  else ja = 1;

  if ((user_c->ksc)) ka = 0;
  else ka = 1;

  VecSet(user_c->Nvert, 0.);
  if (user_c->thislevel > 0) {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
	  if (nvert_h[kh   ][jh   ][ih   ] *
	      nvert_h[kh   ][jh   ][ih-ia] *
	      nvert_h[kh   ][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh   ][ih   ] *
	      nvert_h[kh   ][jh-ja][ih-ia] *
	      nvert_h[kh-ka][jh   ][ih-ia] *
	      nvert_h[kh-ka][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
	    nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
	  }
	}
      }
    }
  }
  else {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
	  if (nvert_h[kh   ][jh   ][ih   ] *
	      nvert_h[kh   ][jh   ][ih-ia] *
	      nvert_h[kh   ][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh   ][ih   ] *
	      nvert_h[kh   ][jh-ja][ih-ia] *
	      nvert_h[kh-ka][jh   ][ih-ia] *
	      nvert_h[kh-ka][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
	    nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
	  }
	}
      }
    }
  }
  DMDAVecRestoreArray(user_h->da, user_h->lNvert, &nvert_h);
  DMDAVecRestoreArray(user_c->da, user_c->Nvert, &nvert);

  DMGlobalToLocalBegin(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DMGlobalToLocalEnd(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DMDAVecGetArray(user_c->da, user_c->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j][i-1] > 1.1 &&
	      nvert[k][j+1][i] + nvert[k][j-1][i] > 1.1 &&
	      nvert[k+1][j][i] + nvert[k-1][j][i] > 1.1) {
	    nvert[k][j][i] = 1.;
	  }
	}
      }
    }
  }

  DMDAVecRestoreArray(user_c->da, user_c->lNvert, &nvert);
  DMLocalToGlobalBegin(user_c->da, user_c->lNvert, INSERT_VALUES, user_c->Nvert);
  DMLocalToGlobalEnd(user_c->da, user_c->lNvert, INSERT_VALUES, user_c->Nvert);
  return 0;
}

PetscErrorCode PoissonSolver_MG(UserMG *usermg)
{
  PetscInt l;
  PetscInt bi;

  MGCtx *mgctx = usermg->mgctx;

  KSP	mgksp, subksp;
  PC	mgpc, subpc;
  UserCtx	*user;

  PetscInt	m_c, m_f, M_c, M_f, flg;

  for (bi=0; bi<block_number; bi++) {

  //  MyNFaceInit(usermg);

  PetscPrintf(PETSC_COMM_WORLD, "PoissonMG NullSpace\n");

  if (immersed) {
    for (l=usermg->mglevels-1; l>0; l--) {
	mgctx[l].user[bi].multinullspace = PETSC_FALSE;
	MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]); 
    }
    /* At the coarsest level, check whether the grid is separated into sevearal 
       blocks by the immersed body */
    l = 0;
    user = mgctx[l].user;
    PetscMalloc(user[bi].info.mx*user[bi].info.my*2*sizeof(PetscInt), &user[bi].KSKE);
    FullyBlocked(&user[bi]);
  }

  l = usermg->mglevels-1;
  user = mgctx[l].user;

  VecDuplicate(user[bi].P, &user[bi].B);
  PetscReal ibm_Flux, ibm_Area;
    
  flg=immersed-1;
  if (MHV || LV|| aneurysm || LVAD) {
    if ((MHV>1 || LV || aneurysm || LVAD) && bi==0 )
      flg=1;
    else
      flg=0;
    VolumeFlux_rev(&user[bi], &ibm_Flux, &ibm_Area, flg);
  }

  VolumeFlux(&user[bi], &ibm_Flux, &ibm_Area, flg);
  PoissonRHS(&(user[bi]), user[bi].B);
   
  for (l=usermg->mglevels-1; l>=0; l--) {
    
  user = mgctx[l].user;
  PoissonLHSNew(&(user[bi]));
  }

  PetscErrorCode ierr;
    /* Create ksp for Multigrid */
  ierr=KSPCreate(PETSC_COMM_WORLD, &mgksp); CHKERRQ(ierr);
  ierr=KSPAppendOptionsPrefix(mgksp, "ps_"); CHKERRQ(ierr);

   /* Use multigrid as preconditioner */
  ierr=KSPGetPC(mgksp, &mgpc); CHKERRQ(ierr);
  ierr=PCSetType(mgpc, PCMG); CHKERRQ(ierr);
  /* Setup MG levels from usergm->mglevels */
  ierr=PCMGSetLevels(mgpc, usermg->mglevels, PETSC_NULL); CHKERRQ(ierr);
  /* V cycle */
  ierr=PCMGSetCycleType(mgpc, PC_MG_CYCLE_V); CHKERRQ(ierr);
  
  ierr=PCMGSetType(mgpc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);

  /*     PCMGSetType(mgpc, PC_MG_FULL); */
  
  /* Create Restriction and Interpolate schemes
     This is needed for all levels other than the coarsest one */
  for (l=usermg->mglevels-1; l>0; l--) {
    user = mgctx[l].user;
    m_c = (usermg->mgctx[l-1].user[bi].info.xm *
	   usermg->mgctx[l-1].user[bi].info.ym *
	   usermg->mgctx[l-1].user[bi].info.zm);
    
    m_f = (usermg->mgctx[l].user[bi].info.xm *
	   usermg->mgctx[l].user[bi].info.ym *
	   usermg->mgctx[l].user[bi].info.zm);
    
    M_c = (usermg->mgctx[l-1].user[bi].info.mx *
	   usermg->mgctx[l-1].user[bi].info.my *
	   usermg->mgctx[l-1].user[bi].info.mz);
    
    M_f = (usermg->mgctx[l].user[bi].info.mx *
	   usermg->mgctx[l].user[bi].info.my *
	   usermg->mgctx[l].user[bi].info.mz);
    
    ierr=MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, (void*)&mgctx[l-1].user[bi], &user[bi].MR); CHKERRQ(ierr);
    ierr=MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, (void*)&user[bi], &user[bi].MP); CHKERRQ(ierr);
    
    ierr=PCMGSetRestriction(mgpc, l, user[bi].MR); CHKERRQ(ierr);
    ierr=PCMGSetInterpolation(mgpc, l, user[bi].MP); CHKERRQ(ierr);
    
    /* Use subroutine MyRestriction and MyInterpolation for
       Mat * Vec operation */
    ierr=MatShellSetOperation(user[bi].MR, MATOP_MULT, (void(*)(void))MyRestriction); CHKERRQ(ierr);
    ierr=MatShellSetOperation(user[bi].MP, MATOP_MULT, (void(*)(void))MyInterpolation); CHKERRQ(ierr);
    
    ierr=MatShellSetOperation(user[bi].MR, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);CHKERRQ(ierr);
    ierr=MatShellSetOperation(user[bi].MP, MATOP_MULT_ADD,(void(*)(void))mymatmultadd); CHKERRQ(ierr);
    
  }
  

  if(visflg)   PetscPrintf(PETSC_COMM_WORLD, "PoissonMG Setup\n");

  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    
    if (l) { /* Not the coarset grid level */
      ierr=PCMGGetSmoother(mgpc, l, &subksp);CHKERRQ(ierr);
      /* Set the left hand side for every KSP at each grid level */
      ierr=KSPSetOperators(subksp, user[bi].A, user[bi].A); CHKERRQ(ierr);
      
      
      ierr=KSPSetFromOptions(subksp);CHKERRQ(ierr);
      KSPSetUp(subksp);
      KSPGetPC(subksp, &subpc);
    
    
      PCSetType(subpc, PCBJACOBI);
      PCSetUp(subpc);

      KSP *subsubksp;
      PC subsubpc;
      PetscInt abi, nlocal;
      
      PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp);
      for (abi = 0; abi<nlocal; abi++) {
	KSPGetPC(subsubksp[abi], &subsubpc);
	PCFactorSetShiftAmount(subsubpc, 1.e-10);
	//	PCFactorSetShiftType(subsubpc,MAT_SHIFT_NONZERO);
      }
      PCFactorSetShiftAmount(subpc, 1.e-10);
      // PCFactorSetShiftType(subpc,MAT_SHIFT_NONZERO);
    }
    else {  /* Coarsest grid */

      /* The default solver for the coarset grid is
	 KSPPreonly and PCLU.
	 One can choose other solvers, such as a KSP solver by changing the
	 following lines. */
      ierr=PCMGGetCoarseSolve(mgpc, &subksp);CHKERRQ(ierr);
     
      KSPGetPC(subksp, &subpc);
      ierr=KSPSetOperators(subksp, user[bi].A, user[bi].A);CHKERRQ(ierr);
    
      PCSetType(subpc, PCBJACOBI);
      // PCSetUp(subpc);
      // KSPSetFromOptions(subksp);
      KSPSetUp(subksp);
     
      
      KSP *subsubksp;
      PC subsubpc;
      PetscInt abi, nlocal;
      PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp);
      for (abi = 0; abi<nlocal; abi++) {
	KSPGetPC(subsubksp[abi], &subsubpc);
	PCFactorSetShiftAmount(subsubpc, 1.e-10);
	// PCFactorSetShiftType(subsubpc,MAT_SHIFT_NONZERO);
      }
      
      ierr=KSPSetTolerances(subksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 40);CHKERRQ(ierr);
      
     
    }

    /* The Poisson equation has Neumann boundary conditions, thus
       need to use NullSpace to solve the equation.
       We use the subroutine PoissonNullSpaceFunction to
       (1) reduce a constant value from the solution
       (2) Set the solution values at boundaries and blanking nodes to
       zero */
    ierr=MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp);CHKERRQ(ierr);
    ierr=MatNullSpaceSetFunction(user[bi].nullsp, PoissonNullSpaceFunction, &user[bi]);CHKERRQ(ierr);
    
    // ierr=KSPSetNullSpace(subksp, user[bi].nullsp);CHKERRQ(ierr);
    MatSetNullSpace( user[bi].A, user[bi].nullsp); //Mohsen Nov 2015
    //    ierr=PCMGSetResidual(mgpc, l, PCMGDefaultResidual, user[bi].A);CHKERRQ(ierr);
    ierr=PCMGSetResidual(mgpc, l, PCMGResidualDefault, user[bi].A);CHKERRQ(ierr);

    ierr=KSPSetUp(subksp);CHKERRQ(ierr);
    
    if (l<usermg->mglevels-1) {
      //  ierr=MatGetVecs(user[bi].A, &user[bi].R, PETSC_NULL);CHKERRQ(ierr);
      ierr=MatCreateVecs(user[bi].A, &user[bi].R, PETSC_NULL);CHKERRQ(ierr);
     
      ierr=PCMGSetRhs(mgpc, l, user[bi].R);CHKERRQ(ierr);
    }
    
  }
  
  
   if(visflg)  PetscPrintf(PETSC_COMM_WORLD, "PC for PoissonMG Setup\n");
  
  l = usermg->mglevels-1;
  user = mgctx[l].user;
  
  ierr=KSPSetOperators(mgksp, user[bi].A, user[bi].A);CHKERRQ(ierr);
  //  ierr=KSPSetNullSpace(mgksp, user[bi].nullsp);CHKERRQ(ierr);
  MatSetNullSpace( user[bi].A, user[bi].nullsp); //Mohsen Nov 2015
  ierr=KSPSetFromOptions(mgksp);CHKERRQ(ierr);
  ierr=KSPSetUp(mgksp);CHKERRQ(ierr);
  if(visflg)  PetscPrintf(PETSC_COMM_WORLD, "KSP Setup & solve for block %d\n",bi);
  
  
  ierr= KSPSolve(mgksp, user[bi].B, user[bi].Phi); CHKERRQ(ierr);

  DMGlobalToLocalBegin(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
  DMGlobalToLocalEnd(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);


  
  /* Release the allocated spaces */
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    MatNullSpaceDestroy(&user[bi].nullsp);

    MatDestroy(&user[bi].A);
    user[bi].assignedA = PETSC_FALSE;
    if (l) { /* Grid level > 0 (coarsest) */
      MatDestroy(&user[bi].MR);
      MatDestroy(&user[bi].MP);
    }
    else {
      if (immersed) PetscFree(user[bi].KSKE);
    }
  }

  
  
  for (l=0; l<usermg->mglevels-1; l++) {
    VecDestroy(&mgctx[l].user[bi].R);
  }
  
  KSPDestroy(&mgksp);
  VecDestroy(&mgctx[usermg->mglevels-1].user[bi].B);
  }
  
  return 0;
}

PetscErrorCode PoissonSolver_MG_old(UserMG *usermg)
{
  PetscInt l;
  PetscInt bi;

  MGCtx *mgctx = usermg->mgctx;

  KSP	*mgksp, subksp;
  PC	*mgpc, subpc;
  UserCtx	*user;

  PetscInt	m_c, m_f, M_c, M_f, flg;

  PetscMalloc(block_number*sizeof(KSP), &mgksp);
  PetscMalloc(block_number*sizeof(PC), &mgpc);

 
  PetscPrintf(PETSC_COMM_WORLD, "PoissonMG NullSpace\n");

  if (immersed) {
    for (l=usermg->mglevels-1; l>0; l--) {
      for (bi=0; bi<block_number; bi++) {
	mgctx[l].user[bi].multinullspace = PETSC_FALSE;
	MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]);
      }
    }
    /* At the corsest level, check whether the grid is separated into sevearal 
       blockes by the immersed body */
    l = 0;
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      PetscMalloc(user[bi].info.mx*user[bi].info.my*2*sizeof(PetscInt), &user[bi].KSKE);
      FullyBlocked(&user[bi]);

    }
  }

 
  l = usermg->mglevels-1;
  user = mgctx[l].user;
  for (bi=0; bi<block_number; bi++) {
    VecDuplicate(user[bi].P, &user[bi].B);
    PetscReal ibm_Flux, ibm_Area;
    
    flg=immersed-1;
    VolumeFlux(&user[bi], &ibm_Flux, &ibm_Area, flg);
    if (MHV || LV) {
    if ((MHV>1 || LV) && bi==0 ) 
      flg=1;
    else 
      flg=0;
    VolumeFlux_rev(&user[bi], &ibm_Flux, &ibm_Area, flg);
    }
    PoissonRHS(&(user[bi]), user[bi].B);
  
  }

  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      if (l==usermg->mglevels-1)
	PoissonLHSNew(&(user[bi]));
      else
	PoissonLHSNew(&(user[bi]));
    }
  }


  for (bi=0; bi<block_number; bi++) {
    /* Create ksp for Multigrid */
    KSPCreate(PETSC_COMM_WORLD, &mgksp[bi]);
    KSPAppendOptionsPrefix(mgksp[bi], "ps_");
    

    
    /* Use multigrid as preconditioner */
    KSPGetPC(mgksp[bi], &mgpc[bi]);
    PCSetType(mgpc[bi], PCMG);
    /* Setup MG levels from usergm->mglevels */
    PCMGSetLevels(mgpc[bi], usermg->mglevels, PETSC_NULL);
    /* V cycle */
    PCMGSetCycleType(mgpc[bi], PC_MG_CYCLE_V);
    
    PCMGSetType(mgpc[bi], PC_MG_MULTIPLICATIVE);


    /* Create Restriction and Interpolate schemes
       This is needed for all levels other than the coarsest one */
    for (l=usermg->mglevels-1; l>0; l--) {
      user = mgctx[l].user;
      m_c = (usermg->mgctx[l-1].user[bi].info.xm *
	     usermg->mgctx[l-1].user[bi].info.ym *
	     usermg->mgctx[l-1].user[bi].info.zm);

      m_f = (usermg->mgctx[l].user[bi].info.xm *
	     usermg->mgctx[l].user[bi].info.ym *
	     usermg->mgctx[l].user[bi].info.zm);

      M_c = (usermg->mgctx[l-1].user[bi].info.mx *
	     usermg->mgctx[l-1].user[bi].info.my *
	     usermg->mgctx[l-1].user[bi].info.mz);

      M_f = (usermg->mgctx[l].user[bi].info.mx *
	     usermg->mgctx[l].user[bi].info.my *
	     usermg->mgctx[l].user[bi].info.mz);

      MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, (void*)&mgctx[l-1].user[bi], &user[bi].MR);
      MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, (void*)&user[bi], &user[bi].MP);

      PCMGSetRestriction(mgpc[bi], l, user[bi].MR);
      PCMGSetInterpolation(mgpc[bi], l, user[bi].MP);

      /* Use subroutine MyRestriction and MyInterpolation for
	 Mat * Vec operation */
      MatShellSetOperation(user[bi].MR, MATOP_MULT, (void(*)(void))MyRestriction);
      MatShellSetOperation(user[bi].MP, MATOP_MULT, (void(*)(void))MyInterpolation);

      MatShellSetOperation(user[bi].MR, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);
      MatShellSetOperation(user[bi].MP, MATOP_MULT_ADD,(void(*)(void))mymatmultadd);

    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "PoissonMG Setup\n");

  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      
      if (l) { /* Not the coarset grid level */
	PCMGGetSmoother(mgpc[bi], l, &subksp);
	/* Set the left hand side for every KSP at each grid level */
	KSPSetOperators(subksp, user[bi].A, user[bi].A);
      
	KSPSetFromOptions(subksp);
	KSPSetUp(subksp);

	KSPGetPC(subksp, &subpc);

	KSP *subsubksp;
	PC subsubpc;
	PetscInt abi, nlocal;
	PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp);
	for (abi = 0; abi<nlocal; abi++) {
	  KSPGetPC(subsubksp[abi], &subsubpc);
	  PCFactorSetShiftType(subsubpc,MAT_SHIFT_NONZERO);
	}	

	PCFactorSetShiftType(subpc,MAT_SHIFT_NONZERO);
      }
      else {  /* Coarsest grid */

	/* The default solver for the coarset grid is
	   KSPPreonly and PCLU.
	   One can choose other solvers, such as a KSP solver by changing the
	   following lines. */
	PCMGGetCoarseSolve(mgpc[bi], &subksp);

	KSPGetPC(subksp, &subpc);
	KSPSetOperators(subksp, user[bi].A, user[bi].A);

	PCSetType(subpc, PCBJACOBI);
	KSPSetUp(subksp);

	KSP *subsubksp;
	PC subsubpc;
	PetscInt abi, nlocal;
	PCBJacobiGetSubKSP(subpc, &nlocal, PETSC_NULL, &subsubksp);
	for (abi = 0; abi<nlocal; abi++) {
	  KSPGetPC(subsubksp[abi], &subsubpc);
	  PCFactorSetShiftType(subsubpc,MAT_SHIFT_NONZERO);
	}

	KSPSetTolerances(subksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 40);
	

      }

      /* The Poisson equation has Neumann boundary conditions, thus
	 need to use NullSpace to solve the equation.
	 We use the subroutine PoissonNullSpaceFunction to
	 (1) reduce a constant value from the solution
	 (2) Set the solution values at boundaries and blanking nodes to
	 zero */
      MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp);
      MatNullSpaceSetFunction(user[bi].nullsp, PoissonNullSpaceFunction, &user[bi]);
    
      // KSPSetNullSpace(subksp, user[bi].nullsp);
      MatSetNullSpace( user[bi].A, user[bi].nullsp); //Mohsen Nov 2015
      //  PCMGSetResidual(mgpc[bi], l, PCMGDefaultResidual, user[bi].A);
      PCMGSetResidual(mgpc[bi], l, PCMGResidualDefault, user[bi].A);

      KSPSetUp(subksp);

      if (l<usermg->mglevels-1) {
	//	MatGetVecs(user[bi].A, &user[bi].R, PETSC_NULL);
	MatCreateVecs(user[bi].A, &user[bi].R,PETSC_NULL);
	PCMGSetRhs(mgpc[bi], l, user[bi].R);
      }

    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "PC for PoissonMG Setup\n");

  l = usermg->mglevels-1;
  user = mgctx[l].user;
  
  for (bi=0; bi<block_number; bi++) {
    KSPSetOperators(mgksp[bi], user[bi].A, user[bi].A);
    // KSPSetNullSpace(mgksp[bi], user[bi].nullsp);
    MatSetNullSpace( user[bi].A, user[bi].nullsp); //Mohsen Nov 2015
    KSPSetFromOptions(mgksp[bi]);
    KSPSetUp(mgksp[bi]);
    PetscPrintf(PETSC_COMM_WORLD, "KSP Setup & sovle for block %d\n",bi);

   
    KSPSolve(mgksp[bi], user[bi].B, user[bi].Phi);
  
  }

  for (bi=0; bi<block_number; bi++) {
    DMGlobalToLocalBegin(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
    DMGlobalToLocalEnd(user[bi].da, user[bi].Phi, INSERT_VALUES, user[bi].lPhi);
  }

 

  /* Release the allocated spaces */
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      MatNullSpaceDestroy(&user[bi].nullsp);

      MatDestroy(&user[bi].A);
      user[bi].assignedA = PETSC_FALSE;
      if (l) { /* Grid level > 0 (coarsest) */
	MatDestroy(&user[bi].MR);
	MatDestroy(&user[bi].MP);
      }
      else {
	PetscFree(user[bi].KSKE);
      }
    }
  }
  for (bi=0; bi<block_number; bi++) {
   
    
    for (l=0; l<usermg->mglevels-1; l++) {
      VecDestroy(&mgctx[l].user[bi].R);
    }
    
    KSPDestroy(&mgksp[bi]);
    VecDestroy(&mgctx[usermg->mglevels-1].user[bi].B);
  }
  PetscFree(mgksp);
  PetscFree(mgpc);

  return 0;
}


PetscErrorCode Projection(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;
  PetscReal	***p;

  Cmpnts	***ucont;
  PetscReal	***nvert;

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

  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lIZet, &izet);

  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lJZet, &jzet);

  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lPhi, &p);
  
  DMDAVecGetArray(fda, user->Ucont, &ucont);


  PetscReal dpdc, dpde, dpdz;
  PetscInt i_end,j_end,k_end;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	// Mohsen Aug 2012
	
	if (user->bctype[0]==7) i_end=mx-1;
	else i_end=mx-2;
	
	if (i<i_end) {
	  
	  dpdc = p[k][j][i+1] - p[k][j][i];
	  
	  dpde = 0.;
	  dpdz = 0.;
	 
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

	  if (!(nvert[k][j][i] + nvert[k][j][i+1])) {
	    ucont[k][j][i].x -= 
	      (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		       icsi[k][j][i].y * icsi[k][j][i].y +
		       icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	       dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		       ieta[k][j][i].y * icsi[k][j][i].y +
		       ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] + 
	       dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		       izet[k][j][i].y * icsi[k][j][i].y +
		       izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i])
	      * user->dt * user->st / COEF_TIME_ACCURACY;
	 
	  }
	}

	if (user->bctype[2]==7) j_end=my-1;
	else j_end=my-2;
	
	if (j<j_end) {
	 
	  dpdc = 0.;
	  dpdz = 0.;

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

	  if (!(nvert[k][j][i] + nvert[k][j+1][i])) {
	    ucont[k][j][i].y -=
	      (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		       jcsi[k][j][i].y * jeta[k][j][i].y +
		       jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	       dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		       jeta[k][j][i].y * jeta[k][j][i].y +
		       jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	       dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		       jzet[k][j][i].y * jeta[k][j][i].y +
		       jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i])
	      * user->dt * user->st / COEF_TIME_ACCURACY;
	  }
	}
	
	if (user->bctype[4]==7) k_end=mz-1;
	else k_end=mz-2;

	if (k < k_end) {
	 
	  dpdc = 0.;
	  dpde = 0.;
	  
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

	  dpdz = p[k+1][j][i] - p[k][j][i];
	  if (!(nvert[k][j][i] + nvert[k+1][j][i])) {
	    
	    ucont[k][j][i].z -=
	      (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		       kcsi[k][j][i].y * kzet[k][j][i].y +
		       kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	       dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		       keta[k][j][i].y * kzet[k][j][i].y +
		       keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	       dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		       kzet[k][j][i].y * kzet[k][j][i].y +
		       kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i])
	      * user->dt *user->st / COEF_TIME_ACCURACY;
	   
	  }
	}
      }
    }
  }

  // Mohsen Sep 2012
  // Update velocity at boundary for periodic BC's//
  // i-direction

  if (user->bctype[0]==7 && xs==0){
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	i=xs;
	
	dpdc = p[k][j][i+1] - p[k][j][i];
	
	dpde = 0.;
	dpdz = 0.;
		
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
	
	
	
	if (!(nvert[k][j][i] + nvert[k][j][i+1])) {
	  ucont[k][j][i].x -=
	    (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		     icsi[k][j][i].y * icsi[k][j][i].y +
		     icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	     dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		     ieta[k][j][i].y * icsi[k][j][i].y +
		     ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	     dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		     izet[k][j][i].y * icsi[k][j][i].y +
		     izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i])
	    * user->dt * user->st / COEF_TIME_ACCURACY;
	  
	}
      }
    }
  }

 
/*     // j-direction */

  if (user->bctype[2]==7 && ys==0){
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	j=ys;
	
	dpdc = 0.;
	dpdz = 0.;
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
	else if ((k == mz-2 || k==1 ) && user->bctype[0]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
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
		
	if (!(nvert[k][j][i] + nvert[k][j+1][i])) {
	  ucont[k][j][i].y -=
	    (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		     jcsi[k][j][i].y * jeta[k][j][i].y +
		     jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	     dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		     jeta[k][j][i].y * jeta[k][j][i].y +
		     jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	     dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		     jzet[k][j][i].y * jeta[k][j][i].y +
		     jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i])
	    * user->dt * user->st / COEF_TIME_ACCURACY;
	}
      }
    }
  }
  //k+direction
  if (user->bctype[4]==7 && zs==0){
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	k=zs;
	dpdc = 0.;
	dpde = 0.;

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
	
	dpdz = p[k+1][j][i] - p[k][j][i];
	
	if (!(nvert[k][j][i] + nvert[k+1][j][i])) {
	  
	  ucont[k][j][i].z -=
	   (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
	  		       kcsi[k][j][i].y * kzet[k][j][i].y +
	  		       kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	  	       dpde * (keta[k][j][i].x * kzet[k][j][i].x +
	  		       keta[k][j][i].y * kzet[k][j][i].y +
	  		       keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	  	       dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
	  		       kzet[k][j][i].y * kzet[k][j][i].y +
	  		       kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i])
	  	      * user->dt *user->st / COEF_TIME_ACCURACY;
	  
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);

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

  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);
  DMDAVecRestoreArray(da, user->lPhi, &p);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
 
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);


 
  Contra2Cart(user);
//  FormBCS(user);
  GhostNodeVelocity(user);
  return(0);
}

PetscErrorCode UpdatePressure(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
 
  PetscInt	lxs, lys, lzs, lxe, lye, lze,gxs, gxe, gys, gye, gzs, gze;
  PetscInt	i, j, k;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  Cmpnts ***cent;
  PetscReal ***p, ***phi, ***nvert, ***lp,***lphi ;

  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(da, user->Phi, &phi);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->Cent, &cent);
 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
		p[k][j][i] += phi[k][j][i];
      }
    }
  }
  
  DMDAVecRestoreArray(da, user->Phi, &phi);
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->Cent, &cent);

  //Mohsen Aug 2012
  //Updating pressure at boundaries for Periodic BC's
  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP,&lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=gzs; k<gze; k++) {
	  for (j=gys; j<gye; j++) {
	    if( k>0 && k<user->KM && j>0 && j<user->JM){
	      lp[k][j][i]=lp[k][j][i-2];
	      lphi[k][j][i]=lphi[k][j][i-2];
	    }
	  }
	}
      }
      if (xe==mx){
	i=mx-1;
	for (k=gzs; k<gze; k++) {
	  for (j=gys; j<gye; j++) {
	    if( k>0 && k<user->KM && j>0 && j<user->JM){
	      lp[k][j][i]=lp[k][j][i+2];
	      lphi[k][j][i]=lphi[k][j][i+2];	
	    }
	  }
	}
      }
    }
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

    DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi);
    DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);
    
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=gzs; k<gze; k++) {
	  for (i=gxs; i<gxe; i++) {
	    if( k>0 && k<user->KM && i>0 && i<user->IM){
	      lp[k][j][i]=lp[k][j-2][i];
	      lphi[k][j][i]=lphi[k][j-2][i];
	    }
	  }
	}
      }
      if (ye==my){
	j=my-1;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    if( k>0 && k<user->KM && i>0 && i<user->IM){
	      lp[k][j][i]=lp[k][j+2][i];
	      lphi[k][j][i]=lphi[k][j+2][i];
	    }
	  }
	}
      }
    }
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

    DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi);
    DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi);

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if( i>0 && i<user->IM && j>0 && j<user->JM){
	      lp[k][j][i]=lp[k-2][j][i];
	      lphi[k][j][i]=lphi[k-2][j][i];
	    }
	  }
	}
      }
      if (ze==mz){
	k=mz-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if( i>0 && i<user->IM && j>0 && j<user->JM){
	      lp[k][j][i]=lp[k+2][j][i];
	      lphi[k][j][i]=lphi[k+2][j][i];
	    }
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);
   
    DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi);
    DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi);

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=gzs; k<gze; k++) {
	  for (j=gys; j<gye; j++) {
	    lp[k][j][i]=lp[k][j][i-2];
	    lphi[k][j][i]=lphi[k][j][i-2];
	  }
	}
      }
      if (xe==mx){
	i=mx-1;
	for (k=gzs; k<gze; k++) {
	  for (j=gys; j<gye; j++) {
	    lp[k][j][i]=lp[k][j][i+2];
	    lphi[k][j][i]=lphi[k][j][i+2];	
	  }
	}
      }
    }
        
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);


    DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

    DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi);
    DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi);

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);
      
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=gzs; k<gze; k++) {
	  for (i=gxs; i<gxe; i++) {
	    lp[k][j][i]=lp[k][j+2][i];
	    lphi[k][j][i]=lphi[k][j+2][i];
	  }
	}
      }
      if (ys==0){
	j=ys;
	for (k=gzs; k<gze; k++) {
	  for (i=gxs; i<gxe; i++) {
	    lp[k][j][i]=lp[k][j-2][i];
	    lphi[k][j][i]=lphi[k][j-2][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

    DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi);
    DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi);

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp); 
    DMDAVecGetArray(da, user->lPhi, &lphi); 
   
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=gys; j<gye; j++) {
	  for (i=gxs; i<gxe; i++) {
	    lp[k][j][i]=lp[k+2][j][i];
	    lphi[k][j][i]=lphi[k+2][j][i];
	  }
	}
      }
      if (zs==0){
	k=zs;
	for (j=gys; j<gye; j++) {
	  for (i=gxs; i<gxe; i++) {
	    lp[k][j][i]=lp[k-2][j][i];
	    lphi[k][j][i]=lphi[k-2][j][i];
	  }
	}
      }
    }
  
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);
    DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi);
    DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi);
  }


  return 0;
}

PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg)
{
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscReal epsilon=1.e-8;
  PetscReal ***nvert, ibmval=1.1;
  
  struct Components {
    PetscReal x;
    PetscReal y;
    PetscReal z;
  }***ucor, ***csi, ***eta, ***zet;
 

  PetscInt xend=mx-2 ,yend=my-2,zend=mz-2;

  if (user->bctype[0]==7) xend=mx-1;
  if (user->bctype[2]==7) yend=my-1;
  if (user->bctype[4]==7) zend=mz-1;

  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area, libm_Flux_abs=0., ibm_Flux_abs;
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] < ibmval && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux += ucor[k][j][i].x;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							  csi[k][j][i].y * csi[k][j][i].y +
							  csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
	    } else 
	      ucor[k][j][i].x=0.;
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux += ucor[k][j][i].y;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/sqrt(eta[k][j][i].x * eta[k][j][i].x +
							  eta[k][j][i].y * eta[k][j][i].y +
							  eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    } else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux += ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux -= ucor[k][j][i].x;
	    if (flg==3)
	    libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							csi[k][j][i].y * csi[k][j][i].y +
							csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	    }else 
	      ucor[k][j][i].x=0.;
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux -= ucor[k][j][i].y;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/ sqrt(eta[k][j][i].x * eta[k][j][i].x +
							   eta[k][j][i].y * eta[k][j][i].y +
							   eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    }else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux -= ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

      }
    }
  }
  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_Flux_abs, &ibm_Flux_abs,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_Flux_abs, &ibm_Flux_abs, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */

  PetscReal correction;

  if (*ibm_Area > 1.e-15) {
    if (flg>1) 
      correction = (*ibm_Flux + user->FluxIntpSum)/ ibm_Flux_abs;
    else if (flg)
      correction = (*ibm_Flux + user->FluxIntpSum) / *ibm_Area;
    else
      correction = *ibm_Flux / *ibm_Area;
  }
  else {
    correction = 0;
  }

  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux %le %le %le\n", *ibm_Flux, *ibm_Area, correction);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] <ibmval && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon){
	    if (flg==3) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x + 
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < 1.1) {
	  if (nvert[k][j][i+1] < 0.1 && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x +
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

      }
    }
  }
  


  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] < ibmval && i < xend) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < yend) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < zend) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < xend) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < yend) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < zend) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux22 %le %le\n", *ibm_Flux, *ibm_Area);



  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xe==mx){
      i=mx-2;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  if(j>0 && k>0 && j<user->JM && k<user->KM){
	    if ((nvert[k][j][i]>1.1 && nvert[k][j][i+1]<0.1) || (nvert[k][j][i]<0.1 && nvert[k][j][i+1]>1.1)) ucor[k][j][i].x=0.0;
	    
	  }
	}
      }
    }
  }

  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ye==my){
      j=my-2;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && k>0 && i<user->IM && k<user->KM){
	    if ((nvert[k][j][i]>1.1 && nvert[k][j+1][i]<0.1) || (nvert[k][j][i]<0.1 && nvert[k][j+1][i]>1.1)) ucor[k][j][i].y=0.0;
	  }
	}
      }
    }
  }

  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (ze==mz){
      k=mz-2;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && j>0 && i<user->IM && j<user->JM){
	    if ((nvert[k][j][i]>1.1 && nvert[k+1][j][i]<0.1) || (nvert[k][j][i]<0.1 && nvert[k+1][j][i]>1.1)) ucor[k][j][i].z=0.0;
	  }
	}
      }
    }
  }


  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  /* periodci boundary condition update corrected flux */

  DMDAVecGetArray(fda, user->lUcont, &ucor);
  
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  if(j>0 && k>0 && j<user->JM && k<user->KM){
	    ucor[k][j][i].x=ucor[k][j][i-2].x;  
	  }
	}
      }
    }
  }
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && k>0 && i<user->IM && k<user->KM){
	    ucor[k][j][i].y=ucor[k][j-2][i].y;
	  }
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && j>0 && i<user->IM && j<user->JM){
	    ucor[k][j][i].z=ucor[k-2][j][i].z;
	  }
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucor);

  DMLocalToGlobalBegin(user->fda, user->lUcont, INSERT_VALUES, user->Ucont);
  DMLocalToGlobalEnd(user->fda, user->lUcont, INSERT_VALUES, user->Ucont);


  return 0;
}

PetscErrorCode VolumeFlux_rev(UserCtx *user, PetscReal *ibm_Flux, 
			      PetscReal *ibm_Area, PetscInt flg)
{
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
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

  PetscReal ***nvert, ibmval=1.5;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area;
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.4 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.4 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > ibmval-0.4 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > ibmval-0.4 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux REV %le %10f\n", *ibm_Flux, *ibm_Area);

  PetscReal correction;

  if (*ibm_Area > 1.e-15) {
    if (flg)
      correction = (*ibm_Flux + user->FluxIntpSum) / *ibm_Area;
    else
      correction = *ibm_Flux / *ibm_Area;
  }
  else {
    correction = 0;
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.4 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.4 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	  }
	  if (nvert[k+1][j][i] > ibmval-0.4 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	  }
	}

	if (nvert[k][j][i] > ibmval-0.4 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	  }
	}

      }
    }
  }
  


  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.4 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.4 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > ibmval-0.4 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > ibmval-0.4 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux22 REV %le %le\n", *ibm_Flux, *ibm_Area);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  return 0;
}


PetscErrorCode MaxPosition(UserCtx *user, PetscInt pos)
{
  
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt i, j, k;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (Gidx(i,j,k,user) == pos) {
	  PetscPrintf(PETSC_COMM_SELF, "Position %i %i %i\n", i, j, k);
	}
      }
    }
  }

  return 0;
}
#define Epsilon_Eq 1.e-6
//#define PartFloat_Eq(a, b) (fabs(a-b)>Epislon_Eq) ? 0:1;
#define Float_Eq(a, b) (a==b) ? PETSC_TRUE : (((a-b)<Epsilon_Eq) && (a-b) > -Epsilon_Eq)

PetscErrorCode MyNFaceFine(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscReal ***nvert;
  Cmpnts ***nface;

  PetscInt      i, j, k;
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

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->lNFace, &nface);

  VecSet(user->lNFace, 0.);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (Float_Eq(nvert[k][j][i], 0)) {
	  if (Float_Eq(nvert[k][j][i+1], 1)) {
	    nface[k][j][i].x = 1.;
	  }
	  if (Float_Eq(nvert[k][j+1][i], 1)) {
	    nface[k][j][i].y = 1.;
	  }
	  if (Float_Eq(nvert[k+1][j][i], 1)) {
	    nface[k][j][i].z = 1.;
	  }
	}
	else {
	  nface[k][j][i].x = 1.;
	  nface[k][j][i].y = 1.;
	  nface[k][j][i].z = 1.;
	}
      }
    }
  }
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lNFace, &nface);

  DMLocalToLocalBegin(fda, user->lNFace, INSERT_VALUES, user->lNFace);
  DMLocalToLocalEnd(fda, user->lNFace, INSERT_VALUES, user->lNFace);
  return 0;
}

PetscErrorCode MyNFaceRestrict(UserCtx *user)
{
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  
  PetscInt i,j,k;
  PetscInt ih, jh, kh, ia, ja, ka;
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

  Cmpnts ***nface_h, ***nface;
  PetscReal ***nvert;

  UserCtx *user_h = user->user_f;
  DM	fda_h = user_h->fda, fda = user->fda, da = user->da;
  DMDAVecGetArray(fda_h, user_h->lNFace, &nface_h);
  DMDAVecGetArray(fda,   user->lNFace, &nface);

  VecSet(user->lNFace, 0.);

  if ((user->isc)) ia = 0;
  else ia = 1;

  if ((user->jsc)) ja = 0;
  else ja = 1;

  if ((user->ksc)) ka = 0;
  else ka = 1;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);
	if (Float_Eq(nface_h[kh   ][jh   ][ih].x, 1) &&
	    Float_Eq(nface_h[kh   ][jh-ja][ih].x, 1) &&
	    Float_Eq(nface_h[kh-ka][jh   ][ih].x, 1) &&
	    Float_Eq(nface_h[kh-ka][jh-ja][ih].x, 1)) {
	  nface[k][j][i].x = 1.;
	}

	if (Float_Eq(nface_h[kh   ][jh][ih   ].y, 1) &&
	    Float_Eq(nface_h[kh   ][jh][ih-ia].y, 1) &&
	    Float_Eq(nface_h[kh-ka][jh][ih   ].y, 1) &&
	    Float_Eq(nface_h[kh-ka][jh][ih-ia].y, 1)) {
	  nface[k][j][i].y = 1.;
	}

	if (Float_Eq(nface_h[kh][jh   ][ih   ].z, 1) &&
	    Float_Eq(nface_h[kh][jh-ja][ih   ].z, 1) &&
	    Float_Eq(nface_h[kh][jh   ][ih-ia].z, 1) &&
	    Float_Eq(nface_h[kh][jh-ja][ih-ia].z, 1)) {
	  nface[k][j][i].z = 1.;
	}
      }
    }
  }

  DMDAVecGetArray(da, user->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (Float_Eq(nvert[k][j][i], 0)) {
	  if (Float_Eq(nvert[k][j][i+1], 1)) {
	    nface[k][j][i].x = 1.;
	  }
	  if (Float_Eq(nvert[k][j+1][i], 1)) {
	    nface[k][j][i].y = 1.;
	  }
	  if (Float_Eq(nvert[k+1][j][i], 1)) {
	    nface[k][j][i].z = 1.;
	  }
	}
	else {
	  nface[k][j][i].x = 1.;
	  nface[k][j][i].y = 1.;
	  nface[k][j][i].z = 1.;
	}
      }
    }
  }


  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  DMDAVecRestoreArray(fda, user->lNFace, &nface);
  DMDAVecRestoreArray(fda_h, user_h->lNFace, &nface_h);

  DMLocalToLocalBegin(fda, user->lNFace, INSERT_VALUES, user->lNFace);
  DMLocalToLocalEnd(fda, user->lNFace, INSERT_VALUES, user->lNFace);

  //  VecSet(user->lNFace, 0.);
  return 0;
}

PetscErrorCode MyNFaceInit(UserMG *usermg)
{
  PetscInt l;

  PetscInt bi;

  MGCtx *mgctx = usermg->mgctx;

  UserCtx *user;
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    for (bi=0; bi<block_number; bi++) {
      VecDuplicate(user[bi].lCsi, &user[bi].lNFace);
      if (l == usermg->mglevels-1) {
	MyNFaceFine(&user[bi]);
      }
      else {
	MyNFaceRestrict(&user[bi]);
      }
    }
  }
  return 0;
}

PetscErrorCode MyNFaceFinalize(UserMG *usermg)
{
  PetscInt l;

  MGCtx *mgctx = usermg->mgctx;

  UserCtx *user;
  for (l=usermg->mglevels-1; l>=0; l--) {
    user = mgctx[l].user;
    VecDestroy(&user->lNFace);
  }
  return 0;
}

PetscErrorCode PoissonLHSNew2(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;

  PetscInt lxs, lxe, lys, lye, lzs, lze;
  PetscInt gxs, gxe, gys, gye, gzs, gze;

  Vec		G11, G12, G13, G21, G22, G23, G31, G32, G33;
  PetscReal	***g11, ***g12, ***g13, ***g21, ***g22, ***g23;
  PetscReal	***g31, ***g32, ***g33;

  PetscReal	***nvert, ***nvert_o;
  PetscScalar	vol[19];
  PetscInt	idx[19], row;

  PetscInt	i, j, k, N;
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (!user->assignedA) {
    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(user->P, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 19, PETSC_NULL, 19, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  MatZeroEntries(user->A);

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

  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lNvert_o, &nvert_o);

  Cmpnts ***nface;
  DMDAVecGetArray(fda, user->lNFace, &nface);

  VecDuplicate(user->lAj, &G11);
  VecDuplicate(user->lAj, &G12);
  VecDuplicate(user->lAj, &G13);
  VecDuplicate(user->lAj, &G21);
  VecDuplicate(user->lAj, &G22);
  VecDuplicate(user->lAj, &G23);
  VecDuplicate(user->lAj, &G31);
  VecDuplicate(user->lAj, &G32);
  VecDuplicate(user->lAj, &G33);

  DMDAVecGetArray(da, G11, &g11);
  DMDAVecGetArray(da, G12, &g12);
  DMDAVecGetArray(da, G13, &g13);
  DMDAVecGetArray(da, G21, &g21);
  DMDAVecGetArray(da, G22, &g22);
  DMDAVecGetArray(da, G23, &g23);
  DMDAVecGetArray(da, G31, &g31);
  DMDAVecGetArray(da, G32, &g32);
  DMDAVecGetArray(da, G33, &g33);

  for (k=gzs; k<gze; k++) {
    for (j=gys; j<gye; j++) {
      for (i=gxs; i<gxe; i++) {
	g11[k][j][i] = (icsi[k][j][i].x * icsi[k][j][i].x +
			icsi[k][j][i].y * icsi[k][j][i].y +
			icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g12[k][j][i] = (ieta[k][j][i].x * icsi[k][j][i].x +
			ieta[k][j][i].y * icsi[k][j][i].y +
			ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g13[k][j][i] = (izet[k][j][i].x * icsi[k][j][i].x +
			izet[k][j][i].y * icsi[k][j][i].y +
			izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
					  			 
	g21[k][j][i] = (jcsi[k][j][i].x * jeta[k][j][i].x +
			jcsi[k][j][i].y * jeta[k][j][i].y +
			jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g22[k][j][i] = (jeta[k][j][i].x * jeta[k][j][i].x +
			jeta[k][j][i].y * jeta[k][j][i].y +
			jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g23[k][j][i] = (jzet[k][j][i].x * jeta[k][j][i].x +
			jzet[k][j][i].y * jeta[k][j][i].y +
			jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
					  			 
	g31[k][j][i] = (kcsi[k][j][i].x * kzet[k][j][i].x +
			kcsi[k][j][i].y * kzet[k][j][i].y +
			kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g32[k][j][i] = (keta[k][j][i].x * kzet[k][j][i].x +
			keta[k][j][i].y * kzet[k][j][i].y +
			keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g33[k][j][i] = (kzet[k][j][i].x * kzet[k][j][i].x +
			kzet[k][j][i].y * kzet[k][j][i].y +
			kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];	
      }
    }
  }

  PetscInt m;
  PetscReal threshold;
  if (user->thislevel == 2) {
    threshold = 0.1;
  }
  else {
    threshold = 0.1;
  }
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  vol[CP] = 1.;	idx[CP] = Gidx(i, j, k, user);
	  MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	}
	else {
	  if (nvert[k][j][i] > 0.1 ||
	      nface[k][j][i].x + nface[k][j][i].y + nface[k][j][i].z +
	      nface[k][j][i-1].x +
	      nface[k][j-1][i].y + nface[k-1][j][i].z > 5.1) { // i, j, k is not a fluid point
	    vol[CP] = 1.; idx[CP] = Gidx(i, j, k, user);
	    MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	  }
	  else { // i, j, k is a fluid point
	    for (m=0; m<19; m++) {
	      vol[m] = 0.;
	    }
	    /* Contribution from i+1 - i */
	    if (i != mx-2 && nface[k][j][i].x < 0.1) { // i+1, j, k is a fluid point
	      /* dpdc{i} = (p_{i+1} - p_{i}) * g11_{i} */
	      vol[CP] -= g11[k][j][i]; //i, j, k
	      vol[EP] += g11[k][j][i]; // i+1, j, k

	      /* dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1})
	       * 0.25 * g12[k][j][i] */
	      if (j == my-2 || nface[k][j][i].y + nface[k][j][i+1].y > 0.1) {
		if (nface[k][j-1][i].y + nface[k][j-1][i+1].y < 0.1 && j!=1) {
		  vol[CP] += g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] += g12[k][j][i] * 0.5; // i+1, j, k
		  vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
		  vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
		}
	      }
	      else if (j == 1 || nface[k][j-1][i].y + nface[k][j-1][i+1].y > 0.1) {
		if (nface[k][j][i].y + nface[k][j][i+1].y < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5;  //i, j+1, k
		  vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
		  vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
		}
	      }
	      else {
		vol[NP] += g12[k][j][i] * 0.25; // i, j+1, k
		vol[NE] += g12[k][j][i] * 0.25; // i+1, j+1, k
		vol[SP] -= g12[k][j][i] * 0.25; // i, j-1, k
		vol[SE] -= g12[k][j][i] * 0.25; // i+1, j-1, k
	      }
	      
	      /* dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if (k == mz-2 || nface[k][j][i].z + nface[k][j][i+1].z > 0.1) {
		if (nface[k-1][j][i].z + nface[k-1][j][i+1].z < 0.1 && k!=1) {
		  vol[CP] += g13[k][j][i] * 0.5; // i, j, k
		  vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
		  vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
		  vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
		}
	      }
	      else if (k == 1 || nface[k-1][j][i].z + nface[k-1][j][i+1].z > 0.1) {
		if (nface[k][j][i].z + nface[k][j][i+1].z < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
		  vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
		  vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
		  vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
		}
	      }
	      else {
		vol[TP] += g13[k][j][i] * 0.25; //i, j, k+1
		vol[TE] += g13[k][j][i] * 0.25; //i+1, j, k+1
		vol[BP] -= g13[k][j][i] * 0.25; //i, j, k-1
		vol[BE] -= g13[k][j][i] * 0.25; //i+1, j, k-1
	      }
	    }  // end of i+1 - i

	    /* Contribution from i - i-1 */
	    if (nface[k][j][i-1].x < threshold && i != 1) { // i-1, j, k is a fluid point
	      /* -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i} */
	      vol[CP] -= g11[k][j][i-1];  //i, j, k
	      vol[WP] += g11[k][j][i-1];  //i-1, j, k

	      /* -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1})
	       * 0.25 * g12[k][j][i-1] */
	      if (j == my-2 || nface[k][j][i].y + nface[k][j][i-1].y > 0.1) {
		if (nface[k][j-1][i].y + nface[k][j-1][i-1].y < 0.1 && j!=1) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
		  vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
		  vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
		}
	      }
	      else if (j == 1 || nface[k][j-1][i].y + nface[k][j-1][i-1].y > 0.1) {
		if (nface[k][j][i].y + nface[k][j][i-1].y < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; // i, j+1, k
		  vol[NW] -= g12[k][j][i-1] * 0.5; // i-1, j+1, k
		  vol[CP] += g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g12[k][j][i-1] * 0.5; // i-1, j, k
		}
	      }
	      else {
		vol[NP] -= g12[k][j][i-1] * 0.25; // i, j+1, k
		vol[NW] -= g12[k][j][i-1] * 0.25; //i-1, j+1, k
		vol[SP] += g12[k][j][i-1] * 0.25; // i, j-1, k
		vol[SW] += g12[k][j][i-1] * 0.25; // i-1, j-1, k
	      }

	      /* -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if (k == mz-2 || nface[k][j][i].z + nface[k][j][i-1].z > 0.1) {
		if (nface[k-1][j][i].z + nface[k-1][j][i-1].z < 0.1 && k!=1) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
		  vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
		  vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
		}
	      }
	      else if (k == 1 || nface[k-1][j][i].z + nface[k-1][j][i-1].z > 0.1) {
		if (nface[k][j][i].z + nface[k][j][i-1].z < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; // i, j, k+1
		  vol[TW] -= g13[k][j][i-1] * 0.5; //i-1, j, k+1
		  vol[CP] += g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g13[k][j][i-1] * 0.5; //i-1, j, k
		}
	      }
	      else {
		vol[TP] -= g13[k][j][i-1] * 0.25;  // i, j, k+1
		vol[TW] -= g13[k][j][i-1] * 0.25; // i-1, j, k+1
		vol[BP] += g13[k][j][i-1] * 0.25;  // i, j, k-1
		vol[BW] += g13[k][j][i-1] * 0.25; // i-1, j, k-1
	      }
	    } // end of i - i-1

	    /* Contribution from j+1 - j */
	    if (nface[k][j][i].y < threshold && j != my-2) {
	      /* dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) *
		 0.25 */
	      if (i == mx-2 || nface[k][j][i].x + nface[k][j+1][i].x > 0.1) {
		if (nface[k][j][i-1].x + nface[k][j+1][i-1].x < 0.1 && i!=1) {
		  vol[CP] += g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
		  vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
		  vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
		}
	      }
	      else if (i == 1 || nface[k][j][i-1].x + nface[k][j+1][i-1].x > 0.1) {
		if (nface[k][j][i].x + nface[k][j+1][i].x < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; // i+1, j, k
		  vol[NE] += g21[k][j][i] * 0.5; // i+1, j+1, k
		  vol[CP] -= g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] -= g21[k][j][i] * 0.5; // i, j+1, k
		}
	      }
	      else {
		vol[EP] += g21[k][j][i] * 0.25; //i+1, j, k
		vol[NE] += g21[k][j][i] * 0.25; //i+1, j+1, k
		vol[WP] -= g21[k][j][i] * 0.25; //i-1, j, k
		vol[NW] -= g21[k][j][i] * 0.25; //i-1, j+1, k
	      }
	     
	      /* dpde{j} = (p{j+1} - p{j}) * g22[k][j][i] */
	      vol[CP] -= g22[k][j][i];
	      vol[NP] += g22[k][j][i];

	      /* dpdz{j} = (p{j, k+1}+p{j+1,k+1} - p{j,k-1}-p{j+1,k-1}) *0.25*/
	      if (k == mz-2 || nface[k][j][i].z + nface[k][j+1][i].z > 0.1) {
		if (nface[k-1][j][i].z + nface[k-1][j+1][i].z < 0.1 && k!=1) {
		  vol[CP] += g23[k][j][i] * 0.5; //i,j,k
		  vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
		  vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
		  vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
		}
	      }
	      else if (k == 1 || nface[k-1][j][i].z + nface[k-1][j+1][i].z > 0.1) {
		if (nface[k][j][i].z + nface[k][j+1][i].z < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; //i, j, k+1
		  vol[TN] += g23[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g23[k][j][i] * 0.5;//i, j, k
		  vol[NP] -= g23[k][j][i] * 0.5;//i, j+1, k
		}
	      }
	      else {
		vol[TP] += g23[k][j][i] * 0.25; // i, j, k+1
		vol[TN] += g23[k][j][i] * 0.25; // i, j+1, k+1
		vol[BP] -= g23[k][j][i] * 0.25; // i, j, k-1
		vol[BN] -= g23[k][j][i] * 0.25; // i, j+1, k-1
	      }
	    } // End of j+1 - j

	    /* Contribution j - j-1 */
	    if (nface[k][j-1][i].y < threshold && j != 1) {
	      /* -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) *
		 0.25 * g21[k][j-1][i] */
	      if (i == mx-2 || nface[k][j][i].x + nface[k][j-1][i].x > 0.1) {
		if (nface[k][j][i-1].x + nface[k][j-1][i-1].x < 0.1 && i!=1) {
		  vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
		  vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
		  vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
		  vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
		}
	      }
	      else if (i == 1 || nface[k][j][i-1].x + nface[k][j-1][i-1].x > 0.1) {
		if (nface[k][j][i].x + nface[k][j-1][i].x < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5;//i+1, j, k
		  vol[SE] -= g21[k][j-1][i] * 0.5;//i+1, j-1, k
		  vol[CP] += g21[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g21[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[EP] -= g21[k][j-1][i] * 0.25;// i+1, j, k
		vol[SE] -= g21[k][j-1][i] * 0.25;// i+1, j-1, k
		vol[WP] += g21[k][j-1][i] * 0.25;// i-1, j, k
		vol[SW] += g21[k][j-1][i] * 0.25;// i-1, j-1, k
	      }
	      
	      /* -dpde{j-1} = -(p{j} - p{j-1}) * g22[k][j-1][i] */
	      vol[CP] -= g22[k][j-1][i];
	      vol[SP] += g22[k][j-1][i];

	      /* -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) *
		 0.25 * g23[k][j-1][i] */
	      if (k == mz-2 || nface[k][j][i].z + nface[k][j-1][i].z > 0.1) {
		if (nface[k-1][j][i].z + nface[k-1][j-1][i].z < 0.1 && k!=1) {
		  vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
		  vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
		  vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if (k == 1 || nface[k-1][j][i].z + nface[k-1][j-1][i].z > 0.1) {
		if (nface[k][j][i].x + nface[k][j-1][i].z < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5;//i, j, k+1
		  vol[TS] -= g23[k][j-1][i] * 0.5;//i, j-1, k+1
		  vol[CP] += g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g23[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[TP] -= g23[k][j-1][i] * 0.25;//i, j, k+1
		vol[TS] -= g23[k][j-1][i] * 0.25;//i, j-1, k+1
		vol[BP] += g23[k][j-1][i] * 0.25;//i, j, k-1
		vol[BS] += g23[k][j-1][i] * 0.25;//i, j-1, k-1
	      }
	    } // End of j - j-1

	    /* contribution from k+1 - k */
	    if (nface[k][j][i].z < threshold && k != mz-2) {
	      /* dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) *
		 0.25 * g31[k][j][i] */
	      if (i == mx-2 || nface[k][j][i].x + nface[k+1][j][i].x > 0.1) {
		if (nface[k][j][i-1].x + nface[k+1][j][i-1].x < 0.1 && i!=1) {
		  vol[CP] += g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
		  vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
		  vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
		}
	      }
	      else if (i == 1 || nface[k][j][i-1].x + nface[k+1][j][i-1].x > 0.1) {
		if (nface[k][j][i].x + nface[k+1][j][i].x < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5;//i+1, j, k
		  vol[TE] += g31[k][j][i] * 0.5;//i+1, j, k+1
		  vol[CP] -= g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g31[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[EP] += g31[k][j][i] * 0.25;//i+1, j, k
		vol[TE] += g31[k][j][i] * 0.25;//i+1, j, k+1
		vol[WP] -= g31[k][j][i] * 0.25;//i-1, j, k
		vol[TW] -= g31[k][j][i] * 0.25;//i-1, j, k+1
	      }

	      /* dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 
		 0.25 * g32[k][j][i] */
	      if (j == my-2 || nface[k][j][i].y + nface[k+1][j][i].y > 0.1) {
		if (nface[k][j-1][i].y + nface[k+1][j-1][i].y < 0.1 && j!=1) {
		  vol[CP] += g32[k][j][i] * 0.5;//i, j,k
		  vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
		  vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
		  vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
		}
	      }
	      else if (j == 1 || nface[k][j-1][i].y + nface[k+1][j-1][i].y > 0.1) {
		if (nface[k][j][i].y + nface[k+1][j][i].y < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5;//i, j+1, k
		  vol[TN] += g32[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g32[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g32[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[NP] += g32[k][j][i] * 0.25;//i, j+1, k
		vol[TN] += g32[k][j][i] * 0.25;//i, j+1, k+1
		vol[SP] -= g32[k][j][i] * 0.25;//i, j-1, k
		vol[TS] -= g32[k][j][i] * 0.25;//i, j-1, k+1
	      }

	      /* dpdz{k} = p{k+1} - p{k} */
	      vol[CP] -= g33[k][j][i]; //i, j, k
	      vol[TP] += g33[k][j][i]; //i, j, k+1
	    } // End of k+1 - k

	    /* Contribution from k - k-1 */
	    if (nface[k-1][j][i].z < threshold && k != 1) {
	      /* -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) *
		 0.25 * g31[k-1][j][i] */
	      if (i == mx-2 || nface[k][j][i].x + nface[k-1][j][i].x > 0.1) {
		if (nface[k][j][i-1].x + nface[k-1][j][i-1].x < 0.1 && i!=1) {
		  vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
		  vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
		  vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
		}
	      }
	      else if (i == 1 || nface[k][j][i-1].x + nface[k-1][j][i-1].x > 0.1) {
		if (nface[k][j][i].x + nface[k-1][j][i].x < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5;//i+1, j, k
		  vol[BE] -= g31[k-1][j][i] * 0.5;//i+1, j, k-1
		  vol[CP] += g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g31[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[EP] -= g31[k-1][j][i] * 0.25;//i+1, j, k
		vol[BE] -= g31[k-1][j][i] * 0.25;//i+1, j, k-1
		vol[WP] += g31[k-1][j][i] * 0.25;//i-1, j, k
		vol[BW] += g31[k-1][j][i] * 0.25;//i-1, j, k-1
	      }
	      
	      /* -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) * 
		 0.25 * g32[k-1][j][i] */
	      if (j == my-2 || nface[k][j][i].y + nface[k-1][j][i].y > 0.1) {
		if (nface[k][j-1][i].y + nface[k-1][j-1][i].y < 0.1 && j!=1) {
		  vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
		  vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
		  vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
		  vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if (j == 1 || nface[k][j-1][i].y + nface[k-1][j-1][i].y > 0.1) {
		if (nface[k][j][i].y + nface[k-1][j][i].y < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5;//i, j+1, k
		  vol[BN] -= g32[k-1][j][i] * 0.5;//i, j+1, k-1
		  vol[CP] += g32[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g32[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[NP] -= g32[k-1][j][i] * 0.25;//i, j+1, k
		vol[BN] -= g32[k-1][j][i] * 0.25;//i, j+1, k-1
		vol[SP] += g32[k-1][j][i] * 0.25;//i, j-1, k
		vol[BS] += g32[k-1][j][i] * 0.25;//i, j-1, k-1
	      }
	      
	      /* -dpdz{k-1} = -(p{k} - p{k-1}) * g33[k-1][j][i] */
	      vol[CP] -= g33[k-1][j][i]; // i, j, k
	      vol[BP] += g33[k-1][j][i]; //i, j, k-1
	    } // End of k - k-1
	    for (m=0; m<19; m++) {
	      vol[m] *= -aj[k][j][i];
	    }

	    idx[CP] = Gidx(i  , j  , k  , user);
	    idx[EP] = Gidx(i+1, j  , k  , user);
	    idx[WP] = Gidx(i-1, j  , k  , user);
	    idx[NP] = Gidx(i  , j+1, k  , user);
	    idx[SP] = Gidx(i  , j-1, k  , user);
	    idx[TP] = Gidx(i  , j  , k+1, user);
	    idx[BP] = Gidx(i  , j  , k-1, user);
	    idx[NE] = Gidx(i+1, j+1, k  , user);
	    idx[SE] = Gidx(i+1, j-1, k  , user);
	    idx[NW] = Gidx(i-1, j+1, k  , user);
	    idx[SW] = Gidx(i-1, j-1, k  , user);
	    idx[TN] = Gidx(i  , j+1, k+1, user);
	    idx[BN] = Gidx(i  , j+1, k-1, user);
	    idx[TS] = Gidx(i  , j-1, k+1, user);
	    idx[BS] = Gidx(i  , j-1, k-1, user);
	    idx[TE] = Gidx(i+1, j  , k+1, user);
	    idx[BE] = Gidx(i+1, j  , k-1, user);
	    idx[TW] = Gidx(i-1, j  , k+1, user);
	    idx[BW] = Gidx(i-1, j  , k-1, user);
	    MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);

/* 	    if (user->thislevel == 2) { */
/* 	      if (row==179161) { */
/* 		PetscInt tpa, tpb, tpc; */
/* 		for (tpa=k-1; tpa<k+2; tpa++) { */
/* 		  for (tpb=j-1; tpb<j+2; tpb++) { */
/* 		    for (tpc=i-1; tpc<i+2; tpc++) { */
/* 		      PetscPrintf(PETSC_COMM_SELF, "TT %i %i %i %e\n", tpa, tpb, tpc, nvert[tpa][tpb][tpc]); */
/* 		    } */
/* 		  } */
/* 		} */
/* 	      } */
/* 	    } */
	  } // End of fluid point
	} // End of interial points
      }
    }
  }
  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  DMDAVecRestoreArray(da, G11, &g11);
  DMDAVecRestoreArray(da, G12, &g12);
  DMDAVecRestoreArray(da, G13, &g13);
  DMDAVecRestoreArray(da, G21, &g21);
  DMDAVecRestoreArray(da, G22, &g22);
  DMDAVecRestoreArray(da, G23, &g23);
  DMDAVecRestoreArray(da, G31, &g31);
  DMDAVecRestoreArray(da, G32, &g32);
  DMDAVecRestoreArray(da, G33, &g33);

  DMDAVecRestoreArray(fda, user->lNFace, &nface);  
  VecDestroy(&G11);
  VecDestroy(&G12);
  VecDestroy(&G13);
  VecDestroy(&G21);
  VecDestroy(&G22);
  VecDestroy(&G23);
  VecDestroy(&G31);
  VecDestroy(&G32);
  VecDestroy(&G33);


  //  VecCopy(user->Phi, user->P);
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

  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);
  return 0;
}
