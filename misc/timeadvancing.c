#include "variables.h"

extern PetscReal FluxInSum, FluxOutSum;
extern PetscInt immersed;

PetscErrorCode MyFieldRestriction(UserCtx *user)
{
  DA	da = user->da, fda = user->fda;

  DA	da_f = *user->da_f;

  UserCtx *user_f = user->user_f;
  DA	fda_f = user_f->fda;

  DALocalInfo	info;
  DAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***ucont, ***ucont_f;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DAVecGetArray(fda, user->Ucont, &ucont);
  DAVecGetArray(fda_f, user_f->lUcont, &ucont_f);

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);
	
	ucont[k][j][i].x = (ucont_f[kh   ][jh   ][ih  ].x +
			    ucont_f[kh-ka][jh   ][ih  ].x +
			    ucont_f[kh   ][jh-ja][ih  ].x +
			    ucont_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	ucont[k][j][i].y = (ucont_f[kh   ][jh  ][ih   ].y +
			    ucont_f[kh-ka][jh  ][ih   ].y +
			    ucont_f[kh   ][jh  ][ih-ia].y +
			    ucont_f[kh-ka][jh  ][ih-ia].y) / (2.-ka) / (2.-ia);
      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	ucont[k][j][i].z = (ucont_f[kh  ][jh   ][ih   ].z +
			    ucont_f[kh  ][jh   ][ih-ia].z +
			    ucont_f[kh  ][jh-ja][ih   ].z +
			    ucont_f[kh  ][jh-ja][ih-ia].z) / (2.-ja) / (2.-ia);
      }
    }
  }

  DAVecRestoreArray(fda, user->Ucont, &ucont);
  DAVecRestoreArray(fda_f, user_f->lUcont, &ucont_f);
  
  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  VecCopy(user->Ucont, user->Ucont_MG);

  PetscReal ***p, ***p_f;

  DAVecGetArray(da, user->P, &p);
  DAVecGetArray(da_f, user_f->lP, &p_f);

  Cmpnts ***ucat, ***ucat_f;
  DAVecGetArray(user_f->fda, user_f->lUcat, &ucat_f);
  DAVecGetArray(user->fda, user->Ucat, &ucat); 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	p[k][j][i] = (p_f[kh   ][jh   ][ih   ] +
		      p_f[kh   ][jh   ][ih-ia] +
		      p_f[kh   ][jh-ja][ih   ] +
		      p_f[kh   ][jh-ja][ih-ia] +
		      p_f[kh-ka][jh   ][ih   ] +
		      p_f[kh-ka][jh   ][ih-ia] +
		      p_f[kh-ka][jh-ja][ih   ] +
		      p_f[kh-ka][jh-ja][ih-ia]) * 0.125;

/* 	ucat[k][j][i].x = (ucat_f[kh   ][jh   ][ih   ].x + */
/* 			   ucat_f[kh   ][jh   ][ih-ia].x + */
/* 			   ucat_f[kh   ][jh-ja][ih   ].x + */
/* 			   ucat_f[kh   ][jh-ja][ih-ia].x + */
/* 			   ucat_f[kh-ka][jh   ][ih   ].x + */
/* 			   ucat_f[kh-ka][jh   ][ih-ia].x + */
/* 			   ucat_f[kh-ka][jh-ja][ih   ].x + */
/* 			   ucat_f[kh-ka][jh-ja][ih-ia].x) * 0.125; */

/* 	ucat[k][j][i].y = (ucat_f[kh   ][jh   ][ih   ].y + */
/* 			   ucat_f[kh   ][jh   ][ih-ia].y + */
/* 			   ucat_f[kh   ][jh-ja][ih   ].y + */
/* 			   ucat_f[kh   ][jh-ja][ih-ia].y + */
/* 			   ucat_f[kh-ka][jh   ][ih   ].y + */
/* 			   ucat_f[kh-ka][jh   ][ih-ia].y + */
/* 			   ucat_f[kh-ka][jh-ja][ih   ].y + */
/* 			   ucat_f[kh-ka][jh-ja][ih-ia].y) * 0.125; */

/* 	ucat[k][j][i].z = (ucat_f[kh   ][jh   ][ih   ].z + */
/* 			   ucat_f[kh   ][jh   ][ih-ia].z + */
/* 			   ucat_f[kh   ][jh-ja][ih   ].z + */
/* 			   ucat_f[kh   ][jh-ja][ih-ia].z + */
/* 			   ucat_f[kh-ka][jh   ][ih   ].z + */
/* 			   ucat_f[kh-ka][jh   ][ih-ia].z + */
/* 			   ucat_f[kh-ka][jh-ja][ih   ].z + */
/* 			   ucat_f[kh-ka][jh-ja][ih-ia].z) * 0.125; */
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user_f->fda, user_f->lUcat, &ucat_f);
  
  DAVecRestoreArray(da, user->P, &p);
  DAVecRestoreArray(da_f, user_f->lP, &p_f);
  
  DAGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DAGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  DAGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DAGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat);
		      
  return 0;
}


PetscErrorCode ComputeRHS(UserCtx *user, PetscInt istage)
{
  Cmpnts ***rhs;
  Vec dUcont;
  VecDuplicate(user->Ucont, &dUcont);

  if (!user->Rhs) VecDuplicate(user->Ucont, &user->Rhs);
  FormFunction1(user->Ucont, user->Rhs, user);
  PetscReal norm;

  VecWAXPY(dUcont, -1., user->Ucont_o, user->Ucont);
  VecScale(dUcont, COEF_TIME_ACCURACY);
/*   VecAXPY(dUcont, -0.5, user->DUold); */
  VecAXPY(user->Rhs, -1./user->dt, dUcont);

  if (user->thislevel != user->mglevels-1) {
    if (!istage) {
      VecAXPY(user->Forcing, -1., user->Rhs);
    }
    VecAXPY(user->Rhs, 1., user->Forcing);
  }

  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;

  DAVecGetArray(fda, user->Rhs, &rhs);
  /* i direction boundary conditions*/
  if (xs ==0) {
    i = 0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  if (xe == mx) {
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	i = mx-2;
	rhs[k][j][i].x = 0;
	i = mx-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  /* j direction boundary conditions */
  if (ys == 0) {
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	j=0;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  if (ye == my) {
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	j=my-2;
	rhs[k][j][i].y = 0;
	j=my-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }
  /* k direction boundary conditions */
  if (zs == 0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

  if (ze == mz) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	k=mz-2;
	rhs[k][j][i].z = 0;
	k=mz-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }
  DAVecRestoreArray(fda, user->Rhs, &rhs);
  if (user->thislevel == user->mglevels-1) {
    VecNorm(user->Rhs, NORM_2, &norm);
    PetscPrintf(PETSC_COMM_WORLD, "TS RK Norm %e\n", norm);

    if (user->thislevel < user->mglevels-1 || PETSC_TRUE) {
      VecMax(user->Rhs, &k, &norm);
      PetscPrintf(PETSC_COMM_WORLD, "TS Forcing Max %i %e\n", k, norm);

      VecMin(user->Rhs, &k, &norm);
      PetscPrintf(PETSC_COMM_WORLD, "TS Forcing Min %i %e\n", k, norm);
    }
  }
  VecDestroy(dUcont);
  return 0;
}

PetscErrorCode RungeKutta_MG(UserCtx *user, PetscInt ps)
{
  PetscReal alfa[4];

  PetscInt i, j, k;
  DA	da = user->da, fda = user->fda;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	istage;


/*   alfa[0] = 1/3.; alfa[1] = 4./15.; alfa[2] = 5./9.; alfa[3] = 1.; */
  alfa[0] = 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

/*   if (user->thislevel == user->mglevels-1) { */
/*     InflowFlux(user); */
/*     OutflowFlux(user); */

/*     FormBCS(user); */
  InflowFlux(user);
/*   } */


  Vec DT;
  VecDuplicate(user->P, &user->Dt);
  VecSet(user->Dt, 0.);
  Spectral(user);
  VecDuplicate(user->Dt, &DT);

  PetscReal ***dt;
  Cmpnts	***rhs;

  Vec pUcont;
  VecDuplicate(user->Ucont, &pUcont);
  VecCopy(user->Ucont, pUcont);
  for (istage=0; istage<4; istage++) {

    Contra2Cart(user);
    GhostNodeVelocity(user);

    PetscInt ttta;
    ttta=istage+ps;
    ComputeRHS(user, ttta);

/*     Vec dUcont; */
/*     VecDuplicate(user->Ucont, &dUcont); */
/*     VecWAXPY(dUcont, -1., user->Ucont_o, user->Ucont); */
/*     VecScale(dUcont, COEF_TIME_ACCURACY); */
/* /\*     VecAXPY(dUcont, -0.5, user->DUold); *\/ */
/*     VecAXPY(user->Rhs, -1./user->dt, dUcont); */
/*     VecDestroy(dUcont); */

    VecCopy(user->Dt, DT);
    VecScale(DT, alfa[istage]);
    DAVecGetArray(fda, user->Rhs, &rhs);
    DAVecGetArray(da, DT, &dt);
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (i==0 || i==mx-1 || j==0 || j==my-1 ||
	      k==0 || k==mz-1) {
	    rhs[k][j][i].x = 0;
	    rhs[k][j][i].y = 0;
	    rhs[k][j][i].z = 0;
	  }
	  else {
	    if (i==mx-2) {
	      rhs[k][j][i].x = 0;
	    }
	    if (j== my-2) {
	      rhs[k][j][i].y = 0;
	    }
	    if (k==mz-2) {
	      rhs[k][j][i].z = 0;
	    }
	    rhs[k][j][i].x = dt[k][j][i] * rhs[k][j][i].x;
	    rhs[k][j][i].y = dt[k][j][i] * rhs[k][j][i].y;
	    rhs[k][j][i].z = dt[k][j][i] * rhs[k][j][i].z;
	  }
	}
      }
    }
    DAVecRestoreArray(fda, user->Rhs, &rhs);
    DAVecRestoreArray(da, DT, &dt);
    VecWAXPY(user->Ucont, 1., user->Rhs, pUcont);

/*     VecWAXPY(user->Ucont, alfa[istage]*user->dt*0.05*user->st, user->Rhs, pUcont); */
    DAGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DAGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

/*     if (user->thislevel==user->mglevels-1) { */
/*       InflowFlux(user); */
/*       OutflowFlux(user); */
/*       FormBCS(user); */
/*     } */
  }

  VecDestroy(pUcont);
  VecDestroy(DT);
  VecDestroy(user->Dt);
  if (user->thislevel == user->mglevels-1) {
    InflowFlux(user);
    OutflowFlux(user);
    
    FormBCS(user);
    if (immersed) {
      ibm_interpolation_advanced(user, user->ibm);
    }

  }

  return 0;
}

PetscErrorCode MyRKRHSRestriction(UserCtx *user)
{
  DA	da = user->da, fda = user->fda;

  UserCtx *user_f = user->user_f;
  DA	fda_f = user_f->fda;

  DALocalInfo	info;
  DAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***f, ***rhs_f;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if (*(user->isc)) ia = 0;
  else ia = 1;

  if (*(user->jsc)) ja = 0;
  else ja = 1;

  if (*(user->ksc)) ka = 0;
  else ka = 1;

  VecSet(user->Forcing, 0.);
  DAVecGetArray(fda, user->Forcing, &f);
  Vec lRhs;
  VecDuplicate(user_f->lUcont, &lRhs);
  DAGlobalToLocalBegin(fda_f, user_f->Rhs, INSERT_VALUES, lRhs);
  DAGlobalToLocalEnd(fda_f, user_f->Rhs, INSERT_VALUES, lRhs);
  DAVecGetArray(fda_f, lRhs, &rhs_f);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	if (i==0 || i==mx-2) {
	  f[k][j][i].x = 0.;
	}
	else {
	  f[k][j][i].x = ((rhs_f[kh   ][jh   ][ih   ].x +
			   rhs_f[kh-ka][jh   ][ih   ].x +
			   rhs_f[kh   ][jh-ja][ih   ].x +
			   rhs_f[kh-ka][jh-ja][ih   ].x) * (1.+ka) * (1.+ja) +
			  (rhs_f[kh   ][jh   ][ih-ia].x +
			   rhs_f[kh-ka][jh   ][ih-ia].x +
			   rhs_f[kh   ][jh-ja][ih-ia].x +
			   rhs_f[kh-ka][jh-ja][ih-ia].x) * (1.+ka) * (1.+ja) / 2.+
			  (rhs_f[kh   ][jh   ][ih+ia].x +
			   rhs_f[kh-ka][jh   ][ih+ia].x +
			   rhs_f[kh   ][jh-ja][ih+ia].x +
			   rhs_f[kh-ka][jh-ja][ih+ia].x) * (1.+ka) * (1.+ja) / 2.)
	    * 0.125;
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	if (j==0 || j==my-2) {
	  f[k][j][i].y = 0.;
	}
	else {
	  f[k][j][i].y = ((rhs_f[kh   ][jh   ][ih   ].y +
			   rhs_f[kh-ka][jh   ][ih   ].y +
			   rhs_f[kh   ][jh   ][ih-ia].y +
			   rhs_f[kh-ka][jh   ][ih-ia].y) * (1.+ka) * (1.+ia) +
			  (rhs_f[kh   ][jh-ja][ih   ].y +
			   rhs_f[kh-ka][jh-ja][ih   ].y +
			   rhs_f[kh   ][jh-ja][ih-ia].y +
			   rhs_f[kh-ka][jh-ja][ih-ia].y) * (1.+ka) * (1.+ia) / 2 +
			  (rhs_f[kh   ][jh+ja][ih   ].y +
			   rhs_f[kh-ka][jh+ja][ih   ].y +
			   rhs_f[kh   ][jh+ja][ih-ia].y +
			   rhs_f[kh-ka][jh+ja][ih-ia].y) * (1.+ka) * (1.+ia) / 2)
	    * 0.125;
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	if (k==0 || k==mz-2) {
	  f[k][j][i].z = 0.;
	}
	else {
	  f[k][j][i].z = ((rhs_f[kh   ][jh   ][ih   ].z +
			   rhs_f[kh   ][jh   ][ih-ia].z +
			   rhs_f[kh   ][jh-ja][ih   ].z +
			   rhs_f[kh   ][jh-ja][ih-ia].z) * (1.+ja) * (1.+ia) +
			  (rhs_f[kh-ka][jh   ][ih   ].z +
			   rhs_f[kh-ka][jh   ][ih-ia].z +
			   rhs_f[kh-ka][jh-ja][ih   ].z +
			   rhs_f[kh-ka][jh-ja][ih-ia].z) * (1.+ja) * (1.+ia) / 2 +
			  (rhs_f[kh+ka][jh   ][ih   ].z +
			   rhs_f[kh+ka][jh   ][ih-ia].z +
			   rhs_f[kh+ka][jh-ja][ih   ].z +
			   rhs_f[kh+ka][jh-ja][ih-ia].z) * (1.+ja) * (1.+ia) / 2)
	    * 0.125;
	}
      }
    }
  }

  DAVecRestoreArray(fda, user->Forcing, &f);
/*   VecScale(user->Forcing, -1); */
  DAVecRestoreArray(fda_f, lRhs, &rhs_f);
  VecDestroy(lRhs);
  return 0;
}

#define GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user) \
  if (*(user->isc)) { \
    ic = i; \
    ia = 0; \
  } \
  else { \
    ic = (i+1) / 2; \
    ia = (i - 2 * (ic)) == 0 ? 1 : -1; \
    if (i==1 || i==mx-2) ia = 0; \
  }\
  if (*(user->jsc)) { \
    jc = j; \
    ja = 0; \
  } \
  else { \
    jc = (j+1) / 2; \
    ja = (j - 2 * (jc)) == 0 ? 1 : -1; \
    if (j==1 || j==my-2) ja = 0; \
  } \
  if (*(user->ksc)) { \
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


PetscErrorCode MyRKRHSInterpolation(UserCtx *user)
{
  DA	da = user->da, fda = user->fda;

  DA	da_c = *user->da_c;

  UserCtx *user_c = user->user_c;

  DA	fda_c = user_c->fda;

  DALocalInfo	info = user->info;
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

  Vec dU;
  Cmpnts	***du, ***f;

  PetscInt	i, j, k, ic, jc, kc, ia, ja, ka;
  PetscReal	***nvert_c;

  PetscInt ii, jj, kk;
  if (*(user->isc)) ii = 0;
  else ii = 1;

  if (*(user->jsc)) jj = 0;
  else jj = 1;

  if (*(user->ksc)) kk = 0;
  else kk = 1;

  VecDuplicate(user->Ucont, &dU);
  VecSet(dU, 0.);

  Vec ldU;
  DAGetLocalVector(fda, &ldU);
  VecSet(ldU, 0.);
  DAVecGetArray(fda, ldU, &du);


  VecWAXPY(user_c->Forcing, -1., user_c->Ucont_MG, user_c->Ucont);

  Vec lForcing;
  DACreateLocalVector(fda_c, &lForcing);
  DAGlobalToLocalBegin(fda_c, user_c->Forcing, INSERT_VALUES, lForcing);
  DAGlobalToLocalEnd(fda_c, user_c->Forcing, INSERT_VALUES, lForcing);

  DAVecGetArray(fda_c, lForcing, &f);

  DAVecGetArray(da_c, *(user->lNvert_c), &nvert_c);

/*   totally wrong! */
/*     1) du is flux not velocity! therefore surface or aj should be */
/*        come into calcutations. */

  /* i component */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	if (!(i%2)) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].x = (f[kc   ][jc   ][ic   ].x * 9 +
			   f[kc   ][jc+ja][ic   ].x * 3 +
			   f[kc+ka][jc   ][ic   ].x * 3 +
			   f[kc+ka][jc+ja][ic   ].x) / 16. / ((1.+kk)*(1.+jj));
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DALocalToLocalBegin(fda, ldU, INSERT_VALUES, ldU);
  DALocalToLocalEnd(fda, ldU, INSERT_VALUES, ldU);

  DAVecGetArray(fda, ldU, &du);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {      
      for (i=lxs; i<lxe; i++) {

	if ((i%2)) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].x = (du[k][j][i-1].x + du[k][j][i+1].x) * 0.5;
	}

      }
    }
  }

  /* j component */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      if (!(j%2)) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].y = (f[kc   ][jc   ][ic   ].y * 9 +
			   f[kc   ][jc   ][ic+ia].y * 3 +
			   f[kc+ka][jc   ][ic   ].y * 3 +
			   f[kc+ka][jc   ][ic+ia].y) / 16. / ((1.+kk)*(1.+ii));
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DALocalToLocalBegin(fda, ldU, INSERT_VALUES, ldU);
  DALocalToLocalEnd(fda, ldU, INSERT_VALUES, ldU);

  DAVecGetArray(fda, ldU, &du);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      if (j%2) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].y = (du[k][j-1][i].y + du[k][j+1][i].y) * 0.5;

	}
      }
    }
  }

  /* k component */
  for (k=lzs; k<lze; k++) {
    if (!(k%2)) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].z = (f[kc   ][jc   ][ic   ].z * 9 +
			   f[kc   ][jc   ][ic+ia].z * 3 +
			   f[kc   ][jc+ja][ic   ].z * 3 +
			   f[kc   ][jc+ja][ic+ia].z) / 16. / ((1.+jj)*(1.+ii));
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DALocalToLocalBegin(fda, ldU, INSERT_VALUES, ldU);
  DALocalToLocalEnd(fda, ldU, INSERT_VALUES, ldU);

  DAVecGetArray(fda, ldU, &du);

  for (k=lzs; k<lze; k++) {
    if ((k%2)) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  du[k][j][i].z = (du[k-1][j][i].z + du[k+1][j][i].z) * 0.5;
	}
      }
    }
  }

  DAVecRestoreArray(fda, ldU, &du);
  DAVecRestoreArray(fda_c, lForcing, &f);

  DALocalToGlobal(fda, ldU, INSERT_VALUES, dU);


  VecAXPY(user->Ucont, 1., dU);

  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DARestoreLocalVector(fda, &ldU);
  VecDestroy(dU);
  VecDestroy(lForcing);

  DAVecRestoreArray(da_c, *(user->lNvert_c), &nvert_c);
  
  return 0;
}

#define MG_MAXIT 12
PetscErrorCode TSMG(UserCtx *user)
{
  PetscInt l;

  PetscInt mglevels = user->mglevels;

  PetscInt mgit, ps, *PseudoSteps;

  PetscMalloc(mglevels*sizeof(PetscInt), &PseudoSteps);
  for (l=0; l<mglevels; l++) {
    PseudoSteps[l] = (l||mglevels==1) ? 1:10;
  }

/*   if (immersed) { */
/*     for (l=usermg->mglevels-1; l>0; l--) { */
/*       for (bi=0; bi<block_number; bi++) { */
/* 	MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]); */
/*       } */
/*     } */
/*   } */
/*   Cmpnts ***zet; */
/*   Cmpnts ***ucont; */
/*   PetscReal temp; */
/*   DAVecGetArray(user->fda, user->lZet, &zet); */
/*   temp = zet[1][1][1].z; */
/*   DAVecGetArray(user->fda, user->Ucont, &ucont); */
/*   PetscInt i,j,k; */
/*   DALocalInfo	info = user->info; */
/*   PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*   PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*   PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*   PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
  
/*   for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	ucont[k][j][i].z = temp; */
/*       } */
/*     } */
/*   } */
/*   DAVecRestoreArray(user->fda, user->Ucont, &ucont); */

/*   DAVecRestoreArray(user->fda, user->lZet, &zet); */


  UserCtx *user_mg;
  user_mg = user;
  for (l=mglevels-1; l>=0; l--) {
    if (user->thislevel == user->mglevels-1) {
      InflowFlux(user);
      OutflowFlux(user);
      
      FormBCS(user);
      if (immersed) {
	ibm_interpolation_advanced(user, user->ibm);
      }
    }

    VecDuplicate(user_mg->Ucont, &user_mg->Ucont_o);
    VecCopy(user_mg->Ucont, user_mg->Ucont_o);
    if (l) { // Not on the coarest grid level
      MyFieldRestriction(user_mg->user_c);
      user_mg = user_mg->user_c;
    }
  }

  user_mg = user;
  PetscTruth MG_Converged = PETSC_FALSE;
  PetscReal norm, norm_old, temp;

  norm_old=0.;
  for (mgit = 0; mgit < MG_MAXIT; mgit++) {
    for (l=mglevels-1; l>=0; l--) {
      for (ps=0; ps<PseudoSteps[l]; ps++) {
	RungeKutta_MG(user_mg, ps);
      }

      if (l==mglevels-1) {
	VecWAXPY(user_mg->Rhs, -1., user_mg->Ucont_o, user_mg->Ucont);
	VecNorm(user_mg->Rhs, NORM_2, &norm);
	temp = norm;
	norm = norm-norm_old;
	norm_old = temp;
	PetscPrintf(PETSC_COMM_WORLD, "TS MG Norm %e\n", norm);
/* 	  PetscInt k; */
/* 	  VecMin(user_mg->Rhs, &k, &norm) */
/* ;	  PetscPrintf(PETSC_COMM_WORLD, "TS Norm %i %e\n", k, norm); */
	if (norm < 1.e-20 && norm > -1.e-20) {
	  PetscPrintf(PETSC_COMM_WORLD, "TSConverged %i\n", mgit);
	  MG_Converged = PETSC_TRUE;
	  break;
	}
      }

      if (l>0) {
	//	user_mg->istage = 0;
	ComputeRHS(user_mg, 1);
	MyFieldRestriction(user_mg->user_c);
	MyRKRHSRestriction(user_mg->user_c);
	user_mg = user_mg->user_c;
      }
      if (MG_Converged) break;
    }

    if (MG_Converged) break;
    if (mglevels>1) user_mg = user_mg->user_f;
    for (l=1; l<mglevels; l++) {
      MyRKRHSInterpolation(user_mg);
/*       for (ps=0; ps<PseudoSteps[l]; ps++) { */
/* 	if (l<mglevels-1) VecSet(user_mg->Forcing, 0.); */
/* 	RungeKutta_MG(user_mg); */
/*       } */
      if (l<mglevels-1) user_mg = user_mg->user_f;
    }
  }
  PetscFree(PseudoSteps);

  user_mg = user;
  for (l=mglevels-1; l>=0; l--) {
    VecDestroy(user_mg->Ucont_o);
    user_mg = user_mg->user_c;
  }

/*   user->mglevels=3; */
  return 0;
}
