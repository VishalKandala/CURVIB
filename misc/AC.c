#include "variables.h"

extern PetscInt ti;
extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed, invicid;
extern PetscInt TwoD;
extern PetscInt moveframe,rotateframe;
extern PetscInt  les, clark, rans;

PetscErrorCode Cart2Contra(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  Vec           Ucat= user->lUcat, Ucont=user->Ucont;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***ucont, ***ucat;

  PetscReal	***nvert;

  PetscReal     innerblank=7.;


  PetscInt	i, j, k;

  Cmpnts	u_s;


  DAVecGetArray(fda, Csi, &csi);
  DAVecGetArray(fda, Eta, &eta);
  DAVecGetArray(fda, Zet, &zet);
  DAVecGetArray(da,  Aj,  &aj);

  DAVecGetArray(fda, Ucont, &ucont);
  DAVecGetArray(fda, Ucat,  &ucat);

  DAVecGetArray(da, user->lNvert, &nvert);

  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe-=2;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ( (nvert[k][j][i]+nvert[k][j][i+1] < 0.1 ) {
	  u_s.x= 0.5 * (ucat[k][j][i].x + ucat[k][j][i+1].x);
	  u_s.y= 0.5 * (ucat[k][j][i].y + ucat[k][j][i+1].y);
	  u_s.z= 0.5 * (ucat[k][j][i].z + ucat[k][j][i+1].z);
	  ucont[k][j][i].x = ( csi[k][j][i].x * u_s.x +
			       csi[k][j][i].y * u_s.y +
			       csi[k][j][i].z * u_s.z );
	}
      }
    }
  }

  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye-=2;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ( (nvert[k][j][i]+nvert[k][j+1][i] < 0.1 ) {
	  u_s.x= 0.5 * (ucat[k][j][i].x + ucat[k][j+1][i].x);
	  u_s.y= 0.5 * (ucat[k][j][i].y + ucat[k][j+1][i].y);
	  u_s.z= 0.5 * (ucat[k][j][i].z + ucat[k][j+1][i].z);
	  ucont[k][j][i].y = ( eta[k][j][i].x * u_s.x +
			       eta[k][j][i].y * u_s.y +
			       eta[k][j][i].z * u_s.z );
	}
      }
    }
  }


  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze-=2;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ( (nvert[k][j][i]+nvert[k+1][j][i] < 0.1 ) {
	  u_s.x= 0.5 * (ucat[k][j][i].x + ucat[k+1][j][i].x);
	  u_s.y= 0.5 * (ucat[k][j][i].y + ucat[k+1][j][i].y);
	  u_s.z= 0.5 * (ucat[k][j][i].z + ucat[k+1][j][i].z);
	  ucont[k][j][i].z = ( zet[k][j][i].x * u_s.x +
			       zet[k][j][i].y * u_s.y +
			       zet[k][j][i].z * u_s.z );
	}
      }
    }
  }

  DAVecRestoreArray(fda, Ucont, &ucont);
  DAVecRestoreArray(fda, Ucat,  &ucat);

  DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DAVecRestoreArray(fda, Csi, &csi);
  DAVecRestoreArray(fda, Eta, &eta);
  DAVecRestoreArray(fda, Zet, &zet);
  DAVecRestoreArray(da,  Aj,  &aj);
  DAGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DAVecRestoreArray(da, user->lNvert, &nvert);
  //  PetscPrintf(PETSC_COMM_WORLD, "Ve End\n");
  return(0);
}

PetscErrorCode FormRHS_BC(UserCtx *user, Vec Rhs)
{

  DA 		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt   lxs = xs, lxe = xe, lys = ys, lye = ye, lzs = zs, lze = ze;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  Cmpnts	***rhs;
  PetscReal	***nvert;
  PetscInt	i, j, k;

  DAVecGetArray(fda, Rhs, &rhs);

  if (TwoD==1)
    rhs[k][j][i].x =0.;
  else if (TwoD==2)
    rhs[k][j][i].y =0.;
  else if (TwoD==3)
    rhs[k][j][i].z =0.;
  
  DAVecGetArray(da, user->lNvert, &nvert);

  //rhs[k][j][i].x = 0.; rhs[k][j][i].y = 0.; rhs[k][j][i].z= 0.;
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
  
  /* i direction boundary conditions*/
/*   if (user->bctype[0]==7 && xs==0) { */
/*     i =  0; */
/*     for (k=lzs; k<lze; k++) { */
/*       for (j=lys; j<lye; j++) { */
/* 	rhs[k][j][i].x = rhs[k][j][mx-2].x; */
/* 	rhs[k][j][i].y = rhs[k][j][mx-2].y; */
/* 	rhs[k][j][i].z = rhs[k][j][mx-2].z; */
/*       } */
/*     } */
/*   } else  */
  if (xs==0) {
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
	if (user->bctype[1]!=7) {
	i = mx-2;
	rhs[k][j][i].x = 0;
	}
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

  DAVecRestoreArray(fda, Rhs, &rhs);

  return(0);
}

PetscErrorCode FormFunctionAC(UserCtx *user, Vec Rhs)
{
  Vec Conv,Visc,dP;
  VecDuplicate(user->lUcont, &Conv);
  VecDuplicate(user->lUcont, &Visc);
  VecDuplicate(user->lUcont, &dP);

  Cart2Contra(user);

  if (moveframe || rotateframe) {
    // a_c is grid speed
    Convection_MV(user, user->lUcont, user->lUcat, Conv);
  } else
  Convection(user, user->lUcont, user->lUcat, Conv);

  // Viscous term
  if (invicid)
    VecSet(Visc,0.);
  else
    Viscous(user, user->lUcont, user->lUcat, Visc);

  // pressure term
  PressureGradient(user, dP);

  // Right hand side terms from convective and viscous terms
  VecWAXPY(Rhs, -1., Conv, Visc);

  FormRHS_BC(user, Rhs);

  VecDestroy(Conv);
  VecDestroy(Visc);
  VecDestroy(dP);

  return(0);
}
