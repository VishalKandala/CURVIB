#include "variables.h"
#include "petsctime.h"

static char help[] = "Testing programming!";
PetscInt ti,tistart=0;
PetscReal	Flux_in = 4.104388e-04, angle = 0;
PetscInt tiout = 10;
PetscInt block_number;
PetscReal FluxInSum, FluxOutSum;
/* immeresed value>1 determines the type of correction
   1      constant velocity correction
   2      proportional (needs FluxOutSum>0)
   3      proportional to flux
   4      proportional to normal velocity (flux/area)
*/
PetscInt immersed = 0; 
PetscInt invicid = 0, channelz = 0;
PetscInt movefsi = 0, rotatefsi=0, sediment=0;
PetscBool rstart_flg;
PetscInt implicit = 0,implicit_type=0;
PetscInt imp_MAX_IT = 50; 
PetscInt radi=10;
PetscInt inletprofile=2, InitialGuessOne=0;
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
PetscInt  mg_MAX_IT=30, mg_idx=1, mg_preItr=1, mg_poItr=1;
PetscReal imp_atol=1e-7, imp_rtol=1e-4, imp_stol=1.e-8;
PetscInt TwoD = 0;
PetscInt STRONG_COUPLING=0;
PetscInt rstart_fsi=0;
PetscInt cop=0, regime=1; // 1 escape regime --- 2 cruise regime
PetscInt fish=0, fish_c=0, eel=0, fishcyl=0, rheology=0, aneurysm=0,turbine=0 , aneu_dom_bn=0;
PetscInt wing=0, hydro=0;
PetscReal St_exp=0.5,wavelength=0.95;
PetscInt MHV=0,LV=0;
PetscInt LVAD=0;
PetscReal max_angle = -54.*3.1415926/180.;// for MHV; min=0
PetscInt thin=0;
PetscInt dgf_z=0,dgf_y=1,dgf_x=0;
PetscInt dgf_az=0,dgf_ay=0,dgf_ax=1 ;
PetscInt NumberOfBodies=1;
PetscInt moveframe=0,rotateframe=0, blank=0;
PetscReal L_dim;
PetscInt les=0, rans=0;
PetscInt wallfunction=0;
PetscInt averaging=0;
int mixed=0;
int clark=0;
int dynamic_freq=1;
int poisson=0;
int periodic=0;
PetscReal max_cs=0.5;
int grid1d=0;
int i_periodic=0;
int j_periodic=0;
int k_periodic=0;
PetscInt blkpbc=10;
int pseudo_periodic=0;
int testfilter_ik=0;
int testfilter_1d=0;
int i_homo_filter=0;
int j_homo_filter=0;
int k_homo_filter=0;
double poisson_tol=5.e-9;	// relative tolerance
///////////////////////////////////////
PetscInt platelet=0;
PetscInt STEPS=0;// for reading from file
PetscInt LVV=0;
///////////////////////////////////////
/* // */
PetscInt orient = 0;
/* // */

//IBMNodes	*ibm_ptr;

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
  
  PetscViewer	viewer;
  
  char filen[90];
  
  
  sprintf(filen, "vfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  PetscInt N;

  VecGetSize(user->Ucont, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad((user->Ucont),viewer);
  
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);

  PetscOptionsClearValue("-vecload_block_size");

  sprintf(filen, "pfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->P),viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert_o),viewer);
  PetscViewerDestroy(&viewer);

  DMGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
  DMGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);

  if(averaging) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat",ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    VecSet(user->Ucat_sum, 0);
    VecSet(user->Ucat_cross_sum, 0);
    VecSet(user->Ucat_square_sum, 0);
    VecSet(user->P_sum, 0);
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and continues the computation ... ***\n\n", filen);
    }
    else {
      fclose(fp);
      PetscBarrier(PETSC_NULL);
      sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su1_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->Ucat_cross_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad((user->Ucat_square_sum),viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
    }
  }
  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat", ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting Cs to 0 and contiues the computation ... ***\n\n", filen);
      VecSet(Cs, 0);
    }
    else {
      fclose(fp);
      
      PetscBarrier(PETSC_NULL);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( Cs,viewer);
      PetscViewerDestroy(&viewer);
    }
    
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    
    VecDestroy(&Cs);
  }

  
  if(rans) {
    // K-Omega
    sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp!=NULL) {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->K_Omega,viewer);
      PetscViewerDestroy(&viewer);
    }
    else {
      K_Omega_IC(user);
    }
    
    VecCopy(user->K_Omega, user->K_Omega_o);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);    
  }

  return 0;
}

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user, PetscInt bi)
{
  PetscViewer	viewer;
  char filen[80];

  int rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  sprintf(filen, "vfield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Ucont, viewer);
  
  PetscViewerDestroy(&viewer);

  
  sprintf(filen, "ufield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Ucat, viewer);
  
  PetscViewerDestroy(&viewer);
  
  sprintf(filen, "pfield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->P, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);
  
  if(averaging && ti!=0) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su0_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "su1_%06d_%1d.dat",  ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_cross_sum, viewer);
    PetscViewerDestroy(&viewer);  
    sprintf(filen, "su1_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_square_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su2_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "sp_%06d_%1d.dat",ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->P_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "sp_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
  }
  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat",  ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(Cs, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(PETSC_NULL);
    VecDestroy(&Cs);
  }
  
  if(rans) {
    sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->K_Omega, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(PETSC_NULL);
  }
  
  return 0;
}

PetscErrorCode Divergence(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj, ***nvert;
  Cmpnts	***ucont;
  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda,user->lUcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);
 /*  for (k=zs; k<lze; k++) { */
/*     for (j=ys; j<lye; j++) { */
/*       for (i=xs; i<lxe; i++) { */
/* 	if (k==2 && j==9 && i==0) PetscPrintf(PETSC_COMM_SELF, "U1 @ i=%d j=%d k=2 is %le\n",i, j,ucont[k][j][i].x); */
/* 	if (k==2 && j==9 && i==1) PetscPrintf(PETSC_COMM_SELF, "U1 @ i=%d j=%d k=2 is %le\n",i, j,ucont[k][j][i].x); */
/* 	if (k==2 && j==8 && i==1) PetscPrintf(PETSC_COMM_SELF, "U2 @ i=%d j=%d k=2 is %le\n",i, j,ucont[k][j][i].y); */
/* 	if (k==2 && j==9 && i==1) PetscPrintf(PETSC_COMM_SELF, "U2 @ i=%d j=%d k=2 is %le\n",i, j,ucont[k][j][i].y); */
/* 	if (k==2 && j==9 && i==1) PetscPrintf(PETSC_COMM_SELF, "U3 @ i=%d j=%d k=2 is %le\n",i, j,ucont[k][j][i].z); */
/* 	if (k==1 && j==9 && i==1) PetscPrintf(PETSC_COMM_SELF, "U3 @ i=%d j=%d k=1 is %le\n",i, j,ucont[k][j][i].z); */
/*       } */
/*     } */
/*   } */

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x +
		       ucont[k][j][i].y - ucont[k][j-1][i].y +
		       ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	    nvert[k][j+1][i] + nvert[k][j-1][i] +
	    nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	div[k][j][i] = maxdiv;
	/* 	if (k==3 && j==4) PetscPrintf(PETSC_COMM_SELF, "div @ i=%d j=%d k=1 is %le\n",i, j,(ucont[k][j][i].x - ucont[k][j][i-1].x + */
/* 		       ucont[k][j][i].y - ucont[k][j-1][i].y + */
/* 		       ucont[k][j][i].z - ucont[k-1][j][i].z)); */

      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, Div, &div);
  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv %d %d %e\n", ti, i, maxdiv);
  PetscInt mi;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==i) {
	  PetscPrintf(PETSC_COMM_SELF, "MMa %d %d %d\n", mi,j, k);
	}
      }
    }
  }
  
  
  /*  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (Gidx(i,j,k,user)==137798) PetscPrintf(PETSC_COMM_WORLD, "Pos %d %d %d %e %e %e %e %e %e\n", i,j,k, ucont[k][j][i].x, ucont[k][j][i-1].x, ucont[k][j][i].y, ucont[k][j-1][i].y, ucont[k][j][i].z, ucont[k-1][j][i].z);
      }
    }
    }*/
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);
  return(0);
}

PetscErrorCode GridDivergence(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj, ***nvert;
  Cmpnts	***ucont;
  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal    norm;
  VecNorm(user->Vcont, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "Grid Flux norm  %le\n", norm);

  DMDAVecGetArray(fda,user->lVcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x +
		       ucont[k][j][i].y - ucont[k][j-1][i].y +
		       ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	    nvert[k][j+1][i] + nvert[k][j-1][i] +
	    nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	div[k][j][i] = maxdiv;
      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, Div, &div);
  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv Grid %d %d %e\n", ti, i, maxdiv);
  PetscInt mi;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==i) {
	  PetscPrintf(PETSC_COMM_SELF, "Maxdiv Grid location %d %d %d\n", mi,j, k);
	}
      }
    }
  }
    
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lVcont, &ucont);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);
  return(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv) {

  UserCtx    *user;
  PetscInt   i, bi, ibi;
  IBMNodes   *ibm;//, *ibm0, *ibm1;
  FSInfo     *fsi;
  PetscBool  DoSCLoop;
  PetscInt   itr_sc;
  Cstart     cstart;
  PetscInt   level;
  UserMG     usermg;
  PetscBool  flg;

  PetscInitialize(&argc, &argv, (char *)0, help);
// ------- INPUT PARAMETERS ----------------------
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-tio", &tiout, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imm", &immersed, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inv", &invicid, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &rstart_flg);
  PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_type", &implicit_type, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_MAX_IT", &imp_MAX_IT, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fsi", &movefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rfsi", &rotatefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-radi", &radi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inlet", &inletprofile, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-str", &STRONG_COUPLING, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rs_fsi", &rstart_fsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-cop", &cop, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-aneurysm", &aneurysm, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-aneu_dom_bn", &aneu_dom_bn, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rheology", &rheology, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-turbine", &turbine, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fishcyl", &fishcyl, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-eel", &eel, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-cstart", &fish_c, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-wing", &wing, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-sediment", &sediment, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mhv", &MHV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-hydro", &hydro, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-lv", &LV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-lvad", &LVAD, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-reg", &regime, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-twoD", &TwoD, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-thin", &thin, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_z", &dgf_z, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_y", &dgf_y, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_x", &dgf_x, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_az", &dgf_az, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_ay", &dgf_ay, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_ax", &dgf_ax, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-body", &NumberOfBodies, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mframe", &moveframe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rframe", &rotateframe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-blk", &blank, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-init1", &InitialGuessOne, PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-x_c", &(CMx_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-y_c", &(CMy_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-z_c", &(CMz_c), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-imp_atol", &(imp_atol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_rtol", &(imp_rtol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_stol", &(imp_stol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-poisson_tol", &poisson_tol, PETSC_NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.

  PetscOptionsGetInt(PETSC_NULL, "-les", &les, PETSC_NULL);				// Seokkoo Kang: if 1 Smagorinsky with Cs=0.1, if 2 Dynamic model
  PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);			// Seokkoo Kang

  PetscOptionsGetInt(PETSC_NULL, "-wallfunction", &wallfunction, PETSC_NULL);	// Seokkoo Kang: 1 or 2
  PetscOptionsGetInt(PETSC_NULL, "-mixed", &mixed, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
  PetscOptionsGetInt(PETSC_NULL, "-clark", &clark, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
  PetscOptionsGetInt(PETSC_NULL, "-dynamic_freq", &dynamic_freq, PETSC_NULL);		// Seokkoo Kang: LES dynamic compute frequency 
  PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL);	// Seokkoo Kang: if 1 do averaging; always begin with -rstart 0

  PetscOptionsGetInt(PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, PETSC_NULL);	
 
  PetscOptionsGetInt(PETSC_NULL, "-pbc_domain", &blkpbc, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-testfilter_ik", &testfilter_ik, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-testfilter_1d", &testfilter_1d, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-poisson", &poisson, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-i_homo_filter", &i_homo_filter, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-j_homo_filter", &j_homo_filter, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-k_homo_filter", &k_homo_filter, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-max_cs", &max_cs, PETSC_NULL);
  
  PetscOptionsGetReal(PETSC_NULL, "-St_exp", &(St_exp), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-wlngth", &(wavelength), PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-orient", &orient, PETSC_NULL);
    
  PetscReal compute_time,time_start,time_end;
  PetscPrintf(PETSC_COMM_WORLD, "tiout %d %le %le thin %d!\n",tiout, imp_atol,imp_rtol,thin);
  PetscTime(&time_start);

// ------------- SETUP PARAMETERS ---------------------------------

  if(LVAD){  // For LVAD; Set Characteristic Length to be 1.0, number of immersed boundaries to be 1 and also set the center of moment of the immersed boundary to be at 0,0,0.
    L_dim=1.;
    CMz_c=0.;CMy_c=0.;CMx_c=0.;
  } 
  else  L_dim=1.; 

  if (immersed) {  // If an immersed boundary exists, allocate memory for IBMNodes and FSInfo (data-types), then call the FsiInitialize() function.
    PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
    PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
    for (ibi=0;ibi<NumberOfBodies;ibi++)
         FsiInitialize(0, &fsi[ibi], ibi);     
  }

  MG_Initial(&usermg, ibm); // Call the Multi-grid initialize function that creates the grid, also initializes the fda, da and other data structures.
  
  if (immersed) {
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      
      PetscMalloc(NumberOfBodies*sizeof(IBMList), &(user[bi].ibmlist)); //For each block, allocate memory for IBMList and also initialize it.
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	InitIBMList(&(user[bi].ibmlist[ibi]));
      }
    }
  }
//------------------ READ IBM DATA -------------------------------------------------
    if(LVAD){
	ibi=0;  // set the number of blocks to be 1 i.e the first block has index 0;
        ibm_read_Icem(&ibm[ibi],ibi);  // This function reads the grid either from grid.dat or from the cartesian grid setup in the control file.
	PetscPrintf(PETSC_COMM_WORLD,"IBM Read LVAD!\n");
        ibm_surface_VTKOut(&ibm[ibi],ibi,0); // This function generates the VTK file that can be used to visualize the immersed boundary surface on ParaView.	
        FsiInitialize(0, &fsi[ibi], ibi);  // This function initializes FSI.
    }
    else {
      for (i=0;i<NumberOfBodies;i++) {
	
	PetscPrintf(PETSC_COMM_WORLD, "Ibm read General!\n");

	if (NumberOfBodies>10 && i>0)
	  ibm_read_ucd(&ibm[i], 0);
	
	else{
	  L_dim=1.;
	  CMz_c=0.0;CMy_c=0.0;CMx_c=0.0;
	  ibm_read_Icem(&ibm[i],i);	 
	  ibm_surface_VTKOut(&ibm[i],i,0);
	}	// init for fsi
	FsiInitialize(0, &fsi[i], i);
	
	if (NumberOfBodies>10 && i>0) {
	  PetscInt Nibm_y=13, Nibm_z=14;
	  PetscReal dYibm=1.5, dZibm=-1.5;
	  ibm_placement(&ibm[i], &fsi[i], Nibm_y, Nibm_z, dYibm, dZibm, i);
	} else
	  fsi[i].z_c +=i*2.;
	PetscBarrier(PETSC_NULL);
      }
    }
  
  if (immersed) {
    ti = 0;
    if (rstart_flg) ti = tistart;
  }
  //--------------- RESTART setup: if Restarting -----------------------------------
  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
  if (rstart_flg) {
    ti = tistart; tistart++;
    
    for (bi=0; bi<block_number; bi++) {
      Ucont_P_Binary_Input(&(user[bi]));
      
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);

      DMGlobalToLocalBegin(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      DMGlobalToLocalEnd(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      
      Contra2Cart(&(user[bi]));
     

      if (rstart_fsi) {
	for (ibi=1;ibi<NumberOfBodies;ibi++) {

	  FSI_DATA_Input(&fsi[ibi],ibi);

	  if (movefsi) {
	    if (!moveframe)
	      Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi],ibi);
	    for (i=0;i<6;i++){
	      fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
	      fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
	    }
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      ibm[ibi].uold[i].x = fsi[ibi].S_real[1];
	      ibm[ibi].uold[i].y = fsi[ibi].S_real[3];
	      ibm[ibi].uold[i].z = fsi[ibi].S_real[5];
	    }
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      ibm[ibi].urm1[i].x = fsi[ibi].S_realm1[1];
	      ibm[ibi].urm1[i].y = fsi[ibi].S_realm1[3];
	      ibm[ibi].urm1[i].z = fsi[ibi].S_realm1[5];
	    }
	  } // for movefsi

	  if ((rotatefsi &&  !rotateframe)|| MHV) {
	    //rotate_ibm(&ibm[ibi],&fsi[ibi]);
	    // calc_ibm_normal(&ibm[ibi]);
	    /* // */
	    Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi],user[bi].dt,ibi);
	    /* // */
	    ibm_surface_VTKOut(&ibm[ibi],ibi,0);
	    // if read ti, then will start for ti+1
	    for (i=0;i<6;i++){
	      fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
	      fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
	    }

	    fsi[ibi].F_x_real=fsi[ibi].F_x;
	    fsi[ibi].F_y_real=fsi[ibi].F_y;
	    fsi[ibi].F_z_real=fsi[ibi].F_z;
 
	    fsi[ibi].M_x_rm3=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm3=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm3=fsi[ibi].M_z;
	    
	    fsi[ibi].M_x_rm2=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm2=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm2=fsi[ibi].M_z;

	    fsi[ibi].M_x_real=fsi[ibi].M_x;
	    fsi[ibi].M_y_real=fsi[ibi].M_y;
	    fsi[ibi].M_z_real=fsi[ibi].M_z;


	    PetscReal rx,ry,rz;
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      rx = ibm[ibi].x_bp[i]-fsi[ibi].x_c;
	      ry = ibm[ibi].y_bp[i]-fsi[ibi].y_c;
	      rz = ibm[ibi].z_bp[i]-fsi[ibi].z_c;

	      ibm[ibi].u[i].x =   0.0  ;
	      ibm[ibi].u[i].y =-( -fsi[ibi].S_ang_n[1]*rz );
	      ibm[ibi].u[i].z =   -fsi[ibi].S_ang_n[1]*ry  ;

	      ibm[ibi].uold[i].x =  0.0  ;
	      ibm[ibi].uold[i].y =-( -fsi[ibi].S_ang_r[1]*rz );
	      ibm[ibi].uold[i].z =   -fsi[ibi].S_ang_r[1]*ry  ;
	      
	      ibm[ibi].urm1[i].x =   0.0  ;
	      ibm[ibi].urm1[i].y =-( -fsi[ibi].S_ang_rm1[1]*rz );
	      ibm[ibi].urm1[i].z =   -fsi[ibi].S_ang_rm1[1]*ry  ;
      
	    } // for nv
	  } // for rotatefsi and not rotateframe

	  if (rotatefsi &&  rotateframe) {
	    // Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt);
	    fsi[ibi].S_ang_n[5]=St_exp;
	    fsi[ibi].S_ang_r[5]=St_exp;

	    // if read ti, then will start for ti+1
	    for (i=0;i<6;i++){
	      fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
	      fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
	    }

	    fsi[ibi].F_x_real=fsi[ibi].F_x;
	    fsi[ibi].F_y_real=fsi[ibi].F_y;
	    fsi[ibi].F_z_real=fsi[ibi].F_z;
 
	    fsi[ibi].M_x_rm3=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm3=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm3=fsi[ibi].M_z;
	    
	    fsi[ibi].M_x_rm2=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm2=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm2=fsi[ibi].M_z;

	    fsi[ibi].M_x_real=fsi[ibi].M_x;
	    fsi[ibi].M_y_real=fsi[ibi].M_y;
	    fsi[ibi].M_z_real=fsi[ibi].M_z;

	    PetscReal rx,ry,rz;
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      rx = ibm[ibi].x_bp[i]-fsi[ibi].x_c;
	      ry = ibm[ibi].y_bp[i]-fsi[ibi].y_c;
	      rz = ibm[ibi].z_bp[i]-fsi[ibi].z_c;

	      ibm[ibi].u[i].x =-( ry*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[3]*rz );
	      ibm[ibi].u[i].y = ( rx*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[1]*rz );
	      ibm[ibi].u[i].z =-( rx*fsi[ibi].S_ang_n[3]-fsi[ibi].S_ang_n[1]*ry );

	      ibm[ibi].uold[i].x =-( ry*fsi[ibi].S_ang_r[5]-fsi[ibi].S_ang_r[3]*rz );
	      ibm[ibi].uold[i].y = ( rx*fsi[ibi/*  ibm->nf_x[i] =  -ibm->nf_x[i]; */
/*     ibm->nf_y[i] =  -ibm->nf_y[i]; */
/*     ibm->nf_z[i] =  -ibm->nf_z[i]; */].S_ang_r[5]-fsi[ibi].S_ang_r[1]*rz );
	      ibm[ibi].uold[i].z =-( rx*fsi[ibi].S_ang_r[3]-fsi[ibi].S_ang_r[1]*ry );

	      ibm[ibi].urm1[i].x =-( ry*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[3]*rz );
	      ibm[ibi].urm1[i].y = ( rx*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[1]*rz );
	      ibm[ibi].urm1[i].z =-( rx*fsi[ibi].S_ang_rm1[3]-fsi[ibi].S_ang_rm1[1]*ry );
      
	    } // for n_v
	  } // for rotateframe and rotatefsi

	}//ibi
      } // if rstart fsi

    }// bi

  } // if rstart
//---------------- SEARCH: Setting Solid, Boundary and Fluid nodes ---------------------------------------
// do the search once if elmt is not moving!
  PetscInt bis=0,bie=block_number;
  if (immersed) {
    for (level = usermg.mglevels-1; level>=usermg.mglevels-1; level--) {
      user = usermg.mgctx[level].user;
      for (bi=bis; bi<bie; bi++) {
 
	user[bi].ibm=ibm;

	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    if(LVAD){
	      ibm_search_advanced_rev(&(user[bi]), &ibm[ibi], ibi);  
  	      PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA REV LVAD  ibi %d bi %d\n", ibi,bi);        
            } 
            else {
	      ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
              PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA ibi %d bi %d\n", ibi, bi); 
	    }
	} //ibi

	PetscBarrier(PETSC_NULL);
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP\n");
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1);
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP End\n");

	} //ibi
      }// bi
    } //mglevels
  } // immersed

// --------------- SOLVER SETUP --------------------------
  ti = 0;
  if (rstart_flg) ti = tistart;

  if (ti==0) {
    if (InitialGuessOne) {  
      for (bi=0; bi<block_number; bi++) {
      SetInitialGuessToOne(&(user[bi]));   // Setup the initial condition for flow
      InflowFlux(&(user[bi]));  // Setup Inflow boundary condition
      OutflowFlux(&(user[bi])); // Setup Outflow boundary condition
      FormBCS(&(user[bi]));     // Setup all boundary conditions with info from bcs.dat
      if (wallfunction) {       // If wall function is used, setup Inflow and BCS again.
	InflowFlux(&(user[bi]));
	FormBCS(&(user[bi]));
      }
      }
      if (block_number>1) {
	Block_Interface_U(user);   // If more than one blocks(meshes) are present, setup interface conditions
      }

      for (bi=0; bi<block_number; bi++) {
      PetscReal normZ, normX, normY;
      VecStrideMax(user[bi].Ucat, 0, PETSC_NULL, &normX);
      VecStrideMax(user[bi].Ucat, 1, PETSC_NULL, &normY);
      VecStrideMax(user[bi].Ucat, 2, PETSC_NULL, &normZ);
      PetscPrintf(PETSC_COMM_WORLD, "Initial max cartesian velocities along: x -  %le; y- %le; z- %le\n",normX, normY, normZ);
      } // for bi
    } // for InitialGuessOne
  } // for ti=0

  
 for (bi=0; bi<block_number; bi++) {
    //VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o));
    VecCopy(user[bi].Ucont, user[bi].Ucont_o);      // Copy Ucont to Ucont_o for the finest level
    VecCopy(user[bi].Ucont, user[bi].Ucont_rm1);    // Copy Ucont to Ucont_rm1 for the finest level
    VecCopy(user[bi].Ucat, user[bi].Ucat_o);        // Copy Ucat to Ucat_o for the finest level
    VecCopy(user[bi].P, user[bi].P_o);              // Copy P to P_o for the finest level

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); // Global to Local for ucont,ucont_o and ucont_rm.
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
  } // for bi
 
  PetscInt tisteps = 5;
  PetscOptionsGetInt(PETSC_NULL, "-totalsteps", &tisteps, &flg); // Get input for real timesteps.
  if (tistart==0) tisteps ++; // if the simulation is starting from zero,run until the final timestep(hence increase timesteps by 1, so loop can be less than timesteps);
/* ================ SOLVE  ==================================================================             */
/*   physical time Step Loop */
  for (ti = tistart; ti<tistart + tisteps; ti++) {
    PetscPrintf(PETSC_COMM_WORLD, "Time level %d\n", ti);   // Print out time-level at the start of every physical time loop.
    /*-----Strong-Coupling (SC) Loop-----*/
    DoSCLoop= PETSC_TRUE ; itr_sc = 0;
    while (DoSCLoop) {
      itr_sc++;
      PetscPrintf(PETSC_COMM_WORLD, "SC LOOP itr # %d\n", itr_sc);
      if (immersed){
      	Struc_Solver(&usermg, ibm, fsi,&cstart, itr_sc,tistart, &DoSCLoop);        //Structral Solver!       
      } else  DoSCLoop = PETSC_FALSE;
      Flow_Solver(&usermg, ibm, fsi);     //Flow Solver!
    }// End of while SC loop
    /*------Save the old values (at ti) for later------*/
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {

      if (immersed) {
    	VecCopy(user[bi].Nvert, user[bi].Nvert_o);   // Saving nvert at ti; storing in a local array.
    	DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
    	DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      }
      
      VecCopy(user[bi].Ucat, user[bi].Ucat_o);           // Copy Ucat to Ucat_o for the finest level
      VecCopy(user[bi].Ucont_o, user[bi].Ucont_rm1);     // Copy Ucont to Ucont_rm1 for the finest level
      VecCopy(user[bi].Ucont, user[bi].Ucont_o);         // Copy Ucont to Ucont_o for the finest level
      VecCopy(user[bi].P, user[bi].P_o);                 // Copy P to P_o for the finest level

      // Save ucont_o and ucont_rm1 to a local array.
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);

      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);

    } // for bi
    /*------Calculate Max Error------- */
    
    Vec Error;
    PetscReal error=0.0;
    for (bi=0; bi<block_number; bi++) {   
      VecDuplicate(user[bi].Ucat,&Error); 
      VecWAXPY(Error,-1.0,user[bi].Ucat,user[bi].Ucat_o);
      VecNorm(Error,NORM_INFINITY,&error);
      // error=error/(user[bi].IM*user[bi].JM*user[bi].KM);
      //  if (error<epsilon) goto nextp;
      PetscPrintf(PETSC_COMM_WORLD, "Ucat l_infinity error is  %le \n", error);
      VecDestroy(&Error);
    } // bi
    /*---save location of immersed boundary and at time ti--*/
   
    if (immersed) {

      for (ibi=0;ibi<NumberOfBodies;ibi++) {

      for (i=0; i<ibm[ibi].n_v; i++) {   // loop through every immersed boundary node
    	ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i]; // copy x,y,z co-ordinates of each node to old co-ordinate arrays.
    	ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
    	ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];

    	ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;   // copy x,y,z velocities(? what is ibm.u??) to old arrays.
    	ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
    	ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;

    	ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
    	ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
    	ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
      }
      // FSI Related data to be stored for next iteration (elaborated later)
      for (i=0;i<6;i++){
    	fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
    	fsi[ibi].S_real[i]=fsi[ibi].S_new[i];

    	fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
    	fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
      }
      

      fsi[ibi].F_x_real=fsi[ibi].F_x;
      fsi[ibi].F_y_real=fsi[ibi].F_y;
      fsi[ibi].F_z_real=fsi[ibi].F_z;
      
      fsi[ibi].M_x_rm3=fsi[ibi].M_x_rm2;
      fsi[ibi].M_y_rm3=fsi[ibi].M_y_rm2;
      fsi[ibi].M_z_rm3=fsi[ibi].M_z_rm2;

      fsi[ibi].M_x_rm2=fsi[ibi].M_x_real;
      fsi[ibi].M_y_rm2=fsi[ibi].M_y_real;
      fsi[ibi].M_z_rm2=fsi[ibi].M_z_real;

      fsi[ibi].M_x_real=fsi[ibi].M_x;
      fsi[ibi].M_y_real=fsi[ibi].M_y;
      fsi[ibi].M_z_real=fsi[ibi].M_z;

      } //ibi
    } // if immersed    
  } // ti (physical time) loop
/* ==================================================================================             */
  PetscTime(&time_end);
  compute_time=time_end-time_start;
  PetscPrintf(PETSC_COMM_WORLD, "start time: %le , end time:%le ,compute time: %le \n",time_start,time_end,compute_time);
  MG_Finalize(&usermg); // Finalize (Destroy) grid fda, da and other data structures.
  PetscFinalize();
  return(0);
  }
