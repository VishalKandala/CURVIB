#include "variables.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscsnes.h"
//#include "petscdmmg.h"
#include "petscdmda.h"             /*I      "petscda.h"     I*/
#include "petscmat.h"
#include "petsctime.h"
#include "petscdmcomposite.h"
//#include "private/dmimpl.h"

extern PetscInt block_number, NumberOfBodies,blank;
extern PetscInt immersed, rans, les;
extern PetscInt ti,tistart,tiout,visflg;
extern PetscInt imp_MAX_IT,implicit_type;
extern PetscReal imp_atol, imp_rtol, imp_stol;
extern PetscInt mg_idx, mg_preItr, mg_poItr, mg_MAX_IT;

//PetscInt lidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);

extern int Contra2Cart(UserCtx *);
extern int  NSNullSpaceFunction(MatNullSpace,Vec, void*);
extern int ModifyJacobian(Mat , UserCtx *);
extern int FormJacobian_FD(SNES,Vec,Mat,Mat,void *);
extern int FormJacobian_MF(SNES,Vec,Mat,Mat,void *);

extern int FormJacobian(SNES,Vec,Mat,Mat,void *);
extern int FormJacobian_Diagonal(SNES,Vec,Mat,Mat,void *);
extern int FormJacobian_Diagonal_MF(Mat,void *);
extern int FormJacobian_Full_Diagonal(SNES,Vec,Mat,Mat,void *);

extern int FormFunction_SNES_Packer(SNES, Vec, Vec, void *);
extern int SNESSetInitialGuess_Packer(Vec,UserMG *);
extern int SNESFormFinal_Packer(Vec ,UserMG *);
extern int FormJacobian_Diagonal_Packer_All(SNES,Vec,Mat,Mat,void *);
extern int FormJacobian_Diagonal_Packer(Mat ,void *);
extern int FormJacobian_Packer(Mat,void *);
extern int ModifyJacobian_packer(Mat,UserCtx *);



typedef struct {
  Vec diag;
} SampleShellPC;
SampleShellPC  *shell;    /* user-defined preconditioner context */
extern PetscErrorCode SampleShellPCCreate(SampleShellPC**);
extern PetscErrorCode SampleShellPCSetUp(PC,Mat,Vec);
extern PetscErrorCode SampleShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode SampleShellPCDestroy(PC);
PetscInt BSize[5],BlSize[5];

PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
PetscInt counter=0;
PetscReal cput=0.0;

PetscErrorCode MySNESMonitor(SNES snes,PetscInt n,PetscReal rnorm,void *dummy)
{
	PetscPrintf(PETSC_COMM_WORLD,"     (%D) SNES Residual norm %14.12e \n",n,rnorm);
	return 0;
}


PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, 
			  FSInfo *fsi)
{
  PetscInt      istage;
  PetscReal     alfa[4];
  PetscInt	bi, ibi;
  
  alfa[0] = 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;
  // Create local working vector for Ucont and Ucat components
  
  /*   DAGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
  /*   DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
  

  
  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    
    FormBCS(&(user[bi]));
    if (immersed) 
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
      }
    
    //    ibm_interpolation_advanced(&user[bi], ibm, fsi);
    //ibm_interpolation_advanced2(&user[bi], ibm);
    /*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    /*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
    /*
      VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
      VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
      VecCopy(user[bi].Ucont, user[bi].dUcont);
      VecCopy(user[bi].Ucont, user[bi].pUcont);
    */
  }
  
  
  //  if (ti>0) {
  
  PetscInt pseudot;
  // pseudo time iteration
  for(pseudot=0; pseudot<1;pseudot++) {
    for (bi=0; bi<block_number; bi++) {
      for (istage=0; istage<4; istage++) {
	
	FormFunction1(&(user[bi]), user[bi].Rhs);
	
	/* 	// Calculate du/dt */
	/* 	  VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont); */
	/* 	  VecScale(user[bi].dUcont, 1.5); */
	/* 	  VecAXPY(user[bi].dUcont, -0.5, user[bi].DUold); */
	/* 	  // Add -du/dt to right hand side */
	/* 	  VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont); */
	// Advanced in time using RK scheme
	VecWAXPY(user[bi].Ucont, alfa[istage] * user[bi].dt * user[bi].st, user[bi].Rhs, user[bi].Ucont_o);
	
	DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	
	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
	PetscPrintf(PETSC_COMM_WORLD, "FluxOutRK %le\n", user->FluxOutSum);
	
	FormBCS(&(user[bi]));
	/* 	FormBCS(user); */
	/* 	if (immersed) { */
	/* 	  ibm_interpolation(ibminfo, user, ibm); */
	/* 	} */
	/* 	  VecCopy(user[bi].Ucont, user[bi].pUcont); */
	
      }//istage
      if (immersed) 
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	}
      
      //ibm_interpolation_advanced(user, ibm, fsi);
      /* 	  ibm_interpolation(ibminfo, user, ibm); */
      
      
      //	VecCopy(user[bi].Ucont, user[bi].pUcont);
    }
    if (block_number>1) {
      Block_Interface_U(user);
    }
  }
  //  }
  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].Rhs);
    /*
      VecWAXPY(user[bi].DUold, -1., user[bi].Ucont_o, user[bi].Ucont);
      //    VecDestroy(user[bi].Ucont_o); 
      VecDestroy(user[bi].dUcont);
      VecDestroy(user[bi].pUcont);
    */
  }
  
  return 0;  
}

PetscErrorCode ImplicitSolverLHSnew03(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ***dtow, ucon;
  PetscReal	***nvert;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.001;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",dir,j,k );

 if (!user->assignedA) {
   //PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

 // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

  MatZeroEntries(user->A);

  if (dir==0) {
    DMDAVecGetArray(fda, user->lICsi, &csi);
    DMDAVecGetArray(fda, user->lIEta, &eta);
    DMDAVecGetArray(fda, user->lIZet, &zet);
    DMDAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecGetArray(fda, user->lJCsi, &csi);
    DMDAVecGetArray(fda, user->lJEta, &eta);
    DMDAVecGetArray(fda, user->lJZet, &zet);
    DMDAVecGetArray(da, user->lJAj, &aj);
  } else {
    DMDAVecGetArray(fda, user->lKCsi, &csi);
    DMDAVecGetArray(fda, user->lKEta, &eta);
    DMDAVecGetArray(fda, user->lKZet, &zet);
    DMDAVecGetArray(da, user->lKAj, &aj);
  }
  //DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->psuedot, &dtow);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 || nvert[k][j][i] + nvert[k][j][i+1] < 0.1) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];
	    
/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    // diagonal 
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*csi0*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*csi0*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*0.25*csi0*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].x)>1e-10)
		val[4] += 0.5*0.25*csi0*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*csi0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*csi0*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].x;			
	      }	  
	    }//csi0

	    if (fabs(csi1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*csi1*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*csi1*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*0.25*csi1*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].y)>1e-10)
		val[4] += 0.5*0.25*csi1*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*csi1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*csi1*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].y;			
	      }	  
	    }//csi1

	    if (fabs(csi2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(csi[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(csi[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(csi[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*0.25*csi2*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(csi[k][j+1][i].z)>1e-10)
		val[4] += 0.5*0.25*csi2*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(csi[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(csi[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*csi2*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].z;			
	      }	  
	    }//csi2

	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 || nvert[k][j][i] + nvert[k][j+1][i] < 0.1) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].x)>1e-10)
		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*0.25*eta0*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].x)>1e-10)
		val[6] += 0.5*0.25*eta0*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].x;			
	      }	  
	    }//eta0

	    if (fabs(eta1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].y)>1e-10)
		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*0.25*eta1*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].y)>1e-10)
		val[6] += 0.5*0.25*eta1*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)	
		  /eta[k+1][j][i].y;			
	      }	  
	    }//eta1

	    if (fabs(eta2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(eta[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(eta[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(eta[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*eta2*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(eta[k][j+1][i].z)>1e-10)
		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(eta[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*0.25*eta2*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(eta[k+1][j][i].z)>1e-10)
		val[6] += 0.5*0.25*eta2*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].z;			
	      }	  
	    }//eta2


	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k<mz-2 || nvert[k][j][i] + nvert[k+1][j][i] < 0.1) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2/Re;

		if (fabs(zet[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2/Re;
		if (fabs(zet[k][j][i+1].x)>1e-10)
		  val[2]+= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2/Re;
		if (fabs(zet[k][j-1][i].x)>1e-10)
		val[3] -= 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2/Re;
		if (fabs(zet[k][j+1][i].x)>1e-10)
		val[4] += 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2/Re;
		if (fabs(zet[k-1][j][i].x)>1e-10)
		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2/Re;
		if (fabs(zet[k+1][j][i].x)>1e-10)
		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].x;			
	      }	  
	    }//zet0

	    if (fabs(zet1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].y)>1e-10)
		  val[1]-= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].y)>1e-10)
		  val[2]+= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].y)>1e-10)
		val[3] -= 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].y)>1e-10)
		val[4] += 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].y)>1e-10)
		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].y)>1e-10)
		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].y;			
	      }	  
	    }//zet1

	    if (fabs(zet2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
		if (fabs(zet[k][j][i-1].z)>1e-10)
		  val[1]-= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
		if (fabs(zet[k][j][i+1].z)>1e-10)
		  val[2]+= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
		if (fabs(zet[k][j-1][i].z)>1e-10)
		val[3] -= 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
		if (fabs(zet[k][j+1][i].z)>1e-10)
		val[4] += 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
		if (fabs(zet[k-1][j][i].z)>1e-10)
		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
		if (fabs(zet[k+1][j][i].z)>1e-10)
		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;			
	      }	  
	    }//zet2


	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	
	
	}	
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DMDAVecRestoreArray(fda, user->lICsi, &csi);
    DMDAVecRestoreArray(fda, user->lIEta, &eta);
    DMDAVecRestoreArray(fda, user->lIZet, &zet);
    DMDAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecRestoreArray(fda, user->lJCsi, &csi);
    DMDAVecRestoreArray(fda, user->lJEta, &eta);
    DMDAVecRestoreArray(fda, user->lJZet, &zet);
    DMDAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DMDAVecRestoreArray(fda, user->lKCsi, &csi);
    DMDAVecRestoreArray(fda, user->lKEta, &eta);
    DMDAVecRestoreArray(fda, user->lKZet, &zet);
    DMDAVecRestoreArray(da, user->lKAj, &aj);
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);  
  DMDAVecRestoreArray(fda, user->psuedot, &dtow);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode ImplicitSolverLHSnew04(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ***dtow, ucon;
  PetscReal	***nvert, ***lnu_t;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren, nu_t=0.;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.0;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  PetscReal     zero=1e-10;
  

 if (!user->assignedA) {

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  MatZeroEntries(user->A);

  if (dir==0) {
    DMDAVecGetArray(fda, user->lICsi, &csi);
    DMDAVecGetArray(fda, user->lIEta, &eta);
    DMDAVecGetArray(fda, user->lIZet, &zet);
    DMDAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecGetArray(fda, user->lJCsi, &csi);
    DMDAVecGetArray(fda, user->lJEta, &eta);
    DMDAVecGetArray(fda, user->lJZet, &zet);
    DMDAVecGetArray(da, user->lJAj, &aj);
  } else {
    DMDAVecGetArray(fda, user->lKCsi, &csi);
    DMDAVecGetArray(fda, user->lKEta, &eta);
    DMDAVecGetArray(fda, user->lKZet, &zet);
    DMDAVecGetArray(da, user->lKAj, &aj);
  }
 
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->psuedot, &dtow);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  if (rans || les) 
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (rans || les) nu_t=lnu_t[k][j][i];

	if (i == 0  || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 && (nvert[k][j][i] + nvert[k][j][i+1] < 0.1)) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	   	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
		    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];
	    
	    // diagonal 
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
   	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i-1].x)>zero)
		  val[1]-= 0.5*aj1*csi0*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3  || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i+1].x)>zero)
		  val[2]+= 0.5*aj1*csi0*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j-1][i].x)>zero)
		val[3] -= 0.5*0.25*csi0*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j+1][i].x)>zero)
		val[4] += 0.5*0.25*csi0*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k-1][j][i].x)>zero)
		val[5] -= 0.5*0.25*csi0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k+1][j][i].x)>zero)
		val[6] += 0.5*0.25*csi0*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].x;			
	      }	  
	    }//csi0

	    if (fabs(csi1)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i-1].y)>zero)
		  val[1]-= 0.5*aj1*csi1*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i+1].y)>zero)
		  val[2]+= 0.5*aj1*csi1*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1|| (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j-1][i].y)>zero)
		val[3] -= 0.5*0.25*csi1*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j+1][i].y)>zero)
		val[4] += 0.5*0.25*csi1*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k-1][j][i].y)>zero)
		val[5] -= 0.5*0.25*csi1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k+1][j][i].y)>zero)
		val[6] += 0.5*0.25*csi1*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].y;			
	      }	  
	    }//csi1

	    if (fabs(csi2)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ( i==1 || (nvert[k][j][i]+nvert[k][j][i-1]>0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i-1].z)>zero)
		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i+1].z)>zero)
		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i]+nvert[k][j-1][i+1]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j-1][i].z)>zero)
		val[3] -= 0.5*0.25*csi2*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i]+nvert[k][j+1][i+1]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j+1][i].z)>zero)
		val[4] += 0.5*0.25*csi2*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j][i+1]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k-1][j][i].z)>zero)
		val[5] -= 0.5*0.25*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j][i+1]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k+1][j][i].z)>zero)
		val[6] += 0.5*0.25*csi2*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].z;			
	      }	  
	    }//csi2

	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 && (nvert[k][j][i] + nvert[k][j+1][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].y+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2.*(nu_t+1/Re)*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2.*(nu_t+1/Re)*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i-1].x)>zero)
		  val[1]-= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i+1].x)>zero)
		  val[2]+= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j-1][i].x)>zero)
		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j+1][i].x)>zero)
		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k-1][j][i].x)>zero)
		val[5] -= 0.5*0.25*eta0*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k+1][j][i].x)>zero)
		val[6] += 0.5*0.25*eta0*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].x;			
	      }	  
	    }//eta0

	    if (fabs(eta1)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i-1].y)>zero)
		  val[1]-= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i+1].y)>zero)
		  val[2]+= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j-1][i].y)>zero)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j+1][i].y)>zero)
		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k-1][j][i].y)>zero)
		val[5] -= 0.5*0.25*eta1*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k+1][j][i].y)>zero)
		val[6] += 0.5*0.25*eta1*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)	
		  /eta[k+1][j][i].y;			
	      }	  
	    }//eta1

	    if (fabs(eta2)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i-1].z)>zero)
		  val[1]-= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i+1].z)>zero)
		  val[2]+= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)		    
		    /eta[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j-1][i].z)>zero)
		val[3] -= 0.5*eta2*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j+1][i].z)>zero)
		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k-1][j][i].z)>zero)
		val[5] -= 0.5*0.25*eta2*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.; 
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k+1][j][i].z)>zero)
		val[6] += 0.5*0.25*eta2*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].z;			
	      }	  
	    }//eta2


	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k<mz-2 && (nvert[k][j][i] + nvert[k+1][j][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].z+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2.*(nu_t+1/Re)*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2.*(nu_t+1/Re)*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2*(nu_t+1/Re);

		if (fabs(zet[k][j][i-1].x)>zero)
		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i+1].x)>zero)
		  val[2]+= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].x;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j-1][i].x)>zero)
		val[3] -= 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j+1][i].x)>zero)
		val[4] += 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].x;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k-1][j][i].x)>zero)
		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k+1][j][i].x)>zero)
		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].x;			
	      }	  
	    }//zet0

	    if (fabs(zet1)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i-1].y)>zero)
		  val[1]-= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].y;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i+1].y)>zero)
		  val[2]+= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].y;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j-1][i].y)>zero)
		val[3] -= 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j+1][i].y)>zero)
		val[4] += 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].y;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k-1][j][i].y)>zero)
		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k+1][j][i].y)>zero)
		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].y;			
	      }	  
	    }//zet1

	    if (fabs(zet2)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i-1].z)>zero)
		  val[1]-= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].z;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i+1].z)>zero)
		  val[2]+= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].z;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j-1][i].z)>zero)
		val[3] -= 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j+1][i].z)>zero)
		val[4] += 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].z;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k-1][j][i].z)>zero)
		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k+1][j][i].z)>zero)
		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;			
	      }	  
	    }//zet2


	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	
	
	}	
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DMDAVecRestoreArray(fda, user->lICsi, &csi);
    DMDAVecRestoreArray(fda, user->lIEta, &eta);
    DMDAVecRestoreArray(fda, user->lIZet, &zet);
    DMDAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecRestoreArray(fda, user->lJCsi, &csi);
    DMDAVecRestoreArray(fda, user->lJEta, &eta);
    DMDAVecRestoreArray(fda, user->lJZet, &zet);
    DMDAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DMDAVecRestoreArray(fda, user->lKCsi, &csi);
    DMDAVecRestoreArray(fda, user->lKEta, &eta);
    DMDAVecRestoreArray(fda, user->lKZet, &zet);
    DMDAVecRestoreArray(da, user->lKAj, &aj);
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);  
  DMDAVecRestoreArray(fda, user->psuedot, &dtow);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  if (rans || les) 
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  
  return(0);
}
PetscErrorCode Update_Metrics_PBC(UserCtx *user)
{// updates metrics  at ghost nodes for periodci bounadry condition
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt      i,j,k;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 

  Cmpnts        ***icsi,***jcsi,***kcsi;
  Cmpnts        ***ieta,***jeta,***keta;
  Cmpnts        ***izet,***jzet,***kzet;
 
  PetscReal     ***iaj,***jaj,***kaj; 

  // PetscPrintf(PETSC_COMM_SELF, " mx=%d my=%d mz=%d \n",mx,my,mz);
 
 
  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lIZet, &izet);
  DMDAVecGetArray(fda, user->lJZet, &jzet);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  if (user->bctype[0]==7){
    if (xs==0){
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  i=0;
	  icsi[k][j][i]= icsi[k][j][i-2];
	  ieta[k][j][i]= ieta[k][j][i-2];
	  izet[k][j][i]= izet[k][j][i-2];
	  jcsi[k][j][i]= jcsi[k][j][i-2];
	  jeta[k][j][i]= jeta[k][j][i-2];
	  jzet[k][j][i]= jzet[k][j][i-2];
	  kcsi[k][j][i]= kcsi[k][j][i-2];
	  keta[k][j][i]= keta[k][j][i-2];
	  kzet[k][j][i]= kzet[k][j][i-2];
	  iaj[k][j][i]= iaj[k][j][i-2];
	  jaj[k][j][i]= jaj[k][j][i-2];
	  kaj[k][j][i]= kaj[k][j][i-2];
	  icsi[k][j][i-1]= icsi[k][j][i-3];
	  ieta[k][j][i-1]= ieta[k][j][i-3];
	  izet[k][j][i-1]= izet[k][j][i-3];
	  jcsi[k][j][i-1]= jcsi[k][j][i-3];
	  jeta[k][j][i-1]= jeta[k][j][i-3];
	  jzet[k][j][i-1]= jzet[k][j][i-3];
	  kcsi[k][j][i-1]= kcsi[k][j][i-3];
	  keta[k][j][i-1]= keta[k][j][i-3];
	  kzet[k][j][i-1]= kzet[k][j][i-3];
	  iaj[k][j][i-1]= iaj[k][j][i-3];
	  jaj[k][j][i-1]= jaj[k][j][i-3];
	  kaj[k][j][i-1]= kaj[k][j][i-3];
	}
      }
    }
  }
  if (user->bctype[0]==7){
    if (xe==mx){
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  i=mx-1;
	  icsi[k][j][i]= icsi[k][j][i+2];
	  ieta[k][j][i]= ieta[k][j][i+2];
	  izet[k][j][i]= izet[k][j][i+2];
	  jcsi[k][j][i]= jcsi[k][j][i+2];
	  jeta[k][j][i]= jeta[k][j][i+2];
	  jzet[k][j][i]= jzet[k][j][i+2];
	  kcsi[k][j][i]= kcsi[k][j][i+2];
	  keta[k][j][i]= keta[k][j][i+2];
	  kzet[k][j][i]= kzet[k][j][i+2];
	  iaj[k][j][i]= iaj[k][j][i+2];
	  jaj[k][j][i]= jaj[k][j][i+2];
	  kaj[k][j][i]= kaj[k][j][i+2];
	  icsi[k][j][i+1]= icsi[k][j][i+3];
	  ieta[k][j][i+1]= ieta[k][j][i+3];
	  izet[k][j][i+1]= izet[k][j][i+3];
	  jcsi[k][j][i+1]= jcsi[k][j][i+3];
	  jeta[k][j][i+1]= jeta[k][j][i+3];
	  jzet[k][j][i+1]= jzet[k][j][i+3];
	  kcsi[k][j][i+1]= kcsi[k][j][i+3];
	  keta[k][j][i+1]= keta[k][j][i+3];
	  kzet[k][j][i+1]= kzet[k][j][i+3];
	  iaj[k][j][i+1]= iaj[k][j][i+3];
	  jaj[k][j][i+1]= jaj[k][j][i+3];
	  kaj[k][j][i+1]= kaj[k][j][i+3];
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);
 
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lIZet, &izet);
  DMDAVecGetArray(fda, user->lJZet, &jzet);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  if (user->bctype[2]==7){
    if (ys==0){
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  j=0;
	  icsi[k][j][i]= icsi[k][j-2][i];
	  ieta[k][j][i]= ieta[k][j-2][i];
	  izet[k][j][i]= izet[k][j-2][i];
	  jcsi[k][j][i]= jcsi[k][j-2][i];
	  jeta[k][j][i]= jeta[k][j-2][i];
	  jzet[k][j][i]= jzet[k][j-2][i];
	  kcsi[k][j][i]= kcsi[k][j-2][i];
	  keta[k][j][i]= keta[k][j-2][i];
	  kzet[k][j][i]= kzet[k][j-2][i];
	  iaj[k][j][i]= iaj[k][j-2][i];
	  jaj[k][j][i]= jaj[k][j-2][i];
	  kaj[k][j][i]= kaj[k][j-2][i];
	  icsi[k][j-1][i]= icsi[k][j-3][i];
	  ieta[k][j-1][i]= ieta[k][j-3][i];
	  izet[k][j-1][i]= izet[k][j-3][i];
	  jcsi[k][j-1][i]= jcsi[k][j-3][i];
	  jeta[k][j-1][i]= jeta[k][j-3][i];
	  jzet[k][j-1][i]= jzet[k][j-3][i];
	  kcsi[k][j-1][i]= kcsi[k][j-3][i];
	  keta[k][j-1][i]= keta[k][j-3][i];
	  kzet[k][j-1][i]= kzet[k][j-3][i];
	  iaj[k][j-1][i]= iaj[k][j-3][i];
	  jaj[k][j-1][i]= jaj[k][j-3][i];
	  kaj[k][j-1][i]= kaj[k][j-3][i];
	}
      }
    }
  }
  if (user->bctype[2]==7){
    if (ye==my){
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  j=my-1;
	  icsi[k][j][i]= icsi[k][j+2][i];
	  ieta[k][j][i]= ieta[k][j+2][i];
	  izet[k][j][i]= izet[k][j+2][i];
	  jcsi[k][j][i]= jcsi[k][j+2][i];
	  jeta[k][j][i]= jeta[k][j+2][i];
	  jzet[k][j][i]= jzet[k][j+2][i];
	  kcsi[k][j][i]= kcsi[k][j+2][i];
	  keta[k][j][i]= keta[k][j+2][i];
	  kzet[k][j][i]= kzet[k][j+2][i];
	  iaj[k][j][i]= iaj[k][j+2][i];
	  jaj[k][j][i]= jaj[k][j+2][i];
	  kaj[k][j][i]= kaj[k][j+2][i];
	  icsi[k][j+1][i]= icsi[k][j+3][i];
	  ieta[k][j+1][i]= ieta[k][j+3][i];
	  izet[k][j+1][i]= izet[k][j+3][i];
	  jcsi[k][j+1][i]= jcsi[k][j+3][i];
	  jeta[k][j+1][i]= jeta[k][j+3][i];
	  jzet[k][j+1][i]= jzet[k][j+3][i];
	  kcsi[k][j+1][i]= kcsi[k][j+3][i];
	  keta[k][j+1][i]= keta[k][j+3][i];
	  kzet[k][j+1][i]= kzet[k][j+3][i];
	  iaj[k][j+1][i]= iaj[k][j+3][i];
	  jaj[k][j+1][i]= jaj[k][j+3][i];
	  kaj[k][j+1][i]= kaj[k][j+3][i];
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);
 
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lIZet, &izet);
  DMDAVecGetArray(fda, user->lJZet, &jzet);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);
  if (user->bctype[4]==7){
    if (zs==0){
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  k=0;
	  icsi[k][j][i]= icsi[k-2][j][i];
	  ieta[k][j][i]= ieta[k-2][j][i];
	  izet[k][j][i]= izet[k-2][j][i];
	  jcsi[k][j][i]= jcsi[k-2][j][i];
	  jeta[k][j][i]= jeta[k-2][j][i];
	  jzet[k][j][i]= jzet[k-2][j][i];
	  kcsi[k][j][i]= kcsi[k-2][j][i];
	  keta[k][j][i]= keta[k-2][j][i];
	  kzet[k][j][i]= kzet[k-2][j][i];
	  iaj[k][j][i]= iaj[k-2][j][i];
	  jaj[k][j][i]= jaj[k-2][j][i];
	  kaj[k][j][i]= kaj[k-2][j][i];
	  icsi[k-1][j][i]= icsi[k-3][j][i];
	  ieta[k-1][j][i]= ieta[k-3][j][i];
	  izet[k-1][j][i]= izet[k-3][j][i];
	  jcsi[k-1][j][i]= jcsi[k-3][j][i];
	  jeta[k-1][j][i]= jeta[k-3][j][i];
	  jzet[k-1][j][i]= jzet[k-3][j][i];
	  kcsi[k-1][j][i]= kcsi[k-3][j][i];
	  keta[k-1][j][i]= keta[k-3][j][i];
	  kzet[k-1][j][i]= kzet[k-3][j][i];
	  iaj[k-1][j][i]= iaj[k-3][j][i];
	  jaj[k-1][j][i]= jaj[k-3][j][i];
	  kaj[k-1][j][i]= kaj[k-3][j][i];
	}
      }
    }
  }
  if (user->bctype[4]==7){
    if (ze==mz){
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  k=mz-1;
	  icsi[k][j][i]= icsi[k+2][j][i];
	  ieta[k][j][i]= ieta[k+2][j][i];
	  izet[k][j][i]= izet[k+2][j][i];
	  jcsi[k][j][i]= jcsi[k+2][j][i];
	  jeta[k][j][i]= jeta[k+2][j][i];
	  jzet[k][j][i]= jzet[k+2][j][i];
	  kcsi[k][j][i]= kcsi[k+2][j][i];
	  keta[k][j][i]= keta[k+2][j][i];
	  kzet[k][j][i]= kzet[k+2][j][i];
	  iaj[k][j][i]= iaj[k+2][j][i];
	  jaj[k][j][i]= jaj[k+2][j][i];
	  kaj[k][j][i]= kaj[k+2][j][i];
	  icsi[k+1][j][i]= icsi[k+3][j][i];
	  ieta[k+1][j][i]= ieta[k+3][j][i];
	  izet[k+1][j][i]= izet[k+3][j][i];
	  jcsi[k+1][j][i]= jcsi[k+3][j][i];
	  jeta[k+1][j][i]= jeta[k+3][j][i];
	  jzet[k+1][j][i]= jzet[k+3][j][i];
	  kcsi[k+1][j][i]= kcsi[k+3][j][i];
	  keta[k+1][j][i]= keta[k+3][j][i];
	  kzet[k+1][j][i]= kzet[k+3][j][i];
	  iaj[k+1][j][i]= iaj[k+3][j][i];
	  jaj[k+1][j][i]= jaj[k+3][j][i];
	  kaj[k+1][j][i]= kaj[k+3][j][i];
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);
 
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);
   return(0);
}
PetscErrorCode Update_U_Cont_PBC(UserCtx *user)
{// updates contravariant velocities at ghost nodes for periodci bounadry condition
  DM		fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt      i,j,k;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 

  Cmpnts        ***ucont;

 
  DMDAVecGetArray(fda, user->lUcont, &ucont);
 
  if (user->bctype[0]==7){
    if (xs==0){
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  i=0;
	  ucont[k][j][i]= ucont[k][j][i-2];
	  ucont[k][j][i-1]= ucont[k][j][i-3];
	}
      }
    }
  }
  
  if (user->bctype[0]==7){
    if (xe==mx){
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  i=mx-1;
	  ucont[k][j][i]= ucont[k][j][i+2];
	  ucont[k][j][i+1]= ucont[k][j][i+3];
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  
  if (user->bctype[2]==7){
    if (ys==0){
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  j=0;
	  ucont[k][j][i]= ucont[k][j-2][i];
	  ucont[k][j-1][i]= ucont[k][j-3][i];
	}
      }
    }
  }
  if (user->bctype[2]==7){
    if (ye==my){
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  j=my-1;
 	  ucont[k][j][i]= ucont[k][j+2][i];
	  ucont[k][j+1][i]= ucont[k][j+3][i];
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->lUcont, &ucont);

  if (user->bctype[4]==7){
    if (zs==0){
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  k=0;
	  ucont[k][j][i]= ucont[k-2][j][i];
	  ucont[k-1][j][i]= ucont[k-3][j][i];
	}
      }
    }
  }
  if (user->bctype[4]==7){
    if (ze==mz){
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  k=mz-1;
	  ucont[k][j][i]= ucont[k+2][j][i];
	  ucont[k+1][j][i]= ucont[k+3][j][i];
	}
      }
    }
  }
 
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  return(0); 
}
PetscErrorCode CopyUContPBC(UserCtx *user)
{//copies contravariant velocities from x=0 to x=mx-2 ...
  DM	        fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt      i,j,k;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 

  Cmpnts        ***ucont;

  //make sure ghost nodes are updated (0->mx->mx-2) 
 DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
 DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);


  DMDAVecGetArray(fda, user->lUcont, &ucont);
   
  if (user->bctype[0]==7){
    if (xe==mx){
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  i=mx-2;
	  ucont[k][j][i].x= ucont[k][j][i+2].x;
	  
	}
      }
    }
  }

 
  if (user->bctype[2]==7){
    if (ye==my){
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  j=my-2;
 	  ucont[k][j][i].y= ucont[k][j+2][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==7){
    if (ze==mz){
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  k=mz-2;
	  ucont[k][j][i].z= ucont[k+2][j][i].z;
	  
	}
      }
    }
  }
 
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);


 DMLocalToGlobalBegin(fda, user->lUcont, INSERT_VALUES, user->Ucont);
 DMLocalToGlobalEnd(fda, user->lUcont, INSERT_VALUES, user->Ucont);



  return(0); 
}
PetscErrorCode ImplicitSolverLHSnew05(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ***dtow, ucon;
  PetscReal	***nvert, ***lnu_t;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren, nu_t=0.;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.0;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  PetscReal     zero=1e-10;
  

 if (!user->assignedA) {

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  MatZeroEntries(user->A);

  if (dir==0) {
    DMDAVecGetArray(fda, user->lICsi, &csi);
    DMDAVecGetArray(fda, user->lIEta, &eta);
    DMDAVecGetArray(fda, user->lIZet, &zet);
    DMDAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecGetArray(fda, user->lJCsi, &csi);
    DMDAVecGetArray(fda, user->lJEta, &eta);
    DMDAVecGetArray(fda, user->lJZet, &zet);
    DMDAVecGetArray(da, user->lJAj, &aj);
  } else {
    DMDAVecGetArray(fda, user->lKCsi, &csi);
    DMDAVecGetArray(fda, user->lKEta, &eta);
    DMDAVecGetArray(fda, user->lKZet, &zet);
    DMDAVecGetArray(da, user->lKAj, &aj);
  }
 
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->psuedot, &dtow);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  if (rans || les) 
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);

  PetscInt xstr,xend,ystr,yend,zstr,zend;
  
  if (user->bctype[0]==7){
    xstr=-1;
    xend=mx-1;
  }else{
    xstr=0;
    xend=mx-2;
  }
  if (user->bctype[2]==7){
    ystr=-1;
    yend=my-1;
  }else{
    ystr=0;
    yend=my-2;
  }
  if (user->bctype[4]==7){
    zstr=-1;
    zend=mz-1;
  }else{
    zstr=0;
    zend=mz-2;
  }
  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (rans || les) nu_t=lnu_t[k][j][i];
	
	//	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0  || k==mz-1) {
	if ((i == 0 &&  user->bctype[0]!=7) || i == mx-1 || (j==0 && user->bctype[2]!=7) || j==my-1 || (k==0 && user->bctype[4]!=7) || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	  if (dir==0) {
	    if ( i> xstr && i<xend && j>0 && k>0 && (nvert[k][j][i] + nvert[k][j][i+1] < 0.1)) {
	      
	      
	      csi0 = csi[k][j][i].x;
	      csi1 = csi[k][j][i].y;
	      csi2 = csi[k][j][i].z;
	      
	      eta0 = eta[k][j][i].x;
	      eta1 = eta[k][j][i].y;
	      eta2 = eta[k][j][i].z;
	      
	      zet0 = zet[k][j][i].x;
	      zet1 = zet[k][j][i].y;
	      zet2 = zet[k][j][i].z;
	      
	      g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	      
	      g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	      
	      g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
	      
	      aj2 = aj[k][j][i]*aj[k][j][i]; 
	      aj1 = aj[k][j][i];
	      
	    // diagonal
	      val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	      val[1] = -eps*eig0;
	      val[2] = -eps*eig0;
	      val[3] = -eps*eig1;
	      val[4] = -eps*eig1;
	      val[5] = -eps*eig2;
	      val[6] = -eps*eig2;

	      if (fabs(csi0)>zero) {
		// diagonal
		val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
		
		// i-1 
		if ((i==1 && user->bctype[0]!=7) || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		  // if ((i==1) || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		  val[1]=0.;
		} else {
		  val[1] -= g11*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i-1].x)>zero)
		    val[1]-= 0.5*aj1*csi0*
		      ucont[k][j][i-1].x/
		      csi[k][j][i-1].x;				
		}
	    
		// i+1
		if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) || (nvert[k][j][i] + nvert[k][j][i+1] > 0.1)){
		  val[2]=0.;
		} else {
		  val[2] -= g11*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i+1].x)>zero)
		    val[2]+= 0.5*aj1*csi0*
		      ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
		}	  
		
		// j-1 
		//  if ((j==1) || (nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1)) {
		if ((j==1 && user->bctype[2]!=7) || (nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1)) {
		  val[3]=0.;
		} else {
		  val[3] -= g22*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j-1][i].x)>zero)
		    val[3] -= 0.5*csi0*aj1*
		      (ucont[k][j-1][i].y+ucont[k][j-1][i+1].y)
		      /csi[k][j-1][i].x;
		}
		
		//j+1
		if ((j==my-2 &&  user->bctype[2]!=7) || (j==my-3 && user->bctype[2]!=7) || (nvert[k][j][i+1] + nvert[k][j][i] > 0.1)){
		  val[4]=0.;
		} else {
		  val[4]  -= g22*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i].x)>zero)
		    val[4] += 0.5*csi0*aj1*
		      (ucont[k][j][i].y+ucont[k][j][i+1].y)
		      /csi[k][j][i].x;
		}	  
		
		// k-1 
		//	if ((k==1 ) || (nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1)) {
		if ((k==1 && user->bctype[4]!=7) || (nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1)) {
		  val[5]=0.;
		} else {
		  val[5] -= g33*aj2*(nu_t+1/Re);
		  if (fabs(csi[k-1][j][i].x)>zero)
		    val[5] -= 0.5*csi0*aj1*
		      (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z)
		      /csi[k-1][j][i].x;			
		}
		
		//k+1
		if ((k==mz-2 && user->bctype[4]!=7) || (k==mz-3 && user->bctype[4]!=7) || (nvert[k][j][i+1] + nvert[k][j][i] > 0.1)){
		  val[6]=0.;
		} else {
		  val[6] -= g11*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i].x)>zero)
		    val[6] += 0.5*csi0*aj1*
		      (ucont[k][j][i].z+ucont[k][j][i+1].z)
		      /csi[k][j][i].x;			
		}
		  
	      }//csi0
	      
	      if (fabs(csi1)>zero) {
	      // diagonal
		val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
		
		// i-1 
		if ((i==1 && user->bctype[0]!=7) || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		  val[1]=0.;
		} else {
		  val[1] -= g11*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i-1].y)>zero)
		    val[1]-= 0.5*aj1*csi1*
		      ucont[k][j][i-1].x/
		      csi[k][j][i-1].y;				
		}
		
		//i+1
		if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		  val[2]=0.;
		} else {
		  val[2] -= g11*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i+1].y)>zero)
		    val[2]+= 0.5*aj1*csi1*
		      ucont[k][j][i+1].x
		      /csi[k][j][i+1].y;
		}	  
		
		// j-1 
		if ((j==1 && user->bctype[2]!=7) || (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		  val[3]=0.;
		} else {
		  val[3] -= g22*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j-1][i].y)>zero)
		    val[3] -= 0.5*csi1*aj1*
		      (ucont[k][j-1][i  ].y+ucont[k][j-1][i+1].y)
		      /csi[k][j-1][i].y;
		}
		
		//j+1
		if ((j==my-2 && user->bctype[2]!=7) || (j==my-3 && user->bctype[2]!=7) || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		  val[4]=0.;
		} else {
		  val[4] -= g22*aj2*(nu_t+1/Re);
		  if (fabs(csi[k][j][i].y)>zero)
		    val[4] += 0.5*csi1*aj1*
		      (ucont[k][j][i].y+ucont[k][j][i+1].y)
		      /csi[k][j][i].y;
		}	  
		
		// k-1 
		if ((k==1 && user->bctype[4]!=7) || (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		  val[5]=0.;
		} else {
		  val[5] -= g33*aj2*(nu_t+1/Re);
		  if (fabs(csi[k-1][j][i].y)>zero)
		    val[5] -= 0.5*csi1*aj1*
		      (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z)
		      /csi[k-1][j][i].y;			
		}
		
		//k+1
		if ((k==mz-2 && user->bctype[4]!=7) || (k==mz-3 && user->bctype[4]!=7) || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		  val[6]=0.;
		} else {
		  val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i].y)>zero)
		  val[6] += 0.5*csi1*aj1*
		    (ucont[k][j][i].z+ucont[k][j][i+1].z)
		    /csi[k][j][i].y;			
		}	  
	      }//csi1

	    if (fabs(csi2)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 

	      if ((i==1 && user->bctype[0]!=7) || (nvert[k][j][i]+nvert[k][j][i-1]>0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i-1].z)>zero)
		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;				
	      }
	      
	      //i+1
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) || (nvert[k][j][i+2]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i+1].z)>zero)
		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }	  
	      
	      // j-1 
	      if ((j==1 && user->bctype[2]!=7) || (nvert[k][j-1][i]+nvert[k][j-1][i+1]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j-1][i].z)>zero)
		val[3] -= 0.5*csi2*aj1*
		  (ucont[k][j-1][i].y+ucont[k][j-1][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if ((j==my-2 && user->bctype[2]!=7) || (j==my-3 && user->bctype[2]!=7) || (nvert[k][j+1][i]+nvert[k][j+1][i+1]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i].z)>zero)
		val[4] += 0.5*csi2*aj1*
		  (ucont[k][j][i].y+ucont[k][j][i+1].y)
		  /csi[k][j][i].z;
	      }


	      // k-1 
	      if ((k==1 && user->bctype[4]!=7) || (nvert[k-1][j][i]+nvert[k-1][j][i+1]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k-1][j][i].z)>zero)
		val[5] -= 0.5*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z)
		  /csi[k-1][j][i].z;			
	      }
	      
	      //k+1
	      if ((k==mz-2 && user->bctype[2]!=7) || (k==mz-3 && user->bctype[2]!=7) || (nvert[k+1][j][i]+nvert[k+1][j][i+1]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(csi[k][j][i].z)>zero)
		val[6] += 0.5*csi2*aj1*
		  (ucont[k][j][i].z+ucont[k][j][i+1].z)
		  /csi[k][j][i].z;			
	      }
	    }//csi2

	  /*   idx[0]=Gidx(i  , j , k,  user); */
/* 	    idx[1]=Gidx(i-1, j , k,  user); */
/* 	    idx[2]=Gidx(i+1, j , k,  user); */
/* 	    idx[3]=Gidx(i  ,j-1, k , user); */
/* 	    idx[4]=Gidx(i  ,j+1, k , user); */
/* 	    idx[5]=Gidx(i  , j ,k-1, user); */
/* 	    idx[6]=Gidx(i  , j ,k+1, user); */
/* 	    idx[0]=Gidx(i  , j , k,  user); */
 
	    idx[0]=Gidx(i  , j , k,  user);

	    if (user->bctype[0]==7 && i==0)  idx[1]=Gidx(mx-3, j , k,  user);
	    else idx[1]=Gidx(i-1, j , k,  user);

	    if (user->bctype[0]==7 && i==mx-2 )	 idx[2]=Gidx(1, j , k,  user);  
	    else idx[2]=Gidx(i+1, j , k,  user);

	    if (user->bctype[2]==7 && j==0)  idx[3]=Gidx(i, my-3 , k,  user);
	    else idx[3]=Gidx(i  ,j-1, k , user);

	    if (user->bctype[2]==7 && j==my-2)  idx[4]=Gidx(i, 1 , k,  user);
	    else idx[4]=Gidx(i  ,j+1, k , user);

	    if (user->bctype[4]==7 && k==0)  idx[5]=Gidx(i, j , mz-3,  user);
	    else idx[5]=Gidx(i  , j ,k-1, user);
	
	    if (user->bctype[4]==7 && k==mz-2)  idx[6]=Gidx(i, j , 1,  user);
	    else idx[6]=Gidx(i  , j ,k+1, user);
	  
	    
	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	    } else {
	      val[0]=1.;
	      idx[0]=row;
	      MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	    }
	   /*  if ((i==2) && j==2 &&  (k==2 || k==1)){ */
/* 		  PetscPrintf(PETSC_COMM_SELF, "val[0] %.15le val[1] %.15le val[2] %.15le val[3] %.15le val[4] %.15le val [5] %.15le val[6] %.15le \n",val[0],val[1],val[2],val[3],val[4],val[5],val[6]); */
/* 		  PetscPrintf(PETSC_COMM_SELF, "idx[0] %d idx[1] %d idx[2] %d idx[3] %d idx[4] %d idx[5] %d idx[6] %d \n",idx[0],idx[1],idx[2],idx[3],idx[4],idx[5],idx[6]); */
/* 		} */
	  } else if (dir==1) {
	    if (j>ystr && j<yend && i>0 && k>0 && (nvert[k][j][i] + nvert[k][j+1][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	   	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

	  /*   ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+ */
/* 			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025; */
/* 	    ucon.y=fabs(ucont[k][j][i  ].y); */
      
/* 	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+ */
/* 			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25; */

/* 	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11)); */
/* 	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22)); */
/* 	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33)); */

	    // diagonal
	    val[0]=1./dtow[k][j][i].y+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(eta0)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ((i==1 && user->bctype[0]!=7) || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i-1].x)>zero)
		  val[1]-= 0.5*aj1*eta0*
		    (ucont[k][j+1][i-1].x+ucont[k][j][i-1].x)
		    /eta[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i].x)>zero)
		  val[2]+= 0.5*aj1*eta0*
		    (ucont[k][j+1][i].x+ucont[k][j][i].x)		    
		    /eta[k][j][i].x;
	      }	  

	      // j-1 
	      if ((j==1 && user->bctype[2]!=7) ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j-1][i].x)>zero)
		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if ((j==my-2 && user->bctype[2]!=7)|| (j==my-3 && user->bctype[2]!=7) || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j+1][i].x)>zero)
		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }	  


	      // k-1 
	      if ((k==1 && user->bctype[4]!=7)|| (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k-1][j][i].x)>zero)
		val[5] -= 0.5*eta0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if ((k==mz-2  && user->bctype[4]!=7) || (k==mz-3  && user->bctype[4]!=7)|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i].x)>zero)
		val[6] += 0.5*eta0*aj1*
		  (ucont[k][j][i].z+ucont[k][j+1][i].z)
		   /eta[k][j][i].x;			
	      }	  
	    }//eta0

	    if (fabs(eta1)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      //if ((i==1 ) || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		if ((i==1 && user->bctype[0]!=7) || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i-1].y)>zero)
		  val[1]-= 0.5*aj1*eta1*
		    (ucont[k][j+1][i-1].x+ucont[k][j][i-1].x)
		    /eta[k][j][i-1].y;				
	      }
	      //i+1
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) ||(nvert[k][j][i]+nvert[k][j+1][i]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i].y)>zero)
		  val[2]+= 0.5*aj1*eta1*
		    (ucont[k][j+1][i].x+ucont[k][j][i].x)		    
		    /eta[k][j][i].y;
	      }	  
	      // j-1 
	      // if ((j==1) ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
	      if ((j==1 && user->bctype[2]!=7) ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j-1][i].y)>zero)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }	
	      //j+1
	     
	      if ((j==my-2 && user->bctype[2]!=7)|| (j==my-3 && user->bctype[2]!=7) || (nvert[k][j][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j+1][i].y)>zero)
		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }	  

	      // k-1 
	      // if ((k==1)|| (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
	      if ((k==1 && user->bctype[4]!=7)|| (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k-1][j][i].y)>zero)
		val[5] -= 0.5*eta1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;			
	      }
	      //k+1
	      
	      if ((k==mz-2  && user->bctype[4]!=7) || (k==mz-3  && user->bctype[4]!=7)|| (nvert[k][j][i]+nvert[k][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i].y)>zero)
		val[6] += 0.5*eta1*aj1*
		  (ucont[k][j][i].z+ucont[k][j+1][i].z)
		   /eta[k][j][i].y;			
	      }	    
	    }//eta1

	    if (fabs(eta2)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ((i==1 && user->bctype[0]!=7) || (nvert[k][j][i-1]+nvert[k][j+1][i]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i-1].z)>zero)
		  val[1]-= 0.5*aj1*eta2*
		    (ucont[k][j+1][i-1].x+ucont[k][j][i-1].x)
		    /eta[k][j][i-1].z;				
	      }
	      //i+1
	     
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i].z)>zero)
		  val[2]+= 0.5*aj1*eta2*
		    (ucont[k][j+1][i].x+ucont[k][j][i].x)		    
		    /eta[k][j][i].z;
	      }	  
	      // j-1 
	    
	      if ((j==1 && user->bctype[2]!=7) ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j-1][i].z)>zero)
		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }	
	      //j+1
	     
	      if ((j==my-2 && user->bctype[2]!=7)|| (j==my-3 && user->bctype[2]!=7) || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j+1][i].z)>zero)
		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }	  

	      // k-1 
	     
	       if ((k==1 && user->bctype[4]!=7)|| (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k-1][j][i].z)>zero)
		val[5] -= 0.5*eta2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;			
	      }
	      //k+1
	     
	      if ((k==mz-2  && user->bctype[4]!=7) || (k==mz-3  && user->bctype[4]!=7)|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2*(nu_t+1/Re);
		if (fabs(eta[k][j][i].z)>zero)
		val[6] += 0.5*eta2*aj1*
		  (ucont[k][j][i].z+ucont[k][j+1][i].z)
		   /eta[k][j][i].z;			
	      }	      
	    }//eta2
	    
	 
	    idx[0]=Gidx(i  , j , k,  user);
 
	    if (user->bctype[0]==7 && i==0)  idx[1]=Gidx(mx-3, j , k,  user);
	    else idx[1]=Gidx(i-1, j , k,  user);

	    if (user->bctype[0]==7 && i==mx-2 )	 idx[2]=Gidx(1, j , k,  user);  
	    else idx[2]=Gidx(i+1, j , k,  user);

	    if (user->bctype[2]==7 && j==0)  idx[3]=Gidx(i, my-3 , k,  user);
	    else idx[3]=Gidx(i  ,j-1, k , user);

	    if (user->bctype[2]==7 && j==my-2)  idx[4]=Gidx(i, 1 , k,  user);
	    else idx[4]=Gidx(i  ,j+1, k , user);

	    if (user->bctype[4]==7 && k==0)  idx[5]=Gidx(i, j , mz-3,  user);
	    else idx[5]=Gidx(i  , j ,k-1, user);
	
	    if (user->bctype[4]==7 && k==mz-2)  idx[6]=Gidx(i, j , 1,  user);
	    else idx[6]=Gidx(i  , j ,k+1, user);
	  /*   idx[0]=Gidx(i  , j , k,  user); */
/* 	    idx[1]=Gidx(i-1, j , k,  user); */
/* 	    idx[2]=Gidx(i+1, j , k,  user); */
/* 	    idx[3]=Gidx(i  ,j-1, k , user); */
/* 	    idx[4]=Gidx(i  ,j+1, k , user); */
/* 	    idx[5]=Gidx(i  , j ,k-1, user); */
/* 	    idx[6]=Gidx(i  , j ,k+1, user); */

	    //  if ((j==0 || j==my-2)&& i==2 &&  k==2 ) PetscPrintf(PETSC_COMM_WORLD, "val[0] %.15le val[1] %.15le val[2] %.15le val[3] %.15le val[4] %.15le val [5] %.15le val[6] %.15le \n",val[0],val[1],val[2],val[3],val[4],val[5],val[6]);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k>zstr && k<zend && j>0 && i>0 && (nvert[k][j][i] + nvert[k+1][j][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	  	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	 	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].z+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(zet0)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ((i==1 && user->bctype[0]!=7)|| (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2*(nu_t+1/Re);
		
		if (fabs(zet[k][j][i-1].x)>zero)
		  val[1]-= 0.5*aj1*zet0*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x)
		    /zet[k][j][i-1].x;				
	      }
	      
	      //i+1
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].x)>zero)
		  val[2]+= 0.5*aj1*zet0*
		    (ucont[k][j][i].x+ucont[k+1][j][i].x)
		    /zet[k][j][i].x;
	      } 

	      // j-1 
	      if ((j==1 && user->bctype[2]!=7)|| (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j-1][i].x)>zero)
		val[3] -= 0.5*zet0*aj1*
		  (ucont[k  ][j-1][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if ((j==my-2 && user->bctype[2]!=7) || (j==my-3 && user->bctype[2]!=7) || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].x)>zero)
		val[4] += 0.5*zet0*aj1*
		  (ucont[k][j][i].y+ucont[k+1][j][i].y)
		  /zet[k][j][i].x;
	      }	  
	      // k-1 
	      if ((k==1 && user->bctype[4]!=7)|| (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k-1][j][i].x)>zero)
		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;			
	      }
	      
	      //k+1
	      if ((k==mz-2 && user->bctype[4]!=7)|| (k==mz-3 && user->bctype[4]!=7) || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].x)>zero)
		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k][j][i].x;			
	      }	  
	    }//zet0

	    if (fabs(zet1)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ((i==1 && user->bctype[0]!=7)|| (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2*(nu_t+1/Re);

		if (fabs(zet[k][j][i-1].y)>zero)
		  val[1]-= 0.5*aj1*zet1*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x)
		     /zet[k][j][i-1].y;				
	      }
	      //i+1
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].y)>zero)
		  val[2]+= 0.5*aj1*zet1*
		    (ucont[k][j][i].x+ucont[k+1][j][i].x)
		    /zet[k][j][i].y;
	      }	  
	      // j-1 
	    if ((j==1 && user->bctype[2]!=7) || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j-1][i].y)>zero)
		val[3] -= 0.5*zet1*aj1*
		  (ucont[k  ][j-1][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }	
	      //j+1
	    if ((j==my-2 && user->bctype[2]!=7) || (j==my-3 && user->bctype[2]!=7) || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].y)>zero)
		val[4] += 0.5*zet1*aj1*
		  (ucont[k][j][i].y+ucont[k+1][j][i].y)
		  /zet[k][j][i].y;
	      }	  
	      // k-1 
	    if ((k==1 && user->bctype[4]!=7)|| (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k-1][j][i].y)>zero)
		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;			
	      }
	      //k+1
	    if ((k==mz-2 && user->bctype[4]!=7)|| (k==mz-3 && user->bctype[4]!=7) || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].y)>zero)
		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k][j][i].y;			
	      }	    
	    }//zet1

	    if (fabs(zet2)>zero) {
	      // diagonal
	      val[0] +=2.*(nu_t+1/Re)*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ((i==1 && user->bctype[0]!=7)|| (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2*(nu_t+1/Re);

		if (fabs(zet[k][j][i-1].z)>zero)
		  val[1]-= 0.5*aj1*zet2*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x)
		     /zet[k][j][i-1].z;				
	      }
	      //i+1
	      if ((i==mx-3 && user->bctype[0]!=7) || (i==mx-2 && user->bctype[0]!=7) || (nvert[k+1][j][i]+nvert[k][j][i]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i+1].z)>zero)
		  val[2]+= 0.5*aj1*zet2*
		    (ucont[k][j][i].x+ucont[k+1][j][i].x)
		    /zet[k][j][i].z;
	      }	  
	      // j-1 
	      //   if ((j==1)|| (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
	      if ((j==1 && user->bctype[2]!=7)|| (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j-1][i].z)>zero)
		val[3] -= 0.5*zet2*aj1*
		  (ucont[k  ][j-1][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }	
	      //j+1
	     
	       if ((j==my-2 && user->bctype[2]!=7) || (j==my-3 && user->bctype[2]!=7) || (nvert[k+1][j][i]+nvert[k][j][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2*(nu_t+1/Re);
		if (fabs(zet[k][j][i].z)>zero)
		val[4] += 0.5*zet2*aj1*
		  (ucont[k][j][i].y+ucont[k+1][j][i].y)
		  /zet[k][j][i].z;
	      }	 

	      // k-1 
	      // if ((k==1)|| (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
	       if ((k==1 && user->bctype[4]!=7)|| (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k-1][j][i].z)>zero)
		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;			
	      }
	      //k+1
	       if ((k==mz-2 && user->bctype[4]!=7)|| (k==mz-3 && user->bctype[4]!=7) || (nvert[k][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2*(nu_t+1/Re);
		if (fabs(zet[k+1][j][i].z)>zero)
		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;			
	       }	     
	    }//zet2


	    idx[0]=Gidx(i  , j , k,  user);

	    if (user->bctype[0]==7 && i==0)  idx[1]=Gidx(mx-3, j , k,  user);
	    else idx[1]=Gidx(i-1, j , k,  user);

	    if (user->bctype[0]==7 && i==mx-2 )	 idx[2]=Gidx(1, j , k,  user);  
	    else idx[2]=Gidx(i+1, j , k,  user);

	    if (user->bctype[2]==7 && j==0)  idx[3]=Gidx(i, my-3 , k,  user);
	    else idx[3]=Gidx(i  ,j-1, k , user);

	    if (user->bctype[2]==7 && j==my-2)  idx[4]=Gidx(i, 1 , k,  user);
	    else idx[4]=Gidx(i  ,j+1, k , user);

	    if (user->bctype[4]==7 && k==0)  idx[5]=Gidx(i, j , mz-3,  user);
	    else idx[5]=Gidx(i  , j ,k-1, user);
	
	    if (user->bctype[4]==7 && k==mz-2)  idx[6]=Gidx(i, j , 1,  user);
	    else idx[6]=Gidx(i  , j ,k+1, user);
/* 	    idx[0]=Gidx(i  , j , k,  user); */
/* 	    idx[1]=Gidx(i-1, j , k,  user); */
/* 	    idx[2]=Gidx(i+1, j , k,  user); */
/* 	    idx[3]=Gidx(i  ,j-1, k , user); */
/* 	    idx[4]=Gidx(i  ,j+1, k , user); */
/* 	    idx[5]=Gidx(i  , j ,k-1, user); */
/* 	    idx[6]=Gidx(i  , j ,k+1, user); */
	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  
	   
	  } else {
	    val[0]=1.;
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	}
      }
    }
  }
  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);
  
  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DMDAVecRestoreArray(fda, user->lICsi, &csi);
    DMDAVecRestoreArray(fda, user->lIEta, &eta);
    DMDAVecRestoreArray(fda, user->lIZet, &zet);
    DMDAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecRestoreArray(fda, user->lJCsi, &csi);
    DMDAVecRestoreArray(fda, user->lJEta, &eta);
    DMDAVecRestoreArray(fda, user->lJZet, &zet);
    DMDAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DMDAVecRestoreArray(fda, user->lKCsi, &csi);
    DMDAVecRestoreArray(fda, user->lKEta, &eta);
    DMDAVecRestoreArray(fda, user->lKZet, &zet);
    DMDAVecRestoreArray(da, user->lKAj, &aj);
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);  
  DMDAVecRestoreArray(fda, user->psuedot, &dtow);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  if (rans || les) 
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  
  return(0);
}
PetscErrorCode GetPsuedoTime(UserCtx *user)
{ 
  DM da = user->da, fda = user->fda; 
  DMDALocalInfo	info = user->info;
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

  PetscInt      i, j, k;
  Cmpnts        ***pdt, ***ucont, ru, absUcont;
  Cmpnts	***csi, ***eta, ***zet;
  PetscReal     ***iaj,***jaj,***kaj;
  PetscReal     g11,g22,g33;
  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal     cfl=user->cfl, vnn=user->vnn*user->ren;

  DMDAVecGetArray(fda, user->psuedot, &pdt);
  DMDAVecGetArray(fda, user->lICsi, &csi);
  DMDAVecGetArray(fda, user->lJEta, &eta);
  DMDAVecGetArray(fda, user->lKZet, &zet);
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	
	csi0 = csi[k][j][i].x*iaj[k][j][i];
	csi1 = csi[k][j][i].y*iaj[k][j][i];
	csi2 = csi[k][j][i].z*iaj[k][j][i];
	
	eta0 = eta[k][j][i].x*jaj[k][j][i];
	eta1 = eta[k][j][i].y*jaj[k][j][i];
	eta2 = eta[k][j][i].z*jaj[k][j][i];
	
	zet0 = zet[k][j][i].x*kaj[k][j][i];
	zet1 = zet[k][j][i].y*kaj[k][j][i];
	zet2 = zet[k][j][i].z*kaj[k][j][i];
	
	g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;	   
	g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;	   	    
	g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
	
	absUcont.x=fabs(ucont[k][j][i].x*iaj[k][j][i]);
	absUcont.y=fabs(ucont[k][j][i].y*jaj[k][j][i]);
	absUcont.z=fabs(ucont[k][j][i].z*kaj[k][j][i]);

	ru.x=absUcont.x+sqrt(absUcont.x*absUcont.x+g11);
	ru.y=absUcont.y+sqrt(absUcont.y*absUcont.y+g22);
	ru.z=absUcont.z+sqrt(absUcont.z*absUcont.z+g33);
	
	if (fabs(ru.x)>1.e-20 && g11>1.e-20)	pdt[k][j][i].x= PetscMin(cfl/ru.x,vnn/g11);
	else 	pdt[k][j][i].x=0.0;
	if (fabs(ru.y)>1.e-20 && g22>1.e-20)	pdt[k][j][i].y= PetscMin(cfl/ru.y,vnn/g22);
	else 	pdt[k][j][i].y=0.0;
	if (fabs(ru.z)>1.e-20 && g33>1.e-20)	pdt[k][j][i].z= PetscMin(cfl/ru.z,vnn/g33);
	else 	pdt[k][j][i].z=0.0;
	//	if (i==0 && j==2 && k==2)	PetscPrintf(PETSC_COMM_WORLD, "csi.x %le \n",csi[k][j][i].x);
       	//if ((i==0 || i==mx-2)&& j==1 && k==1)	PetscPrintf(PETSC_COMM_WORLD, "pdt.x %le pdt.y %le pdt.z %le \n",pdt[k][j][i].x,pdt[k][j][i].y,pdt[k][j][i].z);
      }
    }
  }
  DMDAVecRestoreArray(fda, user->psuedot, &pdt);
  DMDAVecRestoreArray(fda, user->lICsi, &csi);
  DMDAVecRestoreArray(fda, user->lJEta, &eta);
  DMDAVecRestoreArray(fda, user->lKZet, &zet);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);  
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);

 
  return(0);
}

PetscErrorCode ImplicitMomentumSolver(UserCtx *user, IBMNodes *ibm,
				      FSInfo *fsi)
{
  DM            da, fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];  alfa[0]=1.;//0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec           B[block_number];
  Vec           lUcont_i[block_number],Ucont_i[block_number], RB[block_number];

  Cmpnts        ***b , ***ucont;
  PetscReal     ***rb, ***ucont_i,***lucont_i;
  PetscReal     ***nvert;
  PetscReal     ts,te,cput,cfl_i;
 
  PetscInt	bi, ibi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp;
 
  PetscInt      norm,N;


  if (block_number>1) {
    Block_Interface_U(user);
    Block_Interface_Correction(user);
    if(blank) {
     
      Block_Blank_Correction_adv(&user[0],0);
    }
  }
 
  for (bi=0; bi<block_number; bi++) {

    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]));

    if (immersed && ti>0 ) {
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
      }
    }

    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(B[bi]));
    VecDuplicate(user[bi].P, &(RB[bi]));
    VecDuplicate(user[bi].P, &(Ucont_i[bi]));
    VecDuplicate(user[bi].lP, &(lUcont_i[bi]));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);
  
  } // for bi

  PetscInt pseudot;
  PetscReal normdU=10.,normdU1[block_number],reldU=1.,normdT;
  PetscReal normdU_bk[block_number],reldU_bk[block_number], normF_bk[block_number];
  // pseudo time iteration
  
  PetscTime(&ts);
  pseudot=0;
  while ((normdU>imp_atol && reldU>imp_rtol || pseudot<1) && pseudot<imp_MAX_IT) {
    pseudot++;
    for (bi=0; bi<block_number; bi++) { 

      da = user[bi].da;
      fda = user[bi].fda;
      info = user[bi].info;
      ksp = user[bi].ksp;
      cfl_i= user[bi].cfl;
      
      xs = info.xs; xe = info.xs + info.xm;
      ys = info.ys; ye = info.ys + info.ym;
      zs = info.zs; ze = info.zs + info.zm;
      mx = info.mx; my = info.my; mz = info.mz;
      
      lxs = xs; lxe = xe;
      lys = ys; lye = ye;
      lzs = zs; lze = ze;
      
      if (xs==0) lxs = xs+1;
      if (ys==0) lys = ys+1;
      if (zs==0) lzs = zs+1;
      
      if (xe==mx) lxe = xe-1;
      if (ye==my) lye = ye-1;
      if (ze==mz) lze = ze-1;

      for (istage=0; istage<1; istage++) {

	// Calc the RHS
	FormFunction1(&(user[bi]),user[bi].Rhs);


	// Calculate du/dt
	// du = 1.5 u - 2 u_o + 0.5 u_rm1
	if (COEF_TIME_ACCURACY>1.){// && ti!=tistart && ti!=tistart+1 && ti!=1) {
	  VecCopy(user[bi].Ucont, user[bi].dUcont);
	  VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
	  VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
	  VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);
	} else {
	  VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);
	}

	// Add -du/dt to right hand side
	VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont);
	// Calc B=U+dt*RHS 

	// Calc B=RHS
	VecCopy(user[bi].Rhs, B[bi]);

	VecSet(user[bi].dUcont, 0.);


	Update_Metrics_PBC(&(user[bi]));
	
	Update_U_Cont_PBC(&(user[bi]));
	
	GetPsuedoTime(&(user[bi]));

	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){
          if(visflg)  PetscPrintf(PETSC_COMM_WORLD," Implicit Momentum Solver: stage - %d, direction - %d \n",istage,dir);
	  // Set the LHS	 
	  ImplicitSolverLHSnew05(&(user[bi]), ibm,  Ucont_i[bi], dir, alfa[istage]);
	 
	  // set rhs of solver
	  DMDAVecGetArray(fda, B[bi], &b);
	  DMDAVecGetArray(da, RB[bi], &rb);
	  DMDAVecGetArray(da, user[bi].lNvert, &nvert);
	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		} else {
		  
		  rb[k][j][i]=b[k][j][i].z;
		  
		}
	      }
	    }
	  }
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1 || ((i>mx-3 || i<1) && user[bi].bctype[0]!=7) ) {
		    rb[k][j][i]=0.;
		    b[k][j][i].x=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].x;
		  }
		} else if (dir==1) {
		  if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1 || ((j>my-3 || j<1)&& user[bi].bctype[2]!=7) ) {
		    rb[k][j][i]=0.;
		    b[k][j][i].y=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].y;
		  }
		} else {
		  // if ((j==2) && i==2 && k==2) PetscPrintf(PETSC_COMM_WORLD, "@i=%d  j=%d k=%d b.z is %.15le \n",i,j,k,rb[k][j][i]);
		  if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1 || ((k>mz-3  || k<1)&& user[bi].bctype[4]!=7)) {
		    rb[k][j][i]=0.;
		    b[k][j][i].z=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].z;
		    //rb[k][j][i]=0.0;	
		  }
		}
	      }
	    }
	  }

	  DMDAVecRestoreArray(fda, B[bi], &b);
	  DMDAVecRestoreArray(da, RB[bi], &rb);
	  DMDAVecRestoreArray(da, user[bi].lNvert, &nvert);

	  //  VecView(RB[bi],PETSC_VIEWER_STDOUT_WORLD);
	 
	  // Set KSP options
	
	  KSPSetOperators(ksp, user[bi].A, user[bi].A);
	  KSPSetType(ksp, KSPBCGSL);
  
	  // Solve
	 
	  KSPSolve(ksp, RB[bi], Ucont_i[bi]);
	 
	  if(visflg) KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);

	

	  DMGlobalToLocalBegin(da, Ucont_i[bi],INSERT_VALUES,lUcont_i[bi]);
	  DMGlobalToLocalEnd(da, Ucont_i[bi],INSERT_VALUES,lUcont_i[bi]);

	  DMDAVecGetArray(da, Ucont_i[bi], &ucont_i);
	  DMDAVecGetArray(da, lUcont_i[bi], &lucont_i);
	 
	 

	  if (dir==0 && user[bi].bctype[0]==7 && xs==0) {
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		i=0;
		lucont_i[k][j][i]=lucont_i[k][j][i-2];
	      }
	    }
	  }else if (dir==1 && user[bi].bctype[2]==7 && ys==0) {
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		j=0;
		lucont_i[k][j][i]=lucont_i[k][j-2][i];
	      }
	    }
	  }else if (dir==2 && user[bi].bctype[4]==7 && zs==0){
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		k=0;
		lucont_i[k][j][i]=lucont_i[k-2][j][i];
	      }
	    }
	  }

	  DMDAVecRestoreArray(da, Ucont_i[bi], &ucont_i);
	  DMDAVecRestoreArray(da, lUcont_i[bi], &lucont_i);

	  DMLocalToGlobalBegin(da, lUcont_i[bi],INSERT_VALUES,Ucont_i[bi]);
	  DMLocalToGlobalEnd(da, lUcont_i[bi],INSERT_VALUES,Ucont_i[bi]);

	  DMDAVecGetArray(fda, user[bi].dUcont, &ucont);
	  DMDAVecGetArray(da, Ucont_i[bi], &ucont_i);

	if (dir==0) {
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=xs; i<lxe; i++) {
		if((user[bi].bctype[0]==7) || (user[bi].bctype[0]!=7 && i>0 && i<mx-2))
		ucont[k][j][i].x=ucont_i[k][j][i];
	      }
	    }
	  }
	} else if (dir==1) {
	  for (k=lzs; k<lze; k++) {
	    for (j=ys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if((user[bi].bctype[2]==7) || (user[bi].bctype[2]!=7 && j>0 && j<my-2))
		ucont[k][j][i].y=ucont_i[k][j][i];
	      }
	    }
	  }
	} else {
	  for (k=zs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if((user[bi].bctype[4]==7) || (user[bi].bctype[4]!=7 && k>0 && k<mz-2))
		ucont[k][j][i].z=ucont_i[k][j][i];
	      }
	    }
	  }
	}
	 
	  
	  DMDAVecRestoreArray(da, Ucont_i[bi], &ucont_i);
	 
	  DMDAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	  
	  
	/*   MatValid(user[bi].A, &temp); */
/* 	  if (temp) { */
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(&user[bi].A);
	    
	    //	  }
	  
	}//dir
  
	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  


	DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU_bk[bi]);
	if (pseudot==1) normdU1[bi]=normdU_bk[bi];
	if (pseudot==1)
	  reldU_bk[bi]=1.;
	else
	  reldU_bk[bi]=normdU_bk[bi]/normdU1[bi];
	
	VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
	VecNorm(B[bi], NORM_INFINITY, &normF_bk[bi]);
	
	PetscTime(&te);
	cput=te-ts;
      if(visflg) PetscPrintf(PETSC_COMM_WORLD, "Implicit Momentum Solver: stage - %d; Pseudo-timestep - %d; block - %d; Max dU -%le; Relative dU - %le \n", istage,pseudot,bi, normdU_bk[bi], reldU_bk[bi]);

	if (!rank) {
	  FILE *f;
	  char filen[80];
	 
	  sprintf(filen, "results/Converge_dU%1.1d",bi);
	  f = fopen(filen, "a");
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU_bk[bi], reldU_bk[bi], normF_bk[bi],normdT, cput);
	  fclose(f);
	}


	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));


	 if(visflg) // PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user[bi].FluxOutSum);
	 if(visflg) // PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user[bi].FluxInSum);
	
	

	
	// Apply BC

	FormBCS(&(user[bi]));

//	 if(visflg) // PetscPrintf(PETSC_COMM_WORLD, "BCS \n");

      } // istage
            
      if (immersed && ti>0) {
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  CHKMEMQ;
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	 if(visflg)   // PetscPrintf(PETSC_COMM_WORLD, "IBM INTP \n");
	  CHKMEMQ;
	}
	
      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      
    } //bi

    normdU=-100000.;reldU=-100000.;
    for (bi=0; bi<block_number; bi++) {
    normdU=PetscMax(normdU_bk[bi],normdU);
    reldU=PetscMax(reldU_bk[bi],reldU);
    }
  }
  // Destroy   
  for (bi=0; bi<block_number; bi++) {
 
    
    KSPDestroy(&user[bi].ksp);
    VecDestroy(&user[bi].Rhs);
    VecDestroy(&user[bi].dUcont);
    VecDestroy(&user[bi].pUcont);
    
    VecDestroy(&B[bi]);
    VecDestroy(&Ucont_i[bi]);
    VecDestroy(&lUcont_i[bi]);
    VecDestroy(&RB[bi]);
  } //bi
  

  
  return(0);
}

PetscErrorCode ImplicitMomentumSolver1(UserCtx *user, IBMNodes *ibm,
				      FSInfo *fsi)
{
  DM            da, fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];  alfa[0]=1.;//0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec           B;
  Vec           Ucont_i, RB;
  Cmpnts        ***b , ***ucont;
  PetscReal     ***rb, ***ucont_i;
  PetscReal     ***nvert;
  PetscReal     ts,te,cput,cfl_i;
  //Vec    dUcont, pUcont;
  PetscInt	bi, ibi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp;
  //  MatNullSpace  nullsp;
  PetscInt      norm,N;

  if (block_number>1) {
    Block_Interface_U(user);
  }

  for (bi=0; bi<block_number; bi++) {

    da = user[bi].da;
    fda = user[bi].fda;
    info = user[bi].info;
    ksp = user[bi].ksp;
    cfl_i= user[bi].cfl;
    
    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]));

    //ibm_interpolation_advanced2(&user[bi], ibm);
    if (immersed && ti>0 ) {
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
      }
    }

/*     if (!CGSolver) { */
/*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
/*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
/*     //PetscPrintf(PETSC_COMM_WORLD, "CGS!\n" ); */
/*     } */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(B));
    VecDuplicate(user[bi].P, &(RB));
    VecDuplicate(user[bi].P, &(Ucont_i));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);
  

  //if (ti>0) {
    //Get the direction 
    //0==i dir, 1==j dir, 2==k dir
/*     for (dir=0;dir<3;dir++){ */


    PetscInt pseudot,ls_itr;
    PetscReal normdU=10.,normdU_old,normdU1,reldU=1.,normdT, normF;

    // pseudo time iteration
    //for(pseudot=0; pseudot<5;pseudot++) {
    PetscTime(&ts);
    pseudot=0;
    ls_itr=0;
    while ((normdU>imp_atol && reldU>imp_rtol || pseudot<1) && pseudot<imp_MAX_IT) {
      pseudot++;
      GetPsuedoTime(&(user[bi]));
      for (istage=0; istage<1; istage++) {

	// Calc the RHS
	FormFunction1(&(user[bi]),user[bi].Rhs);
	// Calculate du/dt
	// du = 1.5 u - 2 u_o + 0.5 u_rm1
	if (COEF_TIME_ACCURACY>1.) {
	  VecCopy(user[bi].Ucont, user[bi].dUcont);
	  VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
	  VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
	  VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);
	} else {
	  VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);
	}

	// Add -du/dt to right hand side
	VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont);
	// Calc B=U+dt*RHS 
	//VecWAXPY(B, alfa[istage]* user[bi].dt * user[bi].st  , user[bi].Rhs, user[bi].Ucont_o);
	// Calc B=RHS
	VecCopy(user[bi].Rhs, B);
	//VecScale(B,alfa[istage]*user[bi].dt* user[bi].st);
	VecSet(user[bi].dUcont, 0.);
	PetscPrintf(PETSC_COMM_WORLD, "RHS done!\n" );

	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){

	  // Set the LHS	 
	  ImplicitSolverLHSnew04(&(user[bi]), ibm,  Ucont_i, dir, alfa[istage]);

	  // set rhs of solver
	  DMDAVecGetArray(fda, B, &b);
	  DMDAVecGetArray(da, RB, &rb);
	  DMDAVecGetArray(da, user[bi].lNvert, &nvert);
	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		  
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		  
		} else {
		  rb[k][j][i]=b[k][j][i].z;			  
		}
	      }
	    }
	  }
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1 || i>mx-3 || i<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].x=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].x;
		  }
		} else if (dir==1) {
		  if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1 || j>my-3 || j<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].y=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].y;
		  }
		} else {
		  if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1 || k>mz-3 || k<1) {
		    rb[k][j][i]=0.;
		    b[k][j][i].z=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].z;	
		  }
		}
	      }
	    }
	  }

	  DMDAVecRestoreArray(fda, B, &b);
	  DMDAVecRestoreArray(da, RB, &rb);
	  DMDAVecRestoreArray(da, user[bi].lNvert, &nvert);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "options!\n" );
	  // Set KSP options
	  //KSPSetFromOptions(ksp);
	  //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  KSPSetOperators(ksp, user[bi].A, user[bi].A);
	  KSPSetType(ksp, KSPBCGS);
/* 	  KSPSetType(ksp, KSPGMRES); */

/* 	  VecNorm(RB, NORM_INFINITY, &normF); */
/* 	  /\* normF=PetscMin(normF,.9); *\/ */
/* 	  normF=PetscMax(normF,1.e-16); */
/* 	  rtol=normF*normF;//PetscMin(normF*normF,1.e-1); */
/* 	  KSPSetTolerances(ksp,normF*normF,rtol,PETSC_DEFAULT,400); */
	  
	  //MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
	  //KSPSetNullSpace(ksp, nullsp);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "Solve!\n");
	  
	  // Solve
	  KSPSolve(ksp, RB, Ucont_i);
	  //PetscBarrier(PETSC_NULL);
	  //	  KSPTrueMonitor(ksp, N, norm, PETSC_NULL);
	  KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);
/* 	  VecNorm(Ucont_i, NORM_INFINITY, &normdU); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, "norm dU %le dir %d!\n",normdU,dir); */

	  // Restore
	  DMDAVecGetArray(da, Ucont_i, &ucont_i);
	  //DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
	  DMDAVecGetArray(fda, user[bi].dUcont, &ucont);
	  
	  //PetscPrintf(PETSC_COMM_SELF, "du.i=dui!\n");
	  //PetscBarrier(PETSC_NULL);
	  
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if (i<mx-2)
		    /*  PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
		    ucont[k][j][i].x=ucont_i[k][j][i];
		} else if (dir==1) {
		  if (j<my-2)
		    ucont[k][j][i].y=ucont_i[k][j][i];
		} else {
		  if (k<mz-2)
		    ucont[k][j][i].z=ucont_i[k][j][i];	    
/* 		  if ((j==97||j==98)  && (k==33 || k==32 || k==31 ||k==30 || k==29 || k==10) ){ */
/* 		    PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
/* 		  } */
		  
		}
	      }
	    }
	  }
	  //PetscPrintf(PETSC_COMM_SELF, "End Restore!\n");
	  //DMDAVecRestoreArray(fda, B, &b);
	  
	  DMDAVecRestoreArray(da, Ucont_i, &ucont_i);
	  //DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
	  DMDAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	  //DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  //DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  
	  //VecAssemblyBegin(user->Ucont);
	  //VecAssemblyEnd(user->Ucont);
	  
	  //MatNullSpaceDestroy(nullsp);
	 /*  PetscBool temp; */
/* 	  MatValid(user[bi].A, &temp); */
/* 	  if (temp) { */
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(&user[bi].A);
	   
	    //  }
	  
	}//dir
  
	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  
	
	//	PetscBarrier(PETSC_NULL);
	DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	normdU_old=normdU;
	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU1=normdU;
	if (pseudot>1) reldU=normdU/normdU1;
	VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
	VecNorm(B, NORM_INFINITY, &normF);
	PetscTime(&te);
	cput=te-ts;
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le %le %le\n", pseudot, normdU, reldU,normF, normdT, cput);

	if (!rank) {
	  FILE *f;
	  char filen[80];
	  sprintf(filen, "results/Converge_dU%1.1d",dir);
	  //sprintf(filen, "Converge_dU",dir);
	  f = fopen(filen, "a");
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU, reldU, normF,normdT, cput);
	  fclose(f);
	}

/* 	lambda=0.5; */
/* 	if (normdU>normdU_old && ls_itr<10 && pseudot>1) { */
/* 	  ls_itr++ ; */
/* 	  pseudot--; */
/* 	  normdU=normdU_old; */

/* 	  user[bi].cfl *=lambda; */

/* 	  VecCopy(user[bi].pUcont, user[bi].Ucont); */
	  
/* 	  DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); */
/* 	  DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); */
	  	  
/* 	  PetscPrintf(PETSC_COMM_WORLD, "LineSearch  %d |F| %le |Fold| %le %le %le\n",ls_itr,normdU,normdU_old,normF,lambda); */
/* 	} else if (pseudot>15 && normdU<imp_atol*5. && ls_itr>2) { */
/* 	   user[bi].cfl /=lambda; */
/* 	} */

	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
	PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum);
	PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user->FluxInSum);
	
	
	//DMDAVecRestoreArray(fda, Ucont_o, &ucont);
	//VecCopy(Ucont_o, user->Ucont);
	/*   DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	/*   DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	
	// Apply BC
	       
	FormBCS(&(user[bi]));
      } // istage
            
      if (immersed && ti>0) {
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	}

	//ibm_interpolation_advanced(&user[bi], ibm, fsi);
	//ibm_interpolation_advanced2(&user[bi], ibm);
/* 	ibm_interpolation(ibminfo, user, ibm); */
      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      

    } //psuedo time
  //  } // ti>0

/*     }//dir */


    user[bi].cfl=cfl_i;

    // Destroy   
    KSPDestroy(&user[bi].ksp);
    VecDestroy(&user[bi].Rhs);

    VecDestroy(&B);
    VecDestroy(&Ucont_i);
    VecDestroy(&RB);
    VecDestroy(&user[bi].dUcont);
    VecDestroy(&user[bi].pUcont);

  } //bi

  if (block_number>1) {
    Block_Interface_U(user);
  }

  return(0);
}

PetscErrorCode ImplicitSolverLHSnew06(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm;
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm;
  gym = info.gym;
  gzm = info.gzm;


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ***dtow, ucon;
  PetscReal	***nvert;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.01;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;
  PetscReal     zero_val=1e-16;

  // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",dir,j,k );

 if (!user->assignedA) {
   //PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

 // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

  MatZeroEntries(user->A);

  if (dir==0) {
    DMDAVecGetArray(fda, user->lICsi, &csi);
    DMDAVecGetArray(fda, user->lIEta, &eta);
    DMDAVecGetArray(fda, user->lIZet, &zet);
    DMDAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecGetArray(fda, user->lJCsi, &csi);
    DMDAVecGetArray(fda, user->lJEta, &eta);
    DMDAVecGetArray(fda, user->lJZet, &zet);
    DMDAVecGetArray(da, user->lJAj, &aj);
  } else {
    DMDAVecGetArray(fda, user->lKCsi, &csi);
    DMDAVecGetArray(fda, user->lKEta, &eta);
    DMDAVecGetArray(fda, user->lKZet, &zet);
    DMDAVecGetArray(da, user->lKAj, &aj);
  }
  //DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->psuedot, &dtow);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row =Gidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 && (nvert[k][j][i] + nvert[k][j][i+1] < 0.1)) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i];
	    aj1 = aj[k][j][i];
	    
/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    // diagonal
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else if (fabs(csi[k][j][i-1].x)>zero_val) {
		  val[1] -= g11*aj2/Re*csi0
		    /csi[k][j][i-1].x;

		  val[1]-= 0.5*aj1*csi0*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].x;
	      }
	      
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else if (fabs(csi[k][j][i+1].x)>zero_val){
		  val[2] -= g11*aj2/Re*csi0/csi[k][j][i+1].x;

		  val[2]+= 0.5*aj1*csi0*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].x;
	      }

	      // j-1
	      if (j==1 || (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else
		if (fabs(csi[k][j-1][i].x)>zero_val) {
		val[3] -= g22*aj2/Re*csi0/csi[k][j-1][i].x;

		val[3] -= 0.5*0.25*csi0*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else
		if (fabs(csi[k][j+1][i].x)>zero_val){
		val[4]  -= g22*aj2/Re*csi0/csi[k][j+1][i].x;

		val[4] += 0.5*0.25*csi0*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].x;
	      }


	      // k-1
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else
		if (fabs(csi[k-1][j][i].x)>zero_val){
		val[5] -= g33*aj2/Re*csi0/csi[k-1][j][i].x;

		val[5] -= 0.5*0.25*csi0*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].x;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else
		if (fabs(csi[k+1][j][i].x)>zero_val){
		val[6] -= g33*aj2/Re*csi0/csi[k+1][j][i].x;

		val[6] += 0.5*0.25*csi0*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].x;
	      }
	    }//csi0

	    if (fabs(csi1)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else
		if (fabs(csi[k][j][i-1].y)>zero_val){
		val[1] -= g11*aj2/Re*csi1/csi[k][j][i-1].y;

		  val[1]-= 0.5*aj1*csi1*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].y;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else
		if (fabs(csi[k][j][i+1].y)>zero_val){
		  val[2] -= g11*aj2/Re*csi1/csi[k][j][i+1].y;

		  val[2]+= 0.5*aj1*csi1*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].y;
	      }

	      // j-1
	      if (j==1|| (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else
		if (fabs(csi[k][j-1][i].y)>zero_val){
		  val[3] -= g22*aj2/Re*csi1/csi[k][j-1][i].y;

		val[3] -= 0.5*0.25*csi1*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else
		if (fabs(csi[k][j+1][i].y)>zero_val){
		val[4] -= g22*aj2/Re*csi1/csi[k][j+1][i].y;


		val[4] += 0.5*0.25*csi1*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].y;
	      }


	      // k-1
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else
		if (fabs(csi[k-1][j][i].y)>zero_val){
		val[5] -= g33*aj2/Re*csi1/csi[k-1][j][i].y;

		val[5] -= 0.5*0.25*csi1*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].y;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else
		if (fabs(csi[k+1][j][i].y)>zero_val){
		val[6] -= g33*aj2/Re*csi1/csi[k+1][j][i].y;

		val[6] += 0.5*0.25*csi1*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].y;
	      }
	    }//csi1

	    if (fabs(csi2)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if ( i==1 || (nvert[k][j][i]+nvert[k][j][i-1]>0.1)) {
		val[1]=0.;
	      } else
		if (fabs(csi[k][j][i-1].z)>zero_val){
		  val[1] -= g11*aj2/Re*csi2/csi[k][j][i-1].z;

		  val[1]-= 0.5*aj1*csi2*
		    ucont[k][j][i-1].x/
		    csi[k][j][i-1].z;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(csi[k][j][i+1].z)>zero_val){
		  val[2] -= g11*aj2/Re*csi2/csi[k][j][i+1].z;

		  val[2]+= 0.5*aj1*csi2*
		    ucont[k][j][i+1].x
		    /csi[k][j][i+1].z;
	      }

	      // j-1
	      if (j==1 || (nvert[k][j-1][i]+nvert[k][j-1][i+1]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(csi[k][j-1][i].z)>zero_val){
		val[3] -= g22*aj2/Re*csi2/csi[k][j-1][i].z;

		val[3] -= 0.5*0.25*csi2*aj1*
		  (ucont[k][j-1][i  ].y+ucont[k][j][i  ].y+
		   ucont[k][j-1][i+1].y+ucont[k][j][i+1].y)
		  /csi[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i]+nvert[k][j+1][i+1]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(csi[k][j+1][i].z)>zero_val){
		val[4] -= g22*aj2/Re*csi2/csi[k][j+1][i].z;

		val[4] += 0.5*0.25*csi2*aj1*
		  (ucont[k][j+2][i  ].y+ucont[k][j+1][i  ].y+
		   ucont[k][j+2][i+1].y+ucont[k][j+1][i+1].y)
		  /csi[k][j+1][i].z;
	      }


	      // k-1
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j][i+1]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(csi[k-1][j][i].z)>zero_val){
		val[5] -= g33*aj2/Re*csi2/csi[k-1][j][i].z;

		val[5] -= 0.5*0.25*csi2*aj1*
		  (ucont[k-1][j][i].z+ucont[k-1][j][i+1].z+
		   ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)
		  /csi[k-1][j][i].z;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j][i+1]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(csi[k+1][j][i].z)>zero_val){
		val[6] -= g33*aj2/Re*csi2/csi[k+1][j][i].z;

		val[6] += 0.5*0.25*csi2*aj1*
		  (ucont[k+2][j][i].z+ucont[k+2][j][i+1].z+
		   ucont[k+1][j][i].z+ucont[k+1][j][i+1].z)
		  /csi[k+1][j][i].z;
	      }
	    }//csi2

	    idx[0]=Gidx(i  , j , k,  user);
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 && (nvert[k][j][i] + nvert[k][j+1][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i];
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].y+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else
		if (fabs(eta[k][j][i-1].x)>zero_val){
		val[1] -= g11*aj2/Re*eta0/eta[k][j][i-1].x;

		  val[1]-= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].x;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(eta[k][j][i+1].x)>zero_val){
		val[2] -= g11*aj2/Re*eta0/eta[k][j][i+1].x;

		  val[2]+= 0.5*aj1*eta0*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].x;
	      }

	      // j-1
	      if (j==1 ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(eta[k][j-1][i].x)>zero_val){
		val[3] -= g22*aj2/Re*eta0/eta[k][j-1][i].x;

		val[3] -= 0.5*eta0*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(eta[k][j+1][i].x)>zero_val){
		val[4]  -= g22*aj2/Re*eta0/eta[k][j+1][i].x;

		val[4] += 0.5*eta0*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].x;
	      }


	      // k-1
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(eta[k-1][j][i].x)>zero_val){
		val[5] -= g33*aj2/Re*eta0/eta[k-1][j][i].x;

		val[5] -= 0.5*0.25*eta0*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].x;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(eta[k+1][j][i].x)>zero_val){
		val[6] -= g33*aj2/Re*eta0/eta[k+1][j][i].x;

		val[6] += 0.5*0.25*eta0*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].x;
	      }
	    }//eta0

	    if (fabs(eta1)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else
		if (fabs(eta[k][j][i-1].y)>zero_val){
		val[1] -= g11*aj2/Re*eta1/eta[k][j][i-1].y;

		  val[1]-= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].y;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(eta[k][j][i+1].y)>zero_val){
		val[2] -= g11*aj2/Re*eta1/eta[k][j][i+1].y;

		  val[2]+= 0.5*aj1*eta1*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].y;
	      }

	      // j-1
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(eta[k][j-1][i].y)>zero_val){
		val[3] -= g22*aj2/Re*eta1/eta[k][j-1][i].y;

		val[3] -= 0.5*eta1*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(eta[k][j+1][i].y)>zero_val){
		val[4] -= g22*aj2/Re*eta1/eta[k][j+1][i].y;

		val[4] += 0.5*eta1*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].y;
	      }


	      // k-1
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(eta[k-1][j][i].y)>zero_val){
		val[5] -= g33*aj2/Re*eta1/eta[k-1][j][i].y;

		val[5] -= 0.5*0.25*eta1*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].y;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(eta[k+1][j][i].y)>zero_val){
		val[6] -= g33*aj2/Re*eta1/eta[k+1][j][i].y;

		val[6] += 0.5*0.25*eta1*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].y;
	      }
	    }//eta1

	    if (fabs(eta2)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else
		if (fabs(eta[k][j][i-1].z)>zero_val){
		val[1] -= g11*aj2/Re*eta2/eta[k][j][i-1].z;

		  val[1]-= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i-1].x+ucont[k][j+1][i-1].x+
		     ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)
		    /eta[k][j][i-1].z;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(eta[k][j][i+1].z)>zero_val){
		val[2] -= g11*aj2/Re*eta2/eta[k][j][i+1].z;

		  val[2]+= 0.5*aj1*eta2*0.25*
		    (ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
		     ucont[k][j][i+2].x+ucont[k][j+1][i+2].x)
		    /eta[k][j][i+1].z;
	      }

	      // j-1
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(eta[k][j-1][i].z)>zero_val){
		val[3] -= g22*aj2/Re*eta2/eta[k][j-1][i].z;

		val[3] -= 0.5*eta2*aj1*
		  ucont[k][j-1][i].y
		  /eta[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(eta[k][j+1][i].z)>zero_val){
		val[4] -= g22*aj2/Re*eta2/eta[k][j+1][i].z;

		val[4] += 0.5*eta2*aj1*
		  ucont[k][j+1][i].y
		  /eta[k][j+1][i].z;
	      }


	      // k-1
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(eta[k-1][j][i].z)>zero_val){
		val[5] -= g33*aj2/Re*eta2/eta[k-1][j][i].z;

		val[5] -= 0.5*0.25*eta2*aj1*
		  (ucont[k  ][j][i].z+ucont[k  ][j+1][i].z+
		   ucont[k-1][j][i].z+ucont[k-1][j+1][i].z)
		  /eta[k-1][j][i].z;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(eta[k+1][j][i].z)>zero_val){
		val[6] -= g33*aj2/Re*eta2/eta[k+1][j][i].z;

		val[6] += 0.5*0.25*eta2*aj1*
		  (ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
		   ucont[k+2][j][i].z+ucont[k+2][j+1][i].z)
		  /eta[k+1][j][i].z;
	      }
	    }//eta2


	    idx[0]=Gidx(i  , j , k,  user);
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k<mz-2 && (nvert[k][j][i] + nvert[k+1][j][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i];
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].z+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else
		if (fabs(zet[k][j][i-1].x)>zero_val){
		val[1] = -g11*aj2/Re*zet0/zet[k][j][i-1].x;

		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(zet[k][j][i+1].x)>zero_val){
		val[2] = -g11*aj2/Re*zet0/zet[k][j][i+1].x;

		  val[2]+= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].x;
	      }

	      // j-1
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(zet[k][j-1][i].x)>zero_val){
		val[3]  = -g22*aj2/Re*zet0/zet[k][j-1][i].x;

		val[3] -= 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].x;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(zet[k][j+1][i].x)>zero_val){
		val[4]  = -g22*aj2/Re*zet0/zet[k][j+1][i].x;

		val[4] += 0.5*zet0*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].x;
	      }


	      // k-1
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(zet[k-1][j][i].x)>zero_val){
		val[5]  = -g33*aj2/Re*zet0/zet[k-1][j][i].x;

		val[5] -= 0.5*zet0*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].x;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(zet[k+1][j][i].x)>zero_val){
		val[6]  = -g33*aj2/Re*zet0/zet[k+1][j][i].x;

		val[6] += 0.5*zet0*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].x;
	      }
	    }//zet0

	    if (fabs(zet1)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else
		if (fabs(zet[k][j][i-1].y)>zero_val){
		val[1] -= g11*aj2/Re*zet1/zet[k][j][i-1].y;

		  val[1]-= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].y;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(zet[k][j][i+1].y)>zero_val){
		val[2] -= g11*aj2/Re*zet1/zet[k][j][i+1].y;

		  val[2]+= 0.5*aj1*zet1*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].y;
	      }

	      // j-1
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(zet[k][j-1][i].y)>zero_val){
		val[3] -= g22*aj2/Re*zet1/zet[k][j-1][i].y;

		val[3] -= 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].y;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(zet[k][j+1][i].y)>zero_val){
		val[4] -= g22*aj2/Re*zet1
		  /zet[k][j+1][i].y;
		val[4] += 0.5*zet1*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].y;
	      }


	      // k-1
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(zet[k-1][j][i].y)>zero_val){
		val[5] -= g33*aj2/Re*zet1/zet[k-1][j][i].y;

		val[5] -= 0.5*zet1*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].y;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(zet[k+1][j][i].y)>zero_val){
		val[6] -= g33*aj2/Re*zet1/zet[k+1][j][i].y;

		val[6] += 0.5*zet1*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].y;
	      }
	    }//zet1

	    if (fabs(zet2)>zero_val) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else
		if (fabs(zet[k][j][i-1].z)>zero_val){
		val[1] -= g11*aj2/Re*zet2/zet[k][j][i-1].z;

		  val[1]-= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].z;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else
		if (fabs(zet[k][j][i+1].z)>zero_val){
		val[2] -= g11*aj2/Re*zet2/zet[k][j][i+1].z;

		  val[2]+= 0.5*aj1*zet2*0.25*
		    (ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
		     ucont[k][j][i+2].x+ucont[k+1][j][i+2].x)
		    /zet[k][j][i+1].z;
	      }

	      // j-1
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else
		if (fabs(zet[k][j-1][i].z)>zero_val){
		val[3] -= g22*aj2/Re*zet2/zet[k][j-1][i].z;

		val[3] -= 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j][i].y+ucont[k  ][j-1][i].y+
		   ucont[k+1][j][i].y+ucont[k+1][j-1][i].y)
		  /zet[k][j-1][i].z;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else
		if (fabs(zet[k][j+1][i].z)>zero_val){
		val[4] -= g22*aj2/Re*zet2/zet[k][j+1][i].z;

		val[4] += 0.5*zet2*aj1*0.25*
		  (ucont[k  ][j+2][i].y+ucont[k  ][j+1][i].y+
		   ucont[k+1][j+2][i].y+ucont[k+1][j+1][i].y)
		  /zet[k][j+1][i].z;
	      }


	      // k-1
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else
		if (fabs(zet[k-1][j][i].z)>zero_val){
		val[5] -= g33*aj2/Re*zet2/zet[k-1][j][i].z;

		val[5] -= 0.5*zet2*aj1*
		  ucont[k-1][j][i].z
		  /zet[k-1][j][i].z;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else
		if (fabs(zet[k+1][j][i].z)>zero_val){
		val[6] -= g33*aj2/Re*zet2/zet[k+1][j][i].z;

		val[6] += 0.5*zet2*aj1*
		  ucont[k+1][j][i].z
		  /zet[k+1][j][i].z;
	      }
	    }//zet2


	    idx[0]=Gidx(i  , j , k,  user);
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	}
	
	
	}
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DMDAVecRestoreArray(fda, user->lICsi, &csi);
    DMDAVecRestoreArray(fda, user->lIEta, &eta);
    DMDAVecRestoreArray(fda, user->lIZet, &zet);
    DMDAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecRestoreArray(fda, user->lJCsi, &csi);
    DMDAVecRestoreArray(fda, user->lJEta, &eta);
    DMDAVecRestoreArray(fda, user->lJZet, &zet);
    DMDAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DMDAVecRestoreArray(fda, user->lKCsi, &csi);
    DMDAVecRestoreArray(fda, user->lKEta, &eta);
    DMDAVecRestoreArray(fda, user->lKZet, &zet);
    DMDAVecRestoreArray(da, user->lKAj, &aj);
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(fda, user->psuedot, &dtow);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}

PetscErrorCode ImplicitSolverLHS_visc(UserCtx *user, IBMNodes *ibm, Vec Ucont_i, PetscInt dir, PetscReal alfa)
{ // works for Stretched grid also! only viscous terms!
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 


  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts        ***ucont, ***dtow, ucon;
  PetscReal	***nvert;
  PetscReal     ***aj, aj1,aj2;

  PetscReal     dt=user->dt, Re=user->ren;

  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g22, g33;
  PetscReal     eig0,eig1,eig2, eps=0.0;//g13, g23,  g12,

  PetscInt	i, j, k, N;
  PetscScalar	val[7];
  PetscInt	idx[7], row;

  // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",dir,j,k );

 if (!user->assignedA) {
   //PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(Ucont_i, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 7, PETSC_NULL, 7, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

 // PetscPrintf(PETSC_COMM_WORLD, "RHS done %d %d %d!\n",i,j,k );

  MatZeroEntries(user->A);

  if (dir==0) {
    DMDAVecGetArray(fda, user->lICsi, &csi);
    DMDAVecGetArray(fda, user->lIEta, &eta);
    DMDAVecGetArray(fda, user->lIZet, &zet);
    DMDAVecGetArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecGetArray(fda, user->lJCsi, &csi);
    DMDAVecGetArray(fda, user->lJEta, &eta);
    DMDAVecGetArray(fda, user->lJZet, &zet);
    DMDAVecGetArray(da, user->lJAj, &aj);
  } else {
    DMDAVecGetArray(fda, user->lKCsi, &csi);
    DMDAVecGetArray(fda, user->lKEta, &eta);
    DMDAVecGetArray(fda, user->lKZet, &zet);
    DMDAVecGetArray(da, user->lKAj, &aj);
  }
  //DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->psuedot, &dtow);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  // Set Values
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  val[0]=1; idx[0]=row;
	  MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	} else {
	if (dir==0) {
	  if (i<mx-2 && (nvert[k][j][i] + nvert[k][j][i+1] < 0.1)) {
	    

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];
	    
/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    // diagonal 
	    
	    ucon.x=fabs(ucont[k][j][i].x);
	    ucon.y=fabs(ucont[k][j+1][i  ].y+ucont[k][j][i  ].y+
			ucont[k][j+1][i+1].y+ucont[k][j][i+1].y)*0.25;
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j][i+1].z+
			ucont[k  ][j][i].z+ucont[k  ][j][i+1].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    if (fabs(csi0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//csi0

	    if (fabs(csi1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i] + nvert[k][j][i-1] > 0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2] + nvert[k][j][i+1] > 0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1|| (nvert[k][j-1][i+1] + nvert[k][j-1][i] > 0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i+1] + nvert[k][j+1][i] > 0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1|| (nvert[k-1][j][i+1] + nvert[k-1][j][i] > 0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i+1] + nvert[k+1][j][i] > 0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//csi1

	    if (fabs(csi2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if ( i==1 || (nvert[k][j][i]+nvert[k][j][i-1]>0.1)) {
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j][i+2]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j-1][i]+nvert[k][j-1][i+1]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+1][i]+nvert[k][j+1][i+1]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j][i+1]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j][i+1]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//csi2

	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }

	} else if (dir==1) {
	  if (j<my-2 && (nvert[k][j][i] + nvert[k][j+1][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k][j+1][i+1].x+
			ucont[k][j][i  ].x+ucont[k][j+1][i  ].x)*.025;
	    ucon.y=fabs(ucont[k][j][i  ].y);
      
	    ucon.z=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].y+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(eta0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 ||(nvert[k][j][i+1]+nvert[k][j+1][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 ||(nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j][i]+nvert[k-1][j+1][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3|| (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//eta0

	    if (fabs(eta1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//eta1

	    if (fabs(eta2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k][j+1][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k][j+1][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k][j][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k][j+2][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k-1][j+1][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+1][j][i]+nvert[k+1][j+1][i]>0.1)){
		val[6]=0.; 
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//eta2


	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);
	  }
	  
	} else {
	  if (k<mz-2 && (nvert[k][j][i] + nvert[k+1][j][i] < 0.1)) {

	    csi0 = csi[k][j][i].x;
	    csi1 = csi[k][j][i].y;
	    csi2 = csi[k][j][i].z;
	    
	    eta0 = eta[k][j][i].x;
	    eta1 = eta[k][j][i].y;
	    eta2 = eta[k][j][i].z;
	    
	    zet0 = zet[k][j][i].x;
	    zet1 = zet[k][j][i].y;
	    zet2 = zet[k][j][i].z;
	    
	    g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	    //g12 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	    //g13 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	    
	    g22 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	    //g23 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	    
	    g33 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	    aj2 = aj[k][j][i]*aj[k][j][i]; 
	    aj1 = aj[k][j][i];

/* 	    for (ival=0;ival<7;ival++) { */
/* 	      val[ival]=0.; */
/* 	    } */

	    ucon.x=fabs(ucont[k][j][i+1].x+ucont[k+1][j][i+1].x+
			ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)*0.25;
	    ucon.y=fabs(ucont[k+1][j][i].z+ucont[k+1][j+1][i].z+
			ucont[k  ][j][i].z+ucont[k  ][j+1][i].z)*0.25;    
	    ucon.z=fabs(ucont[k][j][i].z);

	    eig0=aj1*(ucon.x+sqrt(ucon.x*ucon.x+g11));
	    eig1=aj1*(ucon.y+sqrt(ucon.y*ucon.y+g22));
	    eig2=aj1*(ucon.z+sqrt(ucon.z*ucon.z+g33));

	    // diagonal
	    val[0]=1./dtow[k][j][i].z+COEF_TIME_ACCURACY/dt+2*eps*(eig0+eig1+eig2);//+
	    //2./Re*aj1*(g11+g12+g13+g12+g22+g23+g13+g23+g33);
	    val[1] = -eps*eig0;
	    val[2] = -eps*eig0;
	    val[3] = -eps*eig1;
	    val[4] = -eps*eig1;
	    val[5] = -eps*eig2;
	    val[6] = -eps*eig2;

	    // diagonal 
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt;
	    //val[0]=1./dtow[k][j][i].x+COEF_TIME_ACCURACY/dt+
	    //2./Re*aj2*(g11+g12+g13+g12+g22+g23+g13+g23+g33);

	    if (fabs(zet0)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] = -g11*aj2/Re;

		if (fabs(zet[k][j][i-1].x)>1e-10)
		  val[1]-= 0.5*aj1*zet0*0.25*
		    (ucont[k][j][i-1].x+ucont[k+1][j][i-1].x+
		     ucont[k][j][i  ].x+ucont[k+1][j][i  ].x)
		    /zet[k][j][i-1].x;				
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] = -g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3]  = -g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4]  = -g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5]  = -g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6]  = -g33*aj2/Re;
	      }	  
	    }//zet0

	    if (fabs(zet1)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//zet1

	    if (fabs(zet2)>1e-10) {
	      // diagonal
	      val[0] +=2./Re*aj2*(g11+g22+g33);
	      
	      // i-1 
	      if (i==1 || (nvert[k][j][i-1]+nvert[k+1][j][i-1]>0.1)){
		val[1]=0.;
	      } else {
		val[1] -= g11*aj2/Re;
	      }
	    
	      //i+1
	      if (i==mx-3 || i==mx-2 || (nvert[k+1][j][i+1]+nvert[k][j][i+1]>0.1)){
		val[2]=0.;
	      } else {
		val[2] -= g11*aj2/Re;
	      }	  

	      // j-1 
	      if (j==1 || (nvert[k+1][j-1][i]+nvert[k][j-1][i]>0.1)) {
		val[3]=0.;
	      } else {
		val[3] -= g22*aj2/Re;
	      }
		
	      //j+1
	      if (j==my-2 || j==my-3 || (nvert[k+1][j+1][i]+nvert[k][j+1][i]>0.1)){
		val[4]=0.;
	      } else {
		val[4] -= g22*aj2/Re;
	      }	  


	      // k-1 
	      if (k==1 || (nvert[k][j][i]+nvert[k-1][j][i]>0.1)) {
		val[5]=0.;
	      } else {
		val[5] -= g33*aj2/Re;
	      }
	      
	      //k+1
	      if (k==mz-2 || k==mz-3 || (nvert[k+2][j][i]+nvert[k+1][j][i]>0.1)){
		val[6]=0.;
	      } else {
		val[6] -= g33*aj2/Re;
	      }	  
	    }//zet2


	    idx[0]=Gidx(i  , j , k,  user);           
	    idx[1]=Gidx(i-1, j , k,  user);
	    idx[2]=Gidx(i+1, j , k,  user);
	    idx[3]=Gidx(i  ,j-1, k , user);
	    idx[4]=Gidx(i  ,j+1, k , user);
	    idx[5]=Gidx(i  , j ,k-1, user);
	    idx[6]=Gidx(i  , j ,k+1, user);

	    MatSetValues(user->A, 1, &row, 7, idx, val, INSERT_VALUES);	  

	  } else {
	    val[0]=1.;///dtow[k][j][i].x+COEF_TIME_ACCURACY/dt; 
	    idx[0]=row;
	    MatSetValues(user->A, 1, &row, 1, idx, val, INSERT_VALUES);	    
	  }
	}
	
	
	}	
      }
    }
  }

  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  //MatView(user->A, PETSC_VIEWER_STDOUT_WORLD);

  //Restore
  if (dir==0) {
    DMDAVecRestoreArray(fda, user->lICsi, &csi);
    DMDAVecRestoreArray(fda, user->lIEta, &eta);
    DMDAVecRestoreArray(fda, user->lIZet, &zet);
    DMDAVecRestoreArray(da, user->lIAj, &aj);
  } else if (dir==1) {
    DMDAVecRestoreArray(fda, user->lJCsi, &csi);
    DMDAVecRestoreArray(fda, user->lJEta, &eta);
    DMDAVecRestoreArray(fda, user->lJZet, &zet);
    DMDAVecRestoreArray(da, user->lJAj, &aj);
  } else {
    DMDAVecRestoreArray(fda, user->lKCsi, &csi);
    DMDAVecRestoreArray(fda, user->lKEta, &eta);
    DMDAVecRestoreArray(fda, user->lKZet, &zet);
    DMDAVecRestoreArray(da, user->lKAj, &aj);
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);  
  DMDAVecRestoreArray(fda, user->psuedot, &dtow);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  return(0);
}



PetscErrorCode CalcRHS(UserCtx *user,FSInfo *fsi)
{
  PetscInt      i, j, k;
  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  Cmpnts        ***rhs,***ucont;
  PetscReal     ***nvert;
  Vec           dUcont;

  DM            da = user->da,fda = user->fda;
  info = user->info;
  
  xs = info.xs; xe = info.xs + info.xm;
  ys = info.ys; ye = info.ys + info.ym;
  zs = info.zs; ze = info.zs + info.zm;
  mx = info.mx; my = info.my; mz = info.mz;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscReal    c=0.8;
  PetscInt     mz_end;

  if (user->bctype[5]==8) 
    mz_end=mz-2;
  else
    mz_end=mz-3;
  

  VecDuplicate(user->Ucont, &dUcont);
  
  // Calc the RHS
  FormFunction1(user,user->Rhs);

  // Calculate du/dt
  // du = 1.5 u - 2 u_o + 0.5 u_rm1  && ti!=tistart+1
   if (COEF_TIME_ACCURACY>1.1 && ti!=tistart && ti!=1) {
 
    VecCopy(user->Ucont, dUcont);
    VecScale(dUcont, COEF_TIME_ACCURACY);
    VecAXPY(dUcont, -2.,user->Ucont_o );
    VecAXPY(dUcont, 0.5, user->Ucont_rm1);
  } else {
    VecWAXPY(dUcont, -1., user->Ucont_o, user->Ucont);
  }
 

  // Add -du/dt to right hand side
  VecAXPY(user->Rhs, -1./user->dt, dUcont);
  
  // set rhs BC
  DMDAVecGetArray(fda, user->Rhs, &rhs);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->Ucont, &ucont);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1 || (i>mx-3 && user->bctype[0]!=7)) {
	  rhs[k][j][i].x=0.;
	} 
	if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1 || (j>my-3 && user->bctype[2]!=7)) {
	  rhs[k][j][i].y=0.;
	} 
 	if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1 || (k>mz_end && user->bctype[4]!=7)) {
	  rhs[k][j][i].z=0.;
	} 
	// RHS of Characteristic BC
	if (mz_end==mz-2 && k==mz-2) {
	  rhs[k][j][i].z += c*(ucont[k][j][i].z-ucont[k-1][j][i].z);
	}	      
      }
    }
  }  

  /* i direction boundary conditions*/
 
  if (xs ==0) {
    i = 0;
    for (k=zs; k<lze; k++) {
      for (j=ys; j<lye; j++) {
	if (user->bctype[0]!=7){
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	  rhs[k][j][i].z = 0;
	}
	else if (user->bctype[0]==7){
	  if (nvert[k][j][i]+nvert[k][j][i+1]>0.1){  //Mohsen Oct 2012
	    rhs[k][j][i].x=0.;
	  }
	  if (nvert[k][j][i]+nvert[k][j+1][i]>0.1){
	    rhs[k][j][i].y=0.;
	  }
	  if (nvert[k][j][i]+nvert[k+1][j][i]>0.1){
	    rhs[k][j][i].z=0.;
	  }
	}
      }
    }
  }

  if (xe == mx) {
    for (k=zs; k<lze; k++) {
      for (j=ys; j<lye; j++) {
	if (user->bctype[0]!=7){
	  i = mx-2;
	  rhs[k][j][i].x = 0;
	}
	else if(user->bctype[0]==7){
	  if (nvert[k][j][i]+nvert[k][j][i+1]>0.1){  //Mohsen Oct 2012
	    rhs[k][j][i].x=0.;
	  }
	  if (nvert[k][j][i]+nvert[k][j+1][i]>0.1){
	    rhs[k][j][i].y=0.;
	  }
	  if (nvert[k][j][i]+nvert[k+1][j][i]>0.1){
	    rhs[k][j][i].z=0.;
	  }
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
    j=0;
    for (k=zs; k<lze; k++) {
      for (i=xs; i<lxe; i++) {
	if(user->bctype[2]!=7){
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	  rhs[k][j][i].z = 0;
	}
	else if(user->bctype[2]==7){
	  if (nvert[k][j][i]+nvert[k][j][i+1]>0.1){  //Mohsen Oct 2012
	    rhs[k][j][i].x=0.;
	  }
	  if (nvert[k][j][i]+nvert[k][j+1][i]>0.1){
	    rhs[k][j][i].y=0.;
	  }
	  if (nvert[k][j][i]+nvert[k+1][j][i]>0.1){
	    rhs[k][j][i].z=0.;
	  }
	}
      }
    }
  }
  if (ye == my) {
    for (k=zs; k<lze; k++) {
      for (i=xs; i<lxe; i++) {
	if(user->bctype[2]!=7){
	  j=my-2;
	  rhs[k][j][i].y = 0;
	}
	else if (user->bctype[2]==7){
	  if (nvert[k][j][i]+nvert[k][j][i+1]>0.1){  //Mohsen Oct 2012
	    rhs[k][j][i].x=0.;
	  }
	  if (nvert[k][j][i]+nvert[k][j+1][i]>0.1){
	    rhs[k][j][i].y=0.;
	  }
	  if (nvert[k][j][i]+nvert[k+1][j][i]>0.1){
	    rhs[k][j][i].z=0.;
	  }
	}
	j=my-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }

/*   k direction boundary conditions */
  if (zs == 0) {
    k=0;
    for (j=ys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	if (user->bctype[4]!=7){
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	  rhs[k][j][i].z = 0;
	}
	else if (user->bctype[4]==7){
	  if (nvert[k][j][i]+nvert[k][j][i+1]>0.1){  //Mohsen Oct 2012
	    rhs[k][j][i].x=0.;
	  }
	  if (nvert[k][j][i]+nvert[k][j+1][i]>0.1){
	    rhs[k][j][i].y=0.;
	  }
	  if (nvert[k][j][i]+nvert[k+1][j][i]>0.1){
	    rhs[k][j][i].z=0.;
	  }
	}
      }
    }
  }
  if (ze == mz) {
    for (j=ys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	if (user->bctype[4]!=7){
	  k=mz-2;
	  rhs[k][j][i].z = 0;
	}
	else if (user->bctype[4]==7){
	  if (nvert[k][j][i]+nvert[k][j][i+1]>0.1){  //Mohsen Oct 2012
	    rhs[k][j][i].x=0.;
	  }
	  if (nvert[k][j][i]+nvert[k][j+1][i]>0.1){
	    rhs[k][j][i].y=0.;
	  }
	  if (nvert[k][j][i]+nvert[k+1][j][i]>0.1){
	    rhs[k][j][i].z=0.;
	  }
	}
	k=mz-1;
	rhs[k][j][i].x = 0;
	rhs[k][j][i].y = 0;
	rhs[k][j][i].z = 0;
      }
    }
  }
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (j==2 && k==2) PetscPrintf(PETSC_COMM_SELF, "at %d %d %d rhs.x is %le rhs.y is %le rhs.z is %le\n", i,j, k,rhs[k][j][i].x,rhs[k][j][i].y,rhs[k][j][i].z); */
/* /\* 	if (Gidx(i,j,k,user)*3+1 ==min_i) { *\/ */
/* /\* 	  PetscPrintf(PETSC_COMM_SELF, "Min Rhs.y  occurs at %d %d %d is %le \n", i,j, k, min_Rhs); *\/ */
/* /\* 	} *\/ */
/* /\* 	if (Gidx(i,j,k,user)*3+1 ==max_i) { *\/ */
/* /\* 	  PetscPrintf(PETSC_COMM_SELF, "Max Rhs.y  occurs at %d %d %d is %le \n", i,j, k, max_Rhs); *\/ */
/* /\* 	} *\/ */
/* /\* 	if (Gidx(i,j,k,user)*3+2 ==min_i) { *\/ */
/* /\* 	  PetscPrintf(PETSC_COMM_SELF, "Min Rhs.z  occurs at %d %d %d is %le \n", i,j, k, min_Rhs); *\/ */
/* /\* 	} *\/ */
/* /\* 	if (Gidx(i,j,k,user)*3+2 ==max_i) { *\/ */
/* /\* 	  PetscPrintf(PETSC_COMM_SELF, "Max Rhs.z  occurs at %d %d %d is %le \n", i,j, k, max_Rhs); *\/ */
/*       } */
/*     } */
/*   } */

  
  DMDAVecRestoreArray(fda, user->Rhs, &rhs);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
 
/*   PetscReal max_Rhs=0.0,min_Rhs=0.0; */
/*   PetscInt  max_i=0,min_i=0; */
  
/*   VecMax(user->Rhs,&max_i,&max_Rhs); */
/*   VecMin(user->Rhs,&min_i,&min_Rhs); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Max Rhs  is %.15le \n",max_Rhs); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Min Rhs  is %.15le \n",min_Rhs); */
  
  
 
  

  PetscBarrier(PETSC_NULL);


  VecDestroy(&dUcont);
  return(0);
}


PetscErrorCode LineSearch_QNewton(UserCtx *user, IBMNodes *ibm,
				  FSInfo *fsi)
{
  DM            da, fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];
  alfa[0] =0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Vec           Ucont_i, RB;
  Cmpnts        ***b,***ucont;
  PetscReal     ***rb, ***ucont_i;
  PetscReal     ts,te,cput;
  //Vec    dUcont, pUcont;
  PetscInt	bi, ibi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp;
  //  MatNullSpace  nullsp;
  PetscInt      norm,N;

  if (block_number>1) {
    Block_Interface_U(user);
  }

  for (bi=0; bi<block_number; bi++) {

    da = user[bi].da;
    fda = user[bi].fda;
    info = user[bi].info;
    ksp = user[bi].ksp;

    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]));

    //ibm_interpolation_advanced2(&user[bi], ibm);
    if (immersed && ti>0) {
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
    }
    }

    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    //VecDuplicate(user[bi].Ucont, &(B));
    VecDuplicate(user[bi].P, &(RB));
    VecDuplicate(user[bi].P, &(Ucont_i));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);


    PetscInt pseudot, ls_itr;
    PetscReal normdU=10.,normdU0,normdU1=1.,reldU=1., normF, normF0, relF=1.,normF1=1.; 
    PetscScalar lambda=.5;
    // pseudo time iteration
    //for(pseudot=0; pseudot<5;pseudot++) {
    PetscTime(&ts);
    pseudot=0;

    // Calc RHS
    CalcRHS(&user[bi],fsi);
    //VecCopy(user[bi].Rhs, B);

    VecNorm(user[bi].Rhs, NORM_INFINITY, &normF0);
    normF1=normF0;

    //while ((normF>imp_atol && relF>imp_rtol && normdU>imp_stol || pseudot<1) && pseudot<imp_MAX_IT) {
    while ((normdU>imp_atol && reldU>imp_rtol || pseudot<1) && pseudot<imp_MAX_IT) {
      pseudot++;
      //      for (istage=0; istage<1; istage++) {

/* 	// Calc the RHS */
/* 	CalcRHS(user, bi); */
/* 	VecCopy(user[bi].Rhs, B); */

	VecSet(user[bi].dUcont, 0.);

	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){

	  // Set the LHS	 
	  //ImplicitSolverLHSnew05(&(user[bi]), ibm,  Ucont_i, dir);
	  ImplicitSolverLHSnew04(&(user[bi]), ibm,  Ucont_i, dir, alfa[istage]);

	  // set RHS of ksp solver
	  DMDAVecGetArray(fda, user[bi].Rhs, &b);
	  DMDAVecGetArray(da, RB, &rb);

	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		  
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		  
		} else {
		  rb[k][j][i]=b[k][j][i].z;
		}
	      }
	    }
	  }

	  DMDAVecRestoreArray(fda, user[bi].Rhs, &b);
	  DMDAVecRestoreArray(da, RB, &rb);

	  //SetSolverRHS(da,fda,&RB,B,dir);

	  //PetscPrintf(PETSC_COMM_WORLD, "options!\n" );
	  // Set KSP options
	  //KSPSetFromOptions(ksp);
	  //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  KSPSetOperators(ksp, user[bi].A, user[bi].A);
	  KSPSetType(ksp, KSPBCGS);
/* 	  KSPSetType(ksp, KSPGMRES); */

	  VecNorm(RB, NORM_INFINITY, &normF);
	  normF=PetscMin(normF,.9);
	  normF=PetscMax(normF,1.e-14);
/* 	  rtol=normF*normF;//PetscMin(normF*normF,1.e-1); */
/* 	  KSPSetTolerances(ksp,normF,PETSC_DEFAULT,PETSC_DEFAULT,400); */
	  KSPSetTolerances(ksp,normF,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	  
	  //MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
	  //KSPSetNullSpace(ksp, nullsp);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "Solve!\n");
	  //	  PetscBarrier(PETSC_NULL);
	  
	  // Solve
	  KSPSolve(ksp, RB, Ucont_i);
	  //	  PetscBarrier(PETSC_NULL);
	  //	  KSPTrueMonitor(ksp, N, norm, PETSC_NULL);
	  KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);
	  /* 	  VecNorm(Ucont_i, NORM_INFINITY, &normdU); */
	  /* 	  PetscPrintf(PETSC_COMM_WORLD, "norm dU %le dir %d!\n",normdU,dir); */

	  // Restore
	  DMDAVecGetArray(da, Ucont_i, &ucont_i);
	  //DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
	  DMDAVecGetArray(fda, user[bi].dUcont, &ucont);
	  	  
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if (i<mx-2)
		    /*  PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
		    ucont[k][j][i].x=ucont_i[k][j][i];
		  else 
		    ucont[k][j][i].x=0.;
		} else if (dir==1) {
		  if (j<my-2)
		    ucont[k][j][i].y=ucont_i[k][j][i];
		  else 
		    ucont[k][j][i].y=0.;
		} else {
		  if (k<mz-2)
		    ucont[k][j][i].z=ucont_i[k][j][i];	    		  
		  else
		    ucont[k][j][i].z=0.;
		}
	      }
	    }
	  }
	  
	  DMDAVecRestoreArray(da, Ucont_i, &ucont_i);
	  //DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
	  DMDAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	  
	 /*  PetscBool temp; */
/* 	  MatValid(user[bi].A, &temp); */
/* 	  if (temp) { */
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(&user[bi].A);
	   
	    //  }
	  
	}//dir
  
	//normF=normF0+10.; // do the 1st itr


	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  
	
	DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
	PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum);
	PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user->FluxInSum);
			
	// Apply BC
	FormBCS(&(user[bi]));

	// Calc the RHS
	CalcRHS(&user[bi],fsi);
	//VecCopy(user[bi].Rhs, B);


	VecNorm(user[bi].Rhs, NORM_INFINITY, &normF);
	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU0=normdU;
	if (pseudot>1) reldU=normdU/normdU0;
	if (pseudot==1) normF0=normF;
	if (pseudot>1) relF=normF/normF0;

	lambda=0.5;
	ls_itr=0;
	if (normdU>normdU1 && ls_itr<7 && pseudot>1) {
	//while (normdU>normdU1 && ls_itr<7 && pseudot>1) {
	//while (normdU>normdU1 && normF>normF1 && ls_itr<10 && pseudot>1) {
	  ls_itr++ ;
	  if (ls_itr==1){
	    VecCopy(user[bi].pUcont, user[bi].Ucont);

	    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	  }	    

	  PetscPrintf(PETSC_COMM_WORLD, "LineSearch  %d |F| %le |Fold| %le %le %le\n",ls_itr,normF,normF1,normdU,lambda);

	  //VecScale(user[bi].dUcont, lambda);
	  //lambda = lambda*lambda;
	  //VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);

/* 	  for (istage=0; istage<4; istage++) { */
/* 	    // Calc the RHS */
/* 	    CalcRHS(&user[bi]); */

/* 	    // Advanced in time using RK scheme */
/* 	    VecWAXPY(user[bi].Ucont, alfa[istage] * user[bi].dt * user[bi].st , user[bi].Rhs, user[bi].pUcont); */
	  	
/* 	    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); */
/* 	    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); */
	    
/* 	    InflowFlux(&(user[bi])); */
/* 	    OutflowFlux(&(user[bi])); */
	    
/* 	    // Apply BC */
/* 	    FormBCS(&(user[bi]),&fsi[0]);	     */
/* 	  }// istage */

/* 	  // Calc the RHS */
/* 	  CalcRHS(&user[bi]); */
	  
/* 	  // Get the norms */
/* 	  //normF=normF1; */
/* 	  //normdU=normdU1; */

/* /\* 	  VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU); *\/ */
/* /\* 	  if (pseudot==1) normdU1=normdU; *\/ */
/* /\* 	  if (pseudot>1) reldU=normdU/normdU0; *\/ */
/* 	  VecNorm(user[bi].Rhs, NORM_INFINITY, &normF); */
	}
	
	normF1=normF;
	normdU1=normdU;

	PetscTime(&te);
	cput=te-ts;
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le %le %le\n", pseudot, normdU, reldU,normF,relF, cput);
	if (!rank) {
	  FILE *f;
	  char filen[80];
	  sprintf(filen, "results/Converge_dU%1.1d",dir);
	  //sprintf(filen, "Converge_dU",dir);
	  f = fopen(filen, "a");
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU, reldU, normF,relF, cput);
	  fclose(f);
	}

	//      } // istage
            
      if (immersed && ti>0) {
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1);
	}
      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      
    } //psuedo time
  //  } // ti>0

/*     }//dir */

  // Destroy   
    KSPDestroy(&user[bi].ksp);
    VecDestroy(&user[bi].Rhs);

    //VecDestroy(B);
    VecDestroy(&Ucont_i);
    VecDestroy(&RB);
    VecDestroy(&user[bi].dUcont);
    VecDestroy(&user[bi].pUcont);

  } //bi

  if (block_number>1) {
    Block_Interface_U(user);
  }

  // Destroy
/*   for (bi=0; bi<block_number; bi++) { */
/*     PetscBool temp; */
/*     MatValid(user[bi].A, &temp); */
/*     if (temp) { */
/*       user[bi].assignedA = PETSC_FALSE; */
/*       MatDestroy(user[bi].A); */
/*     }     */
/*   } */

  return(0);
}

PetscErrorCode ImpRK(UserCtx *user, IBMNodes *ibm, 
		     FSInfo *fsi)
{
  DM da, fda;
  PetscInt istage, i, j, k;
  PetscReal alfa[4],ts=0.0,te=0.0;
 
  PetscInt	bi, ibi,rank;
  PetscInt      pseudot;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  Cmpnts ***ucont, ***ucont_o, ***rhs, ***zet, ***ucat, ***ubcs;
  Cmpnts ***csi, ***eta;

  Cmpnts ducont;
  Vec    dUcont, pUcont;

  PetscReal  ***p, ***aj;
  alfa[0] = 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  PetscReal normdU=10.,normdU1[block_number],normdU0[block_number],reldU=1.,normdT, normF, normF0[block_number], relF=1.,rtol, normF1[block_number]; 
  PetscReal normdU_bk[block_number],reldU_bk[block_number], normF_bk[block_number], relF_bk[block_number];
  PetscReal lambda[block_number];

  PetscReal d_taw;

  PetscOptionsGetReal(PETSC_NULL, "-d_taw", &d_taw, PETSC_NULL);
 
  PetscTime(&ts);
  pseudot=0;

  if (block_number>1) {
    Block_Interface_U(user);
  }

 for (bi=0; bi<block_number; bi++) {
       OutflowFlux(&(user[bi]));
 }

  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]));
    if (immersed)
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
    }
    
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecSet(user[bi].Rhs,0.);
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    //VecCopy(user[bi].Ucont, user[bi].dUcont);
    VecCopy(user[bi].Ucont, user[bi].pUcont);


    // Calc RHS
    CalcRHS(&user[bi],fsi);
    
    VecNorm(user[bi].Rhs, NORM_INFINITY, &normF0[bi]);
    normF1[bi]=1000.;
    normdU1[bi]=1000.;
    lambda[bi]=user[bi].cfl;

  }
    
   while ((normdU>imp_atol && reldU>imp_rtol || pseudot<1) && pseudot<imp_MAX_IT) {
     pseudot++;
     for (bi=0; bi<block_number; bi++) {
       for (istage=0; istage<4; istage++) {

	 // Advanced in time using RK scheme
	 VecWAXPY(user[bi].Ucont,d_taw*alfa[istage], user[bi].Rhs, user[bi].pUcont);
	 DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	 DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	 
	 InflowFlux(&(user[bi]));
	 OutflowFlux(&(user[bi]));
	 
	 FormBCS(&(user[bi]));

	 // Calc the RHS
	 CalcRHS(&user[bi],fsi);
	 
       }//istage
       
       if (immersed)
	 for (ibi=0;ibi<NumberOfBodies;ibi++) {
	   ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	 }
       
       VecWAXPY(user[bi].dUcont, -1.  , user[bi].Ucont, user[bi].pUcont);
       
       normdU1[bi]=normdU_bk[bi];

       VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU_bk[bi]);
     
       if (pseudot==1) {
	 normdU0[bi]=normdU_bk[bi];
	 reldU_bk[bi]=1.;
       } else if (normdU0[bi]>1.e-10)
	 reldU_bk[bi]=normdU_bk[bi]/normdU0[bi];
       VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
       if (pseudot>1) normF1[bi]=normF_bk[bi];
       VecNorm(user[bi].Rhs, NORM_INFINITY, &normF_bk[bi]);
       if (pseudot==1) normF0[bi]=normF_bk[bi];
       relF_bk[bi]=normF_bk[bi]/normF0[bi];
       
       
       PetscBarrier(PETSC_NULL);
       PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d %d  %le %le %le %le %le\n", pseudot,bi, normdU_bk[bi], reldU_bk[bi],normF_bk[bi], normdT, cput);
       
       if (!rank) {
	 FILE *f;
	 char filen[80];
	 
	 sprintf(filen, "results/Converge_dU%1.1d",bi);
	 f = fopen(filen, "a");
	 PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le %le %le\n",ti, pseudot, normdU_bk[bi], reldU_bk[bi], normF_bk[bi],normdT, cput);
	 fclose(f);
       }
     } //bi
     
     normdU=-100000.;reldU=-100000.;relF=-100000.;
     for (bi=0; bi<block_number; bi++) {
       normdU=PetscMax(normdU_bk[bi],normdU);
       reldU=PetscMax(reldU_bk[bi],reldU);
       relF=PetscMax(relF_bk[bi],relF);
     }
     
     for (bi=0; bi<block_number; bi++) {
       if (pseudot>1 && normF_bk[bi]>normF1[bi] && normdU_bk[bi]>normdU1[bi] && lambda[bi]>0.01) {
	 lambda[bi]=0.5*lambda[bi];
	 pseudot--;
	 VecCopy(user[bi].pUcont, user[bi].Ucont);
	 PetscPrintf(PETSC_COMM_WORLD, "!!%d %le\n",pseudot, lambda[bi] );
	 
	 DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	 DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
       }else {
	 VecCopy(user[bi].Ucont, user[bi].pUcont);
	 normF1[bi]=normF_bk[bi];
	 normdU1[bi]=normdU_bk[bi];
       }
     } // for bi
     
     if (block_number>1) {
       
       PetscPrintf(PETSC_COMM_WORLD, "solver B _U \n");
       Block_Interface_U(user);
      
     }
     
   } //while (iteration loop)
  
   if (block_number>1) {
     Block_Interface_Correction(user);
     if(blank) {
       Block_Blank_Correction_adv(&user[0],0);
     }
   }
   for (bi=0; bi<block_number; bi++) {
     VecDestroy(&user[bi].Rhs);
     VecDestroy(&user[bi].dUcont);
     VecDestroy(&user[bi].pUcont);
   }

   PetscTime(&te);
   cput+=te-ts;
   PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------------------------------------------\n");  
   PetscPrintf(PETSC_COMM_WORLD,"TotalNSCputime: %le , CurrentIterationNSCputime: %le\n",cput,te-ts); 
   PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------------------------------------------\n"); 
    
  
  return(0);
}




PetscErrorCode ImplicitSmoother(UserCtx *user, IBMNodes *ibm, FSInfo *fsi,PetscInt ItrNo)
{
  DM da = user->da, fda = user->fda;
  PetscInt      dir, istage,i, j, k;
  PetscReal alfa[4];alfa[0] = 1.;// 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;

  DMDALocalInfo	info = user->info;
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

  Vec           B;
  Vec           Ucont_i, RB;
  Cmpnts        ***b , ***ucont;
  PetscReal     ***rb, ***ucont_i;
  PetscReal     ***nvert;
  //Vec    dUcont, pUcont;
  PetscInt	bi, ibi;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (bi=0; bi<block_number; bi++) {
    KSPCreate(PETSC_COMM_WORLD, &(user[bi].ksp));
  }

  KSP   	ksp = user[0].ksp;
  //  MatNullSpace  nullsp;
  PetscInt      norm,N;


  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));

    FormBCS(&(user[bi]));

    //ibm_interpolation_advanced2(&user[bi], ibm);
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
    }

    //    ibm_interpolation_advanced(&user[bi], ibm, fsi);

/*     if (!CGSolver) { */
/*     VecDuplicate(user[bi].Ucont, &(user[bi].Ucont_o)); */
/*     VecCopy(user[bi].Ucont, user[bi].Ucont_o); */
/*     //PetscPrintf(PETSC_COMM_WORLD, "CGS!\n" ); */
/*     } */
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecDuplicate(user[bi].Ucont, &(B));
    VecDuplicate(user[bi].P, &(RB));
    VecDuplicate(user[bi].P, &(Ucont_i));
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecDuplicate(user[bi].Ucont, &(user[bi].pUcont));
    VecCopy(user[bi].Ucont, user[bi].pUcont);
    VecSet(user[bi].dUcont, 0.);
  }

  //if (ti>0) {

  PetscInt pseudot;
  PetscReal normdU=10.,normdU1,reldU=1.,normdT;
  // pseudo time iteration
  //for(pseudot=0; pseudot<5;pseudot++) {
  pseudot=0;
  for (bi=0; bi<block_number; bi++) {    
    while ((normdU>1e-7 && reldU>1e-3) && pseudot<ItrNo) {
      pseudot++;
      GetPsuedoTime(&(user[bi]));    
      for (istage=0; istage<1; istage++) {

	// Calc the RHS
	FormFunction1(&(user[bi]),user[bi].Rhs);

	// Calculate du/dt
	VecCopy(user[bi].Ucont, user[bi].dUcont);
	VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
	VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
	VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);

	//VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);

	// Add -du/dt to right hand side
	VecAXPY(user[bi].Rhs, -1./user[bi].dt, user[bi].dUcont);
	// Calc B=U+dt*RHS 
	//VecWAXPY(B, alfa[istage]* user[bi].dt * user[bi].st  , user[bi].Rhs, user[bi].Ucont_o);
	// Calc B=RHS
	VecCopy(user[bi].Rhs, B);
	//VecScale(B,alfa[istage]*user[bi].dt* user[bi].st);
	VecSet(user[bi].dUcont, 0.);
	//PetscPrintf(PETSC_COMM_WORLD, "RHS done!\n" );
	//Get the direction 
	//0==i dir, 1==j dir, 2==k dir
	for (dir=0;dir<3;dir++){
	  DMDAVecGetArray(fda, B, &b);
	  DMDAVecGetArray(da, RB, &rb);
	  DMDAVecGetArray(da, user[bi].lNvert, &nvert);
	  for (k=zs; k<ze; k++) {
	    for (j=ys; j<ye; j++) {
	      for (i=xs; i<xe; i++) {
		if (dir==0) {
		  rb[k][j][i]=b[k][j][i].x;
		  
		} else if (dir==1) {
		  rb[k][j][i]=b[k][j][i].y;
		  
		} else {
		  rb[k][j][i]=b[k][j][i].z;			  
		}
	      }
	    }
	  }
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1) {
		    rb[k][j][i]=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].x;
		  }
		} else if (dir==1) {
		  if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1) {
		    rb[k][j][i]=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].y;
		  }
		} else {
		  if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1) {
		    rb[k][j][i]=0.;
		  } else {
		    rb[k][j][i]=b[k][j][i].z;	
		  }
		}
	      }
	    }
	  }

	  DMDAVecRestoreArray(fda, B, &b);
	  DMDAVecRestoreArray(da, RB, &rb);
	  DMDAVecRestoreArray(da, user[bi].lNvert, &nvert);

	  // Set the LHS
	  ImplicitSolverLHSnew04(&(user[bi]), ibm,  Ucont_i, dir, alfa[istage]);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "options!\n" );
	  // Set KSP options
	  //KSPSetFromOptions(ksp);
	  //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  KSPSetOperators(ksp, user[bi].A, user[bi].A);
	  KSPSetType(ksp, KSPBCGS);

	  //MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
	  //KSPSetNullSpace(ksp, nullsp);
	  
	  //PetscPrintf(PETSC_COMM_WORLD, "Solve!\n");
	  
	  // Solve
	  KSPSolve(ksp, RB, Ucont_i);
	  //	  PetscBarrier(PETSC_NULL);
	  //	  KSPTrueMonitor(ksp, N, norm, PETSC_NULL);
	  KSPMonitorTrueResidualNorm(ksp, N, norm, PETSC_NULL);

	  // Restore
	  DMDAVecGetArray(da, Ucont_i, &ucont_i);
	  //DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
	  DMDAVecGetArray(fda, user[bi].dUcont, &ucont);
	  
	  //PetscPrintf(PETSC_COMM_SELF, "du.i=dui!\n");
	  //PetscBarrier(PETSC_NULL);
	  
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		if (dir==0) {
		  if (i<mx-2)
		    /*  PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
		    ucont[k][j][i].x=ucont_i[k][j][i];
		} else if (dir==1) {
		  if (j<my-2)
		    ucont[k][j][i].y=ucont_i[k][j][i];
		} else {
		  if (k<mz-2)
		    ucont[k][j][i].z=ucont_i[k][j][i];	    
/* 		  if ((j==97||j==98)  && (k==33 || k==32 || k==31 ||k==30 || k==29 || k==10) ){ */
/* 		    PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le %le %d %d %d\n", ucont[k][j][i].z,b[k][j][i].z,ucont_i[k][j][i],i,j,k); */
/* 		  } */
		  
		}
	      }
	    }
	  }
	  //PetscPrintf(PETSC_COMM_SELF, "End Restore!\n");
	  //DMDAVecRestoreArray(fda, B, &b);
	  
	  DMDAVecRestoreArray(da, Ucont_i, &ucont_i);
	  //DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
	  DMDAVecRestoreArray(fda, user[bi].dUcont, &ucont);
	 
	 /*  PetscBool temp; */
/* 	  MatValid(user[bi].A, &temp); */
/* 	  if (temp) { */
	    user[bi].assignedA = PETSC_FALSE;
	    MatDestroy(&user[bi].A);
	    
	    //	  }
	  
	}
  
	VecWAXPY(user[bi].Ucont, 1.  , user[bi].dUcont, user[bi].pUcont);  
	
	PetscBarrier(PETSC_NULL);
	DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

	VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU);
	if (pseudot==1) normdU1=normdU;
	if (pseudot>1) reldU=normdU/normdU1;
	VecNorm(user[bi].psuedot, NORM_INFINITY, &normdT);
	PetscPrintf(PETSC_COMM_WORLD, "!!norm of dU %d  %le %le %le\n", pseudot, normdU, reldU, normdT);

/* 	if (!rank) { */
/* 	  FILE *f; */
/* 	  char filen[80]; */
/* 	  sprintf(filen, "Converge_dU"); */
/* 	  f = fopen(filen, "a"); */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %le %le %le\n",ti, pseudot, normdU, reldU, normdT); */
/* 	  fclose(f); */
/* 	} */

	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	
/* 	PetscPrintf(PETSC_COMM_WORLD, "FluxinRK %le\n", user->FluxOutSum); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "Inlet flux %le\n",user->FluxInSum); */
	
	
	//DMDAVecRestoreArray(fda, Ucont_o, &ucont);
	//VecCopy(Ucont_o, user->Ucont);
	/*   DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	/*   DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcong); */
	
	// Apply BC
	       
	FormBCS(&(user[bi]));
      } // istage
            
      if (immersed && ti>0) {
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	}

	//ibm_interpolation_advanced(&user[bi], ibm, fsi);
	//ibm_interpolation_advanced2(&user[bi], ibm);
/* 	ibm_interpolation(ibminfo, user, ibm); */
      }

      VecCopy(user[bi].Ucont, user[bi].pUcont);
      
    } //bi
  } //psuedo time
  //  } // ti>0
  
  // Destroy
  for (bi=0; bi<block_number; bi++) {
   /*  PetscBool temp; */
/*     MatValid(user[bi].A, &temp); */
/*     if (temp) { */
      user[bi].assignedA = PETSC_FALSE;
      MatDestroy(&user[bi].A);
      //  }
    
    KSPDestroy(&user[bi].ksp);
   
    VecDestroy(&user[bi].Rhs);

/*     if (!CGSolver) */
/*       VecDestroy(user[bi].Ucont_o); */

    VecDestroy(&B);
    VecDestroy(&Ucont_i);
    VecDestroy(&RB);
    VecDestroy(&user[bi].dUcont);
    VecDestroy(&user[bi].pUcont);
  }

  return(0);
}



PetscReal LinearInterpolation(PetscReal host[2], PetscReal p, PetscReal val[2])
{
  if (fabs(host[0]-host[1])>1e-9) {

  return(val[0]*(p-host[1])/(host[0]-host[1])
       + val[1]*(p-host[0])/(host[1]-host[0]));
  } else {
    PetscPrintf(PETSC_COMM_SELF, "Linear intrp Failed!!!!!!!!!!!!!\n"); 
  }
  return(0.);
}

  
PetscReal BilinearInterpolation(Cpt2D host[4], Cpt2D p, PetscReal val[4])
{
  /* Needs special ordering of host
             vkji 
     host[0]=v?00
     host[1]=v?10
     host[2]=v?01
     host[3]=v?11

  */
  PetscReal  v,x23,x01,v01,v23;

  // bilinear interpolation
  // in y
  if (fabs(host[0].y-host[1].y)>1e-9) {
  v01=val[0]*(p.y-host[1].y)/(host[0].y-host[1].y)
    + val[1]*(p.y-host[0].y)/(host[1].y-host[0].y);

  x01=host[0].x*(p.y-host[1].y)/(host[0].y-host[1].y)
    + host[1].x*(p.y-host[0].y)/(host[1].y-host[0].y);

  v23=val[2]*(p.y-host[3].y)/(host[2].y-host[3].y)
    + val[3]*(p.y-host[2].y)/(host[3].y-host[2].y);

  x23=host[2].x*(p.y-host[3].y)/(host[2].y-host[3].y)
    + host[3].x*(p.y-host[2].y)/(host[3].y-host[2].y);
  
  // in x
  if (fabs(x01-x23)>1e-9) {
  v = v01*(p.x-x23)/(x01-x23) +
      v23*(p.x-x01)/(x23-x01);
  return(v);
  } else {
   PetscPrintf(PETSC_COMM_SELF, "Bilinear intrp Failed x!!!!!!!!!!!!!\n"); 
  }

  } else {
    PetscPrintf(PETSC_COMM_SELF, "Bilinear intrp Failed y!!!!!!!!!!!!!\n"); 
  }
  return(0.);
}

PetscReal TrilinearInterpolation(Cmpnts host[8], Cmpnts p, PetscReal val[8])
{
  /* Needs special ordering of host
             vkji 
     host[0]=v000
     host[1]=v100
     host[2]=v010
     host[3]=v110
     host[4]=v001
     host[5]=v101
     host[6]=v011
     host[7]=v111
  */

  Cpt2D      bih[4], bip ;
  PetscReal  v, bval[4] ;

  bip.x=p.x;
  bip.y=p.y;
  // in z
  if (fabs(host[0].z-host[1].z)>1e-9) {
    bval[0]=val[0]*(p.z-host[1].z)/(host[0].z-host[1].z)
          + val[1]*(p.z-host[0].z)/(host[1].z-host[0].z);    
    bih[0].x=host[0].x*(p.z-host[1].z)/(host[0].z-host[1].z)
           + host[1].x*(p.z-host[0].z)/(host[1].z-host[0].z);
    bih[0].y=host[0].y*(p.z-host[1].z)/(host[0].z-host[1].z)
           + host[1].y*(p.z-host[0].z)/(host[1].z-host[0].z);
    
    bval[1]=val[2]*(p.z-host[3].z)/(host[2].z-host[3].z)
          + val[3]*(p.z-host[2].z)/(host[3].z-host[2].z);    
    bih[1].x=host[2].x*(p.z-host[3].z)/(host[2].z-host[3].z)
           + host[3].x*(p.z-host[2].z)/(host[3].z-host[2].z);
    bih[1].y=host[2].y*(p.z-host[3].z)/(host[2].z-host[3].z)
           + host[3].y*(p.z-host[2].z)/(host[3].z-host[2].z);

    bval[2]=val[4]*(p.z-host[5].z)/(host[4].z-host[5].z)
          + val[5]*(p.z-host[4].z)/(host[5].z-host[4].z);    
    bih[2].x=host[4].x*(p.z-host[5].z)/(host[4].z-host[5].z)
           + host[5].x*(p.z-host[4].z)/(host[5].z-host[4].z);
    bih[2].y=host[4].y*(p.z-host[5].z)/(host[4].z-host[5].z)
           + host[5].y*(p.z-host[4].z)/(host[5].z-host[4].z);

    bval[3]=val[6]*(p.z-host[7].z)/(host[6].z-host[7].z)
          + val[7]*(p.z-host[6].z)/(host[7].z-host[6].z);    
    bih[3].x=host[6].x*(p.z-host[7].z)/(host[6].z-host[7].z)
          + host[7].x*(p.z-host[6].z)/(host[7].z-host[6].z);
    bih[3].y=host[6].y*(p.z-host[7].z)/(host[6].z-host[7].z)
          + host[7].y*(p.z-host[6].z)/(host[7].z-host[6].z);
    
    // in y,x
    v = BilinearInterpolation(bih,bip,bval);
    return(v);
  } else {
    PetscPrintf(PETSC_COMM_SELF, "Trilinear intrp Failed!!!!!!!!!!!!!\n"); 
  }
 return(0.);
}
 
PetscErrorCode MgRestrictionP(UserCtx *user,PetscInt flevel)
{ 
  /* 6/28/06 Iman
     Restrics Pressure from the fine grid to coarse grid
     flevel   finer grid level
     clevel   coarser grid level
  */

  DM	da = user->da, fda = user->fda;

  UserCtx *user_f = user->user_f;
  DM	da_f =  user_f->da, fda_f = user_f->fda;  	
  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);

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

  PetscReal	***p_c, ***p_f, pRhost[8];
  Cmpnts        ***coor_c, ***coor_f, Rhost[8], R;

  PetscInt      i,j,k;

  // Get the vectors
  DMDAVecGetArray(fda, user->lCent,&coor_f);
  DMDAVecGetArray(fda_f, user_f->lCent,&coor_c);
  DMDAVecGetArray(da, user->lP, &p_f);
  DMDAVecGetArray(da_f, user_f->lP, &p_c);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	R=coor_c[k][j][i];

      /*  trilinear interpolation
	Vxyz = 	V000 (1 - x) (1 - y) (1 - z) +
	V100 x (1 - y) (1 - z) +
	V010 (1 - x) y (1 - z) +
	V001 (1 - x) (1 - y) z +
	V101 x (1 - y) z +
	V011 (1 - x) y z +
	V110 x y (1 - z) +
	V111 x y z
      */ 
       	
	Rhost[0]=coor_f[2*k-1][2*j-1][2*i-1];
	Rhost[1]=coor_f[2*k  ][2*j-1][2*i-1];
	Rhost[2]=coor_f[2*k-1][2*j  ][2*i-1];
	Rhost[3]=coor_f[2*k  ][2*j  ][2*i-1];
	Rhost[4]=coor_f[2*k-1][2*j-1][2*i  ];
	Rhost[5]=coor_f[2*k  ][2*j-1][2*i  ];
	Rhost[6]=coor_f[2*k-1][2*j  ][2*i  ];
	Rhost[7]=coor_f[2*k  ][2*j  ][2*i  ];

	pRhost[0]=p_f[2*k-1][2*j-1][2*i-1];
	pRhost[1]=p_f[2*k  ][2*j-1][2*i-1];
	pRhost[2]=p_f[2*k-1][2*j  ][2*i-1];
	pRhost[3]=p_f[2*k  ][2*j  ][2*i-1];
	pRhost[4]=p_f[2*k-1][2*j-1][2*i  ];
	pRhost[5]=p_f[2*k  ][2*j-1][2*i  ];
	pRhost[6]=p_f[2*k-1][2*j  ][2*i  ];
	pRhost[7]=p_f[2*k  ][2*j  ][2*i  ];

	p_c[k][j][i]= TrilinearInterpolation(Rhost,R,pRhost);

/* 	// Calc Coeff using lagrange formula */
/* 	for (n=0;n<8;n++) { */
/* 	  Rcoeff[n]=1.; */
/* 	  for (l=0;l<8;l++) { */
/* 	    if (l!=n) { */
/* 	      if (fabs(Rhost[n].x-Rhost[l].x)>1e-8)  */
/* 		Rcoeff[n] *=(R.x-Rhost[l].x)/(Rhost[n].x-Rhost[l].x); */
/* 	      if (fabs(Rhost[n].y-Rhost[l].y)>1e-8)  */
/* 		Rcoeff[n] *=(R.y-Rhost[l].y)/(Rhost[n].y-Rhost[l].y); */
/* 	      if (fabs(Rhost[n].z-Rhost[l].z)>1e-8)  */
/* 		Rcoeff[n] *=(R.z-Rhost[l].z)/(Rhost[n].z-Rhost[l].z); */
/* 	    }//if */
/* 	  }//l */
/* 	}//n */	     
/* 	p_c[k][j][i]=0.; */
/* 	for (n=0;n<8;n++) { */
/* 	  Rcoeff[n]=1/8.; */
/* 	  p_c[k][j][i]+=Rcoeff[n]*pRhost[n]; */
/* 	} */

      }
    }
  }

  // Restore vectors
  DMDAVecRestoreArray(fda, user->lCent,&coor_f);
  DMDAVecRestoreArray(fda_f, user_f->lCent,&coor_c);
  DMDAVecRestoreArray(da, user->lP, &p_f);
  DMDAVecRestoreArray(da_f, user_f->lP, &p_c);

  return(0);
}

PetscErrorCode MgGridInterpolation(PetscInt i, PetscInt j, PetscInt k,
				   PetscInt *ih, PetscInt *jh, PetscInt *kh,
				   PetscInt dir, UserCtx *user)
{
  //PetscPrintf(PETSC_COMM_WORLD, "grid intp/n");
  if ((user->isc)) {
    *ih = i;
  }
  else if (dir==0 || i%2==0){
    *ih = i/2;
  }
  else {
    *ih = i/2 + 1;
  }

  if ((user->jsc)) {
    *jh = j;
  }
  else if (dir==1 || j%2==0){
    *jh = j/2;
  }
  else {
    *jh = j/2 + 1;
  }

  if ((user->ksc)) {
    *kh = k;
  }
  else if (dir==2 || k%2==0){
    *kh = k/2;
  }
  else {
    *kh = k/2 + 1;
  }

  return 0;
}

PetscErrorCode MgInterpolationdU_old(UserCtx *user)
{
/*  Input: Coarsegrid dU (user) */
/*  Output: Finegrid dU   */
  DM	da = user->da, fda = user->fda;

  UserCtx *user_f = user->user_f;
  DM	 fda_f = user_f->fda;  	
  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);
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

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka, dir;
  Cmpnts ***dU, ***dU_f;

  if ((user->isc)) ia = 1;
  else ia = 0;

  if ((user->jsc)) ja = 1;
  else ja = 0;

  if ((user->ksc)) ka = 1;
  else ka = 0;

  DMDAVecGetArray(fda, user->dUcont, &dU);
  DMDAVecGetArray(fda_f, user_f->dUcont, &dU_f);

  // x-dir
  dir=0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);

	if (i%2==0) {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih  ].x )
		            / (2.-ka) / (2.-ja);
	} else {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih     ].x +
			   dU[kh   ][jh   ][ih+1-ia].x )
		            / (2.-ka) / (2.-ja) / 2.;
	}
      }
    }
  }

  // y-dir
  dir=1;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (j%2==0) {
	dU_f[k][j][i].y = (dU[kh   ][jh   ][ih  ].y )
		            / (2.-ka) / (2.-ia);
	} else {
	dU_f[k][j][i].y = (dU[kh   ][jh     ][ih  ].y +
			   dU[kh   ][jh+1-ja][ih  ].y )
		            / (2.-ka) / (2.-ia) / 2.;
	}
      }
    }
  }

  // z-dir
  dir=2;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (k%2==0) {
	dU_f[k][j][i].z = (dU[kh   ][jh   ][ih  ].z )
		            / (2.-ia) / (2.-ja);
	} else {
	dU_f[k][j][i].z = (dU[kh     ][jh   ][ih  ].z +
			   dU[kh+1-ka][jh   ][ih  ].z )
		            / (2.-ia) / (2.-ja) / 2.;
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->dUcont, &dU);
  DMDAVecRestoreArray(fda_f, user_f->dUcont, &dU_f);

  return(0);
}

PetscErrorCode MgInterpolationdU(UserCtx *user)
{
/*  Input: Coarsegrid dU (user) */
/*  Output: Finegrid dU   */
  DM	 fda = user->fda;

  UserCtx *user_f = user->user_f;
  DM	da_f =  user_f->da, fda_f = user_f->fda;  	
  DMDALocalInfo	info;
  DMDAGetLocalInfo(da_f, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze, llze,llye,llxe;
  lxs = xs; lxe = xe; llxe=xe;
  lys = ys; lye = ye; llye=ye;
  lzs = zs; lze = ze; llze=ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if (xe==mx) llxe = xe-2;
  if (ye==my) llye = ye-2;
  if (ze==mz) llze = ze-2;

  PetscInt i, j, k, ih, jh, kh, ia, ja, ka, dir;
  Cmpnts ***dU, ***dU_f;
  Cmpnts ***dA, ***dA_f;

  if ((user->isc)) ia = 1;
  else ia = 0;

  if ((user->jsc)) ja = 1;
  else ja = 0;

  if ((user->ksc)) ka = 1;
  else ka = 0;

  DMDAVecGetArray(fda, user->dUcont, &dU);
  DMDAVecGetArray(fda_f, user_f->dUcont, &dU_f);
  DMDAVecGetArray(fda, user->lCsi, &dA);
  DMDAVecGetArray(fda_f, user_f->lCsi, &dA_f);

  // x-dir
  dir=0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<llxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);

	if (i%2==0) {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih  ].x /
			   dA[kh   ][jh   ][ih  ].x)
		         * dA_f[k][j][i].x;
	} else {
	dU_f[k][j][i].x = (dU[kh   ][jh   ][ih     ].x /
			   dA[kh   ][jh   ][ih     ].x +
			   dU[kh   ][jh   ][ih+1-ia].x /
			   dA[kh   ][jh   ][ih+1-ia].x)
		            * dA_f[k][j][i].x / 2.;
	}
       
	if (fabs(dA_f[k][j][i].x)<1e-6)
	  PetscPrintf(PETSC_COMM_SELF, "dA_x is zero!!!! %d %d %d \n", i, j, k);

      }
    }
  }

  // y-dir
  dir=1;
  for (k=lzs; k<lze; k++) {
    for (j=ys; j<llye; j++) {
      for (i=lxs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (j%2==0) {
	dU_f[k][j][i].y = (dU[kh   ][jh   ][ih  ].y /
			   dA[kh   ][jh   ][ih  ].y)
		         * dA_f[k][j][i].y;
	} else {
	dU_f[k][j][i].y = (dU[kh   ][jh     ][ih  ].y /
			   dA[kh   ][jh     ][ih  ].y +
			   dU[kh   ][jh+1-ja][ih  ].y /
			   dA[kh   ][jh+1-ja][ih  ].y)
		         * dA_f[k][j][i].y   / 2.;
	}

	if (fabs(dA_f[k][j][i].y)<1e-6)
	  PetscPrintf(PETSC_COMM_SELF, "dA_y is zero!!!! %d %d %d \n", i, j, k);
      }
    }
  }

  // z-dir
  dir=2;
  for (k=zs; k<llze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	MgGridInterpolation(i, j, k, &ih, &jh, &kh,dir, user);
	
	if (k%2==0) {
	dU_f[k][j][i].z = (dU[kh   ][jh   ][ih  ].z /
			   dA[kh   ][jh   ][ih  ].z)
		          * dA_f[k][j][i].z;
	} else {
	dU_f[k][j][i].z = (dU[kh     ][jh   ][ih  ].z /
			   dA[kh     ][jh   ][ih  ].z +
			   dU[kh+1-ka][jh   ][ih  ].z /
			   dA[kh+1-ka][jh   ][ih  ].z)
		           * dA_f[k][j][i].z / 2.;
	}

	if (fabs(dA_f[k][j][i].z)<1e-6)
	  PetscPrintf(PETSC_COMM_SELF, "dA_z is zero!!!! %d %d %d \n", i, j, k);

      }
    }
  }

  DMDAVecRestoreArray(fda, user->dUcont, &dU);
  DMDAVecRestoreArray(fda_f, user_f->dUcont, &dU_f);
  DMDAVecRestoreArray(fda, user->lCsi, &dA);
  DMDAVecRestoreArray(fda_f, user_f->lCsi, &dA_f);

  return(0);
}

PetscErrorCode MgFieldRestriction(UserCtx *user)
{
  DM	da = user->da, fda = user->fda;

  DM	da_f = *user->da_f;

  UserCtx *user_f = user->user_f;
  DM	fda_f = user_f->fda;

  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***ucont, ***ucont_f;
  Cmpnts ***ucont_o, ***ucont_o_f;
  Cmpnts ***ucont_rm1, ***ucont_rm1_f;

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

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda_f, user_f->lUcont, &ucont_f);
  DMDAVecGetArray(fda, user->Ucont_o, &ucont_o);
  DMDAVecGetArray(fda_f, user_f->lUcont_o, &ucont_o_f);
  DMDAVecGetArray(fda, user->Ucont_rm1, &ucont_rm1);
  DMDAVecGetArray(fda_f, user_f->lUcont_rm1, &ucont_rm1_f);

  if ((user->isc)) ia = 0;
  else ia = 1;

  if ((user->jsc)) ja = 0;
  else ja = 1;

  if ((user->ksc)) ka = 0;
  else ka = 1;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);
	
	ucont[k][j][i].x = (ucont_f[kh   ][jh   ][ih  ].x +
			    ucont_f[kh-ka][jh   ][ih  ].x +
			    ucont_f[kh   ][jh-ja][ih  ].x +
			    ucont_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);

	ucont_o[k][j][i].x = (ucont_o_f[kh   ][jh   ][ih  ].x +
			      ucont_o_f[kh-ka][jh   ][ih  ].x +
			      ucont_o_f[kh   ][jh-ja][ih  ].x +
			      ucont_o_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);

	ucont_rm1[k][j][i].x = (ucont_rm1_f[kh   ][jh   ][ih  ].x +
				ucont_rm1_f[kh-ka][jh   ][ih  ].x +
				ucont_rm1_f[kh   ][jh-ja][ih  ].x +
				ucont_rm1_f[kh-ka][jh-ja][ih  ].x) / (2.-ka) / (2.-ja);
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

	ucont_o[k][j][i].y = (ucont_o_f[kh   ][jh  ][ih   ].y +
			      ucont_o_f[kh-ka][jh  ][ih   ].y +
			      ucont_o_f[kh   ][jh  ][ih-ia].y +
			      ucont_o_f[kh-ka][jh  ][ih-ia].y) / (2.-ka) / (2.-ia);

	ucont_rm1[k][j][i].y = (ucont_rm1_f[kh   ][jh  ][ih   ].y +
				ucont_rm1_f[kh-ka][jh  ][ih   ].y +
				ucont_rm1_f[kh   ][jh  ][ih-ia].y +
				ucont_rm1_f[kh-ka][jh  ][ih-ia].y) / (2.-ka) / (2.-ia);
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

	ucont_o[k][j][i].z = (ucont_o_f[kh  ][jh   ][ih   ].z +
			      ucont_o_f[kh  ][jh   ][ih-ia].z +
			      ucont_o_f[kh  ][jh-ja][ih   ].z +
			      ucont_o_f[kh  ][jh-ja][ih-ia].z) / (2.-ja) / (2.-ia);

	ucont_rm1[k][j][i].z = (ucont_rm1_f[kh  ][jh   ][ih   ].z +
				ucont_rm1_f[kh  ][jh   ][ih-ia].z +
				ucont_rm1_f[kh  ][jh-ja][ih   ].z +
				ucont_rm1_f[kh  ][jh-ja][ih-ia].z) / (2.-ja) / (2.-ia);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda_f, user_f->lUcont, &ucont_f);
  DMDAVecRestoreArray(fda, user->Ucont_o, &ucont_o);
  DMDAVecRestoreArray(fda_f, user_f->lUcont_o, &ucont_o_f);
  DMDAVecRestoreArray(fda, user->Ucont_rm1, &ucont_rm1);
  DMDAVecRestoreArray(fda_f, user_f->lUcont_rm1, &ucont_rm1_f);
  
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DMGlobalToLocalBegin(fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);
  DMGlobalToLocalEnd(fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o);

  DMGlobalToLocalBegin(fda, user->Ucont_rm1, INSERT_VALUES, user->lUcont_rm1);
  DMGlobalToLocalEnd(fda, user->Ucont_rm1, INSERT_VALUES, user->lUcont_rm1);

  VecCopy(user->Ucont, user->Ucont_MG);

  PetscReal ***p, ***p_f, ***v_f, v_tot;

  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(da_f, user_f->lP, &p_f);
  DMDAVecGetArray(da_f, user_f->lAj, &v_f);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	GridRestriction(i, j, k, &ih, &jh, &kh, user);

	v_tot =(      v_f[kh   ][jh   ][ih   ] +		
		      v_f[kh   ][jh   ][ih-ia] +		
		      v_f[kh   ][jh-ja][ih   ] +		
		      v_f[kh   ][jh-ja][ih-ia] +		
		      v_f[kh-ka][jh   ][ih   ] +		
		      v_f[kh-ka][jh   ][ih-ia] +		
		      v_f[kh-ka][jh-ja][ih   ] +		
		      v_f[kh-ka][jh-ja][ih-ia]    );

	p[k][j][i] = (p_f[kh   ][jh   ][ih   ] *
		      v_f[kh   ][jh   ][ih   ] +
		      p_f[kh   ][jh   ][ih-ia] *
		      v_f[kh   ][jh   ][ih-ia] +
		      p_f[kh   ][jh-ja][ih   ] *
		      v_f[kh   ][jh-ja][ih   ] +
		      p_f[kh   ][jh-ja][ih-ia] *
		      v_f[kh   ][jh-ja][ih-ia] +
		      p_f[kh-ka][jh   ][ih   ] *
		      v_f[kh-ka][jh   ][ih   ] +
		      p_f[kh-ka][jh   ][ih-ia] *
		      v_f[kh-ka][jh   ][ih-ia] +
		      p_f[kh-ka][jh-ja][ih   ] *
		      v_f[kh-ka][jh-ja][ih   ] +
		      p_f[kh-ka][jh-ja][ih-ia] *
		      v_f[kh-ka][jh-ja][ih-ia])/v_tot;

      }
    }
  }

  
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(da_f, user_f->lP, &p_f);
  DMDAVecRestoreArray(da_f, user_f->lAj, &v_f);
  
  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
		      
  return 0;
}

/* ==================================================================================             */
/*      Multi-Grid Cycle  */

PetscErrorCode MGCYC(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi,PetscInt level, PetscInt cycl_idx, PetscInt preItr, PetscInt postItr,PetscInt CGItr)
{

  PetscInt i;
  Vec      pUcont, dUcont;

  UserCtx *user, *user_c;

  user   = usermg->mgctx[level].user;
  user_c = usermg->mgctx[level-1].user;

  DM    fda=user->fda;

/* ==================================================================================             */
/*     Presmoothing fG */
  ImplicitSmoother(user,ibm, fsi,preItr);
  PetscPrintf(PETSC_COMM_WORLD, "PRE SMOOTHING!\n");


/* ==================================================================================             */
/*    CGC: */
/*    1) Restrict U,p to CG*/
    
    //MgRestrictionP(user, level, isc, jsc, ksc);
    MgFieldRestriction(user_c);
    //MyNvertRestriction(user, user_c);
    PetscPrintf(PETSC_COMM_WORLD, "RESTRICTION!!!\n");

/* ==================================================================================             */
/*      save Ucont  */
    VecDuplicate(user_c->lUcont, &pUcont);   
    VecDuplicate(user_c->lUcont, &dUcont);   
    VecCopy(user_c->lUcont, pUcont);

/* ==================================================================================             */  
/*      2) Solve Ucont on CG */
    if (level-1==0) { // Coarsest level: Solve to machine zero!
      //ImplicitMomentumSolver(user_c, ibm, fsi);
      ImplicitSmoother(user_c,ibm, fsi, CGItr);
      PetscPrintf(PETSC_COMM_WORLD, "CG Solver!!!\n");

    } else { // not the coarsest level: Solve by MGM
      for (i=0; i<cycl_idx; i++) {
	MGCYC(usermg, ibm, fsi,level-1,cycl_idx, preItr, postItr, CGItr);
      }
    }
  
/* ==================================================================================             */
/*      calc dUcont  */
    VecWAXPY(dUcont, -1., pUcont, user_c->lUcont);
    VecDuplicate(user_c->lUcont, &(user_c->dUcont));
    VecCopy(dUcont, user_c->dUcont);
    VecDuplicate(user->Ucont, &(user->dUcont));

/* ==================================================================================             */      
/*      3) Interpolate dU from CG to finer Grid*/
    MgInterpolationdU(user_c);  
    // Ucont=Ucont+dU on fine Grid
    VecDuplicate(user->Ucont, &(user->pUcont));
    VecCopy(user->Ucont, user->pUcont);
    VecWAXPY(user->Ucont, +1., user->dUcont, user->pUcont);  
    
    DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

    PetscPrintf(PETSC_COMM_WORLD, "INTERPOLATION!!!\n");
    
/* ==================================================================================             */
/*      Destroy */
    VecDestroy(&dUcont);
    VecDestroy(&pUcont);
    VecDestroy(&user->pUcont);
    VecDestroy(&user_c->dUcont);
    VecDestroy(&user->dUcont);


/* ==================================================================================             */
/*    Postsmoothing fG */
  ImplicitSmoother(user,ibm,fsi, postItr);
  PetscPrintf(PETSC_COMM_WORLD, "POST SMOOTHING!\n");

  
  return(0);
}



PetscReal deltf(PetscInt i, PetscInt j)
{
  if (i==j) 
    return 1.;
  else      
    return 0.;
}

PetscErrorCode ConvectiveJacobian(PetscReal ucon, PetscReal u[3],
				  PetscReal csi[3][3], PetscInt j,
				  PetscReal Jac[][3])
{
  PetscInt i,k;
  for (i=0;i<3;i++) {
    for (k=0;k<3;k++) {
      Jac[i][k]=u[i]*deltf(k,j);
      if (fabs(csi[k][i]>1.e-8))
	Jac[i][k]+=ucon/csi[k][i];
      Jac[i][k]*=0.5;
    }
  }
  return(0);
}

PetscErrorCode ConvectiveJacobian_u(PetscReal ucon, Cmpnts u,
				    Cmpnts csi, PetscReal Jac[])
{
  PetscInt i;
  Jac[0]=ucon+u.x*csi.x;
  Jac[1]=     u.x*csi.y;
  Jac[2]=     u.x*csi.z;

  Jac[3]=     u.y*csi.x;
  Jac[4]=ucon+u.y*csi.y;
  Jac[5]=     u.y*csi.z;

  Jac[6]=     u.z*csi.x;
  Jac[7]=     u.z*csi.y;
  Jac[8]=ucon+u.z*csi.z;      
  for (i=0;i<9;i++) Jac[i]*=0.5;
  return(0);
}



PetscErrorCode SetRHSBc(UserCtx *user, Vec Rhs)
{
  PetscInt      i, j, k;
  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  Cmpnts        ***rhs;
  PetscReal     ***nvert;
  DM            da = user->da,fda = user->fda;
  info = user->info;

  xs = info.xs; xe = info.xs + info.xm;
  ys = info.ys; ye = info.ys + info.ym;
  zs = info.zs; ze = info.zs + info.zm;
  mx = info.mx; my = info.my; mz = info.mz;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  // set rhs BC
  DMDAVecGetArray(fda, Rhs, &rhs);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((nvert[k][j][i]+nvert[k][j][i+1])>0.1 || (i>mx-3 && user->bctype[0]!=7) || (i<1 && user->bctype[0]!=7)) {//Mohsen July 2012
	  rhs[k][j][i].x=0.;
	} 
	if ((nvert[k][j][i]+nvert[k][j+1][i])>0.1 || (j>my-3 && user->bctype[2]!=7) || (j<1 && user->bctype[2]!=7)) {
	  rhs[k][j][i].y=0.;
	} 
	if ((nvert[k][j][i]+nvert[k+1][j][i])>0.1 || (k>mz-3 && user->bctype[4]!=7) || (k<1 && user->bctype[4]!=7)) {
	  rhs[k][j][i].z=0.;
	} 
/* 	// RHS of Characteristic BC */

      }
    }
  }  

  /* i direction boundary conditions*/
 
  if (user->bctype[0]!=7) { //Mohsen July 2012
    
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
  }

  /* j direction boundary conditions */
  if (user->bctype[2]!=7) {
 
    if (ys == 0) {
      j=0;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
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
  }

  /* k direction boundary conditions */
  if (user->bctype[4]!=7) {

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
  }

  DMDAVecRestoreArray(fda, Rhs, &rhs);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  return(0);
}

PetscErrorCode FormFunction_SNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr)
{
  Vec    dUcont;
  PetscInt ibi;
  
  UserCtx *user = (UserCtx*)ptr;
  VecCopy(Ucont, user->Ucont);
  
  DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

  //BC
  InflowFlux(user);
  FormBCS(user);
  if (immersed)
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(user, &(user->ibm[ibi]), ibi,1);
    }
  if (block_number>1) 
    Block_Interface_U(user);//this needs work!

  // Calc the RHS
  VecSet(Rhs,0);
  FormFunction1(user, Rhs);

  // Calculate du/dt
  VecDuplicate(user->Ucont, &dUcont);
  // du = 1.5 u - 2 u_o + 0.5 u_rm1  && ti!=tistart+1
  if (COEF_TIME_ACCURACY>1.1 && ti!=tistart && ti!=1) {
    VecCopy(user->Ucont, dUcont);
    VecScale(dUcont, COEF_TIME_ACCURACY);
    VecAXPY(dUcont, -2.,user->Ucont_o );
    VecAXPY(dUcont, 0.5, user->Ucont_rm1);
  } else {
    VecWAXPY(dUcont, -1., user->Ucont_o, user->Ucont);
  }
  
  // Add -du/dt to right hand side
  VecAXPY(Rhs, -1./user->dt, dUcont);

  SetRHSBc(user, Rhs);
  
 

//----Start visualization-------------------------------------------------------
  Cmpnts ***UU;
  PetscInt i,j,k,c,bi;
  bi=0;
  c=0;    
  DMDAVecGetArray(user->fda, user->Ucont, &UU);
  
  for (k=user[bi].info.zs; k<(user[bi].info.zs+user[bi].info.zm); k++) {
    for (j=user[bi].info.ys; j<(user[bi].info.ys+user[bi].info.ym); j++) {
      for (i=user[bi].info.xs; i<(user[bi].info.xs+user[bi].info.xm); i++) {
	if(i==1 && j==3 && k==3){
	  if(fabs(UU[k][j][i].x)>1.e-10 || fabs(UU[k][j][i].y)>1.e-10 || fabs(UU[k][j][i].z)>1.e-10 ){
	    c+=1;
	    // PetscPrintf(PETSC_COMM_WORLD, "%d:Ucont[%d][%d][%d]: {%le , %le , %le}\n",c,k,j,i,UU[k][j][i].x,UU[k][j][i].y,UU[k][j][i].z); 
	  }
	}
	}
    }
  }
  
  DMDAVecRestoreArray(user->fda, user->Ucont, &UU);
  //----End visualization-------------------------------------------------------
  

 VecDestroy(&dUcont);



  return 0;
}










PetscErrorCode PFormFunction_SNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr)
{
  Vec    dUcont,pdUcont;
  PetscInt ibi;
  
  UserCtx *user = (UserCtx*)ptr;

  PetscReal d_taw;
  PetscOptionsGetReal(PETSC_NULL, "-d_taw", &d_taw, PETSC_NULL);


  VecCopy(Ucont, user->Ucont);
  
  DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

  //BC
  InflowFlux(user);
  FormBCS(user);
  if (immersed)
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(user, &(user->ibm[ibi]), ibi,1);
    }
  if (block_number>1)
    Block_Interface_U(user);//this needs work!

  // Calc the RHS
  VecSet(Rhs,0);
  FormFunction1(user, Rhs);

  // Calculate du/dt
  VecDuplicate(user->Ucont, &dUcont);
  VecDuplicate(user->Ucont, &pdUcont);
  // du = 1.5 u - 2 u_o + 0.5 u_rm1  && ti!=tistart+1
  if (COEF_TIME_ACCURACY>1.1 && ti!=tistart && ti!=1) {
    VecCopy(user->pUcont, dUcont);
    VecScale(dUcont, COEF_TIME_ACCURACY);
    VecAXPY(dUcont, -2.,user->Ucont_o );
    VecAXPY(dUcont, 0.5, user->Ucont_rm1);
  } else {
    VecWAXPY(dUcont, -1., user->Ucont_o, user->pUcont);
  }
  
  // Add -du/dt to right hand side
  VecScale(dUcont,1./user->dt);
  VecWAXPY(pdUcont, -1., user->pUcont, user->Ucont);
  VecScale(pdUcont, 1./d_taw);
  VecAXPY(dUcont, 1.,pdUcont);
  VecAXPY(Rhs, -1.,dUcont);
 

  SetRHSBc(user, Rhs);
  
 

//----Start visualization-------------------------------------------------------
  Cmpnts ***UU;
  PetscInt i,j,k,c,bi;
  bi=0;
  c=0;
  DMDAVecGetArray(user->fda, user->Ucont, &UU);
  
  for (k=user[bi].info.zs; k<(user[bi].info.zs+user[bi].info.zm); k++) {
    for (j=user[bi].info.ys; j<(user[bi].info.ys+user[bi].info.ym); j++) {
      for (i=user[bi].info.xs; i<(user[bi].info.xs+user[bi].info.xm); i++) {
	if(i==1 && j==3 && k==3){
	  if(fabs(UU[k][j][i].x)>1.e-10 || fabs(UU[k][j][i].y)>1.e-10 || fabs(UU[k][j][i].z)>1.e-10 ){
	    c+=1;
	    // PetscPrintf(PETSC_COMM_WORLD, "%d:Ucont[%d][%d][%d]: {%le , %le , %le}\n",c,k,j,i,UU[k][j][i].x,UU[k][j][i].y,UU[k][j][i].z);
	  }
	}
	}
    }
  }
  
  DMDAVecRestoreArray(user->fda, user->Ucont, &UU);
  //----End visualization-------------------------------------------------------
  

 VecDestroy(&dUcont);
 VecDestroy(&pdUcont);


  return 0;
}

























int snes_created=0;

PetscErrorCode Implicit_MatrixFree(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
  
/*   DALocalInfo	info = user->info; */
/*   PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*   PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*   PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*   PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
	
  PetscInt	bi,i, ibi;
	

  int rank;
  PetscReal ts,te,cput;
  PetscTime(&ts);
	
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	
  for (bi=0; bi<block_number; bi++) {
		
    KSP ksp;
    PC pc;
		
    VecDuplicate(user[bi].Ucont, &user[bi].Rhs);
    if(!snes_created) {
      SNESCreate(PETSC_COMM_WORLD,&user[bi].snes);
    }
    SNESSetFunction(user[bi].snes,user[bi].Rhs,FormFunction_SNES,(void *)&user[bi]);
		
#if defined(PETSC_HAVE_ADIC___)
    DAGetMatrix(user->da,MATAIJ,&user[bi].J);
    ISColoring             iscoloring;
    DAGetColoring(user->da,IS_COLORING_GHOSTED,&iscoloring);
    MatSetColoring(J[bi],iscoloring);
    ISColoringDestroy(iscoloring);
    SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,SNESDAComputeJacobianWithAdic,(void *)&user[bi]);
#else
    MatCreateSNESMF(user[bi].snes, &user[bi].J);
    SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,MatMFFDComputeJacobian,(void *)&user[bi]);
#endif
		
    SNESSetFromOptions(user[bi].snes);

    //    extern double imp_free_tol;
    SNESSetType(user[bi].snes,SNESNEWTONTR);			//SNESTR,SNESLS	: SNESLS is better for stiff PDEs such as the one including IB but slower

    SNESSetMaxLinearSolveFailures(user[bi].snes,10000);
    SNESSetMaxNonlinearStepFailures(user[bi].snes,10000);
    SNESKSPSetUseEW(user[bi].snes, PETSC_TRUE);
    SNESKSPSetParametersEW(user[bi].snes,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    //SNESSetTolerances(user[bi].snes,imp_free_tol,PETSC_DEFAULT,PETSC_DEFAULT,50,50000);
    SNESSetTolerances(user[bi].snes,PETSC_DEFAULT,imp_rtol,PETSC_DEFAULT,50,50000);
		
    SNESGetKSP(user[bi].snes, &ksp);
    KSPGetPC(ksp,&pc);
		
    if(!snes_created) {
      KSPSetType(ksp, KSPGMRES);
      //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);	//2009.09.22 poor performance
      //KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);	//2009.09.22
			
      //KSPFischerGuess itg;
      //KSPFischerGuessCreate(ksp,1,100,&itg);
      //KSPSetFischerGuess(ksp, itg);	//2009.09.22
			
      //KSPGMRESSetPreAllocateVectors(ksp);	--> crazy thing consumes memory
    }
#if defined(PETSC_HAVE_ADIC___)
    PCSetType(pc,PCBJACOBI);
#else
    PCSetType(pc,PCNONE);
#endif
    int maxits=1000;
    //double rtol=PETSC_DEFAULT, atol=imp_free_tol, dtol=PETSC_DEFAULT;
    double rtol=imp_rtol, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT;

    KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
		
    snes_created=1;
  }
	
  bi =0;
  //int it=0;
		
  Vec U;
  VecDuplicate(user[bi].Ucont, &U);
  //  extern int outflow_scale;
  
  InflowFlux(&(user[bi]));
  OutflowFlux(&(user[bi]));
  
  PetscPrintf(PETSC_COMM_WORLD, "FluxOutRK %le\n", user->FluxOutSum);

  //  outflow_scale=0;
  FormBCS(&(user[bi]));
		
  VecCopy(user[bi].Ucont, U);
		
  SNESSolve(user[bi].snes, PETSC_NULL, U);
		
  PetscPrintf(PETSC_COMM_WORLD, "\nMomentum eqs computed ...\n");
	
  double norm;Vec Norm;
  // SNESGetFunctionNorm(user[bi].snes, &norm);
 SNESGetFunction(user[bi].snes, &Norm,0,0);
  VecNorm(Norm, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "\nSNES residual norm=%.5e\n\n", norm);
		
  VecCopy(U, user[bi].Ucont);
		
  DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
  DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
		
/*   Vec Coor; */
/*   Cmpnts ***ucont, ***lucont, ***coor, ***icsi, ***lucat; */
/*   int i, j, k; */
	
/*   DAGetGhostedCoordinates(user[0].da, &Coor); */
/*   DMDAVecGetArray(user[0].fda, Coor, &coor); */
/*   DMDAVecGetArray(user[0].fda, user[0].Ucont, &ucont); */
/*   DMDAVecGetArray(user[0].fda, user[0].lUcont, &lucont); */
/*   DMDAVecGetArray(user[0].fda, user[0].lUcat, &lucat); */
/*   DMDAVecGetArray(user[0].fda, user[0].lICsi, &icsi); */
/*   for (k=zs; k<ze; k++) */
/*     for (j=ys; j<ye; j++) */
/*       for (i=xs; i<xe; i++) { */
		
/* 	if ( (user[0].bctype[0]==1 || user[0].bctype[0]==-1 || user[0].bctype[0]==10) && i==0) ucont[k][j][i].x = 0; */
/* 	if ( (user[0].bctype[1]==1 || user[0].bctype[1]==-1 || user[0].bctype[1]==10) && i==mx-2) ucont[k][j][i].x = 0; */
/* 	if ( (user[0].bctype[2]==1 || user[0].bctype[2]==-1 || user[0].bctype[2]==10) && j==0) ucont[k][j][i].y = 0; */
/* 	if ( (user[0].bctype[3]==1 || user[0].bctype[3]==-1 || user[0].bctype[3]==10 || user[0].bctype[3]==2) && j==my-2) ucont[k][j][i].y = 0; */
/* 	if ( (user[0].bctype[4]==1 || user[0].bctype[4]==-1 || user[0].bctype[4]==10) && k==0) ucont[k][j][i].z = 0; */
/* 	if ( (user[0].bctype[5]==1 || user[0].bctype[5]==-1 || user[0].bctype[5]==10) && k==mz-2) ucont[k][j][i].z = 0; */
			
/* 	if ( user[0].bctype[5]==4 && k==mz-2) ucont[k][j][i].z = lucont[k-1][j][i].z;		// 07.25.2009	// 09.26.2009 when k is very thin k-1 can be other processor's */
		
/* 	if (user->bctype[0]==11 && i==0 && (j!=0 && j!=my-1 && k!=0 && k!=mz-1) ) { */
/* 	  double zc = ( coor[k][j][i+1].z + coor[k-1][j][i+1].z + coor[k][j-1][i+1].z + coor[k-1][j-1][i+1].z ) * 0.25; */
/* 	  if( zc > 0 ) { */
/* 	    ucont[k][j][i].x = 0 * icsi[k][j][i].x + 0 * icsi[k][j][i].y + lucat[k][j][i].z * icsi[k][j][i].z; */
/* 	  } */
/* 	} */
/*       } */
/*   DMDAVecRestoreArray(user[0].fda, Coor, &coor); */
/*   DMDAVecRestoreArray(user[0].fda, user[0].Ucont, &ucont); */
/*   DMDAVecRestoreArray(user[0].fda, user[0].lUcont, &lucont); */
/*   DMDAVecRestoreArray(user[0].fda, user[0].lUcat, &lucat); */
/*   DMDAVecRestoreArray(user[0].fda, user[0].lICsi, &icsi); */
	
/*   outflow_scale=1; */
  InflowFlux(&(user[bi]));
  OutflowFlux(&(user[bi]));
  FormBCS(&(user[bi]));
  if (immersed && ti>0 ) 
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
    }
  
  if (block_number>1) Block_Interface_U(user);

  PetscTime(&te);
  cput=te-ts;
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "results/Converge_dU");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d(momentum) %.2e(s) %le\n", ti, cput, norm);
    fclose(f);
  }
  VecDestroy(&U);
	
  for (bi=0; bi<block_number; bi++) {
/*     DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); */
/*     DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); */
/*     VecWAXPY(user[bi].DUold, -1., user[bi].Ucont_o, user[bi].Ucont); */
    /*
		
    SNESDestroy(user[bi].snes);
    */
    VecDestroy(&user[bi].Rhs);
    MatDestroy(&user[bi].J);
  }
	
  double max_norm;
  VecMax(user->Ucat, &i, &max_norm);
  PetscPrintf(PETSC_COMM_WORLD, "\n*** Max Ucat = %e \n", max_norm);
      
  return(0);
}



PetscErrorCode  Implicit_SNES(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)  //Solver No.5
{
  PetscInt	bi,ibi;
  KSP ksp;
  PC pc;
  Vec U;
  int rank;
	
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (block_number>1) Block_Interface_U(user);
  for (bi=0; bi<block_number; bi++) {
     OutflowFlux(&(user[bi]));
  }
  //--------------start for each block
  for (bi=0; bi<block_number; bi++) {

    VecDuplicate(user[bi].Ucont, &user[bi].Rhs);
    SNESCreate(PETSC_COMM_WORLD,&user[bi].snes);
    SNESSetFunction(user[bi].snes,user[bi].Rhs,FormFunction_SNES,(void *)&user[bi]);

    if(implicit_type==1){
    //Start--------------------Finite difference jacobian--------------------------------
    DMCreateColoring(user[bi].fda,IS_COLORING_GLOBAL,&(user[bi].iscoloring));
    DMCreateMatrix(user[bi].fda,&user[bi].J);
    MatFDColoringCreate(user[bi].J,user[bi].iscoloring,&(user[bi].matfdcoloring));
    MatFDColoringSetFunction(user[bi].matfdcoloring,(PetscErrorCode (*)(void))FormFunction_SNES,(void *)&user[bi]);
    MatFDColoringSetFromOptions(user[bi].matfdcoloring);
    MatFDColoringSetUp(user[bi].J,user[bi].iscoloring,user[bi].matfdcoloring);
    ISColoringDestroy(&(user[bi].iscoloring));
    SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,FormJacobian_FD,&user[bi]);

    //End--------------------Finite difference jacobian--------------------------------
    } else if(implicit_type==2){
    //Start--------------------Matrix Free---------------------------------------------
      MatCreateSNESMF(user[bi].snes, &user[bi].J);
      PetscInt M,N;
      VecGetSize(user[bi].Rhs,&N);
      VecGetLocalSize(user[bi].Rhs, &M);
      MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N,3, PETSC_NULL,3, PETSC_NULL,&user[bi].PJ);

      SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].PJ,FormJacobian_MF,(void *)&user[bi]);

    //End----------------------Matrix Free---------------------------------------------
    }else if(implicit_type==3){  //Diagonal
       
      PetscInt M,N;
      VecGetSize(user[bi].Rhs,&N);
      VecGetLocalSize(user[bi].Rhs, &M);

      MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N,3, PETSC_NULL,3, PETSC_NULL,&user[bi].J);  //for diagonal 3 , 51 for jacobian
      
      SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,FormJacobian_Diagonal,(void *)&user[bi]);
    }else if(implicit_type==4){  //FUll
 
      PetscInt M,N;
      VecGetSize(user[bi].Rhs,&N);
      VecGetLocalSize(user[bi].Rhs, &M);
      MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N,102, PETSC_NULL,102, PETSC_NULL,&user[bi].J);  //for diagonal 3 , 51 for jacobian
      
      SNESSetJacobian(user[bi].snes,user[bi].J,user[bi].J,FormJacobian,(void *)&user[bi]);

    }
    SNESGetKSP(user[bi].snes, &ksp);
    //---------------NULL space removing from jacobian
    // MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp);
    // MatNullSpaceSetFunction(user[bi].nullsp,NSNullSpaceFunction, &user[bi]);
    // KSPSetNullSpace(ksp, user[bi].nullsp);
   //-------------- End Null space removing

 	
    SNESSetFromOptions(user[bi].snes);

   

    SNESSetMaxLinearSolveFailures(user[bi].snes,10000);
    SNESSetMaxNonlinearStepFailures(user[bi].snes,10000);
    SNESKSPSetUseEW(user[bi].snes, PETSC_TRUE);
    SNESKSPSetParametersEW(user[bi].snes,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
     	
    KSPGetPC(ksp,&pc);

    VecDuplicate(user[bi].Ucont, &U);
  InflowFlux(&(user[bi]));
  OutflowFlux(&(user[bi]));
  FormBCS(&(user[bi]));
		
  VecCopy(user[bi].Ucont, U);


   SNESSolve(user[bi].snes, PETSC_NULL, U);

  MatFDColoringDestroy(&(user[bi].matfdcoloring));
 		
  PetscPrintf(PETSC_COMM_WORLD, "\nMomentum eqs computed ...\n");

  PetscReal norm;
  //SNESGetFunctionNorm(user[bi].snes, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "\nSNES residual norm=%.5e\n\n", norm);
  SNESDestroy(&user[bi].snes);
  VecCopy(U, user[bi].Ucont);
		
  DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
  DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
		

  InflowFlux(&(user[bi]));
  OutflowFlux(&(user[bi]));
  FormBCS(&(user[bi]));
  CopyUContPBC(&(user[bi]));
  //ForceSolidBoundary(&(user[bi]));
 
  //  ForceSolidBoundary(&(user[bi]));
  if (immersed && ti>0 )
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
    }
   
  VecDestroy(&U);

  }  //for loop on through each block number


  //--------------end for each block
  //-------Velocity correction for each block
  if (block_number>1) Block_Interface_U(user);
  if (block_number>1) {
    Block_Interface_Correction(user);
     if(blank) {
       Block_Blank_Correction_adv(&user[0],0);
     }
  }

  //-------Destroying
 for (bi=0; bi<block_number; bi++) {
   //  MatNullSpaceDestroy(&user[bi].nullsp);
   VecDestroy(&user[bi].Rhs);
   //
   MatDestroy(&user[bi].J);
   MatDestroy(&user[bi].PJ);
  }
 
  return(0); 
  	
}



PetscErrorCode ModifyJacobian(Mat PJ, UserCtx *user)
{
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert,Unit=1.0;
  PetscInt	i, j, k,row,kk;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

 /*  if (xs==0) lxs = xs+1; */
/*   if (ys==0) lys = ys+1; */
/*   if (zs==0) lzs = zs+1; */

/*   if (xe==mx) lxe = xe-1; */
/*   if (ye==my) lye = ye-1; */
/*   if (ze==mz) lze = ze-1; */
  


  DMDAVecGetArray(user->da, user->lNvert, &nvert);

 
  if (zs == 0) {
    k = 0;  //physical ghost and boundary node
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	
	row = Gidx(i, j, k, user);
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);

	if (user->bctype[4]!=7) {
	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
      }
    }
  }

  if (ze == mz) {
    k = mz-1;  //physical ghost node
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	row = Gidx(i, j, k, user);
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);	
      }
    }
    k = mz-2; //Boundary node
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	row = Gidx(i, j, k, user);
	if (user->bctype[5]!=7) {
	  kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
      }
    }
  }

  if (ys == 0) {
    j = 0;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	row = Gidx(i, j, k, user);
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);

	if (user->bctype[2]!=7) {
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}

	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);	
      }
    }
  }

  if (ye == my) {
    j = my-1;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	row = Gidx(i, j, k, user);
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);	
      }
    }
    j = my-2;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	row = Gidx(i,j,k,user);
	if (user->bctype[3]!=7) {
       	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	
      }
    }
  }

  if (xs == 0) {
    i = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	row = Gidx(i, j, k, user);
	if (user->bctype[0]!=7) {
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);	
      }
    }
  }

  if (xe == mx) {
    i = mx-1;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	row = Gidx(i, j, k, user);
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);	
      }
    }
    i = mx-2;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	row = Gidx(i, j, k, user);
	if (user->bctype[1]!=7) {
	  kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {   //Immersed Boundary
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe-1; i++) {
   	if (0.5*(nvert[k][j][i]+nvert[k][j][i+1]) > 0.1) {
	 row = Gidx(i, j, k, user);
	kk=3*row;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	
	 
      }
    }
  }
 for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye-1; j++) {
      for (i=lxs; i<lxe; i++) {
	if (0.5*(nvert[k][j][i]+nvert[k][j+1][i]) > 0.1) {
	  row = Gidx(i, j, k, user);
	kk=3*row+1;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
 
      }
    }
  }
 for (k=lzs; k<lze-1; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (0.5*(nvert[k][j][i]+nvert[k+1][j][i]) > 0.1) {
        row = Gidx(i, j, k, user);
	kk=3*row+2;
	MatSetValues(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	 
      }
    }
  }


  DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
  
  


return 0;
}


 PetscErrorCode FormJacobian_FD(SNES snes,Vec Ucont,Mat J,Mat PJ,void *ctx)
{

 //---------------Form FD
  UserCtx *user = (UserCtx*)ctx;
  //*flag = SAME_NONZERO_PATTERN;
  DM da = user[user->_this].da;
  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);
  PetscInt bi=user->_this;

  MatFDColoring  color = user[user->_this].matfdcoloring;
  Vec            f;
  PetscErrorCode (*ff)(void),(*fd)(void);
  void           *fctx;
  
 
 SNESGetFunction(snes,&f,(PetscErrorCode (**)(SNES,Vec,Vec,void*))&ff,0);
 MatFDColoringGetFunction(color,&fd,&fctx);
 MatFDColoringApply(J,color,Ucont,user[bi].snes);
 

 ModifyJacobian(J,&user[bi]);

 MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
 MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);

 //---------------Form AJ
 PetscInt M,N;
 Mat AJ;
 VecGetSize(user[bi].Rhs,&N);
 VecGetLocalSize(user[bi].Rhs, &M);
 MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N,51, PETSC_NULL,51, PETSC_NULL,&AJ); 
 FormJacobian(snes,Ucont,AJ,AJ,&user[bi]);

 MatDestroy(&AJ);


 return 0;
}

 PetscErrorCode FormJacobian_MF(SNES snes,Vec Ucont,Mat J,Mat PJ,void *ctx)
{

  //--------------First Jacobian based on MatMFFDSetFunction------------------------------
  MatMFFDSetFunction(J,(PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction,snes);  
  MatMFFDSetBase(J,Ucont,NULL);

  //-----------Second Jacobian Based on Full-Diagonal as Preconditioner-------------------
  UserCtx *user = (UserCtx*)ctx;
 /*  *flag = SAME_NONZERO_PATTERN; */
  PetscInt bi=user->_this;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert,val[3],c1=1./12.,c2=8./12.;
  PetscInt	i,j,k,m,n,p,row[3],col,kk,dd;
  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;
  PetscReal U0jp,U0kp,U1ip,U1kp,U2ip,U2jp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj,A[5][5];

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else
    dtc=1./user->dt;

  Re=user->ren;
  
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
   
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

   
	  
	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	
	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);
	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);
	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;


	//--Matrix below is needed to implement d/dcsi[Ui dUj]
	A[0][0]=0.125*(aj[k][j][i]*ucont[k][j][i].y); A[0][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y);
	A[0][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].y); A[0][3]=-0.125*(aj[k][j-1][i+1]*ucont[k][j-1][i+1].y);

	A[1][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[1][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z);
	A[1][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].z); A[1][3]=-0.125*(aj[k-1][j][i+1]*ucont[k-1][j][i+1].z);

	A[2][0]=-0.125*(aj[k][j+1][i-1]*ucont[k][j+1][i-1].x); A[2][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x);
	A[2][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].x); A[2][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[3][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[3][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z);
	A[3][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].z); A[3][3]=-0.125*(aj[k-1][j+1][i]*ucont[k-1][j+1][i].z);

	A[4][0]=-0.125*(aj[k+1][j][i-1]*ucont[k+1][j][i-1].x); A[4][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x);
	A[4][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].x); A[4][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[5][0]=-0.125*(aj[k+1][j-1][i]*ucont[k+1][j-1][i].y); A[5][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y);
	A[5][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].y); A[5][3]=0.125*(aj[k][j][i]*ucont[k][j][i].y);

	  kk = Gidx(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=Gidx(i, j, k, user);  //(i,j,k)
	  //--calc needed
	  m=i,n=j,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--

	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re+A[0][0]+A[0][1]+A[0][2]+A[0][3]+A[1][0]+A[1][1]+A[1][2]+A[1][3]);
	  val[1]=-(0.5*AJip*U1ip);
	  val[2]=-(0.5*AJip*U2ip);
	  MatSetValues(PJ,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(0.5*AJjp*U0jp);
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re+A[2][0]+A[2][1]+A[2][2]+A[2][3]+A[3][0]+A[3][1]+A[3][2]+A[3][3]);
	  val[2]=-(0.5*AJjp*U2jp);;
	  MatSetValues(PJ,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(0.5*AJkp*U0kp);;
	  val[1]=-(0.5*AJkp*U1kp);;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re+A[4][0]+A[4][1]+A[4][2]+A[4][3]+A[5][0]+A[5][1]+A[5][2]+A[5][3]);
	  MatSetValues(PJ,3,row,1,&col,val,INSERT_VALUES);
	 
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  

  ModifyJacobian(PJ,&user[bi]);
  
  MatAssemblyBegin(PJ,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(PJ,MAT_FINAL_ASSEMBLY);

  return 0;
}



PetscErrorCode FormJacobian(SNES snes,Vec Ucont,Mat J,Mat PJ,void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;
  // *flag = SAME_NONZERO_PATTERN;
  PetscInt bi=user->_this;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     val[3],c1=1./12.,c2=8./12.;
  PetscInt	i,j,k,m,n,p,row[3],col,kk,dd;
  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;
  PetscReal U0jp,U0kp,U1ip,U1kp,U2ip,U2jp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj,A[6][6]; 

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else 
    dtc=1./user->dt;

  Re=user->ren;
  
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
   
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

  
	 
	  
	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	
	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);
	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);
	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;


	//--Matrix below is needed to implement d/dcsi[Ui dUj]
	A[0][0]=0.125*(aj[k][j][i]*ucont[k][j][i].y); A[0][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y); 
	A[0][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].y); A[0][3]=-0.125*(aj[k][j-1][i+1]*ucont[k][j-1][i+1].y);

	A[1][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[1][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z); 
	A[1][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].z); A[1][3]=-0.125*(aj[k-1][j][i+1]*ucont[k-1][j][i+1].z);

	A[2][0]=-0.125*(aj[k][j+1][i-1]*ucont[k][j+1][i-1].x); A[2][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x); 
	A[2][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].x); A[2][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[3][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[3][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z);
	A[3][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].z); A[3][3]=-0.125*(aj[k-1][j+1][i]*ucont[k-1][j+1][i].z);

	A[4][0]=-0.125*(aj[k+1][j][i-1]*ucont[k+1][j][i-1].x); A[4][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x);
	A[4][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].x); A[4][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[5][0]=-0.125*(aj[k+1][j-1][i]*ucont[k+1][j-1][i].y); A[5][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y); 
	A[5][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].y); A[5][3]=0.125*(aj[k][j][i]*ucont[k][j][i].y);

	  kk = Gidx(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=Gidx(i, j, k, user);  //(i,j,k)
	  //--calc needed
	  m=i,n=j,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--

	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re+A[0][0]+A[0][1]+A[0][2]+A[0][3]+A[1][0]+A[1][1]+A[1][2]+A[1][3]);
	  val[1]=-(0.5*AJip*U1ip);
	  val[2]=-(0.5*AJip*U2ip);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(0.5*AJjp*U0jp);
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re+A[2][0]+A[2][1]+A[2][2]+A[2][3]+A[3][0]+A[3][1]+A[3][2]+A[3][3]);
	  val[2]=-(0.5*AJjp*U2jp);;	  
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(0.5*AJkp*U0kp);;
	  val[1]=-(0.5*AJkp*U1kp);;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re+A[4][0]+A[4][1]+A[4][2]+A[4][3]+A[5][0]+A[5][1]+A[5][2]+A[5][3]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //----------------------------------------------
	  dd=Gidx(i+1, j, k, user);  //(i+1,j,k)
	  //--calc needed
	  m=i+1,n=j,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g11ip/Re)+(c2*ucont[p][n][m].x)+A[0][2]+A[0][3]+A[1][2]+A[1][3]);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(0.5*AJjp*U0jp);
	  val[1]=-((-0.5*AJjp*AJjp*g11jp/Re)+A[2][2]+A[2][3]);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(0.5*AJkp*U0kp);
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g11kp/Re)+A[4][2]+A[4][3]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i-1, j, k, user);  //(i-1,j,k)
	  //--calc needed
	  m=i-1,n=j,p=k;
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  //--
	  col=3*dd;
	  val[0]=-((-0.5*AJip*AJip*g11ip/Re)+(-c2*ucont[p][n][m].x)+A[0][0]+A[0][1]+A[1][0]+A[1][1]);
	  val[1]=-(-0.5*AJip*U1ip);
	  val[2]=-(-0.5*AJip*U2ip);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g11jp/Re)+A[2][0]+A[2][1]);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g11kp/Re)+A[4][0]+A[4][1]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------	
	  dd=Gidx(i, j+1, k, user);  //(i,j+1,k)
	  //--calc needed
	  m=i,n=j+1,p=k;
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g22ip/Re)+A[0][0]+A[0][2]);
	  val[1]=-(0.5*AJip*U1ip);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g22jp/Re)+(c2*ucont[p][n][m].y)+A[2][2]+A[2][0]+A[3][2]+A[3][3]);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=-(0.5*AJkp*U1kp);
	  val[2]=-((-0.5*AJkp*AJkp*g22kp/Re)+A[5][2]+A[5][3]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j-1, k, user);  //(i,j-1,k)
	  //--calc needed
	  m=i,n=j-1,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g22ip/Re)+A[0][1]+A[0][3]);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(-0.5*AJjp*U0jp);
	  val[1]=-((-0.5*AJjp*AJjp*g22jp/Re)+(-c2*ucont[p][n][m].y)+A[2][3]+A[2][1]+A[3][0]+A[3][1]);
	  val[2]=-(-0.5*AJjp*U2jp);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g22kp/Re)+A[5][0]+A[5][1]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j, k+1, user);  //(i,j,k+1)
	  //--calc needed
	  m=i,n=j,p=k+1;
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--
	  col=3*dd;
	  val[0]=-((-0.5*AJip*AJip*g33ip/Re+A[1][0]+A[1][2]));
	  val[1]=0.0;
	  val[2]=-(0.5*AJip*U2ip);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g33jp/Re)+A[3][0]+A[3][2]);
	  val[2]=-(0.5*AJjp*U2jp);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g33kp/Re)+(c2*ucont[p][n][m].z)+A[4][0]+A[4][2]+A[5][0]+A[5][2]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j, k-1, user);  //(i,j,k-1)
	  //--calc needed
	  m=i,n=j,p=k-1;
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g33ip/Re+A[1][1]+A[1][3]));
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g33jp/Re)+A[3][1]+A[3][3]);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(-0.5*AJkp*U0kp);
	  val[1]=-(-0.5*AJkp*U1kp);
	  val[2]=-((-0.5*AJkp*AJkp*g33kp/Re)+(-c2*ucont[p][n][m].z)+A[4][1]+A[4][3]+A[5][1]+A[5][3]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	 
	  //---------------------------------------------
	  dd=Gidx(i-1, j+1, k, user);  //(i-1,j+1,k)
	  //--calc needed
	  m=i-1,n=j+1,p=k;
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  //--
	  col=3*dd;
	  val[0]=-A[0][0];
	  val[1]=-(-0.5*AJip*U1ip);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[2][0];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i-1, j, k+1, user);  //(i-1,j,k+1)
	  //--calc needed
	  m=i-1,n=j,p=k+1;
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  //--
	  col=3*dd;
	  val[0]=-A[1][0];
	  val[1]=0.0;
	  val[2]=-(-0.5*AJip*U2ip);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[4][0];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i+1, j-1, k, user);  //(i+1,j-1,k)
	  //--calc needed
	  m=i+1,n=j-1,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  //--
	  col=3*dd;
	  val[0]=-A[0][3];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=-(-0.5*AJjp*U0jp);
	  val[1]=-A[2][3];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j-1, k+1, user);  //(i,j-1,k+1)
	  //--calc needed
	  m=i,n=j-1,p=k+1;
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][0];
	  val[2]=-(-0.5*AJjp*U2jp);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[5][0];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i+1, j, k-1, user);  //(i+1,j,k-1)
	  //--calc needed
	  m=i+1,n=j,p=k-1;
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  //--
	  col=3*dd;
	  val[0]=-A[1][3];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=-(-0.5*AJkp*U0kp);
	  val[1]=0.0;
	  val[2]=-A[4][3];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j+1, k-1, user);  //(i,j+1,k-1)
	  //--calc needed
	  m=i,n=j+1,p=k-1;
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  //--
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][3];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=-(-0.5*AJkp*U1kp);
	  val[2]=-A[5][3];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i-1, j-1, k, user);  //(i-1,j-1,k)
	 
	  col=3*dd;
	  val[0]=-A[0][1];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[2][1];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i+1, j+1, k, user);  //(i+1,j+1,k)
	 
	  col=3*dd;
	  val[0]=-A[0][2];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[2][2];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i-1, j, k-1, user);  //(i-1,j,k-1)
	 
	  col=3*dd;
	  val[0]=-A[1][1];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[4][1];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i+1, j, k+1, user);  //(i+1,j,k+1)
	 
	  col=3*dd;
	  val[0]=-A[1][2];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[4][2];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j-1, k-1, user);  //(i,j-1,k-1)
	 
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][1];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[5][1];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  dd=Gidx(i, j+1, k+1, user);  //(i,j+1,k+1)
	 
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][2];
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[5][2];
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  //---------------------------------------------
	  if(i != 1){
	    m=i-2,n=j,p=k;
	  dd=Gidx(i-2, j, k, user);  //(i-2,j,k)
	  col=3*dd;
	  val[0]=-(c1*ucont[p][n][m].x);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  }
	  //---------------------------------------------
	  if(j != 1){
	    m=i,n=j-2,p=k;
	  dd=Gidx(i, j-2, k, user);  //(i,j-2,k)
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(c1*ucont[p][n][m].y);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  }
	  //---------------------------------------------
	  if(k != 1){
	    m=i,n=j,p=k-2;
	  dd=Gidx(i, j, k-2, user);  //(i,j,k-2)
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(c1*ucont[p][n][m].z);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  }
	  //---------------------------------------------
	  if(i != mx-2){
	    m=i+2,n=j,p=k;
	  dd=Gidx(i+2, j, k, user);  //(i+2,j,k)
	  col=3*dd;
	  val[0]=-(-c1*ucont[p][n][m].x);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  }
	  //---------------------------------------------
	  if(j != my-2){
	    m=i,n=j+2,p=k;
	  dd=Gidx(i, j+2, k, user);  //(i,j+2,k)
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(-c1*ucont[p][n][m].y);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  }
	  //---------------------------------------------
	  if(k != mz-2){
	    m=i,n=j,p=k+2;
	  dd=Gidx(i, j, k+2, user);  //(i,j,k+2)
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(-c1*ucont[p][n][m].z);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	  }
	  //---------------------------------------------
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

 

  ModifyJacobian(J,&user[bi]);
  
  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);

  //
  MatView(J,PETSC_VIEWER_DRAW_WORLD); 
  //

 return 0;
}

PetscErrorCode FormJacobian_Diagonal(SNES snes,Vec Ucont,Mat J,Mat PJ,void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;
  // *flag = SAME_NONZERO_PATTERN;
  PetscInt bi=user->_this;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     val[3];
  PetscInt	i,j,k,row[3],col,kk,dd;
  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  //Periodic Boundary
  if (user->bctype[0]!=7) {  if (xs==0) lxs = xs+1;}
  if (user->bctype[2]!=7) {  if (ys==0) lys = ys+1;}
  if (user->bctype[4]!=7) {  if (zs==0) lzs = zs+1;}
 /*  if (xs==0) lxs = xs+1; */
/*   if (ys==0) lys = ys+1; */
/*   if (zs==0) lzs = zs+1; */

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
 
  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj; 

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else 
    dtc=1./user->dt;

  Re=user->ren;

  /* update metrics for periodic boundary conditions */
  Update_Metrics_PBC(&(user[bi]));

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
   
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

   	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	

	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);

	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);

	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;
	



	  kk = Gidx(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=Gidx(i, j, k, user);  //(i,j,k)
	 

	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
		   
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
			    
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re); 
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

 

  ModifyJacobian(J,&user[bi]);

  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);



 return 0;
}


PetscErrorCode FormJacobian_Diagonal_MF(Mat J,void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;
   PetscInt bi=user->_this;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     val[3];
  PetscInt	i,j,k,row[3],col,kk,dd;
  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj; 

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else 
    dtc=1./user->dt;

  Re=user->ren;
  
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
   
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

   	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	
	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);
	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);
	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;





	  kk = Gidx(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=Gidx(i, j, k, user);  //(i,j,k)
	 

	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
		   
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re);
	  val[2]=0.0;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
			    
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);

	  
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

 

  ModifyJacobian(J,&user[bi]);

  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);



 return 0;
}

PetscErrorCode FormJacobian_Full_Diagonal(SNES snes,Vec Ucont,Mat J,Mat PJ,void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;
  //*flag = SAME_NONZERO_PATTERN;
  PetscInt bi=user->_this;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     val[3];
  PetscInt	i,j,k,m,n,p,row[3],col,kk,dd;
  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;
  PetscReal U0jp,U0kp,U1ip,U1kp,U2ip,U2jp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj,A[6][6];

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else
    dtc=1./user->dt;

  Re=user->ren;
  
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
   
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

   
	  
	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	
	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);
	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);
	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;


	//--Matrix below is needed to implement d/dcsi[Ui dUj]
	A[0][0]=0.125*(aj[k][j][i]*ucont[k][j][i].y); A[0][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y);
	A[0][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].y); A[0][3]=-0.125*(aj[k][j-1][i+1]*ucont[k][j-1][i+1].y);

	A[1][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[1][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z);
	A[1][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].z); A[1][3]=-0.125*(aj[k-1][j][i+1]*ucont[k-1][j][i+1].z);

	A[2][0]=-0.125*(aj[k][j+1][i-1]*ucont[k][j+1][i-1].x); A[2][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x);
	A[2][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].x); A[2][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[3][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[3][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z);
	A[3][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].z); A[3][3]=-0.125*(aj[k-1][j+1][i]*ucont[k-1][j+1][i].z);

	A[4][0]=-0.125*(aj[k+1][j][i-1]*ucont[k+1][j][i-1].x); A[4][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x);
	A[4][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].x); A[4][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[5][0]=-0.125*(aj[k+1][j-1][i]*ucont[k+1][j-1][i].y); A[5][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y);
	A[5][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].y); A[5][3]=0.125*(aj[k][j][i]*ucont[k][j][i].y);

	  kk = Gidx(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=Gidx(i, j, k, user);  //(i,j,k)
	  //--calc needed
	  m=i,n=j,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--

	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re+A[0][0]+A[0][1]+A[0][2]+A[0][3]+A[1][0]+A[1][1]+A[1][2]+A[1][3]);
	  val[1]=-(0.5*AJip*U1ip);
	  val[2]=-(0.5*AJip*U2ip);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(0.5*AJjp*U0jp);
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re+A[2][0]+A[2][1]+A[2][2]+A[2][3]+A[3][0]+A[3][1]+A[3][2]+A[3][3]);
	  val[2]=-(0.5*AJjp*U2jp);;
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(0.5*AJkp*U0kp);;
	  val[1]=-(0.5*AJkp*U1kp);;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re+A[4][0]+A[4][1]+A[4][2]+A[4][3]+A[5][0]+A[5][1]+A[5][2]+A[5][3]);
	  MatSetValues(J,3,row,1,&col,val,INSERT_VALUES);
	 
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

  

  ModifyJacobian(J,&user[bi]);
  
  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);


 return 0;
}


PetscErrorCode SampleShellPCCreate(SampleShellPC **shell)
{
  SampleShellPC  *newctx;

     PetscNew(&newctx);
     newctx->diag = 0;
     *shell       = newctx;
   return 0;
}

PetscErrorCode SampleShellPCSetUp(PC pc,Mat pmat,Vec x)
{
   SampleShellPC  *shell;
   Vec            diag;
   
   PCShellGetContext(pc,(void**)&shell);
   VecDuplicate(x,&diag);
   MatGetDiagonal(pmat,diag);
   VecReciprocal(diag);
   //VecDuplicate(x,&(shell->diag));
  shell->diag = diag;
  return 0;
 }

PetscErrorCode SampleShellPCApply(PC pc,Vec x,Vec y)
 {
   SampleShellPC  *shell;
   
   PCShellGetContext(pc,(void**)&shell);
   VecPointwiseMult(y,x,shell->diag);
   
   return 0;
 }

PetscErrorCode SampleShellPCDestroy(PC pc)
{
  SampleShellPC  *shell;

  PCShellGetContext(pc,(void**)&shell);
  VecDestroy(&shell->diag);
  PetscFree(shell);

  return 0;
}


// Added by Hafez 06/14  to use snes for multiblock, consider that maximum block(overset grid)numbers is assumed 5.
PetscErrorCode  Implicit_SNES_Packer(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi)  //Solver No.7
{
   
  PetscInt  bi,ibi,level,rank;
  KSP ksp;
  PC pc;
  Vec U,R;
  PetscReal ts=0.0,te=0.0;

  PetscTime(&ts);   
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  //Recover user context
  level = usermg->mglevels-1;
  UserCtx *user = usermg->mgctx[level].user;

  //Create packer 
  DMCompositeCreate(PETSC_COMM_WORLD,&(usermg->packer));
  for (bi=0; bi<block_number; bi++) {
    DMCompositeAddDM(usermg->packer,(DM)user[bi].fda);
  }
 
  
  if (block_number>1) {Block_Interface_U(user);}
  for (bi=0; bi<block_number; bi++) {
    OutflowFlux(&(user[bi]));
    ComputelNface(&(user[bi]));
  }
  

  SNESCreate(PETSC_COMM_WORLD,&(usermg->snespacker));

  //Set Jacobian
  PetscInt N=0,Nl=0;
  Mat JAC;
  for (bi=0; bi<block_number; bi++) {
    
    VecGetSize(user[bi].Ucont,&BSize[bi]);
    N+=BSize[bi];
    VecGetLocalSize(user[bi].Ucont,&BlSize[bi]);
    Nl+=BlSize[bi];
  }

  //--Create Matrix with MatCreatAIJ in Global form and map it to local form(we need to get submatrix from local form)
  ISLocalToGlobalMapping *ltogs,ltog;
  DMCompositeGetISLocalToGlobalMappings(usermg->packer,&ltogs);
  ISLocalToGlobalMappingConcatenate(PETSC_COMM_WORLD,block_number,ltogs,&ltog); 
  if(implicit_type==3){
    MatCreateAIJ(PETSC_COMM_WORLD, Nl, Nl, N, N,3, PETSC_NULL,3, PETSC_NULL,&JAC);  //for diagonal 3 , 51 for jacobian
  }else if(implicit_type==4){
    MatCreateAIJ(PETSC_COMM_WORLD, Nl, Nl, N, N,102, PETSC_NULL,102, PETSC_NULL,&JAC);  //for diagonal 3 , 51 for jacobian
  }
  //  DMCreateMatrix(usermg->packer,MATAIJ,&JAC);

  MatSetLocalToGlobalMapping(JAC,ltog,ltog);
  for (bi=0; bi<block_number; bi++) {
    ISLocalToGlobalMappingDestroy(&ltogs[bi]);
  }
  ISLocalToGlobalMappingDestroy(&ltog);

  PetscFree(ltogs);
  //--End

 
  SNESSetJacobian(usermg->snespacker,JAC,JAC,FormJacobian_Diagonal_Packer_All,(void *)usermg);

   //Creat U and R (JAC*U=R), user[bi].lUcontb is intermediate variable to transfer U to user[bi].Ucont and vise versa 
   DMCreateGlobalVector(usermg->packer,&U);
   VecDuplicate(U,&R);
 
  
   SNESSetInitialGuess_Packer(U,usermg);
   SNESSetFunction(usermg->snespacker,R,FormFunction_SNES_Packer,(void *)usermg);
  
   SNESGetKSP(usermg->snespacker, &ksp);
   SNESSetFromOptions(usermg->snespacker);
   KSPGetPC(ksp,&pc);
  
   SNESSolve(usermg->snespacker, PETSC_NULL,U);
    SNESFormFinal_Packer(U,usermg);
     
  //Boundary and Immersed correction
  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    OutflowFlux(&(user[bi]));
    FormBCS(&(user[bi]));
   
   if (immersed && ti>0 )
     for (ibi=0;ibi<NumberOfBodies;ibi++) {
       ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
     }

  }//Correction
  
  
  //-------Velocity correction for each block
  if (block_number>1) {Block_Interface_U(user);}
  if (block_number>1) {
    Block_Interface_Correction(user);
    if(blank) {
      Block_Blank_Correction_adv(&user[0],0);
    }
  }
 for (bi=0; bi<block_number; bi++) {
   CopyUContPBC(&(user[bi]));
 }
 
  //-------Destroying
  DMDestroy(&(usermg->packer));
  SNESDestroy(&(usermg->snespacker)); 
  MatDestroy(&JAC);
  VecDestroy(&U);
  VecDestroy(&R); 

  PetscTime(&te);
  cput+=te-ts;
  PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------------------------------------------\n");  
  PetscPrintf(PETSC_COMM_WORLD,"TotalNSCputime: %le , CurrentIterationNSCputime: %le\n",cput,te-ts); 
  PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------------------------------------------\n");   
  
  return(0);

}
PetscErrorCode SNESSetInitialGuess_Packer(Vec U,UserMG *usermg)
{
  
  PetscInt bi;
  Vec          Ub[block_number];
  //Recover user context
  PetscInt level = usermg->mglevels-1;
  UserCtx *user = usermg->mgctx[level].user;
 

   //Copy Ucont from user  to U
  DMCompositeGetAccess(usermg->packer,U,&Ub[0],&Ub[1],&Ub[2]);
  for (bi=0; bi<block_number; bi++) {
    VecCopy(user[bi].Ucont, Ub[bi]);
  }
  DMCompositeRestoreAccess(usermg->packer,U,&Ub[0],&Ub[1],&Ub[2]);
 

  return 0;
 }

PetscErrorCode FormFunction_SNES_Packer(SNES snes, Vec U, Vec R, void *usermgv)
{
  PetscInt ibi,bi;
  Vec          Ub[block_number], Rb[block_number];
  UserMG *usermg = (UserMG *)usermgv;
//Recover user context
  PetscInt level = usermg->mglevels-1;
  UserCtx *user = usermg->mgctx[level].user;
 
  //Copy Solution from ksp(U) to user
  DMCompositeGetAccess(usermg->packer, U,&Ub[0],&Ub[1],&Ub[2]);
  for (bi=0; bi<block_number; bi++) {
   VecCopy(Ub[bi], user[bi].Ucont);
 }
 for (bi=0; bi<block_number; bi++) {
   DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
   DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
 }
 DMCompositeRestoreAccess(usermg->packer, U,&Ub[0],&Ub[1],&Ub[2]);

 //----Compute RHS as usual
  if (block_number>1) {
    Block_Interface_U(user);
  }

  for (bi=0; bi<block_number; bi++) {
    InflowFlux(&(user[bi]));
    FormBCS(&(user[bi]));
    if (immersed)
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_advanced(&user[bi], &(user[bi].ibm[ibi]), ibi,1);
    }
    
    VecDuplicate(user[bi].Ucont, &(user[bi].Rhs));
    VecSet(user[bi].Rhs,0.);
    VecDuplicate(user[bi].Ucont, &(user[bi].dUcont));
    VecSet(user[bi].dUcont,0.);

    FormFunction1(&user[bi],user[bi].Rhs);

   
    if (COEF_TIME_ACCURACY>1.1 && ti!=tistart && ti!=1){
      VecCopy(user[bi].Ucont, user[bi].dUcont);
      VecScale(user[bi].dUcont, COEF_TIME_ACCURACY);
      VecAXPY(user[bi].dUcont, -2.,user[bi].Ucont_o );
      VecAXPY(user[bi].dUcont, 0.5, user[bi].Ucont_rm1);
    } else {
      VecWAXPY(user[bi].dUcont, -1., user[bi].Ucont_o, user[bi].Ucont);
    }
    
    // Add -du/dt to right hand side
    VecAXPY(user[bi].Rhs, -1./user->dt, user[bi].dUcont);

    SetRHSBc(&user[bi],user[bi].Rhs);
  } //through each block


  //Copy RHS from FormFunction(user[bi].Rhs) to R(ksp Right hand side)
  VecSet(R,0.);
  DMCompositeGetAccess(usermg->packer, R,&Rb[0],&Rb[1],&Rb[2]);
  for (bi=0; bi<block_number; bi++) {
    VecCopy(user[bi].Rhs, Rb[bi]);
  }
  DMCompositeRestoreAccess(usermg->packer, R,&Rb[0],&Rb[1],&Rb[2]);


 //Copy RHS from FormFunction(user[bi].Rhs) to R(ksp Right hand side)

/*   DMCompositeGetAccess(usermg->packer, U,&Ub[0],&Ub[1]); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(user[bi].Ucont, Ub[bi]); */
/*     } */
/*   DMCompositeRestoreAccess(usermg->packer, U,&Ub[0],&Ub[1]); */
  
  
  // Destroy
  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].Rhs);
    VecDestroy(&user[bi].dUcont);
  }


  return 0;
 }
PetscErrorCode  SNESFormFinal_Packer(Vec U,UserMG *usermg)
{
  
  PetscInt bi;
  Vec          Ub[block_number];
  //Recover user context
  PetscInt level = usermg->mglevels-1;
  UserCtx *user = usermg->mgctx[level].user;
  
  //Copy Solution from ksp(U) to user
  DMCompositeGetAccess(usermg->packer, U,&Ub[0],&Ub[1],&Ub[2]);
  for (bi=0; bi<block_number; bi++) {
    VecCopy(Ub[bi], user[bi].Ucont);
  }
  for (bi=0; bi<block_number; bi++) {
    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
  }
  DMCompositeRestoreAccess(usermg->packer, U,&Ub[0],&Ub[1],&Ub[2]);
  
  
  
  return 0;
}



PetscErrorCode FormJacobian_Diagonal_Packer_All(SNES snes,Vec Ucont,Mat J,Mat PJ,void *usermgv)
{

  
  UserMG *usermg = (UserMG *)usermgv;
//Recover user context
  PetscInt level = usermg->mglevels-1;
  UserCtx *user = usermg->mgctx[level].user;

  PetscInt bi;
  IS  *is;
  Mat JAC1,JAC2,JAC3;
 

  DMCompositeGetLocalISs(usermg->packer,&is);
 


  MatGetLocalSubMatrix(J,is[0],is[0],&JAC1);
  MatGetLocalSubMatrix(J,is[1],is[1],&JAC2);
   MatGetLocalSubMatrix(J,is[2],is[2],&JAC3);

 if(implicit_type==3){
  FormJacobian_Diagonal_Packer(JAC1,&user[0]);
  FormJacobian_Diagonal_Packer(JAC2,&user[1]);
  FormJacobian_Diagonal_Packer(JAC3,&user[2]);

 }else if(implicit_type==4){
  FormJacobian_Packer(JAC1,&user[0]);
  FormJacobian_Packer(JAC2,&user[1]);
  FormJacobian_Packer(JAC3,&user[2]);
 }
  MatRestoreLocalSubMatrix(J,is[0],is[0],&JAC1);
  MatRestoreLocalSubMatrix(J,is[1],is[1],&JAC2);
  MatRestoreLocalSubMatrix(J,is[2],is[2],&JAC3);

  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);


  for (bi=0; bi<block_number; bi++) {
    ISDestroy(&is[bi]);
  }
  PetscFree(is);


 return 0;
}
PetscErrorCode FormJacobian_Diagonal_Packer(Mat J,void *ctx)
{
  

  UserCtx *user = (UserCtx*)ctx;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;
  
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal     val[3];
  PetscInt	i,j,k,row[3],col,kk,dd;
  
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  //Periodic Boundary
  if (user->bctype[0]!=7) {  if (xs==0) lxs = xs+1;}
  if (user->bctype[2]!=7) {  if (ys==0) lys = ys+1;}
  if (user->bctype[4]!=7) {  if (zs==0) lzs = zs+1;}
  /*  if (xs==0) lxs = xs+1; */
  /*   if (ys==0) lys = ys+1; */
  /*   if (zs==0) lzs = zs+1; */
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
 
  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj; 

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else 
    dtc=1./user->dt;

  Re=user->ren;

  /* update metrics for periodic boundary conditions */
    Update_Metrics_PBC(user);

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
   
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

   	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	

	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);
	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);
	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;




	  kk =lidxLocal(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=lidxLocal(i, j, k, user);  //(i,j,k)

	 
	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
		   
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
			    
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  // CheckCol_packer(J,user,i,j,k);

	  if (i==mx-2 && user->bctype[1]!=7) {
	    CheckRowDiagonal_packer(J,user,i,j,k,0);
	  }
	  if (j==my-2 && user->bctype[3]!=7) {
	    CheckRowDiagonal_packer(J,user,i,j,k,1);
	  }
	  if (k==mz-2 && user->bctype[5]!=7) {
	    CheckRowDiagonal_packer(J,user,i,j,k,2);
	  }




	  
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

 

  ModifyJacobian_packer(J,user);

  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);


 return 0;
}
PetscErrorCode FormJacobian_Packer(Mat J,void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;
  DM fda = user->fda;
  DM da = user->da;
  DMDALocalInfo	info = user->info;
  
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert,val[3],c1=1./12.,c2=8./12.;
  PetscInt	i,j,k,m,n,p,row[3],col,kk,dd;
  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  //Compute the variables which is needed(at centers)
  PetscReal dtc,Re,AJip,AJjp,AJkp,g11ip,g22ip,g33ip,g11jp,g22jp,g33jp,g11kp,g22kp,g33kp;
  PetscReal U0jp,U0kp,U1ip,U1kp,U2ip,U2jp;

  Cmpnts	***csi, ***eta, ***zet,***ucont;
  PetscReal ***aj,A[6][6]; 

  if (ti!=tistart && ti!=1)
    dtc=1.5/user->dt;
  else 
    dtc=1./user->dt;

  Re=user->ren;
  
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);
 
  DMDAVecGetArray(da, user->lAj, &aj);
   DMDAVecGetArray(fda, user->lUcont, &ucont);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

  	  
	
	AJip=0.5*(aj[k][j][i]+aj[k][j][i+1]);
	AJjp=0.5*(aj[k][j][i]+aj[k][j+1][i]);
	AJkp=0.5*(aj[k+1][j][i]+aj[k][j][i]);
	
	//----on i+1/2
	g11ip=csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z;
	g22ip=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k][j][i+1].x*eta[k][j][i+1].x+eta[k][j][i+1].y*eta[k][j][i+1].y+eta[k][j][i+1].z*eta[k][j][i+1].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k][j-1][i+1].x*eta[k][j-1][i+1].x+eta[k][j-1][i+1].y*eta[k][j-1][i+1].y+eta[k][j-1][i+1].z*eta[k][j-1][i+1].z);
	g33ip=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j][i+1].x*zet[k][j][i+1].x+zet[k][j][i+1].y*zet[k][j][i+1].y+zet[k][j][i+1].z*zet[k][j][i+1].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j][i+1].x*zet[k-1][j][i+1].x+zet[k-1][j][i+1].y*zet[k-1][j][i+1].y+zet[k-1][j][i+1].z*zet[k-1][j][i+1].z);
	//----on j+1/2
	g11jp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k][j+1][i].x*csi[k][j+1][i].x+csi[k][j+1][i].y*csi[k][j+1][i].y+csi[k][j+1][i].z*csi[k][j+1][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k][j+1][i-1].x*csi[k][j+1][i-1].x+csi[k][j+1][i-1].y*csi[k][j+1][i-1].y+csi[k][j+1][i-1].z*csi[k][j+1][i-1].z);
	
	g22jp=eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z;

	g33jp=0.25*(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z+
		    zet[k][j+1][i].x*zet[k][j+1][i].x+zet[k][j+1][i].y*zet[k][j+1][i].y+zet[k][j+1][i].z*zet[k][j+1][i].z+
		    zet[k-1][j][i].x*zet[k-1][j][i].x+zet[k-1][j][i].y*zet[k-1][j][i].y+zet[k-1][j][i].z*zet[k-1][j][i].z+
		    zet[k-1][j+1][i].x*zet[k-1][j+1][i].x+zet[k-1][j+1][i].y*zet[k-1][j+1][i].y+zet[k-1][j+1][i].z*zet[k-1][j+1][i].z);
	//----on k+1/2
	g11kp=0.25*(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z+
		    csi[k+1][j][i].x*csi[k+1][j][i].x+csi[k+1][j][i].y*csi[k+1][j][i].y+csi[k+1][j][i].z*csi[k+1][j][i].z+
		    csi[k][j][i-1].x*csi[k][j][i-1].x+csi[k][j][i-1].y*csi[k][j][i-1].y+csi[k][j][i-1].z*csi[k][j][i-1].z+
		    csi[k+1][j][i-1].x*csi[k+1][j][i-1].x+csi[k+1][j][i-1].y*csi[k+1][j][i-1].y+csi[k+1][j][i-1].z*csi[k+1][j][i-1].z);

	g22kp=0.25*(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z+
		    eta[k+1][j][i].x*eta[k+1][j][i].x+eta[k+1][j][i].y*eta[k+1][j][i].y+eta[k+1][j][i].z*eta[k+1][j][i].z+
		    eta[k][j-1][i].x*eta[k][j-1][i].x+eta[k][j-1][i].y*eta[k][j-1][i].y+eta[k][j-1][i].z*eta[k][j-1][i].z+
		    eta[k+1][j-1][i].x*eta[k+1][j-1][i].x+eta[k+1][j-1][i].y*eta[k+1][j-1][i].y+eta[k+1][j-1][i].z*eta[k+1][j-1][i].z);

	g33kp=zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z;


	//--Matrix below is needed to implement d/dcsi[Ui dUj]
	A[0][0]=0.125*(aj[k][j][i]*ucont[k][j][i].y); A[0][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y); 
	A[0][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].y); A[0][3]=-0.125*(aj[k][j-1][i+1]*ucont[k][j-1][i+1].y);

	A[1][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[1][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z); 
	A[1][2]=0.125*(aj[k][j][i+1]*ucont[k][j][i+1].z); A[1][3]=-0.125*(aj[k-1][j][i+1]*ucont[k-1][j][i+1].z);

	A[2][0]=-0.125*(aj[k][j+1][i-1]*ucont[k][j+1][i-1].x); A[2][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x); 
	A[2][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].x); A[2][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[3][0]=0.125*(aj[k][j][i]*ucont[k][j][i].z); A[3][1]=-0.125*(aj[k-1][j][i]*ucont[k-1][j][i].z);
	A[3][2]=0.125*(aj[k][j+1][i]*ucont[k][j+1][i].z); A[3][3]=-0.125*(aj[k-1][j+1][i]*ucont[k-1][j+1][i].z);

	A[4][0]=-0.125*(aj[k+1][j][i-1]*ucont[k+1][j][i-1].x); A[4][1]=-0.125*(aj[k][j][i-1]*ucont[k][j][i-1].x);
	A[4][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].x); A[4][3]=0.125*(aj[k][j][i]*ucont[k][j][i].x);

	A[5][0]=-0.125*(aj[k+1][j-1][i]*ucont[k+1][j-1][i].y); A[5][1]=-0.125*(aj[k][j-1][i]*ucont[k][j-1][i].y); 
	A[5][2]=0.125*(aj[k+1][j][i]*ucont[k+1][j][i].y); A[5][3]=0.125*(aj[k][j][i]*ucont[k][j][i].y);

	  kk = lidxLocal(i, j, k, user);
	  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

	  //---------------------------------------------
	  dd=lidxLocal(i, j, k, user);  //(i,j,k)
	  //--calc needed
	  m=i,n=j,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--

	  col=3*dd;
	  val[0]=-(dtc+AJip*AJip*(g11ip+g22ip+g33ip)/Re+A[0][0]+A[0][1]+A[0][2]+A[0][3]+A[1][0]+A[1][1]+A[1][2]+A[1][3]);
	  val[1]=-(0.5*AJip*U1ip);
	  val[2]=-(0.5*AJip*U2ip);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(0.5*AJjp*U0jp);
	  val[1]=-(dtc+AJjp*AJjp*(g11jp+g22jp+g33jp)/Re+A[2][0]+A[2][1]+A[2][2]+A[2][3]+A[3][0]+A[3][1]+A[3][2]+A[3][3]);
	  val[2]=-(0.5*AJjp*U2jp);;	  
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(0.5*AJkp*U0kp);;
	  val[1]=-(0.5*AJkp*U1kp);;
	  val[2]=-(dtc+AJkp*AJkp*(g11kp+g22kp+g33kp)/Re+A[4][0]+A[4][1]+A[4][2]+A[4][3]+A[5][0]+A[5][1]+A[5][2]+A[5][3]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //----------------------------------------------
	  dd=lidxLocal(i+1, j, k, user);  //(i+1,j,k)
	  //--calc needed
	  m=i+1,n=j,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g11ip/Re)+(c2*ucont[p][n][m].x)+A[0][2]+A[0][3]+A[1][2]+A[1][3]);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(0.5*AJjp*U0jp);
	  val[1]=-((-0.5*AJjp*AJjp*g11jp/Re)+A[2][2]+A[2][3]);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(0.5*AJkp*U0kp);
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g11kp/Re)+A[4][2]+A[4][3]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i-1, j, k, user);  //(i-1,j,k)
	  //--calc needed
	  m=i-1,n=j,p=k;
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  //--
	  col=3*dd;
	  val[0]=-((-0.5*AJip*AJip*g11ip/Re)+(-c2*ucont[p][n][m].x)+A[0][0]+A[0][1]+A[1][0]+A[1][1]);
	  val[1]=-(-0.5*AJip*U1ip);
	  val[2]=-(-0.5*AJip*U2ip);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g11jp/Re)+A[2][0]+A[2][1]);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g11kp/Re)+A[4][0]+A[4][1]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------	
	  dd=lidxLocal(i, j+1, k, user);  //(i,j+1,k)
	  //--calc needed
	  m=i,n=j+1,p=k;
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g22ip/Re)+A[0][0]+A[0][2]);
	  val[1]=-(0.5*AJip*U1ip);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g22jp/Re)+(c2*ucont[p][n][m].y)+A[2][2]+A[2][0]+A[3][2]+A[3][3]);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=-(0.5*AJkp*U1kp);
	  val[2]=-((-0.5*AJkp*AJkp*g22kp/Re)+A[5][2]+A[5][3]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j-1, k, user);  //(i,j-1,k)
	  //--calc needed
	  m=i,n=j-1,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g22ip/Re)+A[0][1]+A[0][3]);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=-(-0.5*AJjp*U0jp);
	  val[1]=-((-0.5*AJjp*AJjp*g22jp/Re)+(-c2*ucont[p][n][m].y)+A[2][3]+A[2][1]+A[3][0]+A[3][1]);
	  val[2]=-(-0.5*AJjp*U2jp);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g22kp/Re)+A[5][0]+A[5][1]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j, k+1, user);  //(i,j,k+1)
	  //--calc needed
	  m=i,n=j,p=k+1;
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--
	  col=3*dd;
	  val[0]=-((-0.5*AJip*AJip*g33ip/Re+A[1][0]+A[1][2]));
	  val[1]=0.0;
	  val[2]=-(0.5*AJip*U2ip);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g33jp/Re)+A[3][0]+A[3][2]);
	  val[2]=-(0.5*AJjp*U2jp);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-((-0.5*AJkp*AJkp*g33kp/Re)+(c2*ucont[p][n][m].z)+A[4][0]+A[4][2]+A[5][0]+A[5][2]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j, k-1, user);  //(i,j,k-1)
	  //--calc needed
	  m=i,n=j,p=k-1;
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  //--
	  col=3*dd; 
	  val[0]=-((-0.5*AJip*AJip*g33ip/Re+A[1][1]+A[1][3]));
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-((-0.5*AJjp*AJjp*g33jp/Re)+A[3][1]+A[3][3]);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	
	  col=3*dd+2;
	  val[0]=-(-0.5*AJkp*U0kp);
	  val[1]=-(-0.5*AJkp*U1kp);
	  val[2]=-((-0.5*AJkp*AJkp*g33kp/Re)+(-c2*ucont[p][n][m].z)+A[4][1]+A[4][3]+A[5][1]+A[5][3]);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i-1, j+1, k, user);  //(i-1,j+1,k)
	  //--calc needed
	  m=i-1,n=j+1,p=k;
	  U1ip=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p][n][m+1].y+ucont[p][n-1][m+1].y);
	  //--
	  col=3*dd;
	  val[0]=-A[0][0];
	  val[1]=-(-0.5*AJip*U1ip);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[2][0];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i-1, j, k+1, user);  //(i-1,j,k+1)
	  //--calc needed
	  m=i-1,n=j,p=k+1;
	  U2ip=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n][m+1].z+ucont[p-1][n][m+1].z);
	  //--
	  col=3*dd;
	  val[0]=-A[1][0];
	  val[1]=0.0;
	  val[2]=-(-0.5*AJip*U2ip);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[4][0];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i+1, j-1, k, user);  //(i+1,j-1,k)
	  //--calc needed
	  m=i+1,n=j-1,p=k;
	  U0jp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p][n+1][m].x+ucont[p][n+1][m-1].x);
	  //--
	  col=3*dd;
	  val[0]=-A[0][3];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=-(-0.5*AJjp*U0jp);
	  val[1]=-A[2][3];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j-1, k+1, user);  //(i,j-1,k+1)
	  //--calc needed
	  m=i,n=j-1,p=k+1;
	  U2jp=0.25*(ucont[p][n][m].z+ucont[p-1][n][m].z+ucont[p][n+1][m].z+ucont[p-1][n+1][m].z);
	  //--
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][0];
	  val[2]=-(-0.5*AJjp*U2jp);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[5][0];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i+1, j, k-1, user);  //(i+1,j,k-1)
	  //--calc needed
	  m=i+1,n=j,p=k-1;
	  U0kp=0.25*(ucont[p][n][m].x+ucont[p][n][m-1].x+ucont[p+1][n][m].x+ucont[p+1][n][m-1].x);
	  //--
	  col=3*dd;
	  val[0]=-A[1][3];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=-(-0.5*AJkp*U0kp);
	  val[1]=0.0;
	  val[2]=-A[4][3];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j+1, k-1, user);  //(i,j+1,k-1)
	  //--calc needed
	  m=i,n=j+1,p=k-1;
	  U1kp=0.25*(ucont[p][n][m].y+ucont[p][n-1][m].y+ucont[p+1][n][m].y+ucont[p+1][n-1][m].y);
	  //--
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][3];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=-(-0.5*AJkp*U1kp);
	  val[2]=-A[5][3];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i-1, j-1, k, user);  //(i-1,j-1,k)
	  m=i-1, n=j-1, p=k;
	  col=3*dd;
	  val[0]=-A[0][1];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[2][1];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i+1, j+1, k, user);  //(i+1,j+1,k)
	  m=i+1, n=j+1, p=k;
	  col=3*dd;
	  val[0]=-A[0][2];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[2][2];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i-1, j, k-1, user);  //(i-1,j,k-1)
	  m=i-1, n=j, p=k-1;
	  col=3*dd;
	  val[0]=-A[1][1];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[4][1];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i+1, j, k+1, user);  //(i+1,j,k+1)
	  m=i+1, n=j, p=k+1;
	  col=3*dd;
	  val[0]=-A[1][2];
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[4][2];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j-1, k-1, user);  //(i,j-1,k-1)
	  m=i, n=j-1, p=k-1;
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][1];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[5][1];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  dd=lidxLocal(i, j+1, k+1, user);  //(i,j+1,k+1)
	  m=i, n=j+1, p=k+1;
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-A[3][2];
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);

	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-A[5][2];
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  //---------------------------------------------
	  if(i != 1){
	    m=i-2,n=j,p=k;
	  dd=lidxLocal(i-2, j, k, user);  //(i-2,j,k)
	  col=3*dd;
	  val[0]=-(c1*ucont[p][n][m].x);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  }
	  //---------------------------------------------
	  if(j != 1){
	    m=i,n=j-2,p=k;
	  dd=lidxLocal(i, j-2, k, user);  //(i,j-2,k)
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(c1*ucont[p][n][m].y);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  }
	  //---------------------------------------------
	  if(k != 1){
	    m=i,n=j,p=k-2;
	  dd=lidxLocal(i, j, k-2, user);  //(i,j,k-2)
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(c1*ucont[p][n][m].z);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  }
	  //---------------------------------------------
	  if(i != mx-2){
	    m=i+2,n=j,p=k;
	  dd=lidxLocal(i+2, j, k, user);  //(i+2,j,k)
	  col=3*dd;
	  val[0]=-(-c1*ucont[p][n][m].x);
	  val[1]=0.0;
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  }
	  //---------------------------------------------
	  if(j != my-2){
	    m=i,n=j+2,p=k;
	  dd=lidxLocal(i, j+2, k, user);  //(i,j+2,k)
	  col=3*dd+1;
	  val[0]=0.0;
	  val[1]=-(-c1*ucont[p][n][m].y);
	  val[2]=0.0;
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  }
	  //---------------------------------------------
	  if(k != mz-2){
	    m=i,n=j,p=k+2;
	  dd=lidxLocal(i, j, k+2, user);  //(i,j,k+2)
	  col=3*dd+2;
	  val[0]=0.0;
	  val[1]=0.0;
	  val[2]=-(-c1*ucont[p][n][m].z);
	  MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
	  CheckCol_packer(J,user,m,n,p);
	  }
	  //---------------------------------------------
	 
	  if (i==mx-2) {
	    CheckRow_packer(J,user,i,j,k,0);
	  }
	  if (j==my-2) {
	    CheckRow_packer(J,user,i,j,k,1);
	  }
	  if (k==mz-2) {
	    CheckRow_packer(J,user,i,j,k,2);
	  }
	/*   if (0.5*(nvert[k][j][i]+nvert[k][j][i+1]) > 0.1) { */
/* 	    CheckRow_packer(J,user,i,j,k,0); */
/* 	  } */
/* 	  if (0.5*(nvert[k][j][i]+nvert[k][j+1][i]) > 0.1) { */
/* 	    CheckRow_packer(J,user,i,j,k,1); */
/* 	  } */
/* 	  if (0.5*(nvert[k][j][i]+nvert[k+1][j][i]) > 0.1) { */
/* 	    CheckRow_packer(J,user,i,j,k,2); */
/* 	  } */


      }
    }
  }
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);

 

  ModifyJacobian_packer(J,user);
  
  MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);
  

 return 0;
}
PetscErrorCode CheckRow_packer(Mat J,UserCtx *ctx,PetscInt i,PetscInt j,PetscInt k,PetscInt direction)
{
  PetscReal val[3];
  PetscInt dd,kk,col[3],row,m,n,p;
  UserCtx *user = (UserCtx*)ctx;
  DMDALocalInfo	info = user->info;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  val[0]=0.0;
  val[1]=0.0;
  val[2]=0.0;

  kk = lidxLocal(i, j, k, user);
  row=3*kk+direction;
 //----------------------------------------------------
  m=i,n=j,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i+1,n=j,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i-1,n=j,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i,n=j+1,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i,n=j-1,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i,n=j,p=k+1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES); 
  //----------------------------------------------------
  m=i,n=j,p=k-1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i-1,n=j+1,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i-1,n=j,p=k+1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i+1,n=j-1,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i,n=j-1,p=k+1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i+1,n=j,p=k-1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i,n=j+1,p=k-1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i-1,n=j-1,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i+1, n=j+1,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i-1, n=j, p=k-1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i+1, n=j, p=k+1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i, n=j-1, p=k-1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  m=i, n=j+1, p=k+1;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
 //----------------------------------------------------
  if(i != 1){
    m=i-2, n=j, p=k;
    dd=lidxLocal(m,n,p,user); 
    col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
    MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
  }
 //----------------------------------------------------
  if(j != 1){
    m=i,n=j-2,p=k;
    dd=lidxLocal(m,n,p,user); 
    col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
    MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
  }
 // ----------------------------------------------------
   if(k != 1){
     m=i,n=j,p=k-2;
     dd=lidxLocal(m,n,p,user); 
     col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
     MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
   }
  //----------------------------------------------------
  if(i != mx-2){
    m=i+2,n=j,p=k;
    dd=lidxLocal(m,n,p,user); 
    col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
    MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
  }
 //----------------------------------------------------
  if(j != my-2){
    m=i,n=j+2,p=k;
    dd=lidxLocal(m,n,p,user); 
    col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
    MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
  }
 //----------------------------------------------------
  if(k != mz-2){
    m=i,n=j,p=k+2;
    dd=lidxLocal(m,n,p,user); 
    col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
    MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
  }
 //----------------------------------------------------

 return 0;
}

PetscErrorCode CheckRowDiagonal_packer(Mat J,UserCtx *ctx,PetscInt i,PetscInt j,PetscInt k,PetscInt direction)
{
  PetscReal val[3];
  PetscInt dd,kk,col[3],row,m,n,p;
  UserCtx *user = (UserCtx*)ctx;
  DMDALocalInfo	info = user->info;

  val[0]=0.0;
  val[1]=0.0;
  val[2]=0.0;

  kk = lidxLocal(i, j, k, user);
  row=3*kk+direction;
 //----------------------------------------------------
  m=i,n=j,p=k;
  dd=lidxLocal(m,n,p,user); 
  col[0]=3*dd; col[1]=3*dd+1; col[2]=3*dd+2;
  MatSetValuesLocal(J,1,&row,3,col,val,INSERT_VALUES);
  //----------------------------------------------------

 return 0;
}
PetscErrorCode CheckCol_packer(Mat J,UserCtx *ctx,PetscInt m,PetscInt n,PetscInt p)
{
  PetscReal val[3];
  PetscInt dd,kk,col,row[3],direction;
  UserCtx *user = (UserCtx*)ctx;
  DMDALocalInfo	info = user->info;
  DM fda = user->fda;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  Cmpnts ***nvertF;

  val[0]=0.0;
  val[1]=0.0;
  val[2]=0.0;

 
  kk = lidxLocal(m,n,p,user);
  row[0]=3*kk; row[1]=3*kk+1;   row[2]=3*kk+2;

  DMDAVecGetArray(fda, user->lNFace, &nvertF);


 /*  if ((nvertF[p][n][m].x > 0.1) || (m==mx-2)) { */
/*     dd=lidxLocal(m,n,p,user); */
/*     col=3*dd; */
/*     MatSetValuesLocal(*J,3,row,1,&col,val,INSERT_VALUES); */
/*   }  */
/*   if ((nvertF[p][n][m].y > 0.1) || (n==my-2)) { */
/*     dd=lidxLocal(m,n,p,user); */
/*     col=3*dd+1; */
/*     MatSetValuesLocal(*J,3,row,1,&col,val,INSERT_VALUES); */
/*   }  */
/*   if ((nvertF[p][n][m].z > 0.1) || (p==mz-2)) { */
/*     dd=lidxLocal(m,n,p,user); */
/*     col=3*dd+2; */
/*     MatSetValuesLocal(*J,3,row,1,&col,val,INSERT_VALUES); */
/*   } */

 if ((m==mx-2)) {
    dd=lidxLocal(m,n,p,user);
    col=3*dd;
    MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
  } 
  if ( (n==my-2)) {
    dd=lidxLocal(m,n,p,user);
    col=3*dd+1;
    MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
  } 
  if ((p==mz-2)) {
    dd=lidxLocal(m,n,p,user);
    col=3*dd+2;
    MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
  }

 DMDAVecRestoreArray(fda, user->lNFace, &nvertF);

  if ( (m==0)|| (m==mx-1)||(n==0)|| (n==my-1)||(p==0)|| (p==mz-1) ) {
    for(direction=0;direction<3;direction++){
      dd=lidxLocal(m,n,p,user);
      col=3*dd+direction;
      MatSetValuesLocal(J,3,row,1,&col,val,INSERT_VALUES);
    }
  }
 


 return 0;
}
PetscErrorCode ComputelNface(UserCtx *ctx){


  UserCtx *user = (UserCtx*)ctx;
  DM da = user->da;
  DM fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert;
  Cmpnts ***nvertF;
  PetscInt	i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

 
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->lNFace, &nvertF);
 
 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
       if (i==mx-1) {
	 nvertF[k][j][i].x=nvert[k][j][i];
       } else if (j==my-1 ) {
	 nvertF[k][j][i].y=nvert[k][j][i];
       }else if (k==mz-1) {
	 nvertF[k][j][i].z=nvert[k][j][i];
       } else{
	 nvertF[k][j][i].x=0.5*(nvert[k][j][i]+nvert[k][j][i+1]);
	 nvertF[k][j][i].y=0.5*(nvert[k][j][i]+nvert[k][j+1][i]);
	 nvertF[k][j][i].z=0.5*(nvert[k][j][i]+nvert[k+1][j][i]);
       }
       
     }
    }
  }
  
 DMDAVecRestoreArray(fda, user->lNFace, &nvertF);
 DMDAVecRestoreArray(da, user->lNvert, &nvert);





  return 0;
}
PetscErrorCode ModifyJacobian_packer(Mat PJ,UserCtx *ctx)
{
  UserCtx *user = (UserCtx*)ctx;
  DM fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     Unit=1.0;
  PetscInt	i, j, k,row,kk;
  Cmpnts        ***nvertF;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

/*   if (xs==0) lxs = xs+1; */
/*   if (ys==0) lys = ys+1; */
/*   if (zs==0) lzs = zs+1; */

/*   if (xe==mx) lxe = xe-1; */
/*   if (ye==my) lye = ye-1; */
/*   if (ze==mz) lze = ze-1; */
  

  DMDAVecGetArray(fda, user->lNFace, &nvertF);

 
  if (zs == 0) {
    k = 0;  //physical ghost and boundary node
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	
	row = lidxLocal(i, j, k, user);
	kk=3*row;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);

	if (user->bctype[4]!=7) {
	kk=3*row+2;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
      }
    }
  }

  if (ze == mz) {
    k = mz-1;  //physical ghost node
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	row = lidxLocal(i, j, k, user);
	kk=3*row;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
      }
    }
    k = mz-2; //Boundary node
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	row = lidxLocal(i, j, k, user);
	if (user->bctype[5]!=7) {
	  kk=3*row+2;
	 MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
      }
    }
  }

  if (ys == 0) {
    j = 0;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	row = lidxLocal(i, j, k, user);
	kk=3*row;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);

	if (user->bctype[2]!=7) {
	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}

	kk=3*row+2;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
      }
    }
  }

  if (ye == my) {
    j = my-1;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	row = lidxLocal(i, j, k, user);
	kk=3*row;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
      }
    }
    j = my-2;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	row = lidxLocal(i,j,k,user);
	if (user->bctype[3]!=7) {
       	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	
      }
    }
  }

  if (xs == 0) {
    i = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	row = lidxLocal(i, j, k, user);
	if (user->bctype[0]!=7) {
	kk=3*row;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
      }
    }
  }

  if (xe == mx) {
    i = mx-1;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	row = lidxLocal(i, j, k, user);
	kk=3*row;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+1;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	kk=3*row+2;
	MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
      }
    }
    i = mx-2;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	row = lidxLocal(i, j, k, user);
	if (user->bctype[1]!=7) {
	  kk=3*row;
	  MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
      }
    }
  }

  for (k=lzs; k<lze; k++) {   //Immersed Boundary
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
   	if (nvertF[k][j][i].x > 0.1) {
	 row = lidxLocal(i, j, k, user);
	 kk=3*row;
	 MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	
	 
      }
    }
  }
 for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvertF[k][j][i].y > 0.1) {
	  row = lidxLocal(i, j, k, user);
	  kk=3*row+1;
	  MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
 
      }
    }
  }
 for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvertF[k][j][i].z > 0.1) {
	  row = lidxLocal(i, j, k, user);
	  kk=3*row+2;
	  MatSetValuesLocal(PJ,1,&kk,1,&kk,&Unit,INSERT_VALUES);
	}
	 
      }
    }
  }


  DMDAVecRestoreArray(fda, user->lNFace, &nvertF);
  
  


return 0;
}
PetscInt lidxLocal(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)
{
  DMDALocalInfo	info = user->info;


  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;


 
  return ((k-gzs) * (info.gxm*info.gym) + (j-gys)*(info.gxm) + (i-gxs));
  
}
