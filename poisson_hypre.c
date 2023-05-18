#include "variables.h"
#include "petscksp.h"
#include "petscpc.h"
//#include "petscpcmg.h"
#include <stdlib.h>
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_IJ_mv.h"
#include "petsctime.h"

#define lidx2(i,j,k,user)	(int)(gid[k][j][i])

extern int periodic, i_periodic, j_periodic, k_periodic, pseudo_periodic;
extern PetscInt block_number, freesurface, immersed, ti, tistart, les, tiout;
extern PetscInt tistart;
extern double poisson_tol;
//extern char path[256];
extern double mean_pressure_gradient;
//extern PetscInt movefsi, rotatefsi;

PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg);

HYPRE_Solver *pcg_solver_p, *precon_p;
HYPRE_IJMatrix *Ap;
HYPRE_ParCSRMatrix *par_Ap;
HYPRE_IJVector *Vec_p, *Vec_p_rhs;
HYPRE_ParVector *par_Vec_p, *par_Vec_p_rhs;
PetscBool      hypre_symmetric;


PetscErrorCode Initialize_Hypre_Var(PetscInt nb)
{
  PetscMalloc(nb*sizeof(HYPRE_Solver), &pcg_solver_p);
  PetscMalloc(nb*sizeof(HYPRE_Solver), &precon_p);
  PetscMalloc(nb*sizeof(HYPRE_IJMatrix), &Ap);
  PetscMalloc(nb*sizeof(HYPRE_ParCSRMatrix), &par_Ap);
  PetscMalloc(nb*sizeof(HYPRE_IJVector), &Vec_p);
  PetscMalloc(nb*sizeof(HYPRE_IJVector), &Vec_p_rhs);
  PetscMalloc(nb*sizeof(HYPRE_ParVector), &par_Vec_p);
  PetscMalloc(nb*sizeof(HYPRE_ParVector), &par_Vec_p_rhs);

  return(0);
}

double time_coeff()
{
	if(ti<=tistart+1) return 1;
	else if(ti<=tistart+2) return 1.5;
	else return 11./6.;		// 1.8333
	
};

void Convert_Phi2_Phi(UserCtx *user)
{
	DMDALocalInfo	info = user->info;
	DM		da = user->da;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, ***phi, ***lphi, *phi2;
	PetscReal poisson_threshold=0.1;
	
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	VecSet(user->Phi,0);
	DMDAVecGetArray(da, user->Phi, &phi);
	DMDAVecGetArray(da, user->lNvert, &nvert);
	VecGetArray(user->Phi2, &phi2);
	
	PetscOptionsHasName(PETSC_NULL,"-hypre_symmetric",&hypre_symmetric);
 
	int pos=0;
	for(k=lzs; k<lze; k++)
	for(j=lys; j<lye; j++)
	for(i=lxs; i<lxe; i++) {
		/*if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
			phi[k][j][i]=0;
		}
		else */if( (int)nvert[k][j][i]>poisson_threshold ) {
			if(!hypre_symmetric) {
				phi[k][j][i] = phi2[pos++];
			}
			else phi[k][j][i]=0;
		}
		else {
			phi[k][j][i] = phi2[pos++];
		}
	}
	
	DMDAVecRestoreArray(da, user->Phi, &phi);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	VecRestoreArray(user->Phi2, &phi2);
	
	DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
	DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);
	
	DMDAVecGetArray(da, user->lPhi, &lphi);
	DMDAVecGetArray(da, user->Phi, &phi);
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
		
	/* 	if(ii_periodic && i==0) a=-2, flag=1; */
/* 		else if(ii_periodic && i==mx-1) a=mx+1, flag=1; */
		
/* 		if(jj_periodic && j==0) b=-2, flag=1; */
/* 		else if(jj_periodic && j==my-1) b=my+1, flag=1; */
		
/* 		if(kk_periodic && k==0) c=-2, flag=1; */
/* 		else if(kk_periodic && k==mz-1) c=mz+1, flag=1; */
		
		if(flag) {
			phi[k][j][i] = lphi[c][b][a];
		}
		
		//if(k==1 && i!=0)	printf("%d,%d,%d, %f %f %f\n", i,j,k, lphi[-2][j][i], lphi[k][j][i], lphi[k+1][j][i]);
	}
	DMDAVecRestoreArray(da, user->lPhi, &lphi);
	DMDAVecRestoreArray(da, user->Phi, &phi);
	//exit(0);
}

PetscInt setup_lidx2(UserCtx *user)
{
	DMDALocalInfo	info = user->info;
	DM		da = user->da;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	int i, j, k;
	PetscReal ***nvert, ***gid, ***lid;

	Vec	Lid;

	PetscReal poisson_threshold=0.1;
	
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	VecDuplicate(user->lNvert, &user->Gid);
	VecSet(user->Gid, -1);
	
	VecDuplicate(user->lNvert, &Lid);
	VecSet(Lid, -1);
	
	DMDAVecGetArray(da, user->Gid, &gid);
	DMDAVecGetArray(da, Lid, &lid);
	DMDAVecGetArray(da, user->lNvert, &nvert);

	int r, myrank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
		
	int ndof_node[size], ndof_node_tmp[size];
	//std::vector<int> ndof_node(size), ndof_node_tmp(size);	// # of pressure dof for processors
	
	int ndof_node_accu;
	PetscOptionsHasName(PETSC_NULL,"-hypre_symmetric",&hypre_symmetric);

	for(r=0; r<size; r++) {
		ndof_node_tmp[r] = 0;
	}
	
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
			//lid[k][j][i]=-10;
		}
		else if( (int)nvert[k][j][i]>poisson_threshold ) {
			if(!hypre_symmetric) {
				lid[k][j][i] = (PetscReal)ndof_node_tmp[myrank]++;
				//ndof_node_tmp[myrank] += 1;
			}
			else {}//lid[k][j][i]=-10;
		}
		else {
			lid[k][j][i] = (PetscReal)ndof_node_tmp[myrank]++;
			//ndof_node_tmp[myrank] += 1;
		}
	}
	
	MPI_Allreduce( &ndof_node_tmp[0], &ndof_node[0], size, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
		
	ndof_node_accu = 0;
	for(r=0; r<myrank; r++) ndof_node_accu += ndof_node[r];
		
	int n;
	user->p_global_begin = ndof_node_accu;
		
	VecGetSize(user->Phi,&n);
	if(myrank==size-1) {
		printf("\n\n********* %d %d ***********\n\n", ndof_node_accu + ndof_node[myrank], n);
		user->reduced_p_size = ndof_node_accu + ndof_node[myrank];
	}
		
	MPI_Bcast(&user->reduced_p_size, 1, MPI_INT, size-1, PETSC_COMM_WORLD);
		
	PetscBarrier(PETSC_NULL);
		
	for(k=zs; k<ze; k++)
	for(j=ys; j<ye; j++)
	for(i=xs; i<xe; i++) {
		if((int)(lid[k][j][i])>=0) {
			gid[k][j][i] = lid[k][j][i] + ndof_node_accu;	// gid is double, be careful
		}
	}
		
	user->local_Phi2_size = ndof_node[myrank];
	
	if(ti!=tistart) VecDestroy (&user->Phi2);
	
	VecCreateMPI(PETSC_COMM_WORLD, ndof_node[myrank], PETSC_DETERMINE, &user->Phi2);
		
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->Gid, &gid);
	DMDAVecRestoreArray(da, Lid, &lid);
	
	VecDestroy(&Lid);

        DMDALocalToLocalBegin(da, user->Gid, INSERT_VALUES, user->Gid);
        DMDALocalToLocalEnd(da, user->Gid, INSERT_VALUES, user->Gid);
	
	if(periodic) {
		DMDAVecGetArray(da, user->Gid, &gid);
		for(k=zs; k<ze; k++)
		for(j=ys; j<ye; j++)
		for(i=xs; i<xe; i++) {
			int flag=0, a=i, b=j, c=k;
			
			if(i_periodic && i==0) a=mx-2, flag=1;
			else if(i_periodic && i==mx-1) a=1, flag=1;
			
			if(j_periodic && j==0) b=my-2, flag=1;
			else if(j_periodic && j==my-1) b=1, flag=1;
			
			if(k_periodic && k==0) c=mz-2, flag=1;
			else if(k_periodic && k==mz-1) c=1, flag=1;
			
			/* if(ii_periodic && i==0) a=-2, flag=1; */
/* 			else if(ii_periodic && i==mx-1) a=mx+1, flag=1; */
			
/* 			if(jj_periodic && j==0) b=-2, flag=1; */
/* 			else if(jj_periodic && j==my-1) b=my+1, flag=1; */
			
/* 			if(kk_periodic && k==0) c=-2, flag=1; */
/* 			else if(kk_periodic && k==mz-1) c=mz+1, flag=1; */
			
			if(flag) gid[k][j][i] = gid[c][b][a];
			
			//if(i==1 && k==mz-2) printf("%d, %d, %d, %d <= %d %d %d \n", (int)gid[k][j][i], k,j,i, c,b,a);
		}
		DMDAVecRestoreArray(da, user->Gid, &gid);
		
		DMDALocalToLocalBegin(da, user->Gid, INSERT_VALUES, user->Gid);
		DMDALocalToLocalEnd(da, user->Gid, INSERT_VALUES, user->Gid);
	}	
	
	return 0;
}

#define CP 0
#define EP 1
#define WP 2
#define NP 3
#define SP 4
#define TP 5
#define BP 6
#define NE 7
#define SE 8
#define NW 9
#define SW 10
#define TN 11
#define BN 12
#define TS 13
#define BS 14
#define TE 15
#define BE 16
#define TW 17
#define BW 18

#define AA 19
#define BB 20
#define CC 21
#define DD 22
#define EE 23
#define FF 24
#define GG 25 
#define HH 26

double buffer[5][5][5]; 

#define	SetValue(i,j,k,v) 	buffer[i+2][j+2][k+2]=(v);
#define AddValue(i,j,k,v) 	buffer[i+2][j+2][k+2]+=(v);

PetscErrorCode PoissonLHS_HYPRE(UserCtx *user)
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

	PetscReal	***nvert, ***nvert_o, ***gid;
	PetscScalar	vol[27];
	PetscInt	idx[27], row;
	PetscReal       poisson_threshold=0.1;

	PetscInt	i, j, k, N;
	PetscInt	p,r,q;

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
		user->assignedA = PETSC_TRUE;
		
		setup_lidx2(user);
		N = mx * my * mz;
		PetscInt M;
		
		MatCreate(PETSC_COMM_WORLD, &(user->A));
		VecGetLocalSize(user->Phi2, &M);
		MatSetSizes(user->A,M,M,PETSC_DETERMINE,PETSC_DETERMINE);
		MatSetType(user->A,MATMPIAIJ);
		MatMPIAIJSetPreallocation(user->A, 19, PETSC_NULL, 19, PETSC_NULL);
		MatSetFromOptions(user->A);
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

	VecDuplicate(user->lAj, &G11);
	VecDuplicate(user->lAj, &G12);
	VecDuplicate(user->lAj, &G13);
	VecDuplicate(user->lAj, &G21);
	VecDuplicate(user->lAj, &G22);
	VecDuplicate(user->lAj, &G23);
	VecDuplicate(user->lAj, &G31);
	VecDuplicate(user->lAj, &G32);
	VecDuplicate(user->lAj, &G33);
/*	
	VecSet(G11,1.e10);
	VecSet(G12,1.e10);
	VecSet(G13,1.e10);
	VecSet(G21,1.e10);
	VecSet(G22,1.e10);
	VecSet(G23,1.e10);
	VecSet(G31,1.e10);
	VecSet(G32,1.e10);
	VecSet(G33,1.e10);
*/
	DMDAVecGetArray(da, G11, &g11);
	DMDAVecGetArray(da, G12, &g12);
	DMDAVecGetArray(da, G13, &g13);
	DMDAVecGetArray(da, G21, &g21);
	DMDAVecGetArray(da, G22, &g22);
	DMDAVecGetArray(da, G23, &g23);
	DMDAVecGetArray(da, G31, &g31);
	DMDAVecGetArray(da, G32, &g32);
	DMDAVecGetArray(da, G33, &g33);
	
	/*for (k=gzs; k<gze; k++)
	for (j=gys; j<gye; j++)
	for (i=gxs; i<gxe; i++) */
	for (k=lzs-1; k<lze+1; k++)
	for (j=lys-1; j<lye+1; j++)
	for (i=lxs-1; i<lxe+1; i++) 
	{
		int a=i, b=j, c=k;
		
		
			int i_flag=0, j_flag=0, k_flag=0;
			
			if(i_periodic && i==0) a=mx-2, i_flag=1;
			else if(i_periodic && i==mx-1) a=1, i_flag=1;
			
			if(j_periodic && j==0) b=my-2, j_flag=1;
			else if(j_periodic && j==my-1) b=1, j_flag=1;
			
			if(k_periodic && k==0) c=mz-2, k_flag=1;
			else if(k_periodic && k==mz-1) c=1, k_flag=1;
			
			/* if(ii_periodic && i==0) a=-2, i_flag=1; */
/* 			else if(ii_periodic && i==mx-1) a=mx+1, i_flag=1; */
			
/* 			if(jj_periodic && j==0) b=-2, j_flag=1; */
/* 			else if(jj_periodic && j==my-1) b=my+1, j_flag=1; */
			
/* 			if(kk_periodic && k==0) c=-2, k_flag=1; */
/* 			else if(kk_periodic && k==mz-1) c=mz+1, k_flag=1; */
			
		g11[k][j][i] = (icsi[c][b][a].x * icsi[c][b][a].x + icsi[c][b][a].y * icsi[c][b][a].y + icsi[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
		g12[k][j][i] = (ieta[c][b][a].x * icsi[c][b][a].x + ieta[c][b][a].y * icsi[c][b][a].y + ieta[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
		g13[k][j][i] = (izet[c][b][a].x * icsi[c][b][a].x + izet[c][b][a].y * icsi[c][b][a].y + izet[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
		g21[k][j][i] = (jcsi[c][b][a].x * jeta[c][b][a].x + jcsi[c][b][a].y * jeta[c][b][a].y + jcsi[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
		g22[k][j][i] = (jeta[c][b][a].x * jeta[c][b][a].x + jeta[c][b][a].y * jeta[c][b][a].y + jeta[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
		g23[k][j][i] = (jzet[c][b][a].x * jeta[c][b][a].x + jzet[c][b][a].y * jeta[c][b][a].y + jzet[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
		g31[k][j][i] = (kcsi[c][b][a].x * kzet[c][b][a].x + kcsi[c][b][a].y * kzet[c][b][a].y + kcsi[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
		g32[k][j][i] = (keta[c][b][a].x * kzet[c][b][a].x + keta[c][b][a].y * kzet[c][b][a].y + keta[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
		g33[k][j][i] = (kzet[c][b][a].x * kzet[c][b][a].x + kzet[c][b][a].y * kzet[c][b][a].y + kzet[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];	
		/*
		if(i==0 && j==2) printf("%d,%d,%d, (%d,%d,%d) %.3f %.3f %.3f   %.3f %.3f %.3f   %.3f %.3f %.3f\n", i, j, k, a,b,c,
							g11[k][j][i], g12[k][j][i], g13[k][j][i],
							g21[k][j][i], g22[k][j][i], g23[k][j][i],
							g31[k][j][i], g32[k][j][i], g33[k][j][i]
						);*/
	}
	
	//exit(0);
	/*
	
	if( periodic ) {
		DMDAVecRestoreArray(da, G11, &g11);
		DMDAVecRestoreArray(da, G12, &g12);
		DMDAVecRestoreArray(da, G13, &g13);
		DMDAVecRestoreArray(da, G21, &g21);
		DMDAVecRestoreArray(da, G22, &g22);
		DMDAVecRestoreArray(da, G23, &g23);
		DMDAVecRestoreArray(da, G31, &g31);
		DMDAVecRestoreArray(da, G32, &g32);
		DMDAVecRestoreArray(da, G33, &g33);
		    
		DALocalToLocalBegin(da, G11, INSERT_VALUES, G11);
		DALocalToLocalEnd(da, G11, INSERT_VALUES, G11);
		    
		DALocalToLocalBegin(da, G12, INSERT_VALUES, G12);
		DALocalToLocalEnd(da, G12, INSERT_VALUES, G12);
		
		DALocalToLocalBegin(da, G13, INSERT_VALUES, G13);
		DALocalToLocalEnd(da, G13, INSERT_VALUES, G13);
		    
		DALocalToLocalBegin(da, G21, INSERT_VALUES, G21);
		DALocalToLocalEnd(da, G21, INSERT_VALUES, G21);
		
		DALocalToLocalBegin(da, G22, INSERT_VALUES, G22);
		DALocalToLocalEnd(da, G22, INSERT_VALUES, G22);
		
		DALocalToLocalBegin(da, G23, INSERT_VALUES, G23);
		DALocalToLocalEnd(da, G23, INSERT_VALUES, G23);
		
		DALocalToLocalBegin(da, G31, INSERT_VALUES, G31);
		DALocalToLocalEnd(da, G31, INSERT_VALUES, G31);
		
		DALocalToLocalBegin(da, G32, INSERT_VALUES, G32);
		DALocalToLocalEnd(da, G32, INSERT_VALUES, G32);
		
		DALocalToLocalBegin(da, G33, INSERT_VALUES, G33);
		DALocalToLocalEnd(da, G33, INSERT_VALUES, G33);

		DMDAVecGetArray(da, G11, &g11);
		DMDAVecGetArray(da, G12, &g12);
		DMDAVecGetArray(da, G13, &g13);
		DMDAVecGetArray(da, G21, &g21);
		DMDAVecGetArray(da, G22, &g22);
		DMDAVecGetArray(da, G23, &g23);
		DMDAVecGetArray(da, G31, &g31);
		DMDAVecGetArray(da, G32, &g32);
		DMDAVecGetArray(da, G33, &g33);
		
		for (k=zs; k<ze; k++)
		for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
		      
			int a=i, b=j, c=k;
			int i_flag=0, j_flag=0, k_flag=0;
			
			if(i_periodic && i==0) a=mx-2, i_flag=1;
			else if(i_periodic && i==mx-1) a=1, i_flag=1;
			
			if(j_periodic && j==0) b=my-2, j_flag=1;
			else if(j_periodic && j==my-1) b=1, j_flag=1;
			
			if(k_periodic && k==0) c=mz-2, k_flag=1;
			else if(k_periodic && k==mz-1) c=1, k_flag=1;
			
			if(ii_periodic && i==0) a=-2, i_flag=1;
			else if(ii_periodic && i==mx-1) a=mx+1, i_flag=1;
			
			if(jj_periodic && j==0) b=-2, j_flag=1;
			else if(jj_periodic && j==my-1) b=my+1, j_flag=1;
			
			if(kk_periodic && k==0) c=-2, k_flag=1;
			else if(kk_periodic && k==mz-1) c=mz+1, k_flag=1;
			
			if(i_flag) {
				g11[k][j][i] = g11[c][b][a];
				g12[k][j][i] = g12[c][b][a];
				g13[k][j][i] = g13[c][b][a];
			}
			if(j_flag) {
				g21[k][j][i] = g21[c][b][a];
				g22[k][j][i] = g22[c][b][a];
				g23[k][j][i] = g23[c][b][a];
			}
			if(k_flag) {
				g31[k][j][i] = g31[c][b][a];
				g32[k][j][i] = g32[c][b][a];
				g33[k][j][i] = g33[c][b][a];
			}
		}
	}
	
	*/
	
	PetscInt m;
	DMDAVecGetArray(da, user->Gid, &gid);	// for macro gid2()

	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {
		row = lidx2(i, j, k, user);
		if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
			vol[CP] = 1.;	idx[CP] = lidx2(i, j, k, user);
			continue;
		}
		else if ( nvert[k][j][i]>0.1 && row>=0 ) {	// for fsi
			vol[CP] = 1.;	idx[CP] = lidx2(i, j, k, user);
			MatSetValues(user->A, 1, &row, 1, &idx[CP], &vol[CP], INSERT_VALUES);
		}
		else {
			if (nvert[k][j][i] > poisson_threshold) { // i, j, k is not a fluid point
				vol[CP] = 1.; idx[CP] = lidx2(i, j, k, user);
				continue;
			}
			else { // i, j, k is a fluid point
				for(p=-2; p<=2; p++)
				for(q=-2; q<=2; q++)
				for(r=-2; r<=2; r++) SetValue(p,q,r,0);
			  
				for (m=0; m<19; m++) vol[m] = 0.;
			    
				/* Contribution from i+1 - i */
				if (nvert[k][j][i+1] < poisson_threshold && (i != mx-2 || i_periodic ) ) { // i+1, j, k is a fluid point
					/* dpdc{i} = (p_{i+1} - p_{i}) * g11_{i} */
					vol[CP] -= g11[k][j][i]; //i, j, k
					vol[EP] += g11[k][j][i]; // i+1, j, k
					
					AddValue(0,0,0, -g11[k][j][i]);	//i, j, k
					AddValue(1,0,0, g11[k][j][i]);	// i+1, j, k

					/* dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1}) * 0.25 * g12[k][j][i] */
					if ( (j == my-2 && !j_periodic ) || nvert[k][j+1][i] + nvert[k][j+1][i+1] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < poisson_threshold && j!=1 ) {
							vol[CP] += g12[k][j][i] * 0.5; //i, j, k
							vol[EP] += g12[k][j][i] * 0.5; // i+1, j, k
							vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
							vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
						
							AddValue(0,0,0, 0.5*g12[k][j][i]);	//i, j, k
							AddValue(1,0,0, 0.5*g12[k][j][i]);	// i+1, j, k
							AddValue(0,-1,0, -0.5*g12[k][j][i]);	//i, j-1, k
							AddValue(1,-1,0, -0.5*g12[k][j][i]);	//i+1, j-1, k
						}
					}
					else if ( (j == 1 && !j_periodic ) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < poisson_threshold) {
							vol[NP] += g12[k][j][i] * 0.5;  //i, j+1, k
							vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
							vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
							vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
							
							AddValue(0,1,0, 0.5*g12[k][j][i]);	//i, j+1, k
							AddValue(1,1,0, 0.5*g12[k][j][i]);	//i+1, j+1, k
							AddValue(0,0,0, -0.5*g12[k][j][i]);	//i, j, k
							AddValue(1,0,0, -0.5*g12[k][j][i]);	//i+1, j, k
						}
					}
					else {
						vol[NP] += g12[k][j][i] * 0.25; // i, j+1, k
						vol[NE] += g12[k][j][i] * 0.25; // i+1, j+1, k
						vol[SP] -= g12[k][j][i] * 0.25; // i, j-1, k
						vol[SE] -= g12[k][j][i] * 0.25; // i+1, j-1, k
						
						AddValue(0,1,0, 0.25*g12[k][j][i]);	// i, j+1, k
						AddValue(1,1,0, 0.25*g12[k][j][i]);	// i+1, j+1, k
						AddValue(0,-1,0, -0.25*g12[k][j][i]);	// i, j-1, k
						AddValue(1,-1,0, -0.25*g12[k][j][i]);	// i+1, j-1, k
					}
					/* dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1}) * 0.25 * g13[k][j][i] */
					if ( (k == mz-2 && !k_periodic  ) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < poisson_threshold && k!=1) {
							vol[CP] += g13[k][j][i] * 0.5; // i, j, k
							vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
							vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
							vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
						
							AddValue(0,0,0, 0.5*g13[k][j][i]);	// i, j, k
							AddValue(1,0,0, 0.5*g13[k][j][i]);	// i+1, j, k
							AddValue(0,0,-1, -0.5*g13[k][j][i]);	// i, j, k-1
							AddValue(1,0,-1, -0.5*g13[k][j][i]);	// i+1, j, k-1
						}
					}
					else if ( (k == 1 && !k_periodic ) || nvert[k-1][j][i] + nvert[k-1][j][i+1] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < poisson_threshold) {
							vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
							vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
							vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
							vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
							
							AddValue(0,0,1, 0.5*g13[k][j][i]);	// i, j, k+1
							AddValue(1,0,1, 0.5*g13[k][j][i]);	// i+1, j, k+1
							AddValue(0,0,0, -0.5*g13[k][j][i]);	// i, j, k
							AddValue(1,0,0, -0.5*g13[k][j][i]);	// i+1, j, k
						}
					}
					else {
						vol[TP] += g13[k][j][i] * 0.25; //i, j, k+1
						vol[TE] += g13[k][j][i] * 0.25; //i+1, j, k+1
						vol[BP] -= g13[k][j][i] * 0.25; //i, j, k-1
						vol[BE] -= g13[k][j][i] * 0.25; //i+1, j, k-1
						
						AddValue(0,0,1, 0.25*g13[k][j][i]);	 //i, j, k+1
						AddValue(1,0,1, 0.25*g13[k][j][i]);	//i+1, j, k+1
						AddValue(0,0,-1, -0.25*g13[k][j][i]);	//i, j, k-1
						AddValue(1,0,-1, -0.25*g13[k][j][i]);	//i+1, j, k-1
					}
				}  // end of i+1 - i

				/* Contribution from i - i-1 */
				if (nvert[k][j][i-1] < poisson_threshold && (i != 1 || i_periodic ) ) { // i-1, j, k is a fluid point
					/* -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i} */
					vol[CP] -= g11[k][j][i-1];  //i, j, k
					vol[WP] += g11[k][j][i-1];  //i-1, j, k

					/* -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1}) * 0.25 * g12[k][j][i-1] */
					if ( (j == my-2 && !j_periodic ) || nvert[k][j+1][i] + nvert[k][j+1][i-1] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < poisson_threshold && j!=1 ) {
							vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
							vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
							vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
							vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
						}
					}
					else if ( (j == 1 && !j_periodic ) || nvert[k][j-1][i] + nvert[k][j-1][i-1] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < poisson_threshold) {
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

					/* -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1}) * 0.25 * g13[k][j][i] */
					if ( (k == mz-2 && !k_periodic ) || nvert[k+1][j][i] + nvert[k+1][j][i-1] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < poisson_threshold && k!=1 ) {
							vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
							vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
							vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
							vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
						}
					}
					else if ( (k == 1 && !k_periodic ) || nvert[k-1][j][i] + nvert[k-1][j][i-1] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < poisson_threshold) {
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
				if (nvert[k][j+1][i] < poisson_threshold && (j != my-2 || j_periodic ) ) {
					/* dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) * 0.25 */
					if ( (i == mx-2 && !i_periodic ) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] += g21[k][j][i] * 0.5; // i, j, k
							vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
							vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
							vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
						}
					}
					else if ( (i == 1 && !i_periodic ) || nvert[k][j][i-1] + nvert[k][j+1][i-1] > poisson_threshold) {
						if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < poisson_threshold) {
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
					if ( (k == mz-2 && !k_periodic ) || nvert[k+1][j][i] + nvert[k+1][j+1][i] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < poisson_threshold && k!=1 ) {
							vol[CP] += g23[k][j][i] * 0.5; //i,j,k
							vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
							vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
							vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
						}
					}
					else if ( (k == 1 && !k_periodic ) || nvert[k-1][j][i] + nvert[k-1][j+1][i] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < poisson_threshold) {
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
				if (nvert[k][j-1][i] < poisson_threshold && (j!=1 || j_periodic ) ) {
					/* -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) * 0.25 * g21[k][j-1][i] */
					if ( (i == mx-2 && !i_periodic ) || nvert[k][j][i+1] + nvert[k][j-1][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
							vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
							vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
							vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
						}
					}
					else if ( (i == 1 && !i_periodic ) || nvert[k][j][i-1] + nvert[k][j-1][i-1] > poisson_threshold) {
						if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < poisson_threshold) {
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

					/* -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) * 0.25 * g23[k][j-1][i] */
					if ( (k == mz-2 && !k_periodic ) || nvert[k+1][j][i] + nvert[k+1][j-1][i] > poisson_threshold) {
						if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < poisson_threshold && k!=1 ) {
							vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
							vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
							vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
							vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
						}
					}
					else if ( (k == 1 && !k_periodic ) || nvert[k-1][j][i] + nvert[k-1][j-1][i] > poisson_threshold) {
						if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < poisson_threshold) {
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
				if (nvert[k+1][j][i] < poisson_threshold && (k != mz-2 || k_periodic ) ) {
					/* dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) * 0.25 * g31[k][j][i] */
					if ( (i == mx-2 && !i_periodic ) || nvert[k][j][i+1] + nvert[k+1][j][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] += g31[k][j][i] * 0.5;//i, j, k
							vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
							vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
							vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
						}
					}
					else if ( (i == 1 && !i_periodic ) || nvert[k][j][i-1] + nvert[k+1][j][i-1] > poisson_threshold) {
							if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < poisson_threshold) {
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

					/* dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 0.25 * g32[k][j][i] */
					if ( (j == my-2 && !j_periodic ) || nvert[k][j+1][i] + nvert[k+1][j+1][i] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < poisson_threshold && j!=1 ) {
							vol[CP] += g32[k][j][i] * 0.5;//i, j,k
							vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
							vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
							vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
						}
					}
					else if ( (j == 1 && !j_periodic ) || nvert[k][j-1][i] + nvert[k+1][j-1][i] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < poisson_threshold) {
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
				if (nvert[k-1][j][i] < poisson_threshold && (k != 1 || k_periodic ) ) {
					/* -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) * 0.25 * g31[k-1][j][i] */
					if ( (i == mx-2 && !i_periodic ) || nvert[k][j][i+1] + nvert[k-1][j][i+1] > poisson_threshold) {
						if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < poisson_threshold && i!=1 ) {
							vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
							vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
							vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
							vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
						}
					}
					else if ( (i == 1 && !i_periodic ) || nvert[k][j][i-1] + nvert[k-1][j][i-1] > poisson_threshold) {
						if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < poisson_threshold) {
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
			      
					/* -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) *  0.25 * g32[k-1][j][i] */
					// ( p{i, j+1,k-1/2} - p{i, j-1,k-1/2} ) / (2eta)
					if ( (j == my-2 && !j_periodic ) || nvert[k][j+1][i] + nvert[k-1][j+1][i] > poisson_threshold) {
						if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < poisson_threshold && j!=1 ) {
							vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
							vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
							vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
							vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
						}
					}
					else if ( (j == 1 && !j_periodic ) || nvert[k][j-1][i] + nvert[k-1][j-1][i] > poisson_threshold) {
						if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < poisson_threshold) {
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
					//vol[m] *= aj[k][j][i];
				}

			    
				idx[CP] = lidx2(i  , j  , k  , user);
				idx[EP] = lidx2(i+1, j  , k  , user);
				idx[WP] = lidx2(i-1, j  , k  , user);
				idx[NP] = lidx2(i  , j+1, k  , user);
				idx[SP] = lidx2(i  , j-1, k  , user);
				idx[TP] = lidx2(i  , j  , k+1, user);
				idx[BP] = lidx2(i  , j  , k-1, user);
				idx[NE] = lidx2(i+1, j+1, k  , user);
				idx[SE] = lidx2(i+1, j-1, k  , user);
				idx[NW] = lidx2(i-1, j+1, k  , user);
				idx[SW] = lidx2(i-1, j-1, k  , user);
				idx[TN] = lidx2(i  , j+1, k+1, user);
				idx[BN] = lidx2(i  , j+1, k-1, user);
				idx[TS] = lidx2(i  , j-1, k+1, user);
				idx[BS] = lidx2(i  , j-1, k-1, user);
				idx[TE] = lidx2(i+1, j  , k+1, user);
				idx[BE] = lidx2(i+1, j  , k-1, user);
				idx[TW] = lidx2(i-1, j  , k+1, user);
				idx[BW] = lidx2(i-1, j  , k-1, user);
				
				//MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);
				//for(m=0; m<19; m++) if( fabs(vol[m])>1.e-70 || m==CP ) MatSetValues(user->A, 1, &row, 1, &idx[m], &vol[m], INSERT_VALUES);
				for(m=0; m<19; m++) {
					if( fabs(vol[m])>1.e-10 && idx[m]>=0/*|| m==CP */) MatSetValues(user->A, 1, &row, 1, &idx[m], &vol[m], INSERT_VALUES);
				}
				/*
				if( (i==1 && j==2 && k==1) || (i==1 && j==2 && k==2)  || (i==1 && j==2 && k==3)  || (i==1 && j==2 && k==59) || (i==1 && j==2 && k==60) ) {
					printf("%d %d %d, my_idx = %d\n", i, j, k, idx[0]);
					for(m=0; m<19; m++) printf("%d\t", m);
					printf("\n");
					for(m=0; m<19; m++) printf("%d\t", idx[m]);
					printf("\n");
					for(m=0; m<19; m++) printf("%.2f\t", vol[m]);
					printf("\n\n");
				}*/
			} // End of fluid point
		} // End of interial points
	}
    //exit(0);
	DMDAVecRestoreArray(da, user->Gid, &gid);

	//kangsk  
	PetscPrintf(PETSC_COMM_WORLD, "MatAssemblyBegin...\n");
	MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);
	PetscPrintf(PETSC_COMM_WORLD, "MatAssemblyEnd...\n");

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
	DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);
	
	return 0;
}


PetscErrorCode PoissonRHS2(UserCtx *user, Vec B)
{
	DMDALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
	PetscReal	***nvert, ***aj, ***gid, dt = user->dt;
	Cmpnts	***ucont, ***lucont;
	PetscReal poisson_threshold=0.1;

	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  	
	DMDAVecGetArray(user->fda, user->lUcont, &ucont);
	DMDAVecGetArray(user->da, user->lNvert, &nvert);
	DMDAVecGetArray(user->da, user->lAj, &aj);
  
	DMDAVecGetArray(user->da, user->Gid, &gid);
	
	int lcount=0;
		
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		double val;
		if (nvert[k][j][i] >= poisson_threshold) {
			if( (int) (gid[k][j][i]) >=0 ) val = 0;	// for fsi
			else continue;
		}
		else {
			double coeff=time_coeff();
			#ifdef DIRICHLET
			//if(freesurface && j==user->free_surface_j[i][k]) coeff *= user->vof[i][k];
			#endif

			val=0;
			
			val -= ucont[k][j][i].x;

			if(i==1 && i_periodic) val += ucont[k][j][mx-2].x;
			//	else if(i==1 && ii_periodic)  val += ucont[k][j][-2].x;
			else val += ucont[k][j][i-1].x;

			val -= ucont[k][j][i].y;

			if(j==1 && j_periodic) val += ucont[k][my-2][i].y;
			//else if(j==1 && jj_periodic) val += ucont[k][-2][i].y;
			else val += ucont[k][j-1][i].y;

			val -= ucont[k][j][i].z;

			if(k==1 && k_periodic) val += ucont[mz-2][j][i].z;
			//else if(k==1 && kk_periodic) val += ucont[-2][j][i].z;
			else val += ucont[k-1][j][i].z;
		/*
			if(  (i==1 && j==2) ) {
				printf("%d %d %d, my_idx = %d, val=%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", i, j, k, 
					(int)gid[k][j][i], ucont[k][j][i-1].x, ucont[k][j][i].x, ucont[k][j-1][i].y, ucont[k][j][i].y, ucont[k-1][j][i].z, ucont[k][j][i].z);
				printf("\n");
			}*/
		
			val *=  -1.0 / dt * user->st * coeff;
				
			lcount++;
		}
		VecSetValue(B, (int)gid[k][j][i], val, INSERT_VALUES);
		
		
	}
  //exit(0);
	VecAssemblyBegin(B);
	VecAssemblyEnd(B);

	
  
#ifndef DIRICHLET  
	
	
	double sum, sum1;
	VecSum(B, &sum);
	/*
	int N;
	VecGetSize(B, &N);
	sum  = sum/(-1.0*N);
	VecShift(B,sum);
	user->multinullspace = PETSC_FALSE;
	*/
	MPI_Allreduce( &lcount, &user->rhs_count, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
	
	double val = -sum/(double) (user->rhs_count);	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]<0.1) VecSetValue(B, (int)gid[k][j][i], val, ADD_VALUES);
	}
	
	VecAssemblyBegin(B);
	VecAssemblyEnd(B);
	VecSum(B, &sum1);
#endif 
  
	PetscPrintf(PETSC_COMM_WORLD, "Summation RHS %e %e\n", sum, sum1); 
	
	DMDAVecRestoreArray(user->da, user->Gid, &gid);
	DMDAVecRestoreArray(user->fda, user->lUcont, &ucont);
	DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
	DMDAVecRestoreArray(user->da, user->lAj, &aj);
  

  return 0;
}


void Create_Hypre_Solver(PetscInt bi)
{
	/*
	0  CLJP-coarsening (a parallel coarsening algorithm using independent sets.
	1  classical Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
	3  classical Ruge-Stueben coarsening on each processor, followed by a third pass, which adds coarse
	   points on the boundaries
	6  Falgout coarsening (uses 1 first, followed by CLJP using the interior coarse points
	   generated by 1 as its first independent set)
	7  CLJP-coarsening (using a fixed random vector, for debugging purposes only)
	8  PMIS-coarsening (a parallel coarsening algorithm using independent sets, generating
	   lower complexities than CLJP, might also lead to slower convergence)
	9  PMIS-coarsening (using a fixed random vector, for debugging purposes only)
	10 HMIS-coarsening (uses one pass Ruge-Stueben on each processor independently, followed
	   by PMIS using the interior C-points generated as its first independent set)
	11 one-pass Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
	21 CGC coarsening by M. Griebel, B. Metsch and A. Schweitzer
	22 CGC-E coarsening by M. Griebel, B. Metsch and A.Schweitzer
	*/
	
	
	HYPRE_BoomerAMGCreate (&precon_p[bi]);
	HYPRE_BoomerAMGSetInterpType(precon_p[bi], 6);	// 6:ext, 13: FF1
	HYPRE_BoomerAMGSetAggNumLevels(precon_p[bi], 0);	// FF1+aggresive coarsening is good for > 50mil grids, 0 for better convergence
	//HYPRE_BoomerAMGSetTol (precon_p, 1.e-6);
	HYPRE_BoomerAMGSetTol (precon_p[bi], poisson_tol);
	HYPRE_BoomerAMGSetPrintLevel(precon_p[bi],1);
	HYPRE_BoomerAMGSetCoarsenType(precon_p[bi],8);	// 0:CLJP, 6:Falgout, 8:PMIS, 10:HMIS, 
	HYPRE_BoomerAMGSetCycleType(precon_p[bi],1);// 1
	HYPRE_BoomerAMGSetMaxIter(precon_p[bi],1);
	HYPRE_BoomerAMGSetStrongThreshold(precon_p[bi],0.75);// 0.5 : Cartesian, 0.6 : Distorted
	
	//HYPRE_BoomerAMGSetRelaxType (precon_p, 6);
	//HYPRE_BoomerAMGSetRelaxWt(precon_p,0.8);
	
	//HYPRE_BoomerAMGSetCycleNumSweeps(precon_p,5,3);	// 5 sweep at coarsest level
	//HYPRE_BoomerAMGSetLevelOuterWt(precon_p,0.5,0);
	//HYPRE_BoomerAMGSetSmoothType(precon_p, 7);	// more complex smoother
	//HYPRE_BoomerAMGSetSmoothNumSweeps(precon_p,2);
	
	// Pressure solver
	
	/*
	HYPRE_ParCSRPCGCreate ( PETSC_COMM_WORLD, &pcg_solver_p );
	HYPRE_PCGSetPrecond ( pcg_solver_p, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precon_p );
	HYPRE_PCGSetMaxIter( pcg_solver_p, 50 );
	HYPRE_PCGSetTol( pcg_solver_p, 1.e-09 );
	HYPRE_PCGSetPrintLevel( pcg_solver_p, 3 ); 
	*/
	
	HYPRE_ParCSRGMRESCreate ( PETSC_COMM_WORLD, &pcg_solver_p[bi] );
	HYPRE_ParCSRGMRESSetKDim (pcg_solver_p[bi], 51);
	//HYPRE_ParCSRGMRESSetAbsoluteTol (pcg_solver_p, 1.e-12);
	HYPRE_GMRESSetPrecond ( pcg_solver_p[bi], (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precon_p[bi] );
	PetscInt poisson_it=10;
	PetscOptionsGetInt(PETSC_NULL, "-poisson_it", &poisson_it, PETSC_NULL);
	HYPRE_GMRESSetMaxIter( pcg_solver_p[bi], poisson_it);
	HYPRE_GMRESSetTol( pcg_solver_p[bi], poisson_tol );
	HYPRE_GMRESSetPrintLevel( pcg_solver_p[bi], 3 ); 
};

void myMatHYPRE_IJMatrixCopy(Mat v,HYPRE_IJMatrix ij)
 {
	int i, rstart,rend;//
	const int *cols;
	int ncols;
	const PetscScalar *values;

	//HYPRE_IJMatrixInitialize(ij);	//->> Call just once when initialize matrix. memory leaks
	MatGetOwnershipRange(v,&rstart,&rend);

	for (i=rstart; i<rend; i++) {
		 MatGetRow(v,i,&ncols,&cols,&values);
		 HYPRE_IJMatrixSetValues(ij,1,&ncols,&i,cols,values);
		 MatRestoreRow(v,i,&ncols,&cols,&values);
	}

	HYPRE_IJMatrixAssemble(ij);
 };
 
void Petsc_to_Hypre_Vector(Vec A, HYPRE_IJVector B, int i_lower)
{
  int       localsize, *indices, i;
  PetscReal *a;
  double    *values;
	
  VecGetLocalSize(A, &localsize);
  
  PetscMalloc(localsize*sizeof(double), &values);
  PetscMalloc(localsize*sizeof(int), &indices);
  //	std::vector<double> values (localsize);
  //	std::vector<int> indices (localsize);
  
  //HYPRE_IJVectorInitialize(B);	//->> Call just once when initialize vector. memory leaks
  
  VecGetArray(A, &a);
  for(i=0; i<localsize; i++) 
    values[i] = a[i], indices[i] = i_lower + i;

  VecRestoreArray(A, &a);
  
  HYPRE_IJVectorSetValues(B, localsize, &indices[0], &values[0]);
  HYPRE_IJVectorAssemble(B);

  PetscFree(indices);PetscFree(values);
  //HYPRE_IJVectorGetObject(B, (void **) &par_B);
};

void Hypre_to_Petsc_Vector(HYPRE_IJVector B, Vec A, int i_lower)
{
  int       localsize, *indices, i;
  double    *values;

  VecGetLocalSize(A, &localsize);
  
  //	std::vector<double> values (localsize);
  //	std::vector<int> indices (localsize);	
  PetscMalloc(localsize*sizeof(double), &values);
  PetscMalloc(localsize*sizeof(int), &indices);
  
  for(i=0; i<localsize; i++) indices[i] = i_lower + i;
  
  HYPRE_IJVectorGetValues(B, localsize, &indices[0], &values[0]);
  /*
    VecSetValues(A, localsize, &indices[0], &values[0], INSERT_VALUES);
    VecAssemblyBegin(A);
    VecAssemblyEnd(A);*/
  
  PetscReal *a;
  VecGetArray(A, &a);
  for(i=0; i<localsize; i++) a[i] = values[i];
  VecRestoreArray(A, &a);

  PetscFree(indices);PetscFree(values);
};


void Create_Hypre_Matrix(UserCtx *user)
{
  int localsize_p, *nz_p,i;
  VecGetLocalSize(user->Phi2, &localsize_p);

  PetscMalloc(localsize_p*sizeof(int), &nz_p);
	
  int p_lower = user->p_global_begin;
  int p_upper = p_lower + localsize_p - 1;
	
  //  std::vector<int> nz_p (localsize_p);
  //std::fill ( nz_p.begin(), nz_p.end(), 19);
  for(i=0; i<localsize_p; i++) nz_p[i]=19;
  
  PetscPrintf(PETSC_COMM_WORLD, "\nbegin HYPRE_IJMatrixCreate\n");
  HYPRE_IJMatrixCreate(PETSC_COMM_WORLD, p_lower, p_upper, p_lower, p_upper, &Ap[user->_this]);
  HYPRE_IJMatrixSetObjectType(Ap[user->_this], HYPRE_PARCSR);
  HYPRE_IJMatrixSetRowSizes(Ap[user->_this], &nz_p[0]);
  HYPRE_IJMatrixSetMaxOffProcElmts (Ap[user->_this], 10);
  HYPRE_IJMatrixInitialize(Ap[user->_this]);
  PetscPrintf(PETSC_COMM_WORLD, "end HYPRE_IJMatrixCreate\n\n");
  PetscFree(nz_p);
}

void Create_Hypre_Vector(UserCtx *user)
{
	int localsize_p;
	
	VecGetLocalSize(user->Phi2, &localsize_p);
	
	int p_lower = user->p_global_begin;
	int p_upper = p_lower + localsize_p - 1;
	
	PetscPrintf(PETSC_COMM_WORLD, "\nbegin HYPRE_IJVectorCreate\n");
	
	// p vector
	HYPRE_IJVectorCreate(PETSC_COMM_WORLD, p_lower, p_upper, &Vec_p[user->_this]);
	HYPRE_IJVectorSetObjectType(Vec_p[user->_this], HYPRE_PARCSR);
	HYPRE_IJVectorSetMaxOffProcElmts (Vec_p[user->_this], 10);
		
	HYPRE_IJVectorCreate(PETSC_COMM_WORLD, p_lower, p_upper, &Vec_p_rhs[user->_this]);
	HYPRE_IJVectorSetObjectType(Vec_p_rhs[user->_this], HYPRE_PARCSR);
	HYPRE_IJVectorSetMaxOffProcElmts (Vec_p_rhs[user->_this], 10);
	
	HYPRE_IJVectorInitialize(Vec_p[user->_this]);
	HYPRE_IJVectorInitialize(Vec_p_rhs[user->_this]);
	
	PetscPrintf(PETSC_COMM_WORLD, "end HYPRE_IJVectorCreate\n\n");
};


void PoissonSolver_Hypre(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo)
{
	DMDALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;
	
	int lxs = xs, lxe = xe;
	int lys = ys, lye = ye;
	int lzs = zs, lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
	
	PetscInt l, i,j,k;
	const PetscInt bi=0;

	PetscReal ts,te,cput;
	PetscTime(&ts);
	
	PetscInt	m_c, m_f, M_c, M_f;	

	PetscOptionsHasName(PETSC_NULL,"-hypre_symmetric",&hypre_symmetric);
	PetscPrintf(PETSC_COMM_WORLD, "SYMMETRIC HYPRE!!!!\n");

	if(ti==tistart) {
		PoissonLHS_HYPRE(user);
		// 07.09.2009 Seokkoo
		
		Create_Hypre_Matrix(user);
		Create_Hypre_Vector(user);
		
		myMatHYPRE_IJMatrixCopy(user->A, Ap[user->_this]);
		HYPRE_IJMatrixGetObject(Ap[user->_this], (void**) &par_Ap[user->_this]);
		MatDestroy(&user->A);	user->assignedA=PETSC_FALSE;
		Create_Hypre_Solver(user->_this);
	}
	else if(!hypre_symmetric) {
		VecDestroy(&user->Gid);
		PoissonLHS_HYPRE(user);
		
		HYPRE_IJMatrixDestroy(Ap[user->_this]);
		HYPRE_IJVectorDestroy(Vec_p[user->_this]);
		HYPRE_IJVectorDestroy(Vec_p_rhs[user->_this]);

		Create_Hypre_Vector(user);
		Create_Hypre_Matrix(user);
		
		myMatHYPRE_IJMatrixCopy(user->A, Ap[user->_this]);
		HYPRE_IJMatrixGetObject(Ap[user->_this], (void**) &par_Ap[user->_this]);
		MatDestroy(&user->A);	user->assignedA=PETSC_FALSE;
	} 
	
	VecDuplicate(user->Phi2, &user->B2);
	PetscReal ibm_Flux, ibm_Area;

	
	VolumeFlux(&user[bi], &ibm_Flux, &ibm_Area, 0);

	PoissonRHS2(user, user->B2);
	
	// 07.09.2009 Seokkoo
	Petsc_to_Hypre_Vector(user->B2, Vec_p_rhs[user->_this], user->p_global_begin);
	Petsc_to_Hypre_Vector(user->Phi2, Vec_p[user->_this], user->p_global_begin);
	
	HYPRE_IJVectorGetObject(Vec_p[user->_this], (void **) &par_Vec_p[user->_this]);
	HYPRE_IJVectorGetObject(Vec_p_rhs[user->_this], (void **) &par_Vec_p_rhs[user->_this]);
	if(ti==tistart) {
		//HYPRE_ParCSRPCGSetup (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
		HYPRE_ParCSRGMRESSetup (pcg_solver_p[user->_this], par_Ap[user->_this], par_Vec_p_rhs[user->_this], par_Vec_p[user->_this]);
	}
	else if(!hypre_symmetric) {
		HYPRE_ParCSRGMRESSetup (pcg_solver_p[user->_this], par_Ap[user->_this], par_Vec_p_rhs[user->_this], par_Vec_p[user->_this]);
	}
	
	PetscPrintf(PETSC_COMM_WORLD, "Solving Poisson ...\n");
	
	MPI_Barrier(PETSC_COMM_WORLD);
	
	//HYPRE_ParCSRPCGSolve (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
	HYPRE_ParCSRGMRESSolve (pcg_solver_p[user->_this], par_Ap[user->_this], par_Vec_p_rhs[user->_this], par_Vec_p[user->_this]);
	Hypre_to_Petsc_Vector(Vec_p[user->_this], user->Phi2, user->p_global_begin);
	
	MPI_Barrier(PETSC_COMM_WORLD);
	#ifndef DIRICHLET
	//int N;
	double sum;
	//VecGetSize(user->Phi2,&N);
	VecSum(user->Phi2,&sum);
	//sum  = sum/(-1.0*N);
	//VecShift(user->Phi2,sum);
	double val = -sum/(double) (user->rhs_count);
	
	PetscReal ***gid, ***nvert;
	
	DMDAVecGetArray(user->da, user->Gid, &gid);
	DMDAVecGetArray(user->da, user->lNvert, &nvert);
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]<0.1) {
			VecSetValue(user->Phi2, (int)gid[k][j][i], val, ADD_VALUES);
		}
	}
	DMDAVecRestoreArray(user->da, user->Gid, &gid);
	DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
	
	VecAssemblyBegin(user->Phi2);
	VecAssemblyEnd(user->Phi2);
	#endif
	
	Convert_Phi2_Phi(user);
	/*
	DMDAVecGetArray(user->da, user->Phi, &gid);
	for (int k=lzs; k<lze; k++)
	for (int j=lys; j<lye; j++)
	for (int i=lxs; i<lxe; i++) {
		if(i==1 && j==2) printf("%d,%d,%d   %.5f\n", i,j,k,gid[k][j][i]);
	}
	DMDAVecRestoreArray(user->da, user->Phi, &gid);
	exit(0);
	*/
		
	DMGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi);
	DMGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi);
	
	VecDestroy(&user->B2);
	
	PetscTime(&te);
	cput=te-ts;
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (!rank) {
		FILE *f;
		char filen[80];
		sprintf(filen, "Converge_dU");
		f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, " %d(poisson)  %.2e(s)", ti, cput);
		fclose(f);
	}
	

}

