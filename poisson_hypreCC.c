#include "variables.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscmg.h"
#include <stdlib.h>
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_IJ_mv.h"

extern int i_periodic, j_periodic, k_periodic, pseudo_periodic;
extern PetscInt block_number, freesurface, immersed, ti, tistart, les, tiout;
extern PetscInt tistart;
extern double poisson_tol;
//extern char path[256];
extern double mean_pressure_gradient;
extern PetscInt movefsi, rotatefsi;

PetscErrorCode PoissonLHSNew(UserCtx *user);
PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg);

HYPRE_Solver pcg_solver_p, precon_p;
HYPRE_IJMatrix Ap;
HYPRE_ParCSRMatrix par_Ap;
HYPRE_IJVector Vec_p, Vec_p_rhs;
HYPRE_ParVector par_Vec_p, par_Vec_p_rhs;


double time_coeff()
{
	if(ti<=tistart+1) return 1;
	else if(ti<=tistart+2) return 1.5;
	else return 11./6.;		// 1.8333
	
};

void Convert_Phi2_Phi(UserCtx *user)
{
	DALocalInfo	info = user->info;
	DA		da = user->da;
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
	DAVecGetArray(da, user->Phi, &phi);
	DAVecGetArray(da, user->lNvert, &nvert);
	VecGetArray(user->Phi2, &phi2);
	
	int pos=0;
	for(k=lzs; k<lze; k++)
	for(j=lys; j<lye; j++)
	for(i=lxs; i<lxe; i++) {
		/*if(i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
			phi[k][j][i]=0;
		}
		else */if( (int)nvert[k][j][i]>poisson_threshold ) {
			if(movefsi || rotatefsi) {
				phi[k][j][i] = phi2[pos++];
			}
			else phi[k][j][i]=0;
		}
		else {
			phi[k][j][i] = phi2[pos++];
		}
	}
	
	DAVecRestoreArray(da, user->Phi, &phi);
	DAVecRestoreArray(da, user->lNvert, &nvert);
	VecRestoreArray(user->Phi2, &phi2);
	
	DAGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
	DAGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);
	
	DAVecGetArray(da, user->lPhi, &lphi);
	DAVecGetArray(da, user->Phi, &phi);
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
		
		if(ii_periodic && i==0) a=-2, flag=1;
		else if(ii_periodic && i==mx-1) a=mx+1, flag=1;
		
		if(jj_periodic && j==0) b=-2, flag=1;
		else if(jj_periodic && j==my-1) b=my+1, flag=1;
		
		if(kk_periodic && k==0) c=-2, flag=1;
		else if(kk_periodic && k==mz-1) c=mz+1, flag=1;
		
		if(flag) {
			phi[k][j][i] = lphi[c][b][a];
		}
		
		//if(k==1 && i!=0)	printf("%d,%d,%d, %f %f %f\n", i,j,k, lphi[-2][j][i], lphi[k][j][i], lphi[k+1][j][i]);
	}
	DAVecRestoreArray(da, user->lPhi, &lphi);
	DAVecRestoreArray(da, user->Phi, &phi);
	//exit(0);
}

PetscErrorCode PoissonRHS2(UserCtx *user, Vec B)
{
	DALocalInfo info = user->info;
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
  	
	DAVecGetArray(user->fda, user->lUcont, &ucont);
	DAVecGetArray(user->da, user->lNvert, &nvert);
	DAVecGetArray(user->da, user->lAj, &aj);
  
	DAVecGetArray(user->da, user->Gid, &gid);
	
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
			else if(i==1 && ii_periodic)  val += ucont[k][j][-2].x;
			else val += ucont[k][j][i-1].x;

			val -= ucont[k][j][i].y;

			if(j==1 && j_periodic) val += ucont[k][my-2][i].y;
			else if(j==1 && jj_periodic) val += ucont[k][-2][i].y;
			else val += ucont[k][j-1][i].y;

			val -= ucont[k][j][i].z;

			if(k==1 && k_periodic) val += ucont[mz-2][j][i].z;
			else if(k==1 && kk_periodic) val += ucont[-2][j][i].z;
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
	
	DAVecRestoreArray(user->da, user->Gid, &gid);
	DAVecRestoreArray(user->fda, user->lUcont, &ucont);
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	DAVecRestoreArray(user->da, user->lAj, &aj);
  

  return 0;
}


void Create_Hypre_Solver()
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
	
	
	HYPRE_BoomerAMGCreate (&precon_p);
	HYPRE_BoomerAMGSetInterpType(precon_p, 6);	// 6:ext, 13: FF1
	HYPRE_BoomerAMGSetAggNumLevels(precon_p, 1);	// FF1+aggresive coarsening is good for > 50mil grids
	//HYPRE_BoomerAMGSetTol (precon_p, 1.e-6);
	HYPRE_BoomerAMGSetTol (precon_p, poisson_tol);
	HYPRE_BoomerAMGSetPrintLevel(precon_p,1);
	HYPRE_BoomerAMGSetCoarsenType(precon_p,8);	// 0:CLJP, 6:Falgout, 8:PMIS, 10:HMIS, 
	HYPRE_BoomerAMGSetCycleType(precon_p,1);// 1
	HYPRE_BoomerAMGSetMaxIter(precon_p,1);
	HYPRE_BoomerAMGSetStrongThreshold(precon_p,0.75);// 0.5 : Cartesian, 0.6 : Distorted
	
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
	
	HYPRE_ParCSRGMRESCreate ( PETSC_COMM_WORLD, &pcg_solver_p );
	HYPRE_ParCSRGMRESSetKDim (pcg_solver_p, 51);
	//HYPRE_ParCSRGMRESSetAbsoluteTol (pcg_solver_p, 1.e-12);
	HYPRE_GMRESSetPrecond ( pcg_solver_p, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precon_p );
	PetscInt poisson_it;
	PetscOptionsGetInt(PETSC_NULL, "-poisson_it", &poisson_it, PETSC_NULL);
	HYPRE_GMRESSetMaxIter( pcg_solver_p, poisson_it);
	HYPRE_GMRESSetTol( pcg_solver_p, poisson_tol );
	HYPRE_GMRESSetPrintLevel( pcg_solver_p, 3 ); 
};

void MatHYPRE_IJMatrixCopy(Mat v,HYPRE_IJMatrix &ij)
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
 
void Petsc_to_Hypre_Vector(Vec A, HYPRE_IJVector &B, int i_lower)
{
	int localsize;
	PetscReal *a;
	
	VecGetLocalSize(A, &localsize);

	std::vector<double> values (localsize);
	std::vector<int> indices (localsize);
	
	//HYPRE_IJVectorInitialize(B);	//->> Call just once when initialize vector. memory leaks

	VecGetArray(A, &a);
	for(register int i=0; i<localsize; i++) values[i] = a[i], indices[i] = i_lower + i;
	VecRestoreArray(A, &a);
	
	HYPRE_IJVectorSetValues(B, localsize, &indices[0], &values[0]);
	HYPRE_IJVectorAssemble(B);
	//HYPRE_IJVectorGetObject(B, (void **) &par_B);
};

void Hypre_to_Petsc_Vector(HYPRE_IJVector &B, Vec A, int i_lower)
{
	int localsize;
	VecGetLocalSize(A, &localsize);
	
	std::vector<double> values (localsize);
	std::vector<int> indices (localsize);
	
	for(register int i=0; i<localsize; i++) indices[i] = i_lower + i;
	HYPRE_IJVectorGetValues(B, localsize, &indices[0], &values[0]);
	/*
	VecSetValues(A, localsize, &indices[0], &values[0], INSERT_VALUES);
	VecAssemblyBegin(A);
	VecAssemblyEnd(A);*/
	
	PetscReal *a;
	VecGetArray(A, &a);
	for(register int i=0; i<localsize; i++) a[i] = values[i];
	VecRestoreArray(A, &a);
};


void Create_Hypre_Matrix(UserCtx *user)
{
	int localsize_p;
	VecGetLocalSize(user->Phi2, &localsize_p);
	
	int p_lower = user->p_global_begin;
	int p_upper = p_lower + localsize_p - 1;
	
	std::vector<int> nz_p (localsize_p);
	std::fill ( nz_p.begin(), nz_p.end(), 19);

	PetscPrintf(PETSC_COMM_WORLD, "\nbegin HYPRE_IJMatrixCreate\n");
	HYPRE_IJMatrixCreate(PETSC_COMM_WORLD, p_lower, p_upper, p_lower, p_upper, &Ap);
	HYPRE_IJMatrixSetObjectType(Ap, HYPRE_PARCSR);
	HYPRE_IJMatrixSetRowSizes(Ap, &nz_p[0]);
	HYPRE_IJMatrixSetMaxOffProcElmts (Ap, 10);
	HYPRE_IJMatrixInitialize(Ap);
	PetscPrintf(PETSC_COMM_WORLD, "end HYPRE_IJMatrixCreate\n\n");
}

void Create_Hypre_Vector(UserCtx *user)
{
	int localsize_p;
	
	VecGetLocalSize(user->Phi2, &localsize_p);
	
	int p_lower = user->p_global_begin;
	int p_upper = p_lower + localsize_p - 1;
	
	PetscPrintf(PETSC_COMM_WORLD, "\nbegin HYPRE_IJVectorCreate\n");
	
	// p vector
	HYPRE_IJVectorCreate(PETSC_COMM_WORLD, p_lower, p_upper, &Vec_p);
	HYPRE_IJVectorSetObjectType(Vec_p, HYPRE_PARCSR);
	HYPRE_IJVectorSetMaxOffProcElmts (Vec_p, 10);
		
	HYPRE_IJVectorCreate(PETSC_COMM_WORLD, p_lower, p_upper, &Vec_p_rhs);
	HYPRE_IJVectorSetObjectType(Vec_p_rhs, HYPRE_PARCSR);
	HYPRE_IJVectorSetMaxOffProcElmts (Vec_p_rhs, 10);
	
	HYPRE_IJVectorInitialize(Vec_p);
	HYPRE_IJVectorInitialize(Vec_p_rhs);
	
	PetscPrintf(PETSC_COMM_WORLD, "end HYPRE_IJVectorCreate\n\n");
};


void PoissonSolver_Hypre(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo)
{
	DALocalInfo info = user->info;
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
	
	PetscInt l;
	const PetscInt bi=0;

	PetscReal ts,te,cput;
	PetscGetTime(&ts);
	
	PetscInt	m_c, m_f, M_c, M_f;	

	if(ti==tistart) {
		PoissonLHSNew(user);
		// 07.09.2009 Seokkoo
		
		Create_Hypre_Matrix(user);
		Create_Hypre_Vector(user);
		
		MatHYPRE_IJMatrixCopy(user->A, Ap);
		HYPRE_IJMatrixGetObject(Ap, (void**) &par_Ap);
		MatDestroy(user->A);	user->assignedA=PETSC_FALSE;
		Create_Hypre_Solver();
	}
	else if(movefsi || rotatefsi) {
		VecDestroy(user->Gid);
		PoissonLHSNew(user);
		
		HYPRE_IJMatrixDestroy(Ap);
		HYPRE_IJVectorDestroy(Vec_p);
		HYPRE_IJVectorDestroy(Vec_p_rhs);

		Create_Hypre_Vector(user);
		Create_Hypre_Matrix(user);
		
		MatHYPRE_IJMatrixCopy(user->A, Ap);
		HYPRE_IJMatrixGetObject(Ap, (void**) &par_Ap);
		MatDestroy(user->A);	user->assignedA=PETSC_FALSE;
	} 
	
	VecDuplicate(user->Phi2, &user->B2);
	PetscReal ibm_Flux, ibm_Area;

	
	VolumeFlux(&user[bi], &ibm_Flux, &ibm_Area, 0);

	PoissonRHS2(user, user->B2);
	
	// 07.09.2009 Seokkoo
	Petsc_to_Hypre_Vector(user->B2, Vec_p_rhs, user->p_global_begin);
	Petsc_to_Hypre_Vector(user->Phi2, Vec_p, user->p_global_begin);
	
	HYPRE_IJVectorGetObject(Vec_p, (void **) &par_Vec_p);
	HYPRE_IJVectorGetObject(Vec_p_rhs, (void **) &par_Vec_p_rhs);
	if(ti==tistart) {
		//HYPRE_ParCSRPCGSetup (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
		HYPRE_ParCSRGMRESSetup (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
	}
	else if(movefsi || rotatefsi) {
		HYPRE_ParCSRGMRESSetup (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
	}
	
	PetscPrintf(PETSC_COMM_WORLD, "Solving Poisson ...\n");
	
	MPI_Barrier(PETSC_COMM_WORLD);
	
	//HYPRE_ParCSRPCGSolve (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
	HYPRE_ParCSRGMRESSolve (pcg_solver_p, par_Ap, par_Vec_p_rhs, par_Vec_p);
	Hypre_to_Petsc_Vector(Vec_p, user->Phi2, user->p_global_begin);
	
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
	
	DAVecGetArray(user->da, user->Gid, &gid);
	DAVecGetArray(user->da, user->lNvert, &nvert);
	
	for (int k=lzs; k<lze; k++)
	for (int j=lys; j<lye; j++)
	for (int i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]<0.1) {
			VecSetValue(user->Phi2, (int)gid[k][j][i], val, ADD_VALUES);
		}
	}
	DAVecRestoreArray(user->da, user->Gid, &gid);
	DAVecRestoreArray(user->da, user->lNvert, &nvert);
	
	VecAssemblyBegin(user->Phi2);
	VecAssemblyEnd(user->Phi2);
	#endif
	
	Convert_Phi2_Phi(user);
	/*
	DAVecGetArray(user->da, user->Phi, &gid);
	for (int k=lzs; k<lze; k++)
	for (int j=lys; j<lye; j++)
	for (int i=lxs; i<lxe; i++) {
		if(i==1 && j==2) printf("%d,%d,%d   %.5f\n", i,j,k,gid[k][j][i]);
	}
	DAVecRestoreArray(user->da, user->Phi, &gid);
	exit(0);
	*/
		
	DAGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi);
	DAGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi);
	
	VecDestroy(user->B2);
	
	PetscGetTime(&te);
	cput=te-ts;
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (!rank) {
		FILE *f;
		char filen[80];
		sprintf(filen, "Converge_dU");
		f = fopen(filen, "a");
		PetscFPrintf(PETSC_COMM_WORLD, f, "%d(poisson)  %.2e(s)", ti, cput);
		fclose(f);
	}
	

}

