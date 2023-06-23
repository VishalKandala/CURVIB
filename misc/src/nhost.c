#include "petscvec.h"
#include "petscda.h"
#include "petscksp.h"
static char help[] = "Testing programming!";

typedef struct UserCtx {
  DA da;	/* Data structure for scalars (include the grid geometry
		   informaion, to obtain the grid information, 
		   use DAGetCoordinates) */
  DA fda;	// Data Structure for vectors
  DALocalInfo info;
  Vec         nhostU, Ucat;
  PetscInt	IM, JM, KM; // dimensions of grid
} UserCtx;

PetscInt block_number=2;


#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
  UserCtx	*user;
  PetscInt      bi,ti,tistart=0;
  Vec	hostU;
  VecScatter tolocalall;
  PetscInt tisteps = 1000000;
  PetscTruth	flg;
  PetscErrorCode ierr;
  DAPeriodicType  wrap=DA_NONPERIODIC;
  PetscInt m, n, p, MM,NN,PP;

  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-totalsteps", &tisteps, &flg);

  PetscMalloc(block_number*sizeof(UserCtx), &user);

  CHKMEMQ
  for (bi=0; bi<block_number; bi++) {
    if (bi==0) {
      user[bi].IM=181;
      user[bi].JM=191;
      user[bi].KM=173;
    } else {
      user[bi].IM=141;
      user[bi].JM=141;
      user[bi].KM=133;
    } 

  CHKMEMQ
    ierr=DACreate3d(PETSC_COMM_WORLD, wrap, DA_STENCIL_BOX,
	       user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, PETSC_DECIDE, 
	       PETSC_DECIDE,
	       PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
	       &(user[bi].da));   
    if (ierr) SETERRQ1(1, "problem creating DA %d",ierr);

  CHKMEMQ
    DASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DAGetCoordinateDA(user[bi].da, &(user[bi].fda));
    DAGetLocalInfo(user[bi].da, &(user[bi].info));
    DAGetInfo(user[bi].da, PETSC_NULL, &MM,&NN,&PP,
	      &m, &n, &p, PETSC_NULL, PETSC_NULL,
	      PETSC_NULL, PETSC_NULL);

    PetscPrintf(PETSC_COMM_WORLD, "DA Distribution: %i %i %i for block %i\n", m, n, p,bi);

  CHKMEMQ    
    ierr = DACreateGlobalVector(user[bi].fda, &(user[bi].Ucat));
    VecSet(user[bi].Ucat,0.);
  CHKMEMQ
  }

  for (ti = tistart; ti<tistart + tisteps; ti++) {
    PetscPrintf(PETSC_COMM_WORLD, "Time %d\n", ti);
  CHKMEMQ
    for (bi=0; bi<block_number; bi++) {
      /* hostU is a parallel PETSc vector that will hold vector values
	 in the natural numbering, rather than in the PETSc parallel
	 numbering associated with the DA */
  CHKMEMQ
      DACreateNaturalVector(user[bi].fda, &hostU);
      
  CHKMEMQ
      // put the center node velocities in the natural ordering in hostU
      DAGlobalToNaturalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);
      DAGlobalToNaturalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);
      
  CHKMEMQ
      /* the velocities at cell centers are saved in user[bi].nhostU
	 which is a sequential vector*/
      VecScatterCreateToAll(hostU, &tolocalall, &(user[bi].nhostU));
  CHKMEMQ      
      VecScatterBegin(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		      SCATTER_FORWARD);
      VecScatterEnd(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		    SCATTER_FORWARD);
  CHKMEMQ
      
      VecScatterDestroy(tolocalall);
      VecDestroy(hostU);
    }
  CHKMEMQ

    for (bi=0; bi<block_number; bi++) {
      VecDestroy(user[bi].nhostU);
    }
  }
  PetscFinalize();
}
