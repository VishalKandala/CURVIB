static char help[] = "Testing programming!";

#include "variables.h"
#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

PetscInt NumberOfBodies=1;

PetscErrorCode ibm_read_tecplot(IBMNodes *ibm, PetscInt nt, PetscTruth flg)
{
PetscInt	rank;
PetscInt	n_v , n_elmt ;
PetscReal	*x_bp , *y_bp , *z_bp ;
PetscInt	*nv1 , *nv2 , *nv3 ;
PetscReal	*nf_x, *nf_y, *nf_z;
PetscInt	i;
PetscInt	n1e, n2e, n3e;
PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
PetscReal	t, dr;
PetscInt 	temp;
double          xt;
char            string[128];
char            str[20];

  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;

PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT print out \n ");

  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


 FILE *fd;
 char filen[128];
 sprintf(filen, "stress%1.1du.dat", nt);
 
 fd = fopen(filen, "r");
 if (!fd) 
   {
     PetscPrintf(PETSC_COMM_WORLD, "Cannot open IBM Solid node file %d\n", nt);
     SETERRQ(1, "Cannot open IBM node file IBMSolidData.dat!");
   }
 n_v =0;
 
 if (fd) {
     
   //      fscanf(fd, "%d %d",&n_v,&n_elmt);
   fgets(string, 128, fd);
   fgets(string, 128, fd);
   
   n_v=3270;
   n_elmt=4356;
   
   ibm->n_v = n_v;
   ibm->n_elmt = n_elmt;
   
   PetscPrintf(PETSC_COMM_WORLD," nv = %d ne = %d ...\n",n_v, n_elmt);
   
   if (flg) {
   PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
   PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
   PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
   
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
   
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
   
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
   PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
   
   PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->u));
   PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->uold));
   PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->urm1));
   } else {
     for (i=0; i<n_v; i++) {
       ibm->urm1[i].x=ibm->uold[i].x;
     }
   }
   PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
   PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
   PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
   PetscReal cl = 1.;
   PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);
   
   // Loop around all the elements - node points
   for (i=0; i<n_v; i++) {
     fscanf(fd, "%le\n", &(ibm->x_bp[i]));
   }
   for (i=0; i<n_v; i++) {
     fscanf(fd, "%le\n", &ibm->y_bp[i]);
   }
   for (i=0; i<n_v; i++) {
     fscanf(fd, "%le\n", &ibm->z_bp[i]);
   }
   for (i=0; i<n_elmt; i++) {
     fscanf(fd, "%le\n", &ibm->u[i].x);
   }
   for (i=0; i<n_elmt; i++) {
     fscanf(fd, "%le\n", &ibm->u[i].y);
   }
   for (i=0; i<n_elmt; i++) {
     fscanf(fd, "%le\n", &ibm->u[i].z);
   }
   for (i=0; i<n_elmt; i++) {
     fscanf(fd, "%le\n", &ibm->uold[i].x);
   }
   //Starting to scan all the elements
   for (i=0; i<n_elmt; i++) 
     {
       
       /* if (flg) {//3 */
/* 	 fscanf(fd, "%i %i %i\n", nv2+i, nv1+i, nv3+i); */
/* 	 if (i==0) PetscPrintf(PETSC_COMM_WORLD, "IBM tec READ elemt REVVVV!!!\n"); */
/*        } */
/*        else  */
	 fscanf(fd, "%i %i %i\n", nv1+i, nv2+i, nv3+i);
       
       nv1[i] = nv1[i] - 1; 
       nv2[i] = nv2[i] - 1; 
       nv3[i] = nv3[i] - 1;
       
     }
   
   ibm->nv1 = nv1; 
   ibm->nv2 = nv2; 
   ibm->nv3 = nv3;
   
   
   
   fclose(fd);
 }//End of file IBMSolidData_0 exists
 
 /*    PetscInt ti=0; */
/*     sprintf(filen, "surface_nf%3.3d_%2.2d.dat",ti,nt); */
/*     fd = fopen(filen, "w"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt); */
/*     for (i=0; i<n_v; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->x_bp[i]); */
/*     } */
/*     for (i=0; i<n_v; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->y_bp[i]); */
/*     } */
/*     for (i=0; i<n_v; i++) {	 */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->z_bp[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_x[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_y[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_z[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_x[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_y[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_z[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_x[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_y[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_z[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, fd, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*     } */
/*     fclose(fd); */



     PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT successfully \n ");

  return(0);

}

PetscErrorCode ibm_dsdt_allout(IBMNodes *ibm, PetscInt ti)
{
  PetscInt rank, i;
  // 72bpm->T=60/72 dt=T/1000
  PetscReal dt=60./72./1000.;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  //  PetscBarrier(PETSC_NULL);
  for (i=0;i<ibm->n_v;i++) {
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "dSdt_history%2.2du.dat", i);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le \n", ti, (ibm->uold[i].x-ibm->urm1[i].x)/dt);
    fclose(f);
  }
  }
  return(0);
}

PetscErrorCode ibm_stress_allout(IBMNodes *ibm, PetscInt ti)
{
  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  //  PetscBarrier(PETSC_NULL);
  for (i=0;i<ibm->n_v;i++) {
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Stress_history%2.2dd.dat", i);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le \n", ti, ibm->uold[i].x);
    fclose(f);
  }
  }
  return(0);
}


PetscErrorCode ibm_stress_out(IBMNodes *ibm, PetscInt elmt, PetscInt ti)
{
  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Stress_history%2.2du.dat", elmt);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le \n", ti, ibm->uold[elmt].x);
    fclose(f);
  }
  return(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{

  IBMNodes	*ibm;
  PetscInt      elmt=400;
  
  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscTruth STR = PETSC_FALSE;
  PetscOptionsGetTruth(PETSC_NULL, "-str", &STR, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-elmt", &elmt, PETSC_NULL);

  PetscErrorCode flag,flg=PETSC_TRUE;

  PetscInt ti,tis, tie, tsteps=5;
  PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");    
  }

  PetscOptionsGetInt(PETSC_NULL, "-tie", &tie, &flag);
  if (!flag) {
    tie = tis;
  }

  PetscOptionsGetInt(PETSC_NULL, "-ts", &tsteps, &flag);
  if (!flag) {
    tsteps = 5; /* Default increasement is 5 */
  }

  PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);

  for (ti=tis; ti<=tie; ti+=tsteps) {
    if (ti>tis) flg=PETSC_FALSE;

    ibm_read_tecplot(ibm, ti, flg);

    if (ti>tis) 
      ibm_dsdt_allout(ibm, ti);

/*     if (STR)  */
/*       ibm_stress_allout(ibm, ti); */
/*     else */
/*       ibm_stress_out(ibm, elmt, ti); */
  }

  PetscFinalize();
  return(0);
}
