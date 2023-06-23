static char help[] = "Testing programming!";

#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscksp.h"
#include "petscvec.h"

typedef struct {
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
  PetscReal    M_x,M_y,M_z; // Moments
  PetscReal    x_c,y_c,z_c;
} FSInfo;


PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ibi, PetscInt ti)
{
  PetscInt  i;
  PetscReal t;

  FILE *f;
  char filen[80];  
  sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti, ibi);
  f = fopen(filen, "r");
  if (!f) {
    SETERRQ(PETSC_COMM_WORLD,1, "Cannot open FSI DATA file");
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data cannot open file !!!!!!!!!!!!\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input begin %d %s\n",ti,filen);

  fscanf(f, "%le %le %le", &t, &t, &t);	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp mu %le %le %le \n",FSinf->red_vel,FSinf->damp,FSinf->mu_s);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
 

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  if (rheology){
    fscanf(f, "%le %le %le ", &FSinf->R[0][0],&FSinf->R[0][1],&FSinf->R[0][2]);
    fscanf(f, "%le %le %le ", &FSinf->R[1][0],&FSinf->R[1][1],&FSinf->R[1][2]);
    fscanf(f, "%le %le %le ", &FSinf->R[2][0],&FSinf->R[2][1],&FSinf->R[2][2]);
    fscanf(f, "%le %le %le ", &FSinf->L_n[0],&FSinf->L_n[1],&FSinf->L_n[2]);
  }
  fclose(f);
  if (rheology)  PetscPrintf(PETSC_COMM_WORLD,"%le %le %le %le %le %le \n",FSinf->R[0][0],FSinf->R[0][1],FSinf->R[0][2],FSinf->R[1][1],FSinf->R[1][2],FSinf->R[2][2]); 
  else{
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data input y, dy/dt  %le %le %le %le\n",FSinf->S_new[2],FSinf->S_new[3],FSinf->red_vel,FSinf->damp);
  }
 
  return(0);
}

PetscErrorCode FSI_DATA_output(FSInfo *FSinf, PetscInt ti)
{
  FILE *f;
  char filen[80];
  sprintf(filen, "FSI_data");
  f = fopen(filen, "a");
  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti,FSinf->S_new[4],FSinf->S_new[5],FSinf->F_z, FSinf->S_ang_n[0], FSinf->S_ang_n[1],FSinf->M_x);
  fclose(f);
  return(0);
}


int main(int argc, char **argv)
{
  PetscInt      ti,tistart,tiend;
  FSInfo        FSI;
  
  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscErrorCode flg;

  PetscOptionsGetInt(PETSC_NULL, "-tis", &tistart, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-tie", &tiend, &flg);

  for (ti=tistart;ti<=tiend;ti++) {
    FSI_DATA_Input(&FSI, ti);
    FSI_DATA_output(&FSI, ti);
  }

  PetscFinalize();

  return(0);
}
