#include "variables.h"
#define CROSS(dest, v1, v2) \
	dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
	dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
	dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define SUB(dest, v1, v2) \
	dest[0] = v1[0] - v2[0]; \
	dest[1] = v1[1] - v2[1]; \
	dest[2] = v1[2] - v2[2];

extern PetscInt ti, LV;
/* // */
extern char orient[];
/* // */
extern PetscInt tiout,STRONG_COUPLING,rheology,NumberOfBodies;
extern PetscInt dgf_z,dgf_y,dgf_x;
extern PetscInt dgf_az,dgf_ay,dgf_ax;
extern PetscReal St_exp,wavelength;
extern PetscReal Flux_in,CMx_c,CMy_c,CMz_c;
extern PetscReal FluxInSum;//, FluxOutSum;

/* ==================================================================================             */
PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ibi)
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
  //  fscanf(f, "%le %le %le", &(FSinf->red_vel), &(FSinf->damp), &(FSinf->mu_s));	  
  fscanf(f, "%le %le %le", &t, &t, &t);	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp mu %le %le %le \n",FSinf->red_vel,FSinf->damp,FSinf->mu_s);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
  //PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input y, dy/dt  %le %le %le %le\n",FSinf->S_new[2],FSinf->S_new[3],FSinf->red_vel,FSinf->damp);
/*   PetscPrintf(PETSC_COMM_WORLD, "FSI_data input angle_x, dang_x/dt  %le %le %le %le\n",FSinf->S_ang_n[0],FSinf->S_ang_n[1],FSinf->red_vel,FSinf->damp); */
  return(0);
}
/* ==================================================================================             */



PetscErrorCode FSI_DATA_Input1(FSInfo *FSinf, PetscInt ibi)
{
  PetscInt  i,ibi1=ibi-3,ti1=ti-1;
  PetscReal t;

  FILE *f;
  char filen[80];  
  sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti1, ibi1);
  f = fopen(filen, "r");
  if (!f) {
    SETERRQ(PETSC_COMM_WORLD,1, "Cannot open FSI DATA file");
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data cannot open file !!!!!!!!!!!!\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input begin %d %s\n",ti1,filen);
  //  fscanf(f, "%le %le %le", &(FSinf->red_vel), &(FSinf->damp), &(FSinf->mu_s));	  
  fscanf(f, "%le %le %le", &t, &t, &t);	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp mu %le %le %le \n",FSinf->red_vel,FSinf->damp,FSinf->mu_s);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
  //PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input y, dy/dt  %le %le %le %le\n",FSinf->S_new[2],FSinf->S_new[3],FSinf->red_vel,FSinf->damp);
/*   PetscPrintf(PETSC_COMM_WORLD, "FSI_data input angle_x, dang_x/dt  %le %le %le %le\n",FSinf->S_ang_n[0],FSinf->S_ang_n[1],FSinf->red_vel,FSinf->damp); */
  return(0);
}


/* ==================================================================================             */

/* ==================================================================================             */
PetscErrorCode FSI_DATA_Output(FSInfo *FSinfo, PetscInt ibi)
{
  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "FSI_position%2.2d",ibi);
    f = fopen(filen, "a");

    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_new[4],FSinfo->S_new[5],FSinfo->F_z,FSinfo->S_new[0],FSinfo->S_new[1],FSinfo->F_x, FSinfo->Power);
    fclose(f);



    sprintf(filen, "FSI_Agnle%2.2d",ibi);
    f = fopen(filen, "a");
    if (rheology)  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le \n",ti,FSinfo->S_ang_n[1],FSinfo->M_x,FSinfo->S_ang_n[3],FSinfo->S_ang_n[5],FSinfo->alpha[0],FSinfo->alpha[1],FSinfo->alpha[2]);
    else
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_ang_n[0],FSinfo->S_ang_n[1],FSinfo->M_x,FSinfo->S_ang_n[2],FSinfo->S_ang_n[3],FSinfo->S_ang_n[4],FSinfo->S_ang_n[5]);
    fclose(f);

    sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti, ibi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	  
    for (i=0; i<6; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]);
    }
    if (rheology){
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->R[0][0],FSinfo->R[0][1],FSinfo->R[0][2]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->R[1][0],FSinfo->R[1][1],FSinfo->R[1][2]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->R[2][0],FSinfo->R[2][1],FSinfo->R[2][2]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->L_n[0], FSinfo->L_n[1], FSinfo->L_n[2]);
    }	
    fclose(f);
  }
  return(0);
}

/* ==================================================================================             */
/*    Old subroutine use FSI_calc_pos instead! */
/* ==================================================================================             */
PetscErrorCode FSI_LinMom(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime, PetscReal Re,PetscReal ti) 
{ 
  PetscInt     i;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces andOA Area
  
  // Calc Forces
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    FSinfo->S_old[i]=FSinfo->S_new[i];
    S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = FSinfo->F_x;
  F_y = FSinfo->F_y;
  F_z = FSinfo->F_z;

  // solve lin mom equ
  S_new[0]=S_new[0]; // x
  S_new[1]=S_new[1]; // dx/dt
  S_new[2]=S_new[2]-dt/2./dtime*(3*S_new[3]-4*S_real[2]+S_realm1[2])+S_new[3]*dt; // y
  // dy/dt
  S_new[3]=S_new[3]-dt/2./dtime*(3*S_new[3]-4*S_real[3]+S_realm1[3])
    +dt*(-2.*damp*(red_vel)*S_new[3] 
	 -(red_vel*red_vel)*S_old[2] 
	 + mu_s*F_y);

  S_new[4]=S_new[4]; //z
  S_new[5]=S_new[5]; //dz/dt

  // store results
  for (i=0;i<6;i++){
    FSinfo->S_new[i]=S_new[i];
  }

  return(0);
}
/* ==================================================================================             */

/* ==================================================================================             */
/*    Old subroutine! use FSI_calc_ang instead! */
/* ==================================================================================             */
PetscErrorCode FSI_AngMom(FSInfo *FSinfo,IBMNodes *ibm, PetscReal dt, PetscReal dtime, PetscReal Re) 
{ 
  
  PetscInt     i;
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  
  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_n[i]=FSinfo->S_ang_n[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
    FSinfo->S_ang_o[i]=FSinfo->S_ang_n[i];
    S_ang_o[i]=S_ang_n[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  M_x = FSinfo->M_x;
  M_y = FSinfo->M_y;
  M_z = FSinfo->M_z;

  // solve Ang mom equ
  S_ang_n[0]=S_ang_n[0]-dt/2./dtime*(3*S_ang_n[0]-4*S_ang_r[0]+S_ang_rm1[0])+S_ang_n[1]*dt; // x
  S_ang_n[1]=S_ang_n[1]-dt/2./dtime*(3*S_ang_n[1]-4*S_ang_r[1]+S_ang_rm1[1])
    +dt*(-2.*damp*(red_vel)*S_ang_n[1] 
	 -(red_vel*red_vel)*S_ang_o[0] 
	 + mu_s*M_x); // dx/dt
  S_ang_n[2]=S_ang_n[2];//-dt/2./dtime*(3*S_ang_n[3]-4*S_ang_r[2]+S_ang_rm1[2])+S_ang_n[3]*dt; // y
  // dy/dt
  S_ang_n[3]=S_ang_n[3];/*-dt/2./dtime*(3*S_ang_n[3]-4*S_ang_r[3]+S_ang_rm1[3])
    +dt*(-2.*damp*(red_vel)*S_ang_n[3] 
	 -(red_vel*red_vel)*S_ang_o[2] 
	 + mu_s*M_y);*/

  S_ang_n[4]=S_ang_n[4]; //z
  S_ang_n[5]=S_ang_n[5]; //dz/dt

  // store results
  for (i=0;i<6;i++){
    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  return(0);
}

/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos(FSInfo *FSinfo,IBMNodes *ibm, 
			    PetscReal dt, PetscReal dtime, 
			    PetscReal Re) 
{ 
  PetscInt     i,j;
  PetscInt     itr=13;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces and Area
  
  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    FSinfo->S_old[i]=FSinfo->S_new[i];
    //S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_old);
  F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_old);
  F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_old);

  // solve lin mom equ
  for (i=0; i< itr;i++) { 

    for (j=0;j<6;j++) {
      S_old[j]=S_new[j];
    }

    S_new[0]=S_new[0]; // x
    S_new[1]=S_new[1]; // dx/dt
    S_new[2]=S_new[2]-dt/2./dtime*
      (3.*S_new[2]-4.*S_real[2]+S_realm1[2])+S_new[3]*dt; // y
    // dy/dt
    S_new[3]=S_new[3]-dt/2./dtime*
         (3.*S_new[3]-4.*S_real[3]+S_realm1[3])+
      dt*(-2.*damp*(red_vel)*S_new[3]
	  -(red_vel*red_vel)*S_old[2]
	  + mu_s*F_y);
    
    S_new[4]=S_new[4];//-dt/2./dtime*(3*S_new[4]-4*S_real[4]+S_realm1[4])+S_new[5]*dt; //z
    S_new[5]=S_new[5];//-dt/2./dtime*(3*S_new[5]-4*S_real[5]+S_realm1[5])
    /*     +dt*(-2.*damp*(red_vel)*S_new[5] */
    /* 	 -(red_vel*red_vel)*(S_old[4]) */
    /* 	 + mu_s*F_z); //dz/dt */
  }
  
  // FSI convergence
  //PetscPrintf(PETSC_COMM_SELF, "FSI convergence z: %le  u_z:%le\n", S_new[4]-S_old[4],S_new[5]-S_old[5]);
  PetscPrintf(PETSC_COMM_WORLD, "FSI convergence z: %le  u_z:%le\n", S_new[2]-S_old[2],S_new[3]-S_old[3]);
  
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_realm1[i]=FSinfo->S_real[i];
    FSinfo->S_real[i]=S_new[i];
    FSinfo->S_new[i]=S_new[i];
  }

  // output values
  PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt  %le %le\n",S_new[4],S_new[5]);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "FSI_position");
    f = fopen(filen, "a");
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_new[4],S_new[5],F_z, S_real[4], S_real[5], S_realm1[4], S_realm1[5]); */
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_new[2],S_new[3],F_z, S_real[2], S_real[3], S_realm1[2], S_realm1[3]);
    fclose(f);

    sprintf(filen, "DATA_FSI%5.5d.dat",ti);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	  
    for (i=0; i<6; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]);
    }
    fclose(f);
  }
  return(0);
}

/* ==================================================================================             */
/* For strong-coupling */
/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos_SC(FSInfo *FSinfo,IBMNodes *ibm, 
			       PetscReal dt, PetscReal dtime, 
			       PetscReal Re) 
{ 
  PetscInt     i,j;
  PetscInt     itr=23;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces and Area
  
  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_real[i];//FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    //FSinfo->S_old[i]=FSinfo->S_new[i];
    //S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = FSinfo->F_x;
  F_y = FSinfo->F_y;
  F_z = FSinfo->F_z;
/*   F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_real); */
/*   F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_real); */
/*   F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_real); */

  PetscPrintf(PETSC_COMM_WORLD, "FSI  %le %le %le %le %le %le %le %le %le\n",red_vel,damp,mu_s,dt,dtime,F_y,F_z, S_real[2],S_realm1[2] );

  // solve lin mom equ
  for (i=0; i< itr;i++) { 

    for (j=0;j<6;j++) {
      S_old[j]=S_new[j];
    }

    if (dgf_x) {
    S_new[0]=S_new[0]-dt/2./dtime*
      (3.*S_new[0]-4.*S_real[0]+S_realm1[0])+S_new[1]*dt; // x
    S_new[1]=S_new[1]-dt/2./dtime*
         (3.*S_new[1]-4.*S_real[1]+S_realm1[1])+
      dt*(-2.*damp*(red_vel)*S_new[1]
	  -(red_vel*red_vel)*S_old[0]
	  + mu_s*F_x); // dx/dt
    }
    if (dgf_y) {
    S_new[2]=S_new[2]-dt/2./dtime*
      (3.*S_new[2]-4.*S_real[2]+S_realm1[2])+S_new[3]*dt; // y    
    // dy/dt
    S_new[3]=S_new[3]-dt/2./dtime*
         (3.*S_new[3]-4.*S_real[3]+S_realm1[3])+
      dt*(-2.*damp*(red_vel)*S_new[3]
	  -(red_vel*red_vel)*S_old[2]
	  + mu_s*F_y);
    }
    if (dgf_z) {
    S_new[4]=S_new[4]-dt/2./dtime*
      (3*S_new[4]-4*S_real[4]+S_realm1[4])+S_new[5]*dt; //z
    S_new[5]=S_new[5]-dt/2./dtime*(3*S_new[5]-4*S_real[5]+S_realm1[5])
              +dt*(-2.*damp*(red_vel)*S_new[5]
		   -(red_vel*red_vel)*(S_old[4])
		   + mu_s*F_z); //dz/dt
    }

    // FSI convergence
    //PetscPrintf(PETSC_COMM_SELF, "FSI convergence z: %le  u_z:%le\n", S_new[4]-S_old[4],S_new[5]-S_old[5]);
    //PetscPrintf(PETSC_COMM_WORLD, "FSI convergence y: %le  u_y:%le\n", S_new[2]-S_old[2],S_new[3]-S_old[3]);

  }
    
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_new[i]=S_new[i];
  }

  // output values
  PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt %le %le %le\n",S_new[4],S_new[5], F_z);
  PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt %le %le %le\n",S_new[2],S_new[3], F_y);
  PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt %le %le %le\n",S_new[0],S_new[1], F_x);

  return(0);
}

/* ==================================================================================             */
/*  integral equation For strong/weak-coupling */
/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos_intg(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt)				  
{ 
  PetscInt     i;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, F_mg=0.;//0.123; //Forces and Area
  PetscReal    w=0.75;

  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=0.;//FSinfo->S_real[i];//FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    S_old[i]=FSinfo->S_old[i];    
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_mg = 0.;//8./6. * mu_s;

  F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_real);
  F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_real);
  F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_real);
  
  PetscPrintf(PETSC_COMM_WORLD, "mass %le F_x %le F_y %le F_z %le \n",mu_s,F_x,F_y,F_z);
  
  PetscReal     L_x,L_y,L_z; 
  PetscReal	*x_bp =ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscInt      n_v,mode=0;
  if (rheology){
    
    PetscOptionsGetReal(PETSC_NULL, "-L_x", &L_x, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-L_y", &L_y, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-L_z", &L_z, PETSC_NULL);
    
    PetscPrintf(PETSC_COMM_WORLD,"ksi eta zeta length  are %le %le %le \n",L_x,L_y,L_z);
        
    xbp_min = 1.e23; xbp_max = -1.e23;
    ybp_min = 1.e23; ybp_max = -1.e23;
    zbp_min = 1.e23; zbp_max = -1.e23;
    
    n_v = ibm->n_v;
    for(i=0; i<n_v; i++) {
      
      xbp_min = PetscMin(xbp_min, x_bp[i]);
      xbp_max = PetscMax(xbp_max, x_bp[i]);
      
      ybp_min = PetscMin(ybp_min, y_bp[i]);
      ybp_max = PetscMax(ybp_max, y_bp[i]);
      
      zbp_min = PetscMin(zbp_min, z_bp[i]);
      zbp_max = PetscMax(zbp_max, z_bp[i]);
    }
    
    xbp_min -= 0.05; xbp_max += 0.05;
    ybp_min -= 0.05; ybp_max += 0.05;
    zbp_min -= 0.05; zbp_max += 0.05;

    if ( ybp_min < 0.0 && F_y <0.0) mode++;
    else if ( ybp_max > L_y && F_y > 0.0) mode++;

  }
  if (mode)  PetscPrintf(PETSC_COMM_WORLD,"Wall repulsive acts  \n");

  // solve lin mom equ.
  if (dgf_x) {
    S_new[1]=S_real[1]+ dt*(mu_s*F_x); // u=u_r + int(F/mdt)

    S_new[0]=S_real[0]+0.5*(S_new[1]+S_real[1])*dt; //x=x_r+u_avedt
  }
  if (dgf_y && mode<1) {
    S_new[3]=S_real[3]+ dt*(mu_s*F_y); // u=u_r + int(F/mdt)

    S_new[2]=S_real[2]+0.5*(S_new[3]+S_real[3])*dt; //y=y_r+u_avedt
  }
  if (dgf_z) {
    S_new[5]=S_real[5]+ dt*(mu_s*F_z-F_mg); // u=u_r + int(F/mdt)

    S_new[4]=S_real[4]+0.5*(S_new[5]+S_real[5])*dt; //x=x_r+u_avedt
  }

  // Relaxation
  if (STRONG_COUPLING) {
    FSinfo->atk=0.;
    for (i=1;i<6;i+=2){
      FSinfo->dS[i]=S_old[i]-S_new[i];
    
      if (fabs(FSinfo->dS[i]-FSinfo->dS_o[i])>1e-8  &&
     	  FSinfo->atk_o!=0.3) {
	FSinfo->atk+=(FSinfo->dS[i])/
	  (FSinfo->dS_o[i]-FSinfo->dS[i]);
      }
    }
    FSinfo->atk=FSinfo->atk_o+(FSinfo->atk_o-1)*FSinfo->atk;
    if (FSinfo->atk>.9) FSinfo->atk=.9;
    if (FSinfo->atk<-.2) FSinfo->atk=-0.2;
    
    w=(1.-FSinfo->atk);
    if (rheology)  PetscOptionsGetReal(PETSC_NULL, "-w_str", &w, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "under relaxation coefficient %le \n",w);
    
    for (i=1;i<6;i+=2){
      S_new[i]=w*S_new[i]+(1.-w)*S_old[i];
      S_new[i-1]=S_real[i-1]+0.5*(S_new[i]+S_real[i])*dt;
    }
  }

  if (rheology){
    FSinfo->x_c=FSinfo->a_c[0]+S_new[0];
    FSinfo->y_c=FSinfo->a_c[1]+S_new[2];
    FSinfo->z_c=FSinfo->a_c[2]+S_new[4];
  }
  
  // store results
  for (i=0;i<6;i++){
/*     // Relaxation  */
/*     S_new[i]=w*S_new[i]+(1-w)*S_old[i]; */
    FSinfo->S_new[i]=S_new[i];
  }
  
 for (i=0;i<3;i++){
   FSinfo->acc[i]=(S_new[2*i+1]-S_real[2*i+1])/dt;
  }

  // output values
 PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt a_z %le %le %le \n",S_new[4],S_new[5],FSinfo->acc[2]);
 PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt a_y %le %le %le \n",S_new[2],S_new[3],FSinfo->acc[1]);
 PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt a_x %le %le %le \n",S_new[0],S_new[1],FSinfo->acc[0]);

  return(0);
}
/* ==================================================================================             */
PetscErrorCode Calc_FSI_vel_particle(FSInfo *FSinfo,IBMNodes *ibm,PetscReal dt)				  
{ 
  PetscInt     i;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    mu_s; //  mass coeff
  PetscReal    F_x,F_y,F_z;//Forces 
  PetscReal    w=0.5;

   
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=0.;//FSinfo->S_real[i];//FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    S_old[i]=FSinfo->S_old[i];    
  }
  
  mu_s   =FSinfo->mu_s   ;

  F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_real);
  F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_real);
  F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_real);
  
  PetscPrintf(PETSC_COMM_WORLD, "mass %le F_x %le F_y %le F_z %le \n",mu_s,F_x,F_y,F_z);
  
  PetscReal     L_x,L_y,L_z; 
  PetscReal	*x_bp =ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscInt      n_v,mode=0;

    
  PetscOptionsGetReal(PETSC_NULL, "-L_x", &L_x, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_y", &L_y, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_z", &L_z, PETSC_NULL);
    
  PetscPrintf(PETSC_COMM_WORLD,"ksi eta zeta length  are %le %le %le \n",L_x,L_y,L_z);
        
  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;
  
  n_v = ibm->n_v;
  for(i=0; i<n_v; i++) {
    
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);
    
    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);
    
    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }
    
  xbp_min -= 0.1; xbp_max += 0.1;
  ybp_min -= 0.1; ybp_max += 0.1;
  zbp_min -= 0.1; zbp_max += 0.1;
  
  if ( ybp_min < 0.0 && F_y <0.0) mode++;
  else if ( ybp_max > L_y && F_y > 0.0) mode++;
  

  if (mode)  PetscPrintf(PETSC_COMM_WORLD,"Wall repulsive acts  \n");

  // solve lin mom equ.
  if (dgf_x)  S_new[1]=S_real[1]+ dt*(mu_s*F_x); // u=u_r + int(F/mdt)
  if (dgf_y && mode<1) S_new[3]=S_real[3]+ dt*(mu_s*F_y); // u=u_r + int(F/mdt)
  if (dgf_z)  S_new[5]=S_real[5]+ dt*(mu_s*F_z); // u=u_r + int(F/mdt)

  // Relaxation
  if (STRONG_COUPLING) {
  /*   FSinfo->atk=0.; */
/*     for (i=1;i<6;i+=2){ */
/*       FSinfo->dS[i]=S_old[i]-S_new[i]; */
    
/*       if (fabs(FSinfo->dS[i]-FSinfo->dS_o[i])>1e-8  && */
/*      	  FSinfo->atk_o!=0.3) { */
/* 	FSinfo->atk+=(FSinfo->dS[i])/ */
/* 	  (FSinfo->dS_o[i]-FSinfo->dS[i]); */
/*       } */
/*     } */
/*     FSinfo->atk=FSinfo->atk_o+(FSinfo->atk_o-1)*FSinfo->atk; */
/*     if (FSinfo->atk>.9) FSinfo->atk=.9; */
/*     if (FSinfo->atk<-.2) FSinfo->atk=-0.2; */
    
/*     w=(1.-FSinfo->atk); */
    PetscOptionsGetReal(PETSC_NULL, "-w_str", &w, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "under relaxation coefficient %le \n",w);
    
    for (i=1;i<6;i+=2){
      S_new[i]=w*S_new[i]+(1.-w)*S_old[i];
    }
  }

  
 /*  FSinfo->x_c=FSinfo->a_c[0]+S_new[0]; */
/*   FSinfo->y_c=FSinfo->a_c[1]+S_new[2]; */
/*   FSinfo->z_c=FSinfo->a_c[2]+S_new[4]; */
  
  
  // store results
  for (i=0;i<6;i++){
/*     // Relaxation  */
/*     S_new[i]=w*S_new[i]+(1-w)*S_old[i]; */
    FSinfo->S_new[i]=S_new[i];
  }
  
 for (i=0;i<3;i++){
   FSinfo->S_new[2*i+1]=S_new[2*i+1];
   FSinfo->acc[i]=(S_new[2*i+1]-S_real[2*i+1])/dt;
  }

  // output values
 PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt a_z %le %le %le \n",S_new[4],S_new[5],FSinfo->acc[2]);
 PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt a_y %le %le %le \n",S_new[2],S_new[3],FSinfo->acc[1]);
 PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt a_x %le %le %le \n",S_new[0],S_new[1],FSinfo->acc[0]);

  return(0);
}
/* ==================================================================================             */

PetscErrorCode Forced_Motion(FSInfo *fsi,
			     PetscReal dt)
{
  PetscReal  pi = 3.141592653589793;
  PetscReal t,KC;
  t = (ti)*dt;
  KC=fsi->red_vel;

/*   fsi->S_new[2]= A*sin(w*t); */
/*   fsi->S_new[3]= A*w*cos(w*t); */
/*   PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt %le %le %le\n",t,fsi->S_new[2],fsi->S_new[3]); */
  fsi->S_new[0]= KC/2./pi*sin(2*pi/KC*t);
  fsi->S_new[1]= cos(2*pi/KC*t);
  fsi->S_new[2]= KC/2./pi*cos(2*pi/KC*t);
  fsi->S_new[3]= sin(2*pi/KC*t);
/*   fsi->S_new[4]= KC/2./pi*sin(2*pi/KC*t); */
/*   fsi->S_new[5]= cos(2*pi/KC*t); */
  PetscPrintf(PETSC_COMM_WORLD, "t z, dz/dt %le %le %le\n",t,fsi->S_new[4],fsi->S_new[5]);

  return(0);
}


/* ==================================================================================             */
PetscErrorCode Calc_FSI_Ang(FSInfo *FSinfo,IBMNodes *ibm, 
			    PetscReal dt, PetscReal dtime,
			    PetscInt ibi, UserCtx *user)
{  
  PetscInt     i,itr=12,j;
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  

  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_n[i]=FSinfo->S_ang_r[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
    //FSinfo->S_ang_o[i]=FSinfo->S_ang_n[i];
    //S_ang_o[i]=S_ang_n[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  M_x = FSinfo->M_x;
  M_y = FSinfo->M_y;
  M_z = FSinfo->M_z;

  // solve Ang mom equw
  for (i=0; i< itr;i++) { 

  for (j=0;j<6;j++) {
    S_ang_o[j]=S_ang_n[j];
  }

  if (dgf_ax) {
  S_ang_n[0]=S_ang_n[0]-dt/2./dtime*
    (3*S_ang_n[0]-4*S_ang_r[0]+S_ang_rm1[0])+S_ang_n[1]*dt; // x
  S_ang_n[1]=S_ang_n[1]-dt/2./dtime*(3*S_ang_n[1]-4*S_ang_r[1]+S_ang_rm1[1])
    +dt*(-2.*damp*(red_vel)*S_ang_n[1] 
	 //-(red_vel*red_vel)*S_ang_o[0] 
	 + mu_s*M_x); // dx/dt
  }
  if (dgf_ay) {
  S_ang_n[2]=S_ang_n[2]-dt/2./dtime*(3*S_ang_n[2]-4*S_ang_r[2]+S_ang_rm1[2])+S_ang_n[3]*dt; // y
  // dy/dt
  S_ang_n[3]=S_ang_n[3]-dt/2./dtime*(3*S_ang_n[3]-4*S_ang_r[3]+S_ang_rm1[3])
    +dt*(-2.*damp*(red_vel)*S_ang_n[3] 
	 -(red_vel*red_vel)*S_ang_o[2] 
	 + mu_s*M_y);
  }
  if (dgf_az) {
  S_ang_n[4]=S_ang_n[4]-dt/2./dtime*(3*S_ang_n[4]-4*S_ang_r[4]+S_ang_rm1[4])+S_ang_n[5]*dt; //z
  S_ang_n[5]=S_ang_n[5]-dt/2./dtime*(3*S_ang_n[5]-4*S_ang_r[5]+S_ang_rm1[5])
    +dt*(-2.*damp*(red_vel)*S_ang_n[5] 
	 -(red_vel*red_vel)*S_ang_o[4] 
	 + mu_s*M_z); //dz/dt
  }

/* ==================================================================================             */
/*    making Moment from dUn/dt implicit !!! */
/* ==================================================================================             */

/*   // 1) Calc new vel of all surf nodes */
/*   wx= S_ang_n[1]; */
/*   for (nv=0; nv<ibm->n_v; i++) { */
/*     rx = ibm->x_bp[nv]-x_c; */
/*     ry = ibm->y_bp[nv]-y_c; */
/*     rz = ibm->z_bp[nv]-z_c;       */
/*     ibm->u[nv].x =   ry*wz-wy*rz  ; */
/*     ibm->u[nv].y =-( rx*wz-wx*rz ); */
/*     ibm->u[nv].z =   rx*wy-wx*ry  ;       */
/*   } */

/*   // 2) ibm interpolation to calc new P+dU/dt on all ibm nodes */
/*   PetscPrintf(PETSC_COMM_SELF, "IBM_INP\n"); */
/*   ibm_interpolation_advanced(user, ibm, FSinfo, ibi); */


/*   // 3) calc new moments */
/*   PetscPrintf(PETSC_COMM_SELF, "Calc_f\n"); */
/*   Calc_forces_SI(FSinfo,user,ibm, ti, ibi, 0); */

/*   M_z = FSinfo->M_z; */

/* ==================================================================================             */
/* ==================================================================================             */

  }

  // FSI convergence
  PetscPrintf(PETSC_COMM_WORLD, "FSI convergence ang_x: %le  w_x:%le\n", S_ang_n[0]-S_ang_o[0],S_ang_n[1]-S_ang_o[1]);

  // store results
  for (i=0;i<6;i++){
/*     FSinfo->S_ang_rm1[i]=FSinfo->S_ang_r[i]; */
/*     FSinfo->S_ang_r[i]=S_ang_n[i]; */

    // Relaxation
    //    S_ang_n[i]=w*S_ang_n[i]+(1-w)*FSinfo->S_ang_o[i];

    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  // output values
  PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x  %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_rm1[0]);
/*   PetscPrintf(PETSC_COMM_WORLD, "Ang_y, dAng_y/dt M_y  %le %le %le %le %le\n",S_ang_n[2],S_ang_n[3], M_y, S_ang_r[2],S_ang_rm1[2]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Ang_z, dAng_z/dt M_z  %le %le %le %le %le\n",S_ang_n[4],S_ang_n[5], M_z, S_ang_r[4],S_ang_rm1[4]); */
/*   PetscInt rank; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   PetscBarrier(PETSC_NULL); */
/*   if (!rank) { */
/*     FILE *f; */
/*     char filen[80]; */
/*     sprintf(filen, "FSI_Agnle"); */
/*     f = fopen(filen, "a"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_ang_n[0],S_ang_n[1],M_x,S_ang_r[0],S_ang_r[1],S_ang_rm1[0],S_ang_rm1[1]); */
/*     fclose(f); */

/*     sprintf(filen, "DATA_FSI%5.5d.dat",ti); */
/*     f = fopen(filen, "w"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	   */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	   */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	   */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	   */
/*     for (i=0; i<6; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]); */
/*     } */
/*     fclose(f); */
/*   } */

  return(0);
}

/* ==================================================================================             */
PetscErrorCode Calc_FSI_Ang_intg(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt,PetscInt itrSC,
				 PetscInt ibi, UserCtx *user, 
				 PetscReal S_init)
{  
PetscInt     i,itr=12,j,nv;
  PetscReal    pi=3.141592654;
  PetscReal    S_ang_n[6],S_ang_r[6],S_ang_rm1[6],S_ang_o[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
  PetscReal    rx,ry,rz;
  PetscReal    x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  PetscReal    wx=0., wy=0., wz=0.;
  PetscReal    w=.5, wf=1.;//0.82;
  PetscReal    Mdpdn_x;

  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_o[i]=FSinfo->S_ang_o[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;
  Mdpdn_x=FSinfo->Mdpdn_x;

  // Trapezoidal rule
/*   if (STRONG_COUPLING) */
/*   M_x = 0.5*(wf*FSinfo->M_x+(1-wf)*FSinfo->M_x_old */
/* 	     + FSinfo->M_x_real); */
/*     M_x = FSinfo->M_x; // 1st order */
/*   else */
/*     M_x = (55.*FSinfo->M_x - 59.*FSinfo->M_x_real + 37.*FSinfo->M_x_rm2 -9.*FSinfo->M_x_rm3)/24.; */

  M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_real);
  M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real);
  M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real);

/*   M_x = (2.*FSinfo->M_x - FSinfo->M_x_real); */
/*   M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real); */
/*   M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real); */

/*   if (itrSC==1) { */
/*     M_x = FSinfo->M_x; */
 /*     M_y = FSinfo->M_y; */
/*     M_z = FSinfo->M_z; */
/*   } else { */
/* /\*     M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_old); *\/ */
/* /\*     M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_old); *\/ */
/* /\*     M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_old); *\/ */
/*   M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_real); */
/*   M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real); */
/*   M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real); */
/*   } */

  // solve Ang mom equ using trapozoidal rule u(n+1)=u(n)+dt*F(n+1/2)
  /* if (dgf_ax) { */
  /* // */
  if (strcmp(orient,"xx00")==0) {
    S_ang_n[1]=(1-damp*dt)/(1.+damp*dt)*S_ang_r[1]
      + dt/(1.+damp*dt)*(mu_s*M_x);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt)
    S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt
  } else if (strcmp(orient,"xy45")==0) {
    S_ang_n[1]=(1-damp*dt)/(1.+damp*dt)*S_ang_r[1]
      + dt/(1.+damp*dt)*(mu_s*(M_x*cos(pi/4) + M_y*cos(pi/4)));///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt)
    S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt
  } else if (strcmp(orient,"yy00")==0) {
    S_ang_n[1]=(1-damp*dt)/(1.+damp*dt)*S_ang_r[1]
      + dt/(1.+damp*dt)*(mu_s*M_y);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt)
    S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt
  }
  /* // */
/*     S_ang_n[1]=((1-damp*dt-0.25*dt*dt*red_vel*red_vel)*S_ang_r[1] */
/* 		+ dt*(mu_s*M_x) */
/* 		- dt*red_vel*red_vel*(S_ang_r[0]-S_init))/ */
/*       (1.+damp*dt+0.25*dt*dt*red_vel*red_vel);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt) */
/*     S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt */

/*     S_ang_n[1]=S_ang_rm1[1]-damp*2.*dt*S_ang_r[1] */
/*                       + 2*dt*(mu_s*M_x);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt) */
/*     S_ang_n[0]=S_ang_rm1[0]+2*S_ang_r[1]*dt; //ang=ang_r+w_avedt */
  /* } */
  /* if (dgf_ay) { */
  /*   S_ang_n[3]=0.0;//S_ang_r[3]+ dt*(mu_s*M_y); // w=w_r + int(M/Idt) */

  /*   S_ang_n[2]=0.0;//S_ang_r[2]+0.5*(S_ang_n[3]+S_ang_r[3])*dt; //ang=ang_r+w_avedt */
  /* } */
  /* if (dgf_az) { */
  /*   S_ang_n[5]=0.0;//S_ang_r[5]+ dt*(mu_s*M_z); // w=w_r + int(M/Idt) */

  /*   S_ang_n[4]=0.0;//S_ang_r[4]+0.5*(S_ang_n[5]+S_ang_r[5])*dt; //ang=ang_r+w_avedt */
  /* } */
  /* // */
    // Relaxation
  PetscInt m, n;
  m = 1; n = 2;
  /* // */
  if (STRONG_COUPLING) {
    FSinfo->atk=0.;
    for (i=m;i<n;i+=2){
      FSinfo->dS[i]=S_ang_o[i]-S_ang_n[i];
    
      if (fabs(FSinfo->dS[i]-FSinfo->dS_o[i])>1e-8  &&
     	  FSinfo->atk_o!=0.3) {
	FSinfo->atk+=(FSinfo->dS[i])/
	  (FSinfo->dS_o[i]-FSinfo->dS[i]);
      }
    }
    FSinfo->atk=FSinfo->atk_o+(FSinfo->atk_o-1)*FSinfo->atk;
    if (FSinfo->atk>.9) FSinfo->atk=.9;
    if (FSinfo->atk<-.2) FSinfo->atk=-0.2;
    
    w=1.-FSinfo->atk;
    for (i=m;i<n;i+=2){
      S_ang_n[i]=w*S_ang_n[i]+(1.-w)*S_ang_o[i];
      S_ang_n[i-1]=S_ang_r[i-1]+0.5*(S_ang_n[i]+S_ang_r[i])*dt;
    }
  }

  // store results
  for (i=0;i<6;i++){
    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  // output values
  //PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x w %le %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_rm1[0],w);
  /* // */
  if (strcmp(orient,"xx00")==0){
    PetscPrintf(PETSC_COMM_WORLD, "%d Ang_x, dAng_x/dt M_x  %le %le %le %le %le %le %le %le\n",ibi,S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_r[1], w,1.-FSinfo->atk_o,FSinfo->dS[1],FSinfo->dS_o[1]);
  } else if (strcmp(orient,"xy45")==0) {
    PetscPrintf(PETSC_COMM_WORLD, "%d Ang_x, dAng_x/dt M_x  %le %le %le %le %le %le %le %le\n",ibi,S_ang_n[0],S_ang_n[1], (M_x*cos(pi/4) + M_y*cos(pi/4)), S_ang_r[0],S_ang_r[1], w,1.-FSinfo->atk_o,FSinfo->dS[1],FSinfo->dS_o[1]);
  } else if (strcmp(orient,"yy00")==0) {
    PetscPrintf(PETSC_COMM_WORLD, "%d Ang_y, dAng_y/dt M_y  %le %le %le %le %le %le %le %le\n",ibi,S_ang_n[0],S_ang_n[1], M_y, S_ang_r[0],S_ang_r[1], w,1.-FSinfo->atk_o,FSinfo->dS[1],FSinfo->dS_o[1]);
  }
  /* // */
  return(0);
}


/* ==================================================================================             */
PetscErrorCode Calc_FSI_Ang_staggered(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt, 
				 PetscInt ibi, UserCtx *user)
{  
  PetscInt     i;
  PetscReal    S_ang_n[6],S_ang_r[6],S_ang_rm1[6],S_ang_o[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    M_x,M_y,M_z; //Forces and Area
 
  PetscReal    Mdpdn_x;

  // Calc Forces
  //Calc_Moments(&FSinfo,&ibm,Re);
  
  // init values
  for (i=0;i<6;i++) {
    S_ang_o[i]=FSinfo->S_ang_o[i];
    S_ang_r[i]=FSinfo->S_ang_r[i];
    S_ang_rm1[i]=FSinfo->S_ang_rm1[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;
  Mdpdn_x=FSinfo->Mdpdn_x;

  M_x = FSinfo->M_x;
  M_y = FSinfo->M_y;
  M_z = FSinfo->M_z;

  // solve Ang mom equ
  if (dgf_ax) {
    S_ang_n[1]=(1-damp*dt)/(1.+damp*dt)*S_ang_r[1]
                      + dt/(1.+damp*dt)*(mu_s*M_x);///(1.+mu_s*Mdpdn_x); // w=w_r + int(M/Idt)

    S_ang_n[0]=S_ang_r[0]+0.5*(S_ang_n[1]+S_ang_r[1])*dt; //ang=ang_r+w_avedt
  }
  if (dgf_ay) {
    S_ang_n[3]=S_ang_r[3]+ dt*(mu_s*M_y); // w=w_r + int(M/Idt)

    S_ang_n[2]=S_ang_r[2]+0.5*(S_ang_n[3]+S_ang_r[3])*dt; //ang=ang_r+w_avedt
  }
  if (dgf_az) {
    S_ang_n[5]=S_ang_r[5]+ dt*(mu_s*M_z); // w=w_r + int(M/Idt)

    S_ang_n[4]=S_ang_r[4]+0.5*(S_ang_n[5]+S_ang_r[5])*dt; //ang=ang_r+w_avedt
  }

  
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_ang_n[i]=S_ang_n[i];
  }
  
  // output values
  PetscPrintf(PETSC_COMM_WORLD, "Ang_x, dAng_x/dt M_x  %le %le %le %le %le\n",S_ang_n[0],S_ang_n[1], M_x, S_ang_r[0],S_ang_rm1[0]);

  return(0);
}


/* ==================================================================================             */
PetscErrorCode Elmt_Move_FSI_TRANS(FSInfo *FSinfo, IBMNodes *ibm
				   , PetscInt ibi)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscInt i;

/*   for (i=0; i<n_v; i++) { */
/*     ibm->x_bp_o[i] = ibm->x_bp[i]; */
/*     ibm->y_bp_o[i] = ibm->y_bp[i]; */
/*     ibm->z_bp_o[i] = ibm->z_bp[i]; */
/*   } */

  PetscPrintf(PETSC_COMM_WORLD, "MOVE BODY x: %le  y:%le z:%le\n", FSinfo->S_new[0],FSinfo->S_new[2],FSinfo->S_new[4]);

  for (i=0; i<n_v; i++) {
    // change for stat case 4/9/06 iman
    ibm->x_bp[i] = ibm->x_bp0[i]+(FSinfo->S_new[0]);//-FSinfo->S_real[0]);
    ibm->y_bp[i] = ibm->y_bp0[i]+(FSinfo->S_new[2]);//-FSinfo->S_real[2]);
    ibm->z_bp[i] = ibm->z_bp0[i]+(FSinfo->S_new[4]);//-FSinfo->S_real[4]);
  }
  
  for (i=0; i<n_v; i++) {
    ibm->u[i].x = FSinfo->S_new[1];
    ibm->u[i].y = FSinfo->S_new[3];
    ibm->u[i].z = FSinfo->S_new[5];
  }

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d_%2.2d.dat",ti,ibi);
      //      sprintf(filen, "surface%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }

  return(0);
}
/* ==================================================================================             */

/* ==================================================================================             */
PetscErrorCode Elmt_Move_FSI_ROT(FSInfo *FSinfo, IBMNodes *ibm, 
				 PetscReal dt, PetscInt ibi)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  PetscReal rot;

  /* // */
  rot=FSinfo->S_ang_n[0];//-FSinfo->S_ang_o[0];
  /* // */
  PetscInt i;
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal rx,ry,rz;
  /* // */
  PetscReal wx,wy,wz;  // used to set velocities of immersed boundary points.
  /* // */
  PetscPrintf(PETSC_COMM_WORLD, "Immesed Body Rotation: Axis- %-s, Angle: %le, Center %le %le %le\n",orient,rot, x_c,y_c,z_c);
    
/* //Rotation of MHV around hinge */
  PetscReal pi = 3.14159265359, angle = 0.0, R[4][4], x[4], xn[4], u, v, w, a, b, c;
  angle = (-rot*pi)/180.0; 
  PetscInt m,n;

  // For more info check out: https://en.wikipedia.org/wiki/Rotation_matrix.

  if (strcmp(orient,"xx00")==0) {  // axis of rotation parallel to x-axis.
    u = 1. ; v = 0.; w = 0.;
  } else if (strcmp(orient,"yy00")==0){  // axis of rotation parallel to y-axis.
    u = 0.; v = 1.; w = 0.; 
  } else if (strcmp(orient,"zz00")==0){  // axis of rotation parallel to z-axis.
    u = 0.; v = 0.; w = 1.;
  } else if (strcmp(orient,"xy45")==0){
    u = cos(pi/4); v = cos(pi/4); w=0;  // axis of  rotation parallel to a line 45 degs between x and y axis.
  } else if (strcmp(orient,"xz45")==0){
    u = cos(pi/4); v =  0; w = cos(pi/4);
  } else if (strcmp(orient,"yz45")==0){
    u = 0; v = cos(pi/4); w = cos(pi/4);
  }
  
//  PetscPrintf(PETSC_COMM_WORLD, "Rotation - u,v,w: %le,%le,%le \n",u,v,w);
  a = x_c; b = y_c; c = z_c; //Rotation origin
  for (i=0; i<n_v; i++) { 
    x[0] = ibm->x_bp0[i]; x[1] = ibm->y_bp0[i]; x[2] = ibm->z_bp0[i]; x[3] = 1.;
    
    R[0][0] = u*u + (v*v + w*w)*cos(angle);
    R[0][1] = u*v*(1 - cos(angle)) - w*sin(angle);
    R[0][2] = u*w*(1 - cos(angle)) + v*sin(angle);
    R[0][3] = (a*(v*v + w*w) - u*(b*v + c*w))*(1 - cos(angle)) + (b*w -c*v)*sin(angle);
    
    R[1][0] = u*v*(1 - cos(angle)) + w*sin(angle);
    R[1][1] = v*v + (u*u + w*w)*cos(angle);
    R[1][2] = v*w*(1 - cos(angle)) - u*sin(angle);
    R[1][3] = (b*(u*u + w*w) - v*(a*u + c*w))*(1 - cos(angle)) + (c*u - a*w)*sin(angle);
    
    R[2][0] = u*w*(a - cos(angle)) - v*sin(angle);
    R[2][1] = v*w*(1 - cos(angle)) + u*sin(angle);
    R[2][2] = w*w + (u*u + v*v)*cos(angle);
    R[2][3] = (c*(u*u + v*v) - w*(a*u + b*v))*(1 - cos(angle)) + (a*v - b*u)*sin(angle);
    
    R[3][0] = 0.;
    R[3][1] = 0.;
    R[3][2] = 0.;
    R[3][3] = 1.;
    
    for (m=0; m<4; m++) {
      xn[m] = 0.;
      for (n=0; n<4; n++) {
	xn[m] += R[m][n]*x[n];
      }
    }
    ibm->x_bp[i] = xn[0]; ibm->y_bp[i] = xn[1]; ibm->z_bp[i] = xn[2];
  }
  /* // */

  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
    
    /*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
    //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);
    // Addedd 5/7/06 iman
    // ns = nf x k
    if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
	(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
      ibm->ns_x[i] = 1.;
      ibm->ns_y[i] = 0.;
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = 0.;
      ibm->nt_y[i] = 1.;
      ibm->nt_z[i] = 0.;
      } else {
      ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
    }

  }
  /* // */
  if (strcmp(orient,"xx00")==0){
  wx = FSinfo->S_ang_n[1], wy = 0., wz = 0.; 
  } else if (strcmp(orient,"yy00")==0){
    wy = FSinfo->S_ang_n[1];
    wx = 0.;
    wz = 0.;
  } else if (strcmp(orient,"zz00")==0){
    wz = FSinfo->S_ang_n[1];
    wx = 0.;
    wy = 0.;
  }else if (strcmp(orient,"xy45")==0){
    wy = FSinfo->S_ang_n[1]*cos(pi/4);
    wx = wy;
    wz = 0.;
  }else if (strcmp(orient,"xz45")==0){
    wx = FSinfo->S_ang_n[1]*cos(pi/4);
    wz = wx;
    wy = 0.;
  }else if (strcmp(orient,"yz45")==0){
    wy = FSinfo->S_ang_n[1]*cos(pi/4);
    wz = wy;
    wx = 0.;
  }
  /* // */
  // 1st order approx. not good!
  for (i=0; i<n_v; i++) {
      /* ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / dt; */
      /* ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / dt; */
      /* ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / dt; */
    rx = ibm->x_bp[i]-x_c;
    ry = ibm->y_bp[i]-y_c;
    rz = ibm->z_bp[i]-z_c;
    ibm->u[i].x =   ry*wz-wy*rz  ;
    ibm->u[i].y =-( rx*wz-wx*rz );
    ibm->u[i].z =   rx*wy-wx*ry  ;
  }
  
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d_%2.2d.dat",ti,ibi);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }
  return(0);
}
/* ==================================================================================             */

/* ==================================================================================             */

PetscErrorCode CollisionDetectionOfCylinders(FSInfo *fsi,IBMNodes *ibm, PetscInt NumberOfBodies)

{
  PetscReal     x_c,y_c,z_c;
  PetscReal     x_c2,y_c2,z_c2;
  PetscReal     l_c;
  PetscReal     n_x,n_y,n_z; //collision direction
  PetscReal     v_x=0.,v_y=0.,v_z=0.;
  PetscReal     v_x2=0.,v_y2=0.,v_z2=0.;
  PetscReal     v_n1, v_t1; //vel in collision direction
  PetscReal     v_n2, v_t2;
  PetscBool    collision;
  PetscInt	ibi,ibi2,i;
  PetscReal     eps=0.05,y1_wall=0.,y2_wall=10.;
  collision= PETSC_FALSE;

  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    x_c=fsi[ibi].x_c+fsi[ibi].S_new[0]; 
    y_c=fsi[ibi].y_c+fsi[ibi].S_new[2]; 
    z_c=fsi[ibi].z_c+fsi[ibi].S_new[4];

    // collision with walls
    if (y_c<y1_wall+0.5+eps || y_c>y2_wall-0.5-eps) {
      PetscPrintf(PETSC_COMM_WORLD, "Collision Detected!!!! cylinder %d with wall 1\n", ibi);
      PetscPrintf(PETSC_COMM_WORLD,"          Position!!!! cyl1 %le %le %le  wall1 %le\n", x_c,y_c,z_c,y1_wall);
      
      x_c2=x_c; 
      if (y_c>y2_wall-0.5-eps) 
	 y_c2=2*y2_wall-y_c; 
      else
	y_c2=2*y1_wall-y_c;
      z_c2=z_c;     

      l_c=sqrt((x_c-x_c2)*(x_c-x_c2) + 
	       (y_c-y_c2)*(y_c-y_c2) +
	       (z_c-z_c2)*(z_c-z_c2));

      // Collision Direction
      n_x = 0.;//(x_c-x_c2)/l_c;
      n_y = (y_c-y_c2)/l_c;
      n_z = (z_c-z_c2)/l_c;

      fsi[ibi].S_new[2]+= (1.0-l_c)*n_y;

      v_x=fsi[ibi].S_new[1];
      v_y=fsi[ibi].S_new[3];
      v_z=fsi[ibi].S_new[5];
      v_n1 = v_x*n_x + v_y*n_y + v_z*n_z;
      v_t1 = v_x*n_x + v_y*n_z - v_z*n_y;

      v_x = v_n2*n_x + v_t1 *n_x;
      v_y = v_n2*n_y + v_t1 *n_z;
      v_z = v_n2*n_z - v_t1 *n_y;
            
      fsi[ibi].S_new[1]=v_x;
      fsi[ibi].S_new[3]=v_y;
      fsi[ibi].S_new[5]=v_z;      
    }


    for (ibi2=ibi+1;ibi2<NumberOfBodies;ibi2++){
            
      x_c2=fsi[ibi2].x_c+fsi[ibi2].S_new[0]; 
      y_c2=fsi[ibi2].y_c+fsi[ibi2].S_new[2]; 
      z_c2=fsi[ibi2].z_c+fsi[ibi2].S_new[4];

      l_c=sqrt((x_c-x_c2)*(x_c-x_c2) + 
	       (y_c-y_c2)*(y_c-y_c2) +
	       (z_c-z_c2)*(z_c-z_c2));

      if (l_c < 1.0) { // collission !!!
	collision= PETSC_TRUE;

	PetscPrintf(PETSC_COMM_WORLD, "Collision Detected!!!! cylinder %d with %d\n", ibi, ibi2);
	PetscPrintf(PETSC_COMM_WORLD,"          Position!!!! cyl1 %le %le %le  cyl2 %le %le %le\n", x_c,y_c,z_c,x_c2,y_c2,z_c2);

	// Collision Direction
	n_x = 0.;//(x_c-x_c2)/l_c;
	n_y = (y_c-y_c2)/l_c;
	n_z = (z_c-z_c2)/l_c;
	
	/* Move the 2nd cyl to the 1D distance of 1st Cyl */
	fsi[ibi2].S_new[0]-= (1.0-l_c)*n_x;
	fsi[ibi2].S_new[2]-= (1.0-l_c)*n_y;
	fsi[ibi2].S_new[4]-= (1.0-l_c)*n_z;

	PetscPrintf(PETSC_COMM_WORLD,"  New     Position!!!! cyl1 %le %le %le  cyl2 %le %le %le\n", x_c,y_c,z_c,fsi[ibi2].x_c+fsi[ibi2].S_new[0],fsi[ibi2].y_c+fsi[ibi2].S_new[2],fsi[ibi2].z_c+fsi[ibi2].S_new[4]);

	/* Change the Vel to the collsion Vel! */
	/*

	v'.t=v.t
	v'2.t=v2.t
        
	v'.n=v2.n
        v'2.n=v.n

        ' is the new vel (after collision)
        n is the collision direction vector
        t is the bi-normal of collision direction vector

	*/

	v_x=fsi[ibi].S_new[1];
	v_y=fsi[ibi].S_new[3];
	v_z=fsi[ibi].S_new[5];
	v_n1 = v_x*n_x + v_y*n_y + v_z*n_z;
	v_t1 = v_x*n_x + v_y*n_z - v_z*n_y;
	  
	v_x2=fsi[ibi2].S_new[1];
	v_y2=fsi[ibi2].S_new[3];
	v_z2=fsi[ibi2].S_new[5];
	v_n2 = v_x2*n_x + v_y2*n_y + v_z2*n_z;
	v_t2 = v_x2*n_x + v_y2*n_z - v_z2*n_y;

	PetscPrintf(PETSC_COMM_WORLD,"          Velocity!!!! cyl1 %le %le %le  cyl2 %le %le %le\n", v_x,v_y,v_z,v_x2,v_y2,v_z2);
/* 	PetscBarrier(PETSC_NULL); */

	v_x = v_n2*n_x + v_t1 *n_x;
	v_y = v_n2*n_y + v_t1 *n_z;
	v_z = v_n2*n_z - v_t1 *n_y;

	v_x2 = v_n1*n_x + v_t2 *n_x;
	v_y2 = v_n1*n_y + v_t2 *n_z;
	v_z2 = v_n1*n_z - v_t2 *n_y;

	fsi[ibi].S_new[1]=v_x;
	fsi[ibi].S_new[3]=v_y;
	fsi[ibi].S_new[5]=v_z;

	fsi[ibi2].S_new[1]=v_x2;
	fsi[ibi2].S_new[3]=v_y2;
	fsi[ibi2].S_new[5]=v_z2;

	PetscPrintf(PETSC_COMM_WORLD, "Collision Velocity!!!! cyl1 %le %le %le  cyl2 %le %le %le\n", v_x,v_y,v_z,v_x2,v_y2,v_z2);

	
      }
    }
  } 

  if (collision){
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi],ibi);
    }
    
    for (ibi=0;ibi<NumberOfBodies;ibi++) {	  
      
      for (i=0; i<ibm[ibi].n_v; i++) {
	ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
	ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
	ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];
	
	ibm[ibi].urm1[i].x = ibm[ibi].u[i].x;
	ibm[ibi].urm1[i].y = ibm[ibi].u[i].y;
	ibm[ibi].urm1[i].z = ibm[ibi].u[i].z;
	
	ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
	ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
	ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
      }
      
      for (i=0;i<6;i++){
	fsi[ibi].S_realm1[i]=fsi[ibi].S_new[i];
	fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
	
	fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_n[i];
	fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
      }
    }
  }

  return(0);
}

/* ================================================================================== */
PetscReal Particle_Distance(IBMNodes *ibm1,IBMNodes *ibm2,FSInfo *fsi1,FSInfo *fsi2,PetscReal *dir_min,PetscReal *R1,PetscReal *R2)
{
/*   DM	        da = user->da; */
/*   IBMListNode   *current; */
/*   IBMInfo       *ibminfo; */
  PetscReal     D_min,Dist,D_MIN=0.0;
  PetscInt      i,j;
  PetscInt      cell1,cell2;
  struct POINT{
    PetscReal x;
    PetscReal y;
    PetscReal z;
  } Point1,Point2,p1,p2;
  struct {
    PetscReal value;
    PetscInt  index;
  }in,out;
  Point1.x=0.0;
  Point1.y=0.0;
  Point1.z=0.0;
  Point2.x=0.0;
  Point2.y=0.0;
  Point2.z=0.0;
  p1.x=0.0;
  p1.y=0.0;
  p1.z=0.0;
  p2.x=0.0;
  p1.y=0.0;
  p1.z=0.0;
  PetscInt  rank;
  //  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  //  DMDAVecGetArray(da, user->lNvert, &nvert);

  D_min=1.e6;
  Dist=0.0;

  for (i=0;i<ibm1->n_v;i++){
    for (j=0;j<ibm2->n_v;j++){
      p1.x=ibm1->x_bp[i];
      p1.y=ibm1->y_bp[i];
      p1.z=ibm1->z_bp[i];
      p2.x=ibm2->x_bp[j];
      p2.y=ibm2->y_bp[j];
      p2.z=ibm2->z_bp[j];
      Dist=sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
      if (Dist<D_min){
	D_min=Dist;
	Point1.x=p1.x;
	Point1.y=p1.y;
	Point1.z=p1.z;
	cell1=i;
	cell2=j;
	Point2.x=p2.x;
	Point2.y=p2.y;
	Point2.z=p2.z;
      }
    }
  }


  //  current = user->ibmlist[P1].head;

	
  //  particle_found=1;
  
/*   while (current) { */
/*     ibminfo = &current->ibm_intp; */
/*     current = current->next; */
/*     i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk; */
/*    /\*  if (fabs(nvert[k][j][i]-1.0-P1/1000.0)>1.e-5) { *\/ */
/* /\*       PetscPrintf(PETSC_COMM_WORLD, "Partciles %d amd %d are inside each other i %d j %d k %d \n",P1,P2,i,j,k); *\/ */
/* /\*       goto nextp1; *\/ */
/* /\*     } *\/ */
    
/*     for (kk=0; kk<5; kk++){ */
/*       for (jj=0; jj<5; jj++){ */
/* 	for (ii=0; ii<5; ii++){ */
	  
/* 	  if (P2==(int)((nvert[k+kk-2][j+jj-2][i+ii-2]-1.000)*1001) && fabs(nvert[k][j][i]-nvert[k+kk-2][j+jj-2][i+ii-2])>1.e-5) { */
	    
/* 	    cell1=ibminfo->cell; */
	    
/* 	  /\*   p1.x=user->ibm[P1].cent_x[cell1]; *\/ */
/* /\* 	    p1.y=user->ibm[P1].cent_y[cell1]; *\/ */
/* /\* 	    p1.z=user->ibm[P1].cent_z[cell1]; *\/ */
/* 	    //  m=ibm[P2].n_v-1; */
/* 	    //   m=3868; */
/* 	    p1.x=ibminfo->pmin.x; */
/* 	    p1.y=ibminfo->pmin.y; */
/* 	    p1.z=ibminfo->pmin.z; */
/* 	    /\*   p2.x=ibm[P2].x_bp[m]; *\/ */
/* 	    /\* 	    p2.y=ibm[P2].y_bp[m]; *\/ */
/* 	    /\* 	    p2.z=ibm[P2].z_bp[m]; *\/ */
/* 	    for (n=0;n<ibm[P2].n_v;n++){ */
/* 	      p2.x=ibm[P2].x_bp[n]; */
/* 	      p2.y=ibm[P2].y_bp[n]; */
/* 	      p2.z=ibm[P2].z_bp[n]; */
/* 	      Dist=sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z)); */
/* 	      if (Dist<D_min){ */
/* 		D_min=Dist; */
/* 		Point1.x=p1.x; */
/* 		Point1.y=p1.y; */
/* 		Point1.z=p1.z; */
/* 		cell2=n; */
/* 		Point2.x=p2.x; */
/* 		Point2.y=p2.y; */
/* 		Point2.z=p2.z; */
/* 	      } */
/* 	    } */
	  
/* 	    *particle_found ++; */
	    
/* 	  } */
/* 	} */
/*       } */
/*     } */
    
/*   } */
//  PetscPrintf(PETSC_COMM_WORLD, "partcile found  %d \n",*particle_found);
  //  PetscPrintf(PETSC_COMM_WORLD, "partciles closest distance : %le cell-1 %d p1-x %le p1-y %le p1-z %le cell-2 %d p2-x %le p2-y %le p2-z %le \n",D_min,cell1,Point1.x,Point1.y,Point1.z,cell2,Point2.x,Point2.y,Point2.z);

  
  /*   PetscPrintf(PETSC_COMM_WORLD, "cell 1  %d x %le y %le z %le cell 2  %d  x %le y %le z %le \n",cell1,Point1.x,Point1.y,Point1.z,cell2,Point2.x,Point2.y,Point2.z); */

  in.value=D_min;
  in.index=rank;
  MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,0,PETSC_COMM_WORLD);
  
  MPI_Allreduce(&D_min, &D_MIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);
  //  MPI_Allreduce(&particle_found, &particle_found,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
  
  PetscPrintf(PETSC_COMM_WORLD, "closest distance between paricles is  : %le \n",D_MIN);
  // else PetscPrintf(PETSC_COMM_WORLD, "particles are out of threshold \n");
  PetscReal dir_mag=0.0;

  
  dir_min[0]=(Point2.x-Point1.x);
  dir_min[1]=(Point2.y-Point1.y);
  dir_min[2]=(Point2.z-Point1.z);
  dir_mag=sqrt((Point1.x-Point2.x)*(Point1.x-Point2.x)+(Point1.y-Point2.y)*(Point1.y-Point2.y)+(Point1.z-Point2.z)*(Point1.z-Point2.z));
  dir_min[0] = dir_min[0]/dir_mag;
  dir_min[1] = dir_min[1]/dir_mag;
  dir_min[2] = dir_min[2]/dir_mag;
  
  R1[0]=(Point1.x-fsi1->x_c);
  R1[1]=(Point1.y-fsi1->y_c);
  R1[2]=(Point1.z-fsi1->z_c);
/*   R1[0]=(Point1.x-fsi[P1].x_c); */
/*   R1[1]=(Point1.y-fsi[P1].y_c); */
/*   R1[2]=(Point1.z-fsi[P1].z_c); */
  
  R2[0]=(Point2.x-fsi2->x_c);
  R2[1]=(Point2.y-fsi2->y_c);
  R2[2]=(Point2.z-fsi2->z_c);
/*   R2[0]=(Point2.x-fsi[P2].x_c); */
/*   R2[1]=(Point2.y-fsi[P2].y_c); */
/*   R2[2]=(Point2.z-fsi[P2].z_c); */
  
  int root=out.index;
  MPI_Bcast(&root,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  
  MPI_Bcast(dir_min,3,MPI_DOUBLE,root,PETSC_COMM_WORLD);
  MPI_Bcast(R1,3,MPI_DOUBLE,root,PETSC_COMM_WORLD);
  MPI_Bcast(R2,3,MPI_DOUBLE,root,PETSC_COMM_WORLD);
  
  PetscPrintf(PETSC_COMM_WORLD, "min_dir------n_x %le n_y %le n_z %le \n",dir_min[0],dir_min[1],dir_min[2]);

  
  
 
  // nextp1: DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  return D_MIN;
}
/* ==================================================================================             */

PetscInt Bounding_Box(IBMNodes *ibm1,IBMNodes *ibm2)
{

  PetscInt mode=0;
  PetscReal	  *x_bp1 =ibm1->x_bp, *y_bp1 = ibm1->y_bp, *z_bp1 = ibm1->z_bp;
  PetscReal	xbp_min_1, ybp_min_1, zbp_min_1, xbp_max_1, ybp_max_1, zbp_max_1;
  PetscReal	  *x_bp2 =ibm2->x_bp, *y_bp2 = ibm2->y_bp, *z_bp2 = ibm2->z_bp;
  PetscReal	xbp_min_2, ybp_min_2, zbp_min_2, xbp_max_2, ybp_max_2, zbp_max_2;
  PetscInt      i,n_v;

  xbp_min_1 = 1.e23; xbp_max_1 = -1.e23;
  ybp_min_1 = 1.e23; ybp_max_1 = -1.e23;
  zbp_min_1 = 1.e23; zbp_max_1 = -1.e23;
  
  xbp_min_2 = 1.e23; xbp_max_2 = -1.e23;
  ybp_min_2 = 1.e23; ybp_max_2 = -1.e23;
  zbp_min_2 = 1.e23; zbp_max_2 = -1.e23;

  n_v = ibm1->n_v;
  for(i=0; i<n_v; i++) {
      
      xbp_min_1 = PetscMin(xbp_min_1, x_bp1[i]);
      xbp_max_1 = PetscMax(xbp_max_1, x_bp1[i]);
      
      ybp_min_1 = PetscMin(ybp_min_1, y_bp1[i]);
      ybp_max_1 = PetscMax(ybp_max_1, y_bp1[i]);
      
      zbp_min_1 = PetscMin(zbp_min_1, z_bp1[i]);
      zbp_max_1 = PetscMax(zbp_max_1, z_bp1[i]);
    }

    xbp_min_1 -= 0.05; xbp_max_1 += 0.05;
    ybp_min_1 -= 0.05; ybp_max_1 += 0.05;
    zbp_min_1 -= 0.05; zbp_max_1 += 0.05;

    n_v = ibm2->n_v;
    for(i=0; i<n_v; i++) {
      
      xbp_min_2 = PetscMin(xbp_min_2, x_bp2[i]);
      xbp_max_2 = PetscMax(xbp_max_2, x_bp2[i]);
      
      ybp_min_2 = PetscMin(ybp_min_2, y_bp2[i]);
      ybp_max_2 = PetscMax(ybp_max_2, y_bp2[i]);
      
      zbp_min_2 = PetscMin(zbp_min_2, z_bp2[i]);
      zbp_max_2 = PetscMax(zbp_max_2, z_bp2[i]);
    }

    xbp_min_2 -= 0.05; xbp_max_2 += 0.05;
    ybp_min_2 -= 0.05; ybp_max_2 += 0.05;
    zbp_min_2 -= 0.05; zbp_max_2 += 0.05;

    PetscPrintf(PETSC_COMM_WORLD, "Bounding Box1 --- x_min %le x_max %le y_min %le y_max %le z_min %le z_max %le \n", xbp_min_1, xbp_max_1,ybp_min_1,ybp_max_1, zbp_min_1, zbp_max_1);
    PetscPrintf(PETSC_COMM_WORLD, "Bounding Box2 --- x_min %le x_max %le y_min %le y_max %le z_min %le z_max %le \n", xbp_min_2, xbp_max_2,ybp_min_2,ybp_max_2, zbp_min_2, zbp_max_2);
   
    if ((xbp_min_1 >= xbp_min_2 && xbp_min_1 <= xbp_max_2) || (xbp_max_1 >= xbp_min_2 && xbp_max_1 <= xbp_max_2) ||
	(xbp_min_2 >= xbp_min_1 && xbp_min_2 <= xbp_max_1) || (xbp_max_2 >= xbp_min_1 && xbp_max_2 <= xbp_max_1)) mode++;
    if ((ybp_min_1 >= ybp_min_2 && ybp_min_1 <= ybp_max_2) || (ybp_max_1 >= ybp_min_2 && ybp_max_1 <= ybp_max_2) ||
	(ybp_min_2 >= ybp_min_1 && ybp_min_2 <= ybp_max_1) || (ybp_max_2 >= ybp_min_1 && ybp_max_2 <= ybp_max_1)) mode++;
    if ((zbp_min_1 >= zbp_min_2 && zbp_min_1 <= zbp_max_2) || (zbp_max_1 >= zbp_min_2 && zbp_max_1 <= zbp_max_2) ||
	(zbp_min_2 >= zbp_min_1 && zbp_min_2 <= zbp_max_1) || (zbp_max_2 >= zbp_min_1 && zbp_max_2 <= zbp_max_1)) mode++;
 
  return mode;
}


/* ==================================================================================== */
PetscErrorCode Particle_Rot_Matrix(FSInfo *fsi,IBMNodes *ibm, PetscReal *Dw, PetscReal dt)

{
  PetscInt i;
 /*  PetscReal R[3][3],RT[3][3],X[3][3],I[3][3],I_inv[3][3],R1[3][3],R2[3][3],q_n[4],Inn,det; */
/*   PetscReal Dw_mag,theta_n,n[3]; */
/*   PetscInt  n_v=ibm->n_v; */
/*   Cmpnts a_c,p; */

 /*  for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       R1[i][j]=fsi->R[i][j]; */
/*     } */
/*   } */
  
 /*  for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       if (i==j) RT[i][j]=R1[i][j]; */
/*       else      RT[i][j]=R1[j][i]; */
/*     } */
/*   } */
  
/*   for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       X[i][j]=fsi->I_inv[i][0]*RT[0][j]+fsi->I_inv[i][1]*RT[1][j]+fsi->I_inv[i][2]*RT[2][j]; */
/*     } */
/*   } */
  
/*   for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       I_inv[i][j]=R1[i][0]*X[0][j]+R1[i][1]*X[1][j]+R1[i][2]*X[2][j]; */
/*     } */
/*   } */
  
/*   det=I_inv[0][0]*I_inv[1][1]*I_inv[2][2]+2*I_inv[0][1]*I_inv[0][2]*I_inv[1][2]-I_inv[0][0]*I_inv[1][2]*I_inv[1][2]-I_inv[1][1]*I_inv[0][2]*I_inv[0][2]-I_inv[2][2]*I_inv[0][1]*I_inv[0][1]; */
  
/*   I[0][0]=(I_inv[1][1]*I_inv[2][2]-I_inv[1][2]*I_inv[1][2])/det; */
/*   I[0][1]=(I_inv[0][2]*I_inv[1][2]-I_inv[0][1]*I_inv[2][2])/det; */
/*   I[0][2]=(I_inv[0][1]*I_inv[1][2]-I_inv[0][2]*I_inv[1][1])/det; */
/*   I[1][0]=I[0][1]; */
/*   I[1][1]=(I_inv[0][0]*I_inv[2][2]-I_inv[0][2]*I_inv[0][2])/det; */
/*   I[1][2]=(I_inv[0][1]*I_inv[0][2]-I_inv[0][0]*I_inv[1][2])/det; */
/*   I[2][0]=I[0][2]; */
/*   I[2][1]=I[1][2]; */
/*   I[2][2]=(I_inv[0][0]*I_inv[1][1]-I_inv[0][1]*I_inv[0][1])/det; */

 /*  Dw_mag=sqrt(Dw[0]*Dw[0]+Dw[1]*Dw[1]+Dw[2]*Dw[2]); */
 
/*   n[0]=Dw[0]/Dw_mag; */
/*   n[1]=Dw[1]/Dw_mag; */
/*   n[2]=Dw[2]/Dw_mag; */


  // Inn=I[0][0]*n[0]*n[0]+2*I[0][1]*n[0]*n[1]+2*I[0][2]*n[0]*n[2]+I[1][1]*n[1]*n[1]+2*I[1][2]*n[1]*n[2]+I[2][2]*n[2]*n[2];

  //  *w_nn=M_mag*dt/Inn;
 
 /*  theta_n=0.5*dt*Dw_mag; */
/*   PetscPrintf(PETSC_COMM_WORLD, "Dw ------%le --------theta_n-------%le \n",Dw_mag,theta_n); */

/*   R2[0][0]=cos(theta_n)+n[0]*n[0]*(1-cos(theta_n)); */
/*   R2[0][1]=n[0]*n[1]*(1-cos(theta_n))-n[2]*sin(theta_n); */
/*   R2[0][2]=n[0]*n[2]*(1-cos(theta_n))+n[1]*sin(theta_n); */

/*   R2[1][0]=n[0]*n[1]*(1-cos(theta_n))+n[2]*sin(theta_n); */
/*   R2[1][1]=cos(theta_n)+n[1]*n[1]*(1-cos(theta_n)); */
/*   R2[1][2]=n[1]*n[2]*(1-cos(theta_n))-n[0]*sin(theta_n); */

/*   R2[2][0]=n[0]*n[2]*(1-cos(theta_n))-n[1]*sin(theta_n); */
/*   R2[2][1]=n[1]*n[2]*(1-cos(theta_n))+n[0]*sin(theta_n); */
/*   R2[2][2]=cos(theta_n)+n[2]*n[2]*(1-cos(theta_n)); */
  PetscReal A[4],r[4];
  PetscReal    w_r[3],w_n[3],q_n[4],q_r[4];

  w_r[0]=fsi->S_ang_r[1];
  w_r[1]=fsi->S_ang_r[3];
  w_r[2]=fsi->S_ang_r[5];

  w_n[0]=fsi->S_ang_n[1]+Dw[0];
  w_n[1]=fsi->S_ang_n[3]+Dw[1];
  w_n[2]=fsi->S_ang_n[5]+Dw[2];
 
  for (i=0;i<4;i++){
    q_r[i]=fsi->q_r[i];
  }

  A[0]=1/dt;
  A[1]=0.25*w_n[0];
  A[2]=0.25*w_n[1];
  A[3]=0.25*w_n[2];

 /*  A[1][0]=-A[0][1]; */
/*   A[1][1]=A[0][0]; */
/*   A[1][2]=A[0][3]; */
/*   A[1][3]=-A[0][2]; */

/*   A[2][0]=-A[0][2]; */
/*   A[2][1]=-A[0][3]; */
/*   A[2][2]=A[0][0]; */
/*   A[2][3]=A[0][1]; */

/*   A[3][0]=-A[0][3]; */
/*   A[3][1]=A[0][2]; */
/*   A[3][2]=-A[0][1]; */
/*   A[3][3]=A[0][0]; */

  r[0]=q_r[0]/dt+0.25*(-w_r[0]*q_r[1]-w_r[1]*q_r[2]-w_r[2]*q_r[3]);
  r[1]=q_r[1]/dt+0.25*( w_r[0]*q_r[0]+w_r[1]*q_r[3]-w_r[2]*q_r[2]);
  r[2]=q_r[2]/dt+0.25*(-w_r[0]*q_r[3]+w_r[1]*q_r[0]+w_r[2]*q_r[1]);
  r[3]=q_r[3]/dt+0.25*( w_r[0]*q_r[2]-w_r[1]*q_r[1]+w_r[2]*q_r[0]);
  
 /*  for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       R[i][j]=R2[i][0]*R1[0][j]+R2[i][1]*R1[1][j]+R2[i][2]*R1[2][j]; */
/*     } */
/*   } */

/*   q_n[0]=sqrt((R[0][0]+R[1][1]+R[2][2]+1.0)/4.0); */
/*   q_n[1]=(R[2][1]-R[1][2])/(4*q_n[0]); */
/*   q_n[2]=(R[1][0]+R[0][1])/(4*q_n[1]); */
/*   q_n[3]=(R[0][2]+R[2][0])/(4*q_n[1]); */

  q_n[0]=(A[0]*r[0]-A[1]*r[1]-A[2]*r[2]-A[3]*r[3])/(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);
  q_n[1]=(A[1]*r[0]+A[0]*r[1]-A[3]*r[2]+A[2]*r[3])/(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);
  q_n[2]=(A[2]*r[0]+A[3]*r[1]+A[0]*r[2]-A[1]*r[3])/(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);
  q_n[3]=(A[3]*r[0]-A[2]*r[1]+A[1]*r[2]+A[0]*r[3])/(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);

  PetscPrintf(PETSC_COMM_WORLD, "Dq[0] %le Dq[1] %le Dq[2] %le Dq[3] %le \n",q_n[0]-fsi->q[0],q_n[1]-fsi->q[1],q_n[2]-fsi->q[2],q_n[3]-fsi->q[3]);

  for (i=0;i<4;i++){
    fsi->q[i]=q_n[i];
  }

  fsi->R[0][0]=1.0-2.0*q_n[2]*q_n[2]-2.0*q_n[3]*q_n[3];
  fsi->R[0][1]=2.0*(q_n[1]*q_n[2]-q_n[0]*q_n[3]);
  fsi->R[0][2]=2.0*(q_n[1]*q_n[3]+q_n[0]*q_n[2]);
  fsi->R[1][0]=2.0*(q_n[1]*q_n[2]+q_n[0]*q_n[3]);
  fsi->R[1][1]=1.0-2.0*q_n[1]*q_n[1]-2.0*q_n[3]*q_n[3];
  fsi->R[1][2]=2.0*(q_n[2]*q_n[3]-q_n[0]*q_n[1]);
  fsi->R[2][0]=2.0*(q_n[1]*q_n[3]-q_n[0]*q_n[2]);
  fsi->R[2][1]=2.0*(q_n[2]*q_n[3]+q_n[0]*q_n[1]);
  fsi->R[2][2]=1.0-2.0*q_n[1]*q_n[1]-2.0*q_n[2]*q_n[2];

 /*  for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       fsi->R[i][j]=R[i][j]; */
/*     } */
/*   } */

  return 0;
}

/* ==================================================================================== */
PetscErrorCode CollisionDetectionOfParticles(UserCtx *user,IBMNodes *ibm,FSInfo *fsi)

{
 
/*   DM	da = user->da; */
  PetscInt	P1,P2,ibi=0,ibj=0;
  PetscInt      Indx;
 
  PetscReal     dir_min[3],R1[3],R2[3],U_w1[3],U_w2[3],Dw1[3],Dw2[3];
  PetscReal     U_tr1[3],U_tr2[3],uw1d_o=0.0,uw2d_o=0.0,uw1d=0.0,uw2d=0.0,Duwd=0.0,u1d_o=0.0,u2d_o=0.0,u1d=0.0,u2d=0.0,dt=user->dt,Dist=0.0;
  Cmpnts        A,B,C;
  for (ibi=0;ibi<NumberOfBodies-1;ibi++){
    for (ibj=ibi+1; ibj<NumberOfBodies;ibj++){
      
      Indx=Bounding_Box(&ibm[ibi],&ibm[ibj]);
      
      if (Indx==3){
	PetscPrintf(PETSC_COMM_WORLD, "mode of partciles  %d and %d is %d \n",ibi,ibj,Indx);
	P1=ibi;
	P2=ibj;
	
	//	Particle_Distance(user,ibm,fsi,P2,P1,dir_min,R1,R2,&particle);
      
	Dist=Particle_Distance(&ibm[P1],&ibm[P2],&fsi[P1],&fsi[P2],dir_min,R1,R2); 

	if (Dist<0.05){
	  PetscReal omega[3];
	  
	  omega[0]=fsi[P1].S_ang_n[1];
	  omega[1]=fsi[P1].S_ang_n[3];
	  omega[2]=fsi[P1].S_ang_n[5];
	  //	PetscPrintf(PETSC_COMM_WORLD, "w1------ %le  %le %le \n",omega[0],omega[1],omega[2]);
	  CROSS(U_w1,omega,R1);
	  
	  U_tr1[0]=fsi[P1].S_new[1];
	  U_tr1[1]=fsi[P1].S_new[3];
	  U_tr1[2]=fsi[P1].S_new[5];
	  
	  omega[0]=fsi[P2].S_ang_n[1];
	  omega[1]=fsi[P2].S_ang_n[3];
	  omega[2]=fsi[P2].S_ang_n[5];
	  
	  CROSS(U_w2,omega,R2);
	  
	  U_tr2[0]=fsi[P2].S_new[1];
	  U_tr2[1]=fsi[P2].S_new[3];
	  U_tr2[2]=fsi[P2].S_new[5];
	  
	  //	PetscPrintf(PETSC_COMM_WORLD, "U_w1------ %le  %le  %le U_w2------ %le  %le  %le \n",U_w1[0],U_w1[1],U_w1[2],U_w2[0],U_w2[1],U_w2[2]);
	  
	  u1d_o=DOT(U_w1,dir_min)+DOT(U_tr1,dir_min);
	  u2d_o=DOT(U_w2,dir_min)+DOT(U_tr2,dir_min);
	  
	  if((u2d_o-u1d_o)<0.0) {
	    
	    PetscPrintf(PETSC_COMM_WORLD, "partciles are approaching: u1d %le u2d %le \n",u1d_o,u2d_o);
	   
	    uw1d_o=DOT(U_w1,dir_min);
	    uw2d_o=DOT(U_w2,dir_min);
	    
	    Duwd=uw1d_o-uw2d_o;
	   
	    if (Duwd > 1.e-3){

	      PetscPrintf(PETSC_COMM_WORLD,"-----repulsive moment is added \n");
	      
	      B.x=R1[0];B.y=R1[1];B.z=R1[2];
	     
	      if (uw1d_o*uw2d_o < 0.0){

		C.x= -uw1d_o*dir_min[0]; C.y= -uw1d_o*dir_min[1]; C.z= -uw1d_o*dir_min[2];

	      }else{
	     
		C.x=-0.5*Duwd*dir_min[0]; C.y=-0.5*Duwd*dir_min[1]; C.z=-0.5*Duwd*dir_min[2];
	      }

	      A.x=(B.y*C.z-B.z*C.y)/(B.x*B.x+B.y*B.y+B.z*B.z);
	      A.y=(B.z*C.x-B.x*C.z)/(B.x*B.x+B.y*B.y+B.z*B.z);
	      A.z=(B.x*C.y-B.y*C.x)/(B.x*B.x+B.y*B.y+B.z*B.z);
	      
	      Dw1[0]=A.x; Dw1[1]=A.y; Dw1[2]=A.z;
	    

	    /*   if (Dw1[0]*Dw1[0]+Dw1[1]*Dw1[1]+Dw1[2]*Dw1[2]> 0.25*(fsi[P1].S_ang_r[1]* fsi[P1].S_ang_r[1]+ fsi[P1].S_ang_r[3]* fsi[P1].S_ang_r[3]+ fsi[P1].S_ang_r[5]* fsi[P1].S_ang_r[5])){ */

/* 		Dw1[0] =A.x*0.5*sqrt(fsi[P1].S_ang_r[1]* fsi[P1].S_ang_r[1]+ fsi[P1].S_ang_r[3]* fsi[P1].S_ang_r[3]+ fsi[P1].S_ang_r[5]* fsi[P1].S_ang_r[5])/sqrt(A.x*A.x+A.y*A.y+A.z*A.z); */
/* 		Dw1[1] =A.y*0.5*sqrt(fsi[P1].S_ang_r[1]* fsi[P1].S_ang_r[1]+ fsi[P1].S_ang_r[3]* fsi[P1].S_ang_r[3]+ fsi[P1].S_ang_r[5]* fsi[P1].S_ang_r[5])/sqrt(A.x*A.x+A.y*A.y+A.z*A.z); */
/* 		Dw1[2] =A.z*0.5*sqrt(fsi[P1].S_ang_r[1]* fsi[P1].S_ang_r[1]+ fsi[P1].S_ang_r[3]* fsi[P1].S_ang_r[3]+ fsi[P1].S_ang_r[5]* fsi[P1].S_ang_r[5])/sqrt(A.x*A.x+A.y*A.y+A.z*A.z); */
/* 		PetscPrintf(PETSC_COMM_WORLD,"Dw1 is modified %le %le  \n",sqrt(A.x*A.x+A.y*A.y+A.z*A.z),0.5*sqrt(fsi[P1].S_ang_r[1]* fsi[P1].S_ang_r[1]+ fsi[P1].S_ang_r[3]* fsi[P1].S_ang_r[3]+ fsi[P1].S_ang_r[5]* fsi[P1].S_ang_r[5])); */
/* 	      } */

	      Particle_Rot_Matrix(&fsi[P1], &ibm[P1], Dw1,dt);
	    
	      fsi[P1].S_ang_n[1] += Dw1[0];
	      fsi[P1].S_ang_n[3] += Dw1[1];
	      fsi[P1].S_ang_n[5] += Dw1[2];
	 	    
	      omega[0]=fsi[P1].S_ang_n[1];
	      omega[1]=fsi[P1].S_ang_n[3];
	      omega[2]=fsi[P1].S_ang_n[5];
	      
	      CROSS(U_w1,omega,R1);
	      uw1d=DOT(U_w1,dir_min);
	    
	      B.x=R2[0];B.y=R2[1];B.z=R2[2];

	      if (uw1d_o*uw2d_o < 0.0){
		
		C.x= -uw2d_o*dir_min[0]; C.y= -uw2d_o*dir_min[1]; C.z= -uw2d_o*dir_min[2];
		
	      }else{
		
		C.x=0.5*Duwd*dir_min[0]; C.y=0.5*Duwd*dir_min[1]; C.z=0.5*Duwd*dir_min[2];
		
	      }
	     
	      A.x=(B.y*C.z-B.z*C.y)/(B.x*B.x+B.y*B.y+B.z*B.z);
	      A.y=(B.z*C.x-B.x*C.z)/(B.x*B.x+B.y*B.y+B.z*B.z);
	      A.z=(B.x*C.y-B.y*C.x)/(B.x*B.x+B.y*B.y+B.z*B.z);
	    	    
	      Dw2[0]=A.x; Dw2[1]=A.y; Dw2[2]=A.z;
	    	     
	  /*     if (Dw2[0]*Dw2[0]+Dw2[1]*Dw2[1]+Dw2[2]*Dw2[2]> 0.25*(fsi[P2].S_ang_r[1]* fsi[P2].S_ang_r[1]+ fsi[P2].S_ang_r[3]* fsi[P2].S_ang_r[3]+ fsi[P2].S_ang_r[5]* fsi[P2].S_ang_r[5])){ */

/* 		Dw2[0] =A.x*0.5*sqrt(fsi[P2].S_ang_r[1]* fsi[P2].S_ang_r[1]+ fsi[P2].S_ang_r[3]* fsi[P2].S_ang_r[3]+ fsi[P2].S_ang_r[5]* fsi[P2].S_ang_r[5])/sqrt(A.x*A.x+A.y*A.y+A.z*A.z); */
/* 		Dw2[1] =A.y*0.5*sqrt(fsi[P2].S_ang_r[1]* fsi[P2].S_ang_r[1]+ fsi[P2].S_ang_r[3]* fsi[P2].S_ang_r[3]+ fsi[P2].S_ang_r[5]* fsi[P2].S_ang_r[5])/sqrt(A.x*A.x+A.y*A.y+A.z*A.z); */
/* 		Dw2[2] =A.z*0.5*sqrt(fsi[P2].S_ang_r[1]* fsi[P2].S_ang_r[1]+ fsi[P2].S_ang_r[3]* fsi[P2].S_ang_r[3]+ fsi[P2].S_ang_r[5]* fsi[P2].S_ang_r[5])/sqrt(A.x*A.x+A.y*A.y+A.z*A.z); */
/* 		PetscPrintf(PETSC_COMM_WORLD,"Dw2 is modified %le %le \n",sqrt(A.x*A.x+A.y*A.y+A.z*A.z),0.5*sqrt(fsi[P2].S_ang_r[1]* fsi[P2].S_ang_r[1]+ fsi[P2].S_ang_r[3]* fsi[P2].S_ang_r[3]+ fsi[P2].S_ang_r[5]* fsi[P2].S_ang_r[5])); */
/* 	      } */
	    
	      Particle_Rot_Matrix(&fsi[P2], &ibm[P2], Dw2,dt);
	    
	      fsi[P2].S_ang_n[1] += Dw2[0];
	      fsi[P2].S_ang_n[3] += Dw2[1];
	      fsi[P2].S_ang_n[5] += Dw2[2];
	    
	      omega[0]=fsi[P2].S_ang_n[1];
	      omega[1]=fsi[P2].S_ang_n[3];
	      omega[2]=fsi[P2].S_ang_n[5];
	      
	      CROSS(U_w2,omega,R2);
	      uw2d=DOT(U_w2,dir_min);
	    
	      u1d=DOT(U_w1,dir_min)+DOT(U_tr1,dir_min);
	      u2d=DOT(U_w2,dir_min)+DOT(U_tr2,dir_min);
	     

	      if((u2d-u1d)<0.0) {
		PetscPrintf(PETSC_COMM_WORLD, "u1d %le u2d %le -----translational repulsive is added \n",u1d,u2d);   
	      
		if (u1d*u2d<0.0) {
	     
		  fsi[P1].S_new[1] -=u1d*dir_min[0];
		  fsi[P1].S_new[3] -=u1d*dir_min[1];
		  fsi[P1].S_new[5] -=u1d*dir_min[2];
		  
		  fsi[P2].S_new[1] -=u2d*dir_min[0];
		  fsi[P2].S_new[3] -=u2d*dir_min[1];
		  fsi[P2].S_new[5] -=u2d*dir_min[2];
		}else{
		  fsi[P1].S_new[1] -=0.5*(u1d-u2d)*dir_min[0];
		  fsi[P1].S_new[3] -=0.5*(u1d-u2d)*dir_min[1];
		  fsi[P1].S_new[5] -=0.5*(u1d-u2d)*dir_min[2];
		  
		  fsi[P2].S_new[1] +=0.5*(u1d-u2d)*dir_min[0];
		  fsi[P2].S_new[3] +=0.5*(u1d-u2d)*dir_min[1];
		  fsi[P2].S_new[5] +=0.5*(u1d-u2d)*dir_min[2];
		}
		
	      } 
	    }else{
	      PetscPrintf(PETSC_COMM_WORLD, "Duwd %le -----translational repulsive is added \n",Duwd);

	      if (u1d_o*u2d_o<0.0) {
	      
		fsi[P1].S_new[1] -=u1d_o*dir_min[0];
		fsi[P1].S_new[3] -=u1d_o*dir_min[1];
		fsi[P1].S_new[5] -=u1d_o*dir_min[2];
		
		fsi[P2].S_new[1] -=u2d_o*dir_min[0];
		fsi[P2].S_new[3] -=u2d_o*dir_min[1];
		fsi[P2].S_new[5] -=u2d_o*dir_min[2];
	      }else{
		fsi[P1].S_new[1] -=0.5*(u1d_o-u2d_o)*dir_min[0];
		fsi[P1].S_new[3] -=0.5*(u1d_o-u2d_o)*dir_min[1];
		fsi[P1].S_new[5] -=0.5*(u1d_o-u2d_o)*dir_min[2];
		
		fsi[P2].S_new[1] +=0.5*(u1d_o-u2d_o)*dir_min[0];
		fsi[P2].S_new[3] +=0.5*(u1d_o-u2d_o)*dir_min[1];
		fsi[P2].S_new[5] +=0.5*(u1d_o-u2d_o)*dir_min[2];
	      }
	    }
	    
	    U_tr1[0]=fsi[P1].S_new[1];
	    U_tr1[1]=fsi[P1].S_new[3];
	    U_tr1[2]=fsi[P1].S_new[5];
	    
	    U_tr2[0]=fsi[P2].S_new[1];
	    U_tr2[1]=fsi[P2].S_new[3];
	    U_tr2[2]=fsi[P2].S_new[5];
	    
	    u1d=DOT(U_w1,dir_min)+DOT(U_tr1,dir_min);
	    u2d=DOT(U_w2,dir_min)+DOT(U_tr2,dir_min);
	    
	    PetscPrintf(PETSC_COMM_WORLD, " after collision stretegy --------u1d %le u2d %le \n",u1d,u2d);
	    
	  }
	}else 	PetscPrintf(PETSC_COMM_WORLD, "particles are out of influential distance %le \n",Dist);
      }
    }
  }
  
  return(0);
}


/* /\* ==================================================================================             *\/ */
/* PetscErrorCode Repulsive(FSInfo *fsi1,IBMNodes *ibm1,FSInfo *fsi2,IBMNodes *ibm2,PetscReal dt,PetscReal *dir_min,PetscReal *R1,PetscReal *R2,PetscReal *w1,PetscReal *w2,PetscReal *q1,PetscReal *q2,PetscReal uw1d_o,PetscReal uw2d_o, PetscReal epsilon) */

/* { */
 
/*   PetscReal U1[3],U2[3],C12[3],w_dir[3];  */
/*   U1[0]=fsi1->S_new[1];U1[1]=fsi1->S_new[3];U1[2]=fsi1->S_new[5]; */
/*   U2[0]=fsi2->S_new[1];U2[1]=fsi2->S_new[3];U2[2]=fsi2->S_new[5]; */
/*  /\*  PetscPrintf(PETSC_COMM_WORLD, "U1_x: %le U1_y %le U1_z %le U2_x  %le U2_y %le U2_z %le \n",U1[0],U1[1],U1[2],U2[0],U2[1],U2[2]);  *\/ */
/* /\*   PetscPrintf(PETSC_COMM_WORLD, "U_w1_x: %le U_w1_y %le U_w1_z %le U_w2_x  %le U_w2_y %le U_w2_z %le \n",U_w1[0],U_w1[1],U_w1[2],U_w2[0],U_w2[1],U_w2[2]); *\/  */
/*   PetscReal F1[3],F2[3],M1[3],M2[3],Fw[3],n[3],q_n[4],w_n=0.0; */
/*   PetscReal u1d,u2d,C12d; */
/*   PetscInt i; */

/*  /\*  u1d=DOT(U1,dir_min); *\/ */
/* /\*   u2d=DOT(U2,dir_min); *\/ */
/* /\*   if ((u2d-u1d)<0.0) { *\/ */
/* /\*     PetscPrintf(PETSC_COMM_WORLD, "repulsive force is added \n"); *\/ */
/* /\*     C12[0]=fsi1->x_c-fsi2->x_c; *\/ */
/* /\*     C12[1]=fsi1->y_c-fsi2->y_c; *\/ */
/* /\*     C12[2]=fsi1->z_c-fsi2->z_c; *\/ */
/* /\*     C12d=DOT(C12,dir_min); *\/ */
/* /\*     C12d=C12d/fabs(C12d); *\/ */
    
/* /\*     F1[0] =-C12d*fabs(u2d-u1d)*dir_min[0]*epsilon; *\/ */
/* /\*     F1[1] =-C12d*fabs(u2d-u1d)*dir_min[0]*epsilon; *\/ */
/* /\*     F1[2] =-C12d*fabs(u2d-u1d)*dir_min[0]*epsilon; *\/ */
    
/* /\*     F2[0] =C12d*fabs(u2d-u1d)*dir_min[0]*epsilon; *\/ */
/* /\*     F2[1] =C12d*fabs(u2d-u1d)*dir_min[0]*epsilon; *\/ */
/* /\*     F2[2] =C12d*fabs(u2d-u1d)*dir_min[0]*epsilon; *\/ */
/* /\*   } *\/ */

  
/*   C12[0]=fsi2->x_c-fsi1->x_c; */
/*   C12[1]=fsi2->y_c-fsi1->y_c; */
/*   C12[2]=fsi2->z_c-fsi1->z_c; */
/*   C12d=DOT(C12,dir_min); */
/*   C12d=C12d/fabs(C12d); */

/*   Fw[0] =-C12d*fabs(uw2d_o-uw1d_o)*dir_min[0]*epsilon; */
/*   Fw[1] =-C12d*fabs(uw2d_o-uw1d_o)*dir_min[1]*epsilon; */
/*   Fw[2] =-C12d*fabs(uw2d_o-uw1d_o)*dir_min[2]*epsilon; */
  
/*   CROSS(M1,R1,Fw); */
  
/*   n[0]=M1[0]/sqrt(M1[0]*M1[0]+M1[1]*M1[1]+M1[2]*M1[2]); */
/*   n[1]=M1[1]/sqrt(M1[0]*M1[0]+M1[1]*M1[1]+M1[2]*M1[2]); */
/*   n[2]=M1[2]/sqrt(M1[0]*M1[0]+M1[1]*M1[1]+M1[2]*M1[2]); */


/*   Particle_Rot_Matrix(fsi1,ibm1,ibmv1,M1,dt,q_n,&w_n); */
/*   PetscPrintf(PETSC_COMM_WORLD, "w_n for partcle 1 is  %le \n",w_n); */
/*   for (i=0;i<4;i++) q1[i]=q_n[i]; */
/*   for (i=0;i<3;i++) w1[i]=w_n*n[i]; */
  
/*   Fw[0] =C12d*fabs(uw2d_o-uw1d_o)*dir_min[0]*epsilon; */
/*   Fw[1] =C12d*fabs(uw2d_o-uw1d_o)*dir_min[1]*epsilon; */
/*   Fw[2] =C12d*fabs(uw2d_o-uw1d_o)*dir_min[2]*epsilon; */
  
/*   CROSS(M2,R2,Fw); */
  
/*   n[0]=M2[0]/sqrt(M2[0]*M2[0]+M2[1]*M2[1]+M2[2]*M2[2]); */
/*   n[1]=M2[1]/sqrt(M2[0]*M2[0]+M2[1]*M2[1]+M2[2]*M2[2]); */
/*   n[2]=M2[2]/sqrt(M2[0]*M2[0]+M2[1]*M2[1]+M2[2]*M2[2]); */

/*   Particle_Rot_Matrix(fsi2,ibm2,ibmv2,M2,dt,q_n,&w_n); */
/*   PetscPrintf(PETSC_COMM_WORLD, "w_n for partcle 2 is %le  \n",w_n); */
/*   for (i=0;i<4;i++) q2[i]=q_n[i]; */
/*   for (i=0;i<3;i++) w2[i]=w_n*n[i]; */
  



/*   return 0; */
/* } */
/* ==================================================================================             */


PetscErrorCode SwingCylinder(FSInfo *fsi, IBMNodes *ibm)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscReal  x_c=fsi->x_c, y_c=fsi->y_c, z_c=fsi->z_c;
  PetscReal  rot_y=0, rot_z=fsi->S_ang_n[4];//fsi->S_ang_n[2]
  PetscInt i;
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal rx,ry,rz;
  PetscReal wx=fsi->S_ang_n[1], wy=fsi->S_ang_n[3], wz=fsi->S_ang_n[5];
  PetscReal lrot,denom;

  denom=sqrt(1.+tan(rot_z)*tan(rot_z)+tan(rot_y)*tan(rot_y));

  PetscPrintf(PETSC_COMM_WORLD, "Body Rotate y,z: %le %le Center %le %le %le %le max z y %le %le\n",rot_y,rot_z, x_c,y_c,z_c, denom, 31/denom*tan(rot_y),-31/denom*tan(rot_z));

  for (i=0; i<n_v; i++) {          
    lrot= ibm->x_bp0[i] / denom; 
    ibm->x_bp[i] = lrot;
    ibm->z_bp[i] = ibm->z_bp0[i]+ lrot*tan(rot_y);
    ibm->y_bp[i] = ibm->y_bp0[i]- lrot*tan(rot_z);
  }

  /*   calc the new normal and assign vel */
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
    
    // ns = nf x k
    if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
	(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
      ibm->ns_x[i] = 1.;     
      ibm->ns_y[i] = 0.;     
      ibm->ns_z[i] = 0 ;

      ibm->nt_x[i] = 0.;
      ibm->nt_y[i] = 1.;
      ibm->nt_z[i] = 0.;
      } else {
      ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);      
      ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);     
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      }

  }

  // 2nd order approx. 
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      rx = ibm->x_bp[i]-x_c;
      ry = ibm->y_bp[i]-y_c;
      rz = ibm->z_bp[i]-z_c;      
      ibm->u[i].x =   ry*wz-wy*rz  ;
      ibm->u[i].y =-( rx*wz-wx*rz );
      ibm->u[i].z =   rx*wy-wx*ry  ;      
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    }
  }
  return(0);
}

/* ==================================================================================             */
/* Left Ventricle Subroutines! */
/* ==================================================================================             */

PetscErrorCode LV_beat(IBMNodes *ibm, PetscReal time
		       ,PetscReal delti, PetscInt ibi)
{
  PetscReal  pi = 3.141592653589793;
  PetscInt   i,ii, n_v=ibm->n_v, n_elmt=ibm->n_elmt;
  PetscReal  x0,x1,y0,y1,z0,z1;
  PetscReal   kwave= 2*pi/wavelength;
  PetscReal  omega=2*pi*St_exp;  // non-dim  frequency St=fL/U
  PetscReal   c2;

  Cmpnts          Zenit, Apex;
  PetscReal       Width,Length;
  Cmpnts          normal;
  PetscReal       Cylinder_Radius,tmp,distance_to_axis,angle,distance_to_Zenit;
  Cmpnts          r,d,displacement,Top;
  PetscReal       Magnitude=0.8; 

  PetscPrintf(PETSC_COMM_WORLD, "LV: %d %le %le %le %le %le %le\n",ti,St_exp, omega, wavelength, kwave, time, delti);

  y0=1.e+5;
  y1=-1.e+5;
  x0=y0;
  x1=y1;
  z0=y0;
  z1=y1;
  for (i=0; i<n_v; i++) {
    x0=PetscMin(x0,(ibm->x_bp0[i]));
    x1=PetscMax(x1,(ibm->x_bp0[i]));
    y0=PetscMin(y0,(ibm->y_bp0[i]));
    y1=PetscMax(y1,(ibm->y_bp0[i]));
    z0=PetscMin(z0,(ibm->z_bp0[i]));
    z1=PetscMax(z1,(ibm->z_bp0[i]));
  }

  PetscPrintf(PETSC_COMM_WORLD, "MAX LV: %d %le %le %le %le %le %le\n",ti,z1,z0,y1,y0,x1,x0);

/*   if (LV) { */
    Top.x = -0.6024+CMx_c;
    Top.y = -0.8051+CMy_c;
    Top.z = -0.052 +CMz_c;
    //Set the axis of the contraction
    
    Apex.x = 1.86095 +CMx_c;
    Apex.y = 1.32064 +CMy_c;
    Apex.z = -1.84186+CMz_c;
    
    //Size of the box
    Width = 3;
    Length = 2.4;
    
    //normal 
    
    normal.x = Top.x - Apex.x;
    normal.y = Top.y - Apex.y;
    normal.z = Top.z - Apex.z;
    
    tmp = sqrt( normal.x * normal.x + 
		normal.y * normal.y + 
		normal.z * normal.z);
    
    normal.x /= tmp;
    normal.y /= tmp;
    normal.z /= tmp;
    
    
    Zenit.x = Apex.x + Length*normal.x;
    Zenit.y = Apex.y + Length*normal.y;
    Zenit.z = Apex.z + Length*normal.z;
    
    
    //Define the cylinder radius
    Cylinder_Radius = 2.5;
    
    //Remember the old position
 /*    for (i=0; i<libm->n_v; i++)  */
/*       { */
/* 	libm->x_bp_o[i] = libm->x_bp[i]; */
/* 	libm->y_bp_o[i] = libm->y_bp[i]; */
/* 	libm->z_bp_o[i] = libm->z_bp[i]; */
/*       } */
    
    
    //Move it    
    for (i=0; i<ibm->n_v; i++) 
      {
	
	//Test inside the cylinder
	r.x = ibm->x_bp0[i] - Apex.x;
	r.y = ibm->y_bp0[i] - Apex.y;
	r.z = ibm->z_bp0[i] - Apex.z;
	
	
	tmp = r.x * normal.x + r.y * normal.y + r.z * normal.z;
	
	d.x = tmp * normal.x;
	d.y = tmp * normal.y;
	d.z = tmp * normal.z;
	
	distance_to_axis = sqrt(r.x * r.x + r.y * r.y + r.z *r.z - 
				d.x * d.x - d.y * d.y - d.z * d.z);
	
	//Zenit computation
	r.x = ibm->x_bp0[i] - Zenit.x;
	r.y = ibm->y_bp0[i] - Zenit.y;
	r.z = ibm->z_bp0[i] - Zenit.z;
	
	angle = r.x * normal.x +
	  r.y * normal.y + 
	  r.z * normal.z;
	
	distance_to_Zenit = fabs(r.x * normal.x + r.y * normal.y + r.z * normal.z);
	
	
	if (distance_to_axis < Cylinder_Radius && angle < 0)
	  {
	    r.x = ibm->x_bp0[i] - Apex.x;
	    r.y = ibm->y_bp0[i] - Apex.y;
	    r.z = ibm->z_bp0[i] - Apex.z; 

	    //Displacement direction
	    displacement.x = d.x - r.x;
	    displacement.y = d.y - r.y;
	    displacement.z = d.z - r.z;
	    
	    tmp=sqrt(displacement.x * displacement.x + 
		     displacement.y * displacement.y + 
		     displacement.z * displacement.z);
	    
	    displacement.x /= tmp;
	    displacement.y /= tmp;
	    displacement.z /= tmp;

/* 	    distance_to_Zenit=1.; */
/* 	    distance_to_axis=1.; */
	    c2=0.5*omega*time*sin(omega*time);
	    
/* 	    ibm->x_bp[i] = ibm->x_bp0[i] + c2*displacement.x * sin(omega*time)*distance_to_axis*distance_to_Zenit*Magnitude; */
	    
/* 	    ibm->y_bp[i] = ibm->y_bp0[i] + c2*displacement.y * sin(omega*time)*distance_to_axis*distance_to_Zenit*Magnitude; */
	   
/* 	    ibm->z_bp[i] = ibm->z_bp0[i] + c2*displacement.z * sin(omega*time)*distance_to_axis*distance_to_Zenit*Magnitude; */

	    ibm->x_bp[i] = ibm->x_bp0[i] + c2*displacement.x * distance_to_axis*distance_to_Zenit*Magnitude;
	    
	    ibm->y_bp[i] = ibm->y_bp0[i] + c2*displacement.y * distance_to_axis*distance_to_Zenit*Magnitude;
	   
	    ibm->z_bp[i] = ibm->z_bp0[i] + c2*displacement.z * distance_to_axis*distance_to_Zenit*Magnitude;

	  }
	
	
      } 
/*   } else {  */
/*     /\*   LV beat for baloon LV *\/ */
/*     for (i=0; i<n_v; i++) { */
/*       zc=ibm->z_bp0[i]; */
/*       yc=ibm->y_bp0[i]; */
/*       xc=ibm->x_bp0[i]; */
      
/*       c2=7.*St_exp*time; */
/*       dz=(zc-z1)/(z0-z1); */
/*       rad = sqrt(xc * xc + yc * yc);  */
/*       az0=c1*c2*dz*(1.-dz)*rad; */
/*       //az0=c1*r; */
      
/*       if (rad>1.e-10){ */
/* 	rat[0]=xc/rad; */
/* 	rat[1]=yc/rad; */
/*       }else { */
/* 	rat[0]=0.; */
/* 	rat[1]=0.; */
/*       }     */
      
/*       ibm->x_bp[i]=ibm->x_bp0[i] - */
/* 	az0*sin(omega*time)*rat[0]; */
/*       //      az0*sin(kwave*dz-omega*time)*rat[0]; */
/*       ibm->y_bp[i]=ibm->y_bp0[i] - */
/* 	az0*sin(omega*time)*rat[1]; */
/*       //az0*sin(kwave*dz-omega*time)*rat[1]; */
/*     } */
/*   } */
  
  
  /*    Calculate the new normal & velcity */
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal  Vol=0.,p1[3],p2[3],p3[3],p4[3],p5[3],p6[3],sign,nf[3];
  PetscReal  edge1[3],edge2[3],edge3[3],edge2c3[3],volTH,cent[3],cent_o[3],dir[3];

  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;
    
    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;


    // calculate volume flux
    nf[0]=ibm->nf_x[i];nf[1]=ibm->nf_y[i];nf[2]=ibm->nf_z[i];

    p1[0]=ibm->x_bp[n1e];p1[1]=ibm->y_bp[n1e];p1[2]=ibm->z_bp[n1e];
    p2[0]=ibm->x_bp[n2e];p2[1]=ibm->y_bp[n2e];p2[2]=ibm->z_bp[n2e];
    p3[0]=ibm->x_bp[n3e];p3[1]=ibm->y_bp[n3e];p3[2]=ibm->z_bp[n3e];

    p4[0]=ibm->x_bp_o[n1e];p4[1]=ibm->y_bp_o[n1e];p4[2]=ibm->z_bp_o[n1e];
    p5[0]=ibm->x_bp_o[n2e];p5[1]=ibm->y_bp_o[n2e];p5[2]=ibm->z_bp_o[n2e];
    p6[0]=ibm->x_bp_o[n3e];p6[1]=ibm->y_bp_o[n3e];p6[2]=ibm->z_bp_o[n3e];

    for (ii=0; ii<3; ii++) {
      cent[ii]  =(p1[ii]+p2[ii]+p3[ii])/3.;
      cent_o[ii]=(p4[ii]+p5[ii]+p6[ii])/3.;
    }
    SUB(dir,cent,cent_o);
    sign=DOT(dir,nf);
    if (fabs(sign)>1e-15) 
      sign /=fabs(sign);
    else
      sign =0.;

    SUB(edge1,p4,p1);
    SUB(edge2,p4,p2);
    SUB(edge3,p4,p3);
    CROSS(edge2c3,edge2,edge3);
    volTH=DOT(edge1,edge2c3);
    
    Vol +=sign*fabs(volTH/6.)/delti;

    SUB(edge1,p5,p4);
    SUB(edge2,p5,p2);
    SUB(edge3,p5,p3);
    CROSS(edge2c3,edge2,edge3);
    volTH=DOT(edge1,edge2c3);
    
    Vol +=sign*fabs(volTH/6.)/delti;

    SUB(edge1,p6,p5);
    SUB(edge2,p6,p4);
    SUB(edge3,p6,p3);
    CROSS(edge2c3,edge2,edge3);
    volTH=DOT(edge1,edge2c3);
    
    Vol +=sign*fabs(volTH/6.)/delti;
  }
  
  PetscPrintf(PETSC_COMM_WORLD, "Volume LV Flux %e\n", Vol);

  FluxInSum = Vol;

  PetscReal  v_max=0.;
  PetscInt   i_vmax=0;
  for (i=0; i<n_v; i++) {
    if (ti>0) {
      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / delti;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / delti;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / delti;
    } else {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }

    if (v_max<fabs(ibm->u[i].z)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].z);
    }
    if (v_max<fabs(ibm->u[i].y)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].y);
    }
    if (v_max<fabs(ibm->u[i].x)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].x);
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "MAX LV Velocity: %d %le %le %le %le\n",ti, v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax]);
  
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  FILE *f;
  char filen[80];

  if (!rank) {
    if (ti == (ti/tiout)*tiout) {
      sprintf(filen, "surface%3.3d_%2.2d.dat",ti,ibi);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
      for (i=0; i<n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
      }
      for (i=0; i<n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
      }
      for (i=0; i<n_v; i++) {	
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);

    }
    if (!rank) {
      sprintf(filen, "InflowFlux.dat");
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le\n",ti, FluxInSum, delti);
      fclose(f);
    }

  }

  return(0);
}

PetscErrorCode LV_Volume_Flux(UserCtx *user) 
{
  DM		da = user->da;//, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
 

 PetscReal   Volume=0., Volume_o=0.;
 PetscReal   VolumeSum=0., Volume_oSum=0.;
 PetscReal   ***aj, ***nvert,***nvert_o;
 
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lNvert_o, &nvert_o);
  DMDAVecGetArray(da, user->lAj, &aj);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i]<0.1)
	  Volume += 1/aj[k][j][i];

	if (nvert_o[k][j][i]<0.1)
	  Volume_o += 1/aj[k][j][i];
      }
    }
  }

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);
  DMDAVecRestoreArray(da, user->lAj, &aj);

  MPI_Allreduce(&Volume, &VolumeSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Volume_o, &Volume_oSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

/*   PetscGlobalSum(&Volume, &VolumeSum, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&Volume_o, &Volume_oSum, PETSC_COMM_WORLD); */

  if (ti>0)
    FluxInSum = (Volume_oSum-VolumeSum)/user->dt;
  else
    FluxInSum = 0.;


  PetscPrintf(PETSC_COMM_WORLD, "Volume LV %e %e Flux %e %e\n", VolumeSum, Volume_oSum, FluxInSum, Flux_in);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
    FILE *f;
    char filen[80];
    //sprintf(filen, "Converge_dU%1.1d",dir);
    sprintf(filen, "InflowFlux.dat");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le\n",ti, FluxInSum,VolumeSum, Volume_oSum, user->dt);
    fclose(f);
  }

  FluxInSum = Flux_in;
  return 0;
}

PetscErrorCode Bed_Change(IBMNodes *ibm, PetscReal delti, PetscInt tstep) 
{
  PetscInt vert;
  PetscReal  PI = 3.141592653589793;
  for( vert = 0; vert < ibm->n_v; vert++ )
    {
      if(ibm->z_bp[ vert ]<1.e-5)
	{
	}
      else
	{
	  if(ibm->x_bp[vert]<= 20.5 && ibm->x_bp[vert]>=9.5)
	    { 
              ibm->z_bp[vert] = 1.+0.38*0.5*(1.-cos(2*PI*delti*tstep*0.1*0.037));
	    }
	  else
	    {
	    }
	}			
    }
  return(0);
}

/* ==================================================================================             */
Cmpnts rotate_y(Cmpnts v0, Cmpnts v_c, Cmpnts rot)
{
  Cmpnts v1;

  // only rot around y
  v1.x =  v_c.x + (v0.z-v_c.z)*sin(rot.y) + (v0.x-v_c.x)*cos(rot.y);
  v1.y =  v0.y;
  v1.z =  v_c.z + (v0.z-v_c.z)*cos(rot.y) - (v0.x-v_c.x)*sin(rot.y);
 
  return(v1);
}

Cmpnts rotateT_y(Cmpnts v0, Cmpnts v_c, Cmpnts rot)
{
  Cmpnts v1; rot.y=-rot.y;

  // only rot around y
  v1.x =  v_c.x + (v0.z-v_c.z)*sin(rot.y) + (v0.x-v_c.x)*cos(rot.y);
  v1.y =  v0.y;
  v1.z =  v_c.z + (v0.z-v_c.z)*cos(rot.y) - (v0.x-v_c.x)*sin(rot.y);
 
  return(v1);
}

Cmpnts rotate_xyz(Cmpnts v, Cmpnts v_c, Cmpnts rot)
{
  Cmpnts v0, v1, v2;

  // only rot around x
  v0.x  =    v.x;
  v0.y  = v_c.y - (v.z-v_c.z)*sin(rot.x) + (v.y-v_c.y)*cos(rot.x);
  v0.z  = v_c.z + (v.z-v_c.z)*cos(rot.x) + (v.y-v_c.y)*sin(rot.x);


  // only rot around y
  v1.x =  v_c.x + (v0.z-v_c.z)*sin(rot.y) + (v0.x-v_c.x)*cos(rot.y);
  v1.y =  v0.y;
  v1.z =  v_c.z + (v0.z-v_c.z)*cos(rot.y) - (v0.x-v_c.x)*sin(rot.y);

    // around z
  v2.x =  v_c.x - (v1.y-v_c.y)*sin(rot.z) + (v1.x-v_c.x)*cos(rot.z);
  v2.y =  v_c.y + (v1.y-v_c.y)*cos(rot.z) + (v1.x-v_c.x)*sin(rot.z);
  v2.z =     v1.z;
  
  return(v2);
}

Cmpnts rotateT_xyz(Cmpnts v, Cmpnts rot)
{
  Cmpnts v0, v1, v2;
  rot.z=-rot.z;
  rot.y=-rot.y;
  rot.x=-rot.x;

    // around z
  v0.x =  - (v.y)*sin(rot.z) + (v.x)*cos(rot.z);
  v0.y =  + (v.y)*cos(rot.z) + (v.x)*sin(rot.z);
  v0.z =     v.z;

  // only rot around y
  v1.x =  (v0.z)*sin(rot.y) + (v0.x)*cos(rot.y);
  v1.y =   v0.y;
  v1.z =  (v0.z)*cos(rot.y) - (v0.x)*sin(rot.y);

  // only rot around x
  v2.x  =    v1.x;
  v2.y  = - (v1.z)*sin(rot.x) + (v1.y)*cos(rot.x);
  v2.z  = + (v1.z)*cos(rot.x) + (v1.y)*sin(rot.x);
  
  return(v2);
}

/* ==================================================================================             */
PetscErrorCode rotate_force(FSInfo *fsi, PetscInt ibi, PetscInt bi)
{  
  Cmpnts    F, rot, a_c;
  PetscReal lift1, lift2, drag1, drag2, radial2, lift3, drag3;

  rot.x=fsi->S_ang_n[0];//-FSinfo->S_ang_o[0];
  rot.y=fsi->S_ang_n[2];//-FSinfo->S_ang_o[0];
  rot.z=fsi->S_ang_n[4];//-FSinfo->S_ang_o[0];

  F.x = fsi->F_x;
  F.y = fsi->F_y;
  F.z = fsi->F_z;

  lift1 =  F.z*cos(rot.y)-F.x*sin(rot.y);
  drag1 =  F.z*sin(rot.y)+F.x*cos(rot.y);

  lift3 =  F.z*cos(rot.y);
  drag3 =  F.z*sin(rot.y);

  a_c.x=a_c.y=a_c.z=0.;

  PetscPrintf(PETSC_COMM_WORLD, "Force Rotate: %le %le %le \n",rot.x,rot.y,rot.z);

  F = rotate_xyz(F, a_c, rot);
  //F = rotateT_xyz(F, rot);
  
  lift2 =  F.z;
  drag2 =  F.x*cos(-rot.z)-F.y*sin(-rot.z);
  radial2= F.x*sin(-rot.z)+F.y*cos(-rot.z);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);                                                                                                                                                                                                              
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Inertial%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le %le\n",ti, F.x,F.y,F.z, lift1, lift2,lift3, drag1, drag2, drag3, fsi->F_y, radial2);
    fclose(f);
  }

  return(0);
}

PetscErrorCode rotate_ibm(IBMNodes *ibm, FSInfo *fsi)
{
  PetscInt n_v = ibm->n_v;
  PetscInt i;
  Cmpnts alfa_c, a_c, rr;

  // axis of frame rotation
  a_c.x=fsi->x_c;
  a_c.y=fsi->y_c;
  a_c.z=fsi->z_c;
  
  // angle of rotation
  alfa_c.x=fsi->S_ang_n[0];
  alfa_c.y=fsi->S_ang_n[2];
  alfa_c.z=fsi->S_ang_n[4];

  //  PetscPrintf(PETSC_COMM_WORLD, "Body Rotate: %le %le %le Center %le %le %le\n",rot_x,rot_y,rot_z,x_c,y_c,z_c);

  for (i=0; i<n_v; i++) {
    rr.x = ibm->x_bp0[i];
    rr.y = ibm->y_bp0[i];
    rr.z = ibm->z_bp0[i];

    rr   = rotate_y(rr, a_c, alfa_c);

    ibm->x_bp[i] = rr.x;
    ibm->y_bp[i] = rr.y;
    ibm->z_bp[i] = rr.z;   
  }
  return(0);
}

PetscErrorCode wing_motion(IBMNodes *ibm, FSInfo *fsi, PetscInt ti, PetscReal *dt){

  /*===================  model constants  ============== */
  PetscReal  pi = 3.141592653589793, twopi=2*pi;
  PetscInt   N_period=3000;
  PetscOptionsGetInt(PETSC_NULL, "-N_period", &N_period, PETSC_NULL);
  
  PetscReal  T_period=5.2;//4.133674;  
  *dt=T_period/N_period; 
  //  PetscReal  pr_dtau_d=0.03 ,pr_dtau_t1=0.01, pr_dtau_r1=0.01;
  PetscReal  pr_dtau_d=0.24 ,pr_dtau_t1=0.12, pr_dtau_r1=0.24;
  PetscOptionsGetReal(PETSC_NULL, "-pr_dtau_r", &pr_dtau_r1, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-pr_dtau_t", &pr_dtau_t1, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-pr_dtau_d", &pr_dtau_d, PETSC_NULL);
  PetscReal  alfa_0=pi/4.;
  PetscReal  L0=0.19, Ltip=0.25 , rtil_0=Ltip/L0;
  PetscReal  Utip=0.26, U0=0.166, upl_0=Utip/U0;

  /*===================      variables     ============== */
  PetscReal  time_T, delt_t;
  PetscReal  tau_r2,tau_r3,tau_r4;
  PetscReal  tau_t2,tau_t3,tau_t4;
  PetscReal  al_0, alfa_1, dalfa_0, alfa_r, dalfa_r;
  PetscReal  upl_t, psi_t, dpsi_t;
  PetscReal  tau_r1, dtau_r1, dtau_r2;
  PetscReal  tau_t1, dtau_t1, dtau_t2, dtau_d;

  dtau_d =T_period*pr_dtau_d;
  dtau_t1=T_period*pr_dtau_t1;
  dtau_r1=T_period*pr_dtau_r1;

  if(dtau_t1==0.0)
    dtau_t2=T_period;
  else
    dtau_t2=T_period/2.0 - dtau_t1;
    
  if(dtau_r1==0.0) 
    dtau_r2=T_period;
  else
    dtau_r2=T_period/2.0 - dtau_r1;
  
  tau_t1=dtau_t2    ;
  tau_r1=tau_t1+0.5*dtau_t1-dtau_d; //tau_t1-0.5*(dtau_r1-dtau_t1)-dtau_d;
    
  /*===================   rotation  alfa  ============== */
  time_T=ti - (ti/N_period*N_period);
  time_T=time_T*( *dt);
  tau_r2=tau_r1+dtau_r1;
  tau_r3=tau_r2+dtau_r2;
  tau_r4=tau_r3+dtau_r1;
  
  alfa_1=pi - alfa_0 ;
  
  if(alfa_0<pi/2.) 
    al_0=alfa_0;
  else if(alfa_0>pi)  
    al_0=alfa_0-pi;
  
  dalfa_0=(pi - 2.*al_0)/(0.5*dtau_r1*(1.-sin(twopi)/twopi));
  
/*   PetscPrintf(PETSC_COMM_WORLD, "time, time_T, T_period %le %le %le\n", ti*( *dt), time_T, T_period); */
/*   PetscPrintf(PETSC_COMM_WORLD, "tau_r1-4 %le %le %le %le\n",tau_r1,tau_r2,tau_r3,tau_r4); */
  
  /*  !write(3,*) 'time,time_T,T_period',time,time_T,T_period */
  /*  !write(3,*) 'tau_r1,tau_r2,tau_r3,tau_r4',tau_r1,tau_r2,tau_r3,tau_r4 */
  
  if(tau_r1<=time_T && time_T<=tau_r2) {
    delt_t=(time_T-tau_r1);
    alfa_r=alfa_0+0.5*dalfa_0*(delt_t-(dtau_r1/twopi)*sin(twopi*delt_t/dtau_r1));
    dalfa_r=0.5*dalfa_0*(1-cos(twopi*delt_t/dtau_r1));  
    //    !write(3,*) 'tau_r1,time_T,tau_r2,alfa_r',tau_r1,time_T,tau_r2,alfa_r
  } else if(tau_r2<= time_T && time_T<=tau_r3) {
    alfa_r=alfa_1;
    dalfa_r=0.;
    //    !write(3,*) 'tau_r2,time_T,tau_r3,alfa_r',tau_r2,time_T,tau_r3,alfa_r
  } else if(tau_r3<= time_T && time_T<=tau_r4) {
    delt_t=(time_T-tau_r3);
    alfa_r=alfa_1-0.5*dalfa_0*(delt_t-(dtau_r1/twopi)*sin(twopi*delt_t/dtau_r1));
    dalfa_r=-0.5*dalfa_0*(1-cos(twopi*delt_t/dtau_r1));
    //    !write(3,*) 'tau_r3,time_T,tau_r4,alfa_r',tau_r3,time_T,tau_r4,alfa_r
  } else {
    alfa_r=alfa_0  ;
    dalfa_r= 0.;
  }
  
  if(alfa_0<0.0){
    alfa_r=0.0;
    dalfa_r=0.0;
  }
  //  !write(3,*) 'alfa_r',alfa_r
  
  /*===================   translation motion:  psi  ==============*/
  
  tau_t2=tau_t1+dtau_t1;
  tau_t3=tau_t2+dtau_t2;
  tau_t4=tau_t3+dtau_t1;

  // PetscPrintf(PETSC_COMM_WORLD, "tau_t1-4 %le %le %le %le\n",tau_t1,tau_t2,tau_t3,tau_t4);
  
  /*  !write(3,*)  */
  /*  !write(3,*) 'tau_t1-t4',tau_t1,tau_t2,tau_t3,tau_t4 */
  
  if(0.0<=time_T && time_T<=tau_t1) {
    delt_t=time_T;
    upl_t=upl_0;
    //    !write(3,*) '01: time_T,psi_t,u_pl',time_T,psi_t*180./pi,upl_t
  } else if(tau_t1<=time_T && time_T<=tau_t2) {
    delt_t=(time_T-tau_t1);
    upl_t=upl_0*cos(pi*delt_t/dtau_t1);
    //    !write(3,*) '12: time_T,psi_t,u_pl',time_T,psi_t*180./pi,upl_t
  } else if(tau_t2<= time_T && time_T<=tau_t3) {
    delt_t=(time_T-tau_t2);
    upl_t=-upl_0;
    //    !write(3,*) '23: time_T,psi_t,u_pl',time_T,psi_t*180./pi,upl_t
  } else if(tau_t3<= time_T && time_T<=tau_t4) {
    delt_t=(time_T-tau_t3);
    upl_t=-upl_0*cos(pi*delt_t/dtau_t1);
    //    !write(3,*) '34: time_T,psi_t,u_pl',time_T,psi_t*180./pi,upl_t
  } else if(tau_t4<= time_T && time_T<=T_period) {
    delt_t=(time_T-tau_t4);
    upl_t=upl_0;
    //    !write(3,*) '40: time_T,psi_t,u_pl',time_T,psi_t*180./pi,upl_t
  }
  psi_t=fsi->S_ang_n[4]-(upl_t/rtil_0)*( *dt);
  dpsi_t=-upl_t/rtil_0;
  
  //x-dir
  fsi->S_ang_n[0]=0.;
  fsi->S_ang_n[1]=0.;
  //y-dir
  fsi->S_ang_n[2]=-alfa_r;//pi*sin(2*pi*t);//0.;
  fsi->S_ang_n[3]=-dalfa_r;//pi*cos(2*pi*t);//0.;
  //z-dir
  fsi->S_ang_n[4]= psi_t;//pi*sin(2*pi*t);
  fsi->S_ang_n[5]= dpsi_t;//pi*cos(2*pi*t);

  PetscPrintf(PETSC_COMM_WORLD, "ti %d %le  psi %le %le alfa %le %le\n",ti,*dt,fsi->S_ang_n[4],fsi->S_ang_n[5],fsi->S_ang_n[2],fsi->S_ang_n[3]);
    
  return(0);
}

/* 	current = user->ibmlist[P1].head; */
/* 	D_min=1.e6; */
/* 	Dist=0.0; */
	
/* 	particle_found=0; */

/* 	while (current) { */
/* 	  ibminfo = &current->ibm_intp; */
/* 	  current = current->next; */
/* 	  i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk; */
/* 	  if (fabs(nvert[k][j][i]-1.0-ibi/1000.0)>1.e-5) { */
/* 	    P2=(int)((nvert[k][j][i]-1.0)*1001); */
/* 	    PetscPrintf(PETSC_COMM_WORLD, "Partciles %d amd %d are inside each other i %d j %d k %d \n",P1,P2,i,j,k); */
/* 	    goto nextp1; */
/* 	  } */
	  
/* 	  for (kk=0; kk<5; kk++){ */
/* 	    for (jj=0; jj<5; jj++){ */
/* 	      for (ii=0; ii<5; ii++){ */
		
/* 		if (P2==(int)((nvert[k+kk-2][j+jj-2][i+ii-2]-1.000)*1001) && fabs(nvert[k][j][i]-nvert[k+kk-2][j+jj-2][i+ii-2])>1.e-5) { */
	
/* 		  cell1=ibminfo->cell; */

/* 		  p1_x=user->ibm[P1].cent_x[cell1]; */
/* 		  p1_y=user->ibm[P1].cent_y[cell1]; */
/* 		  p1_z=user->ibm[P1].cent_z[cell1]; */

/* 		  for (elmt=0;elmt<ibm[P2].n_elmt;elmt++){  */
/* 		    p2_x=ibm[P2].cent_x[elmt]; */
/* 		    p2_y=ibm[P2].cent_y[elmt]; */
/* 		    p2_z=ibm[P2].cent_z[elmt]; */
/* 		    Dist=sqrt((p1_x-p2_x)*(p1_x-p2_x)+(p1_y-p2_y)*(p1_y-p2_y)+(p1_z-p2_z)*(p1_z-p2_z)); */
/* 		    if (Dist<D_min){ */
/* 		      D_min=Dist; */
/* 		      Point1.x=p1_x; */
/* 		      Point1.y=p1_y; */
/* 		      Point1.z=p1_z; */
/* 		      // cell_1=cell1; */
/* 		      //cell_2=cell2; */
/* 		      Point2.x=p2_x; */
/* 		      Point2.y=p2_y; */
/* 		      Point2.z=p2_z; */
/* 		    } */
/* 		  } */
		  		  
/* 		  particle_found++; */
				
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
	  
/* 	} */
/* 	in.value=D_min; */
/* 	in.index=rank; */
/* 	MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,0,PETSC_COMM_WORLD); */
	 
/* 	MPI_Allreduce(&D_min, &D_MIN,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD); */
/* 	MPI_Allreduce(&particle_found, &particle_found,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD); */
/* 	//	PetscPrintf(PETSC_COMM_WORLD, "partciles closest distance : %le cell-1 %d p1-x %le p1-y %le p1-z %le cell-2 %d p2-x %le p2-y %le p2-z %le \n",D_min,cell_1,Point1.x,Point1.y,Point1.z,cell_2,Point2.x,Point2.y,Point2.z); */
/* 	if (particle_found) PetscPrintf(PETSC_COMM_WORLD, "closest distance between paricle %d and partcile %d is  : %le \n",P1,P2,D_MIN); */

/* 	PetscReal dir_mag; */
/* 	PetscReal omega[3]; */

/* 	dir_min[0]=(Point2.x-Point1.x); */
/* 	dir_min[1]=(Point2.y-Point1.y); */
/* 	dir_min[2]=(Point2.z-Point1.z); */
/* 	dir_mag=sqrt((Point1.x-Point2.x)*(Point1.x-Point2.x)+(Point1.y-Point2.y)*(Point1.y-Point2.y)+(Point1.z-Point2.z)*(Point1.z-Point2.z)); */
/* 	dir_min[0] = dir_min[0]/dir_mag; */
/* 	dir_min[1] = dir_min[1]/dir_mag; */
/* 	dir_min[2] = dir_min[2]/dir_mag; */

/* 	R1[0]=(Point1.x-fsi[P1].x_c); */
/* 	R1[1]=(Point1.y-fsi[P1].y_c); */
/* 	R1[2]=(Point1.z-fsi[P1].z_c); */

/* 	omega[0]=fsi[P1].S_ang_n[1]; */
/* 	omega[1]=fsi[P1].S_ang_n[3]; */
/* 	omega[2]=fsi[P1].S_ang_n[5]; */

/* 	PetscPrintf(PETSC_COMM_WORLD, "omega1------w_x %le w_y %le w_z %le \n",omega[0],omega[1],omega[2]); */
/* 	CROSS(U_w1,omega,R1); */

/* 	dir_mag=sqrt(R1[0]*R1[0]+R1[1]*R1[1]+R1[2]*R1[2]); */
/* 	R1[0]=R1[0]/dir_mag; */
/* 	R1[1]=R1[1]/dir_mag; */
/* 	R1[2]=R1[2]/dir_mag; */

/* 	R2[0]=(Point2.x-fsi[P2].x_c); */
/* 	R2[1]=(Point2.y-fsi[P2].y_c); */
/* 	R2[2]=(Point2.z-fsi[P2].z_c); */

/* 	omega[0]=fsi[P2].S_ang_n[1]; */
/* 	omega[1]=fsi[P2].S_ang_n[3]; */
/* 	omega[2]=fsi[P2].S_ang_n[5]; */

/* 	PetscPrintf(PETSC_COMM_WORLD, "omega2------w_x %le w_y %le w_z %le \n",omega[0],omega[1],omega[2]); */
/* 	CROSS(U_w2,omega,R2); */

/* 	dir_mag=sqrt(R2[0]*R2[0]+R2[1]*R2[1]+R2[2]*R2[2]); */

/* 	R2[0]=R2[0]/dir_mag; */
/* 	R2[1]=R2[1]/dir_mag; */
/* 	R2[2]=R2[2]/dir_mag; */
     
/* 	int root=out.index; */
/* 	MPI_Bcast(&root,1,MPI_DOUBLE,0,PETSC_COMM_WORLD); */

/* 	MPI_Bcast(dir_min,3,MPI_DOUBLE,root,PETSC_COMM_WORLD); */
/* 	MPI_Bcast(R1,3,MPI_DOUBLE,root,PETSC_COMM_WORLD); */
/* 	MPI_Bcast(R2,3,MPI_DOUBLE,root,PETSC_COMM_WORLD); */
/* 	MPI_Bcast(U_w1,3,MPI_DOUBLE,root,PETSC_COMM_WORLD); */
/* 	MPI_Bcast(U_w2,3,MPI_DOUBLE,root,PETSC_COMM_WORLD); */

/* 	PetscPrintf(PETSC_COMM_WORLD, "min_dir------n_x %le n_y %le n_z %le \n",dir_min[0],dir_min[1],dir_min[2]); */

/* 	PetscPrintf(PETSC_COMM_WORLD, "R1------n_x %le n_y %le n_z %le \n",R1[0],R1[1],R1[2]); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "R2------n_x %le n_y %le n_z %le \n",R2[0],R2[1],R2[2]); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "x1_c %le y1_c %le z1_c %le \n",fsi[P1].x_c,fsi[P1].y_c,fsi[P1].z_c); */
/* 	PetscPrintf(PETSC_COMM_WORLD, "x2_c %le y2_c %le z2_c %le \n",fsi[P2].x_c,fsi[P2].y_c,fsi[P2].z_c); */

/*       } */
/*     } */
/*   } */
/* 	while (current) { */
/* 	  ibminfo = &current->ibm_intp; */
/* 	  current = current->next; */
/* 	  i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk; */
/* 	  if (fabs(nvert[k][j][i]-1.0-ibi/1000.0)>1.e-5) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD, "Partciles inside each other i %d j %d k %d \n",i,j,k); */
/* 	    /\*       goto nextp1; *\/ */
/* 	  } */
/* 	  cell1=ibminfo->cell; */
/* 	  p1_x=user->ibm[P1].cent_x[cell1]; */
/* 	  p1_y=user->ibm[P1].cent_y[cell1]; */
/* 	  p1_z=user->ibm[P1].cent_z[cell1]; */
/* 	  //PetscPrintf(PETSC_COMM_SELF, "first partcile cell %d ---x :%le y :%le z :%le  \n",cell1,p1_x,p1_y,p1_z); */
/* 	  current2 = user->ibmlist[P2].head; */
/* 	  while (current2) { */
/* 	    ibminfo2=&current2->ibm_intp; */
/* 	    current2=current2->next; */
/* 	    cell2=ibminfo2->cell; */
/* 	    p2_x=user->ibm[P2].cent_x[cell2]; */
/* 	    p2_y=user->ibm[P2].cent_y[cell2]; */
/* 	    p2_z=user->ibm[P2].cent_z[cell2]; */
	    
/* 	    //PetscPrintf(PETSC_COMM_SELF, "second partcile cell %d ---x :%le y :%le z :%le  \n",cell2,p2_x,p2_y,p2_z); */
/* 	    Dist=sqrt((p1_x-p2_x)*(p1_x-p2_x)+(p1_y-p2_y)*(p1_y-p2_y)+(p1_z-p2_z)*(p1_z-p2_z)); */
/* 	    if (Dist<D_min){ */
/* 	      D_min=Dist; */
/* 	      Point1.x=p1_x; */
/* 	      Point1.y=p1_y; */
/* 	      Point1.z=p1_z; */
/* 	      cell_1=cell1; */
/* 	      cell_2=cell2; */
/* 	      Point2.x=p2_x; */
/* 	      Point2.y=p2_y; */
/* 	      Point2.z=p2_z; */
/* 	    } */
	    
/* 	  } */
/* 	} */
/* 	for (i=0;i<ibm[P1].n_elmt;i++){ */
/* 	  for (j=0;j<ibm[P2].n_elmt;j++){ */
/* 	    p1_x=ibm[P1].cent_x[i]; */
/* 	    p1_y=ibm[P1].cent_y[i]; */
/* 	    p1_z=ibm[P1].cent_z[i]; */
/* 	    p2_x=ibm[P2].cent_x[j]; */
/* 	    p2_y=ibm[P2].cent_y[j]; */
/* 	    p2_z=ibm[P2].cent_z[j]; */
/* 	    Dist=sqrt((p1_x-p2_x)*(p1_x-p2_x)+(p1_y-p2_y)*(p1_y-p2_y)+(p1_z-p2_z)*(p1_z-p2_z)); */
/* 	    if (Dist<D_min){ */
/* 	      D_min=Dist; */
/* 	      Point1.x=p1_x; */
/* 	      Point1.y=p1_y; */
/* 	      Point1.z=p1_z; */
/* 	      cell_1=i; */
/* 	      cell_2=j; */
/* 	      Point2.x=p2_x; */
/* 	      Point2.y=p2_y; */
/* 	      Point2.z=p2_z; */
/* 	    } */
/* 	  } */
/* 	} */
/* 	PetscPrintf(PETSC_COMM_WORLD, "partciles closest distance : %le cell-1 %d p1-x %le p1-y %le p1-z %le cell-2 %d p2-x %le p2-y %le p2-z %le \n",D_min,cell_1,Point1.x,Point1.y,Point1.z,cell_2,Point2.x,Point2.y,Point2.z); */


/* ==================================================================================             */
/* Left Ventricle(second) Subroutines! */
/* ==================================================================================             */

PetscErrorCode LV_beat2(IBMNodes *ibm, PetscInt time, PetscReal delti,PetscInt ibi) {
  
  PetscInt   i, ii, n_v=ibm->n_v, n_elmt=ibm->n_elmt;
  PetscReal D0[80]={4.4416,4.4416,4.4416,4.4416,4.4164,4.3674,4.2976,4.2109,4.1116,4.0047,3.8951,3.7877,3.6881,3.6104,3.5501,3.5035,3.4676,3.4398,3.4166,3.4126,3.4126,3.4126,3.4126,3.4126,3.4126,3.4126,3.4126,3.4319,3.4614,3.5020,3.5514,3.6070,3.6667,3.7285,3.7908,3.8521,3.9113,3.9673,4.0193,4.0664,4.1082,4.1440,4.1737,4.1971,4.2143,4.2255,4.2313,4.2323,4.2334,4.2386,4.2501,4.2686,4.2933,4.3221,4.3524,4.3812,4.4059,4.4247,4.4363,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416,4.4416};
  
  //for (i=0; i<80; i++) D0[i]=D0[i]/1.628;
  
  PetscInt  c, c1,c2,c3;
  PetscReal position1[n_v], position11[n_v], position12[n_v];
  PetscReal position2[n_v], position21[n_v], position22[n_v];
  PetscReal position3[n_v], position31[n_v], position32[n_v];
  PetscReal pointer[n_v];
  PetscReal zm[n_v];
  PetscReal zt0[n_v];
  PetscReal zm0[n_v];
  PetscReal xt0[n_v];
  PetscReal yt0[n_v];
  PetscReal rm[n_v];
  PetscReal rm0[n_v];
  PetscReal H_scale=3.35, HoverD=2., H_tot=D0[0]*HoverD, H_cyl=D0[0]*0.5, z_scale=1.;
  PetscOptionsGetReal(PETSC_NULL, "-z_scale", &z_scale, PETSC_NULL);
  
  //H_scale = 10.8; //HoverD=2*3.22; H_cyl = H_cyl*3.22;
  //HoverD=1.;  

  PetscInt j=0;
  PetscInt period;
  PetscOptionsGetInt(PETSC_NULL, "-period", &period, PETSC_NULL);
  
  PetscInt n;
  c=1;
  n=80;
  PetscReal co1[80], co2[80], co3[80];  
  PetscInt  u[80];
  PetscReal  s;

  for (i=0;i<80;i++) co1[i]=0; 
  for (i=0;i<80;i++) co2[i]=0; 
  for (i=0;i<80;i++) co3[i]=0; 
  for (i=0;i<80;i++) u[i]=i; 
  spli(n,u,D0,co1,co2,co3);
 
  s=(time-(time/period)*period)*80.00/period;

  PetscReal D;
  sev(n,s,u,D0,co1,co2,co3,&D);
  if ((time-time/period*period)>=3876) D=D0[0];
  /*  position1x is the initial coor of the LV */
  for(i=0;i<n_v;i++) {
    position1[i]=ibm->x_bp0[i];
    position11[i]=ibm->y_bp0[i];
    position12[i]=ibm->z_bp0[i];
  }

  j=0;
  /* get ride of the useless points and keep them in new vectors */
  for(i=0;i<n_v;i++) {
    if(position12[i]>0.5){
      pointer[j]=i;
      position3[j]=position1[i];
      position1[i]=0;
      position31[j]=position11[i];
      position11[i]=0;
      position32[j]=position12[i];
      position12[i]=0;
      j=j+1;
    }
  }
  
  /* move the original point to 0,0,0*/
  for(i=0;i<n_v;i++)   {
    position1[i]=position1[i]+0.3;
    position11[i]=position11[i]+0.3;
    position12[i]=position12[i]-0.4;
  }
  
  /* rotate coordinate and rotated coor is postion2x*/
  for(i=0;i<n_v;i++) {
    position2[i]=0.507*position1[i]-0.87*position11[i]-0*position12[i];
    position21[i]=0.62*position1[i]+0.36*position11[i]+0.75*position12[i];
    position22[i]=0.627*position1[i]+0.36*position11[i]-0.68*position12[i];
  }
  
  /* new zt*/
  for(i=0;i<n_v;i++)  {
    zt0[i]=position22[i];
    xt0[i]=position2[i];
    yt0[i]=position21[i];
    // zt=z0*D(c)/D0
    position22[i]=zt0[i]*D/D0[0]; 
  }
  
  /* z mapped in model,3.35 is the length of long axis of real heart*/
  /*   zm= zt* H/H_scale ==> zm=zt* H/D*D/H_scale */
  for(i=0;i<n_v;i++)  {
    zm[i]=position22[i]*HoverD*D0[0]/H_scale;
    zm0[i]=zt0[i]*HoverD*D0[0]/H_scale;
  }
  
  /* calculate new xt,yt*/
  /*   if z<H_cyl then cylinder part of the model, i.e., r=D/2 
       if z>H_cyl then round part of the model*/
  for(i=0;i<n_v;i++){
    if (zm[i]<H_cyl*D/D0[0]) {
      rm[i]=0.5*D;
    } else {
      rm[i]=D*(HoverD*D-zm[i])/(4.*(HoverD-H_cyl/D0[0]));
      if (rm[i]>=0)
	{
	  rm[i]=sqrt(rm[i]);
	}     else       {
	rm[i]=sqrt(-rm[i]);
      }
    }
    
    if (zm0[i]<H_cyl){
      rm0[i]=0.5*D0[0];
    } else  {
      rm0[i]=D0[0]*(HoverD*D0[0]-zm0[i])/(4.*(HoverD-H_cyl/D0[0]));
      if (rm0[i]>=0) {
	rm0[i]=sqrt(rm0[i]);
      }     else	{
	rm0[i]=sqrt(-rm0[i]);
      }
    }              
    // z_nomove=1, 1<zm<2 contrained, 3.35 is lenght of the long axis of heart
    PetscReal z_nomove=1.0/(2*4.4416)*HoverD*D0[0]; 
    if (zm[i]<=z_nomove)
      { // don't move
	position2[i]=xt0[i];
	position21[i]=yt0[i];
	position22[i]=zt0[i];
      }
    else if (zm[i]<2*z_nomove && zm[i]>z_nomove )
      { // move relative
	rm[i]=(2*z_nomove-zm[i])/z_nomove*rm0[i]+(zm[i]-z_nomove)/z_nomove*rm[i]*sqrt(z_scale);
	position22[i]=(2*z_nomove-zm[i])/z_nomove*zt0[i]+(zm[i]-z_nomove)/z_nomove*position22[i]/z_scale;
      }
    else
      {
	rm[i]=rm[i]*sqrt(z_scale);
	position22[i]=position22[i]/z_scale;
      }
    
    /*     move only points with z>1 */
    if (zm[i]>z_nomove)
      { // the locations can be scaled differnetly in differnt directions
	position2[i]=xt0[i]*rm[i]/rm0[i];
	position21[i]=yt0[i]*rm[i]/rm0[i];
	position22[i]=position22[i];
      }
  }
  
  j=0;
  for(i=0;i<n_v;i++) {                        /* rotate back*/
    
    if (i==pointer[j])
      {
	ibm->x_bp[i]=position3[j];
	ibm->y_bp[i]=position31[j];
	ibm->z_bp[i]=position32[j];
	j=j+1;
      }
    else
      {
	ibm->x_bp[i]=0.4965*position2[i]+0.5705*position21[i]+0.6293*position22[i]-0.3;
	ibm->y_bp[i]=-0.8601*position2[i]+0.3325*position21[i]+0.3667*position22[i]-0.3;
	ibm->z_bp[i]=0.0024*position2[i]+0.7021*position21[i]-0.6962*position22[i]+0.4;
      }
  }
  
  /*    Calcula1te the new normal & velcity */
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal  Vol=0.,p1[3],p2[3],p3[3],p4[3],p5[3],p6[3],sign,nf[3];
  PetscReal  edge1[3],edge2[3],edge3[3],edge2c3[3],volTH,cent[3],cent_o[3],dir[3];


  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;
    
    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
    
    
    // calculate volume flux
    if (time>0) {
      nf[0]=ibm->nf_x[i];nf[1]=ibm->nf_y[i];nf[2]=ibm->nf_z[i];
      
      p1[0]=ibm->x_bp[n1e];p1[1]=ibm->y_bp[n1e];p1[2]=ibm->z_bp[n1e];
      p2[0]=ibm->x_bp[n2e];p2[1]=ibm->y_bp[n2e];p2[2]=ibm->z_bp[n2e];
      p3[0]=ibm->x_bp[n3e];p3[1]=ibm->y_bp[n3e];p3[2]=ibm->z_bp[n3e];
      
      p4[0]=ibm->x_bp_o[n1e];p4[1]=ibm->y_bp_o[n1e];p4[2]=ibm->z_bp_o[n1e];
      p5[0]=ibm->x_bp_o[n2e];p5[1]=ibm->y_bp_o[n2e];p5[2]=ibm->z_bp_o[n2e];
      p6[0]=ibm->x_bp_o[n3e];p6[1]=ibm->y_bp_o[n3e];p6[2]=ibm->z_bp_o[n3e];
      
      for (ii=0; ii<3; ii++) {
	cent[ii]  =(p1[ii]+p2[ii]+p3[ii])/3.;
	cent_o[ii]=(p4[ii]+p5[ii]+p6[ii])/3.;
      }
      SUB(dir,cent,cent_o);
      sign=DOT(dir,nf);
      if (fabs(sign)>1e-15) 
	sign /=fabs(sign);
      else
	sign =0.;
      
      SUB(edge1,p4,p1);
      SUB(edge2,p4,p2);
      SUB(edge3,p4,p3);
      CROSS(edge2c3,edge2,edge3);
      volTH=DOT(edge1,edge2c3);
    
      Vol +=sign*fabs(volTH/6.)/delti;
      
      SUB(edge1,p5,p4);
      SUB(edge2,p5,p2);
      SUB(edge3,p5,p3);
      CROSS(edge2c3,edge2,edge3);
      volTH=DOT(edge1,edge2c3);
      
      Vol +=sign*fabs(volTH/6.)/delti;
      
      SUB(edge1,p6,p5);
      SUB(edge2,p6,p4);
      SUB(edge3,p6,p3);
      CROSS(edge2c3,edge2,edge3);
      volTH=DOT(edge1,edge2c3);
      
      Vol +=sign*fabs(volTH/6.)/delti;
      
    } else { //ti>0
      Vol = 0.;
    }  
  }
  
  PetscPrintf(PETSC_COMM_WORLD, "Volume LV Flux %e\n", Vol);
  
  FluxInSum = Vol;
  
  PetscReal  v_max=0.;
  PetscInt   i_vmax=0;
  for (i=0; i<n_v; i++) {
    if (time>0) {
      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / delti;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / delti;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / delti;
    } else {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
    
    if (v_max<fabs(ibm->u[i].z)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].z);
    }
    if (v_max<fabs(ibm->u[i].y)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].y);
    }
    if (v_max<fabs(ibm->u[i].x)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].x);
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "MAX LV Velocity: %d %le %le %le %le\n",time, v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax]);
  
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  FILE *f;
  char filen[80];
  
  if (!rank) {
    sprintf(filen, "InflowFlux.dat");
    f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le\n",time, FluxInSum, delti);
      fclose(f);
  }
  
  return(0);
}


PetscErrorCode spli (PetscInt n,PetscInt x[],PetscReal  y[],PetscReal b[],PetscReal c[],PetscReal d[])
  
{  /* begin procedure spline() */
  
  PetscInt    nm1, ib, i;
  PetscReal   t;
  nm1    = n - 1;
  
  /* ---- Set up the symmetric tri-diagonal system
     b = diagonal
     d = offdiagonal
     c = right-hand-side  */
  d[0] = x[1] - x[0];
  c[1] = (y[1] - y[0]) / d[0];
  for (i = 1; i < nm1; i++)
    {
      d[i]   = x[i+1] - x[i];
      b[i]   = 2.0 * (d[i-1] + d[i]);
      c[i+1] = (y[i+1] - y[i]) / d[i];
      c[i]   = c[i+1] - c[i];
    }
  
  /* ---- Default End conditions
     Third derivatives at x[0] and x[n-1] obtained
     from divided differences  */
  b[0]   = -d[0];
  b[nm1] = -d[n-2];
  c[0]   = 0.0;
  c[nm1] = 0.0;
  if (n != 3)
    {
      c[0]   = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
      c[nm1] = c[n-2] / (x[nm1] - x[n-3]) - c[n-3] / (x[n-2] - x[n-4]);
      c[0]   = c[0] * d[0] * d[0] / (x[3] - x[0]);
      c[nm1] = -c[nm1] * d[n-2] * d[n-2] / (x[nm1] - x[n-4]);
    }
  
  /* Forward elimination */
  for (i = 1; i < n; ++i)
    {
      t    = d[i-1] / b[i-1];
      b[i] = b[i] - t * d[i-1];
      c[i] = c[i] - t * c[i-1];
    }
  
  /* Back substitution */
  c[nm1] = c[nm1] / b[nm1];
  for (ib = 0; ib < nm1; ++ib)
    {
      i    = n - ib - 2;
      c[i] = (c[i] - d[i] * c[i+1]) / b[i];
    }
  
  /* c[i] is now the sigma[i] of the text */
  
  /* Compute the polynomial coefficients */
  b[nm1] = (y[nm1] - y[n-2]) / d[n-2] + d[n-2] * (c[n-2] + 2.0 * c[nm1]);
  for (i = 0; i < nm1; ++i)
    {
      b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2.0 * c[i]);
      d[i] = (c[i+1] - c[i]) / d[i];
      c[i] = 3.0 * c[i];
    }
  c[nm1] = 3.0 * c[nm1];
  d[nm1] = d[n-2];
  
}  /* end of spline() */
					     



PetscErrorCode sev (PetscInt n, PetscReal u,  PetscInt  x[],  PetscReal  y[], PetscReal b[], 
		    PetscReal c[], PetscReal d[], PetscReal *w)
  
/* int    n; */
/* double u; */
/* double x[], y[], b[], c[], d[]; */
/* int    *last; */
  
/* #endif */
  
/* PetscReal seval(PetscInt n, PetscReal u, PetscReal *x, PetscReal *y, */
/* 		PetscReal *b, PetscReal *c, PetscReal *d, PetscInt *last) */
  
/*Purpose ...
  -------
  Evaluate the cubic spline function
  
  S(xx) = y[i] + b[i] * w + c[i] * w**2 + d[i] * w**3
  where w = u - x[i]
  and   x[i] <= u <= x[i+1]
  Note that Horner's rule is used.
  If u < x[0]   then i = 0 is used.
  If u > x[n-1] then i = n-1 is used.

  Input :
  -------
  n       : The number of data points or knots (n >= 2)
  u       : the abscissa at which the spline is to be evaluated
  Last    : the segment that was last used to evaluate U
  x[]     : the abscissas of the knots in strictly increasing order
  y[]     : the ordinates of the knots
  b, c, d : arrays of spline coefficients computed by spline().

  Output :
  --------
  seval   : the value of the spline function at u
  Last    : the segment in which u lies

  Notes ...
  -----
  (1) If u is not in the same interval as the previous call then a
      binary search is performed to determine the proper interval.

*/
/*-------------------------------------------------------------------*/

{  /* begin function seval() */
  
  PetscInt   i, j, k;
   
  for (i=0;i<n;i++)
    {
      if ((i <=u) && (i >=u-1))
	{ 
	  *w = u - x[i];
	  *w = y[i] + *w * (b[i] + *w * (c[i] + *w * d[i]));
	  break;
	}
    }
  
  return(0);
}
