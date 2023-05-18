#include "variables.h"
#include "petsc.h"

extern PetscReal CMx_c,CMy_c,CMz_c, L_dim;
extern PetscInt  NumberOfBodies, tiout;
extern PetscInt  cop, ti;


PetscReal seval(PetscInt n, PetscReal u, PetscReal *x, PetscReal *y,
		PetscReal *b, PetscReal *c, PetscReal *d, PetscInt *last);

PetscReal Sval(PetscInt n, PetscReal u, PetscReal *x, PetscReal *y,
	       PetscReal *b, PetscReal *c, PetscReal *d, PetscInt *last);

PetscReal deriv(PetscInt n, PetscReal u, PetscReal *x, 
		PetscReal *b, PetscReal *c, PetscReal *d, PetscInt *last);

PetscErrorCode spline(PetscInt n, PetscInt end1, PetscInt end2, PetscReal slope1,
		      PetscReal slope2, PetscReal *x, PetscReal *y, PetscReal *b,
		      PetscReal *c, PetscReal *d, PetscInt *iflag);

PetscErrorCode read_midpoint(Cstart *cstart)
{
  PetscInt   rank, i;
  PetscInt   n_tstep, n_midp,n_subit, ts;
  PetscReal  *x_midp, *y_midp;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("midpoints.dat", "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open C-start node file!!! (midpoints.dat)");
    fscanf(fd, "%d %d %d\n", &n_tstep, &n_midp, &n_subit);

    //  PetscReal  x_midp[n_tstep][n_midp], y_midp[n_tstep][n_midp];
    
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &x_midp);
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &y_midp);
    PetscPrintf(PETSC_COMM_WORLD, "malloc %d %d %d! \n",n_tstep, n_midp, n_subit);

    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_midp, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_subit, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Bcast! \n");
    for (ts=0; ts<n_tstep; ts++) {
      for (i=0; i<n_midp; i++) {
	fscanf(fd, "%le %le\n", &x_midp[ts*n_midp+i], &y_midp[ts*n_midp+i]);      
      }
    }
    fclose(fd);
    PetscPrintf(PETSC_COMM_WORLD, "READ! \n");
    MPI_Bcast(x_midp, n_tstep*n_midp, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_midp, n_tstep*n_midp, MPIU_REAL, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Bcast! \n");

    cstart->n_time=n_tstep;
    cstart->n_subit=n_subit;
    cstart->n_midp=n_midp;
    cstart->x_midp=x_midp;
    cstart->y_midp=y_midp;
    PetscMalloc(n_midp*sizeof(PetscReal), &cstart->s1);
    PetscMalloc(n_midp*sizeof(PetscReal), &cstart->s2);
    PetscMalloc(n_midp*sizeof(PetscReal), &cstart->s3);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->st1);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->st2);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->st3);

    fd = fopen("com_motion.dat", "r"); 
    if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open C-start node file!!! (com_motion.dat)");

    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->x_com);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->y_com);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->head_ang);
    
    for (ts=0; ts<n_tstep; ts++) {
      fscanf(fd, "%le %le %le\n", &cstart->x_com[ts], &cstart->y_com[ts],
	     &cstart->head_ang[ts]);      
    }
    MPI_Bcast(cstart->x_com, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->y_com, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->head_ang, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
  } else {
    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_midp, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_subit, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_SELF, "Bcast Rank 1! \n");
    //PetscReal  x_midp[n_tstep][n_midp], y_midp[n_tstep][n_midp];
    // PetscMalloc(n_midp*sizeof(PetscReal), &x_midp);
    //PetscMalloc(n_midp*sizeof(PetscReal), &y_midp);
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &x_midp);
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &y_midp);

    MPI_Bcast(x_midp,n_tstep*n_midp, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_midp,n_tstep*n_midp, MPIU_REAL, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_SELF, "malloc Rank1 %d %d %d! \n",n_tstep, n_midp, n_subit);
    //    PetscPrintf(PETSC_COMM_SELF, "Bcast Rank1! \n");

    cstart->n_time=n_tstep;
    cstart->n_subit=n_subit;
    cstart->n_midp=n_midp;
    cstart->x_midp=x_midp;
    cstart->y_midp=y_midp;
    PetscMalloc(n_midp*sizeof(PetscReal), &cstart->s1);
    PetscMalloc(n_midp*sizeof(PetscReal), &cstart->s2);
    PetscMalloc(n_midp*sizeof(PetscReal), &cstart->s3);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->st1);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->st2);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->st3);

    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->x_com);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->y_com);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->head_ang);

    MPI_Bcast(cstart->x_com, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->y_com, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->head_ang, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  return(0); 
}


PetscErrorCode form_cstart(Cstart *cstart,IBMNodes *ibm, FSInfo *fsi,
			   PetscReal dt, PetscInt ibi, PetscInt nonIntertial)
{  
  PetscInt   flg, i,seg,last,ts,subts,tss;
  PetscReal  u[1001],v[1001],w[1001],*x_midp, *y_midp,*ubar,*dx,*dy;
  PetscReal  *x_midt, *y_midt,*ut,cl;
  PetscInt   rank;
  PetscBool print_spline=PETSC_FALSE;


  PetscMalloc(ibm->n_v*sizeof(PetscReal), &dx);
  PetscMalloc(ibm->n_v*sizeof(PetscReal), &dy);

  PetscMalloc(cstart->n_midp*sizeof(PetscReal), &ubar);
  PetscMalloc(cstart->n_midp*sizeof(PetscReal), &x_midp);
  PetscMalloc(cstart->n_midp*sizeof(PetscReal), &y_midp);
  PetscMalloc(cstart->n_time*sizeof(PetscReal), &ut);
  PetscMalloc(cstart->n_time*sizeof(PetscReal), &x_midt); 
  PetscMalloc(cstart->n_time*sizeof(PetscReal), &y_midt);
  
  for (tss=0; tss<cstart->n_time; tss++) {
    ut[tss]=tss*1.;
  }

  ts=ti/(1+cstart->n_subit);
  subts=ti-ts*(1+cstart->n_subit);
  PetscPrintf(PETSC_COMM_WORLD, "ti ts subts %d %d %d\n",ti,ts, subts);
  
  /* spline interpolation in time for finding midline points at
     the time between experimental time steps */
  for (seg=0; seg<cstart->n_midp; seg++) {
    ubar[seg]=seg*1.;

    if (subts) {
      for (tss=0; tss<cstart->n_time; tss++) {
	x_midt[tss]=cstart->x_midp[tss*cstart->n_midp+seg];
	y_midt[tss]=cstart->y_midp[tss*cstart->n_midp+seg];
	//      PetscPrintf(PETSC_COMM_WORLD, "Sval %d %le %le %le\n",seg,ut[tss],x_midt[tss],y_midt[tss]);
      }
      
      u[0]=ts+(subts*1./(double)(cstart->n_subit+1.));
      spline(cstart->n_time, 0, 0, 0., 0., ut, x_midt,
	     cstart->st1, cstart->st2, cstart->st3, &flg);
      //PetscPrintf(PETSC_COMM_WORLD, "Spline time intp done! %le %le %le %d\n",cstart->st1[3],cstart->st2[3],cstart->st3[3], flg);
      
      x_midp[seg]=Sval(cstart->n_midp, u[0]
			 , ut, x_midt, cstart->st1, cstart->st2, cstart->st3, &ts);
      
      spline(cstart->n_time, 0, 0, 0., 0., ut, y_midt,
	     cstart->st1, cstart->st2, cstart->st3, &flg);
      y_midp[seg]=Sval(cstart->n_midp,  u[0]
			 , ut, y_midt, cstart->st1, cstart->st2, cstart->st3, &ts);
      
      // PetscPrintf(PETSC_COMM_WORLD, "Sval %d %le %le %le\n",ts,u[0],x_midp[seg],y_midp[seg]);
    } else {
      //      x_midp[seg]=cstart->x_midp[seg];
      //      y_midp[seg]=cstart->y_midp[seg];
      x_midp[seg]=cstart->x_midp[ts*cstart->n_midp+seg];
      y_midp[seg]=cstart->y_midp[ts*cstart->n_midp+seg];
      //      PetscPrintf(PETSC_COMM_WORLD, "Sval %d %le %le %le \n",ts,ubar[seg],x_midp[seg],y_midp[seg]);
     }
  }
 

  /* spline for spatial interpolation on midpoints and move
     the mesh points similar to midline */

  // 1) in the y direction S(u)=Y(u) 
  spline(cstart->n_midp, 0, 0, 0., 0., ubar, y_midp, 
	 cstart->s1, cstart->s2, cstart->s3, &flg); 
  PetscPrintf(PETSC_COMM_WORLD, "Spline done! %le %le %le %d\n",cstart->s1[3],cstart->s2[3],cstart->s3[3], flg);

  if (print_spline) {
  for (i=0;i<1001; i++){
    u[i]=i*cstart->n_midp*0.001;
    v[i]=seval(cstart->n_midp, u[i], ubar, y_midp,
	       cstart->s1, cstart->s2, cstart->s3, &last);
    //PetscPrintf(PETSC_COMM_WORLD, "Sval %le %le \n",u[i],v[i]);
  }
  }
  
  for (i=0; i<ibm->n_v; i++) {
    u[0]=(ibm->z_bp0[i]-cstart->zmin)
      /  (  cstart->zmax    -cstart->zmin)*(cstart->n_midp-1.);
    seg=(int) (u[0]);
    ibm->z_bp[i]= //-cstart->y_midp[0]+
      Sval(cstart->n_midp, u[0]
	   , ubar, y_midp, cstart->s1, cstart->s2, cstart->s3, &seg); 
    dy[i]=deriv(cstart->n_midp, u[0] // this is tanget to the spline
		, ubar, cstart->s1, cstart->s2, cstart->s3, &seg); 
    /*       if (i==1936) */
    /* 	PetscPrintf(PETSC_COMM_WORLD, "Sval %d %d %le %le\n",i,seg,u[1],deriv(cstart->n_midp, u[0] */
    /* 									      , ubar, cstart->s1, cstart->s2, cstart->s3, &seg)); */
  }
  
  // 2) in the x direction S(u)=X(u) 
  spline(cstart->n_midp, 0, 0, 0., 0., ubar, x_midp, 
	 cstart->s1, cstart->s2, cstart->s3, &flg); 
  PetscPrintf(PETSC_COMM_WORLD, "Spline done! %le %le %le %d\n",cstart->s1[3],cstart->s2[3],cstart->s3[3], flg);

  if (print_spline) {
  for (i=0;i<1001; i++){
    u[i]=i*cstart->n_midp*0.001;
    w[i]=seval(cstart->n_midp, u[i], ubar, x_midp,
	       cstart->s1, cstart->s2, cstart->s3, &last);
    //PetscPrintf(PETSC_COMM_WORLD, "Sval %le %le \n",u[i],v[i]);
  }
  }

  for (i=0; i<ibm->n_v; i++) {
    u[0]=(ibm->z_bp0[i]-cstart->zmin)
      /  (  cstart->zmax    -cstart->zmin)*(cstart->n_midp-1.);
    seg=(int) (u[0]);
    ibm->x_bp[i]=ibm->x_bp0[i]+
      cstart->x_midp[0]-0.5*(cstart->xmin+cstart->xmax)+      
      Sval(cstart->n_midp, u[0]
	   , ubar, x_midp, cstart->s1, cstart->s2, cstart->s3, &seg)-cstart->x_midp[0]; 
    dx[i]=deriv(cstart->n_midp, u[0] // this is tangent to the spline 
		, ubar, cstart->s1, cstart->s2, cstart->s3, &seg); 
    /*  if (i==1936) */
    /* 	PetscPrintf(PETSC_COMM_WORLD, "Sval %d %d %le %le\n",i,seg,u[1],deriv(cstart->n_midp, u[0] */
    /* 	     , ubar, cstart->s1, cstart->s2, cstart->s3, &seg)); */
  }
  

  /* spline interpolation in time for finding the center of mass & angle at
     the time between experimental time steps */
  //PetscReal x_com, y_com, dxdt_com, dydt_com, rot_x, drot_x;

  if (nonIntertial) {
  u[0]=ts+(subts*1./(double)(cstart->n_subit+1.));
/*   spline(cstart->n_time, 0, 0, 0., 0., ut, cstart->x_com, */
/* 	 cstart->st1, cstart->st2, cstart->st3, &flg); */
/*   fsi->S_new[0]=Sval(cstart->n_time, u[0] */
/* 		     , ut, cstart->x_com, cstart->st1, cstart->st2, cstart->st3, &ts); */
/*   fsi->S_new[1]=deriv(cstart->n_time, u[0] */
/* 		      , ut, cstart->st1, cstart->st2, cstart->st3, &ts)/((cstart->n_subit+1.)*dt); */
  
/*   spline(cstart->n_time, 0, 0, 0., 0., ut, cstart->y_com, */
/* 	 cstart->st1, cstart->st2, cstart->st3, &flg); */
/*   fsi->S_new[4]=Sval(cstart->n_time,  u[0] */
/* 		     , ut, cstart->y_com, cstart->st1, cstart->st2, cstart->st3, &ts); */
/*   fsi->S_new[5]=deriv(cstart->n_time, u[0] */
/* 		      , ut, cstart->st1, cstart->st2, cstart->st3, &ts)/((cstart->n_subit+1.)*dt); */

  fsi->S_new[0]=x_midp[0];
  fsi->S_new[1]=(fsi->S_new[0]-fsi->S_real[0])/dt;
  
  fsi->S_new[4]=y_midp[0];
  fsi->S_new[5]=(fsi->S_new[4]-fsi->S_real[4])/dt;

  spline(cstart->n_time, 0, 0, 0., 0., ut, cstart->head_ang,
	 cstart->st1, cstart->st2, cstart->st3, &flg);
  fsi->S_ang_n[2]=Sval(cstart->n_time,  u[0]
		       , ut, cstart->head_ang, cstart->st1, cstart->st2, cstart->st3, &ts);
  fsi->S_ang_n[3]=(fsi->S_ang_n[2]-fsi->S_ang_r[2])/dt;

    //deriv(cstart->n_time, u[0]
    //			, ut, cstart->st1, cstart->st2, cstart->st3, &ts)/((cstart->n_subit+1.)*dt);
  PetscPrintf(PETSC_COMM_WORLD, "%d %d COM! x  %le %le %le z %le %le %le ang %le %le %le\n", ti, ts, fsi->S_new[0],fsi->S_new[1],(fsi->S_new[0]-fsi->S_real[0])/dt, fsi->S_new[4],fsi->S_new[5],(fsi->S_new[4]-fsi->S_real[4])/dt, fsi->S_ang_n[2],fsi->S_ang_n[3],(fsi->S_ang_n[2]-fsi->S_ang_r[2])/dt);
  }

  // correct for distance from midline
  for (i=0; i<ibm->n_v; i++) {
    u[1]=-(ibm->x_bp0[i]-0.5*(cstart->xmin+cstart->xmax));
    cl=sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
    dx[i] /= cl; //dx is n_y compnent of normal 
    dy[i] /= cl; //dy is n_z componet of normal
    ibm->x_bp[i] -=-u[1]+ u[1]*dy[i];
    ibm->z_bp[i] += u[1]*dx[i];
  }	

  PetscReal x_p, y_p, x_c=0.2, y_c=-0.6, rot_x=fsi->S_ang_n[2];
  if (nonIntertial) {
    for (i=0; i<ibm->n_v; i++) {
      /* R is rotation matrix R=[cos sin; -sin cos] 
	 position:   X= x_c + R(x_init - x_head - x_c)  
	 velocity: R dx_init = dX + R dx_head + dR R^T(X-x_c) 
	 dX        : velocity in the non-inertial frame
	 R dx_head : rotated velocity of the head in the intertial frame 
	 dR R^T(X-x_c): rotational velocity of the non-inertial frame 
      */

      // subtract position of head
      ibm->x_bp[i] -= fsi->S_new[0]-x_c;//cstart->x_com[ts];//-cstart->x_com[0] ;
      ibm->z_bp[i] -= fsi->S_new[4]-y_c;//cstart->y_com[ts];//-cstart->y_com[0] ;  

      // rotate base on -head angle around head at 0, 0
      x_p = ibm->x_bp[i];
      y_p = ibm->z_bp[i];
      
      ibm->x_bp[i] = x_c + (x_p-x_c)*cos(rot_x) + (y_p-y_c)*sin(rot_x);
      ibm->z_bp[i] = y_c - (x_p-x_c)*sin(rot_x) + (y_p-y_c)*cos(rot_x);

       
    }

    // rotate the com velcoity
    x_p =  fsi->S_new[1];
    y_p =  fsi->S_new[5];
    
    fsi->S_new[1] = (x_p)*cos(rot_x) + (y_p)*sin(rot_x);
    fsi->S_new[5] =-(x_p)*sin(rot_x) + (y_p)*cos(rot_x);
  }

  calc_ibm_normal(ibm);
  
  if (ti>0) calc_ibm_velocity(ibm, dt);
    
  //  if (ti == (ti/tiout)*tiout) 
  //  ibm_surface_out(&ibm[ibi],ti,ibi);

  if (print_spline) {
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(!rank) { // root processor write the data
      FILE *fd;
      char filen[80];
      sprintf(filen, "spline%5.5d.dat",ti);
      fd = fopen(filen, "w"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open spline.dat to write!!");
      for (i=0;i<1000; i++){
	PetscFPrintf(PETSC_COMM_SELF, fd, "%le %le\n", w[i],v[i]);
      }
      fclose(fd);
    }
  }

  PetscFree(dx);      PetscFree(dy); 
  PetscFree(ubar);    PetscFree(ut);
  PetscFree(x_midp);  PetscFree(y_midp);
  PetscFree(x_midt);  PetscFree(y_midt);

  return(0);
}
