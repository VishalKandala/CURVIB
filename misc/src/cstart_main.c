static char help[] = "Testing programming!";

#include "variables.h"
#include "petsc.h"

typedef struct {
  PetscInt   n_time, n_midp, n_subit; //number of time steps and number of control points
  PetscReal  *x_midp, *y_midp; //control points x/y_midp[time][n_midp]
  PetscReal  *x_com, *y_com, *head_ang; //position of center of mass and head angle
  PetscReal  *s1,*s2,*s3; /* spline coefficients for spatial interpolation
			     S(xx) = Y_mipd[i] + s1[i] * w + s2[i] * w**2 + s3[i] * w**3
			     where w = xx - x[i]
			  */
  PetscReal  *st1,*st2,*st3; /* spline coefficients for the temporal interpolation    
				S(xx) = Y_mipd[i] + s1[i] * w + s2[i] * w**2 + s3[i] * w**3
				where w = xx - x[i]
			     */
  PetscReal  xmin,xmax,ymin,ymax,zmin,zmax;//min max of ibm
} Cstart;

PetscReal CMx_c=0.,CMy_c=0.,CMz_c=0., L_dim=1.;
PetscInt  NumberOfBodies=1, tiout=1;
PetscInt  cop=0;


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
  PetscInt	rank, i;
  PetscInt   n_tstep, n_midp,n_subit, ts;
  PetscReal  *x_midp, *y_midp;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("midpoints.dat", "r"); if (!fd) SETERRQ(1, "Cannot open C-start node file!!! (midpoints.dat)");
    fscanf(fd, "%d %d %d\n", &n_tstep, &n_midp, &n_subit);

    //  PetscReal  x_midp[n_tstep][n_midp], y_midp[n_tstep][n_midp];
    
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &x_midp);
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &y_midp);
    PetscPrintf(PETSC_COMM_WORLD, "malloc %d %d! \n",n_tstep, n_midp);

    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_midp, 1, MPI_INT, 0, PETSC_COMM_WORLD);
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
    if (!fd) SETERRQ(1, "Cannot open C-start node file!!! (com_motion.dat)");

    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->x_com);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->y_com);
    PetscMalloc(n_tstep*sizeof(PetscReal), &cstart->head_ang);
    
    for (ts=0; ts<n_tstep; ts++) {
      fscanf(fd, "%le %le %le\n", &cstart->x_com[ts], &cstart->y_com[ts],
	     &cstart->head_ang[ts]);      
    }
    MPI_Bcast(cstart->x_midp, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->y_midp, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->head_ang, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
  } else {
    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_midp, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_SELF, "Bcast Rank 1! \n");
    //PetscReal  x_midp[n_tstep][n_midp], y_midp[n_tstep][n_midp];
    // PetscMalloc(n_midp*sizeof(PetscReal), &x_midp);
    //PetscMalloc(n_midp*sizeof(PetscReal), &y_midp);
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &x_midp);
    PetscMalloc(n_midp*n_tstep*sizeof(PetscReal), &y_midp);

    MPI_Bcast(x_midp,n_tstep*n_midp, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_midp,n_tstep*n_midp, MPIU_REAL, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_SELF, "Bcast Rank1! \n");

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

    MPI_Bcast(cstart->x_midp, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->y_midp, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(cstart->head_ang, n_tstep, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  return(0); 
}


int main(int argc, char **argv)
{  
  PetscInt   flg, i,seg,last,ts,subts,tss,ti,ibi;
  PetscReal  u[1001],v[1001],w[1001],*x_midp, *y_midp,*ubar,*dx,*dy;
  PetscReal  *x_midt, *y_midt,*ut,scale[7],cl, dt=0.001;
  PetscInt   rank;
  Cstart     cstart;
  IBMNodes   *ibm;
  PetscTruth print_spline=PETSC_FALSE;


  PetscInitialize(&argc, &argv, (char *)0, help);
  PetscPrintf(PETSC_COMM_WORLD, "Init \n");
  PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);

  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    ibm_read_Ansys(&ibm[ibi],ibi);

    ibm_scale(&ibm[ibi],&scale);

    cstart.zmin=scale[0];cstart.zmax=scale[1];
    cstart.ymin=scale[2];cstart.ymax=scale[3];
    cstart.xmin=scale[4];cstart.xmax=scale[5];
    PetscMalloc(ibm[ibi].n_v*sizeof(PetscReal), &dx);
    PetscMalloc(ibm[ibi].n_v*sizeof(PetscReal), &dy);
    PetscPrintf(PETSC_COMM_WORLD, "zmin/max %le %le y %le %le x %le %le\n",cstart.zmin,cstart.zmax,cstart.ymin,cstart.ymax,cstart.xmin,cstart.xmax);
  }

  read_midpoint(&cstart);

  PetscPrintf(PETSC_COMM_WORLD, "READ midpoints done! %d %d %d %le %le\n",cstart.n_time,cstart.n_midp, cstart.n_subit,cstart.x_midp[1],cstart.y_midp[1]);
  PetscMalloc(cstart.n_midp*sizeof(PetscReal), &ubar);
  PetscMalloc(cstart.n_midp*sizeof(PetscReal), &x_midp);
  PetscMalloc(cstart.n_midp*sizeof(PetscReal), &y_midp);
  PetscMalloc(cstart.n_time*sizeof(PetscReal), &ut);
  PetscMalloc(cstart.n_time*sizeof(PetscReal), &x_midt); 
  PetscMalloc(cstart.n_time*sizeof(PetscReal), &y_midt);
  
  for (tss=0; tss<cstart.n_time; tss++) {
    ut[tss]=tss*1.;
  }

  for (ti=0; ti<cstart.n_time+ (cstart.n_time-1)*(cstart.n_subit); ti++) {

    ts=ti/(1+cstart.n_subit);
    subts=ti-ts*(1+cstart.n_subit);
    PetscPrintf(PETSC_COMM_WORLD, "ti ts subts %d %d %d\n",ti,ts, subts);

    /* spline interpolation in time for finding midline points at
       the time between experimental time steps */

  for (seg=0; seg<cstart.n_midp; seg++) {
    ubar[seg]=seg*1.;

    if (subts) {
      for (tss=0; tss<cstart.n_time; tss++) {
	//ut[tss]=tss*1.;
	x_midt[tss]=cstart.x_midp[tss*cstart.n_midp+seg];
	y_midt[tss]=cstart.y_midp[tss*cstart.n_midp+seg];      
	//      PetscPrintf(PETSC_COMM_WORLD, "Sval %d %le %le %le\n",seg,ut[tss],x_midt[tss],y_midt[tss]);
      }
      
      u[0]=ts+(subts*1./(double)(cstart.n_subit+1.));
      spline(cstart.n_time, 0, 0, 0., 0., ut, x_midt, 
	     cstart.st1, cstart.st2, cstart.st3, &flg); 
      PetscPrintf(PETSC_COMM_WORLD, "Spline time intp done! %le %le %le %d\n",cstart.st1[3],cstart.st2[3],cstart.st3[3], flg);
      
      x_midp[seg]=Sval(cstart.n_midp, u[0]
			 , ut, x_midt, cstart.st1, cstart.st2, cstart.st3, &ts);   
      
      spline(cstart.n_time, 0, 0, 0., 0., ut, y_midt, 
	     cstart.st1, cstart.st2, cstart.st3, &flg); 
      y_midp[seg]=Sval(cstart.n_midp,  u[0]
			 , ut, y_midt, cstart.st1, cstart.st2, cstart.st3, &ts);
      
      // PetscPrintf(PETSC_COMM_WORLD, "Sval %d %le %le %le\n",ts,u[0],x_midp[seg],y_midp[seg]);
    } else {
      u[0]=ts+(subts*1./(double)(cstart.n_subit+1.));
      x_midp[seg]=cstart.x_midp[ts*cstart.n_midp+seg];
      y_midp[seg]=cstart.y_midp[ts*cstart.n_midp+seg];
      //PetscPrintf(PETSC_COMM_WORLD, "Sval %d %le %le %le \n",ts,ubar[seg],x_midp[seg],y_midp[seg]);
    }
  }
 
  /* spline for spatial interpolation on midpoints and move
     the mesh points similar to midline */

  // 1) in the y direction S(u)=Y(u) 
  spline(cstart.n_midp, 0, 0, 0., 0., ubar, y_midp, 
	 cstart.s1, cstart.s2, cstart.s3, &flg); 
  PetscPrintf(PETSC_COMM_WORLD, "Spline done! %le %le %le %d\n",cstart.s1[3],cstart.s2[3],cstart.s3[3], flg);

  if (print_spline) {
  for (i=0;i<1001; i++){
    u[i]=i*cstart.n_midp*0.001;
    v[i]=seval(cstart.n_midp, u[i], ubar, y_midp,
	       cstart.s1, cstart.s2, cstart.s3, &last);
    //PetscPrintf(PETSC_COMM_WORLD, "Sval %le %le \n",u[i],v[i]);
  }
  }

  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    for (i=0; i<ibm[ibi].n_v; i++) {
      u[0]=(ibm[ibi].z_bp0[i]-cstart.zmin)
	/  (  cstart.zmax    -cstart.zmin)*(cstart.n_midp-1.);
      seg=(int) (u[0]);
      ibm[ibi].z_bp[i]= //-cstart.y_midp[0]+
	Sval(cstart.n_midp, u[0]
	     , ubar, y_midp, cstart.s1, cstart.s2, cstart.s3, &seg); 
      dy[i]=deriv(cstart.n_midp, u[0]
		  , ubar, cstart.s1, cstart.s2, cstart.s3, &seg); 
      /*       if (i==1936) */
      /* 	PetscPrintf(PETSC_COMM_WORLD, "Sval %d %d %le %le\n",i,seg,u[1],deriv(cstart.n_midp, u[0] */
      /* 									      , ubar, cstart.s1, cstart.s2, cstart.s3, &seg)); */
    }
  }

  // 2) in the x direction S(u)=X(u) 
  spline(cstart.n_midp, 0, 0, 0., 0., ubar, x_midp, 
	 cstart.s1, cstart.s2, cstart.s3, &flg); 
  PetscPrintf(PETSC_COMM_WORLD, "Spline done! %le %le %le %d\n",cstart.s1[3],cstart.s2[3],cstart.s3[3], flg);

  if (print_spline) {
  for (i=0;i<1001; i++){
    u[i]=i*cstart.n_midp*0.001;
    w[i]=seval(cstart.n_midp, u[i], ubar, x_midp,
	       cstart.s1, cstart.s2, cstart.s3, &last);
    //PetscPrintf(PETSC_COMM_WORLD, "Sval %le %le \n",u[i],v[i]);
  }
  }
  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    for (i=0; i<ibm[ibi].n_v; i++) {
      u[0]=(ibm[ibi].z_bp0[i]-cstart.zmin)
	/  (  cstart.zmax    -cstart.zmin)*(cstart.n_midp-1.);
      seg=(int) (u[0]);
      ibm[ibi].x_bp[i]=ibm[ibi].x_bp0[i]+
	cstart.x_midp[0]-0.5*(cstart.xmin+cstart.xmax)+      
	Sval(cstart.n_midp, u[0]
	     , ubar, x_midp, cstart.s1, cstart.s2, cstart.s3, &seg)-cstart.x_midp[0]; 
      dx[i]=deriv(cstart.n_midp, u[0]
		  , ubar, cstart.s1, cstart.s2, cstart.s3, &seg); 
      /*  if (i==1936) */
      /* 	PetscPrintf(PETSC_COMM_WORLD, "Sval %d %d %le %le\n",i,seg,u[1],deriv(cstart.n_midp, u[0] */
      /* 	     , ubar, cstart.s1, cstart.s2, cstart.s3, &seg)); */
    }
  }
  

  /* spline interpolation in time for finding the center of mass & angle at
     the time between experimental time steps */
  PetscReal x_p, y_p, x_c=0., y_c=0., rot_x;
  PetscReal x_com, y_com, dxdt_com, dydt_com, drot_x;

  u[0]=ts+(subts*1./(double)(cstart.n_subit+1.));
  spline(cstart.n_time, 0, 0, 0., 0., ut, cstart.x_com, 
	 cstart.st1, cstart.st2, cstart.st3, &flg); 
  x_com=Sval(cstart.n_midp, u[0]
	     , ut, cstart.x_com, cstart.st1, cstart.st2, cstart.st3, &ts);   
  dxdt_com=deriv(cstart.n_midp, u[0]
		 , ut, cstart.st1, cstart.st2, cstart.st3, &ts)/((cstart.n_subit+1.)*dt);
  
  spline(cstart.n_time, 0, 0, 0., 0., ut, cstart.y_com, 
	 cstart.st1, cstart.st2, cstart.st3, &flg); 
  y_com=Sval(cstart.n_midp,  u[0]
	     , ut, cstart.y_com, cstart.st1, cstart.st2, cstart.st3, &ts);
  dydt_com=deriv(cstart.n_midp, u[0]
		 , ut, cstart.st1, cstart.st2, cstart.st3, &ts)/((cstart.n_subit+1.)*dt);

  spline(cstart.n_time, 0, 0, 0., 0., ut, cstart.head_ang, 
	 cstart.st1, cstart.st2, cstart.st3, &flg); 
  rot_x=Sval(cstart.n_midp,  u[0]
	     , ut, cstart.head_ang, cstart.st1, cstart.st2, cstart.st3, &ts);
  drot_x=deriv(cstart.n_midp, u[0]
	       , ut, cstart.st1, cstart.st2, cstart.st3, &ts)/((cstart.n_subit+1.)*dt);
  PetscPrintf(PETSC_COMM_WORLD, "COM! x  %le %le y %le %le ang %le %le\n",x_com,dxdt_com,y_com,dydt_com, rot_x,drot_x);

  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    // correct for distance from midline
    for (i=0; i<ibm[ibi].n_v; i++) {
      u[1]=-(ibm[ibi].x_bp0[i]-0.5*(cstart.xmin+cstart.xmax));
      cl=sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
      dx[i] /= cl;
      dy[i] /= cl;
      ibm[ibi].x_bp[i] -=-u[1]+ u[1]*dy[i];
      ibm[ibi].z_bp[i] += u[1]*dx[i];

      // rotate base on -head angle
/*       x_p = ibm[ibi].x_bp[i]; */
/*       y_p = ibm[ibi].z_bp[i]; */

/*       ibm[ibi].x_bp[i] = x_c + (x_p-x_c)*cos(rot_x) + (y_p-y_c)*sin(rot_x); */
/*       ibm[ibi].z_bp[i] = y_c - (x_p-x_c)*sin(rot_x) + (y_p-y_c)*cos(rot_x); */

/*       // add position of com */
/*       ibm[ibi].x_bp[i] += cstart.x_com[ts];//-cstart.x_com[0] ; */
/*       ibm[ibi].z_bp[i] += cstart.y_com[ts];//-cstart.y_com[0] ;  */   
    }	
    
    calc_ibm_normal(&ibm[ibi]);
    calc_ibm_velocity(&ibm[ibi], dt);

    for (i=0; i<ibm[ibi].n_v; i++) {
      ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
      ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
      ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];
      
      ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;
      ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
      ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;
      
      ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
      ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
      ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
    }
    
    ibm_surface_out(&ibm[ibi],ti,ibi);
  }

  if (print_spline) {
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor write the data
    FILE *fd;
    char filen[80];
    sprintf(filen, "spline%5.5d.dat",ti);
    fd = fopen(filen, "w"); if (!fd) SETERRQ(1, "Cannot open spline.dat to write!!");
    for (i=0;i<1000; i++){
      PetscFPrintf(PETSC_COMM_SELF, fd, "%le %le\n", w[i],v[i]);
    }
    fclose(fd);
  }
  }

  }//ti  
  PetscFinalize();
  return(0);
}
