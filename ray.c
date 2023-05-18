#include "variables.h"
#include "petsc.h"

PetscErrorCode read_raykinematics(Cstart *ray)
{

  PetscInt   rank, i;
  PetscInt   n_tstep, n_p,n_subit, ts;
  PetscReal  *x_p, *y_p, *z_p;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("stingraykin.dat", "r"); if (!fd) SETERRQ(1, "Cannot open C-start node file!!! (midpoints.dat)");
    fscanf(fd, "%d %d %d\n", &n_p, &n_tstep, &n_subit);
  
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &x_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &y_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &z_p);
    PetscPrintf(PETSC_COMM_WORLD, "malloc %d %d %d! \n",n_tstep, n_p, n_subit);

    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_p, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_subit, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Bcast! \n");
    for (ts=0; ts<n_tstep; ts++) {
      for (i=0; i<n_p; i++) {
	fscanf(fd, "%le %le %le ", &x_p[ts*n_p+i], &y_p[ts*n_p+i],&z_p[ts*n_p+i]);      
      }
      fscanf(fd, "\n");      
    }
    fclose(fd);
    ts=30; i=24;
    PetscPrintf(PETSC_COMM_WORLD, "READ! time %d point %d xyz %le %le %le\n",ts+1,i+1,x_p[ts*n_p+i], y_p[ts*n_p+i],z_p[ts*n_p+i]);
    MPI_Bcast(x_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(z_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Bcast! \n");


  } else {
    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_p, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_subit, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &x_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &y_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &z_p);

    MPI_Bcast(x_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(z_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

  ray->n_time=n_tstep;
  ray->n_subit=n_subit;
  ray->n_midp=n_p;
  ray->x_midp=x_p;
  ray->y_midp=y_p;
  ray->z_midp=z_p;
  PetscMalloc(n_p*sizeof(PetscReal), &ray->s1);
  PetscMalloc(n_p*sizeof(PetscReal), &ray->s2);
  PetscMalloc(n_p*sizeof(PetscReal), &ray->s3);

  return(0);
}


PetscErrorCode scale_raykinematics(Cstart *ray)
{
  PetscInt   i, ts, cl=135, xm,ym,zm;
  ts =0 ;
  xm=ray->x_midp[ts*ray->n_midp+13];
  ym=ray->y_midp[ts*ray->n_midp+13];
  zm=ray->z_midp[ts*ray->n_midp+13];

  //  for (ts=ray->n_time-1;ts<0; ts--) {
  for (ts=0;ts<ray->n_time; ts++) {
    for (i=0; i<ray->n_midp; i++) {
      ray->x_midp[ts*ray->n_midp+i] -= xm ;
      ray->y_midp[ts*ray->n_midp+i] -= ym ;
      ray->z_midp[ts*ray->n_midp+i] -= zm ;//;
    }
  }

  for (ts=0; ts<ray->n_time; ts++) {
    for (i=0; i<ray->n_midp; i++) {
      ray->x_midp[ts*ray->n_midp+i] /=cl;
      ray->y_midp[ts*ray->n_midp+i] /=cl;
      ray->z_midp[ts*ray->n_midp+i] /=cl;
    }
  }
  return(0);
}

PetscErrorCode test_raykinematics(Cstart *ray)
{

  PetscInt   rank, i,j;
  PetscInt   n_tstep, n_p,n_subit, ts;
  PetscReal  *x_p, *y_p, *z_p, r, tetha, pi=3.14159265;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    n_tstep=100;
    n_p=32;
    n_subit=0;
  
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &x_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &y_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &z_p);
    PetscPrintf(PETSC_COMM_WORLD, "malloc %d %d %d! \n",n_tstep, n_p, n_subit);
    
    // set up the test kinematics
    for (ts=0; ts<n_tstep; ts++) {
      for (j=0;j<4;j++) {
	for (i=0; i<8; i++) {
	  r=.5/(j+1);
	  tetha=2.*pi/8.*i;
	  x_p[ts*n_p+j*8+i]=r*cos(tetha);
	  y_p[ts*n_p+j*8+i]=r*sin(tetha);
	  z_p[ts*n_p+j*8+i]=0.2*sin(2*pi*x_p[ts*n_p+j*8+i]-2*pi*ts/50.);
	  if (ts==0) PetscPrintf(PETSC_COMM_WORLD, "ti %d n_p %d x_p %f y_p %f z_p %f! \n", ts,8*j+i, x_p[ts*n_p+j*4+i],y_p[ts*n_p+j*4+i],z_p[ts*n_p+j*4+i]);
	}
      }
    }
    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_p, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_subit, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Bcast! \n");


    MPI_Bcast(x_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(z_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Bcast! \n");


  } else {
    MPI_Bcast(&n_tstep, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_p, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&n_subit, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &x_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &y_p);
    PetscMalloc(n_p*n_tstep*sizeof(PetscReal), &z_p);

    MPI_Bcast(x_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(y_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(z_p, n_tstep*n_p, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

  ray->n_time=n_tstep;
  ray->n_subit=n_subit;
  ray->n_midp=n_p;
  ray->x_midp=x_p;
  ray->y_midp=y_p;
  ray->z_midp=z_p;
  PetscMalloc(n_p*sizeof(PetscReal), &ray->s1);
  PetscMalloc(n_p*sizeof(PetscReal), &ray->s2);
  PetscMalloc(n_p*sizeof(PetscReal), &ray->s3);

  return(0);
}


PetscReal phi_basis(Cmpnts p, Cmpnts pj)
{
  PetscReal r;
  r= ((p.x-pj.x)*(p.x-pj.x) +
      (p.y-pj.y)*(p.y-pj.y));
/*   // Thin-plate spline */
/*   if (r<1e-10)  */
/*     return(0); */
/*   else */
/*     return(r*log(sqrt(r))); */

  // cubic spline
  return(pow(sqrt(r),3));
}


PetscErrorCode Create_Interpolation_Matrix(Cstart *ray)
{
  PetscInt n_solu=ray->n_midp;

  MatCreateSeqDense(PETSC_COMM_SELF, n_solu, n_solu, PETSC_NULL, &(ray->Mphi));
  MatZeroEntries(ray->Mphi);
  return(0);
}
  
PetscErrorCode Setup_Interpolation_Matrix(Cstart *ray, PetscInt ts)
{
  PetscInt      i,j, n_solu;
  n_solu=ray->n_midp;
  PetscReal     val[ray->n_midp];
  PetscInt      idx[ray->n_midp],row;
  Cmpnts        p,pj;

  for (j=0; j<ray->n_midp; j++) {
    idx[j]=j;
  }

  for (i=0; i<ray->n_midp; i++) {
    p.x=ray->x_midp[ts*ray->n_midp+i];
    p.y=ray->y_midp[ts*ray->n_midp+i];

    for (j=0; j<ray->n_midp; j++) {
      pj.x=ray->x_midp[ts*ray->n_midp+j];
      pj.y=ray->y_midp[ts*ray->n_midp+j];

      val[j]=phi_basis(p,pj);
      if (fabs(val[j])<1e-3) PetscPrintf(PETSC_COMM_SELF, "val !!!! %d %d p_i %e %e p_j %e %e\n",i,j,p.x,p.y,pj.x,pj.y);
    }
    row=i;
    
    MatSetValues(ray->Mphi, 1, &row, n_solu, idx, val, INSERT_VALUES);
  }
  MatAssemblyBegin(ray->Mphi, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(ray->Mphi, MAT_FINAL_ASSEMBLY);
  // MatView(ray->Mphi,0);

  IS rperm, cperm;
  MatFactorInfo info;
  MatFactorInfoInitialize(&info);
  info.dtcol=1;
  MatGetOrdering(ray->Mphi,MATORDERING_NATURAL,&rperm, &cperm);
  MatLUFactor(ray->Mphi, rperm, cperm, &info);
  ISDestroy(rperm);
  ISDestroy(cperm);
  PetscPrintf(PETSC_COMM_SELF, "LU Factorization'!!!!\n");
  
  return(0);
}

PetscErrorCode Get_Interpolation_Coefficients(Cstart *ray, PetscInt ts)
{
  PetscInt i;
  Vec solu, rhs;
  PetscReal *z_p;
  // s1 is the coefficients
  // head_ang is the rhs to get the coefficients
  VecCreateSeq(PETSC_COMM_SELF,ray->n_midp,&solu);
  VecDuplicate(solu, &rhs);

  // z-direction ray->s1
  VecGetArray(rhs, &z_p);
  for (i=0; i<ray->n_midp; i++) {
    z_p[i]=ray->z_midp[ts*ray->n_midp+i];//-ray->z_midp[(ts-1)*ray->n_midp+i];
    //    PetscPrintf(PETSC_COMM_SELF, "Intepolation rhs %f !!!!\n", z_p[i]);
  }
  VecRestoreArray(rhs, &z_p);
  MatSolve(ray->Mphi, rhs, solu); 

  VecGetArray(solu, &z_p);
  for (i=0; i<ray->n_midp; i++) {
    ray->s1[i]=z_p[i];
    //    PetscPrintf(PETSC_COMM_SELF, "Intepolation coeff %f !!!!\n", ray->s1[i]);
  }
  VecRestoreArray(solu, &z_p);

  // y-direction ray->s2
  VecGetArray(rhs, &z_p);
  for (i=0; i<ray->n_midp; i++) {
    z_p[i]=ray->y_midp[ts*ray->n_midp+i];//-ray->y_midp[(ts-1)*ray->n_midp+i];
    //    PetscPrintf(PETSC_COMM_SELF, "Intepolation rhs %f !!!!\n", z_p[i]);
  }
  VecRestoreArray(rhs, &z_p);
  MatSolve(ray->Mphi, rhs, solu); 

  VecGetArray(solu, &z_p);
  for (i=0; i<ray->n_midp; i++) {
    ray->s2[i]=z_p[i];
    //    PetscPrintf(PETSC_COMM_SELF, "Intepolation coeff %f !!!!\n", ray->s1[i]);
  }
  VecRestoreArray(solu, &z_p);

  // x-direction ray->s3
  VecGetArray(rhs, &z_p);
  for (i=0; i<ray->n_midp; i++) {
    z_p[i]=ray->x_midp[ts*ray->n_midp+i];//-ray->x_midp[(ts-1)*ray->n_midp+i];
    //    PetscPrintf(PETSC_COMM_SELF, "Intepolation rhs %f !!!!\n", z_p[i]);
  }
  VecRestoreArray(rhs, &z_p);
  MatSolve(ray->Mphi, rhs, solu); 

  VecGetArray(solu, &z_p);
  for (i=0; i<ray->n_midp; i++) {
    ray->s3[i]=z_p[i];
    //    PetscPrintf(PETSC_COMM_SELF, "Intepolation coeff %f !!!!\n", ray->s1[i]);
  }
  VecRestoreArray(solu, &z_p);

  VecDestroy(solu);
  VecDestroy(rhs);
  return(0);
}

PetscReal Surface_Interpolation_point(Cstart *ray, PetscReal s1[],Cmpnts p)
{
  Cmpnts        pi;
  PetscInt      i;
  PetscReal     val=0.,phi;

  for (i=0; i<ray->n_midp; i++) {
    pi.x=ray->x_midp[i];
    pi.y=ray->y_midp[i];
    phi  = phi_basis(p,pi);
    val+=s1[i]*phi;
    //   PetscPrintf(PETSC_COMM_SELF, "Intepolation coeff %f phi %f val %f!!!!\n", ray->s1[i], phi,val);
  }
  return(val);
}

PetscErrorCode Interpolate_ibm(IBMNodes *ibm, Cstart *ray)
{
  PetscInt      i;
  Cmpnts p;

  for (i=0; i<ibm->n_v; i++) {
    p.x=ibm->x_bp[i];
    p.y=ibm->y_bp[i];

    ibm->z_bp[i] = ibm->z_bp0[i]+Surface_Interpolation_point(ray,ray->s1,p);
    //    PetscPrintf(PETSC_COMM_SELF, " initial %f interpolated %f !!!!\n", ibm->z_bp[i],ibm->z_bp0[i]);
  }
  return(0);
}

PetscErrorCode Interpolate_ray_ibm(IBMNodes *ibm, Cstart *ray)
{
  PetscInt      i;
  Cmpnts p;

  for (i=0; i<ibm->n_v; i++) {
    p.x=ibm->x_bp[i];
    p.y=fabs(ibm->y_bp[i]);

    ibm->z_bp[i] = ibm->z_bp0[i]+Surface_Interpolation_point(ray,ray->s1,p);
    // ibm->y_bp[i] = ibm->y_bp0[i]+Surface_Interpolation_point(ray,ray->s2,p);
    // ibm->x_bp[i] = ibm->x_bp0[i]+Surface_Interpolation_point(ray,ray->s3,p);
    //    PetscPrintf(PETSC_COMM_SELF, " initial %f interpolated %f !!!!\n", ibm->z_bp[i],ibm->z_bp0[i]);
  }
  return(0);
}

