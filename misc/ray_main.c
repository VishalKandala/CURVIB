static char help[] = "Testing programming!";

#include "variables.h"
#include "petsc.h"

PetscReal CMx_c=0.,CMy_c=0.,CMz_c=0., L_dim=1.;
PetscInt  NumberOfBodies=1, tiout=1;
PetscInt  cop=0, wing=0;

int main(int argc, char **argv)
{  
  PetscInt   flg, i,seg,last,ts,subts,tss,ti,ibi;
  PetscReal  u[1001],v[1001],w[1001],*x_midp, *y_midp,*ubar,*dx,*dy;
  PetscReal  *x_midt, *y_midt,*ut,scale[7],cl, dt=0.001,val;
  PetscInt   rank;
  Cstart     ray;
  IBMNodes   *ibm;
  PetscTruth print_spline=PETSC_FALSE;
  Cmpnts        pi;

  PetscInitialize(&argc, &argv, (char *)0, help);
  PetscPrintf(PETSC_COMM_WORLD, "Init \n");
  PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);

  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    ibm_read_Ansys(&ibm[ibi],ibi);
  }

  
  read_raykinematics(&ray);  
  scale_raykinematics(&ray);
  
  //  test_raykinematics(&ray); 

  Create_Interpolation_Matrix(&ray);
  Setup_Interpolation_Matrix(&ray, 0);

  for (ti=1;ti<ray.n_time;ti++) {
    PetscPrintf(PETSC_COMM_WORLD, "ti %d \n",ti);
    Get_Interpolation_Coefficients(&ray,ti);
/*     for (i=0; i<ray.n_midp; i++) { */
/*       pi.x=ray.x_midp[i]; */
/*       pi.y=ray.y_midp[i]; */
/*       extern PetscReal Surface_Interpolation_point(Cstart *ray, Cmpnts p); */
/*       val=Surface_Interpolation_point(&ray,pi); */
/*       PetscPrintf(PETSC_COMM_SELF, " initial %f interpolated %f !!!!\n", ray.z_midp[i],val); */
/*     } */
    Interpolate_ray_ibm(&ibm[0], &ray);
    calc_ibm_normal(&ibm[0]);
  
    ibm_surface_VTKOut(&ibm[0],0,ti);
  }
  MatDestroy(ray.Mphi);
  PetscPrintf(PETSC_COMM_WORLD, "done! ti %d \n",ti);
  PetscFinalize();
  return(0);
}
