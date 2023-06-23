#include "variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c, L_dim;
extern PetscInt  cop, wing,NumberOfBodies,rheology;
extern char orient[];
PetscReal l_x,l_y,l_z;

PetscErrorCode ibm_surface_out(IBMNodes *ibm, PetscInt ti,
			       PetscInt ibi)
{
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == ti) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d_%2.2d.dat",ti,ibi);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,u_x,u_y,u_z,n_x,n_y,n_z,nt_x,nt_y,nt_z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-6]=NODAL,[7-12]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);
      for (i=0; i<ibm->n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
      }
      for (i=0; i<ibm->n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
      }
      for (i=0; i<ibm->n_v; i++) {	
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
      }
      for (i=0; i<ibm->n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->u[i].x);
      }
      for (i=0; i<ibm->n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->u[i].y);
      }
      for (i=0; i<ibm->n_v; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->u[i].z);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]);
      }
      for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }
      fclose(f);
    } 
  }
  return(0);
}


PetscErrorCode ibm_surface_VTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti)
{
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surface%2.2d_%5.5d.vtk", ibi,ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
//    PetscPrintf(PETSC_COMM_WORLD,"Surface file headers done!\n"); 
    PetscPrintf(PETSC_COMM_WORLD,"Number of nodes is %d \n",ibm->n_v);
    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibm->n_v);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->x_bp[i],ibm->y_bp[i],ibm->z_bp[i]);
    }
//    PetscPrintf(PETSC_COMM_WORLD,"Surface File no.of nodes and co-ordinates printed!\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibm->n_elmt, (n_cells+1)*ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);
    }
//    PetscPrintf(PETSC_COMM_WORLD,"Surface File no.of elements and nodes printed!\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
    }
//    PetscPrintf(PETSC_COMM_WORLD,"Surface File no.of nodes and co-ordinates printed!\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", ibm->n_v);
    PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->u[i].x,ibm->u[i].y,ibm->u[i].z);
    }
//    PetscPrintf(PETSC_COMM_WORLD,"Surface File point data printed!\n");

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n", ibm->n_elmt);
    PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nf float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nf_x[i],ibm->nf_y[i],ibm->nf_z[i]);
    }
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS nt float\n"); */
/*     for (i=0; i<ibm->n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nt_x[i],ibm->nt_y[i],ibm->nt_z[i]); */
/*     } */
    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS p double\n");    
    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
//    PetscPrintf(PETSC_COMM_WORLD,"Surface File cell data printed!\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", ibm->nt_y[i]);
    }
    fclose(f);
  }
  return(0);
}

PetscErrorCode ibm_volume_VTKOut(IBMVNodes *ibmv, PetscInt ibi, PetscInt ti)//added by Qiping
{
    // vtk file name
  PetscInt n_cells=4;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "volume%2.2d_%5.5d.vtk", ibi,ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Volume Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
  
    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibmv->n_v);
    for (i=0; i<ibmv->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibmv->x_bp[i],ibmv->y_bp[i],ibmv->z_bp[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibmv->n_elmt, (n_cells+1)*ibmv->n_elmt);
    for (i=0; i<ibmv->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d %d\n",n_cells, ibmv->nv1[i], ibmv->nv2[i], ibmv->nv3[i], ibmv->nv4[i]);
     
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibmv->n_elmt);
    for (i=0; i<ibmv->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",10);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", ibmv->n_v);
    PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
    for (i=0; i<ibmv->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibmv->u[i].x,ibmv->u[i].y,ibmv->u[i].z);
    }
    fclose(f);
  }
  return(0);
}

PetscErrorCode ibm_read(IBMNodes *ibm1, IBMNodes *ibm2)
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
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscReal     cl=30.;

  double xt;
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata0", "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file");
    n_v =0;
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%le", &xt);
    ibm1->n_v = n_v/2;
    ibm2->n_v = n_v/2;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->urm1));

    MPI_Bcast(&(ibm2->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->urm1));

    for (i=0; i<n_v; i++) {
      fscanf(fd, "%le %le %le %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i], &t, &t, &t);
      x_bp[i] = x_bp[i] / cl;// 28.
      y_bp[i] = y_bp[i] / cl;
      z_bp[i] = z_bp[i] / cl;
    }
/*     ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */   

    PetscReal temp;
    for (i=0; i<n_v/2; i++) {
      temp = y_bp[i];
      ibm1->x_bp0[i] = z_bp[i];
      ibm1->y_bp0[i] = x_bp[i]-CMy_c;
      ibm1->z_bp0[i] = -temp+CMz_c;

      ibm1->x_bp[i] = z_bp[i];
      ibm1->y_bp[i] = x_bp[i]-CMy_c;
      ibm1->z_bp[i] = -temp+CMz_c;

      ibm1->x_bp_o[i] = z_bp[i];
      ibm1->y_bp_o[i] = x_bp[i]-CMy_c;
      ibm1->z_bp_o[i] = -temp+CMz_c;   

      ibm1->u[i].x = 0.;
      ibm1->u[i].y = 0.;
      ibm1->u[i].z = 0.;
      
      ibm1->uold[i].x = 0.;
      ibm1->uold[i].y = 0.;
      ibm1->uold[i].z = 0.;
      
      ibm1->urm1[i].x = 0.;
      ibm1->urm1[i].y = 0.;
      ibm1->urm1[i].z = 0.;   
    }

    for (i=0; i<n_v/2; i++) {
      temp = y_bp[i+n_v/2];
      ibm2->x_bp0[i] = z_bp[i+n_v/2];
      ibm2->y_bp0[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp0[i] = -temp+CMz_c;

      ibm2->x_bp[i] = z_bp[i+n_v/2];
      ibm2->y_bp[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp[i] = -temp+CMz_c;

      ibm2->x_bp_o[i] = z_bp[i+n_v/2];
      ibm2->y_bp_o[i] = x_bp[i+n_v/2]+CMy_c;
      ibm2->z_bp_o[i] = -temp+CMz_c;

      ibm2->u[i].x = 0.;
      ibm2->u[i].y = 0.;
      ibm2->u[i].z = 0.;
      
      ibm2->uold[i].x = 0.;
      ibm2->uold[i].y = 0.;
      ibm2->uold[i].z = 0.;
      
      ibm2->urm1[i].x = 0.;
      ibm2->urm1[i].y = 0.;
      ibm2->urm1[i].z = 0.;   
    }

/*     ibm1->x_bp=ibm1->x_bp0;ibm1->y_bp=ibm1->y_bp0;ibm1->z_bp=ibm1->z_bp0; */
/*     ibm2->x_bp=ibm2->x_bp0;ibm2->y_bp=ibm2->y_bp0;ibm2->z_bp=ibm2->z_bp0; */
/*     ibm1->x_bp_o=ibm1->x_bp0;ibm1->y_bp_o=ibm1->y_bp0;ibm1->z_bp_o=ibm1->z_bp0; */
/*     ibm2->x_bp_o=ibm2->x_bp0;ibm2->y_bp_o=ibm2->y_bp0;ibm2->z_bp_o=ibm2->z_bp0; */

    MPI_Bcast(ibm1->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    fscanf(fd, "%i\n", &n_elmt);
    ibm1->n_elmt = n_elmt/2;
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm2->n_elmt = n_elmt/2;
    //MPI_Bcast(&(ibm2->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "elmts , %d \n", n_elmt);

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    // Added 4/1/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area
    
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
    
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_z));

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_z));

    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i %i %i\n", nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      //      PetscPrintf(PETSC_COMM_WORLD, "I %d %d %d\n", nv1[i], nv2[i], nv3[i]);
    }

    //ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

    fclose(fd);

    for (i=0; i<n_elmt; i++) {
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;

      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) 
      {
	if (nf_z[i]>0) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
	} else {
	ns_x[i] = -1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = -1.;
	nt_z[i] = 0.;
	}	  
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
/*       ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.; */
/*       ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.; */
/*       ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	 */
    }
   
    //ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    for (i=0; i<n_elmt/2; i++) {
      ibm1->nv1[i]=nv1[i];
      ibm1->nv2[i]=nv2[i];
      ibm1->nv3[i]=nv3[i];

      ibm1->nf_x[i]=  nf_z[i];
      ibm1->nf_y[i]=  nf_x[i];
      ibm1->nf_z[i]= -nf_y[i];

      ibm1->ns_x[i]=  ns_z[i];
      ibm1->ns_y[i]=  ns_x[i];
      ibm1->ns_z[i]= -ns_y[i];

      ibm1->nt_x[i]=  nt_z[i];
      ibm1->nt_y[i]=  nt_x[i];
      ibm1->nt_z[i]= -nt_y[i];
    }

    for (i=0; i<n_elmt/2; i++) {
      ibm2->nv1[i]=nv1[i+n_elmt/2]-n_v/2;
      ibm2->nv2[i]=nv2[i+n_elmt/2]-n_v/2;
      ibm2->nv3[i]=nv3[i+n_elmt/2]-n_v/2;

      ibm2->nf_x[i]=  nf_z[i+n_elmt/2];
      ibm2->nf_y[i]=  nf_x[i+n_elmt/2];
      ibm2->nf_z[i]= -nf_y[i+n_elmt/2];

      ibm2->ns_x[i]=  ns_z[i+n_elmt/2];
      ibm2->ns_y[i]=  ns_x[i+n_elmt/2];
      ibm2->ns_z[i]= -ns_y[i+n_elmt/2];

      ibm2->nt_x[i]=  nt_z[i+n_elmt/2];
      ibm2->nt_y[i]=  nt_x[i+n_elmt/2];
      ibm2->nt_z[i]= -nt_y[i+n_elmt/2];

    }

/*     for (i=0; i<n_elmt; i++) { */
/*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
/*     } */

    MPI_Bcast(ibm1->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm1->n_v = n_v/2;

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm1->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm1->urm1));

    MPI_Bcast(&(ibm2->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp0));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp0));

    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->x_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->y_bp_o));
    PetscMalloc(n_v/2*sizeof(PetscReal), &(ibm2->z_bp_o));

    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->u));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->uold));
    PetscMalloc(n_v/2*sizeof(Cmpnts), &(ibm2->urm1));

    for (i=0; i<n_v/2; i++) {

      ibm1->u[i].x = 0.;
      ibm1->u[i].y = 0.;
      ibm1->u[i].z = 0.;
      
      ibm1->uold[i].x = 0.;
      ibm1->uold[i].y = 0.;
      ibm1->uold[i].z = 0.;
      
      ibm1->urm1[i].x = 0.;
      ibm1->urm1[i].y = 0.;
      ibm1->urm1[i].z = 0.;   

      ibm2->u[i].x = 0.;
      ibm2->u[i].y = 0.;
      ibm2->u[i].z = 0.;
      
      ibm2->uold[i].x = 0.;
      ibm2->uold[i].y = 0.;
      ibm2->uold[i].z = 0.;
      
      ibm2->urm1[i].x = 0.;
      ibm2->urm1[i].y = 0.;
      ibm2->urm1[i].z = 0.;   

    }

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm1->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp0, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->x_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->y_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->z_bp_o, n_v/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     ibm1->x_bp=ibm1->x_bp0;ibm1->y_bp=ibm1->y_bp0;ibm1->z_bp=ibm1->z_bp0; */
/*     ibm2->x_bp=ibm2->x_bp0;ibm2->y_bp=ibm2->y_bp0;ibm2->z_bp=ibm2->z_bp0; */
/*     ibm1->x_bp_o=ibm1->x_bp0;ibm1->y_bp_o=ibm1->y_bp0;ibm1->z_bp_o=ibm1->z_bp0; */
/*     ibm2->x_bp_o=ibm2->x_bp0;ibm2->y_bp_o=ibm2->y_bp0;ibm2->z_bp_o=ibm2->z_bp0; */


    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm1->n_elmt = n_elmt/2;
    ibm2->n_elmt = n_elmt/2;
    
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm1->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm1->nt_z));

    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv1));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv2));
    PetscMalloc(n_elmt/2*sizeof(PetscInt), &(ibm2->nv3));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nf_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->ns_z));

    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_x));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_y));
    PetscMalloc(n_elmt/2*sizeof(PetscReal), &(ibm2->nt_z));

/*     ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3; */
/*     ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z; */


/*     for (i=0; i<n_elmt/2; i++) { */
/*       ibm1->nv1[i]=nv1[i]; */
/*       ibm1->nv2[i]=nv2[i]; */
/*       ibm1->nv3[i]=nv3[i]; */

/*       ibm1->nf_x[i]= nf_x[i]; */
/*       ibm1->nf_y[i]= nf_y[i]; *//*       ibm1->nf_z[i]= nf_z[i]; */
/*     } */

/*     for (i=0; i<n_elmt/2; i++) { */
/*       ibm2->nv1[i]=nv1[i+n_elmt/2]; */
/*       ibm2->nv2[i]=nv2[i+n_elmt/2]; */
/*       ibm2->nv3[i]=nv3[i+n_elmt/2]; */

/*       ibm2->nf_x[i]= nf_x[i+n_elmt/2]; */
/*       ibm2->nf_y[i]= nf_y[i+n_elmt/2]; */
/*       ibm2->nf_z[i]= nf_z[i+n_elmt/2]; */
/*     } */

/*     for (i=0; i<n_elmt; i++) { */
    i=10;
      PetscPrintf(PETSC_COMM_SELF, " ibm1 xbp %d %le %le %le\n", i, ibm1->z_bp0[i], ibm1->z_bp_o[i], ibm1->z_bp[i]);
/*     } */

    MPI_Bcast(ibm1->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm1->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm1->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nv1, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv2, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nv3, n_elmt/2, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nf_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nf_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->nt_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->nt_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm2->ns_x, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_y, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm2->ns_z, n_elmt/2, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}

PetscErrorCode ibm_read_thin(IBMNodes *ibm)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
/*   PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270]; */
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscReal	t;
 
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("NodeData_A1_m", "r");
    
    n_v = 3270;
    ibm->n_v = n_v;

    MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);

/*     PetscMalloc(n_v*101*sizeof(PetscReal), &x_bp_in); */
/*     PetscMalloc(n_v*101*sizeof(PetscReal), &y_bp_in); */
/*     PetscMalloc(n_v*101*sizeof(PetscReal), &z_bp_in); */

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscInt t1, j;
    for (j=0; j<101; j++)
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%d %le %le %le %le", &t1, &(ibm->x_bp_in[j][i]),
	       &(ibm->y_bp_in[j][i]), &(ibm->z_bp_in[j][i]), &t);
	ibm->x_bp_in[j][i] = ibm->x_bp_in[j][i] / 24.;
	ibm->y_bp_in[j][i] = ibm->y_bp_in[j][i] / 24.;
	ibm->z_bp_in[j][i] = ibm->z_bp_in[j][i] / 24.;
/*       PetscPrintf(PETSC_COMM_WORLD, "%le %le %le\n", ibm->x_bp_in[j][i], ibm->y_bp_in[j][i], ibm->z_bp_in[j][i]); */
    }

    for (j=45; j<50; j++) {
      for (i=0; i<n_v; i++) {
	ibm->x_bp_in[j][i] = ibm->x_bp_in[50-j][i];
	ibm->y_bp_in[j][i] = ibm->y_bp_in[50-j][i];
	ibm->z_bp_in[j][i] = ibm->z_bp_in[50-j][i];
      }
    }
	
/*     ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

/*     for (i=0; i<n_v; i++) { */
/*       PetscReal temp; */
/*       temp = ibm->y_bp0[i]; */
/*       ibm->y_bp0[i] = ibm->z_bp0[i]; */
/*       ibm->z_bp0[i] = -temp; */
/*     } */


    MPI_Bcast(ibm->x_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */

/*     fscanf(fd, "%i\n", &n_elmt); */
    fclose(fd);

    n_elmt = 6228;
    ibm->n_elmt = n_elmt;
    MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_z);
    fd = fopen("Elements_id_A1", "r");
    
    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i, %i, %i, %i\n", &t1, nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      //      PetscPrintf(PETSC_COMM_WORLD, "I %d %d %d\n", nv1[i], nv2[i], nv3[i]);
    }
    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

    fclose(fd);

    for (i=0; i<n_elmt; i++) {
      nv1[i+n_elmt] = nv1[i];
      nv2[i+n_elmt] = nv3[i];
      nv3[i+n_elmt] = nv2[i];
    }

    n_elmt = n_elmt*2;
    ibm->n_elmt = n_elmt;
/*     for (i=0; i<n_elmt; i++) { */
/*       n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i]; */
/*       dx12 = x_bp[n2e] - x_bp[n1e]; */
/*       dy12 = y_bp[n2e] - y_bp[n1e]; */
/*       dz12 = z_bp[n2e] - z_bp[n1e]; */

/*       dx13 = x_bp[n3e] - x_bp[n1e]; */
/*       dy13 = y_bp[n3e] - y_bp[n1e]; */
/*       dz13 = z_bp[n3e] - z_bp[n1e]; */

/*       nf_x[i] = dy12 * dz13 - dz12 * dy13; */
/*       nf_y[i] = -dx12 * dz13 + dz12 * dx13; */
/*       nf_z[i] = dx12 * dy13 - dy12 * dx13; */

/*       dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]); */

/*       nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr; */

      
/*     } */
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

/*     for (i=0; i<n_elmt; i++) { */
/*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
/*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    n_v = ibm->n_v;
/*     PetscMalloc(n_v*sizeof(PetscReal), &x_bp); */
/*     PetscMalloc(n_v*sizeof(PetscReal), &y_bp); */
/*     PetscMalloc(n_v*sizeof(PetscReal), &z_bp); */

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

/*     PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0)); */
/*     PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0)); */
/*     PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0)); */

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

/*     MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */

    MPI_Bcast(ibm->x_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_in, n_v*101, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    n_elmt = ibm->n_elmt;
/*     ibm->n_elmt = n_elmt; */

    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*2*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*2*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    n_elmt = 2 * n_elmt;
    ibm->n_elmt = n_elmt;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}

PetscErrorCode ibm_read_tecplot(IBMNodes *ibm, PetscInt nt, PetscBool flg)
{
PetscInt	rank;
PetscInt	n_v , n_elmt ;
PetscReal	*x_bp , *y_bp , *z_bp ;
PetscInt	*nv1 , *nv2 , *nv3 ;
PetscReal	*nf_x, *nf_y, *nf_z;
PetscInt	i;
PetscInt	n1e, n2e, n3e;
PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
PetscReal	dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;

PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT print out \n ");

  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


 if(!rank) 
  { // root processor read in the data
    FILE *fd;
    char filen[128];
    sprintf(filen, "IBMSolidData_%1.1d.dat", nt);

    fd = fopen(filen, "r");
    if (!fd) 
    {
      PetscPrintf(PETSC_COMM_WORLD, "Cannot open IBM Solid node file %d\n", nt);
      SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file IBMSolidData.dat!");
     }
    n_v =0;

    if (fd) 
     {

      fscanf(fd, "%d %d",&n_v,&n_elmt);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;

      PetscPrintf(PETSC_COMM_WORLD," nv = %d ne = %d ...\n",n_v, n_elmt);

      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
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

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
      
      PetscReal cl = 1.;
      PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);
      //cl=1./L_dim;

      // Loop around all the elements - node points
      for (i=0; i<n_v; i++) 
      {
	fscanf(fd, "%le %le %le", &x_bp[i], &y_bp[i], &z_bp[i]);
      	
	x_bp[i] = x_bp[i] / cl + CMx_c;
	y_bp[i] = y_bp[i] / cl + CMy_c;
	z_bp[i] = z_bp[i] / cl + CMz_c;
	
	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];


	ibm->x_bp0[i] = x_bp[i];
	ibm->y_bp0[i] = y_bp[i];
	ibm->z_bp0[i] = z_bp[i];

	ibm->x_bp_o[i] = x_bp[i];
	ibm->y_bp_o[i] = y_bp[i];
	ibm->z_bp_o[i] = z_bp[i];

	ibm->u[i].x = 0.;
	ibm->u[i].y = 0.;
	ibm->u[i].z = 0.;

	ibm->uold[i].x = 0.;
	ibm->uold[i].y = 0.;
	ibm->uold[i].z = 0.;

	ibm->urm1[i].x = 0.;
	ibm->urm1[i].y = 0.;
	ibm->urm1[i].z = 0.;
      }//End of i loop
      
      PetscPrintf(PETSC_COMM_WORLD, "IBM tec READ\n");
    
      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

 /*      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*       MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*       MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD); */
      
      MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

      // Added 4/1/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      // Added 6/4/06 iman
      //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added

      //Starting to scan all the elements
      for (i=0; i<n_elmt; i++) 
      {

	if (flg) {//3
	  fscanf(fd, "%i %i %i\n", nv2+i, nv1+i, nv3+i);
	  if (i==0) PetscPrintf(PETSC_COMM_WORLD, "IBM tec READ elemt REVVVV!!!\n");
	}
	else 
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
  
  
    // Now loop around all the elements
    // Find the normal vector to the surface
    for (i=0; i<n_elmt; i++) 
     {

      n1e = nv1[i]; 
      n2e = nv2[i]; 
      n3e = nv3[i];


      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;

            // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
    }
    PetscPrintf(PETSC_COMM_WORLD, "IBM tec clac n\n");

    ibm->nf_x = nf_x; 
    ibm->nf_y = nf_y;  
    ibm->nf_z = nf_z;

    //Added 4/1/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

    PetscPrintf(PETSC_COMM_WORLD, "IBM tec clac n\n");

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    // Added 4/1/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/2/06 iman
    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 6/4/06 iman
    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    PetscInt ti=0;
    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",ti,nt);
    fd = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<n_v; i++) {	
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->z_bp[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
    }
    fclose(fd);

  }
 else if (rank) //If it is root processor
   {

    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    for (i=0; i<ibm->n_v; i++) 
    {
      ibm->x_bp[i] = ibm->x_bp0[i];
      ibm->y_bp[i] = ibm->y_bp0[i];
      ibm->z_bp[i] = ibm->z_bp0[i];

      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;

      ibm->uold[i].x = 0.;
      ibm->uold[i].y = 0.;
      ibm->uold[i].z = 0.;
      
      ibm->urm1[i].x = 0.;
      ibm->urm1[i].y = 0.;
      ibm->urm1[i].z = 0.;     
    }

    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    //Added 4/1/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

    //Added 4/2/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
    // Added 4/2/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

    // Added 6/4/06
    //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

    ibm->nv1 = nv1; 
    ibm->nv2 = nv2; 
    ibm->nv3 = nv3;

    ibm->nf_x = nf_x; 
    ibm->nf_y = nf_y; 
    ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    //Added 4/2/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
   }//If it is root 


     PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT successfully \n ");

  return(0);

}

 
PetscErrorCode ibm_placement(IBMNodes *ibm, FSInfo *fsi,
			     PetscInt Nibm_y, PetscInt Nibm_z,
			     PetscReal dYibm, PetscReal dZibm, 
			     PetscInt ibi)
{
  PetscInt	n_v=ibm->n_v , n_elmt=ibm->n_elmt, n_y, n_z, i ;

  PetscInt ti=0,rank;
  FILE *f;
  char filen[80];
  

  
  n_z= ibi/Nibm_y;
  n_y= ibi - (int) (n_z*Nibm_y);
  
  fsi->y_c += n_y*dYibm;
  fsi->z_c += n_z*dZibm;
 
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) {
 
    PetscPrintf(PETSC_COMM_WORLD, "IBM Placement ibi %d %d %d, n_y %d, n_z %d, y_c %le, z_c %le %le %le\n ", ibi, n_v, n_elmt,n_y, n_z, fsi->y_c, fsi->z_c,  n_y*dYibm,  n_z*dZibm);
    
    for (i=0; i<n_v; i++) {      
      ibm->x_bp[i] += 0.;//0.25 ;// 24.;	
      ibm->y_bp[i] += n_y*dYibm;//+ ibi*2.;//5.;//2 .;//8.;//6.;//2.   ;// 24.;
      ibm->z_bp[i] += n_z*dZibm;//2.;//8.;//15.;//2.   ;// 24.;
      
      ibm->x_bp0[i] = ibm->x_bp[i];
      ibm->y_bp0[i] = ibm->y_bp[i];
      ibm->z_bp0[i] = ibm->z_bp[i];
      
      ibm->x_bp_o[i] = ibm->x_bp[i];
      ibm->y_bp_o[i] = ibm->y_bp[i];
      ibm->z_bp_o[i] = ibm->z_bp[i];
    }
    
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",ti,ibi);
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
    
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt); */
/*     for (i=0; i<n_v; i++) { */
      
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*     } */
    fclose(f);
  } else {
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
  }
  return(0);
}

PetscErrorCode ibm_read_ucd(IBMNodes *ibm, PetscInt ibi)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;

  char   ss[20];
  //double xt;
  char string[128];

  PetscReal cl = 1.;
  cl=1./L_dim;
  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, PETSC_NULL);     


  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ ibmdata\n");
    char filen[80];  
    sprintf(filen,"ibmdata%2.2d" , ibi);
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file");
    n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      
      fscanf(fd, "%i %i %i %i %i",&n_v,&n_elmt,&ii,&ii,&ii);
      PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d\n",n_v, n_elmt);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;      
      
      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
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
      
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
      /* // */
      PetscPrintf(PETSC_COMM_WORLD, "Orienttation %-s !\n",orient);
      /* // */
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	
/* //Rotation of MHV */
	PetscReal pi = 3.14159265359, angle, R[4][4], x[4], xn[4], u, v, w, a, b, c, xd[3];
	PetscInt m,n;

	xd[0] = 0.02;  xd[1] = -0.05;  xd[2] = -0.12; //Translation
	
	if (strcmp(orient,"xx00") == 0) {
	  angle =0.;
	} else if (strcmp(orient,"xy45") == 0) {
	  angle = pi/4.;
	} else if (strcmp(orient,"yy00") == 0) {
	  angle = pi/2.;
	}

	u = 0.; v = 0.; w = 1.; //Rotation vector
	a = 0.04; b = 0.04; c = 0.134; //Rotation origin
	x[0] = x_bp[i]; x[1] = y_bp[i]; x[2] = z_bp[i]; x[3] = 1.;
	
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
	x_bp[i] = xn[0]; y_bp[i] = xn[1]; z_bp[i] = xn[2];
       
        x_bp[i] = x_bp[i]/cl + CMx_c + xd[0];//0.25 ;// 24.;
        y_bp[i] = y_bp[i]/cl + CMy_c + xd[1];//+ ibi*2.;//5.;//2.;//8.;//6.;//2.   ;// 24.;
        z_bp[i] = z_bp[i]/cl + CMz_c + xd[2];//+ ibi*2.;//2.;//8.;//15.;//2.   ;// 24.;
/* // */
	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

	ibm->x_bp0[i] = x_bp[i];
	ibm->y_bp0[i] = y_bp[i];
	ibm->z_bp0[i] = z_bp[i];

	ibm->x_bp_o[i] = x_bp[i];
	ibm->y_bp_o[i] = y_bp[i];
	ibm->z_bp_o[i] = z_bp[i];

	ibm->u[i].x = 0.;
	ibm->u[i].y = 0.;
	ibm->u[i].z = 0.;

	ibm->uold[i].x = 0.;
	ibm->uold[i].y = 0.;
	ibm->uold[i].z = 0.;

	ibm->urm1[i].x = 0.;
	ibm->urm1[i].y = 0.;
	ibm->urm1[i].z = 0.;
      }
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

/*       ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */
/*       ibm->x_bp_o = x_bp; ibm->y_bp_o = y_bp; ibm->z_bp_o = z_bp; */

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      // Added 4/1/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area

      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      // Added 6/4/06 iman
      //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added

      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, ss, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;

      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

      fclose(fd);
    }
      
    for (i=0; i<n_elmt; i++) {
      
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];
      
      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];
      
      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;
      
      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
      // Temp sol. 2D
/*       if (fabs(nf_x[i])<.5) */
/* 	nf_x[i]=0.; */

      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
    }
    
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
    //Added 4/1/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    
    
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/1/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/2/06 iman
    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 6/4/06 iman
    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

 /*    PetscFree(dA); */
/*     PetscFree(nf_x);PetscFree(nf_y);PetscFree(nf_z); */
/*     PetscFree(nt_x);PetscFree(nt_y);PetscFree(nt_z); */
/*     PetscFree(ns_x);PetscFree(ns_y);PetscFree(ns_z); */
/*     PetscFree(nv1);PetscFree(nv2);PetscFree(nv3); */
/*     PetscFree(x_bp);PetscFree(y_bp);PetscFree(z_bp); */
    PetscInt ti=0;
    FILE *f;
    //char filen[80];
    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",ti,ibi);
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
    
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt); */
/*     for (i=0; i<n_v; i++) { */
      
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*     } */
    fclose(f);

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;
    
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
    
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;

      ibm->uold[i].x = 0.;
      ibm->uold[i].y = 0.;
      ibm->uold[i].z = 0.;
      
      ibm->urm1[i].x = 0.;
      ibm->urm1[i].y = 0.;
      ibm->urm1[i].z = 0.;      
    }
        
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    //Added 4/1/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

    //Added 4/2/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
    // Added 4/2/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

    // Added 6/4/06
    //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    //Added 4/2/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

  return(0);
}

PetscErrorCode ibm_read_ucd_old(IBMNodes *ibm)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal	dr;
  PetscInt 	temp;
  
  char string[128];
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata", "r");
    if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file");
    n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);

      fscanf(fd, "%i %i %i %i %i\n", &n_v, &n_elmt, &temp, &temp, &temp);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;

      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

      
      PetscReal cl = 1.;

      PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);
      
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &temp, &x_bp[i], &y_bp[i], &z_bp[i]);
	x_bp[i] = (x_bp[i]    ) / cl;
	y_bp[i] = (y_bp[i]+6. ) / cl ;
	z_bp[i] = (z_bp[i]+15.) / cl;
	
/* 	ibm->x_bp[i] = x_bp[i]; */
/* 	ibm->y_bp[i] = y_bp[i]; */
/* 	ibm->z_bp[i] = z_bp[i]; */

	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

	ibm->u[i].x = 0.;
	ibm->u[i].y = 0.;
	ibm->u[i].z = 0.;
      }
      ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp;

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

	
      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      char str[20];
      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &temp, &temp, str, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;

      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      fclose(fd);
    }
    for (i=0; i<n_elmt; i++) {
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;

      
    }
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

    /*     for (i=0; i<n_elmt; i++) { */
    /*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
    /*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    for (i=0; i<ibm->n_v; i++) {
      ibm->x_bp[i] = ibm->x_bp0[i];
      ibm->y_bp[i] = ibm->y_bp0[i];
      ibm->z_bp[i] = ibm->z_bp0[i];

      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}
PetscErrorCode ibm_duplicate(IBMNodes *ibm,IBMNodes *ibm0,PetscReal r_x, PetscReal r_y, PetscReal r_z)
{
 
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;

  PetscReal     *dA ;
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
 
 
  ibm->n_v = ibm0->n_v;
  n_v=ibm0->n_v;           
       
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
  
  PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
  PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
  PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
  
  
  for (i=0;i<n_v;i++){
    ibm->x_bp[i] = ibm0->x_bp0[i]+r_x;
    ibm->y_bp[i] = ibm0->y_bp0[i]+r_y;
    ibm->z_bp[i] = ibm0->z_bp0[i]+r_z;
    
    ibm->x_bp0[i] = ibm0->x_bp0[i]+r_x;
    ibm->y_bp0[i] = ibm0->y_bp0[i]+r_y;
    ibm->z_bp0[i] = ibm0->z_bp0[i]+r_z;
    
    ibm->x_bp_o[i] = ibm0->x_bp0[i]+r_x;
    ibm->y_bp_o[i] = ibm0->x_bp0[i]+r_y;
    ibm->z_bp_o[i] = ibm0->x_bp0[i]+r_z;
  
    
    ibm->u[i].x = 0.;
    ibm->u[i].y = 0.;
    ibm->u[i].z = 0.;
    
    ibm->uold[i].x = 0.;
    ibm->uold[i].y = 0.;
    ibm->uold[i].z = 0.;
    
    ibm->urm1[i].x = 0.;
    ibm->urm1[i].y = 0.;
    ibm->urm1[i].z = 0.;      
  }
  
    
  ibm->n_elmt = ibm0->n_elmt;      
  n_elmt= ibm0->n_elmt;
  
  
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
  
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
  
  PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area
  
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);
  
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
  PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
  
  
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
  
  
    
  for (i=0;i<n_elmt;i++){
    nv1[i] = ibm0->nv1[i];
    nv2[i] = ibm0->nv2[i];
    nv3[i] = ibm0->nv3[i];
  }
  ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3; 
  x_bp=ibm->x_bp;  y_bp=ibm->y_bp ; z_bp = ibm->z_bp ;
  
    for (i=0; i<n_elmt; i++) {
      
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];
      
      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];
      
      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;
      
      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
      
      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	
      
    }
    
   
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
    //Added 4/1/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;  
    
    return(0);
}

PetscErrorCode  Random_Particle(UserCtx *user, IBMNodes *ibm,  FSInfo *fsi)
{
  PetscOptionsGetReal(PETSC_NULL, "-L_x", &l_x, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_y", &l_y, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_z", &l_z, PETSC_NULL);
  PetscInt N_x=2,N_y=2,N_z=2;
  PetscOptionsGetInt(PETSC_NULL, "-No_x", &N_x, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-No_y", &N_y, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-No_z", &N_z, PETSC_NULL);
  PetscReal s1,s2,s3;
  PetscReal P_x,P_y,P_z,O_x,O_y,O_z,D_x,D_y,D_z;
  time_t t;
  PetscReal	*x_bp =ibm[0].x_bp, *y_bp = ibm[0].y_bp, *z_bp = ibm[0].z_bp;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscInt      i,j,k,n_v,n,n_elmt;

  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;
  
 
  n_v = ibm[0].n_v;
  for(i=0; i<n_v; i++) {
    
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);
    
    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);
    
    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }
 
  D_x= xbp_max- xbp_min;
  D_y= ybp_max- ybp_min;
  D_z= zbp_max- zbp_min;
  D_x +=0.01;
  D_y +=0.01;
  D_z +=0.01;

  PetscPrintf(PETSC_COMM_WORLD, "Random particle --- xbp_min %le ybp_min %le zbp_min %le \n",xbp_min,ybp_min,zbp_min);
  srand((unsigned) time(&t));
 
  s1 = rand()/((double) RAND_MAX+1);
  s2 = rand()/((double) RAND_MAX+1);
  s3 = rand()/((double) RAND_MAX+1);
  P_x=s1*(l_x/N_x-D_x)-xbp_min;
  P_y=s2*(l_y/N_y-D_y)-ybp_min;
  P_z=s3*(l_z/N_z-D_z)-zbp_min;
 /*  P_x=0.10; */
/*   P_y=0.10; */
/*   P_z=0.10; */
  // PetscPrintf(PETSC_COMM_WORLD, "Random particle --- P_x %le P_y %le P_z %le \n",P_x,P_y,P_z);

  for (i=0; i<n_v;i++) {

    ibm[0].x_bp[i] =  ibm[0].x_bp0[i]+P_x;
    ibm[0].y_bp[i] =  ibm[0].y_bp0[i]+P_y;
    ibm[0].z_bp[i] =  ibm[0].z_bp0[i]+P_z;
       
    ibm[0].x_bp_o[i] = ibm[0].x_bp0[i]+P_x;
    ibm[0].y_bp_o[i] = ibm[0].y_bp0[i]+P_y;
    ibm[0].z_bp_o[i] = ibm[0].z_bp0[i]+P_z;
  }
  n_elmt = ibm[0].n_elmt;
  
  for (i=0; i<n_elmt; i++) {
    ibm[0].cent_x[i] += P_x;
    ibm[0].cent_y[i] += P_y;
    ibm[0].cent_z[i] += P_z;	
  }

  for(k=0;k<N_z;k++){
    for(j=0;j<N_y;j++){
      for (i=0;i<N_x;i++){
  
	s1 = rand()/((double) RAND_MAX+1);
	s2 = rand()/((double) RAND_MAX+1);
	s3 = rand()/((double) RAND_MAX+1);
	//	PetscPrintf(PETSC_COMM_WORLD,"random value %le %le %le \n",s1,s2,s3);
	O_x=i*l_x/N_x;
	O_y=j*l_y/N_y;
	O_z=k*l_z/N_z;
	P_x=O_x+s1*(l_x/N_x-D_x)-xbp_min;
	P_y=O_y+s2*(l_y/N_y-D_y)-ybp_min;
	P_z=O_z+s3*(l_z/N_z-D_z)-zbp_min;

	n=i+N_x*j+N_x*N_y*k;
/* 	if (n==1) { */
/* 	  P_x=1.01; */
/* 	  P_y=0.09; */
/* 	  P_z=0.09; */
/* 	} */
/* 	if (n==2) { */
/* 	  P_x=0.1; */
/* 	  P_y=1.00; */
/* 	  P_z=0.11; */
/* 	} */
/* 	if (n==3) { */
/* 	  P_x=1.01; */
/* 	  P_y=1.01; */
/* 	  P_z=0.09; */
/* 	} */
/* 	if (n==4) { */
/* 	  P_x=0.09; */
/* 	  P_y=0.11; */
/* 	  P_z=1.00; */
/* 	} */
/* 	if (n==5) { */
/* 	  P_x=0.09; */
/* 	  P_y=1.02; */
/* 	  P_z=1.01; */
/* 	} */
/* 	if (n==6) { */
/* 	  P_x=1.01; */
/* 	  P_y=0.08; */
/* 	  P_z=1.00; */
/* 	} */
/* 	if (n==7) { */
/* 	  P_x=1.00; */
/* 	  P_y=1.01; */
/* 	  P_z=1.02; */
/* 	} */
	if (n>0) ibm_duplicate(&ibm[n],&ibm[0],P_x,P_y,P_z);	  
      }
    }
  }
 

  return 0;
}

PetscErrorCode  Duplicate_Particle(UserCtx *user, IBMNodes *ibm,  FSInfo *fsi)
{
 

  PetscInt i=0;
  PetscOptionsGetReal(PETSC_NULL, "-L_x", &l_x, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_y", &l_y, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_z", &l_z, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD,"ksi eta zeta length  are %le %le %le \n",l_x,l_y,l_z);

  PetscPrintf(PETSC_COMM_WORLD,"original partcile  %d  , n_v %d , n_elmt %d \n",i,ibm[i].n_v,ibm[i].n_elmt); 

  i=1;
  ibm_duplicate(&ibm[i],&ibm[0], l_x, 0.0, 0.0);
  fsi[i].x_c=fsi[0].x_c+l_x;
  fsi[i].y_c=fsi[0].y_c;
  fsi[i].z_c=fsi[0].z_c;

  fsi[i].a_c[0]= fsi[i].x_c;
  fsi[i].a_c[1]= fsi[i].y_c;
  fsi[i].a_c[2]= fsi[i].z_c;

  PetscPrintf(PETSC_COMM_WORLD,"partcile copy %d, center of mass  x_c %le y_c %le z_c %le \n",i,fsi[i].x_c,fsi[i].y_c,fsi[i].z_c);

  i=2;
  ibm_duplicate(&ibm[i],&ibm[0], -l_x, 0.0, 0.0);
  fsi[i].x_c=fsi[0].x_c-l_x;
  fsi[i].y_c=fsi[0].y_c;
  fsi[i].z_c=fsi[0].z_c;


  fsi[i].a_c[0]= fsi[i].x_c;
  fsi[i].a_c[1]= fsi[i].y_c;
  fsi[i].a_c[2]= fsi[i].z_c;

  PetscPrintf(PETSC_COMM_WORLD,"partcile copy %d, center of mass  x_c %le y_c %le z_c %le \n",i,fsi[i].x_c,fsi[i].y_c,fsi[i].z_c);

  i=3;
  ibm_duplicate(&ibm[i],&ibm[0], 0.0, 0.0, l_z);
  fsi[i].x_c=fsi[0].x_c;
  fsi[i].y_c=fsi[0].y_c;
  fsi[i].z_c=fsi[0].z_c+l_z;

  fsi[i].a_c[0]= fsi[i].x_c;
  fsi[i].a_c[1]= fsi[i].y_c;
  fsi[i].a_c[2]= fsi[i].z_c;

  PetscPrintf(PETSC_COMM_WORLD,"partcile copy %d, center of mass  x_c %le y_c %le z_c %le \n",i,fsi[i].x_c,fsi[i].y_c,fsi[i].z_c);

  i=4;
  ibm_duplicate(&ibm[i],&ibm[0], 0.0, 0.0, -l_z);
  fsi[i].x_c=fsi[0].x_c;
  fsi[i].y_c=fsi[0].y_c;
  fsi[i].z_c=fsi[0].z_c-l_z;

  fsi[i].a_c[0]= fsi[i].x_c;
  fsi[i].a_c[1]= fsi[i].y_c;
  fsi[i].a_c[2]= fsi[i].z_c;

  PetscPrintf(PETSC_COMM_WORLD,"partcile copy %d, canter of mass  x_c %le y_c %le z_c %le \n",i,fsi[i].x_c,fsi[i].y_c,fsi[i].z_c);
 
 

  return(0);
}


PetscErrorCode ibm_read_Ansys(IBMNodes *ibm, PetscInt ibi)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscInt      itr;
 
  PetscReal     cl=1.,cl1=1.,cl2=1,cl3=1;//168. for copepod;
  PetscReal X,Y,Z;
  PetscReal theta,u_x,u_y,u_z;

  theta=0.0;u_x=0.0;u_y=0.0;u_z=0.0;
  if (rheology){
    PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm1", &cl1, PETSC_NULL);  
    PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm2", &cl2, PETSC_NULL); 
    PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm3", &cl3, PETSC_NULL);    
    PetscOptionsGetReal(PETSC_NULL, "-theta", &theta, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-u_x", &u_x, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-u_y", &u_y, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-u_z", &u_z, PETSC_NULL);
    
    if (ibi==1){
      cl1=1.75;
      cl2=1.75;
      cl3=1.75;
      u_x=1.0;
      theta=3.1415/3;
    }else if (ibi==2){
      cl1=3.8;
      cl2=4.3;
      cl3=4.7;
      u_x=1.0;
      u_y=0.0;
      u_z=0.0;
      theta=3.1415/1.8; 
    }
    PetscPrintf(PETSC_COMM_WORLD, "cl1 %le cl2 %le cl3 %le theta %le n_x %le n_y %le n_z %le CMx %le CMy %le CMz %le \n",cl1,cl2,cl3,theta,u_x,u_y,u_z,CMx_c,CMy_c,CMz_c);
  }else{
    PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "cl %le CMx %le CMy %le CMz %le \n",cl1,cl2,cl3,theta,u_x,u_y,u_z,CMx_c,CMy_c,CMz_c);
  }

  char string[128];
  

  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // rot processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ nlist, %le\n", L_dim);
    char filen[80];  
    sprintf(filen,"nlist%2.2d" , ibi);

    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file");
 
   
    if (fd) {
      
      fscanf(fd, "%i",&n_v);
      PetscPrintf(PETSC_COMM_SELF, "number of nodes of list %d %d \n",ibi, n_v);
      
      ibm->n_v = n_v;
            
      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
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
      
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

      i=-1;
      fgets(string, 128, fd);// miss 2 lines
      while(fgets(string, 128, fd)) {

	itr = 0;
	
	// read the data on 20 lines
	while (i+1<n_v && itr<20) {
	itr++;	i++;

	fscanf(fd, "%d %d %d %le %le %le\n", &ii, &ii, &ii,&(x_bp[i]), &(y_bp[i]), &(z_bp[i]));

	X=x_bp[i];
	Y=y_bp[i];
	Z=z_bp[i];
	if (rheology){
	  x_bp[i] = ((cos(theta)+u_x*u_x*(1-cos(theta)))*X+(u_x*u_y*(1-cos(theta))-u_z*sin(theta))*Y+(u_x*u_z*(1-cos(theta))+u_y*sin(theta))*Z)/cl1+CMx_c;
	  y_bp[i] = ((u_x*u_y*(1-cos(theta))+u_z*sin(theta))*X+(cos(theta)+u_y*u_y*(1-cos(theta)))*Y+(u_y*u_z*(1-cos(theta))-u_x*sin(theta))*Z)/cl2+CMy_c;
	  z_bp[i] = ((u_x*u_z*(1-cos(theta))-u_y*sin(theta))*X+(u_y*u_z*(1-cos(theta))+u_x*sin(theta))*Y+(cos(theta)+u_z*u_z*(1-cos(theta)))*Z)/cl3+CMz_c;
	}else{

	  x_bp[i] = X/cl+CMx_c;
	  y_bp[i] = Y/cl+CMy_c;
	  z_bp[i] = Z/cl+CMz_c;
	}
	
	if (cop) {
	ibm->x_bp[i] = y_bp[i]*L_dim;
	ibm->y_bp[i] = z_bp[i]*L_dim;
	ibm->z_bp[i] = x_bp[i]*L_dim;

	ibm->x_bp0[i] = x_bp[i]*L_dim;
	ibm->y_bp0[i] = y_bp[i]*L_dim;
	ibm->z_bp0[i] = z_bp[i]*L_dim;

	ibm->x_bp_o[i] = y_bp[i]*L_dim;
	ibm->y_bp_o[i] = z_bp[i]*L_dim;
	ibm->z_bp_o[i] = x_bp[i]*L_dim;
	} else {
	ibm->x_bp[i] = x_bp[i]*L_dim;
	ibm->y_bp[i] = y_bp[i]*L_dim;
	ibm->z_bp[i] = z_bp[i]*L_dim;

	ibm->x_bp0[i] = x_bp[i]*L_dim;
	ibm->y_bp0[i] = y_bp[i]*L_dim;
	ibm->z_bp0[i] = z_bp[i]*L_dim;

	ibm->x_bp_o[i] = x_bp[i]*L_dim;
	ibm->y_bp_o[i] = y_bp[i]*L_dim;
	ibm->z_bp_o[i] = z_bp[i]*L_dim;
	}

	ibm->u[i].x = 0.;
	ibm->u[i].y = 0.;
	ibm->u[i].z = 0.;

	ibm->uold[i].x = 0.;
	ibm->uold[i].y = 0.;
	ibm->uold[i].z = 0.;
      
	ibm->urm1[i].x = 0.;
	ibm->urm1[i].y = 0.;
	ibm->urm1[i].z = 0.;      
	}
      }
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);
      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);
      i=n_v-1;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      fclose(fd);
    }

    PetscPrintf(PETSC_COMM_SELF, "READ elist\n");

    sprintf(filen,"elist%2.2d" , ibi);
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file");

    if (fd) {
           
      fscanf(fd, "%i",&n_elmt);
      PetscPrintf(PETSC_COMM_SELF, "number of element of list %d %d \n",ibi, n_elmt);
      
      ibm->n_elmt = n_elmt;      
         
      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      // Added 4/1/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      // Added 6/4/06 iman
    
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added
      
      
      i=0;
      fgets(string, 128, fd);// miss 2 lines
      while(fgets(string, 128, fd)) {
	itr = 0;	
	// read the data on 20 lines
	while (i<n_elmt && itr<20) {
	itr++;	i++;

	fscanf(fd, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", &ii,&ii, &ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii, &nv1[i-1], &nv2[i-1], &nv3[i-1],&ii);
	nv1[i-1] = nv1[i-1] - 1; nv2[i-1] = nv2[i-1]-1; nv3[i-1] = nv3[i-1] - 1;

	}
      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);
      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);
      i=n_elmt-1;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

      fclose(fd);
    }

    x_bp=ibm->x_bp;  y_bp=ibm->y_bp ; z_bp = ibm->z_bp ;
    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n");      
    for (i=0; i<n_elmt; i++) {
           
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];
      
      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];
      
      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;
      
      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
      if (i==10484)  PetscPrintf(PETSC_COMM_WORLD, "n1e %d n2e %d n3e %d\n",n1e,n2e,n3e);

      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
      
      //Added 4/1/06 iman
      dA[i] = dr/2.; 
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;
      // if (i==10484)  PetscPrintf(PETSC_COMM_WORLD, "centx %le centy %le centz %le\n",ibm->cent_x[i],ibm->cent_y[i],ibm->cent_z[i]);	
    }
          
    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n"); 

    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
    //Added 4/1/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;  

    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n"); 
    
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/1/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/2/06 iman
    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 6/4/06 iman
    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n"); 

// write the surface
    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",0,ibi);
    fd = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->z_bp[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
    }


    fclose(fd);


  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;
    
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
    
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;

      ibm->uold[i].x = 0.;
      ibm->uold[i].y = 0.;
      ibm->uold[i].z = 0.;
         
      ibm->urm1[i].x = 0.;
      ibm->urm1[i].y = 0.;
      ibm->urm1[i].z = 0.;      
    }
        
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    //Added 4/1/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

    //Added 4/2/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
    // Added 4/2/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

    // Added 6/4/06
   
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    //Added 4/2/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

  return(0);
}


PetscErrorCode ibm_read_Ansys_vol(IBMVNodes *ibmv, PetscInt ibi)// added by Qiping
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ,*nv4;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e, n4e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13, dx14, dy14, dz14;
  PetscReal     dr;

  PetscInt      itr;

  PetscReal     cl1=1.,cl2=1,cl3=1;//168. for copepod;
  
  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm1", &cl1, PETSC_NULL);  
  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm2", &cl2, PETSC_NULL); 
  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm3", &cl3, PETSC_NULL);    

  PetscReal X,Y,Z;
  PetscReal theta=0.0,u_x=0.0,u_y=0.0,u_z=0.0;

  PetscOptionsGetReal(PETSC_NULL, "-theta", &theta, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-u_x", &u_x, PETSC_NULL);  
  PetscOptionsGetReal(PETSC_NULL, "-u_y", &u_y, PETSC_NULL); 
  PetscOptionsGetReal(PETSC_NULL, "-u_z", &u_z, PETSC_NULL);

 if (ibi==1){
    cl1=1.75;
    cl2=1.75;
    cl3=1.75;
    u_x=1.0;
    theta=3.1415/3;
  }else if (ibi==2){
    cl1=3.8;
    cl2=4.3;
    cl3=4.7;
    u_x=1.0;
    u_y=0.0;
    u_z=0.0;
    theta=3.1415/1.8;
 }

/*   PetscReal     cl=1.;//168. for copepod; */

/*   PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, PETSC_NULL);     */ 
  //char   ss[20];
  //double xt;
  char string[128];

  // Temp Sol. Change based on file !!!

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ nlist, %le\n", L_dim);
    char filen[80];  
    // sprintf(filen,"NLIST.lis" );
    sprintf(filen,"NLIST%2.2d.lis", ibi);
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file");

    if (fd) {
      
      fscanf(fd, "%i",&n_v);
      PetscPrintf(PETSC_COMM_SELF, "number of nodes of list %d %d \n",ibi, n_v);
      
      ibmv->n_v = n_v;
      MPI_Bcast(&(ibmv->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->z_bp));
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->z_bp_o));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->z_bp0));
      
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibmv->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibmv->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibmv->urm1));

      for (i=0; i<n_v; i++) {
	ibmv->u[i].x = 0.;
	ibmv->u[i].y = 0.;
	ibmv->u[i].z = 0.;
	
	ibmv->uold[i].x = 0.;
	ibmv->uold[i].y = 0.;
	ibmv->uold[i].z = 0.;
	
	ibmv->urm1[i].x = 0.;
	ibmv->urm1[i].y = 0.;
	ibmv->urm1[i].z = 0.;      
      }

      i=-1;
      fgets(string, 128, fd);// miss 2 lines
      while(fgets(string, 128, fd)) {
	itr = 0;
	
	// read the data on 20 lines
	while (i+1<n_v && itr<20) {
	itr++;	i++;
	fscanf(fd, "%d %le %le %le\n", &ii, &(x_bp[i]), &(y_bp[i]), &(z_bp[i]));


	X=x_bp[i];
	Y=y_bp[i];
	Z=z_bp[i];


/* 	x_bp[i] = X/cl1 + CMx_c; */
/* 	y_bp[i] = Y/cl2 + CMy_c; */
/* 	z_bp[i] = Z/cl3 + CMz_c; */


	x_bp[i] = ((cos(theta)+u_x*u_x*(1-cos(theta)))*X+(u_x*u_y*(1-cos(theta))-u_z*sin(theta))*Y+(u_x*u_z*(1-cos(theta))+u_y*sin(theta))*Z)/cl1+CMx_c;
	y_bp[i] = ((u_x*u_y*(1-cos(theta))+u_z*sin(theta))*X+(cos(theta)+u_y*u_y*(1-cos(theta)))*Y+(u_y*u_z*(1-cos(theta))-u_x*sin(theta))*Z)/cl2+CMy_c;
	z_bp[i] = ((u_x*u_z*(1-cos(theta))-u_y*sin(theta))*X+(u_y*u_z*(1-cos(theta))+u_x*sin(theta))*Y+(cos(theta)+u_z*u_z*(1-cos(theta)))*Z)/cl3+CMz_c;

	ibmv->x_bp[i] = x_bp[i]*L_dim;
	ibmv->y_bp[i] = y_bp[i]*L_dim;
	ibmv->z_bp[i] = z_bp[i]*L_dim;

	ibmv->x_bp0[i] = x_bp[i]*L_dim;
	ibmv->y_bp0[i] = y_bp[i]*L_dim;
	ibmv->z_bp0[i] = z_bp[i]*L_dim;

	ibmv->x_bp_o[i] = x_bp[i]*L_dim;
	ibmv->y_bp_o[i] = y_bp[i]*L_dim;
	ibmv->z_bp_o[i] = z_bp[i]*L_dim;

    	}
      }

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);
      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);
      i=n_v-1;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

      MPI_Bcast(ibmv->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibmv->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibmv->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      MPI_Bcast(ibmv->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibmv->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibmv->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      MPI_Bcast(ibmv->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibmv->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibmv->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      fclose(fd);
    }

    PetscPrintf(PETSC_COMM_SELF, "READ elist\n");
    sprintf(filen,"ELIST%2.2d.lis" , ibi);   
 //  sprintf(filen,"ELIST.lis" , ibi);
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file");

    if (fd) {
      fscanf(fd, "%i",&n_elmt);
      PetscPrintf(PETSC_COMM_SELF, "number of element of list %d %d \n",ibi, n_elmt);
      
      ibmv->n_elmt = n_elmt;      
      // n_v =0;
      
      MPI_Bcast(&(ibmv->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv4);
      
     
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->dV0)); 
     
     
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->cent_z));
      // end added
      
      
      i=0;
      fgets(string, 128, fd);// miss 2 lines
   
      while(fgets(string, 128, fd)) {
	itr = 0;	
	// read the data on 20 lines
	while (i<n_elmt && itr<20) {
	itr++;	i++;

	fscanf(fd, "%d %d %d %d %d %d %d %d %d %d\n", &ii,&ii, &ii,&ii,&ii,&ii, &nv1[i-1], &nv2[i-1], &nv3[i-1],&nv4[i-1]);
	nv1[i-1] = nv1[i-1] - 1; nv2[i-1] = nv2[i-1]-1; nv3[i-1] = nv3[i-1] - 1; nv4[i-1]=nv4[i-1]-1;
	}
      }
      ibmv->nv1 = nv1; ibmv->nv2 = nv2; ibmv->nv3 = nv3; ibmv->nv4 = nv4;

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d %d\n", nv1[i], nv2[i], nv3[i], nv4[i]);
      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d %d\n", nv1[i], nv2[i], nv3[i], nv4[i]);
      i=n_elmt-1;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d %d\n", nv1[i], nv2[i], nv3[i], nv4[i]);

      fclose(fd);
    }

    x_bp=ibmv->x_bp;  y_bp=ibmv->y_bp ; z_bp = ibmv->z_bp ;
    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n");      
    for (i=0; i<n_elmt; i++) {
      
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i]; n4e=nv4[i];
      
	  
      dx12 = x_bp[n1e] - x_bp[n2e];
      dy12 = y_bp[n1e] - y_bp[n2e];
      dz12 = z_bp[n1e] - z_bp[n2e];
      
      dx13 = x_bp[n1e] - x_bp[n3e];
      dy13 = y_bp[n1e] - y_bp[n3e];
      dz13 = z_bp[n1e] - z_bp[n3e];

      dx14 = x_bp[n1e] - x_bp[n4e];
      dy14 = y_bp[n1e] - y_bp[n4e];
      dz14 = z_bp[n1e] - z_bp[n4e];
      
      dr = dx12*(dy13*dz14-dz13*dy14)-dy12*(dx13*dz14-dz13*dx14)+dz12*(dx13*dy14-dy13*dx14);
     
      ibmv->dV0[i]=fabs(dr/6.0);

      // Added 6/4/06 iman
      // Calc the center of the element
      ibmv->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e]+x_bp[n4e])/4.;
      ibmv->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e]+y_bp[n4e])/4.;
      ibmv->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e]+z_bp[n4e])/4.;	
    }
          
      
    MPI_Bcast(ibmv->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->nv4, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
  
    MPI_Bcast(ibmv->dV0, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    PetscPrintf(PETSC_COMM_WORLD, "cop mmmmmnf!\n"); //in order to determine location

// write the surface
    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",0,ibi);
    fd = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables=x,y,z\n");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibmv->x_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibmv->y_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibmv->z_bp[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%d %d %d %d\n", ibmv->nv1[i]+1, ibmv->nv2[i]+1, ibmv->nv3[i]+1, ibmv->nv4[i]+1);
    }


    fclose(fd);


  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibmv->n_v = n_v;
    
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->z_bp));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->z_bp0));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibmv->z_bp_o));
    
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibmv->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibmv->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibmv->urm1));

    for (i=0; i<n_v; i++) {
      ibmv->u[i].x = 0.;
      ibmv->u[i].y = 0.;
      ibmv->u[i].z = 0.;

      ibmv->uold[i].x = 0.;
      ibmv->uold[i].y = 0.;
      ibmv->uold[i].z = 0.;
         
      ibmv->urm1[i].x = 0.;
      ibmv->urm1[i].y = 0.;
      ibmv->urm1[i].z = 0.;      
    }
        
    MPI_Bcast(ibmv->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
    MPI_Bcast(ibmv->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibmv->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibmv->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibmv->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv4);
  
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->dV0));

    ibmv->nv1 = nv1; ibmv->nv2 = nv2; ibmv->nv3 = nv3; ibmv->nv4 = nv4;
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->cent_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->cent_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibmv->cent_z));

    MPI_Bcast(ibmv->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->nv4, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
  
    MPI_Bcast(ibmv->dV0, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibmv->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibmv->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

  return(0);
}

//------------------------------------------------------------------
//ibm_read_Ansys is updated for icem Hafez 1/20/14
PetscErrorCode ibm_read_Icem(IBMNodes *ibm, PetscInt ibi)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscInt      itr;
  PetscReal     tt;
  PetscReal     cl=1.;//168. for copepod;

  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, PETSC_NULL);     
  
  char string[128];

 

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ nlist begins\n");
    char filen[80];  
    sprintf(filen,"nlist%2.2d" , ibi);


    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file");
 
    PetscPrintf(PETSC_COMM_WORLD, "Geometric Center of Immersed Boundary:  CMx_c = %le,CMy_c = %le,CMz_c = %le \n",CMx_c,CMy_c,CMz_c);
    
    if (fd) {
      
      fscanf(fd, "%i",&n_v);  // Read the first line of nlist, to get number of nodes
      PetscPrintf(PETSC_COMM_SELF, "number of nodes of list %d %d \n",ibi, n_v);
      
      ibm->n_v = n_v;   // No of nodes is designated as n_v. 
          
      
      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);  // Total number of nodes is broadcast to all processors (since the file reading is done by only the master processor).
    
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);  // Memory allocation for x,y,z co-ordinate values.
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));  // Memory allocation for x,y,z co-ordinate values in the datastructure ibm.
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o)); // Memory allocation for x,y,z co-ordinate values in ibm, for co-ordinates at previous timestep.
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0)); // Memory allocation for x,y,z co-ordinate values in ibm, for co-ordinates at initial timestep.
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
      
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));       // u,u_old and urm-1 velocities with three components will be created (memory allocated), these are used for euler-lagrange modelling(FSI).
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));

   
      i=-1;
         fgets(string,128, fd);// skip line one

    	while (i+1<n_v) {  // as number of lines in nlist is equal to n_v+1 and we skipped one line, we start i from -1  and go through each line to get co-ordinates for each node.

         i++;
	 fscanf(fd, "%d %d %d %le %le %le\n", &ii, &ii,&ii,&(x_bp[i]), &(y_bp[i]), &(z_bp[i]));

	 x_bp[i] = x_bp[i]/cl + CMx_c;
	 y_bp[i] = y_bp[i]/cl + CMy_c;
	 z_bp[i] = z_bp[i]/cl+ CMz_c;

	 ibm->x_bp[i] = x_bp[i]*L_dim;
	 ibm->y_bp[i] = y_bp[i]*L_dim;
	 ibm->z_bp[i] = z_bp[i]*L_dim;

	 ibm->x_bp0[i] = x_bp[i]*L_dim;
	 ibm->y_bp0[i] = y_bp[i]*L_dim;
	 ibm->z_bp0[i] = z_bp[i]*L_dim;

	 ibm->x_bp_o[i] = x_bp[i]*L_dim;
	 ibm->y_bp_o[i] = y_bp[i]*L_dim;
	 ibm->z_bp_o[i] = z_bp[i]*L_dim;

	 ibm->u[i].x = 0.;
	 ibm->u[i].y = 0.;
	 ibm->u[i].z = 0.;

	 ibm->uold[i].x = 0.;
	 ibm->uold[i].y = 0.;
	 ibm->uold[i].z = 0.;
      
	 ibm->urm1[i].x = 0.;
	 ibm->urm1[i].y = 0.;
	 ibm->urm1[i].z = 0.;      
	
	} // i(n_v) or lines in nlist.
      /*------- Print out a few lines */
      PetscPrintf(PETSC_COMM_WORLD," nlist READ Complete \n");
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "co-ordinates(xyz) for node n=%d: %le %le %le\n",i, x_bp[i], y_bp[i], z_bp[i]);
      i=50;
      PetscPrintf(PETSC_COMM_WORLD, "co-ordinates(xyz) for node n=%d: %le %le %le\n",i, x_bp[i], y_bp[i], z_bp[i]);
      i=n_v-1;
      PetscPrintf(PETSC_COMM_WORLD, "co-ordinates(xyz) for node n=%d: %le %le %le\n",i, x_bp[i], y_bp[i], z_bp[i]);
      //-------------------------------

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);   // Broadcast co-ordinate values to all processors.
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

      fclose(fd);
    } // for fd: nlist file open check.

    //Reading elements list
    PetscPrintf(PETSC_COMM_SELF, "READ elist begins\n");

    sprintf(filen,"elist%2.2d" , ibi);
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file");

    if (fd) {
     
      
      fscanf(fd, "%i",&n_elmt);  // read first line of elist file
      PetscPrintf(PETSC_COMM_SELF, "number of element of list %d %d \n",ibi, n_elmt);
      
      ibm->n_elmt = n_elmt;      // set number of elements to n_elmt.
     
      
      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);  // broadcast number of elements to all processors.
      
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);    // Allocate memory for integer lists nv1,nv2 and nv3 each of which are the size of number of elements as nv1, nv2 and nv3 represent the indices of  first,second and third node that make up an element.
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x); // similary to nv1,nv2 and nv3, nfx,nfy and nfz represent the x,y and z components of normals of elements and hence are each real valued lists of length n_elmt.
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      // Added 4/1/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &dA); //Area of each element.
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x); // similar to nf_x,nf_y and nf_z but for tangent.
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);  // what is ns_x?? add to questions.txt
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);
      
      // Added 6/4/06 iman
 
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));  // cent_x,cent_y and cent_z represent geometric centres of each element.
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added
      
      
      i=0;
      fgets(string, 128, fd);//skip one line 
      //    while(fgets(string, 128, fd)) {
//	itr = 0;	

	while (i<n_elmt) {
	i++;

	fscanf(fd, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", &ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&nv1[i-1], &nv2[i-1], &nv3[i-1],&ii);
	nv1[i-1] = nv1[i-1] - 1; nv2[i-1] = nv2[i-1]-1; nv3[i-1] = nv3[i-1] - 1; // Adjusting for the first line of elist00 being the number of elements.

	} // elements loop (or number of lines in elist)
	//  } closing of first while loop
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;  // The nodes(1,2,3) corresponding to each element are stored. 
       /*------- Print out a few lines */
      PetscPrintf(PETSC_COMM_WORLD," elist READ Complete \n");      
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "Elemtn n_elmt = %d is comprised of nodes nv = %d,%d and %d\n", i,nv1[i], nv2[i], nv3[i]);
      i=30;
       PetscPrintf(PETSC_COMM_WORLD, "Elemtn n_elmt = %d is comprised of nodes nv = %d,%d and %d\n", i,nv1[i], nv2[i], nv3[i]);
      i=n_elmt-1;
       PetscPrintf(PETSC_COMM_WORLD, "Elemtn n_elmt = %d is comprised of nodes nv = %d,%d and %d\n", i,nv1[i], nv2[i], nv3[i]); 

      fclose(fd);
    } // file elist open condition.

    x_bp=ibm->x_bp;  y_bp=ibm->y_bp ; z_bp = ibm->z_bp ; // The co-ordinates of every node are taken into local memory.  
    for (i=0; i<n_elmt; i++) {
      // going through each element, the nodes of the element are considered, the x,y,z distances between each consecutive nodes of the element are measured.   
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i]; 
      dx12 = x_bp[n2e] - x_bp[n1e];  // x displacement between node 2 and 1
      dy12 = y_bp[n2e] - y_bp[n1e];  // y """""""""""""""""""""""""""""""""
      dz12 = z_bp[n2e] - z_bp[n1e];  // z """""""""""""""""""""""""""""""""
      
      dx13 = x_bp[n3e] - x_bp[n1e];  // x displacement between node 2 and 1 
      dy13 = y_bp[n3e] - y_bp[n1e];  // y displacement between node 3 and 1
      dz13 = z_bp[n3e] - z_bp[n1e];  // z displacement between node 3 and 1

      
      nf_x[i] = dy12 * dz13 - dz12 * dy13;  // normal component(projected area) along x is calculated (un-normalized)
      nf_y[i] = -dx12 * dz13 + dz12 * dx13; //""""""""""""""""""""""" y """""""""""""""""""""""""""""
      nf_z[i] = dx12 * dy13 - dy12 * dx13;  //""""""""""""""""""""""" z """""""""""""""""""""""""""""
      
      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);  // area of parallelogram created by normals calculated  
   
      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr; // normals are normalized using the area.
      
      // Addedd 4/2/06 iman
      if ((((1.-nf_z[i])<=1e-6 )&((-1.+nf_z[i])<1e-6))|
	  (((nf_z[i]+1.)<=1e-6 )&((-1.-nf_z[i])<1e-6))) {
	ns_x[i] = 1.;     
	ns_y[i] = 0.;     
	ns_z[i] = 0. ;
	
	nt_x[i] = 0.;
	nt_y[i] = 1.;
	nt_z[i] = 0.;
      } else {
	ns_x[i] =  nf_y[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);      
	ns_y[i] = -nf_x[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);     
	ns_z[i] = 0. ;
	
        // tangents  are calculated likewise.
	nt_x[i] = -nf_x[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_y[i] = -nf_y[i]*nf_z[i]/ sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
	nt_z[i] = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i]);
      }
 
      //Added 4/1/06 iman
      dA[i] = dr/2.;  // Area of element is also calculated.
      
      // Added 6/4/06 iman
      // Calc the center of the element
      ibm->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e])/3.;
      ibm->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e])/3.;
      ibm->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e])/3.;	

 
    }
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z; // Normals stored in ibm datastructure (global memory)
     PetscPrintf(PETSC_COMM_WORLD,"Element Normals(nf_x,y,z) calculated\n");   
    //Added 4/1/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;  
    PetscPrintf(PETSC_COMM_WORLD,"Element tangents(nt_x,y,z) and ns_x,y,z calculated\n");     
    PetscPrintf(PETSC_COMM_WORLD,"Element centers(centt_x,y,z) and areas(dA) calculated\n");

    // The calculated normals,tangents,centers,areas and ns are all broadcast to all processors.
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/1/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 4/2/06 iman
    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    // Added 6/4/06 iman
    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/***********************************************************************/
// write the surface file.
    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",0,ibi);
    fd = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    PetscFPrintf(PETSC_COMM_WORLD, fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->z_bp[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nf_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->nt_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_x[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_y[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%e\n", ibm->ns_z[i]);
    }
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
    }
    fclose(fd);
    }// if(!rank) i.e for master processor
  else if (rank) {
    /* All other processors just allocate memory for all the data being read in by master processor */
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;
    
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);  
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));
    
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
    /* All the velocities are initialized to zero*/
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;

      ibm->uold[i].x = 0.;
      ibm->uold[i].y = 0.;
      ibm->uold[i].z = 0.;
         
      ibm->urm1[i].x = 0.;
      ibm->urm1[i].y = 0.;
      ibm->urm1[i].z = 0.;      
    }
        
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);    
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->x_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp_o, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    //Added 4/1/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &dA);

    //Added 4/2/06 iman
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nt_z);

    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
    
    // Added 4/2/06 iman
    ibm->dA = dA;
    ibm->nt_x = nt_x; ibm->nt_y = nt_y;  ibm->nt_z = nt_z;
    ibm->ns_x = ns_x; ibm->ns_y = ns_y;  ibm->ns_z = ns_z;    

    // Added 6/4/06
    //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
    PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    //Added 4/2/06 iman
    MPI_Bcast(ibm->dA, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nt_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nt_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->ns_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->ns_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->cent_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->cent_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }

  return(0);
}
//----------------------------------------------------------------------

PetscErrorCode ibm_scale(IBMNodes *ibm, PetscReal scale[])
{
  PetscInt   i;
  PetscReal  xmin,xmax,ymin,ymax,zmin,zmax, cl;

  ymin=1.e+5;
  ymax=-1.e+5;
  xmin=ymin;
  xmax=ymax;
  zmin=ymin;
  zmax=ymax;

  for (i=0; i<ibm->n_v; i++) {
    xmin=PetscMin(xmin,(ibm->x_bp0[i]));
    xmax=PetscMax(xmax,(ibm->x_bp0[i]));
    ymin=PetscMin(ymin,(ibm->y_bp0[i]));
    ymax=PetscMax(ymax,(ibm->y_bp0[i]));
    zmin=PetscMin(zmin,(ibm->z_bp0[i]));
    zmax=PetscMax(zmax,(ibm->z_bp0[i]));
  }

  PetscBool flg;
  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, &flg);     

  if (flg)
    cl=1.;
  else
    cl=(zmax-zmin);
  if (wing)
    cl=(ymax-ymin);

  zmin /= cl;
  zmax /= cl;
  ymin /= cl;
  ymax /= cl;
  xmin /= cl;
  xmax /= cl;

  for (i=0; i<ibm->n_v; i++) {
    ibm->x_bp0[i]/=cl; 
    ibm->y_bp0[i]/=cl; 
    ibm->z_bp0[i]/=cl; 
    
    ibm->x_bp[i]/=cl; 
    ibm->y_bp[i]/=cl; 
    ibm->z_bp[i]/=cl; 
    
    ibm->x_bp_o[i]/=cl; 
    ibm->y_bp_o[i]/=cl; 
    ibm->z_bp_o[i]/=cl;

    if (!wing) {
      ibm->x_bp0[i]-=0.5*(xmin+xmax);
      ibm->y_bp0[i]-=0.5*(ymin+ymax);
      ibm->z_bp0[i]-=zmin;

      ibm->x_bp[i]-=0.5*(xmin+xmax);
      ibm->y_bp[i]-=0.5*(ymin+ymax);
      ibm->z_bp[i]-=zmin;

      ibm->x_bp_o[i]-=0.5*(xmin+xmax);
      ibm->y_bp_o[i]-=0.5*(ymin+ymax);
      ibm->z_bp_o[i]-=zmin;
    }
  }

  cl=zmin;           zmin -= cl;  zmax -= cl;
  cl=0.5*(ymin+ymax);ymin -= cl;  ymax -= cl;
  cl=0.5*(xmin+xmax);xmin -= cl;  xmax -= cl;

  scale[0]=zmin;
  scale[1]=zmax;
  scale[2]=ymin;
  scale[3]=ymax;
  scale[4]=xmin;
  scale[5]=xmax;
  scale[6]=cl;
  return(0);
}

PetscErrorCode calc_ibm_normal(IBMNodes *ibm)
{
  PetscInt   n1e, n2e, n3e, i;
  PetscReal  dx12, dy12, dz12, dx13, dy13, dz13, dr;
  
  for (i=0; i<ibm->n_elmt; i++) {
    //PetscPrintf(PETSC_COMM_WORLD, "cop nf %d !\n",i);       
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e = ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e];
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e];
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e];
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;
    
    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + 
	      ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;

    if ((((1.-ibm->nf_z[i])<=1e-6 )&((-1.+ibm->nf_z[i])<1e-6))|
	(((ibm->nf_z[i]+1.)<=1e-6 )&((-1.-ibm->nf_z[i])<1e-6))) {
      ibm->ns_x[i] = 1.;     
      ibm->ns_y[i] = 0.;     
      ibm->ns_z[i] = 0. ;
      
      ibm->nt_x[i] = 0.;
      ibm->nt_y[i] = 1.;
      ibm->nt_z[i] = 0.;
    } else {
      ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + 
					 ibm->nf_y[i]*ibm->nf_y[i]);      
      ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + 
					 ibm->nf_y[i]*ibm->nf_y[i]);     
      ibm->ns_z[i] = 0. ;
      
      ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + 
						      ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + 
						      ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
    }
    
    //Added 4/1/06 iman
    ibm->dA[i] = dr/2.; 
    
    // Added 6/4/06 iman
    // Calc the center of the element
    ibm->cent_x[i]= (ibm->x_bp[n1e]+ibm->x_bp[n2e]+ibm->x_bp[n3e])/3.;
    ibm->cent_y[i]= (ibm->y_bp[n1e]+ibm->y_bp[n2e]+ibm->y_bp[n3e])/3.;
    ibm->cent_z[i]= (ibm->z_bp[n1e]+ibm->z_bp[n2e]+ibm->z_bp[n3e])/3.;
  }
  return(0);
}

PetscErrorCode calc_ibm_velocity(IBMNodes *ibm, PetscReal delti)
{
  PetscReal  v_max=0.;
  PetscInt   i, i_vmax=0;
  for (i=0; i<ibm->n_v; i++) {
    ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / delti;
    ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / delti;
    ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / delti;
    
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
  PetscPrintf(PETSC_COMM_WORLD, "MAX fish Velocity:%le %le %le %le\n", v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax]);
  
  return(0);
}


#define CROSS(dest, v1, v2) \
	dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
	dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
	dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define SUB(dest, v1, v2) \
	dest[0] = v1[0] - v2[0]; \
	dest[1] = v1[1] - v2[1]; \
	dest[2] = v1[2] - v2[2];

PetscErrorCode calc_ibm_volumeFlux(IBMNodes *ibm, PetscReal delti, PetscReal *VolumeFlux)
{
  PetscInt   i,ii, n_elmt=ibm->n_elmt;
  PetscInt n1e, n2e, n3e;
  PetscReal  Vol=0.,p1[3],p2[3],p3[3],p4[3],p5[3],p6[3],sign,nf[3];
  PetscReal  edge1[3],edge2[3],edge3[3],edge2c3[3],volTH,cent[3],cent_o[3],dir[3];

  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];

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
    
    // calculate volume flux
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
  *VolumeFlux = Vol;
  PetscPrintf(PETSC_COMM_WORLD, "Volume Flux %e\n", Vol);

  return(0);
}
