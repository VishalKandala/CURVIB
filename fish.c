#include "variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c;
extern PetscInt ti;
extern PetscInt tiout,eel,pizza;
extern PetscReal St_exp,wavelength;
PetscReal ampl,V_slip,omega, kwave,alpha;
PetscReal T_period_fish;
PetscReal c_0,c_1,c_2;
PetscInt  N_period_fish,pizza;

PetscErrorCode ibm_read_fish(IBMNodes *ibm,PetscReal cl)
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
 
  //PetscReal     cl=2.*15.;
  //double xt;
  char string[128];

  // Temp Sol. Change based on file !!!
/*   n_v= 3397; */
/*   n_elmt = 6676; */

//  cl = 1.;
  n_v= 1694; //1682;//1618;
  n_elmt = 3384;//3312;//3152;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ nlist\n");
/*     fd = fopen("f_nlist_6676.txt", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file nlist"); */
    fd = fopen("nlist_Mack.txt", "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file nlist");
    //   n_v =0;

    if (fd) {

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

      //for (i=0; i<n_v; i++) {
      i=-1;
      fgets(string, 128, fd);// miss 2 lines
      while(fgets(string, 128, fd)) {
	//fgets(string, 128, fd);// miss 2 lines
	//fscanf(fd, "\n");
	//fscanf(fd, "%s %s %s %s %s %s %s\n", ss, ss, ss, ss, ss, ss, ss);
	itr = 0;
	
	// read the data on 20 lines
	while (i+1<n_v && itr<20) {
	itr++;	i++;
	//fscanf(fd, "%d %3.3e %3.4e %3.3e %1.2e %1.2e %1.2e\n", &ii, &(x_bp[i]), &(y_bp[i]), &(z_bp[i]), &tt, &tt, &tt);
      	//fscanf(fd, "%d %le %le %le %le %le %le\n", &ii, &(x_bp[i]), &(y_bp[i]), &(z_bp[i]), &tt, &tt, &tt);
	fscanf(fd, "%d %le %le %le\n",&ii,&(x_bp[i]),&(y_bp[i]),&(z_bp[i]));

	x_bp[i] = x_bp[i]/cl + CMx_c;
	y_bp[i] = y_bp[i]/cl + CMy_c;
	z_bp[i] = z_bp[i]/cl + CMz_c;
	
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
      }
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %d %le %le %le\n",i, x_bp[i], y_bp[i], z_bp[i]);
      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %d %le %le %le\n",i, x_bp[i], y_bp[i], z_bp[i]);
      i=100;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %d %le %le %le\n",i, x_bp[i], y_bp[i], z_bp[i]);

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

      fclose(fd);
    }

    PetscPrintf(PETSC_COMM_SELF, "READ elist\n");
/*     fd = fopen("f_elist_6676.txt", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file elist"); */
    fd = fopen("elist_Mack.txt", "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file elist");
						  // n_v =0;
    if (fd) {
      //for (i=0; i<n_elmt; i++) {
      i=0;
      fgets(string, 128, fd);// miss 2 lines
      while(fgets(string, 128, fd)) {
	itr = 0;	
	// read the data on 20 lines
	while (i<n_elmt && itr<20) {
	itr++;	i++;

	fscanf(fd, "%d %d %d %d %d %d %d %d %d %d\n", &ii,&ii, &ii,&ii,&ii,&ii, &nv1[i-1], &nv2[i-1], &nv3[i-1],&ii);
	nv1[i-1] = nv1[i-1] - 1; nv2[i-1] = nv2[i-1]-1; nv3[i-1] = nv3[i-1] - 1;
	//PetscPrintf(PETSC_COMM_WORLD, "ii %d nv %d %d %d\n",i, nv1[i-1], nv2[i-1], nv3[i-1]);
	}
      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

      fclose(fd);
    }

    x_bp=ibm->x_bp;  y_bp=ibm->y_bp ; z_bp = ibm->z_bp ;
    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n");      
    for (i=0; i<n_elmt; i++) {
      //PetscPrintf(PETSC_COMM_WORLD, "cop nf %d !\n",i);       
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

 /*      ibm->nf_x[i] = nf_x[i]; ibm->nf_y[i] = nf_y[i];  ibm->nf_z[i] = nf_z[i]; */
      
/*       //Added 4/1/06 iman */
/*       ibm->dA[i] = dA[i]; */
/*       ibm->nt_x[i] = nt_x[i]; ibm->nt_y[i] = nt_y[i];  ibm->nt_z[i] = nt_z[i]; */
/*       ibm->ns_x[i] = ns_x[i]; ibm->ns_y[i] = ns_y[i];  ibm->ns_z[i] = ns_z[i];   */
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

 /*    PetscFree(dA); */
/*     PetscFree(nf_x);PetscFree(nf_y);PetscFree(nf_z); */
/*     PetscFree(nt_x);PetscFree(nt_y);PetscFree(nt_z); */
/*     PetscFree(ns_x);PetscFree(ns_y);PetscFree(ns_z); */
/*     PetscFree(nv1);PetscFree(nv2);PetscFree(nv3); */
/*     PetscFree(x_bp);PetscFree(y_bp);PetscFree(z_bp); */

    char filen[80];
   
    sprintf(filen, "surface%3.3d.stl",ti);
    fd = fopen(filen, "w");
    PetscPrintf(PETSC_COMM_WORLD, "fish nf!\n"); 

    PetscFPrintf(PETSC_COMM_WORLD,fd, "%s %s\n", "solid", "fish.stl");
    PetscPrintf(PETSC_COMM_WORLD, "fish nf!\n"); 

    for (i=0; i<ibm->n_elmt; i++) {
      //PetscPrintf(PETSC_COMM_WORLD, "reading i,ii,j,  %d %d %d %s \n", i,ii,j,ss);
      PetscFPrintf(PETSC_COMM_WORLD,fd, "%s %s %le %le %le\n", "facet", "normal", ibm->nf_x[i], ibm->nf_y[i], ibm->nf_z[i]);
      PetscFPrintf(PETSC_COMM_WORLD,fd, "  %s %s\n","outer", "loop");

      n1e=ibm->nv1[i]; n2e=ibm->nv2[i]; n3e=ibm->nv3[i];

      //PetscPrintf(PETSC_COMM_WORLD, "reading i,ii,j,  %d %d %d %s \n", i,ii,j,ss);
      PetscFPrintf(PETSC_COMM_WORLD,fd, "    %s %le %le %le\n", "vertex", ibm->x_bp[n1e], ibm->y_bp[n1e], ibm->z_bp[n1e]);
      PetscFPrintf(PETSC_COMM_WORLD,fd, "    %s %le %le %le\n", "vertex", ibm->x_bp[n2e], ibm->y_bp[n2e], ibm->z_bp[n2e]);
      PetscFPrintf(PETSC_COMM_WORLD,fd, "    %s %le %le %le\n", "vertex", ibm->x_bp[n3e], ibm->y_bp[n3e], ibm->z_bp[n3e]);
 
      PetscFPrintf(PETSC_COMM_WORLD,fd, "  %s\n", "endloop");
      PetscFPrintf(PETSC_COMM_WORLD,fd, "%s\n", "endfacet");
     //PetscPrintf(PETSC_COMM_WORLD, "reading i,ii,j,  %d %d %d %s \n", i,ii,j,ss);
    }
    PetscFPrintf(PETSC_COMM_WORLD,fd, "%s %s\n", "endsolid", "fish.stl");
    fclose(fd);

    FILE *f;
    sprintf(filen, "surface_nf%3.3d.dat",ti);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL)\n", n_v, n_elmt);
    //PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    //    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt);
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<n_v; i++) {	
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
    }
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]); */
/*     } */
/*     for (i=0; i<n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]); */
/*     } */
    for (i=0; i<n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
    }
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

PetscErrorCode fish_init(PetscReal *delti)
{
  PetscReal  pi = 3.141592653589793, St;
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-N_period_fish", &N_period_fish, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "N_period_fish is : %d \n",N_period_fish);
  V_slip=.7;

  alpha = 2.18;
  if (eel) 
    ampl=0.1;//0.089;
  else
    ampl=0.1;//0.09 // non-dim ampl a/L
  c_0=0.02;
  c_1=-0.08;
  c_2=0.16;
/*   c_0=0.02; */
/*   c_2=(0.3*(ampl-c_0)+c_0)/0.291;  */
/*   c_1=ampl - c_0 - c_2 ;//  ! Ampl. is defined at the end of tail */

  St=St_exp/(2.*ampl);  // non-dim  frequency St=fL/U
  omega=2*pi*St;//2.*pi; //non-dimensionalized w
  //kwave=2.*pi*St*V_slip; //non-dimensionalized k
  //wavelength = 0.95;
  kwave= 2*pi/wavelength;//8.307767;
  V_slip= kwave*0.5/pi/St;

  T_period_fish=2.*pi/omega;
  *delti = T_period_fish/N_period_fish;  

  PetscPrintf(PETSC_COMM_WORLD, "fish init: St_exp %le Amplitude %le coeff a0 a1 a2 %le %le %le\n",St_exp,ampl,c_0,c_1,c_2);
  PetscPrintf(PETSC_COMM_WORLD, "fish init: dt %le w %le f %le k %le lambda %le T %le N %d  V_slip %le\n",*delti,omega,St,kwave,wavelength,T_period_fish,N_period_fish,V_slip);
 
  return(0);
}

PetscErrorCode fish_swim(IBMNodes *ibm, PetscReal time
			,PetscReal delti)
{

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-pizza", &pizza, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "N_period_fish is : %d \n",N_period_fish);
 
  PetscReal  h0,h1,h2,h3,h4;

  PetscInt   i, n_v=ibm->n_v, n_elmt=ibm->n_elmt;
  PetscReal  x0,x1,y0,y1,z0,z1;
  PetscReal  az0,dz;

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

  PetscPrintf(PETSC_COMM_WORLD, "MAX fish: %d %le %le %le %le %le %le\n",ti,z1,z0,y1,y0,x1,x0);
  if(pizza){
    kwave=0.0;
    h0=exp(-1.0)*ampl*sin(-omega*time);
    h1=exp(-0.75)*ampl*sin(kwave*0.25-omega*time);
    h2=exp(-0.5)*ampl*sin(kwave*0.5-omega*time);
    h3=exp(-0.25)*ampl*sin(kwave*0.75-omega*time);
    h4=ampl*sin(kwave*1.0-omega*time);
    PetscPrintf(PETSC_COMM_WORLD, " t %le  h0 %le  h1 %le h2 %le h3 %le h4 %le\n",time,h0,h1,h2,h3,h4);
  }
  /*   oscillation for Mackerel in x-dir */
  for (i=0; i<n_v; i++) {
    dz=(ibm->z_bp0[i]-z0);//+(13.3086/15.0);//for the tail simulation only
    if(pizza){
      if(dz <=0.25)        ibm->x_bp[i]=ibm->x_bp0[i]+h1+(dz-0.25)*(h1-h0)/(0.25);
      else  if(dz <=0.50)  ibm->x_bp[i]=ibm->x_bp0[i]+h2+(dz-0.5)*(h2-h1)/(0.25);
      else  if(dz <=0.75)  ibm->x_bp[i]=ibm->x_bp0[i]+h3+(dz-0.75)*(h3-h2)/(0.25);
      else  if(dz <=1.0)   ibm->x_bp[i]=ibm->x_bp0[i]+h4+(dz-1.0)*(h4-h3)/(0.25);
      else                 ibm->x_bp[i]=ibm->x_bp0[i]+0.0;
    }else{
      if (eel) // Anguiliform
	az0=ampl*exp(alpha*(dz-1.));
      else    // Carangiform:
	az0=c_0+c_1*dz+c_2*dz*dz;

      ibm->x_bp[i]=ibm->x_bp0[i] +

	az0*sin(kwave*dz-omega*time);
    
    }
  }
  
/*    Calculate the new normal & velcity */
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  
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
  }
  
  PetscReal  v_max=0.;
  PetscInt   i_vmax=0;
  for (i=0; i<n_v; i++) {
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
  PetscPrintf(PETSC_COMM_WORLD, "MAX fish Velocity: %d %le %le %le %le\n",ti, v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax]);
  
/*   PetscInt rank; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   if (!rank) { */
/*     if (ti == (ti/tiout)*tiout) { */
/*       FILE *f; */
/*       char filen[80]; */
/*       sprintf(filen, "surface%3.3d.dat",ti); */
/*       f = fopen(filen, "w"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt); */
/*       for (i=0; i<n_v; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]); */
/*       } */
/*       for (i=0; i<n_v; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]); */
/*       } */
/*       for (i=0; i<n_v; i++) {	 */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*       } */
/*       fclose(f); */

/* /\*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n"); *\/ */
/* /\*       PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt); *\/ */
/* /\*       for (i=0; i<n_v; i++) { *\/ */
/* /\* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]); *\/ */
/* /\*       } *\/ */
/* /\*       for (i=0; i<n_elmt; i++) { *\/ */
/* /\* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); *\/ */
/* /\*       } *\/ */
/* /\*       fclose(f); *\/ */

/*     } */
/*   } */

  return(0);
}
