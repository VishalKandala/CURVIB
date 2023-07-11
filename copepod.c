#include "variables.h"
#include "cop_variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c;
extern PetscInt ti;
extern PetscInt tiout, NumberOfBodies;
extern PetscReal L_dim;

PetscErrorCode ibm_read_cop(IBMNodes *ibm)
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
  PetscReal     cl=168.;
  //char   ss[20];
  //double xt;
  char string[128];

  // Temp Sol. Change based on file !!!
/*   n_v= 4682; */
/*   n_elmt = 9284; */
  n_v= 4554;
  n_elmt = 9028;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ nlist, %le\n", L_dim);
/*     fd = fopen("nlist_0210.txt", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file nlist"); */
    fd = fopen("nlist_cop01.txt", "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file nlist");
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
	fscanf(fd, "%d %le %le %le %le %le %le\n", &ii, &(x_bp[i]), &(y_bp[i]), &(z_bp[i]), &tt, &tt, &tt);
	//PetscPrintf(PETSC_COMM_WORLD, "i,ii, xyz_bp %d %d %e %e %e\n",i,ii, x_bp[i], y_bp[i], z_bp[i]);
	x_bp[i] = x_bp[i]/cl + CMx_c;//0.25 ;// 24.;
	y_bp[i] = y_bp[i]/cl + CMy_c;//2.;//8.;//6.;//2.   ;// 24.;
	z_bp[i] = z_bp[i]/cl + CMz_c;//2.;//8.;//15.;//2.   ;// 24.;
	
/* 	ibm->x_bp[i] = x_bp[i]; */
/* 	ibm->y_bp[i] = y_bp[i]; */
/* 	ibm->z_bp[i] = z_bp[i]; */

	ibm->x_bp[i] = y_bp[i]*L_dim;
	ibm->y_bp[i] = z_bp[i]*L_dim;
	ibm->z_bp[i] = x_bp[i]*L_dim;

	ibm->x_bp0[i] = x_bp[i]*L_dim;
	ibm->y_bp0[i] = y_bp[i]*L_dim;
	ibm->z_bp0[i] = z_bp[i]*L_dim;

	ibm->x_bp_o[i] = y_bp[i]*L_dim;
	ibm->y_bp_o[i] = z_bp[i]*L_dim;
	ibm->z_bp_o[i] = x_bp[i]*L_dim;

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

      fclose(fd);
    }

    PetscPrintf(PETSC_COMM_SELF, "READ elist\n");
/*     fd = fopen("elist_0210.txt", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file elist"); */
    fd = fopen("elist_cop01.txt", "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file elist");
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

    PetscPrintf(PETSC_COMM_WORLD, "cop nf!\n"); 
 /*    PetscFree(dA); */
/*     PetscFree(nf_x);PetscFree(nf_y);PetscFree(nf_z); */
/*     PetscFree(nt_x);PetscFree(nt_y);PetscFree(nt_z); */
/*     PetscFree(ns_x);PetscFree(ns_y);PetscFree(ns_z); */
/*     PetscFree(nv1);PetscFree(nv2);PetscFree(nv3); */
/*     PetscFree(x_bp);PetscFree(y_bp);PetscFree(z_bp); */

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


PetscErrorCode cop_init(IBMNodes *ibm)
{
  nib_tail=1;
  nib_antene=2;
  nib_leg12=3;
  nib_leg34=4;
  nib_leg56=5;
  nib_leg78=6;
/*   x_bp=ibm->x_bp0; y_bp=ibm->y_bp0; z_bp=ibm->z_bp0; */

  //L_dim = 1.;//0.005;

/*   n1_min=1    ; */
/*   n1_max=203; */

/*   n2_min=204; */
/*   n2_max=406; */

/*   n3_min=407; */
/*   n3_max=609; */

/*   n4_min=610; */
/*   n4_max=812; */

/*   n5_min=813; */
/*   n5_max=1015; */

/*   n6_min=1016; */
/*   n6_max=1218; */

/*   n7_min=1219; */
/*   n7_max=1565; */

/*   n8_min=1566; */
/*   n8_max=1912; */

/*   n9_min=1913; */
/*   n9_max=2039; */

/*   n10_min=2040; */
/*   n10_max=2166; */

/*   n11_min=2167; */
/*   n11_max=2342; */

/*   n12_min=2343; */
/*   n12_max=2518; */

/*   n_antene_mn=2519; */
/*   n_antene_mx=3882; */

/*   n_antene_0=3817;  // z_max */
/*   n_antene_1=3841;  // z_min */
/*   n_antene_2=3818;  // x_min  */
/*   n_antene_3=3883;  // point of rot. */

/*   n_tail_0  =5508; */
    
/*   n_bot_0=123; */
/*   n_bot_1=2488; //!!!!this is not correct */
/*   n_bot_2=2488; */

/*   n1_min=1-1   ;  */
/*   n1_max=94-1; */

/*   n2_min=95-1; */
/*   n2_max=188-1; */

/*   n3_min=189-1; */
/*   n3_max=282-1; */

/*   n4_min=283-1; */
/*   n4_max=376-1; */
    
/*   n5_min=377-1; */
/*   n5_max=470-1; */

/*   n6_min=471-1; */
/*   n6_max=564-1; */

/*   n7_min=565-1; */
/*   n7_max=719-1; */

/*   n8_min=720-1; */
/*   n8_max=874-1; */

/*   n9_min=875-1; */
/*   n9_max=975-1; */

/*   n10_min=976-1; */
/*   n10_max=1076-1; */

/*   n11_min=1077-1; */
/*   n11_max=1209-1; */

/*   n12_min=1210-1; */
/*   n12_max=1342-1; */

/*   n_antene_mn=1343-1; */
/*   n_antene_mx=2088-1; */

/*   n_antene_0=2048-1  ;//! z_max */
/*   n_antene_1=2065-1  ;//! z_min */
/*   n_antene_2=2053-1  ;//! x_min  */
/*   n_antene_3=1672-1  ;//! point of rot. */

/*   n_tail_0  =2856-1; */
    
/*   n_bot_0=4011-1; */
/*   n_bot_1=4300-1; */
/*   n_bot_2=4300-1; */

/*   x_tail_0=x_bp[n_tail_0]; */

/*   z_antene_0=z_bp[n_antene_0]; */
/*   z_antene_1=z_bp[n_antene_1]; */
/*   x_antene_2=x_bp[n_antene_2]; */

/* /\*   write(*,*) 'z_antene_0,z_antene_1,x_antene_2',z_antene_0,z_antene_1,x_antene_2 *\/ */
/* /\*   write(3,*) 'z_antene_0,z_antene_1,x_antene_2',z_antene_0,z_antene_1,x_antene_2 *\/ */
  
/*   z_bot_0=z_bp[n_bot_0]; */
/*   z_bot_1=z_bp[n_bot_1]; */
/*   z_bot_2=z_bp[n_bot_2]; */

/* /\*   write(*,*) 'z_bot_0,z_bot_1,z_bot_2',z_bot_0,z_bot_1,z_bot_2 *\/ */
/* /\*   write(3,*) 'z_bot_0,z_bot_1,z_bot_2',z_bot_0,z_bot_1,z_bot_2 *\/ */

/*   n_tail=0; */
/*   n_antene=0; */
/*   n_leg1=0;n_leg2=0;n_leg3=0;n_leg4=0; */
/*   n_leg5=0;n_leg6=0;n_leg7=0;n_leg8=0; */
/*   n_leg9=0;n_leg10=0;n_leg11=0;n_leg12=0; */


  //def. of points of rotation
  n_tail_0=7-1;
  sx_tail=ibm[nib_tail].x_bp0[n_tail_0];
  sy_tail=ibm[nib_tail].y_bp0[n_tail_0];
  sz_tail=ibm[nib_tail].z_bp0[n_tail_0];

/*   write(*,*) 'n_tail_0,sx_tail,sy_tail,sz_tail',n_tail_0,sx_tail,sy_tail,sz_tail */
/*   write(3,*) 'n_tail_0,sx_tail,sy_tail,sz_tail',n_tail_0,sx_tail,sy_tail,sz_tail */

  PetscPrintf(PETSC_COMM_WORLD, "sx tail %le %le %le\n",sx_tail,sy_tail,sz_tail);

  n_antene_3=353-1;//339-1;
/*   sx_antene=ibm[nib_antene].x_bp0[n_antene_3]; */
/*   sy_antene=ibm[nib_antene].y_bp0[n_antene_3]; */
/*   sz_antene=ibm[nib_antene].z_bp0[n_antene_3]; */
  sx_antene=ibm[nib_antene].x_bp[n_antene_3];
  sy_antene=ibm[nib_antene].y_bp[n_antene_3];
  sz_antene=ibm[nib_antene].z_bp[n_antene_3];

  PetscPrintf(PETSC_COMM_WORLD, "sx antene %le %le %le\n",sx_antene,sy_antene,sz_antene);

  n_leg_1_rot=188-1;
  n_leg_3_rot=188-1;
  n_leg_5_rot=94-1;
  n_leg_7_rot=154-1;
/*     //! for regime 2 */
/*   n_leg_9_rot=938-1; */
/*   n_leg_10_rot=1036-1; */
/*   n_leg_11_rot=1162-1; */
/*   n_leg_12_rot=1295-1; */
    
  sx_l1 = ibm[nib_leg12].x_bp0[n_leg_1_rot];   sy_l1= ibm[nib_leg12].y_bp0[n_leg_1_rot];    sz_l1= ibm[nib_leg12].z_bp0[n_leg_1_rot];
  sx_l3 = ibm[nib_leg34].x_bp0[n_leg_3_rot];   sy_l3= ibm[nib_leg34].y_bp0[n_leg_3_rot];    sz_l3= ibm[nib_leg34].z_bp0[n_leg_3_rot];
  sx_l5 = ibm[nib_leg56].x_bp0[n_leg_5_rot];   sy_l5= ibm[nib_leg56].y_bp0[n_leg_5_rot];    sz_l5= ibm[nib_leg56].z_bp0[n_leg_5_rot];
  sx_l7 = ibm[nib_leg78].x_bp0[n_leg_7_rot];   sy_l7= ibm[nib_leg78].y_bp0[n_leg_7_rot];    sz_l7= ibm[nib_leg78].z_bp0[n_leg_7_rot];
/*   sx_l9 = x_bp[n_leg_9_rot];   sy_l9=y_bp[n_leg_9_rot];    sz_l9=z_bp[n_leg_9_rot]; */
/*   sx_l10= x_bp[n_leg_10_rot];  sy_l10=y_bp[n_leg_10_rot];  sz_l10=z_bp[n_leg_10_rot]; */
/*   sx_l11= x_bp[n_leg_11_rot];  sy_l11=y_bp[n_leg_11_rot];  sz_l11=z_bp[n_leg_11_rot]; */
/*   sx_l12= x_bp[n_leg_12_rot];  sy_l12=y_bp[n_leg_12_rot];  sz_l12=z_bp[n_leg_12_rot]; */

  PetscPrintf(PETSC_COMM_WORLD, "sx legs 12 %le %le %le\n",sx_l1,sy_l1,sz_l1);
  PetscPrintf(PETSC_COMM_WORLD, "sx legs 34 %le %le %le\n",sx_l3,sy_l3,sz_l3);
  PetscPrintf(PETSC_COMM_WORLD, "sx legs 56 %le %le %le\n",sx_l5,sy_l5,sz_l5);
  PetscPrintf(PETSC_COMM_WORLD, "sx legs 78 %le %le %le\n",sx_l7,sy_l7,sz_l7);

/*   write(*,*) 'antenna:sx_antene,sy_antene,sz_antene',sx_antene,sy_antene,sz_antene */
/*   write(3,*) 'antenna:sx_antene,sy_antene,sz_antene',sx_antene,sy_antene,sz_antene */

/*   !stop */

/*   n_leg_1_rot=46; */
/*   n_leg_2_rot=219; */
/*   n_leg_3_rot=422; */
/*   n_leg_4_rot=625; */
/*   n_leg_5_rot=828; */
/*   n_leg_6_rot=1031; */
/*   n_leg_7_rot=1225; */
/*   n_leg_8_rot=1572; */
/*   n_leg_9_rot=1913; */
/*   n_leg_10_rot=2087; */
/*   n_leg_11_rot=2174; */
/*   n_leg_12_rot=2350; */

/*   n_leg_1_rot=63-1; */
/*   n_leg_2_rot=157-1; */
/*   n_leg_3_rot=251-1; */
/*   n_leg_4_rot=345-1; */
/*   n_leg_5_rot=439-1; */
/*   n_leg_6_rot=533-1; */
/*   n_leg_7_rot=569-1; */
/*   n_leg_8_rot=724-1; */
/*     //! for regime 2 */
/*   n_leg_9_rot=938-1; */
/*   n_leg_10_rot=1036-1; */
/*   n_leg_11_rot=1162-1; */
/*   n_leg_12_rot=1295-1; */
    
/*   sx_l1 = x_bp[n_leg_1_rot];   sy_l1=y_bp[n_leg_1_rot];    sz_l1=z_bp[n_leg_1_rot]; */
/*   sx_l2 = x_bp[n_leg_2_rot];   sy_l2=y_bp[n_leg_2_rot];    sz_l2=z_bp[n_leg_2_rot]; */
/*   sx_l3 = x_bp[n_leg_3_rot];   sy_l3=y_bp[n_leg_3_rot];    sz_l3=z_bp[n_leg_3_rot]; */
/*   sx_l4 = x_bp[n_leg_4_rot];   sy_l4=y_bp[n_leg_4_rot];    sz_l4=z_bp[n_leg_4_rot]; */
/*   sx_l5 = x_bp[n_leg_5_rot];   sy_l5=y_bp[n_leg_5_rot];    sz_l5=z_bp[n_leg_5_rot]; */
/*   sx_l6 = x_bp[n_leg_6_rot];   sy_l6=y_bp[n_leg_6_rot];    sz_l6=z_bp[n_leg_6_rot]; */
/*   sx_l7 = x_bp[n_leg_7_rot];   sy_l7=y_bp[n_leg_7_rot];    sz_l7=z_bp[n_leg_7_rot]; */
/*   sx_l8 = x_bp[n_leg_8_rot];   sy_l8=y_bp[n_leg_8_rot];    sz_l8=z_bp[n_leg_8_rot]; */
/*   sx_l9 = x_bp[n_leg_9_rot];   sy_l9=y_bp[n_leg_9_rot];    sz_l9=z_bp[n_leg_9_rot]; */
/*   sx_l10= x_bp[n_leg_10_rot];  sy_l10=y_bp[n_leg_10_rot];  sz_l10=z_bp[n_leg_10_rot]; */
/*   sx_l11= x_bp[n_leg_11_rot];  sy_l11=y_bp[n_leg_11_rot];  sz_l11=z_bp[n_leg_11_rot]; */
/*   sx_l12= x_bp[n_leg_12_rot];  sy_l12=y_bp[n_leg_12_rot];  sz_l12=z_bp[n_leg_12_rot]; */


/*   for (n=0; n<n_v; n++) { */
/*     //!write(3,*) 'n,xf[n],yf[n],zf[n]',n,xf[n],yf[n],zf[n] */
/*     if(x_bp[n]>=x_tail_0) { */
/*       n_tail=n_tail+1; */
/*       //nv_tail(n_tail)=n; */
/*     } */

/*     if(n>=n_antene_mn && n<=n_antene_mx) { */
/*       if(n>875) { */
/* 	n_antene=n_antene+1; */
/* 	//nv_antene(n_antene)=n; */
/*       } */
/*     } */

/*     if(z_bp[n]<z_bot_0) { */
/*       if(n>=n1_min && n<=n1_max) { */
/* 	n_leg1=n_leg1+1; */
/* 	//nv_leg1(n_leg1)=n; */
/*       } */
/*       if(n>=n2_min && n<=n2_max) { */
/* 	n_leg2=n_leg2+1; */
/* 	//nv_leg2(n_leg2)=n; */
/*       } */
/*       if(n>=n3_min && n<=n3_max) { */
/* 	n_leg3=n_leg3+1; */
/* 	//nv_leg3(n_leg3)=n; */
/*       } */
/*       if(n>=n4_min && n<=n4_max) { */
/* 	n_leg4=n_leg4+1; */
/* 	//nv_leg4(n_leg4)=n; */
/*       } */
/*       if(n>=n5_min && n<=n5_max) { */
/* 	n_leg5=n_leg5+1; */
/* 	//nv_leg5(n_leg5)=n; */
/*       } */
/*       if(n>=n6_min && n<=n6_max) { */
/* 	n_leg6=n_leg6+1; */
/* 	//nv_leg6(n_leg6)=n; */
/*       } */
/*     } */

/*     if(z_bp[n]<z_bot_1) { */
/*       if(n>=n7_min && n<=n7_max) { */
/* 	n_leg7=n_leg7+1; */
/* 	//nv_leg7(n_leg7)=n; */
/*       } */
/*       if(n>=n8_min && n<=n8_max) { */
/* 	n_leg8=n_leg8+1; */
/* 	//nv_leg8(n_leg8)=n; */
/*       } */
/*     } */

/*     //!write(*,*) 'n,zf,z_bot_2',n,zf[n],z_bot_2 */
/*     if(n>=n9_min && n<=n9_max) { */
/*       //!if(n/=5664)  */
/*       n_leg9=n_leg9+1; */
/*       //nv_leg9(n_leg9)=n; */
/*       //} */
/*     } */
/*     if(n>=n10_min && n<=n10_max) { */
/*       n_leg10=n_leg10+1; */
/*       //nv_leg10(n_leg10)=n; */
/*     } */
   
/*     if(z_bp[n]<sz_l11) { */
/*       if(n>=n11_min && n<=n11_max) { */
/* 	n_leg11=n_leg11+1; */
/* 	//nv_leg11(n_leg11)=n; */
/*       } */
/*       if(n>=n12_min && n<=n12_max) { */
/* 	n_leg12=n_leg12+1; */
/* 	//nv_leg12(n_leg12)=n; */
/*       } */
/*     } */

/*   } */

/*   // mem allocation */
/*   PetscMalloc(n_tail*sizeof(PetscInt), &nv_tail); */
/*   PetscMalloc(n_antene*sizeof(PetscInt), &nv_antene); */
/*   PetscMalloc(n_leg1*sizeof(PetscInt), &nv_leg1); */
/*   PetscMalloc(n_leg2*sizeof(PetscInt), &nv_leg2); */
/*   PetscMalloc(n_leg3*sizeof(PetscInt), &nv_leg3); */
/*   PetscMalloc(n_leg4*sizeof(PetscInt), &nv_leg4); */
/*   PetscMalloc(n_leg5*sizeof(PetscInt), &nv_leg5); */
/*   PetscMalloc(n_leg6*sizeof(PetscInt), &nv_leg6); */
/*   PetscMalloc(n_leg7*sizeof(PetscInt), &nv_leg7); */
/*   PetscMalloc(n_leg8*sizeof(PetscInt), &nv_leg8); */
/*   PetscMalloc(n_leg9*sizeof(PetscInt), &nv_leg9); */
/*   PetscMalloc(n_leg10*sizeof(PetscInt), &nv_leg10); */
/*   PetscMalloc(n_leg11*sizeof(PetscInt), &nv_leg11); */
/*   PetscMalloc(n_leg12*sizeof(PetscInt), &nv_leg12); */

/*   n_tail=0; */
/*   n_antene=0; */
/*   n_leg1=0;n_leg2=0;n_leg3=0;n_leg4=0; */
/*   n_leg5=0;n_leg6=0;n_leg7=0;n_leg8=0; */
/*   n_leg9=0;n_leg10=0;n_leg11=0;n_leg12=0; */


/*   for (n=0; n<n_v; n++) { */
/*     //!write(3,*) 'n,xf[n],yf[n],zf[n]',n,xf[n],yf[n],zf[n] */
/*     if(x_bp[n]>=x_tail_0) { */
/*       n_tail=n_tail+1; */
/*       nv_tail[n_tail-1]=n; */
/*     } */

/*     if(n>=n_antene_mn && n<=n_antene_mx) { */
/*       if(n>875) { */
/* 	n_antene=n_antene+1; */
/* 	nv_antene[n_antene-1]=n; */
/*       } */
/*     } */

/*     if(z_bp[n]<z_bot_0) { */
/*       if(n>=n1_min && n<=n1_max) { */
/* 	n_leg1=n_leg1+1; */
/* 	nv_leg1[n_leg1-1]=n; */
/*       } */
/*       if(n>=n2_min && n<=n2_max) { */
/* 	n_leg2=n_leg2+1; */
/* 	nv_leg2[n_leg2-1]=n; */
/*       } */
/*       if(n>=n3_min && n<=n3_max) { */
/* 	n_leg3=n_leg3+1; */
/* 	nv_leg3[n_leg3-1]=n; */
/*       } */
/*       if(n>=n4_min && n<=n4_max) { */
/* 	n_leg4=n_leg4+1; */
/* 	nv_leg4[n_leg4-1]=n; */
/*       } */
/*       if(n>=n5_min && n<=n5_max) { */
/* 	n_leg5=n_leg5+1; */
/* 	nv_leg5[n_leg5-1]=n; */
/*       } */
/*       if(n>=n6_min && n<=n6_max) { */
/* 	n_leg6=n_leg6+1; */
/* 	nv_leg6[n_leg6-1]=n; */
/*       } */
/*     } */

/*     if(z_bp[n]<z_bot_1) { */
/*       if(n>=n7_min && n<=n7_max) { */
/* 	n_leg7=n_leg7+1; */
/* 	nv_leg7[n_leg7-1]=n; */
/*       } */
/*       if(n>=n8_min && n<=n8_max) { */
/* 	n_leg8=n_leg8+1; */
/* 	nv_leg8[n_leg8-1]=n; */
/*       } */
/*     } */

/*     //!write(*,*) 'n,zf,z_bot_2',n,zf[n],z_bot_2 */
/*     if(n>=n9_min && n<=n9_max) { */
/*       //!if(n/=5664)  */
/*       n_leg9=n_leg9+1; */
/*       nv_leg9[n_leg9-1]=n; */
/*       //} */
/*     } */
/*     if(n>=n10_min && n<=n10_max) { */
/*       n_leg10=n_leg10+1; */
/*       nv_leg10[n_leg10-1]=n; */
/*     } */
   
/*     if(z_bp[n]<sz_l11) { */
/*       if(n>=n11_min && n<=n11_max) { */
/* 	n_leg11=n_leg11+1; */
/* 	nv_leg11[n_leg11-1]=n; */
/*       } */
/*       if(n>=n12_min && n<=n12_max) { */
/* 	n_leg12=n_leg12+1; */
/* 	nv_leg12[n_leg12-1]=n; */
/*       } */
/*     } */

/*   } */

  return(0);
}

    
PetscErrorCode cop_init_time(PetscInt regime, PetscReal *delti)
{

/* ==================================================================================             */
/*   cop_0210 */
/* ==================================================================================             */

  nfrq=1; N_period=500;
  tet_tmin=-20.0; tet_tmax=60.0; tet_amin=0.;//5.0; 
  tet_amax=70.0 ;                                                // escape regime
  tet_l12max=90.0; tet_l34max=90.0; tet_l56max=90.0; tet_l78max=90.0; tet_l910max=50.0; tet_l1112max=120.0;  // escape regime
  //tet_tmin=-0.0; tet_tmax=0.0; tet_amin=0.0; tet_amax=0.0                                              !cruising regime
  //tet_l12max=0.0; tet_l34max=0.0; tet_l56max=0.0; tet_l78max=0.0; tet_l910max=20.0  ;                   !cruising regime
  tet_l12min=0.0; tet_l34min=0.0; tet_l56min=0.0; tet_l78min=0.0; tet_l910min=0.0; tet_l1112min=0.0;
  tet_adif=15;
  t1_t=0.0;        t2_t=0.075;       t3_t=0.18;       t4_t=0.225;
  t1_a=0.0;        t2_a=0.075;       t3_a=0.145;      t4_a=0.3;
  t3d1_a=0.2; t3d2_a=0.25; t3d3_a=0.29;
  t3Rs_a=0.2; t3Re_a=0.285;
  t2d_a=0.1; t2R_a=0.1;

/*   t1_a=0.0;        t2_a=0.075;       t3_a=0.18;      t4_a=0.225; */
  t1_l12=0.101;    t2_l12=0.145;     t3_l12=0.145;    t4_l12=0.225;
  t1_l34=0.0925;   t2_l34=0.127;     t3_l34=0.145;    t4_l34=0.225;
  t1_l56=0.084;    t2_l56=0.110;     t3_l56=0.145;    t4_l56=0.225;
  t1_l78=0.075;    t2_l78=0.0925;    t3_l78=0.145;    t4_l78=0.225;
/*   t1_l910 =0.0;    t2_l910 =0.01;    t3_l910 =0.014;   t4_l910 =0.014; */
/*   t1_l1112=0.0;    t2_l1112=0.01;    t3_l1112=0.014;   t4_l1112=0.014; */
  t1_l910 =0.0;    t2_l910 =0.068;    t3_l910 =0.070;   t4_l910 =0.078;
  t1_l1112=0.0;    t2_l1112=0.068;    t3_l1112=0.070;   t4_l1112=0.078;

/* ==================================================================================             */


  if(*delti < 10.0)  {    //!Problem is nonstationary

    //! this is temporary change
    //!nfrq=10
    //write(*,*) 'Frequencies of copepod are change nfrq = ',nfrq
    //write(3,*) 'Frequencies of copepod are change nfrq = ',nfrq

    t1_t=t1_t*nfrq;t2_t=t2_t*nfrq;t3_t=t3_t*nfrq;t4_t=t4_t*nfrq;
    t1_a=t1_a*nfrq;t2_a=t2_a*nfrq;t3_a=t3_a*nfrq;t4_a=t4_a*nfrq;
    t3d1_a=t3d1_a*nfrq;t3Rs_a=t3Rs_a*nfrq;t3Re_a=t3Re_a*nfrq;
    t3d2_a=t3d2_a*nfrq;t3d3_a=t3d3_a*nfrq;
    t2d_a=t2d_a*nfrq;t2R_a=t2R_a*nfrq;

    t1_l12=t1_l12*nfrq;  t2_l12=t2_l12*nfrq;  t3_l12=t3_l12*nfrq;  t4_l12=t4_l12*nfrq;
    t1_l34=t1_l34*nfrq;  t2_l34=t2_l34*nfrq;  t3_l34=t3_l34*nfrq;  t4_l34=t4_l34*nfrq;
    t1_l56=t1_l56*nfrq;  t2_l56=t2_l56*nfrq;  t3_l56=t3_l56*nfrq;  t4_l56=t4_l56*nfrq;
    t1_l78=t1_l78*nfrq;  t2_l78=t2_l78*nfrq;  t3_l78=t3_l78*nfrq;  t4_l78=t4_l78*nfrq;
    t1_l910=t1_l910*nfrq;t2_l910=t2_l910*nfrq;t3_l910=t3_l910*nfrq;t4_l910=t4_l910*nfrq;

    if(regime==0 || regime==1){
      T_period=t4_a-t1_a;   //    ! for  escape regime
/*       write(*,*) '********* Copepode. Escape regime ********: T_period=', T_period,' f=',1./T_period */
/*       write(3,*) '********* Copepode. Escape regime ********: T_period=', T_period,' f=',1./T_period */

      // ! to exclude of chilla moving
      t1_l910 =-1.0;
      t1_l1112=-1.0;
    }
    else if(regime==2) {
      T_period=t4_l910-t1_l910; //     ! for cruise regime
/*       write(*,*) '********* Copepode. Cruise regime ********: T_period=', T_period,' f=',1./T_period */
/*       write(3,*) '********* Copepode. Cruise regime ********: T_period=', T_period,' f=',1./T_period */

      // ! to exclude of anntenne, legs, tail moving
      t1_t=-1.0;
      t1_a=-1.0;
      t1_l12=-1.0;
      t1_l34=-1.0;
      t1_l56=-1.0;
      t1_l78=-1.0;
    }

    *delti = T_period/N_period;
    //s_work = 0.0 ;
    tet_t=0.0;
    tet_a=tet_amin;
    tet_l12=tet_l12min;
    tet_l34=tet_l34min;
    tet_l56=tet_l56min ;
    tet_l78=tet_l78min;
    tet_l910=tet_l910min;

/*     write(*,*) 'Nonstationary problem:' */
/*     write(*,*) 'delti, T_period=',delti, T_period,'N steps=',T_period/delti */
/*     write(3,*)  */
/*     write(3,*) 'Nonstationary problem:' */
/*     write(3,*) 'delti, T_period=',delti, T_period,'N steps=',T_period/delti */
/*     write(3,*)  */
  }
  else {
/*     write(*,*) 'Copepode. Stationary problem: delti', delti */
/*     write(3,*) 'Copepode. Stationary problem: delti', delti */
    T_period=*delti;
    //ijk_direction=-1;
  }

  PetscPrintf(PETSC_COMM_WORLD, "Perid T %le , N Period %le , dt %le \n",T_period, N_period, *delti );

  return(0);
}

PetscErrorCode rotate_cop_app_old(IBMNodes *ibm,PetscInt n_cop, 
				  PetscInt *nv_cop, PetscReal dphix,
				  PetscReal dphiy,PetscReal dphiz,
				  PetscReal sx1,PetscReal sy1,
				  PetscReal sz1,PetscReal delti)
{
/* program of a simple rotation on an angle: dphix,dphiy,dphiz */
  PetscReal *xf,*yf,*zf;
  PetscReal v_max;
  PetscInt  n,l;
  PetscReal xflm,yflm,zflm;
  //PetscReal uflm,vflm,wflm;
  PetscReal x1,x2,x3,x4;
  PetscReal y1,y2,y3,y4;
  PetscReal z1,z2,z3,z4;
  //Cmpnts    *uf;

  xf = ibm->x_bp0; yf=ibm->y_bp0; zf=ibm->z_bp0;
  //  uf = ibm->u;

  v_max=0.0;
  for (n=0;n<n_cop; n++) {
    l=nv_cop[n];
    xflm=xf[l]; yflm=yf[l]; zflm=zf[l];

    x1=xf[l]-sx1; y1=yf[l]-sy1; z1=zf[l]-sz1;
    //! rot. of x
    x2= x1;
    y2= y1*cos(dphix)+z1*sin(dphix);
    z2=-y1*sin(dphix)+z1*cos(dphix);
    //! rot. of y
    x3= x2*cos(dphiy)-z2*sin(dphiy);
    y3= y2;
    z3= x2*sin(dphiy)+z2*cos(dphiy);
    //! rot. of z
    if(y1>=0.0) {
      x4= x3*cos(dphiz)+y3*sin(dphiz);
      y4=-x3*sin(dphiz)+y3*cos(dphiz);
      z4= z3;
    }
    else {
      x4= x3*cos(-dphiz)+y3*sin(-dphiz);
      y4=-x3*sin(-dphiz)+y3*cos(-dphiz);
      z4= z3;
    }
    xf[l]=x4+sx1;
    yf[l]=y4+sy1;
    zf[l]=z4+sz1;

/*     uflm=uf[l].x; */
/*     vflm=uf[l].y; */
/*     wflm=uf[l].z; */

/*     uf[l].x=(xf[l]-xflm)/delti; */
/*     uf[l].y=(yf[l]-yflm)/delti; */
/*     uf[l].z=(zf[l]-zflm)/delti; */

/*     dudtf[l]=(uf[l]-uflm)/delti; */
/*     dvdtf[l]=(vf[l]-vflm)/delti; */
/*     dwdtf[l]=(wf[l]-wflm)/delti; */

/*     v_max=PetscMax(fabs(uf[l].x),fabs(uf[l].y)); */
/*     v_max=PetscMax(v_max,fabs(uf[l].z)); */

/*     PetscPrintf(PETSC_COMM_WORLD, "Cop Body Rotate: %le vel %le %le %le\n",dphi); */

  }
/*     write(*,*)' v_max=',v_max */

/* if(.false.) then */
/*     ! changing of surface */
/*     sign= dphiy/abs(dphiy)*0.1 */
/*     yflm=0.0 */
/*     do n=1,n_cop */
/*      l=nv_cop[n] */
/*      yflm=yflm+yf(l) */
/*     end do */
/*     yflm=yflm/n_cop */
    
/*     write(*,*)' mean value of y: yflm=',yflm,'sign=',sign */

/*     do n=1,n_cop */
/*      l=nv_cop[n] */
/*      yf(l)=yf(l) - sign*(yflm-yf(l)) */
/*     end do */
/*     yflm=yflm/n_cop */
/* end if */

  return(0);
}

PetscErrorCode deform_antene(IBMNodes *ibm, PetscReal A,
			     PetscReal L,PetscReal R,
			     PetscReal sx1,PetscReal sy1,
			     PetscReal sz1)

{
  PetscReal  pi = 3.141592653589793;
  PetscReal *xf,*yf,*zf;
  PetscInt  l,n_app  ;
  PetscReal x1,y1,z1;

  n_app= ibm->n_v;
  xf = ibm->x_bp; yf=ibm->y_bp; zf=ibm->z_bp;

  for (l=0;l<n_app; l++) {

    x1=xf[l]-sx1; y1=yf[l]-sy1; z1=zf[l]-sz1;

    x1=R*x1;
    if (x1>0)
      z1=z1-A*sin(pi*x1/L);
    else
      z1=z1+A*sin(pi*x1/L);

    xf[l]=x1+sx1;
    yf[l]=y1+sy1;
    zf[l]=z1+sz1;    
  }
  return 0;
}
    

PetscErrorCode rotate_antene(IBMNodes *ibm,
			     PetscReal dphix,
			     PetscReal dphiy,PetscReal dphiz,
			     PetscReal sx1,PetscReal sy1,
			     PetscReal sz1,PetscReal delti)
{
/* program of a simple rotation on an angle: dphix,dphiy,dphiz */
  PetscReal *xf,*yf,*zf;
  //  PetscReal v_max;
  PetscInt  l,n_app  ;
  PetscReal xflm,yflm,zflm;
  //PetscReal uflm,vflm,wflm;
  PetscReal x1,x2,x3,x4;
  PetscReal y1,y2,y3,y4;
  PetscReal z1,z2,z3,z4;
  //Cmpnts    *uf;

  n_app= ibm->n_v;
  xf = ibm->x_bp; yf=ibm->y_bp; zf=ibm->z_bp;
  //  uf = ibm->u;
  
  for (l=0;l<n_app; l++) {   
    xflm=xf[l]; yflm=yf[l]; zflm=zf[l];

    x1=xf[l]-sx1; y1=yf[l]-sy1; z1=zf[l]-sz1;
    //! rot. of x
    x2= x1;
    y2= y1*cos(dphix)+z1*sin(dphix);
    z2=-y1*sin(dphix)+z1*cos(dphix);
    //! rot. of y
    if (x1>=0.) {
      x3= x2*cos(dphiy)-z2*sin(dphiy);
      y3= y2;
      z3= x2*sin(dphiy)+z2*cos(dphiy);
    } else {
      x3= x2*cos(-dphiy)-z2*sin(-dphiy);
      y3= y2;
      z3= x2*sin(-dphiy)+z2*cos(-dphiy);
    }
    //! rot. of z
    if(y1>=0.0) {
      x4= x3*cos(dphiz)+y3*sin(dphiz);
      y4=-x3*sin(dphiz)+y3*cos(dphiz);
      z4= z3;
    }
    else {
      x4= x3*cos(-dphiz)+y3*sin(-dphiz);
      y4=-x3*sin(-dphiz)+y3*cos(-dphiz);
      z4= z3;
    }

/*     ibm->x_bp[l] = y4+sy1;//yf[i]; */
/*     ibm->y_bp[l] = z4+sz1;//zf[i]; */
/*     ibm->z_bp[l] = x4+sx1;//xf[i]; */
 
    xf[l]=x4+sx1;
    yf[l]=y4+sy1;
    zf[l]=z4+sz1;
  }

  /*    for (i=0;i<n_app;i++) { */
  /*      ibm->x_bp[i] = yf[i]; */
  /*      ibm->y_bp[i] = zf[i]; */
  /*      ibm->z_bp[i] = xf[i]; */
  /*    } */

  return(0);
}


PetscErrorCode rotate_cop_app(IBMNodes *ibm,
			      PetscReal dphix,
			      PetscReal dphiy,PetscReal dphiz,
			      PetscReal sx1,PetscReal sy1,
			      PetscReal sz1,PetscReal delti)
{
/* program of a simple rotation on an angle: dphix,dphiy,dphiz */
  PetscReal *xf,*yf,*zf;
  //  PetscReal v_max;
  PetscInt  i,l,n_app  ;
  PetscReal xflm,yflm,zflm;
  //PetscReal uflm,vflm,wflm;
  PetscReal x1,x2,x3,x4;
  PetscReal y1,y2,y3,y4;
  PetscReal z1,z2,z3,z4;
  //Cmpnts    *uf;

  n_app= ibm->n_v;
  xf = ibm->x_bp0; yf=ibm->y_bp0; zf=ibm->z_bp0;
  //  uf = ibm->u;
  
  for (l=0;l<n_app; l++) {   
    xflm=xf[l]; yflm=yf[l]; zflm=zf[l];

    x1=xf[l]-sx1; y1=yf[l]-sy1; z1=zf[l]-sz1;
    //! rot. of x
    x2= x1;
    y2= y1*cos(dphix)+z1*sin(dphix);
    z2=-y1*sin(dphix)+z1*cos(dphix);
    //! rot. of y
    x3= x2*cos(dphiy)-z2*sin(dphiy);
    y3= y2;
    z3= x2*sin(dphiy)+z2*cos(dphiy);
    //! rot. of z
    if(y1>=0.0) {
      x4= x3*cos(dphiz)+y3*sin(dphiz);
      y4=-x3*sin(dphiz)+y3*cos(dphiz);
      z4= z3;
    }
    else {
      x4= x3*cos(-dphiz)+y3*sin(-dphiz);
      y4=-x3*sin(-dphiz)+y3*cos(-dphiz);
      z4= z3;
    }
    xf[l]=x4+sx1;
    yf[l]=y4+sy1;
    zf[l]=z4+sz1;
  }

   for (i=0;i<n_app;i++) {
     ibm->x_bp[i] = yf[i];
     ibm->y_bp[i] = zf[i];
     ibm->z_bp[i] = xf[i];
   }

  return(0);
}


PetscErrorCode cop_swim(IBMNodes *ibm, PetscReal time
			,PetscReal delti)
{
  PetscReal  pi = 3.141592653589793;
  PetscReal  pi0180 = pi/180.0;
  PetscReal  dphix,dphiy,dphiz;
  PetscReal  time_i, tet_ti, tet_ai;
  PetscReal  tet_l12i,tet_l34i,tet_l56i,tet_l78i;
  PetscInt   i,  ibi;
  PetscReal  v_max;
  PetscInt   i_vmax;

/*   x_bp=ibm->x_bp0;y_bp=ibm->y_bp0;z_bp=ibm->z_bp0; */

/*  write(3,*) */
/*  write(3,*) 'time',time */
/*  write(3,*) */
  //ang_vel=0.0;
  //!T_period=t4_a-t1_a      ! for  cruise regime

/*    PetscPrintf(PETSC_COMM_WORLD, "MAX Cop Velocity: \n"); */

  if(t1_t>=0.0) {
/*     !tail */
/*     !T_period=t4_t-t1_t */
    time_i=time/T_period-(int)(time/T_period);
    if(time>0.0 && time_i==0.0) time_i=1.0;
    time_i=time_i*T_period;
    //!time_i=time
    //write(3,*) 'tail:t1_t,time_i,t4_t',t1_t,time_i,t4_t
    if(t1_t<=time_i && time_i<=t2_t)
      tet_ti=((tet_tmax-0.0)/(t2_t-t1_t))*(time_i-t1_t);
    else if(t2_t<=time_i && time_i<=t3_t)
      tet_ti=tet_tmax+((tet_tmin-tet_tmax)/(t3_t-t2_t))*(time_i-t2_t);
    else if(t3_t<=time_i && time_i<=t4_t)
      tet_ti=tet_tmin+((0.0-tet_tmin)/(t4_t-t3_t))*(time_i-t3_t);
     
    if(t1_t<=time_i && time_i<=t4_t) {
      dphix=0.0;
      dphiy=(tet_ti-tet_t)*pi0180;
      dphiz=0.0;
/*        write(*,*) 'tet_ti,tet_t,delt',tet_ti,tet_t,tet_ti-tet_t */
      rotate_cop_app(&ibm[1],
	  dphix,dphiy,dphiz,sx_tail,sy_tail,sz_tail,delti);
      tet_t=tet_ti;
      //ang_vel(1)=dphiy;
    }
  }

  //antene
  if(t1_a>=0.0) {

    PetscReal L0,L,R,A;
    L0=0.75;    
    // !T_period=t4_a-t1_a
    time_i=time/T_period-(int)(time/T_period);
    if(time>=0.0 && time_i==0.0) time_i=1.0;
    time_i=time_i*T_period;

    //!time_i=time
    if(t1_a<=time_i && time_i<=t2_a)
      tet_ai=((tet_amax-tet_amin)/t2_t)*time_i+tet_amin;
    
    else if(t2_a<=time_i && time_i<=t3_a){
      if (time_i<=t2R_a)      
	tet_ai=tet_amax+tet_adif/(t2R_a-t2_a)*(time_i-t2_a);
      else
	tet_ai=tet_amax+tet_adif;

      if (time_i<=t2d_a) {
	R=1.;
	L=R*L0;
	A=L0/8./(t2d_a-t2_a)*(time_i-t2_a);
      } else {
	R=1.;
	L=R*L0;
	A=L0/8.;
      }
    }
    
    else if(t3_a<time_i && time_i<=t4_a) {
      if (time_i <= t3d1_a) {
	R=((0.4-1.)/(t3d1_a-t3_a))*(time_i-t3_a)+1.;
	L=R*L0;
	A=L0/(8.*R);
	//	PetscPrintf(PETSC_COMM_WORLD, "Deform R %le A %le L %le time_i %le time_3d \n",R,A,L,time_i,t3d_a);
      }	else if (time_i>t3d2_a && time_i<t3d3_a){    
	R=((0.4-1.)/(t3d2_a-t3d3_a))*(time_i-t3d3_a)+1.;
	L=R*L0;
	A=L0/(8.*R);
      } else if (time_i>=t3d3_a){
	R=1.;
	L=R*L0;
	A=L0/8./(t3d3_a-t4_a)*(time_i-t4_a);	
      } else {
	R=0.4;//((0.2-1.)/(t3d_a-t3_a))*(time_i-t3_a)+1.;
	L=R*L0;
	A=L0/(8.*R);
      }

      if (time_i>=t3Rs_a && time_i<t3Re_a) 
	tet_ai=tet_amax+tet_adif+
	  ((-tet_amax-tet_adif+tet_amin)/(t3Re_a-t3Rs_a))*(time_i-t3Rs_a);
      else if (time_i<t3Rs_a)
	tet_ai=tet_amax+tet_adif;
      else if (time_i>=t3Re_a)
	tet_ai=tet_amin;
    }

    for (i=0; i<ibm[2].n_v; i++) {    
      ibm[2].x_bp[i] = ibm[2].y_bp0[i];
      ibm[2].y_bp[i] = ibm[2].z_bp0[i];
      ibm[2].z_bp[i] = ibm[2].x_bp0[i];
    }
    
    if(t2_a<time_i && time_i<=t4_a) {
      deform_antene(&ibm[2],A,L,R,
		    sx_antene,sy_antene,sz_antene);
      PetscPrintf(PETSC_COMM_WORLD, "Deform R %le A %le L %le tetha %le  t3d_a %le %le\n",R,A,L,tet_ai,t3d1_a,time_i);
    }

    if(t1_a<=time_i && time_i<=t4_a) {
      dphix=0.0;
      //dphiy=0.0;
      dphiz=0.0;
      //dphiz=(tet_ai-tet_a)*pi0180;
      dphiy=(tet_ai)*pi0180;

      rotate_antene(&ibm[2],dphix,dphiy,dphiz,
		    sx_antene,sy_antene,sz_antene,delti);
      /*       rotate_cop_app(&ibm[2], */
      /* 	  dphix,dphiy,dphiz,sx_antene,sy_antene,sz_antene,delti); */
      /*       tet_a=tet_ai; */
    }
  }
  
  /*  legs ============================================================================= */
/*  l_12 */
  if(t1_l12>=0.0) {
/*     sx_l1 = x_bp[n_leg_1_rot];   sy_l1=y_bp[n_leg_1_rot];    sz_l1=z_bp[n_leg_1_rot]; */
/*     sx_l2 = x_bp[n_leg_2_rot];   sy_l2=y_bp[n_leg_2_rot];    sz_l2=z_bp[n_leg_2_rot]; */

    // !T_period=t4_l12-t1_l12
    time_i=time/T_period-(int)(time/T_period);
    if(time>0.0 && time_i==0.0) time_i=1.0;
    time_i=time_i*T_period;
    //!time_i=time
    //write(3,*) 'l12:t1_l12,time_i,t4_l12',t1_l12,time_i,t4_l12
    if(t1_l12<=time_i && time_i<=t2_l12)
      tet_l12i=tet_l12min+((tet_l12max-tet_l12min)/(t2_l12-t1_l12))*(time_i-t1_l12);
    else if(t2_l12<=time_i && time_i<=t3_l12)
      tet_l12i=tet_l12max;
    else if(t3_l12<=time_i && time_i<=t4_l12)
      tet_l12i=tet_l12max+((tet_l12min-tet_l12max)/(t4_l12-t3_l12))*(time_i-t3_l12);
    
    if(t1_l12<=time_i && time_i<=t4_l12) {
      dphix=0.0;
      dphiy=(tet_l12i-tet_l12)*pi0180;
      dphiz=0.0;
/*        write(*,*) 'tet_l12i,tet_l12,delt',tet_l12i,tet_l12,tet_l12i-tet_l12 */
      rotate_cop_app(&ibm[3],
	  dphix,dphiy,dphiz,sx_l1,sy_l1,sz_l1,delti);
/*       rotate_cop_app(ibm,n_leg2,nv_leg2, */
/* 	  dphix,dphiy,dphiz,sx_l2,sy_l2,sz_l2,delti); */
	tet_l12=tet_l12i;
	//ang_vel(3)=dphiy;
    }
  }

  /*  l_34 */
  if(t1_l34>=0.0) {
/*     sx_l3 = x_bp[n_leg_3_rot];   sy_l3=y_bp[n_leg_3_rot];    sz_l3=z_bp[n_leg_3_rot]; */
/*     sx_l4 = x_bp[n_leg_4_rot];   sy_l4=y_bp[n_leg_4_rot];    sz_l4=z_bp[n_leg_4_rot]; */

    //!T_period=t4_l34-t1_l34
    time_i=time/T_period-(int)(time/T_period);
    if(time>0.0 && time_i==0.0) time_i=1.0;
    time_i=time_i*T_period;
    //!time_i=time
    //write(3,*) 'l34:t1_l34,time_i,t4_l34',t1_l34,time_i,t4_l34
    if(t1_l34<=time_i && time_i<=t2_l34)
      tet_l34i=tet_l34min+((tet_l34max-tet_l34min)/(t2_l34-t1_l34))*(time_i-t1_l34);
    else if(t2_l34<=time_i && time_i<=t3_l34)
      tet_l34i=tet_l34max;
    else if(t3_l34<=time_i && time_i<=t4_l34)
      tet_l34i=tet_l34max+((tet_l34min-tet_l34max)/(t4_l34-t3_l34))*(time_i-t3_l34);
    
    if(t1_l34<=time_i && time_i<=t4_l34) {
      /*        write(*,*) 'tet_l34i,tet_l34,delt',tet_l34i,tet_l34,tet_l34i-tet_l34 */
      dphix=0.0;
      dphiy=(tet_l34i-tet_l34)*pi0180;
      dphiz=0.0;
      rotate_cop_app(&ibm[4],
	  dphix,dphiy,dphiz,sx_l3,sy_l3,sz_l3,delti);
/*       rotate_cop_app(ibm,n_leg4,nv_leg4, */
/* 	  dphix,dphiy,dphiz,sx_l4,sy_l4,sz_l4,delti); */
      tet_l34=tet_l34i;
      //ang_vel(4)=dphiy;
    }
  }
  
/*  !l_56 */
  if(t1_l56>=0.0) {
/*     sx_l5 = x_bp[n_leg_5_rot];   sy_l5=y_bp[n_leg_5_rot];    sz_l5=z_bp[n_leg_5_rot]; */
/*     sx_l6 = x_bp[n_leg_6_rot];   sy_l6=y_bp[n_leg_6_rot];    sz_l6=z_bp[n_leg_6_rot]; */

    //!T_period=t4_l56-t1_l56
    time_i=time/T_period-(int)(time/T_period);
    if(time>0.0 && time_i==0.0) time_i=1.0;
    time_i=time_i*T_period;
    //!time_i=time
/*     write(3,*) 'l56:t1_l56,time_i,t4_l56',t1_l56,time_i,t4_l56; */
    if(t1_l56<=time_i && time_i<=t2_l56)
      tet_l56i=tet_l56min+((tet_l56max-tet_l56min)/(t2_l56-t1_l56))*(time_i-t1_l56);
    else if(t2_l56<=time_i && time_i<=t3_l56)
      tet_l56i=tet_l56max;
    else if(t3_l56<=time_i && time_i<=t4_l56)
      tet_l56i=tet_l56max+((tet_l56min-tet_l56max)/(t4_l56-t3_l56))*(time_i-t3_l56);
    
    if(t1_l56<=time_i && time_i<=t4_l56) {
      /*        write(*,*) 'tet_l56i,tet_l56,delt',tet_l56i,tet_l56,tet_l56i-tet_l56; */
      dphix=0.0;
      dphiy=(tet_l56i-tet_l56)*pi0180;
      dphiz=0.0;
      rotate_cop_app(&ibm[5],
	  dphix,dphiy,dphiz,sx_l5,sy_l5,sz_l5,delti);
/*       rotate_cop_app(ibm,n_leg6,nv_leg6, */
/* 	  dphix,dphiy,dphiz,sx_l6,sy_l6,sz_l6,delti); */
      tet_l56=tet_l56i;
      //ang_vel(5)=dphiy;
    }
  }

/*  !l_78 */
  if(t1_l78>=0.0) {
/*     sx_l7 = x_bp[n_leg_7_rot];   sy_l7=y_bp[n_leg_7_rot];    sz_l7=z_bp[n_leg_7_rot]; */
/*     sx_l8 = x_bp[n_leg_8_rot];   sy_l8=y_bp[n_leg_8_rot];    sz_l8=z_bp[n_leg_8_rot]; */

    //!T_period=t4_l78-t1_l78
    time_i=time/T_period-(int)(time/T_period);
    if(time>0.0 && time_i==0.0) time_i=1.0;
    time_i=time_i*T_period;
    //!time_i=time
    //write(3,*) 'l78:t1_l78,time_i,t4_l78',t1_l78,time_i,t4_l78
    if(t1_l78<=time_i && time_i<=t2_l78)
      tet_l78i=tet_l78min+((tet_l78max-tet_l78min)/(t2_l78-t1_l78))*(time_i-t1_l78);
    else if(t2_l78<=time_i && time_i<=t3_l78)
      tet_l78i=tet_l78max;
    else if(t3_l78<=time_i && time_i<=t4_l78)
      tet_l78i=tet_l78max+((tet_l78min-tet_l78max)/(t4_l78-t3_l78))*(time_i-t3_l78);
    
    if(t1_l78<=time_i && time_i<=t4_l78) {
/*        write(*,*) 'tet_l78i,tet_l78,delt',tet_l78i,tet_l78,tet_l78i-tet_l78 */
      dphix=0.0;
      dphiy=(tet_l78i-tet_l78)*pi0180;
      dphiz=0.0;
      rotate_cop_app(&ibm[6],
	  dphix,dphiy,dphiz,sx_l7,sy_l7,sz_l7,delti);
/*       rotate_cop_app(ibm,n_leg8,nv_leg8, */
/* 	  dphix,dphiy,dphiz,sx_l8,sy_l8,sz_l8,delti); */
      tet_l78=tet_l78i;
      //ang_vel(6)=dphiy;
    }
  }

/* /\* ==============================  Cruise regime   ==================================== *\/ */
/* /\*    por=0.0; *\/ */

/* /\* ------------------------------------------------ 9--10 ---------------------------------------- *\/ */

/*    if(t1_l910>=0.0) { */
/*      sx_l9 = x_bp[n_leg_9_rot];   sy_l9=y_bp[n_leg_9_rot];    sz_l9=z_bp[n_leg_9_rot]; */
/*      sx_l10= x_bp[n_leg_10_rot];  sy_l10=y_bp[n_leg_10_rot];  sz_l10=z_bp[n_leg_10_rot]; */

/*      //!if(t2_l910<=time_i && time_i<=t4_l910) por=1.0 */

/*      //!l_910 */
/*      //!T_period=t4_l910-t1_l910     ! for cruise regime */
/*      time_i=time/T_period-(int)(time/T_period); */
/*      if(time>0.0 && time_i==0.0) time_i=1.0; */
/*      time_i=time_i*T_period; */
/* /\*      write(3,*) 'l910:t1_l910,time_i,t4_l910',t1_l910,time_i,t4_l910 *\/ */
/* /\*      write(*,*) 'l910:t1_l910,time_i,t4_l910',t1_l910,time_i,t4_l910 *\/ */

/* /\*      write(*,*) 'l910:t1_l910-dt,time_i,t1_l910+dt',t1_l910-delti,time_i,t1_l910+delti *\/ */
/* /\*      !n_arise=0  *\/ */

/*      if(0)  {  //! I don't use deminishing of legs */
/*                    //!if(t1_l910-delti<=time_i && time_i<=t1_l910+delti) then */
/*        if(t1_l910<=time_i && time_i<=t1_l910+1.1*delti) { */
/* 	 //!arise of legs */
/* /\* 	 n_arise=1 ; *\/ */
/* /\* 	 por=1.0; //  ! to exclude arising of perturbation *\/ */
/*                   //   !cop_arise(n_leg9, nv_leg9, lv_max,xf_nd,yf_nd,zf_nd,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf) */
/*                   //   !cop_arise(n_leg10,nv_leg10,lv_max,xf_nd,yf_nd,zf_nd,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf) */
/*        } */

/*         /\*               write(*,*) 'l910:t2_l910-dt,time_i,t2_l910+dt',t2_l910-delti,time_i,t2_l910+delti *\/ */
/* /\*                      !if(t2_l910-delti<time_i && time_i<t2_l910+delti) then  *\/ */
/* /\*                      !if(t2_l910-delti<time_i && time_i<t2_l910) then  *\/ */
/*        if(t2_l910-delti<time_i && time_i<t4_l910){ */
/* 	 /\*        !diminish of legs *\/ */
/* 	 /\* 	 !cop_diminish(n_leg9, nv_leg9, lv_max,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf,sx_l9,sy_l9,sz_l9) *\/ */
/* 	 /\* 	 !cop_diminish(n_leg10,nv_leg10,lv_max,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf,sx_l10,sy_l10,sz_l10) *\/ */
/*        } */
/*      } */

/*      //!if(n_arise==0) then */
/*          if(t1_l910<=time_i && time_i<=t2_l910) */
/*            tet_l910i=tet_l910min+((tet_l910max-tet_l910min)/(t2_l910-t1_l910))*(time_i-t1_l910); */
/*          else if(t2_l910<=time_i && time_i<=t3_l910) */
/*            tet_l910i=tet_l910max; */
/*          else if(t3_l910<=time_i && time_i<=t4_l910) */
/*            tet_l910i=tet_l910max+((tet_l910min-tet_l910max)/(t4_l910-t3_l910))*(time_i-t3_l910); */
         

/*          if(t1_l910<=time_i && time_i<=t4_l910) { */
/* /\*            write(3,*) 'tet_l910i,tet_l910,delt',tet_l910i,tet_l910,tet_l910i-tet_l910 *\/ */
/* /\*            write(*,*) 'tet_l910i,tet_l910,delt',tet_l910i,tet_l910,tet_l910i-tet_l910 *\/ */
/*            dphix=0.0; */
/*            dphiy=0.0; */
/*            dphiz=(tet_l910i-tet_l910)*pi0180; */
/* /\*            if(dphiz<=0.0) por=1.0; *\/ */
/*            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!if(t1_l910<=time_i && time_i<=t2_l910)& */
/*            if(t1_l910<=time_i && time_i<=t4_l910) */
/*              rotate_cop_app(ibm,n_leg9,nv_leg9,dphix,dphiy,dphiz,sx_l9,sy_l9,sz_l9,delti); */
/*            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!if(t1_l910<=time_i && time_i<=t2_l910)& */
/*            if(t1_l910<=time_i && time_i<=t4_l910) */
/*              rotate_cop_app(ibm,n_leg10,nv_leg10,dphix,dphiy,dphiz,sx_l10,sy_l10,sz_l10,delti); */
/*            tet_l910=tet_l910i; */
/*            //ang_vel(7)=dphiy; */
/*          } */
/* 	 //!end if */
/*    } */
  
/* /\* ------------------------------------------------ 11--12 ---------------------------------------- *\/ */
/* /\*  l_11-12 *\/ */
/*    if(t1_l1112>=0.0) { */
/*      sx_l11= x_bp[n_leg_11_rot];  sy_l11=y_bp[n_leg_11_rot];  sz_l11=z_bp[n_leg_11_rot]; */
/*      sx_l12= x_bp[n_leg_12_rot];  sy_l12=y_bp[n_leg_12_rot];  sz_l12=z_bp[n_leg_12_rot]; */

/*      time_i=time/T_period-(int)(time/T_period); */
/*      if(time>0.0 && time_i==0.0) time_i=1.0; */
/*      time_i=time_i*T_period; */
/* /\*      write(3,*) 'l11-12:t1_l1112,time_i,t4_l1112',t1_l1112,time_i,t4_l1112 *\/ */
/* /\*      write(*,*) 'l11-12:t1_l1112,time_i,t4_l1112',t1_l1112,time_i,t4_l1112 *\/ */

/* /\*                 if(.false.) then   *\/ */
/* /\*                  write(*,*) 'l1112:t1_l1112-dt,time_i,t1_l1112+dt',t1_l1112-delti,time_i,t1_l1112+delti *\/ */
/* /\*                  !if(t1_l1112-delti<=time_i && time_i<=t1_l1112+delti) then  *\/ */
/* /\*                  if(t1_l1112-delti<=time_i && time_i<=t1_l1112+1.1*delti) then  *\/ */
/* /\*                    !arise of legs *\/ */
/* /\*                    por=1.0  ! to exclude arising of perturbation *\/ */
/* /\*                    !cop_arise(n_leg11,nv_leg11,lv_max,xf_nd,yf_nd,zf_nd,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf) *\/ */
/* /\*                    !cop_arise(n_leg12,nv_leg12,lv_max,xf_nd,yf_nd,zf_nd,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf) *\/ */
/* /\*                  end if *\/ */

/* /\*                  write(*,*) 'l1112:t2_l1112-dt,time_i,t2_l1112+dt',t2_l1112-delti,time_i,t2_l1112+delti *\/ */
/* /\*                  if(t2_l1112-delti<time_i && time_i<t4_l1112) then  *\/ */
/* /\*                    !diminish of legs *\/ */
/* /\*                    !cop_diminish(n_leg11,nv_leg11,lv_max,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf,sx_l11,sy_l11,sz_l11) *\/ */
/* /\*                    !cop_diminish(n_leg12,nv_leg12,lv_max,xf,yf,zf,uf,vf,wf,dudtf,dvdtf,dwdtf,sx_l12,sy_l12,sz_l12) *\/ */
/* /\*                  end if *\/ */
/* /\*                 end if *\/ */

/*      //!if(n_arise==0) then */
/*         if(t1_l1112<=time_i && time_i<=t2_l1112) */
/* 	  tet_l1112i=tet_l1112min+((tet_l1112max-tet_l1112min)/(t2_l1112-t1_l1112))*(time_i-t1_l1112); */
/* 	else if(t2_l1112<=time_i && time_i<=t3_l1112) */
/* 	  tet_l1112i=tet_l1112max; */
/* 	else if(t3_l1112<=time_i && time_i<=t4_l1112) */
/* 	  tet_l1112i=tet_l1112max+((tet_l1112min-tet_l1112max)/(t4_l1112-t3_l1112))*(time_i-t3_l1112); */
	

/* 	if(t1_l1112<=time_i && time_i<=t4_l1112) { */
/* /\*            write(3,*) 'tet_l1112i,tet_l1112,delt',tet_l1112i,tet_l1112,tet_l1112i-tet_l1112 *\/ */
/* /\*            write(*,*) 'tet_l1112i,tet_l1112,delt',tet_l1112i,tet_l1112,tet_l1112i-tet_l1112 *\/ */
/* 	  dphix=0.0; */
/* 	  dphiy=(tet_l1112i-tet_l1112)*pi0180; */
/* 	  dphiz=0.0; */
/* 	  //!!!!!!!!!!!!!!!!!!!!!!!!if(t1_l1112<=time_i && time_i<=t2_l1112)& */
/*            if(t1_l1112<=time_i && time_i<=t4_l1112) */
/*              rotate_cop_app(ibm,n_leg11,nv_leg11,dphix,dphiy,dphiz,sx_l11,sy_l11,sz_l11,delti); */
/*            //!!!!!!!!!!!!!!!!!!!!!!!!!if(t1_l1112<=time_i && time_i<=t2_l1112) */
/*            if(t1_l1112<=time_i && time_i<=t4_l1112) */
/*              rotate_cop_app(ibm,n_leg12,nv_leg12,dphix,dphiy,dphiz,sx_l12,sy_l12,sz_l12,delti); */
/*            tet_l1112=tet_l1112i; */
/*            //ang_vel(7)=dphiy; */
/* 	} */
/* 	//!end if */
/*    } */


   PetscInt n1e, n2e, n3e;
   PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
   FILE *f;
   char filen[80];
   PetscInt rank;

   for (ibi=0;ibi<NumberOfBodies;ibi++) {
     
     /*    Calculate the new normal & velcity */
     for (i=0; i<ibm[ibi].n_elmt; i++) {
       n1e = ibm[ibi].nv1[i]; n2e =ibm[ibi].nv2[i]; n3e =ibm[ibi].nv3[i];
       dx12 = ibm[ibi].x_bp[n2e] - ibm[ibi].x_bp[n1e]; 
       dy12 = ibm[ibi].y_bp[n2e] - ibm[ibi].y_bp[n1e]; 
       dz12 = ibm[ibi].z_bp[n2e] - ibm[ibi].z_bp[n1e]; 
       
       dx13 = ibm[ibi].x_bp[n3e] - ibm[ibi].x_bp[n1e]; 
       dy13 = ibm[ibi].y_bp[n3e] - ibm[ibi].y_bp[n1e]; 
       dz13 = ibm[ibi].z_bp[n3e] - ibm[ibi].z_bp[n1e]; 
       
       ibm[ibi].nf_x[i] = dy12 * dz13 - dz12 * dy13;
       ibm[ibi].nf_y[i] = -dx12 * dz13 + dz12 * dx13;
       ibm[ibi].nf_z[i] = dx12 * dy13 - dy12 * dx13;
       
       dr = sqrt(ibm[ibi].nf_x[i]*ibm[ibi].nf_x[i] + ibm[ibi].nf_y[i]*ibm[ibi].nf_y[i] + 
		 ibm[ibi].nf_z[i]*ibm[ibi].nf_z[i]);
       
       ibm[ibi].nf_x[i] /=dr; ibm[ibi].nf_y[i]/=dr; ibm[ibi].nf_z[i]/=dr;
     }
     
     for (i=0; i<ibm[ibi].n_v; i++) {
       ibm[ibi].u[i].x = (ibm[ibi].x_bp[i] - ibm[ibi].x_bp_o[i]) / delti;
       ibm[ibi].u[i].y = (ibm[ibi].y_bp[i] - ibm[ibi].y_bp_o[i]) / delti;
       ibm[ibi].u[i].z = (ibm[ibi].z_bp[i] - ibm[ibi].z_bp_o[i]) / delti;
       
       /*      v_max=PetscMax(v_max,fabs(ibm[ibi].u[i].z)); */
       /*      v_max=PetscMax(v_max,fabs(ibm[ibi].u[i].y));      */
       /*      v_max=PetscMax(v_max,fabs(ibm[ibi].u[i].x)); */

       v_max=0.;
       i_vmax=0;
       
       if (v_max<fabs(ibm[ibi].u[i].z)) {
	 i_vmax=i;
	 v_max= ibm[ibi].u[i].z;
       }
       if (v_max<fabs(ibm[ibi].u[i].y)) {
	 i_vmax=i;
	 v_max= ibm[ibi].u[i].y;
       }
       if (v_max<fabs(ibm[ibi].u[i].x)) {
	 i_vmax=i;
	 v_max= ibm[ibi].u[i].x;
       }
     }
     PetscPrintf(PETSC_COMM_WORLD, "MAX Cop Velocity: %d %le %le %le %le\n",ti, v_max, ibm[ibi].x_bp[i_vmax],ibm[ibi].y_bp[i_vmax],ibm[ibi].z_bp[i_vmax]);
     
     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
     if (!rank) {
       if (ti == (ti/tiout)*tiout) {

	 sprintf(filen, "results/surface%3.3d_%2.2d.dat",ti,ibi);
	 f  = fopen(filen, "w");
	 PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
	 PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", ibm[ibi].n_v, ibm[ibi].n_elmt);
	 for (i=0; i<ibm[ibi].n_v; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].x_bp[i]);
	 }
	 for (i=0; i<ibm[ibi].n_v; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].y_bp[i]);
	 }
	 for (i=0; i<ibm[ibi].n_v; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].z_bp[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].nf_x[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].nf_y[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].nf_z[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].nt_x[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].nt_y[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].nt_z[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].ns_x[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].ns_y[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm[ibi].ns_z[i]);
	 }
	 for (i=0; i<ibm[ibi].n_elmt; i++) {
	   PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm[ibi].nv1[i]+1, ibm[ibi].nv2[i]+1, ibm[ibi].nv3[i]+1);
	 }

/* 	 sprintf(filen, "surface%3.3d_%2.2d.dat",ti,ibi);     */
/* 	 //sprintf(filen, "surface%3.3d.dat",ti); */
/* 	 f = fopen(filen, "w"); */
/* 	 PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n"); */
/* 	 PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", ibm[ibi].n_v, ibm[ibi].n_elmt); */
/* 	 for (i=0; i<ibm[ibi].n_v; i++) { */
/* 	   PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm[ibi].x_bp[i], ibm[ibi].y_bp[i], ibm[ibi].z_bp[i]); */
/* 	 } */
/* 	 for (i=0; i<ibm[ibi].n_elmt; i++) { */
/* 	   PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm[ibi].nv1[i]+1, ibm[ibi].nv2[i]+1, ibm[ibi].nv3[i]+1); */
/* 	 } */

	 fclose(f);
       }
     }     
   }     

/*  ! this is trying of body motion */

/*  u_bodn=u_bodnm1+cff_m*delti*F_x */
/*  v_bodn=v_bodnm1+cff_m*delti*F_y */
/*  w_bodn=w_bodnm1+cff_m*delti*F_z */

/*  write(3,*) 'u_bodn,dx=delti*u_bodn',u_bodn,delti*u_bodn */
/*  write(*,*) 'u_bodn,dx=delti*u_bodn',u_bodn,delti*u_bodn */

/*  do l=1,lv_max! */
/*   xf(l)=xf(l)+delti*u_bodn */
/*   yf(l)=yf(l)+delti*v_bodn */
/*   zf(l)=zf(l)+delti*w_bodn */
/*  end do */

/*  u_bodnm1=u_bodn */
/*  v_bodnm1=v_bodn */
/*  w_bodnm1=w_bodn */

   return(0);
}

