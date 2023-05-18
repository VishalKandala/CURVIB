static char help[] = "Interface Searching!";

#include "petscdmda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "stdlib.h"

typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct node {
  PetscInt Node;
  struct node *next;
} node;

typedef struct list{
  node *head;
} List;

typedef struct {
  DM	da, fda;
  PetscInt IM, JM, KM;
  Vec	Coor, Cent, BCS, Nvert;
  /* bctype is used to specify boundary conditionsim
     if bctype == 0 then the whole side is interface points */
  PetscInt bctype[6];
  PetscInt itfcptsnumber;
  PetscInt *itfcI, *itfcJ, *itfcK;
  PetscInt *itfchostI, *itfchostJ, *itfchostK, *itfchostB;
  PetscReal *itfchostx, *itfchosty, *itfchostz;
  DMDALocalInfo info;
  PetscReal Max_X,Max_Y,Max_Z,Min_X,Min_Y,Min_Z;
  PetscReal ddx, ddy, ddz;
  
} UserCtx;

typedef struct {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;	// number of vertices and number of elements
  PetscInt	*nv1, *nv2, *nv3;	// Node index of each triangle
  PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
  PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes

  // for radius check
  Cmpnts *qvec;
  PetscReal *radvec; 
} IBMNodes;


typedef struct {
  PetscInt *is, *ie, *js, *je, *ks, *ke;
} SearchRange;

struct list_node {
  PetscInt i, j, k, bi;
  struct list_node *next;
};

typedef struct list_node llnode;

llnode *list_add(llnode **p, int i, int j, int k, int bi) {
  llnode *n = malloc(sizeof(llnode));
  n->next = *p;
  *p = n;
  n->i = i;
  n->j = j;
  n->k = k;
  n->bi = bi;
  return n;  
}

PetscErrorCode distance(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
{
  PetscReal xn1, yn1, zn1;
  PetscReal xc, yc, zc;
  
  PetscReal dx1, dy1, dz1, dx2, dy2, dz2, r;

  dx1 = p3.x - p1.x;
  dy1 = p3.y - p1.y;
  dz1 = p3.z - p1.z;

  dx2 = p4.x - p2.x;
  dy2 = p4.y - p2.y;
  dz2 = p4.z - p2.z;

  xn1 = dy1 * dz2 - dz1 * dy2;
  yn1 = - (dx1 * dz2 - dz1 * dx2);
  zn1 = dx1 * dy2 - dy1 * dx2;

  r = sqrt(xn1 * xn1 + yn1 * yn1 + zn1 * zn1);
  xn1 /= r; yn1 /= r; zn1 /= r;

  xc = 0.25 * (p1.x + p2.x + p3.x + p4.x);
  yc = 0.25 * (p1.y + p2.y + p3.y + p4.y);
  zc = 0.25 * (p1.z + p2.z + p3.z + p4.z);

  *d = (p.x - xc) * xn1 + (p.y - yc) * yn1 + (p.z - zc) * zn1;
  if (PetscAbsReal(*d)<1.e-6) *d=0.;
  return (0);
}

PetscBool ISInsideCell(Cmpnts p, Cmpnts cell[8], PetscReal d[6])
{
  // k direction
  distance(cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  if (d[4]<0) return(PETSC_FALSE);
  distance(cell[4], cell[7], cell[6], cell[5], p, &(d[5]));
  if (d[5]<0) return(PETSC_FALSE);

  // j direction
  distance(cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  if (d[2]<0) return(PETSC_FALSE);

  distance(cell[3], cell[2], cell[6], cell[7], p, &(d[3]));
  if (d[3]<0) return(PETSC_FALSE);

  // i direction
  distance(cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  if (d[0]<0) return(PETSC_FALSE);
  
  distance(cell[1], cell[5], cell[6], cell[2], p, &(d[1]));
  if (d[1]<0) return(PETSC_FALSE);
  return(PETSC_TRUE);
}

PetscErrorCode randomdirection(Cmpnts p, PetscInt ip, PetscInt jp,
			       PetscReal xbp_min, PetscReal ybp_min,
			       PetscReal zbp_max, PetscReal dcx, PetscReal dcy,
			       PetscReal dir[3],PetscInt seed)
{
  Cmpnts endpoint;
  PetscReal s;

  PetscReal xpc, ypc; 
  
  xpc = dcx * (ip+0.5) + xbp_min;
  ypc = dcy * (jp+0.5) + ybp_min;
    
  // init rand()
  //  srand(time(NULL)+seed);
  srand(seed);
  // Generate a random number [-0.5, 0.5)
  s = rand() / ((double)RAND_MAX + 1) - 0.5;
  endpoint.x = xpc + s * dcx;
  endpoint.y = ypc + s * dcy;
  endpoint.z = zbp_max + 0.2;

  dir[0] = endpoint.x - p.x;
  dir[1] = endpoint.y - p.y;
  dir[2] = endpoint.z - p.z;

  s = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  dir[0] /= s;
  dir[1] /= s;
  dir[2] /= s;
  return 0;
}

PetscInt point_cell_advanced(Cmpnts p, PetscInt ip, PetscInt jp, PetscInt kp,
			     IBMNodes *ibm, PetscInt ncx, PetscInt ncy,
			     PetscInt ncz, PetscReal dcx, PetscReal dcy,
			     PetscReal xbp_min, PetscReal ybp_min,
			     PetscReal zbp_max, List *cell_trg,
			     PetscInt flg)
{
  PetscInt	*nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3 = ibm->nv3;
  //  PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt	i, j, k, ln_v, n1e, n2e, n3e, nintp;
  
  PetscInt	nvert_l;
  PetscReal	dt[1000], ndotn=0., dirdotn;
  //  Cmpnts        dnn[1000],nn;

  PetscReal	epsilon = 1.e-8;
  PetscReal     eps_tangent=1.e-10;

  PetscBool	*Element_Searched;
  j = jp; i = ip;

  PetscBool NotDecided = PETSC_TRUE, Singularity = PETSC_FALSE;
  PetscReal t, u, v;
  PetscReal orig[3], dir[3], vert0[3], vert1[3], vert2[3];

  node *current;
  PetscInt searchtimes=0;
  PetscMalloc(ibm->n_elmt*sizeof(PetscBool), &Element_Searched);
  if (flg) 
    PetscPrintf(PETSC_COMM_SELF, " serch itr\n");
  
  while (NotDecided) {
    
    searchtimes++;
    nintp = 0 ;
    randomdirection(p, ip, jp, xbp_min, ybp_min, zbp_max, dcx, dcy, dir, searchtimes);
    Singularity = PETSC_FALSE;
    if (flg) 
      PetscPrintf(PETSC_COMM_SELF, " serch itr, dir %d %le %le %le\n", searchtimes,dir[0],dir[1],dir[2]);
    
    for (ln_v=0; ln_v<ibm->n_elmt; ln_v++) {
      Element_Searched[ln_v] = PETSC_FALSE;
    }
    
    for (k=kp; k<ncz; k++) {
      current = cell_trg[k*ncx*ncy+j*ncx+i].head;
      while (current) {
	ln_v = current->Node;
	if (!Element_Searched[ln_v]) {
	  Element_Searched[ln_v] = PETSC_TRUE;
	  n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];
	 
	  orig[0] = p.x; orig[1] = p.y, orig[2] = p.z;

	  vert0[0] = x_bp[n1e]; vert0[1] = y_bp[n1e]; vert0[2] = z_bp[n1e];
	  vert1[0] = x_bp[n2e]; vert1[1] = y_bp[n2e]; vert1[2] = z_bp[n2e];
	  vert2[0] = x_bp[n3e]; vert2[1] = y_bp[n3e]; vert2[2] = z_bp[n3e];
            
	  //	  dirdotn=dir[0]*nn.x+dir[1]*nn.y+dir[2]*nn.z;

	  nvert_l = intsect_triangle(orig, dir, vert0, vert1, vert2, &t, &u, &v);
	  
	  if (flg) 
	    PetscPrintf(PETSC_COMM_SELF, "elm, %d %d %le %le %le %d %d %d %le\n",ln_v,nvert_l,t,u,v,n1e,n2e,n3e,dirdotn);
	  
	  if (nvert_l > 0 && t>0) {
	    dt[nintp] = t;
	    //    dnn[nintp].x=nn.x;dnn[nintp].y=nn.y;dnn[nintp].z=nn.z;

	    nintp ++;
	    PetscInt temp;
	    for (temp = 0; temp < nintp-1; temp++) {
	      // Two interception points are the same, this leads to huge
	      // trouble for crossing number test
	      // Rather to program for all cases, we use a new line to
	      // repeat the test
	      //	      ndotn=dnn[temp].x*nn.x+dnn[temp].y*nn.y+dnn[temp].z*nn.z;	      
	      
	      if ((fabs(t-dt[temp]) < epsilon && ndotn>-0.97)){ 
		
		Singularity = PETSC_TRUE;
	      }
	    }
	    if (Singularity) break;
	  }
	}
	if (Singularity) {
	  break;
	}
	else {
	  current = current->next;
	}
	} // Search through the list
      if (Singularity) {
	break;
      }
    } // for k
    if (flg) 
      PetscPrintf(PETSC_COMM_SELF, " serch itr, %d %le \n",nintp,dirdotn);
    
    if (!Singularity) {
      NotDecided = PETSC_TRUE;
      if (nintp%2) { // The interception point number is odd, inside body
	PetscFree(Element_Searched);
	return 4;
      }
      else {
	PetscFree(Element_Searched);
	return 0;
      }
    }
  }
  PetscFree(Element_Searched);
  return 0;
}

PetscErrorCode blank_search_advanced(UserCtx *user, IBMNodes *ibm, 
				     PetscInt ibi)

/*      Note : Always go from ibi (immersed body number) 0 -> NumberOfBodies  */
/*             Nvert should be set to zero before any new search and this */
/*             happens if ibi==0--not anymore! set nvert=0 manually!*/
{
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt	ncx = 40, ncy = 40, ncz = 40;
  List          *cell_trg;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt 	ln_v, n_v = ibm->n_v;

  PetscInt	i, j, k;

  PetscReal	dcx, dcy, dcz;
  PetscInt	n1e, n2e, n3e;
  PetscReal	xv_min, yv_min, zv_min, xv_max, yv_max, zv_max;
  PetscInt	iv_min, iv_max, jv_min, jv_max, kv_min, kv_max;
  PetscReal	***nvert;
  PetscInt	ic, jc, kc;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;

  for(i=0; i<n_v; i++) {
    
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);
    
    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);

    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }

  xbp_min -= 0.01; xbp_max += 0.01;
  ybp_min -= 0.01; ybp_max += 0.01;
  zbp_min -= 0.01; zbp_max += 0.01;
  
  dcx = (xbp_max - xbp_min) / (ncx - 1.);
  dcy = (ybp_max - ybp_min) / (ncy - 1.);
  dcz = (zbp_max - zbp_min) / (ncz - 1.);

  PetscMalloc(ncz * ncy * ncx * sizeof(List), &cell_trg);
 
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);
      }
    }
  }
  
  for (ln_v=0; ln_v < ibm->n_elmt; ln_v++) {

    n1e = ibm->nv1[ln_v]; n2e = ibm->nv2[ln_v]; n3e = ibm->nv3[ln_v];
    
    xv_min = PetscMin(PetscMin(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);
    xv_max = PetscMax(PetscMax(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);

    yv_min = PetscMin(PetscMin(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    yv_max = PetscMax(PetscMax(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    
    zv_min = PetscMin(PetscMin(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    zv_max = PetscMax(PetscMax(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    
    iv_min = floor((xv_min - xbp_min) / dcx); //  +1???
    iv_max = floor((xv_max - xbp_min) / dcx) +1;

    jv_min = floor((yv_min - ybp_min) / dcy); //  +1???
    jv_max = floor((yv_max - ybp_min) / dcy) +1;

    kv_min = floor((zv_min - zbp_min) / dcz); //  +1???
    kv_max = floor((zv_max - zbp_min) / dcz) +1;

    iv_min = (iv_min<0) ? 0:iv_min;
    iv_max = (iv_max>ncx) ? ncx:iv_max;

    jv_min = (jv_min<0) ? 0:jv_min;
    jv_max = (jv_max>ncx) ? ncy:jv_max;

    kv_min = (kv_min<0) ? 0:kv_min;
    kv_max = (kv_max>ncz) ? ncz:kv_max;

    // Insert IBM node information into a list
    for (k=kv_min; k<kv_max; k++) {
      for (j=jv_min; j<jv_max; j++) {
	for (i=iv_min; i<iv_max; i++) {
	  insertnode(&(cell_trg[k *ncx*ncy + j*ncx +i]), ln_v);
	}
      }
    }
  }

  PetscInt rank, flg=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  Cmpnts ***coor;
  DMDAVecGetArray(fda, user->Cent, &coor);
  DMDAVecGetArray(da, user->Nvert, &nvert);

  PetscBarrier(PETSC_NULL);
  // for this body nvert 4 is inside, 2 is near bndry
  // for previous bodies nvert sb*10 inside
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
	    coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	    coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	  ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	  jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	  kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	  
	  nvert[k][j][i] =  PetscMax(nvert[k][j][i],
				     point_cell_advanced(coor[k][j][i], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg));
	  if (nvert[k][j][i] < 0) nvert[k][j][i] = 0;
	}
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->Cent, &coor);

  // Back to the old nvert 3 and 1 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 4) nvert[k][j][i]=ibi*10;
      }
    }
  }

  PetscInt ip, im, jp, jm, kp, km;
  PetscInt ii, jj, kk;
  // Near Boundary?
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ip = (i<mx-2?(i+1):(i));
	im = (i>1   ?(i-1):(i));
	
	jp = (j<my-2?(j+1):(j));
	jm = (j>1   ?(j-1):(j));
	
	kp = (k<mz-2?(k+1):(k));
	km = (k>1   ?(k-1):(k));
	
	if ((int)(nvert[k][j][i]+0.5) == ibi*10) {
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(nvert[kk][jj][ii] +0.5) == 0) {
		  nvert[k][j][i] = PetscMin(ibi*10-1, nvert[k][j][i]);
		  
		}
	      }
	    }
	  }
	}
      }
    }
  }


/*   for (k=lzs; k<lze; k++) { */
/*     for (j=lys; j<lye; j++) { */
/*       for (i=lxs; i<lxe; i++) { */
/* 	ip = (i<mx-2?(i+1):(i)); */
/* 	im = (i>1   ?(i-1):(i)); */
	
/* 	jp = (j<my-2?(j+1):(j)); */
/* 	jm = (j>1   ?(j-1):(j)); */
	
/* 	kp = (k<mz-2?(k+1):(k)); */
/* 	km = (k>1   ?(k-1):(k)); */
	
/* 	if ((int)(nvert[k][j][i]+0.5) == ibi*10) { */
/* 	  for (kk=km; kk<kp+1; kk++) { */
/* 	    for (jj=jm; jj<jp+1; jj++) { */
/* 	      for (ii=im; ii<ip+1; ii++) { */
/* 		if ((int)(nvert[kk][jj][ii] +0.5) == ibi*10-1) { */
/* 		  nvert[k][j][i] = PetscMin(ibi*10-2, nvert[k][j][i]); */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */

  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	destroy(&cell_trg[k*ncx*ncy+j*ncx+i]);
      }
    }
  }

  PetscFree(cell_trg);
  
  PetscViewer	viewer;
  char filen[80];
  sprintf(filen, "nvfield_blank.dat");
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);

  return 0;
}

PetscErrorCode CenterNodeCoor(UserCtx *user)
{

  PetscInt       i, j, k, IM, JM, KM;
  Cmpnts ***coor, ***cent;

  DMDAVecGetArray(user->fda, user->Coor, &coor);
  DMDAVecGetArray(user->fda, user->Cent, &cent);

  
  for (k=1; k<user->KM; k++) {
    for (j=1; j<user->JM; j++) {
      for (i=1; i<user->IM; i++) {
	cent[k][j][i].x = 0.125 * (coor[k  ][j  ][i  ].x +
				   coor[k  ][j  ][i-1].x +
				   coor[k  ][j-1][i  ].x +
				   coor[k  ][j-1][i-1].x +
				   coor[k-1][j  ][i  ].x +
				   coor[k-1][j  ][i-1].x +
				   coor[k-1][j-1][i  ].x +
				   coor[k-1][j-1][i-1].x);
	
	cent[k][j][i].y = 0.125 * (coor[k  ][j  ][i  ].y +
				   coor[k  ][j  ][i-1].y +
				   coor[k  ][j-1][i  ].y +
				   coor[k  ][j-1][i-1].y +
				   coor[k-1][j  ][i  ].y +
				   coor[k-1][j  ][i-1].y +
				   coor[k-1][j-1][i  ].y +
				   coor[k-1][j-1][i-1].y);
	
	cent[k][j][i].z = 0.125 * (coor[k  ][j  ][i  ].z +
				   coor[k  ][j  ][i-1].z +
				   coor[k  ][j-1][i  ].z +
				   coor[k  ][j-1][i-1].z +
				   coor[k-1][j  ][i  ].z +
				   coor[k-1][j  ][i-1].z +
				   coor[k-1][j-1][i  ].z +
				   coor[k-1][j-1][i-1].z);
	
      }
    }
  }
  
  // Ghost nodes
  k=0;
  for (j=1; j<user->JM; j++) {
    for (i=1; i<user->IM; i++) {
      cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i  ].x +
				   coor[k  ][j  ][i-1].x +
				   coor[k  ][j-1][i  ].x +
				   coor[k  ][j-1][i-1].x)-
	cent[k+1][j][i].x;
      
      cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i  ].y +
				   coor[k  ][j  ][i-1].y +
				   coor[k  ][j-1][i  ].y +
				   coor[k  ][j-1][i-1].y )-
	cent[k+1][j][i].y;
      
      cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i  ].z +
				   coor[k  ][j  ][i-1].z +
				   coor[k  ][j-1][i  ].z +
				   coor[k  ][j-1][i-1].z)-
	cent[k+1][j][i].z;      
    }
  }
  
  k=user->KM;
  for (j=1; j<user->JM; j++) {
    for (i=1; i<user->IM; i++) {
      cent[k][j][i].x = 2.*0.25 * (coor[k-1][j  ][i  ].x +
				   coor[k-1][j  ][i-1].x +
				   coor[k-1][j-1][i  ].x +
				   coor[k-1][j-1][i-1].x)-
	cent[k-1][j][i].x;
      
      cent[k][j][i].y = 2.*0.25 * (coor[k-1][j  ][i  ].y +
				   coor[k-1][j  ][i-1].y +
				   coor[k-1][j-1][i  ].y +
				   coor[k-1][j-1][i-1].y )-
	cent[k-1][j][i].y;
      
      cent[k][j][i].z = 2.*0.25 * (coor[k-1][j  ][i  ].z +
				   coor[k-1][j  ][i-1].z +
				   coor[k-1][j-1][i  ].z +
				   coor[k-1][j-1][i-1].z)-
	cent[k-1][j][i].z;      
    }
  }

  j=0;
  for (k=1; k<user->KM; k++) {
    for (i=1; i<user->IM; i++) {
      cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i  ].x +
				   coor[k  ][j  ][i-1].x +
				   coor[k-1][j  ][i  ].x +
				   coor[k-1][j  ][i-1].x)-
	cent[k][j+1][i].x;      
	
      cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i  ].y +
				   coor[k  ][j  ][i-1].y +
				   coor[k-1][j  ][i  ].y +
				   coor[k-1][j  ][i-1].y)-
	cent[k][j+1][i].y;      
      
      cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i  ].z +
				    coor[k  ][j  ][i-1].z +
				    coor[k-1][j  ][i  ].z +
				    coor[k-1][j  ][i-1].z) -
	cent[k][j+1][i].z;                  
    }
  }

  j=user->JM;
  for (k=1; k<user->KM; k++) {
    for (i=1; i<user->IM; i++) {
      cent[k][j][i].x = 2.*0.25 * (coor[k  ][j-1][i  ].x +
				   coor[k  ][j-1][i-1].x +
				   coor[k-1][j-1][i  ].x +
				   coor[k-1][j-1][i-1].x)-
	cent[k][j-1][i].x;      
	
      cent[k][j][i].y = 2.*0.25 * (coor[k  ][j-1][i  ].y +
				   coor[k  ][j-1][i-1].y +
				   coor[k-1][j-1][i  ].y +
				   coor[k-1][j-1][i-1].y)-
	cent[k][j-1][i].y;      
      
      cent[k][j][i].z = 2.*0.25 * (coor[k  ][j-1][i  ].z +
				   coor[k  ][j-1][i-1].z +
				   coor[k-1][j-1][i  ].z +
				   coor[k-1][j-1][i-1].z) -
	cent[k][j-1][i].z;                  
    }
  }

  i=0;
  for (k=1; k<user->KM; k++) {
    for (j=1; j<user->JM; j++) {
	cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i  ].x +
				     coor[k  ][j-1][i  ].x +
				     coor[k-1][j  ][i  ].x +
				     coor[k-1][j-1][i  ].x)-
	  cent[k][j][i+1].x;  
	
	cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i  ].y +
				     coor[k  ][j-1][i  ].y +
				     coor[k-1][j  ][i  ].y +
				     coor[k-1][j-1][i  ].y)-
	  cent[k][j][i+1].y;  

		
	cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i  ].z +
				     coor[k  ][j-1][i  ].z +
				     coor[k-1][j  ][i  ].z +
				     coor[k-1][j-1][i  ].z)-
	  cent[k][j][i+1].z;  	
    }
  }

  i=user->IM;
  for (k=1; k<user->KM; k++) {
    for (j=1; j<user->JM; j++) {
	cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i-1].x +
				     coor[k  ][j-1][i-1].x +
				     coor[k-1][j  ][i-1].x +
				     coor[k-1][j-1][i-1].x)-
	  cent[k][j][i-1].x;  
	
	cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i-1].y +
				     coor[k  ][j-1][i-1].y +
				     coor[k-1][j  ][i-1].y +
				     coor[k-1][j-1][i-1].y)-
	  cent[k][j][i-1].y;  
		
	cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i-1].z +
				     coor[k  ][j-1][i-1].z +
				     coor[k-1][j  ][i-1].z +
				     coor[k-1][j-1][i-1].z)-
	  cent[k][j][i-1].z;  	
    }
  }

  // lines 

  // kmin plane lines
  k=0;

  j=0;
  for (i=1; i<user->IM; i++) {
    cent[k][j][i].x = 2.*coor[k][j][i].x - cent[k+1][j+1][i].x;
    cent[k][j][i].y = 2.*coor[k][j][i].y - cent[k+1][j+1][i].y;
    cent[k][j][i].z = 2.*coor[k][j][i].z - cent[k+1][j+1][i].z;
  }

  j=user->JM;
  for (i=1; i<user->IM; i++) {
    cent[k][j][i].x = 2.*coor[k][j-1][i].x - cent[k+1][j-1][i].x;
    cent[k][j][i].y = 2.*coor[k][j-1][i].y - cent[k+1][j-1][i].y;
    cent[k][j][i].z = 2.*coor[k][j-1][i].z - cent[k+1][j-1][i].z;
  }

  i=0;
  for (j=1; j<user->JM; j++) {
    cent[k][j][i].x = 2.*coor[k][j][i].x - cent[k+1][j][i+1].x;
    cent[k][j][i].y = 2.*coor[k][j][i].y - cent[k+1][j][i+1].y;
    cent[k][j][i].z = 2.*coor[k][j][i].z - cent[k+1][j][i+1].z;
  }

  i=user->IM;
  for (j=1; j<user->JM; j++) {
    cent[k][j][i].x = 2.*coor[k][j][i-1].x - cent[k+1][j][i-1].x;
    cent[k][j][i].y = 2.*coor[k][j][i-1].y - cent[k+1][j][i-1].y;
    cent[k][j][i].z = 2.*coor[k][j][i-1].z - cent[k+1][j][i-1].z;
  }

  // kmax plane lines
  k=user->KM;

  j=0;
  for (i=1; i<user->IM; i++) {
    cent[k][j][i].x = 2.*coor[k-1][j][i].x - cent[k-1][j+1][i].x;
    cent[k][j][i].y = 2.*coor[k-1][j][i].y - cent[k-1][j+1][i].y;
    cent[k][j][i].z = 2.*coor[k-1][j][i].z - cent[k-1][j+1][i].z;
  }

  j=user->JM;
  for (i=1; i<user->IM; i++) {
    cent[k][j][i].x = 2.*coor[k-1][j-1][i].x - cent[k-1][j-1][i].x;
    cent[k][j][i].y = 2.*coor[k-1][j-1][i].y - cent[k-1][j-1][i].y;
    cent[k][j][i].z = 2.*coor[k-1][j-1][i].z - cent[k-1][j-1][i].z;
  }

  i=0;
  for (j=1; j<user->JM; j++) {
    cent[k][j][i].x = 2.*coor[k-1][j][i].x - cent[k-1][j][i+1].x;
    cent[k][j][i].y = 2.*coor[k-1][j][i].y - cent[k-1][j][i+1].y;
    cent[k][j][i].z = 2.*coor[k-1][j][i].z - cent[k-1][j][i+1].z;
  }

  i=user->IM;
  for (j=1; j<user->JM; j++) {
    cent[k][j][i].x = 2.*coor[k-1][j][i-1].x - cent[k-1][j][i-1].x;
    cent[k][j][i].y = 2.*coor[k-1][j][i-1].y - cent[k-1][j][i-1].y;
    cent[k][j][i].z = 2.*coor[k-1][j][i-1].z - cent[k-1][j][i-1].z;
  }

  // jmin lines
  j=0;

  i=0;
  for (k=1; k<user->KM; k++) {
    cent[k][j][i].x = 2.*coor[k][j][i].x - cent[k][j+1][i+1].x;
    cent[k][j][i].y = 2.*coor[k][j][i].y - cent[k][j+1][i+1].y;
    cent[k][j][i].z = 2.*coor[k][j][i].z - cent[k][j+1][i+1].z;
  }

  i=user->IM;
  for (k=1; k<user->KM; k++) {
    cent[k][j][i].x = 2.*coor[k][j][i-1].x - cent[k][j+1][i-1].x;
    cent[k][j][i].y = 2.*coor[k][j][i-1].y - cent[k][j+1][i-1].y;
    cent[k][j][i].z = 2.*coor[k][j][i-1].z - cent[k][j+1][i-1].z;	

  }

  // jmax lines
  j=user->JM;

  i=0;
  for (k=1; k<user->KM; k++) {
    cent[k][j][i].x = 2.*coor[k][j-1][i].x - cent[k][j-1][i+1].x;
    cent[k][j][i].y = 2.*coor[k][j-1][i].y - cent[k][j-1][i+1].y;
    cent[k][j][i].z = 2.*coor[k][j-1][i].z - cent[k][j-1][i+1].z;
  }

  i=user->IM;
  for (k=1; k<user->KM; k++) {
    cent[k][j][i].x = 2.*coor[k][j-1][i-1].x - cent[k][j-1][i-1].x;
    cent[k][j][i].y = 2.*coor[k][j-1][i-1].y - cent[k][j-1][i-1].y;
    cent[k][j][i].z = 2.*coor[k][j-1][i-1].z - cent[k][j-1][i-1].z;
  }

  // Corners
  i=0;j=0;k=0;
  cent[k][j][i].x = 2.*coor[k][j][i].x - cent[k+1][j+1][i+1].x;
  cent[k][j][i].y = 2.*coor[k][j][i].y - cent[k+1][j+1][i+1].y;
  cent[k][j][i].z = 2.*coor[k][j][i].z - cent[k+1][j+1][i+1].z;

  i=0;j=0;k=user->KM;
  cent[k][j][i].x = 2.*coor[k-1][j][i].x - cent[k-1][j+1][i+1].x;
  cent[k][j][i].y = 2.*coor[k-1][j][i].y - cent[k-1][j+1][i+1].y;
  cent[k][j][i].z = 2.*coor[k-1][j][i].z - cent[k-1][j+1][i+1].z;

  i=0;j=user->JM;k=0;
  cent[k][j][i].x = 2.*coor[k][j-1][i].x - cent[k+1][j-1][i+1].x;
  cent[k][j][i].y = 2.*coor[k][j-1][i].y - cent[k+1][j-1][i+1].y;
  cent[k][j][i].z = 2.*coor[k][j-1][i].z - cent[k+1][j-1][i+1].z;

  i=0;j=user->JM;k=user->KM;
  cent[k][j][i].x = 2.*coor[k-1][j-1][i].x - cent[k-1][j-1][i+1].x;
  cent[k][j][i].y = 2.*coor[k-1][j-1][i].y - cent[k-1][j-1][i+1].y;
  cent[k][j][i].z = 2.*coor[k-1][j-1][i].z - cent[k-1][j-1][i+1].z;

  i=user->IM;j=0;k=0;
  cent[k][j][i].x = 2.*coor[k][j][i-1].x - cent[k+1][j+1][i-1].x;
  cent[k][j][i].y = 2.*coor[k][j][i-1].y - cent[k+1][j+1][i-1].y;
  cent[k][j][i].z = 2.*coor[k][j][i-1].z - cent[k+1][j+1][i-1].z;

  i=user->IM;j=0;k=user->KM;
  cent[k][j][i].x = 2.*coor[k-1][j][i-1].x - cent[k-1][j+1][i-1].x;
  cent[k][j][i].y = 2.*coor[k-1][j][i-1].y - cent[k-1][j+1][i-1].y;
  cent[k][j][i].z = 2.*coor[k-1][j][i-1].z - cent[k-1][j+1][i-1].z;

  i=user->IM;j=user->JM;k=0;
  cent[k][j][i].x = 2.*coor[k][j-1][i-1].x - cent[k+1][j-1][i-1].x;
  cent[k][j][i].y = 2.*coor[k][j-1][i-1].y - cent[k+1][j-1][i-1].y;
  cent[k][j][i].z = 2.*coor[k][j-1][i-1].z - cent[k+1][j-1][i-1].z;

  i=user->IM;j=user->JM;k=user->KM;
  cent[k][j][i].x = 2.*coor[k-1][j-1][i-1].x - cent[k-1][j-1][i-1].x;
  cent[k][j][i].y = 2.*coor[k-1][j-1][i-1].y - cent[k-1][j-1][i-1].y;
  cent[k][j][i].z = 2.*coor[k-1][j-1][i-1].z - cent[k-1][j-1][i-1].z;

  DMDAVecRestoreArray(user->fda, user->Cent, &cent);
  DMDAVecRestoreArray(user->fda, user->Coor, &coor);

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
  //  cl=1./L_dim;
  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, PETSC_NULL);     


  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ ibmdata\n");
    char filen[80];  
    sprintf(filen,"blank_grid%2.2d.dat",ibi);
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file for Blanking");
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
      
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	
	x_bp[i] = x_bp[i]/cl;// + CMx_c;//0.25 ;// 24.;	
	y_bp[i] = y_bp[i]/cl;// + CMy_c ;//+ ibi*2.;//5.;//2.;//8.;//6.;//2.   ;// 24.;
	z_bp[i] = z_bp[i]/cl;// + CMz_c ;//+ ibi*2.;//2.;//8.;//15.;//2.   ;// 24.;
	
	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];

      }
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);

      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
      
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, &ss, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

      fclose(fd);
    }
      
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    FILE *f;
    sprintf(filen, "surface_trigrid.dat");
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL)\n", ibm->n_v, ibm->n_elmt);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<ibm->n_v; i++) {	
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
    }
    for (i=0; i<ibm->n_elmt; i++) {
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
    
    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;
  }
  return(0);
}

PetscErrorCode CreateTriGrid(UserCtx *user,PetscInt ld[6], IBMNodes *ibm)
{
  PetscInt    i, j, k, n_v=0, n_elmt=0,n_vp[6];
  PetscInt    is,ie, js,je,ks,ke; 
  Cmpnts      ***coor;

  is=ld[0];ie=ld[1]+1;
  js=ld[2];je=ld[3]+1;
  ks=ld[4];ke=ld[5]+1;
  
  DMDAVecGetArray(user->fda, user->Coor, &coor);
  ibm->n_v=2*(ie-is-2)*(je-js-2)+2*(ie-is-2)*(ke-ks)+2*(ke-ks)*(je-js);
  n_vp[0]=  (ke-ks)*(je-js); //imin
  n_vp[1]=2*(ke-ks)*(je-js);//imin+imax
  n_vp[2]=2*(ke-ks)*(je-js)+  (ie-is-2)*(ke-ks);//imin+imax+jmin
  n_vp[3]=2*(ke-ks)*(je-js)+2*(ie-is-2)*(ke-ks);//imin+imax+jmin+jmax
  n_vp[4]=2*(ke-ks)*(je-js)+2*(ie-is-2)*(ke-ks)+(ie-is-2)*(je-js-2);//imin+imax+jmin+jamx+kmin

  ibm->n_elmt=4*(ie-is-1)*(je-js-1)+4*(ie-is-1)*(ke-ks-1)+4*(ke-ks-1)*(je-js-1);
  PetscMalloc(ibm->n_v*sizeof(PetscReal), &(ibm->x_bp));
  PetscMalloc(ibm->n_v*sizeof(PetscReal), &(ibm->y_bp));
  PetscMalloc(ibm->n_v*sizeof(PetscReal), &(ibm->z_bp));
  
  //i_min plane  
  i=is;
  for (k=ks;k<ke;k++) {
    for (j=js;j<je;j++) {
      ibm->x_bp[n_v]=coor[k][j][i].x;
      ibm->y_bp[n_v]=coor[k][j][i].y;
      ibm->z_bp[n_v]=coor[k][j][i].z;      
      n_v++;
    }
  }

  //i_max plane  
  i=ie-1;
  for (k=ks;k<ke;k++) {
    for (j=js;j<je;j++) {
      ibm->x_bp[n_v]=coor[k][j][i].x;
      ibm->y_bp[n_v]=coor[k][j][i].y;
      ibm->z_bp[n_v]=coor[k][j][i].z;      
      n_v++;
    }
  }

  //j_min plane  
  j=js;
  for (k=ks;k<ke;k++) {
    for (i=is+1;i<ie-1;i++) {
      ibm->x_bp[n_v]=coor[k][j][i].x;
      ibm->y_bp[n_v]=coor[k][j][i].y;
      ibm->z_bp[n_v]=coor[k][j][i].z;      
      n_v++;
    }
  }

  //j_max plane  
  j=je-1;
  for (k=ks;k<ke;k++) {
    for (i=is+1;i<ie-1;i++) {
      ibm->x_bp[n_v]=coor[k][j][i].x;
      ibm->y_bp[n_v]=coor[k][j][i].y;
      ibm->z_bp[n_v]=coor[k][j][i].z;      
      n_v++;
    }
  }

  //k_min plane  
  k=ks;
  for (j=js+1;j<je-1;j++) {
    for (i=is+1;i<ie-1;i++) {
      ibm->x_bp[n_v]=coor[k][j][i].x;
      ibm->y_bp[n_v]=coor[k][j][i].y;
      ibm->z_bp[n_v]=coor[k][j][i].z;      
      n_v++;
    }
  }

  //k_max plane  
  k=ke-1;
  for (j=js+1;j<je-1;j++) {
    for (i=is+1;i<ie-1;i++) {
      ibm->x_bp[n_v]=coor[k][j][i].x;
      ibm->y_bp[n_v]=coor[k][j][i].y;
      ibm->z_bp[n_v]=coor[k][j][i].z;      
      n_v++;
    }
  }

  DMDAVecRestoreArray(user->fda, user->Coor, &coor);

  // Elements
  PetscMalloc(ibm->n_elmt*sizeof(PetscInt), &(ibm->nv1));
  PetscMalloc(ibm->n_elmt*sizeof(PetscInt), &(ibm->nv2));
  PetscMalloc(ibm->n_elmt*sizeof(PetscInt), &(ibm->nv3));

  //i_min plane  
  i=is;
  for (k=ks;k<ke-1;k++) {
    for (j=js;j<je-1;j++) {
      ibm->nv1[n_elmt]=(j  -js)+(je-js)*(k  -ks);
      ibm->nv2[n_elmt]=(j+1-js)+(je-js)*(k  -ks);
      ibm->nv3[n_elmt]=(j  -js)+(je-js)*(k+1-ks);
      n_elmt++;
      ibm->nv1[n_elmt]=(j+1-js)+(je-js)*(k+1-ks);
      ibm->nv2[n_elmt]=(j+1-js)+(je-js)*(k  -ks);
      ibm->nv3[n_elmt]=(j  -js)+(je-js)*(k+1-ks);
      n_elmt++; 
    }
  }

  //i_max plane  
  i=ie-1;
  for (k=ks;k<ke-1;k++) {
    for (j=js;j<je-1;j++) {
      ibm->nv1[n_elmt]=(j  -js)+(je-js)*(k  -ks)+n_vp[0];
      ibm->nv2[n_elmt]=(j+1-js)+(je-js)*(k  -ks)+n_vp[0];
      ibm->nv3[n_elmt]=(j  -js)+(je-js)*(k+1-ks)+n_vp[0];
      n_elmt++;
      ibm->nv1[n_elmt]=(j+1-js)+(je-js)*(k+1-ks)+n_vp[0];
      ibm->nv2[n_elmt]=(j+1-js)+(je-js)*(k  -ks)+n_vp[0];
      ibm->nv3[n_elmt]=(j  -js)+(je-js)*(k+1-ks)+n_vp[0];
      n_elmt++;
    }
  }

  //j_min plane  
  j=js;
  for (k=ks;k<ke-1;k++) {
    for (i=is;i<ie-1;i++) {
      // i k
      if (i==is)        ibm->nv1[n_elmt]=(j  -js  )+(je-js  )*(k  -ks);
      else      	ibm->nv1[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks)+n_vp[1];
      // i+1 k
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks)+n_vp[0];
      else	        ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks)+n_vp[1];
      // i  k+1
      if (i==is) 	ibm->nv3[n_elmt]=(j  -js  )+(je-js  )*(k+1-ks);
      else      	ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k+1-ks)+n_vp[1];
      n_elmt++;

      // i+1, k+1
      if (i==ie-2)	ibm->nv1[n_elmt]=(j  -js  )+(je-js  )*(k+1-ks)+n_vp[0];
      else      	ibm->nv1[n_elmt]=(i+1-is-1)+(ie-is-2)*(k+1-ks)+n_vp[1];
      // i+1 k
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks)+n_vp[0];
      else      	ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks)+n_vp[1];
      // i k+1     
      if (i==is)  	ibm->nv3[n_elmt]=(j  -js  )+(je-js  )*(k+1-ks);
      else      	ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k+1-ks)+n_vp[1];
      n_elmt++; 
    }
  }

  //j_max plane  
  j=je-1;
  for (k=ks;k<ke-1;k++) {
    for (i=is;i<ie-1;i++) {
      // i k
      if (i==is)	ibm->nv1[n_elmt]=(j  -js  )+(je-js  )*(k  -ks);
      else              ibm->nv1[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks)+n_vp[2];
      // i+1 k
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks)+n_vp[0];
      else	        ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks)+n_vp[2];
      // i   k+1
      if (i==is) 	ibm->nv3[n_elmt]=(j  -js  )+(je-js  )*(k+1-ks);
      else              ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k+1-ks)+n_vp[2];
      n_elmt++;

      // i+1 k+1
      if (i==ie-2)	ibm->nv1[n_elmt]=(j  -js  )+(je-js  )*(k+1-ks)+n_vp[0];
      else              ibm->nv1[n_elmt]=(i+1-is-1)+(ie-is-2)*(k+1-ks)+n_vp[2];
      // i+1 k
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks)+n_vp[0];
      else              ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks)+n_vp[2];
      // i  k+1
      if (i==is)  	ibm->nv3[n_elmt]=(j  -js  )+(je-js  )*(k+1-ks);
      else              ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k+1-ks)+n_vp[2];
      n_elmt++; 
    }
  }

  //k_min plane  
  k=ks;
  for (j=js;j<je-1;j++) {
    for (i=is;i<ie-1;i++) {
      // i j 
      if (i==is)	ibm->nv1[n_elmt]=(j  -js  )+(je-js  )*(k  -ks  );
      else if (j==js)   ibm->nv1[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks  )+n_vp[1];
      else              ibm->nv1[n_elmt]=(i  -is-1)+(ie-is-2)*(j  -js-1)+n_vp[3];
      // i+1 j
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks  )+n_vp[0];
      else if (j==js)   ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks  )+n_vp[1];
      else              ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(j  -js-1)+n_vp[3];
      // i  j+1
      if (i==is) 	ibm->nv3[n_elmt]=(j+1-js  )+(je-js  )*(k  -ks  );
      else if (j==je-2) ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks  )+n_vp[2];
      else              ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(j+1-js-1)+n_vp[3];
      n_elmt++;
      // i+1 j+1
      if (i==ie-2)	ibm->nv1[n_elmt]=(j+1-js  )+(je-js  )*(k  -ks  )+n_vp[0];
      else if (j==je-2) ibm->nv1[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks  )+n_vp[2];      
      else              ibm->nv1[n_elmt]=(i+1-is-1)+(ie-is-2)*(j+1-js-1)+n_vp[3];
      // i+1 j
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks  )+n_vp[0];
      else if (j==js)   ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks  )+n_vp[1];
      else              ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(j  -js-1)+n_vp[3];
      // i  j+1
      if (i==is)        ibm->nv3[n_elmt]=(j+1-js  )+(je-js  )*(k  -ks  );
      else if (j==je-2) ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks  )+n_vp[2];
      else              ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(j+1-js-1)+n_vp[3];
      n_elmt++; 
    }
  }

  //k_max plane  
  k=ke-1;
  for (j=js;j<je-1;j++) {
    for (i=is;i<ie-1;i++) {
      // i j 
      if (i==is)	ibm->nv1[n_elmt]=(j  -js  )+(je-js  )*(k  -ks  );
      else if (j==js)   ibm->nv1[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks  )+n_vp[1];
      else              ibm->nv1[n_elmt]=(i  -is-1)+(ie-is-2)*(j  -js-1)+n_vp[4];
      // i+1 j
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks  )+n_vp[0];
      else if (j==js)   ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks  )+n_vp[1];
      else              ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(j  -js-1)+n_vp[4];
      // i  j+1
      if (i==is) 	ibm->nv3[n_elmt]=(j+1-js  )+(je-js  )*(k  -ks  );
      else if (j==je-2) ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks  )+n_vp[2];
      else              ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(j+1-js-1)+n_vp[4];
      n_elmt++;
      // i+1 j+1
      if (i==ie-2)	ibm->nv1[n_elmt]=(j+1-js  )+(je-js  )*(k  -ks  )+n_vp[0];
      else if (j==je-2) ibm->nv1[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks  )+n_vp[2];      
      else              ibm->nv1[n_elmt]=(i+1-is-1)+(ie-is-2)*(j+1-js-1)+n_vp[4];
      // i+1 j
      if (i==ie-2)	ibm->nv2[n_elmt]=(j  -js  )+(je-js  )*(k  -ks  )+n_vp[0];
      else if (j==js)   ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(k  -ks  )+n_vp[1];
      else              ibm->nv2[n_elmt]=(i+1-is-1)+(ie-is-2)*(j  -js-1)+n_vp[4];
      // i  j+1
      if (i==is)        ibm->nv3[n_elmt]=(j+1-js  )+(je-js  )*(k  -ks  );
      else if (j==je-2) ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(k  -ks  )+n_vp[2];
      else              ibm->nv3[n_elmt]=(i  -is-1)+(ie-is-2)*(j+1-js-1)+n_vp[4];
      n_elmt++; 
    }
  }
  
  PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d %d %d\n",ibm->n_v, n_v,ibm->n_elmt,n_elmt);
 
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surface_trigrid.dat");
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL)\n", ibm->n_v, ibm->n_elmt);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]);
    }
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]);
    }
    for (i=0; i<ibm->n_v; i++) {	
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]);
    }
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
    }
    fclose(f);    
  }
  return(0);
}

int main(int argc, char**argv)
{
  PetscInt block_number;
  UserCtx *user;

  PetscInt      bi, i, j, k,nn,indx, IM, JM, KM;
  PetscReal	***bcs, ***nvert;
  PetscReal	cl = 1.;
  //  Vec           Coor;
  PetscReal     *gc;

  PetscReal  Max_X,Max_Y,Max_Z;
  PetscReal  Min_X,Min_Y,Min_Z;
  PetscReal  x,y,z;
  PetscInt   ip[10],jp[10],kp[10],sb;
  PetscInt   ip1,jp1,kp1, kp0;
  PetscReal  xpmin,xpmax,ypmin,ypmax,zpmin,zpmax;
  PetscInt   dispx[10],dispy[10],dispz[10],dispz0;

  PetscInt dIM=40, dJM=40, dKM=40;
  PetscInt dI, dJ, dK;
  PetscReal xmin=1.e10, xmax=-1.e10, ymin=1.e10, ymax=-1.e10, 
    zmin=1.e10, zmax=-1.e10, tmp;
  PetscReal ddx, ddy, ddz;

  llnode *ptsincell[dKM][dJM][dIM], *head[dKM][dJM][dIM], *current;

  Cmpnts ***coor, ***cent;
  IBMNodes ibm;

  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscInt	blank = 0,TwoD=0, blank_surface=0;
  //  PetscOptionsInsertFile("control.dat");
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_NULL);
  
  PetscOptionsGetInt(PETSC_NULL, "-blk_grd", &blank_surface, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-blk", &blank, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-twoD", &TwoD, PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-cl", &cl, PETSC_NULL);

  if (TwoD==1) dIM=1;
  if (TwoD==2) dJM=1;

  PetscInt ids[dKM][dJM][dIM][10], ide[dKM][dJM][dIM][10],
    jds[dKM][dJM][dIM][10], jde[dKM][dJM][dIM][10],
    kds[dKM][dJM][dIM][10], kde[dKM][dJM][dIM][10];


  PetscInt generate_grid=0, grid1d=0;
  PetscOptionsGetInt(PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);

  // Read in Grid data
  FILE *fd;
  fd = fopen("grid.dat", "r");
  fscanf(fd, "%i\n", &block_number);
  PetscPrintf(PETSC_COMM_WORLD, "%i\n", block_number);
  PetscMalloc(block_number*sizeof(UserCtx), &user);
  
  for (bi=0; bi<block_number; bi++) {

    if (bi==0)  cl = 1.;
    if (bi==1)  cl = 0.85;


    fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
    IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;

    DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1, 1, PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &(user[bi].da));

    DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DMGetCoordinateDM(user[bi].da, &(user[bi].fda));
    /* Commented out BY IMAN Feb 5 2009*/

    /* End Commented out BY IMAN Feb 5 2009*/
    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i\n", user[bi].IM, user[bi].JM, user[bi].KM);

    /* ADDED BY IMAN */
    DMDAGetLocalInfo(user[bi].da, &(user[bi].info));

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    

    //DMDAGetGhostedCoordinates(user[bi].da, &user[bi].Coor);
    DMGetCoordinatesLocal(user[bi].da, &user[bi].Coor);
    DMDAVecGetArray(user[bi].fda, user[bi].Coor, &coor);

    if (grid1d) {
      PetscReal xx;

      PetscMalloc((IM+JM+KM)*sizeof(PetscReal), &gc);
      PetscPrintf(PETSC_COMM_SELF, "Malloc0\n");
      // read i
      for (i=0; i<IM; i++) 
	fscanf(fd, "%le %le %le\n",&gc[i],&xx,&xx);
      // read j
      for (j=0; j<JM; j++) 
	fscanf(fd, "%le %le %le\n",&xx,&gc[IM+j],&xx);
      // read k
      for (i=0; i<KM; i++) 
	fscanf(fd, "%le %le %le\n",&xx,&xx,&gc[IM+JM+i]);

      MPI_Bcast(gc, (IM+JM+KM), MPIU_REAL, 0, PETSC_COMM_WORLD);

      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if (k<KM && j<JM && i<IM) {
	      coor[k][j][i].x = *(gc + i)/cl;
	      coor[k][j][i].y = *(gc + IM + j)/cl;
	      coor[k][j][i].z = *(gc + IM + JM + k)/cl;
	    }
	  }
	}
      }
       
   } else { // if 3d gridgen file
      PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);
      for (k=0; k<KM; k++) {
	for (j=0; j<JM; j++) {
	  for (i=0; i<IM; i++) {
	    fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
/* 	    *(gc+(k*JM*IM + j*IM + i)*3) = 1./(IM-1.) * i;	 */
	  }
	}
      }
      
      for (k=0; k<KM; k++) {
	for (j=0; j<JM; j++) {
	  for (i=0; i<IM; i++) {
	    fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
/* 	    *(gc+(k*JM*IM + j*IM + i)*3+1) = 1./(JM-1.) * j; */
	  }
	}
      }

      for (k=0; k<KM; k++) {
	for (j=0; j<JM; j++) {
	  for (i=0; i<IM; i++) {
	    fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
/* 	    *(gc+(k*JM*IM + j*IM + i)*3+2) = 1./(KM-1.) * k; */
	  }
	}
      }
        
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if (k<KM && j<JM && i<IM) {
	      coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl;
	      coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl;
	      coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl;
	    }
	  }
	}
      }
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if ((i==0 || i==mx-2) && (j==1 || j==my-2) && (k==1 || k==mz-2))
	      PetscPrintf(PETSC_COMM_SELF, "@i=%d j=%d k=%d coor.y is %le coor.z is %le \n",i,j,k,coor[k][j][i].y,coor[k][j][i].z);
	  }
	}
      }
    }
    PetscFree(gc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Coor, &coor);


    Vec	gCoor;
    //DMDAGetCoordinates(user[bi].da, &gCoor);
    DMGetCoordinates(user[bi].da, &gCoor);
    DMLocalToGlobalBegin(user[bi].fda, user[bi].Coor, INSERT_VALUES, gCoor);
    DMLocalToGlobalEnd(user[bi].fda, user[bi].Coor, INSERT_VALUES, gCoor);

    DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, user[bi].Coor);
    DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, user[bi].Coor);

    VecDestroy(&gCoor);
    
    // read grid ends
    xmin=1.e10; xmax=-1.e10; ymin=1.e10; ymax=-1.e10; 
    zmin=1.e10; zmax=-1.e10;
/*     /\* Find the Min and Max of grid *\/ */
    
    VecStrideMin(user[bi].Coor, 0, PETSC_NULL, &tmp);
    xmin = PetscMin(xmin, tmp);
    VecStrideMax(user[bi].Coor, 0, PETSC_NULL, &tmp);
    xmax = PetscMax(xmax, tmp);

    VecStrideMin(user[bi].Coor, 1, PETSC_NULL, &tmp);
    ymin = PetscMin(ymin, tmp);
    VecStrideMax(user[bi].Coor, 1, PETSC_NULL, &tmp);
    ymax = PetscMax(ymax, tmp);

    VecStrideMin(user[bi].Coor, 2, PETSC_NULL, &tmp);
    zmin = PetscMin(zmin, tmp);
    VecStrideMax(user[bi].Coor, 2, PETSC_NULL, &tmp);
    zmax = PetscMax(zmax, tmp);

    PetscPrintf(PETSC_COMM_WORLD, "%i xmin %le %le ymin %le %le zmin %le %le\n",bi,xmin,xmax, ymin,ymax, zmin,zmax);

  } //bi

  fclose(fd);

  /**************************************************************************************************************************/
  /* Read boundary conditions  */
  /**************************************************************************************************************************/

  fd = fopen("bcs.dat", "r");
  for (bi=0; bi<block_number; bi++) {
    fscanf(fd, "%i %i %i %i %i %i\n", &(user[bi].bctype[0]),
	   &(user[bi].bctype[1]), &(user[bi].bctype[2]), &(user[bi].bctype[3]),
	   &(user[bi].bctype[4]), &(user[bi].bctype[5]));
    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i %i %i %i\n",bi,user[bi].bctype[0],user[bi].bctype[1],user[bi].bctype[2],user[bi].bctype[3],user[bi].bctype[4],user[bi].bctype[5] );
  }
  fclose(fd);

  /**************************************************************************************************************************/  
  /* Calculate the locations of grid center and
     find Min Max of grid
  */
  /**************************************************************************************************************************/
  PetscReal  xc[block_number],yc[block_number],zc[block_number];

  for (bi=0; bi<block_number; bi++) {
    VecDuplicate(user[bi].Coor, &(user[bi].Cent));
    VecSet(user[bi].Cent, 0.);
    CenterNodeCoor(&user[bi]);

    // Find Min Max of grid centers
    VecSetBlockSize(user[bi].Cent, 3);
    xmin=1.e10; xmax=-1.e10; ymin=1.e10; ymax=-1.e10; 
    zmin=1.e10; zmax=-1.e10;
    /* Find the Min and Max of grid */
    VecStrideMin(user[bi].Cent, 0, PETSC_NULL, &tmp);
    xmin = PetscMin(xmin, tmp);
    VecStrideMax(user[bi].Cent, 0, PETSC_NULL, &tmp);
    xmax = PetscMax(xmax, tmp);

    VecStrideMin(user[bi].Cent, 1, PETSC_NULL, &tmp);
    ymin = PetscMin(ymin, tmp);
    VecStrideMax(user[bi].Cent, 1, PETSC_NULL, &tmp);
    ymax = PetscMax(ymax, tmp);

    VecStrideMin(user[bi].Cent, 2, &indx, &tmp);
    zmin = PetscMin(zmin, tmp);
    VecStrideMax(user[bi].Cent, 2, PETSC_NULL, &tmp);
    zmax = PetscMax(zmax, tmp);

    xmin -= 0.01; xmax += 0.01;
    ymin -= 0.01; ymax += 0.01;
    zmin -= 0.01; zmax += 0.01;

    user[bi].Max_X = xmax;
    user[bi].Max_Y = ymax;
    user[bi].Max_Z = zmax;
    user[bi].Min_X = xmin;
    user[bi].Min_Y = ymin;
    user[bi].Min_Z = zmin;

    // find the center of inner domain
    xc[bi]=(user[bi].Min_X+user[bi].Max_X)/2.0;
    yc[bi]=(user[bi].Min_Y+user[bi].Max_Y)/2.0;
    zc[bi]=(user[bi].Min_Z+user[bi].Max_Z)/2.0;         

    VecGetBlockSize(user[bi].Cent, &nn);
    PetscPrintf(PETSC_COMM_WORLD, "%i %i xmin %le %le ymin %le %le zmin %le %le %i\n",nn,bi,xmin,xmax, ymin,ymax, zmin,zmax,indx);
    PetscPrintf(PETSC_COMM_WORLD, "%i %i xc %le %le %le\n",nn,bi,xc[bi], yc[bi], zc[bi]);
  } //bi

  /**************************************************************************************************************************/
  /* Partition the grid  */
  /**************************************************************************************************************************/
  for (bi=0; bi<block_number; bi++) {	 
    for (k=0; k<dKM; k++) {
      for (j=0; j<dJM; j++) {
	for (i=0; i<dIM; i++) {
	  ids[k][j][i][bi] = 1000;
	  ide[k][j][i][bi] = 0;
	  
	  jds[k][j][i][bi] = 1000;
	  jde[k][j][i][bi] = 0;
	  
	  kds[k][j][i][bi] = 1000;
	  kde[k][j][i][bi] = 0;
	}
      }
    }
  }//bi loop

  /**************************************************************************************************************************/  
  /* 
     locate the center points into a dIM * dJM * dKM cell system
  */
  /**************************************************************************************************************************/
    
  for (bi=0; bi<block_number; bi++) {
    DMDAVecGetArray(user[bi].fda, user[bi].Cent, &coor);

    ddx = (user[bi].Max_X-user[bi].Min_X)/(double)(dIM);
    ddy = (user[bi].Max_Y-user[bi].Min_Y)/(double)(dJM);
    ddz = (user[bi].Max_Z-user[bi].Min_Z)/(double)(dKM);

    user[bi].ddx =ddx;
    user[bi].ddy =ddy;
    user[bi].ddz =ddz;

    PetscPrintf(PETSC_COMM_WORLD, "bi %i ddxyz %le %le %le\n",bi,ddx,ddy,ddz);
    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
	for (i=0; i<user[bi].IM; i++) {
	  dI = floor((coor[k][j][i].x - user[bi].Min_X) / ddx);
	  dJ = floor((coor[k][j][i].y - user[bi].Min_Y) / ddy);
	  dK = floor((coor[k][j][i].z - user[bi].Min_Z) / ddz);
	  //head[dK][dJ][dI] = list_add(&(head[dK][dJ][dI]), i, j, k, bi);

	  if (i<=ids[dK][dJ][dI][bi]) ids[dK][dJ][dI][bi]=i;
	  if (i>=ide[dK][dJ][dI][bi]) ide[dK][dJ][dI][bi]=i;
	  
	  if (j<=jds[dK][dJ][dI][bi]) jds[dK][dJ][dI][bi]=j;
	  if (j>=jde[dK][dJ][dI][bi]) jde[dK][dJ][dI][bi]=j;
	  
	  if (k<=kds[dK][dJ][dI][bi]) kds[dK][dJ][dI][bi]=k;
	  if (k>=kde[dK][dJ][dI][bi]) kde[dK][dJ][dI][bi]=k;	  
	}
      }  
    }
    DMDAVecRestoreArray(user[bi].fda, user[bi].Cent, &coor);
    nn=4;
    for (dK=0; dK<dKM; dK++) {
      for (dJ=0; dJ<dJM; dJ++) {
	for (dI=0; dI<dIM; dI++) {
	  if(ids[dK][dJ][dI][bi]>nn && ids[dK][dJ][dI][bi]<999) 
	    ids[dK][dJ][dI][bi]-= nn;
	  else
	    ids[dK][dJ][dI][bi]=0;
	  if(ide[dK][dJ][dI][bi]<user[bi].IM-nn && ide[dK][dJ][dI][bi]>0)
	    ide[dK][dJ][dI][bi]+= nn;
	  else 
	    ide[dK][dJ][dI][bi]=user[bi].IM;
	  if(jds[dK][dJ][dI][bi]> nn && jds[dK][dJ][dI][bi]<999)
	    jds[dK][dJ][dI][bi]-= nn;
	  else
	    jds[dK][dJ][dI][bi]=0;
	  if(jde[dK][dJ][dI][bi]<user[bi].JM-nn && jde[dK][dJ][dI][bi]>0)
	    jde[dK][dJ][dI][bi]+=nn;
	  else 
	    jde[dK][dJ][dI][bi]=user[bi].JM;
	  if(kds[dK][dJ][dI][bi]>nn && kds[dK][dJ][dI][bi]<999)
	    kds[dK][dJ][dI][bi] -= nn;
	  else
	    kds[dK][dJ][dI][bi]=0;
	  if(kde[dK][dJ][dI][bi]<user[bi].KM-nn && kde[dK][dJ][dI][bi]>0)
	    kde[dK][dJ][dI][bi]+=nn;
	  else
	    kde[dK][dJ][dI][bi]=user[bi].KM;
	}
      }
    }  
  } //bi

  /**************************************************************************************************************************/
  /* Decide wheter a point is interface point */
  /**************************************************************************************************************************/
  PetscInt    is,ie, js,je,ks,ke, ld[6];

  for (bi = 0; bi<block_number; bi++) {
    DMCreateGlobalVector(user[bi].da, &(user[bi].BCS));
    VecSet(user[bi].BCS, 1.);
    //    DMDAVecGetArray(user[bi].fda, user[bi].Cent, &cent);
    user[bi].itfcptsnumber=0;

    /**************************************************************************************************************************/
    /*  inside another domain for blanking
     change Iman */
    /**************************************************************************************************************************/
    if (blank) {
      if (bi==0) {
    	VecDuplicate(user[bi].BCS, &(user[bi].Nvert));
	for (sb=0; sb<block_number; sb++) {
	  if (sb!=bi) {
	    // find the is, ie, js, je, ks, ke of grid sb inside which 
	    // the background mesh will be blanked
	    DMDALocalInfo info = user[sb].info;
	    PetscInt    mx = info.mx, my = info.my, mz = info.mz;	  
	    
	    is=mx/4;//20;// 
	    ie=mx-2-mx/4;//121;//
	    js=my/4;//15;//101-71;//
	    je=my-2-my/4;//126;//101+71;//
	    ks=mz/4;//20; //101-71;//35;//
	    ke=mz-2-mz/4;//105;//101+71;//
	    PetscOptionsGetInt(PETSC_NULL, "-is", &is, PETSC_NULL);
	    PetscOptionsGetInt(PETSC_NULL, "-ie", &ie, PETSC_NULL);
	    PetscOptionsGetInt(PETSC_NULL, "-js", &js, PETSC_NULL);
	    PetscOptionsGetInt(PETSC_NULL, "-je", &je, PETSC_NULL);
	    PetscOptionsGetInt(PETSC_NULL, "-ks", &ks, PETSC_NULL);
	    PetscOptionsGetInt(PETSC_NULL, "-ke", &ke, PETSC_NULL);
	    
	    PetscPrintf(PETSC_COMM_WORLD, "Interface blak is-e %d %d  js-e %d %d  ks-e %d %d\n", is,ie,js,je,ks,ke);
	    
	    ld[0]=is;ld[1]=ie;
	    ld[2]=js;ld[3]=je;
	    ld[4]=ks;ld[5]=ke;
	    
	    // create a triangular grid for sb
	    if (blank_surface) 
	      ibm_read_ucd(&ibm, sb);
	    else
	      CreateTriGrid(&user[sb],ld, &ibm);
	    
	  // find the background grid nodes of block bi inside 
	  // the triangular grid for sb
	    blank_search_advanced(&user[bi], &ibm, sb);

	  // Find the buffer nodes that (3 layers) that will be
	  // interpolated
	    DMDAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
	    DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
	    
	    info = user[bi].info;
	    mx = info.mx, my = info.my, mz = info.mz;
	    PetscInt    xs = info.xs, xe = info.xs + info.xm;
	    PetscInt    ys = info.ys, ye = info.ys + info.ym;
	    PetscInt    zs = info.zs, ze = info.zs + info.zm;
	    
	    PetscInt ip, im, jp, jm, kp, km;
	    PetscInt ii, jj, kk;
	    
	    for (k=zs; k<ze; k++) {
	      for (j=ys; j<ye; j++) {
		for (i=xs; i<xe; i++) {
		  ip = (i<mx-1?(i+1):(i));
		  im = (i>0   ?(i-1):(i));
		  
		  jp = (j<my-1?(j+1):(j));
		  jm = (j>0   ?(j-1):(j));
		  
		  kp = (k<mz-1?(k+1):(k));
		  km = (k>0   ?(k-1):(k));
		  
		  if (((int)(nvert[k][j][i]+0.1) < sb*10) && 
		      ((int)(nvert[k][j][i]+0.1) > sb*10-3) ) {
		    for (kk=km; kk<kp; kk++) {
		      for (jj=jm; jj<jp; jj++) {
			for (ii=im; ii<ip; ii++) {
			  bcs[kk][jj][ii] = PetscMin(-sb*10, bcs[kk][jj][ii]);
			}
		      }
		    }
		  }
		}
	      }
	    }
	    
	    DMDAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
	    DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); 
	  } // if sb
	}// for sb
	VecDestroy(&user[bi].Nvert);
      } //if bi==0
    } //if blank

    /**************************************************************************************************************************/
    /*  Check if outer boundary are interfaces
     */
    /**************************************************************************************************************************/    
    DMDAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    if (user[bi].bctype[0] == 0) { // imin boundary
      for (k=0; k<user[bi].KM; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  //	  for (i=0; i<2; i++) {
	  for (i=0; i<1; i++) {
	    bcs[k][j][i] = -0.1;
	    // user[bi].itfcptsnumber++;
	  }
	}
      }
    }

    if (user[bi].bctype[1] == 0) { // imax boundary
      for (k=0; k<user[bi].KM; k++) {
	for (j=0; j<user[bi].JM; j++) {
/* 	  for (i=user[bi].IM-2; i<user[bi].IM; i++) { */
	  for (i=user[bi].IM-1; i<user[bi].IM; i++) {
	    bcs[k][j][i] = -0.1;
	    // user[bi].itfcptsnumber++;
	  }
	}
      }
    }

    if (user[bi].bctype[2] == 0) { // jmin boundary
      for (k=0; k<user[bi].KM; k++) {
	//	for (j=0; j<2; j++) {
	for (j=0; j<1; j++) {
	  for (i=0; i<user[bi].IM; i++) {
	    bcs[k][j][i] = -1.;
	    // user[bi].itfcptsnumber++;
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0) { // jmax boundary
      for (k=0; k<user[bi].KM; k++) {
	//	for (j=user[bi].JM-2; j<user[bi].JM; j++) {
	for (j=user[bi].JM-1; j<user[bi].JM; j++) {
	  for (i=0; i<user[bi].IM; i++) {
	    bcs[k][j][i] = -1.;
	    // user[bi].itfcptsnumber++;
	  }
	}
      }
    }

    if (user[bi].bctype[4] == 0) { // kmin boundary
      //      for (k=0; k<2; k++) {
      for (k=0; k<1; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  for (i=0; i<user[bi].IM; i++) {
	    bcs[k][j][i] = -2.;
	    //user[bi].itfcptsnumber++;
	  }
	}
      }
    }

    if (user[bi].bctype[5] == 0) { // kmax boundary
      //      for (k=user[bi].KM-2; k<user[bi].KM; k++) {
      for (k=user[bi].KM-1; k<user[bi].KM; k++) {
	for (j=0; j<user[bi].JM; j++) {
	  for (i=0; i<user[bi].IM; i++) {
	    bcs[k][j][i] = -2.;
	    //user[bi].itfcptsnumber++;
	  }
	}
      }
    }

    DMDAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
  }//for bi
  /**************************************************************************************************************************/
  /*   Get the center node coordinates for host nodes */
  /**************************************************************************************************************************/
  Cmpnts ***host[10];
  for (bi=0; bi<block_number; bi++) {
    //DMDAVecGetArray(user[bi].fda, user[bi].Coor, &(host[bi]));
      DMDAVecGetArray(user[bi].fda, user[bi].Cent, &(host[bi]));// change Iman
  }

  /**************************************************************************************************************************/
  /* Find the number of interface points
    i.e. itfcptsnumber and allocate 
    memory */
  /**************************************************************************************************************************/

  PetscReal d[6];
  Cmpnts pc;
  PetscBool found;
  PetscReal  xh[8], yh[8], zh[8];
  PetscInt bh, si, sj, sk;
  Cmpnts cell[8];

  for (bi=0; bi<block_number; bi++) {
    DMDAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    user[bi].itfcptsnumber=0;
    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
        for (i=0; i<user[bi].IM; i++) {
          if (bcs[k][j][i]<1.e-6) {
            user[bi].itfcptsnumber++;
          }
        }
      }
    }
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfcI);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfcJ);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfcK);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostI);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostJ);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostK);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscInt), &user[bi].itfchostB);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscReal), &user[bi].itfchostx);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscReal), &user[bi].itfchosty);
    PetscMalloc(user[bi].itfcptsnumber*sizeof(PetscReal), &user[bi].itfchostz);
    DMDAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
  }

  /**************************************************************************************************************************/
  /* Find host cells for the interface nodes  */
  /**************************************************************************************************************************/
  PetscInt number,nNotFound;//, is, ie, js, je, ks, ke;
  PetscInt dKs,dKe,dJs,dJe,dIs,dIe;

  for (bi=0; bi<block_number; bi++) {
    DMDAVecGetArray(user[bi].fda, user[bi].Coor, &coor);
   
    DMDAVecGetArray(user[bi].da, user[bi].BCS, &bcs);
    number = -1;
    nNotFound =0 ;
    PetscPrintf(PETSC_COMM_WORLD, "%i %i %i %i\n", user[bi].IM, user[bi].JM, user[bi].KM, bi);
    for (k=0; k<user[bi].KM; k++) {
      for (j=0; j<user[bi].JM; j++) {
	for (i=0; i<user[bi].IM; i++) {
	  if(bcs[k][j][i]<-1.e-6) {
	 
	    number++;
	    	    pc.x = coor[k][j][i].x;
	    	    pc.y = coor[k][j][i].y;
	    	    pc.z = coor[k][j][i].z;

	 	    found = PETSC_FALSE;
	   	    for (sb=0; sb<block_number; sb++) {
		      if (sb!=bi) {
			if (bi>-1 && sb>-1) {
		
			  ddx = (user[sb].Max_X-user[sb].Min_X)/(double)(dIM);
			  ddy = (user[sb].Max_Y-user[sb].Min_Y)/(double)(dJM);
			  ddz = (user[sb].Max_Z-user[sb].Min_Z)/(double)(dKM);
		  
			  // if in bounding box
			  if (pc.x >= user[sb].Min_X && pc.x <= user[sb].Max_X &&
			      pc.y >= user[sb].Min_Y && pc.y <= user[sb].Max_Y &&
			      pc.z >= user[sb].Min_Z && pc.z <= user[sb].Max_Z){
			    dI = floor((pc.x - user[sb].Min_X) / ddx);
			    dJ = floor((pc.y - user[sb].Min_Y) / ddy);
			    dK = floor((pc.z - user[sb].Min_Z) / ddz);

	
			    for (sk=kds[dK][dJ][dI][sb]; sk<kde[dK][dJ][dI][sb]; sk++) {
			      for (sj=jds[dK][dJ][dI][sb]; sj<jde[dK][dJ][dI][sb]; sj++) {
				for (si=ids[dK][dJ][dI][sb]; si<ide[dK][dJ][dI][sb]; si++) {
				  //		      PetscPrintf(PETSC_COMM_WORLD, "sk %i %i %i\n", sk, sj, si);
				  cell[0] = host[sb][sk  ][sj  ][si  ];
				  cell[1] = host[sb][sk  ][sj  ][si+1];
				  cell[2] = host[sb][sk  ][sj+1][si+1];
				  cell[3] = host[sb][sk  ][sj+1][si  ];
				  
				  cell[4] = host[sb][sk+1][sj  ][si  ];
				  cell[5] = host[sb][sk+1][sj  ][si+1];
				  cell[6] = host[sb][sk+1][sj+1][si+1];
				  cell[7] = host[sb][sk+1][sj+1][si  ];
				  
				  if(ISInsideCell(pc, cell, d)) {
				    found = PETSC_TRUE;
				    user[bi].itfcI[number] = i;
				    user[bi].itfcJ[number] = j;
				    user[bi].itfcK[number] = k;
				    user[bi].itfchostI[number] = si;
				    user[bi].itfchostJ[number] = sj;
				    user[bi].itfchostK[number] = sk;
				    user[bi].itfchostB[number] = sb;
				    user[bi].itfchostx[number] = d[0] / (d[0] + d[1]);
				    user[bi].itfchosty[number] = d[2] / (d[2] + d[3]);
				    user[bi].itfchostz[number] = d[4] / (d[4] + d[5]);
				    
				    goto nextp;
				  }
				}
			      }
			    }
			  }// if in bounding box
			  else
			    PetscPrintf(PETSC_COMM_WORLD, "not in bounding box %d %d %d bi %d %le\n", i,j,k,bi,bcs[k][j][i]);
			} // if bi>sb
		      } //if sb
		    } //sb
	  nextp: number=number;
		    if (!found) {
		      nNotFound++;
		      user[bi].itfcI[number] = i;
		      user[bi].itfcJ[number] = j;
		      user[bi].itfcK[number] = k;
		      user[bi].itfchostI[number] = 0;
		      user[bi].itfchostJ[number] = 0;
		      user[bi].itfchostK[number] = 0;
		      user[bi].itfchostB[number] = sb;
		      user[bi].itfchostx[number] = -2.;//d[0] / (d[0] + d[1]);
		      user[bi].itfchosty[number] = 0.;//d[2] / (d[2] + d[3]);
		      user[bi].itfchostz[number] = 0.;//d[4] / (d[4] + d[5]);
		    }	   	  
	  }	  
        }
      }
    }
    
    
    PetscPrintf(PETSC_COMM_WORLD, "%i number %i itfc number %i  not found %i\n", bi, number+1,user[bi].itfcptsnumber, nNotFound);
    
    DMDAVecRestoreArray(user[bi].da, user[bi].BCS, &bcs);
   
    DMDAVecRestoreArray(user[bi].fda, user[bi].Coor, &coor);
    
  } //bi

  fd = fopen("interface.dat", "w");
  for (bi=0; bi<block_number; bi++) {
    if (blank && bi==0) {     
      for (sb=1; sb<block_number; sb++) {
	PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %i %i %i %i %i\n", is,ie,js,je,ks,ke,nn);
      }
    }    
    PetscFPrintf(PETSC_COMM_WORLD, fd, "%i\n", user[bi].itfcptsnumber);
    PetscPrintf(PETSC_COMM_WORLD, "bi=%i  ifcp#=%i\n", bi,user[bi].itfcptsnumber);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %i\n", user[bi].itfcI[i], user[bi].itfcJ[i], user[bi].itfcK[i]);
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%i %i %i %i\n", user[bi].itfchostI[i], user[bi].itfchostJ[i], user[bi].itfchostK[i], user[bi].itfchostB[i]);
      PetscFPrintf(PETSC_COMM_WORLD, fd, "%le %le %le\n", user[bi].itfchostx[i], user[bi].itfchosty[i], user[bi].itfchostz[i]);
    }
  }
  fclose(fd);
  PetscFinalize();
  return(0);
}

