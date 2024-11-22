#include "variables.h"
#include "stdlib.h"
#include "time.h"
extern PetscInt thin, block_number,invicid;
extern PetscInt NumberOfBodies, moveframe,ti;

#define Dist(p1, p2) sqrt(( p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z))

#define Cross(Resu, v1, v2) \
        Resu.x = v1.y * v2.z - v1.z * v2.y; \
        Resu.y = v1.z * v2.x - v1.x * v2.z; \
        Resu.z = v1.x * v2.y - v1.y * v2.x;

#define VecAMinusB(C, A, B) \
        C.x = A.x - B.x; \
        C.y = A.y - B.y; \
        C.z = A.z - B.z;
 

PetscInt flagprint = 0;

PetscReal detmnt(PetscReal a[3][3])
{
  PetscReal tmp;
  tmp = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
	 a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
	 a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
  return tmp;
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

PetscErrorCode ibm_search_advanced(UserCtx *user, IBMNodes *ibm, 
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
 

  dcx = (xbp_max - xbp_min) / (ncx - 1.);//double (ncx) //
  dcy = (ybp_max - ybp_min) / (ncy - 1.);//double (ncy) //
  dcz = (zbp_max - zbp_min) / (ncz - 1.);//double (ncz) //


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

 //Mohsen 

  PetscInt rank, flg=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscPrintf(PETSC_COMM_SELF, "test002\n");
  Cmpnts ***coor;
  DMDAVecGetArray(fda, user->Cent, &coor);
  DMDAVecGetArray(da, user->Nvert, &nvert);

  
  // for this body nvert 4 is inside, 2 is near bndry
  // for previous bodies nvert 3 inside, 1 near bndry
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
	    coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	    coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	  ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	  jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	  kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	  
	  /* search only if the node is not inside another body 
	     already! */

	  nvert[k][j][i] =  PetscMax(nvert[k][j][i],
				     point_cell_advanced(coor[k][j][i], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg));
	  if (nvert[k][j][i] < 0) nvert[k][j][i] = 0;
	}
      }
    }
  }
  
  
  DMDAVecRestoreArray(da, user->Nvert, &nvert);

  DMDAVecGetArray(da, user->Nvert, &nvert);

  if (thin) {
  PetscPrintf(PETSC_COMM_WORLD, "IBM thin  %d %d %le %le %le %le %le %le\n", ibm->n_v, ibm->n_elmt, xbp_max, xbp_min, ybp_max, ybp_min, zbp_max, zbp_min);
    PetscInt cutthrough;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
	      coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	      coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	    ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	    jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	    kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	    
	    cutthrough = point_cell_thin(coor[k][j][i],coor[k][j][i+1],coor[k][j+1][i],coor[k+1][j][i],coor[k+1][j+1][i+1], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg);

	    if (cutthrough) {
	      if (nvert[k  ][j  ][i  ] < 0.5) nvert[k  ][j  ][i  ]=2.;
	      if (nvert[k  ][j  ][i+1] < 0.5) nvert[k  ][j  ][i+1]=2.;
	      if (nvert[k  ][j+1][i  ] < 0.5) nvert[k  ][j+1][i  ]=2.;
	      if (nvert[k  ][j+1][i+1] < 0.5) nvert[k  ][j+1][i+1]=2.;

	      if (nvert[k+1][j  ][i  ] < 0.5) nvert[k+1][j  ][i  ]=2.;
	      if (nvert[k+1][j  ][i+1] < 0.5) nvert[k+1][j  ][i+1]=2.;
	      if (nvert[k+1][j+1][i  ] < 0.5) nvert[k+1][j+1][i  ]=2.;
	      if (nvert[k+1][j+1][i+1] < 0.5) nvert[k+1][j+1][i+1]=2.;
	    }
	  }
	}
      }
    }    
  }

  PetscInt ip, im, jp, jm, kp, km;
  PetscInt ii, jj, kk;

  // Near boundary?
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert[k][j][i] <0) nvert[k][j][i] = 0;
	ip = (i<mx-1?(i+1):(i));
	im = (i>0   ?(i-1):(i));

	jp = (j<my-1?(j+1):(j));
	jm = (j>0   ?(j-1):(j));

	kp = (k<mz-1?(k+1):(k));
	km = (k>0   ?(k-1):(k));

	if ((int)(nvert[k][j][i]+0.5) != 4) {
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(nvert[kk][jj][ii] +0.5) == 4) {
		  nvert[k][j][i] = PetscMax(2, nvert[k][j][i]);
		}
	      }
	    }
	  }
	}
      }
    }
  }

 

  DMDAVecRestoreArray(da, user->Nvert, &nvert);
 
  if (user->ibmlist[ibi].head) DestroyIBMList(&(user->ibmlist[ibi]));
 
  InitIBMList(&user->ibmlist[ibi]);
  PetscInt number;
  number = 0;

  IBMInfo ibm_intp;

  DMDAVecGetArray(da, user->Nvert, &nvert);

  BoundingSphere(ibm); 

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) {
	  number ++;
	  ic = (int)((coor[k][j][i].x - xbp_min) / dcx);
	  jc = (int)((coor[k][j][i].y - ybp_min) / dcy);
	  kc = (int)((coor[k][j][i].z - zbp_min) / dcz);

	  if (ic<0) ic=0;
	  else if (ic>=ncx) ic=ncx-1;

	  if (jc<0) jc=0;
	  else if (jc>=ncy) jc = ncy-1;

	  if (kc<0) kc=0;
	  else if (kc>=ncz) kc = ncz-1;

	  ibm_intp.ni = i;
	  ibm_intp.nj = j;
	  ibm_intp.nk = k;

	  nearestcell1(coor[k][j][i], ibm, &ibm_intp);
	  
	  InterceptionPoint(coor[k][j][i], i, j, k, &ibm_intp, user);
	 
	  if (ibm_intp.imode<0) {
	    PetscInt cell;
	    Cmpnts ptmp;
	    if (i==1 || i==mx-2 ||
		j==1 || j==my-2) {

	      cell = ibm_intp.cell;
	      if (ibm->nf_z[cell] > 0) {
		ptmp = coor[k+1][j][i];
		ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				    (coor[k][j][i].x - ptmp.x) +
				    (coor[k][j][i].y - ptmp.y) *
				    (coor[k][j][i].y - ptmp.y) +
				    (coor[k][j][i].z - ptmp.z) *
				    (coor[k][j][i].z - ptmp.z));
		ibm_intp.cr1 = 1.;
		ibm_intp.cr2 = 0.;
		ibm_intp.cr3 = 0.;
		ibm_intp.i1 = i;
		ibm_intp.j1 = j;
		ibm_intp.k1 = k+1;

		ibm_intp.i2 = i;
		ibm_intp.j2 = j;
		ibm_intp.k2 = k+1;
		ibm_intp.i3 = i;
		ibm_intp.j3 = j;
		ibm_intp.k3 = k+1;
	      }
	      else {
		ptmp = coor[k-1][j][i];
		ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				    (coor[k][j][i].x - ptmp.x) +
				    (coor[k][j][i].y - ptmp.y) *
				    (coor[k][j][i].y - ptmp.y) +
				    (coor[k][j][i].z - ptmp.z) *
				    (coor[k][j][i].z - ptmp.z));
		ibm_intp.cr1 = 1.;
		ibm_intp.cr2 = 0.;
		ibm_intp.cr3 = 0.;
		ibm_intp.i1 = i;
		ibm_intp.j1 = j;
		ibm_intp.k1 = k-1;

		ibm_intp.i2 = i;
		ibm_intp.j2 = j;
		ibm_intp.k2 = k-1;
		ibm_intp.i3 = i;
		ibm_intp.j3 = j;
		ibm_intp.k3 = k-1;
	      }
	    }
	    else if (k==1 || k==mz-2) {
	      cell = ibm_intp.cell;
	      ptmp = coor[k][j+1][i];
	      ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				  (coor[k][j][i].x - ptmp.x) +
				  (coor[k][j][i].y - ptmp.y) *
				  (coor[k][j][i].y - ptmp.y) +
				  (coor[k][j][i].z - ptmp.z) *
				  (coor[k][j][i].z - ptmp.z));
	      ibm_intp.cr1 = 1.;
	      ibm_intp.cr2 = 0.;
	      ibm_intp.cr3 = 0.;
	      ibm_intp.i1 = i;
	      ibm_intp.j1 = j+1;
	      ibm_intp.k1 = k;
	      
	      ibm_intp.i2 = i;
	      ibm_intp.j2 = j+1;
	      ibm_intp.k2 = k;
	      ibm_intp.i3 = i;
	      ibm_intp.j3 = j+1;
	      ibm_intp.k3 = k;
	    }
	    else {
	      PetscPrintf(PETSC_COMM_SELF, "%%%%IBM Searching Fail! %d %d %d\n", i, j, k);
	    }
	  }

	  AddIBMNode(&user->ibmlist[ibi], ibm_intp);
	}
      }
    }
  }

  PetscBarrier(PETSC_NULL);

 

  // Back to the old nvert 3 and 1 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) nvert[k][j][i]=1;
	if ((int)(nvert[k][j][i]+0.5) == 4) nvert[k][j][i]=3;
      }
    }
  }


  DMDAVecRestoreArray(fda, user->Cent,&coor);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	destroy(&cell_trg[k*ncx*ncy+j*ncx+i]);
      }
    }
  }

  PetscFree(cell_trg);
  PetscFree(ibm->qvec);
  PetscFree(ibm->radvec); 

  return 0;
}

PetscErrorCode ibm_search_advanced_rev(UserCtx *user, IBMNodes *ibm, 
				       PetscInt ibi)

/*      Note : Always go from ibi (immersed body number) 0 -> NumberOfBodies  */
/*             Nvert should be set to zero before any new search  */
/*              */
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
  DMDAVecGetArray(fda, user->lCent, &coor);
  DMDAVecGetArray(da, user->Nvert, &nvert);

  
  // for this body nvert 4 is inside, 2 is near bndry
  // for previous bodies nvert 3 inside, 1 near bndry
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
	    coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	    coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	  ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	  jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	  kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	 
	  /* search only if the node is not inside another body 
	     already! */

	  nvert[k][j][i] =  PetscMax(nvert[k][j][i],
				     point_cell_advanced(coor[k][j][i], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg));
	  nvert[k][j][i] -=4;
	  if (nvert[k][j][i] < 0) nvert[k][j][i] = 4;
	} else 
	  nvert[k][j][i] = 4;
      }
    }
  }
  
 
  DMDAVecRestoreArray(da, user->Nvert, &nvert);

  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  if (thin) {
  PetscPrintf(PETSC_COMM_WORLD, "IBM thin  %d %d %le %le %le %le %le %le\n", ibm->n_v, ibm->n_elmt, xbp_max, xbp_min, ybp_max, ybp_min, zbp_max, zbp_min);
    PetscInt cutthrough;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
	      coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	      coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	    ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	    jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	    kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	    
	    cutthrough = point_cell_thin(coor[k][j][i],coor[k][j][i+1],coor[k][j+1][i],coor[k+1][j][i],coor[k+1][j+1][i+1], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg);

	    if (cutthrough) {
	      if (nvert[k  ][j  ][i  ] < 0.5) nvert[k  ][j  ][i  ]=2.;
	      if (nvert[k  ][j  ][i+1] < 0.5) nvert[k  ][j  ][i+1]=2.;
	      if (nvert[k  ][j+1][i  ] < 0.5) nvert[k  ][j+1][i  ]=2.;
	      if (nvert[k  ][j+1][i+1] < 0.5) nvert[k  ][j+1][i+1]=2.;

	      if (nvert[k+1][j  ][i  ] < 0.5) nvert[k+1][j  ][i  ]=2.;
	      if (nvert[k+1][j  ][i+1] < 0.5) nvert[k+1][j  ][i+1]=2.;
	      if (nvert[k+1][j+1][i  ] < 0.5) nvert[k+1][j+1][i  ]=2.;
	      if (nvert[k+1][j+1][i+1] < 0.5) nvert[k+1][j+1][i+1]=2.;
	    }
	  }
	}
      }
    }    
  }

  
  PetscInt ip, im, jp, jm, kp, km;
  PetscInt ii, jj, kk;

  // Near boundary?
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert[k][j][i] <0) nvert[k][j][i] = 0;
	ip = (i<mx-1?(i+1):(i));
	im = (i>0   ?(i-1):(i));

	jp = (j<my-1?(j+1):(j));
	jm = (j>0   ?(j-1):(j));

	kp = (k<mz-1?(k+1):(k));
	km = (k>0   ?(k-1):(k));

	if ((int)(nvert[k][j][i]+0.5) != 4) {
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(nvert[kk][jj][ii] +0.5) == 4) {
		  nvert[k][j][i] = PetscMax(2, nvert[k][j][i]);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  PetscBarrier(PETSC_NULL);
 
  PetscReal	***nvert_o;
  DMDAVecGetArray(da, user->lNvert_o, &nvert_o);
  if (ibi==NumberOfBodies-1)
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert_o[k][j][i] >2.5 && nvert[k][j][i] < 0.5) {
	  PetscPrintf(PETSC_COMM_SELF, "Phase Change at %d, %d, %d!\n", i, j, k);
	  nvert[k][j][i]=2;
	}
      }
    }
  }
  DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);
  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert);
  DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert);

  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  if (user->ibmlist[ibi].head) DestroyIBMList(&(user->ibmlist[ibi]));

  InitIBMList(&user->ibmlist[ibi]);
  PetscInt number;
  number = 0;

  IBMInfo ibm_intp;

  DMDAVecGetArray(da, user->lNvert, &nvert);

  BoundingSphere(ibm); 

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) {
	  number ++;
	  ic = (int)((coor[k][j][i].x - xbp_min) / dcx);
	  jc = (int)((coor[k][j][i].y - ybp_min) / dcy);
	  kc = (int)((coor[k][j][i].z - zbp_min) / dcz);

	  if (ic<0) ic=0;
	  else if (ic>=ncx) ic=ncx-1;

	  if (jc<0) jc=0;
	  else if (jc>=ncy) jc = ncy-1;

	  if (kc<0) kc=0;
	  else if (kc>=ncz) kc = ncz-1;

	  ibm_intp.ni = i;
	  ibm_intp.nj = j;
	  ibm_intp.nk = k;

	  nearestcell1(coor[k][j][i], ibm, &ibm_intp);
	  
	  InterceptionPoint(coor[k][j][i], i, j, k, &ibm_intp, user);

	  if (ibm_intp.imode<0) {
	    PetscInt cell;
	    Cmpnts ptmp;
	    if (i==1 || i==mx-2 ||
		j==1 || j==my-2) {

	      cell = ibm_intp.cell;
	      if (ibm->nf_z[cell] > 0) {
		ptmp = coor[k+1][j][i];
		ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				    (coor[k][j][i].x - ptmp.x) +
				    (coor[k][j][i].y - ptmp.y) *
				    (coor[k][j][i].y - ptmp.y) +
				    (coor[k][j][i].z - ptmp.z) *
				    (coor[k][j][i].z - ptmp.z));
		ibm_intp.cr1 = 1.;
		ibm_intp.cr2 = 0.;
		ibm_intp.cr3 = 0.;
		ibm_intp.i1 = i;
		ibm_intp.j1 = j;
		ibm_intp.k1 = k+1;

		ibm_intp.i2 = i;
		ibm_intp.j2 = j;
		ibm_intp.k2 = k+1;
		ibm_intp.i3 = i;
		ibm_intp.j3 = j;
		ibm_intp.k3 = k+1;
	      }
	      else {
		ptmp = coor[k-1][j][i];
		ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				    (coor[k][j][i].x - ptmp.x) +
				    (coor[k][j][i].y - ptmp.y) *
				    (coor[k][j][i].y - ptmp.y) +
				    (coor[k][j][i].z - ptmp.z) *
				    (coor[k][j][i].z - ptmp.z));
		ibm_intp.cr1 = 1.;
		ibm_intp.cr2 = 0.;
		ibm_intp.cr3 = 0.;
		ibm_intp.i1 = i;
		ibm_intp.j1 = j;
		ibm_intp.k1 = k-1;

		ibm_intp.i2 = i;
		ibm_intp.j2 = j;
		ibm_intp.k2 = k-1;
		ibm_intp.i3 = i;
		ibm_intp.j3 = j;
		ibm_intp.k3 = k-1;
	      }
	    }
	    else if (k==1 || k==mz-2) {
	      cell = ibm_intp.cell;
	      ptmp = coor[k][j+1][i];
	      ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				  (coor[k][j][i].x - ptmp.x) +
				  (coor[k][j][i].y - ptmp.y) *
				  (coor[k][j][i].y - ptmp.y) +
				  (coor[k][j][i].z - ptmp.z) *
				  (coor[k][j][i].z - ptmp.z));
	      ibm_intp.cr1 = 1.;
	      ibm_intp.cr2 = 0.;
	      ibm_intp.cr3 = 0.;
	      ibm_intp.i1 = i;
	      ibm_intp.j1 = j+1;
	      ibm_intp.k1 = k;
	      
	      ibm_intp.i2 = i;
	      ibm_intp.j2 = j+1;
	      ibm_intp.k2 = k;
	      ibm_intp.i3 = i;
	      ibm_intp.j3 = j+1;
	      ibm_intp.k3 = k;
	    }
	    else {
	      PetscPrintf(PETSC_COMM_SELF, "%%%%IBM Searching Fail! %d %d %d\n", i, j, k);
	    }
	  }

	  AddIBMNode(&user->ibmlist[ibi], ibm_intp);
	}
      }
    }
  }

  PetscBarrier(PETSC_NULL);

   // Back to the old nvert 3 and 1 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) nvert[k][j][i]=1.4;
	if ((int)(nvert[k][j][i]+0.5) == 4) nvert[k][j][i]=3.4;
      }
    }
  }


  DMDAVecRestoreArray(fda, user->lCent,&coor);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert);
  DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert);

  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	destroy(&cell_trg[k*ncx*ncy+j*ncx+i]);
      }
    }
  }

  PetscFree(cell_trg);
  PetscFree(ibm->qvec);
  PetscFree(ibm->radvec); 

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
  PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt	i, j, k, ln_v, n1e, n2e, n3e, nintp;
  
  PetscInt	nvert_l;
  PetscReal	dt[1000], ndotn, dirdotn;
  Cmpnts        dnn[1000],nn;

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

    //    for (ln_v=0; ln_v<ibm->n_elmt; ln_v++) {
    for (k=kp; k<ncz; k++) {
      current = cell_trg[k*ncx*ncy+j*ncx+i].head;
      while (current) {
	ln_v = current->Node;
	if (!Element_Searched[ln_v]) {
	  Element_Searched[ln_v] = PETSC_TRUE;
	  n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];
	  nn.x=nf_x[ln_v]; nn.y=nf_y[ln_v]; nn.z=nf_z[ln_v];

	  orig[0] = p.x; orig[1] = p.y, orig[2] = p.z;

	  vert0[0] = x_bp[n1e]; vert0[1] = y_bp[n1e]; vert0[2] = z_bp[n1e];
	  vert1[0] = x_bp[n2e]; vert1[1] = y_bp[n2e]; vert1[2] = z_bp[n2e];
	  vert2[0] = x_bp[n3e]; vert2[1] = y_bp[n3e]; vert2[2] = z_bp[n3e];
            
	  dirdotn=dir[0]*nn.x+dir[1]*nn.y+dir[2]*nn.z;

	  nvert_l = intsect_triangle(orig, dir, vert0, vert1, vert2, &t, &u, &v);
	
	  if (flg) 
	    PetscPrintf(PETSC_COMM_SELF, "elm, %d %d %le %le %le %d %d %d %le\n",ln_v,nvert_l,t,u,v,n1e,n2e,n3e,dirdotn);
	  
	  if (nvert_l > 0 && t>0) {
	    dt[nintp] = t;
	    dnn[nintp].x=nn.x;dnn[nintp].y=nn.y;dnn[nintp].z=nn.z;

	    nintp ++;
	    PetscInt temp;
	    for (temp = 0; temp < nintp-1; temp++) {
	      // Two interception points are the same, this leads to huge
	      // trouble for crossing number test
	      // Rather to program for all cases, we use a new line to
	      // repeat the test
	      ndotn=dnn[temp].x*nn.x+dnn[temp].y*nn.y+dnn[temp].z*nn.z;	      
	      
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

PetscInt point_cell_thin(Cmpnts p,Cmpnts p1,Cmpnts p2,Cmpnts p3,Cmpnts p4,
			 PetscInt ip, PetscInt jp, PetscInt kp,
			 IBMNodes *ibm, PetscInt ncx, PetscInt ncy,
			 PetscInt ncz, PetscReal dcx, PetscReal dcy,
			 PetscReal xbp_min, PetscReal ybp_min,
			 PetscReal zbp_max, List *cell_trg,
			 PetscInt flg)
{


  PetscInt	i, j, k, ln_v;//, n1e, n2e, n3e, nintp;
  PetscInt    cut;
 
  j = jp; i = ip; k = kp;

  node *current;
  
  if (flg)
    PetscPrintf(PETSC_COMM_SELF, "ip,jp,kp,nc_xyz %d %d %d %d %d %d\n",i,j,k,ncx,ncy,ncz);
  
  for (ln_v=0;ln_v<ibm->n_v;ln_v++) {
    
    if (flg) PetscPrintf(PETSC_COMM_SELF, "test010 ln_v %d \n",ln_v);
      
	cut = ISLineTriangleIntp(p,p1,ibm,ln_v);
	if (cut) return(2);
	
	cut = ISLineTriangleIntp(p,p2,ibm,ln_v);
	if (cut) return(2);
	
	cut = ISLineTriangleIntp(p,p3,ibm,ln_v);
	if (cut) return(2);

	cut = ISLineTriangleIntp(p,p4,ibm,ln_v);
	if (cut) return(2);

      current = current->next;
   
  }//for
  
  return(0);
}


/* Implementing the closest triangle algorithm described as the attached
   point-pairs.pdf */

PetscErrorCode BoundingSphere(IBMNodes *ibm)
{
  PetscInt    *nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3 = ibm->nv3;
  PetscReal    *x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt    n_elmt = ibm->n_elmt;
  
  PetscInt    ln_v;
  
  Cmpnts    p1, p2, p3, p0;
  
  PetscInt    n1e, n2e, n3e;
  
  p0.x = 0; p0.y = 0; p0.z = 0;

  PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->qvec));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->radvec));

  Cmpnts *qvec = ibm->qvec;
  PetscReal *radvec = ibm->radvec;

  Cmpnts pa, pb, pc, pu, pv, pf, pd, pt;
  PetscReal l12, l23, l31;
  PetscReal gama, lamda;
  for (ln_v = 0; ln_v < n_elmt; ln_v ++) {
    n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];

    p1.x = x_bp[n1e]; p1.y = y_bp[n1e]; p1.z = z_bp[n1e];
    p2.x = x_bp[n2e]; p2.y = y_bp[n2e]; p2.z = z_bp[n2e];
    p3.x = x_bp[n3e]; p3.y = y_bp[n3e]; p3.z = z_bp[n3e];

    l12 = Dist(p1, p2); l23 = Dist(p2, p3); l31 = Dist(p3, p1);

    /* Find the longest edge and assign the corresponding two vertices
       to pa and pb */
    if (l12 > l23) {
      if (l12 > l31) {
    pa = p1; pb = p2; pc = p3;
      }
      else {
    pa = p3; pb = p1; pc = p2;
      }
    }
    else {
      if (l31 < l23) {
    pa = p2; pb = p3; pc = p1;
      }
      else {
    pa = p3; pb = p1; pc = p2;
      }
    }

    pf.x = 0.5 * (pa.x + pb.x);
    pf.y = 0.5 * ( pa.y + pb.y);
    pf.z = 0.5 * (pa.z + pb.z);

    // u = a - f; v = c - f;
    VecAMinusB(pu, pa, pf);
    VecAMinusB(pv, pc, pf);

    // d = (u X v) X u;
    Cross(pt, pu, pv);
    Cross(pd, pt, pu);

    // gama = (u^2 - v^2) / (2 d \dot (v - u));
    gama = (Dist(pu, p0)*Dist(pu, p0) - Dist(pv, p0) * Dist(pv, p0));

    VecAMinusB(pt, pv, pu);
    lamda = 2 * (pd.x * pt.x + pd.y * pt.y + pd.z * pt.z);

    gama /= lamda;
    
    if (gama <0) {
      lamda = 0;
    }
    else {
      lamda = gama;
    }
    
    qvec[ln_v].x = pf.x + lamda * pd.x;
    qvec[ln_v].y = pf.y + lamda * pd.y;
    qvec[ln_v].z = pf.z + lamda * pd.z;

    radvec[ln_v] = Dist(qvec[ln_v], pa);
  }

  return 0;
}


PetscErrorCode nearestcell1(Cmpnts p, IBMNodes *ibm, IBMInfo *ibminfo)
{
  PetscInt    *nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3 = ibm->nv3;
  PetscReal    *nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;
  PetscReal    *x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt    n_elmt = ibm->n_elmt;
  
  PetscInt    ln_v;

  Cmpnts    p1, p2, p3;
  PetscReal    tf;

  PetscInt    n1e, n2e, n3e;
  PetscReal    nfx, nfy, nfz;

  Cmpnts    pj; // projection point
  Cmpnts    pmin, po;
  PetscReal    dmin, d;
  PetscInt    cell_min;
  dmin = 1.e20;

  cell_min = -100;

  PetscReal d_center;

  for (ln_v=0; ln_v<n_elmt; ln_v++) {
    d_center = Dist(p, ibm->qvec[ln_v]);
    if (d_center - ibm->radvec[ln_v] < dmin) {
      n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];
      nfx = nf_x[ln_v];
      nfy = nf_y[ln_v];
      nfz = nf_z[ln_v];

      p1.x = x_bp[n1e]; p1.y = y_bp[n1e]; p1.z = z_bp[n1e];
      p2.x = x_bp[n2e]; p2.y = y_bp[n2e]; p2.z = z_bp[n2e];
      p3.x = x_bp[n3e]; p3.y = y_bp[n3e]; p3.z = z_bp[n3e];

      tf = ((p.x - x_bp[n1e]) * nfx +
        (p.y - y_bp[n1e]) * nfy +
        (p.z - z_bp[n1e]) * nfz);
    
      if (fabs(tf) < 1.e-10) tf = 1.e-10;
      if (tf>=0) { // Point p locates on the positive side of surface triangle
  
    pj.x = p.x - tf * nfx;
    pj.y = p.y - tf * nfy;
    pj.z = p.z - tf * nfz;


  
    if (ISPointInTriangle(pj, p1, p2, p3, nfx, nfy, nfz) == 1) { /* The projected point
                                    is inside the
                                    triangle */
     
      if (tf < dmin) {
        dmin = tf;
        pmin.x = pj.x;
        pmin.y = pj.y;
        pmin.z = pj.z;
        cell_min = ln_v;
      }
    }
    else {
      Dis_P_Line(p, p1, p2, &po, &d);
      if (d < dmin) {
        dmin = d;
        pmin.x = po.x;
        pmin.y = po.y;
        pmin.z = po.z;

        cell_min = ln_v;
      }
      Dis_P_Line(p, p2, p3, &po, &d);
      if (d < dmin) {
        dmin = d;
        pmin.x = po.x;
        pmin.y = po.y;
        pmin.z = po.z;

        cell_min = ln_v;
      }
      Dis_P_Line(p, p3, p1, &po, &d);
      if (d < dmin) {
        dmin = d;
        pmin.x = po.x;
        pmin.y = po.y;
        pmin.z = po.z;
        cell_min = ln_v;
      }      
    }
      }
    }
  }

  if (cell_min == -100) {
    PetscPrintf(PETSC_COMM_SELF, "Nearest Cell Searching Error!\n");
    exit(0);
  }
   
  ibminfo->cell = cell_min;
  ibminfo->pmin = pmin;
  ibminfo->d_s = dmin;

  Cpt2D pjp, pj1, pj2, pj3;
  nfx = nf_x[cell_min]; nfy = nf_y[cell_min]; nfz=nf_z[cell_min];

  n1e = nv1[cell_min]; n2e = nv2[cell_min]; n3e = nv3[cell_min];
  p1.x = x_bp[n1e]; p1.y = y_bp[n1e]; p1.z = z_bp[n1e];
  p2.x = x_bp[n2e]; p2.y = y_bp[n2e]; p2.z = z_bp[n2e];
  p3.x = x_bp[n3e]; p3.y = y_bp[n3e]; p3.z = z_bp[n3e];

  if (fabs(nfx) >= fabs(nfy) && fabs(nfx)>= fabs(nfz)) {
    pjp.x = pmin.y; pjp.y = pmin.z;
    pj1.x = p1.y;   pj1.y = p1.z;
    pj2.x = p2.y;   pj2.y = p2.z;
    pj3.x = p3.y;   pj3.y = p3.z;
    triangle_intp2(pjp, pj1, pj2, pj3, ibminfo);
  }
  else if (fabs(nfy) >= fabs(nfx) && fabs(nfy)>= fabs(nfz)) {
    pjp.x = pmin.x; pjp.y = pmin.z;
    pj1.x = p1.x;   pj1.y = p1.z;
    pj2.x = p2.x;   pj2.y = p2.z;
    pj3.x = p3.x;   pj3.y = p3.z;
    triangle_intp2(pjp, pj1, pj2, pj3, ibminfo);
  }
  else if (fabs(nfz) >= fabs(nfy) && fabs(nfz)>= fabs(nfx)) {
    pjp.x = pmin.y; pjp.y = pmin.x;
    pj1.x = p1.y;   pj1.y = p1.x;
    pj2.x = p2.y;   pj2.y = p2.x;
    pj3.x = p3.y;   pj3.y = p3.x;
    triangle_intp2(pjp, pj1, pj2, pj3, ibminfo);
  }
  if (ibminfo->cs1 != ibminfo->cs1) {
    PetscPrintf(PETSC_COMM_SELF, "INTP2 %e %e %e %i %i %i\n", nfx, nfy, nfz, n1e, n2e, n3e);
  }
}



PetscErrorCode nearestcell(Cmpnts p, IBMNodes *ibm, IBMInfo *ibminfo)
{
  PetscInt	*nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3 = ibm->nv3;
  PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt	n_elmt = ibm->n_elmt;
  
  PetscInt	ln_v;

  Cmpnts	p1, p2, p3;
  PetscReal	tf;

  PetscInt	n1e, n2e, n3e;
  PetscReal	nfx, nfy, nfz;

  Cmpnts	pj; // projection point
  Cmpnts	pmin, po;
  PetscReal	dmin, d;
  PetscInt	cell_min;
  dmin = 1.e20;

  cell_min = -100;
  for (ln_v=0; ln_v<n_elmt; ln_v++) {
    n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];
    nfx = nf_x[ln_v];
    nfy = nf_y[ln_v];
    nfz = nf_z[ln_v];

    p1.x = x_bp[n1e]; p1.y = y_bp[n1e]; p1.z = z_bp[n1e];
    p2.x = x_bp[n2e]; p2.y = y_bp[n2e]; p2.z = z_bp[n2e];
    p3.x = x_bp[n3e]; p3.y = y_bp[n3e]; p3.z = z_bp[n3e];

    tf = ((p.x - x_bp[n1e]) * nfx +
	  (p.y - y_bp[n1e]) * nfy +
	  (p.z - z_bp[n1e]) * nfz);
   
    if (fabs(tf) < 1.e-10) tf = 1.e-10;
    if (tf>=0) { // Point p locates on the positive side of surface triangle
    
      pj.x = p.x - tf * nfx;
      pj.y = p.y - tf * nfy;
      pj.z = p.z - tf * nfz;


     
      if (ISPointInTriangle(pj, p1, p2, p3, nfx, nfy, nfz) == 1) { /* The projected point
							is inside the
							triangle */
	if (tf < dmin) {
	  dmin = tf;
	  pmin.x = pj.x;
	  pmin.y = pj.y;
	  pmin.z = pj.z;
	  cell_min = ln_v;
	}
      }
      else {
	Dis_P_Line(p, p1, p2, &po, &d);
	if (d < dmin) {
	  dmin = d;
	  pmin.x = po.x;
	  pmin.y = po.y;
	  pmin.z = po.z;

	  cell_min = ln_v;
	}
	Dis_P_Line(p, p2, p3, &po, &d);
	if (d < dmin) {
	  dmin = d;
	  pmin.x = po.x;
	  pmin.y = po.y;
	  pmin.z = po.z;

	  cell_min = ln_v;
	}
	Dis_P_Line(p, p3, p1, &po, &d);
	if (d < dmin) {
	  dmin = d;
	  pmin.x = po.x;
	  pmin.y = po.y;
	  pmin.z = po.z;
	  cell_min = ln_v;
	}      
      }
    }
    if (flagprint == 1) {
      PetscPrintf(PETSC_COMM_SELF, "TF %d %d %d %d %le %le %le %le %le %le %le %le\n", n1e, n2e, n3e, ln_v, tf,dmin,nfx,nfy,nfz,p1.x,p1.y,p1.z);
    }

  }

  if (cell_min == -100) {
    PetscPrintf(PETSC_COMM_SELF, "Nearest Cell Searching Error! %le %le %le %le\n",p.x,p.y,p.z,tf);
    exit(0);
  }
   if (flagprint==1)
      PetscPrintf(PETSC_COMM_SELF, "L %d %e %e %e %d p. %le %le %le\n", cell_min, pmin.x, pmin.y, pmin.z, n1e,p.x,p.y,p.z);
  
  ibminfo->cell = cell_min;
  ibminfo->pmin = pmin;
  ibminfo->d_s = dmin;

  Cpt2D pjp, pj1, pj2, pj3;
  nfx = nf_x[cell_min]; nfy = nf_y[cell_min]; nfz=nf_z[cell_min];

  n1e = nv1[cell_min]; n2e = nv2[cell_min]; n3e = nv3[cell_min];
  p1.x = x_bp[n1e]; p1.y = y_bp[n1e]; p1.z = z_bp[n1e];
  p2.x = x_bp[n2e]; p2.y = y_bp[n2e]; p2.z = z_bp[n2e];
  p3.x = x_bp[n3e]; p3.y = y_bp[n3e]; p3.z = z_bp[n3e];

  if (fabs(nfx) >= fabs(nfy) && fabs(nfx)>= fabs(nfz)) {
    pjp.x = pmin.y; pjp.y = pmin.z;
    pj1.x = p1.y;   pj1.y = p1.z;
    pj2.x = p2.y;   pj2.y = p2.z;
    pj3.x = p3.y;   pj3.y = p3.z;
    triangle_intp2(pjp, pj1, pj2, pj3, ibminfo);
  }
  else if (fabs(nfy) >= fabs(nfx) && fabs(nfy)>= fabs(nfz)) {
    pjp.x = pmin.x; pjp.y = pmin.z;
    pj1.x = p1.x;   pj1.y = p1.z;
    pj2.x = p2.x;   pj2.y = p2.z;
    pj3.x = p3.x;   pj3.y = p3.z;
    triangle_intp2(pjp, pj1, pj2, pj3, ibminfo);
  }
  else if (fabs(nfz) >= fabs(nfy) && fabs(nfz)>= fabs(nfx)) {
    pjp.x = pmin.y; pjp.y = pmin.x;
    pj1.x = p1.y;   pj1.y = p1.x;
    pj2.x = p2.y;   pj2.y = p2.x;
    pj3.x = p3.y;   pj3.y = p3.x;
    triangle_intp2(pjp, pj1, pj2, pj3, ibminfo);
  }
  if (ibminfo->cs1 != ibminfo->cs1) {
    PetscPrintf(PETSC_COMM_SELF, "INTP2 %e %e %e %i %i %i\n", nfx, nfy, nfz, n1e, n2e, n3e);
  }
}


PetscErrorCode ICP(Cmpnts p, Cmpnts pc[9], PetscReal nfx, PetscReal nfy,
	       PetscReal nfz, IBMInfo *ibminfo, PetscInt *ip,
	       PetscInt *jp  , PetscInt *kp)
{
  PetscInt 	triangles[3][8];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo->imode = -100;

  triangles[0][0] = 0; triangles[1][0] = 1; triangles[2][0] = 4;
  triangles[0][1] = 1; triangles[1][1] = 2; triangles[2][1] = 4;
  triangles[0][2] = 2; triangles[1][2] = 4; triangles[2][2] = 5;
  triangles[0][3] = 4; triangles[1][3] = 5; triangles[2][3] = 8;
  triangles[0][4] = 4; triangles[1][4] = 7; triangles[2][4] = 8;
  triangles[0][5] = 4; triangles[1][5] = 6; triangles[2][5] = 7;
  triangles[0][6] = 3; triangles[1][6] = 4; triangles[2][6] = 6;
  triangles[0][7] = 3; triangles[1][7] = 4; triangles[2][7] = 0;

  for (i=0; i<8; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;

    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      

      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
      
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;

	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;

	  if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	    pjp.x = pint.y; pjp.y = pint.z;
	    pj1.x = p1.y;   pj1.y = p1.z;
	    pj2.x = p2.y;   pj2.y = p2.z;
	    pj3.x = p3.y;   pj3.y = p3.z;
	    triangle_intp(pjp, pj1, pj2, pj3, ibminfo);
	  }
	  else if	(fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	    pjp.x = pint.x; pjp.y = pint.z;
	    pj1.x = p1.x;   pj1.y = p1.z;
	    pj2.x = p2.x;   pj2.y = p2.z;
	    pj3.x = p3.x;   pj3.y = p3.z;
	    triangle_intp(pjp, pj1, pj2, pj3, ibminfo);
	  }
	  else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	    pjp.x = pint.y; pjp.y = pint.x;
	    pj1.x = p1.y;   pj1.y = p1.x;
	    pj2.x = p2.y;   pj2.y = p2.x;
	    pj3.x = p3.y;   pj3.y = p3.x;
	    triangle_intp(pjp, pj1, pj2, pj3, ibminfo);
	  }

	  ibminfo->d_i = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
				     (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo->imode = cell;

	  return (0);
	}
      }
    }
  }
  return 0;
}

PetscErrorCode InterceptionPoint(Cmpnts p, PetscInt i, PetscInt j, PetscInt k,
			     IBMInfo *ibminfo, UserCtx *user)
{
  DM		da =user->da, fda=user->fda;
  PetscInt	ip[9], jp[9], kp[9];
  Cmpnts	pc[9];

  PetscInt	nif;

  PetscReal	nfx, nfy, nfz;
  PetscReal 	dr;
  Vec		Coor;
  Cmpnts	***coor;


  DMDAVecGetArray(fda, user->Cent, &coor);

  nfx = p.x - ibminfo->pmin.x;
  nfy = p.y - ibminfo->pmin.y;
  nfz = p.z - ibminfo->pmin.z;

  dr = sqrt(nfx*nfx + nfy*nfy + nfz*nfz);
  nfx /= dr; nfy /= dr; nfz /=dr;

  flagprint = 0;
  
  ip[0] = i-1; ip[1] = i-1; ip[2] = i-1;
  ip[3] = i-1; ip[4] = i-1; ip[5] = i-1;
  ip[6] = i-1; ip[7] = i-1; ip[8] = i-1;
  
  jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
  jp[1] = j;   jp[4] = j;   jp[7] = j;
  jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;
  
  kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
  kp[3] = k;   kp[4] = k;   kp[5] = k;
  kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;
  
  for (nif=0; nif<9; nif++) {
    pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
    pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
    pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
  }

  ICP(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);

  switch (ibminfo->imode) {
    
  case(0): {
    ibminfo->i1=ip[0]; ibminfo->j1 = jp[0]; ibminfo->k1 = kp[0];
    ibminfo->i2=ip[1]; ibminfo->j2 = jp[1]; ibminfo->k2 = kp[1];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (1): {
    ibminfo->i1=ip[1]; ibminfo->j1 = jp[1]; ibminfo->k1 = kp[1];
    ibminfo->i2=ip[2]; ibminfo->j2 = jp[2]; ibminfo->k2 = kp[2];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (2): {
    ibminfo->i1=ip[2]; ibminfo->j1 = jp[2]; ibminfo->k1 = kp[2];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[5]; ibminfo->j3 = jp[5]; ibminfo->k3 = kp[5];
    break;
  }
  case (3): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[5]; ibminfo->j2 = jp[5]; ibminfo->k2 = kp[5];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (4): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[7]; ibminfo->j2 = jp[7]; ibminfo->k2 = kp[7];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (5): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[6]; ibminfo->j2 = jp[6]; ibminfo->k2 = kp[6];
    ibminfo->i3=ip[7]; ibminfo->j3 = jp[7]; ibminfo->k3 = kp[7];
    break;
  }
  case (6): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[6]; ibminfo->j3 = jp[6]; ibminfo->k3 = kp[6];
      break;
  }
  case (7): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[0]; ibminfo->j3 = jp[0]; ibminfo->k3 = kp[0];
    break;
  }
  }
  if (ibminfo->imode >=0) {
    DMDAVecRestoreArray(fda, user->Cent, &coor);return(0);
  }
    
  ip[0] = i+1; ip[1] = i+1; ip[2] = i+1;
  ip[3] = i+1; ip[4] = i+1; ip[5] = i+1;
  ip[6] = i+1; ip[7] = i+1; ip[8] = i+1;
  
  jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
  jp[1] = j;   jp[4] = j;   jp[7] = j;
  jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;
  
  kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
  kp[3] = k;   kp[4] = k;   kp[5] = k;
  kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;
  
  for (nif=0; nif<9; nif++) {
    pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
    pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
    pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
  }
  
  ICP(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
  
  switch (ibminfo->imode) {
  case(0): {
    ibminfo->i1=ip[0]; ibminfo->j1 = jp[0]; ibminfo->k1 = kp[0];
    ibminfo->i2=ip[1]; ibminfo->j2 = jp[1]; ibminfo->k2 = kp[1];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (1): {
    ibminfo->i1=ip[1]; ibminfo->j1 = jp[1]; ibminfo->k1 = kp[1];
    ibminfo->i2=ip[2]; ibminfo->j2 = jp[2]; ibminfo->k2 = kp[2];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (2): {
    ibminfo->i1=ip[2]; ibminfo->j1 = jp[2]; ibminfo->k1 = kp[2];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[5]; ibminfo->j3 = jp[5]; ibminfo->k3 = kp[5];
    break;
  }
  case (3): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[5]; ibminfo->j2 = jp[5]; ibminfo->k2 = kp[5];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (4): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[7]; ibminfo->j2 = jp[7]; ibminfo->k2 = kp[7];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (5): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[6]; ibminfo->j2 = jp[6]; ibminfo->k2 = kp[6];
    ibminfo->i3=ip[7]; ibminfo->j3 = jp[7]; ibminfo->k3 = kp[7];
    break;
  }
  case (6): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[6]; ibminfo->j3 = jp[6]; ibminfo->k3 = kp[6];
    break;
  }
  case (7): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[0]; ibminfo->j3 = jp[0]; ibminfo->k3 = kp[0];
    break;
  }
  }
  if (ibminfo->imode >=0) {
    DMDAVecRestoreArray(fda, user->Cent, &coor);return(0);
  }
  
  ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
  ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
  ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;
  
  jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
  jp[1] = j-1; jp[4] = j-1; jp[7] = j-1;
  jp[2] = j-1; jp[5] = j-1; jp[8] = j-1;
  
  kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
  kp[3] = k;   kp[4] = k;   kp[5] = k;
  kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;
  
  for (nif=0; nif<9; nif++) {
    pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
    pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
    pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
  }
  
  ICP(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
  switch (ibminfo->imode) {
  case(0): {
    ibminfo->i1=ip[0]; ibminfo->j1 = jp[0]; ibminfo->k1 = kp[0];
    ibminfo->i2=ip[1]; ibminfo->j2 = jp[1]; ibminfo->k2 = kp[1];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (1): {
    ibminfo->i1=ip[1]; ibminfo->j1 = jp[1]; ibminfo->k1 = kp[1];
    ibminfo->i2=ip[2]; ibminfo->j2 = jp[2]; ibminfo->k2 = kp[2];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (2): {
    ibminfo->i1=ip[2]; ibminfo->j1 = jp[2]; ibminfo->k1 = kp[2];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[5]; ibminfo->j3 = jp[5]; ibminfo->k3 = kp[5];
    break;
  }
  case (3): {
      ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
      ibminfo->i2=ip[5]; ibminfo->j2 = jp[5]; ibminfo->k2 = kp[5];
      ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
      break;
  }
  case (4): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[7]; ibminfo->j2 = jp[7]; ibminfo->k2 = kp[7];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (5): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[6]; ibminfo->j2 = jp[6]; ibminfo->k2 = kp[6];
    ibminfo->i3=ip[7]; ibminfo->j3 = jp[7]; ibminfo->k3 = kp[7];
    break;
  }
  case (6): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[6]; ibminfo->j3 = jp[6]; ibminfo->k3 = kp[6];
    break;
  }
  case (7): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[0]; ibminfo->j3 = jp[0]; ibminfo->k3 = kp[0];
    break;
  }
  }
  if (ibminfo->imode >=0) {
    DMDAVecRestoreArray(fda, user->Cent, &coor);return(0);
  }
  
  ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
  ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
  ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;
  
  jp[0] = j+1; jp[3] = j+1; jp[6] = j+1;
  jp[1] = j+1; jp[4] = j+1; jp[7] = j+1;
  jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;
  
  kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
  kp[3] = k;   kp[4] = k;   kp[5] = k;
  kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;
  
  for (nif=0; nif<9; nif++) {
    pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
    pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
    pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
  }
  
  ICP(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
  switch (ibminfo->imode) {
  case(0): {
    ibminfo->i1=ip[0]; ibminfo->j1 = jp[0]; ibminfo->k1 = kp[0];
    ibminfo->i2=ip[1]; ibminfo->j2 = jp[1]; ibminfo->k2 = kp[1];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (1): {
    ibminfo->i1=ip[1]; ibminfo->j1 = jp[1]; ibminfo->k1 = kp[1];
    ibminfo->i2=ip[2]; ibminfo->j2 = jp[2]; ibminfo->k2 = kp[2];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (2): {
    ibminfo->i1=ip[2]; ibminfo->j1 = jp[2]; ibminfo->k1 = kp[2];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[5]; ibminfo->j3 = jp[5]; ibminfo->k3 = kp[5];
    break;
  }
  case (3): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[5]; ibminfo->j2 = jp[5]; ibminfo->k2 = kp[5];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (4): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[7]; ibminfo->j2 = jp[7]; ibminfo->k2 = kp[7];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (5): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[6]; ibminfo->j2 = jp[6]; ibminfo->k2 = kp[6];
    ibminfo->i3=ip[7]; ibminfo->j3 = jp[7]; ibminfo->k3 = kp[7];
    break;
  }
  case (6): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[6]; ibminfo->j3 = jp[6]; ibminfo->k3 = kp[6];
    break;
  }
  case (7): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[0]; ibminfo->j3 = jp[0]; ibminfo->k3 = kp[0];
    break;
  }
  }
  if (ibminfo->imode >=0) {
    DMDAVecRestoreArray(fda, user->Cent, &coor);return(0);
  }
  ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
  ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
  ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;
  
  jp[0] = j-1; jp[3] = j  ; jp[6] = j+1;
  jp[1] = j-1; jp[4] = j  ; jp[7] = j+1;
  jp[2] = j-1; jp[5] = j  ; jp[8] = j+1;
  
  kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
  kp[3] = k-1; kp[4] = k-1; kp[5] = k-1;
  kp[6] = k-1; kp[7] = k-1; kp[8] = k-1;
  
  for (nif=0; nif<9; nif++) {
    pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
    pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
    pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
  }
  
  ICP(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
  switch (ibminfo->imode) {
  case(0): {
    ibminfo->i1=ip[0]; ibminfo->j1 = jp[0]; ibminfo->k1 = kp[0];
    ibminfo->i2=ip[1]; ibminfo->j2 = jp[1]; ibminfo->k2 = kp[1];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (1): {
    ibminfo->i1=ip[1]; ibminfo->j1 = jp[1]; ibminfo->k1 = kp[1];
    ibminfo->i2=ip[2]; ibminfo->j2 = jp[2]; ibminfo->k2 = kp[2];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (2): {
    ibminfo->i1=ip[2]; ibminfo->j1 = jp[2]; ibminfo->k1 = kp[2];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[5]; ibminfo->j3 = jp[5]; ibminfo->k3 = kp[5];
    break;
  }
  case (3): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[5]; ibminfo->j2 = jp[5]; ibminfo->k2 = kp[5];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (4): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[7]; ibminfo->j2 = jp[7]; ibminfo->k2 = kp[7];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (5): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[6]; ibminfo->j2 = jp[6]; ibminfo->k2 = kp[6];
    ibminfo->i3=ip[7]; ibminfo->j3 = jp[7]; ibminfo->k3 = kp[7];
    break;
  }
  case (6): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[6]; ibminfo->j3 = jp[6]; ibminfo->k3 = kp[6];
    break;
  }
  case (7): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[0]; ibminfo->j3 = jp[0]; ibminfo->k3 = kp[0];
    break;
  }
  }
  if (ibminfo->imode >=0) {
    DMDAVecRestoreArray(fda, user->Cent, &coor);return(0);
  }
  
  ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
  ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
  ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;
  
  jp[0] = j-1; jp[3] = j  ; jp[6] = j+1;
  jp[1] = j-1; jp[4] = j  ; jp[7] = j+1;
  jp[2] = j-1; jp[5] = j  ; jp[8] = j+1;
  
  kp[0] = k+1; kp[1] = k+1; kp[2] = k+1;
  kp[3] = k+1; kp[4] = k+1; kp[5] = k+1;
  kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;
  
  for (nif=0; nif<9; nif++) {
    pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
    pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
    pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
  }
  
  ICP(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
  switch (ibminfo->imode) {
  case(0): {
    ibminfo->i1=ip[0]; ibminfo->j1 = jp[0]; ibminfo->k1 = kp[0];
    ibminfo->i2=ip[1]; ibminfo->j2 = jp[1]; ibminfo->k2 = kp[1];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (1): {
    ibminfo->i1=ip[1]; ibminfo->j1 = jp[1]; ibminfo->k1 = kp[1];
    ibminfo->i2=ip[2]; ibminfo->j2 = jp[2]; ibminfo->k2 = kp[2];
    ibminfo->i3=ip[4]; ibminfo->j3 = jp[4]; ibminfo->k3 = kp[4];
    break;
  }
  case (2): {
    ibminfo->i1=ip[2]; ibminfo->j1 = jp[2]; ibminfo->k1 = kp[2];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[5]; ibminfo->j3 = jp[5]; ibminfo->k3 = kp[5];
    break;
  }
  case (3): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[5]; ibminfo->j2 = jp[5]; ibminfo->k2 = kp[5];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (4): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[7]; ibminfo->j2 = jp[7]; ibminfo->k2 = kp[7];
    ibminfo->i3=ip[8]; ibminfo->j3 = jp[8]; ibminfo->k3 = kp[8];
    break;
  }
  case (5): {
    ibminfo->i1=ip[4]; ibminfo->j1 = jp[4]; ibminfo->k1 = kp[4];
    ibminfo->i2=ip[6]; ibminfo->j2 = jp[6]; ibminfo->k2 = kp[6];
    ibminfo->i3=ip[7]; ibminfo->j3 = jp[7]; ibminfo->k3 = kp[7];
    break;
  }
  case (6): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[6]; ibminfo->j3 = jp[6]; ibminfo->k3 = kp[6];
    break;
  }
  case (7): {
    ibminfo->i1=ip[3]; ibminfo->j1 = jp[3]; ibminfo->k1 = kp[3];
    ibminfo->i2=ip[4]; ibminfo->j2 = jp[4]; ibminfo->k2 = kp[4];
    ibminfo->i3=ip[0]; ibminfo->j3 = jp[0]; ibminfo->k3 = kp[0];
    break;
  }
  }
  
  if (ibminfo->imode >=0) {
    DMDAVecRestoreArray(fda, user->Cent, &coor);return(0);
  }
  
  
  DMDAVecRestoreArray(fda, user->Cent, &coor);
  
  return(0);
}

PetscErrorCode ICP2(Cmpnts p, Cmpnts pc[9], PetscReal nfx, PetscReal nfy,
	       PetscReal nfz, IBMInfo *ibminfo, PetscInt *ip,
	       PetscInt *jp  , PetscInt *kp)
{
  PetscInt 	triangles[3][8];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3;

  Cpt2D		pj1, pj2, pj3, pjp;
  PetscInt	cell, flag;

  PetscInt	i;
  Cmpnts	pint; // Interception point
  PetscReal	nfxt, nfyt, nfzt;

  ibminfo->iimode = -100;

  triangles[0][0] = 0; triangles[1][0] = 1; triangles[2][0] = 3;
  triangles[0][1] = 1; triangles[1][1] = 2; triangles[2][1] = 5;
  triangles[0][2] = 1; triangles[1][2] = 4; triangles[2][2] = 5;
  triangles[0][3] = 4; triangles[1][3] = 5; triangles[2][3] = 7;
  triangles[0][4] = 5; triangles[1][4] = 7; triangles[2][4] = 8;
  triangles[0][5] = 4; triangles[1][5] = 3; triangles[2][5] = 7;
  triangles[0][6] = 3; triangles[1][6] = 7; triangles[2][6] = 6;
  triangles[0][7] = 3; triangles[1][7] = 4; triangles[2][7] = 1;

  for (i=0; i<8; i++) {
    p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];
    
    dx1 = p.x - p1.x; dy1 = p.y - p1.y; dz1 = p.z - p1.z;
    dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;
    dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;
    
    d = (nfx * (dy2 * dz3 - dz2 * dy3) - 
	 nfy * (dx2 * dz3 - dz2 * dx3) + 
	 nfz * (dx2 * dy3 - dy2 * dx3));
    if (fabs(d) > 1.e-10) {
      d = -(dx1 * (dy2 * dz3 - dz2 * dy3) - 
	    dy1 * (dx2 * dz3 - dz2 * dx3) + 
	    dz1 * (dx2 * dy3 - dy2 * dx3)) / d;
      
      
      if (d>0) {
	pint.x = p.x + d * nfx;
	pint.y = p.y + d * nfy;
	pint.z = p.z + d * nfz;

	rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
	rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;
	
	nfxt = ry1 * rz2 - rz1 * ry2;
	nfyt = -rx1 * rz2 + rz1 * rx2;
	nfzt = rx1 * ry2 - ry1 * rx2;
	
	flag = ISPointInTriangle(pint, p1, p2, p3, nfxt, nfyt, nfzt);
	if (flag >= 0) {
	  cell = i;
	  
	  if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt)) {
	    pjp.x = pint.y; pjp.y = pint.z;
	    pj1.x = p1.y;   pj1.y = p1.z;
	    pj2.x = p2.y;   pj2.y = p2.z;
	    pj3.x = p3.y;   pj3.y = p3.z;
	    triangle_intpp(pjp, pj1, pj2, pj3, ibminfo);
	  }
	  else if	(fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt)) {
	    pjp.x = pint.x; pjp.y = pint.z;
	    pj1.x = p1.x;   pj1.y = p1.z;
	    pj2.x = p2.x;   pj2.y = p2.z;
	    pj3.x = p3.x;   pj3.y = p3.z;
	    triangle_intpp(pjp, pj1, pj2, pj3, ibminfo);
	  }
	  else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt)) {
	    pjp.x = pint.y; pjp.y = pint.x;
	    pj1.x = p1.y;   pj1.y = p1.x;
	    pj2.x = p2.y;   pj2.y = p2.x;
	    pj3.x = p3.y;   pj3.y = p3.x;
	    triangle_intpp(pjp, pj1, pj2, pj3, ibminfo);
	  }

	  ibminfo->d_ii = sqrt((pint.x-p.x)*(pint.x - p.x) + (pint.y-p.y) *
			       (pint.y-p.y) + (pint.z - p.z)* (pint.z - p.z));
	  ibminfo->iimode = cell;

	  return (0);
	}
      }
    }
  }
  return 0;
}

PetscErrorCode InterceptionPoint2(Cmpnts p, PetscInt i, PetscInt j, PetscInt k,
			     IBMInfo *ibminfo, UserCtx *user)
{
  DM		da =user->da, fda=user->fda;
  PetscInt	ip[9], jp[9], kp[9];
  Cmpnts	pc[9];

  PetscInt	nif;

  PetscReal	nfx, nfy, nfz;
  PetscReal 	dr;
  Vec		Coor;
  Cmpnts	***coor;

  DMDAVecGetArray(fda, user->lCent, &coor);

  nfx = p.x - ibminfo->pmin.x;
  nfy = p.y - ibminfo->pmin.y;
  nfz = p.z - ibminfo->pmin.z;

  dr = sqrt(nfx*nfx + nfy*nfy + nfz*nfz);
  nfx /= dr; nfy /= dr; nfz /=dr;

  flagprint = 0;
 
    ip[0] = i-1; ip[1] = i-1; ip[2] = i-1;
    ip[3] = i-1; ip[4] = i-1; ip[5] = i-1;
    ip[6] = i-1; ip[7] = i-1; ip[8] = i-1;

    jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
    jp[1] = j;   jp[4] = j;   jp[7] = j;
    jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;

    kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
    kp[3] = k;   kp[4] = k;   kp[5] = k;
    kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

    for (nif=0; nif<9; nif++) {
      pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
      pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
      pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
      
    }

    ICP2(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
    switch (ibminfo->iimode) {
    case(0): {
      ibminfo->i11=ip[0]; ibminfo->j11 = jp[0]; ibminfo->k11 = kp[0];
      ibminfo->i22=ip[1]; ibminfo->j22 = jp[1]; ibminfo->k22 = kp[1];
      ibminfo->i33=ip[3]; ibminfo->j33 = jp[3]; ibminfo->k33 = kp[3];
      break;
    }
    case (1): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[2]; ibminfo->j22 = jp[2]; ibminfo->k22 = kp[2];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (2): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (3): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[5]; ibminfo->j22 = jp[5]; ibminfo->k22 = kp[5];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (4): {
      ibminfo->i11=ip[5]; ibminfo->j11 = jp[5]; ibminfo->k11 = kp[5];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[8]; ibminfo->j33 = jp[8]; ibminfo->k33 = kp[8];
      break;
    }
    case (5): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[3]; ibminfo->j22 = jp[3]; ibminfo->k22 = kp[3];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (6): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[6]; ibminfo->j33 = jp[6]; ibminfo->k33 = kp[6];
      break;
    }
    case (7): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[1]; ibminfo->j33 = jp[1]; ibminfo->k33 = kp[1];
      break;
    }
    }
    if (ibminfo->iimode >=0) {
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    ip[0] = i+1; ip[1] = i+1; ip[2] = i+1;
    ip[3] = i+1; ip[4] = i+1; ip[5] = i+1;
    ip[6] = i+1; ip[7] = i+1; ip[8] = i+1;

    jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
    jp[1] = j;   jp[4] = j;   jp[7] = j;
    jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;

    kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
    kp[3] = k;   kp[4] = k;   kp[5] = k;
    kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

    for (nif=0; nif<9; nif++) {
      pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
      pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
      pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
    }

    ICP2(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
    switch (ibminfo->iimode) {
    case(0): {
      ibminfo->i11=ip[0]; ibminfo->j11 = jp[0]; ibminfo->k11 = kp[0];
      ibminfo->i22=ip[1]; ibminfo->j22 = jp[1]; ibminfo->k22 = kp[1];
      ibminfo->i33=ip[3]; ibminfo->j33 = jp[3]; ibminfo->k33 = kp[3];
      break;
    }
    case (1): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[2]; ibminfo->j22 = jp[2]; ibminfo->k22 = kp[2];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (2): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (3): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[5]; ibminfo->j22 = jp[5]; ibminfo->k22 = kp[5];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (4): {
      ibminfo->i11=ip[5]; ibminfo->j11 = jp[5]; ibminfo->k11 = kp[5];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[8]; ibminfo->j33 = jp[8]; ibminfo->k33 = kp[8];
      break;
    }
    case (5): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[3]; ibminfo->j22 = jp[3]; ibminfo->k22 = kp[3];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (6): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[6]; ibminfo->j33 = jp[6]; ibminfo->k33 = kp[6];
      break;
    }
    case (7): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[1]; ibminfo->j33 = jp[1]; ibminfo->k33 = kp[1];
      break;
    }
    }
    if (ibminfo->iimode >=0) {
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
    ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
    ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

    jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
    jp[1] = j-1; jp[4] = j-1; jp[7] = j-1;
    jp[2] = j-1; jp[5] = j-1; jp[8] = j-1;

    kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
    kp[3] = k;   kp[4] = k;   kp[5] = k;
    kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

    for (nif=0; nif<9; nif++) {
      pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
      pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
      pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
    }

    ICP2(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
    switch (ibminfo->iimode) {
    case(0): {
      ibminfo->i11=ip[0]; ibminfo->j11 = jp[0]; ibminfo->k11 = kp[0];
      ibminfo->i22=ip[1]; ibminfo->j22 = jp[1]; ibminfo->k22 = kp[1];
      ibminfo->i33=ip[3]; ibminfo->j33 = jp[3]; ibminfo->k33 = kp[3];
      break;
    }
    case (1): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[2]; ibminfo->j22 = jp[2]; ibminfo->k22 = kp[2];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (2): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (3): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[5]; ibminfo->j22 = jp[5]; ibminfo->k22 = kp[5];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (4): {
      ibminfo->i11=ip[5]; ibminfo->j11 = jp[5]; ibminfo->k11 = kp[5];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[8]; ibminfo->j33 = jp[8]; ibminfo->k33 = kp[8];
      break;
    }
    case (5): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[3]; ibminfo->j22 = jp[3]; ibminfo->k22 = kp[3];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (6): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[6]; ibminfo->j33 = jp[6]; ibminfo->k33 = kp[6];
      break;
    }
    case (7): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[1]; ibminfo->j33 = jp[1]; ibminfo->k33 = kp[1];
      break;
    }
    }
    if (ibminfo->iimode >=0) {
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
    ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
    ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

    jp[0] = j+1; jp[3] = j+1; jp[6] = j+1;
    jp[1] = j+1; jp[4] = j+1; jp[7] = j+1;
    jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;

    kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
    kp[3] = k;   kp[4] = k;   kp[5] = k;
    kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

    for (nif=0; nif<9; nif++) {
      pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
      pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
      pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
    }

    ICP2(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
    switch (ibminfo->iimode) {
    case(0): {
      ibminfo->i11=ip[0]; ibminfo->j11 = jp[0]; ibminfo->k11 = kp[0];
      ibminfo->i22=ip[1]; ibminfo->j22 = jp[1]; ibminfo->k22 = kp[1];
      ibminfo->i33=ip[3]; ibminfo->j33 = jp[3]; ibminfo->k33 = kp[3];
      break;
    }
    case (1): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[2]; ibminfo->j22 = jp[2]; ibminfo->k22 = kp[2];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (2): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (3): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[5]; ibminfo->j22 = jp[5]; ibminfo->k22 = kp[5];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (4): {
      ibminfo->i11=ip[5]; ibminfo->j11 = jp[5]; ibminfo->k11 = kp[5];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[8]; ibminfo->j33 = jp[8]; ibminfo->k33 = kp[8];
      break;
    }
    case (5): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[3]; ibminfo->j22 = jp[3]; ibminfo->k22 = kp[3];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (6): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[6]; ibminfo->j33 = jp[6]; ibminfo->k33 = kp[6];
      break;
    }
    case (7): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[1]; ibminfo->j33 = jp[1]; ibminfo->k33 = kp[1];
      break;
    }
    }
    if (ibminfo->iimode >=0) {
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
    ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
    ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

    jp[0] = j-1; jp[3] = j  ; jp[6] = j+1;
    jp[1] = j-1; jp[4] = j  ; jp[7] = j+1;
    jp[2] = j-1; jp[5] = j  ; jp[8] = j+1;

    kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
    kp[3] = k-1; kp[4] = k-1; kp[5] = k-1;
    kp[6] = k-1; kp[7] = k-1; kp[8] = k-1;

    for (nif=0; nif<9; nif++) {
      pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
      pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
      pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
    }

    ICP2(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
    switch (ibminfo->iimode) {
    case(0): {
      ibminfo->i11=ip[0]; ibminfo->j11 = jp[0]; ibminfo->k11 = kp[0];
      ibminfo->i22=ip[1]; ibminfo->j22 = jp[1]; ibminfo->k22 = kp[1];
      ibminfo->i33=ip[3]; ibminfo->j33 = jp[3]; ibminfo->k33 = kp[3];
      break;
    }
    case (1): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[2]; ibminfo->j22 = jp[2]; ibminfo->k22 = kp[2];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (2): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (3): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[5]; ibminfo->j22 = jp[5]; ibminfo->k22 = kp[5];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (4): {
      ibminfo->i11=ip[5]; ibminfo->j11 = jp[5]; ibminfo->k11 = kp[5];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[8]; ibminfo->j33 = jp[8]; ibminfo->k33 = kp[8];
      break;
    }
    case (5): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[3]; ibminfo->j22 = jp[3]; ibminfo->k22 = kp[3];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (6): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[6]; ibminfo->j33 = jp[6]; ibminfo->k33 = kp[6];
      break;
    }
    case (7): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[1]; ibminfo->j33 = jp[1]; ibminfo->k33 = kp[1];
      break;
    }
    }
    if (ibminfo->iimode >=0) {
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }
    
    ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
    ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
    ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

    jp[0] = j-1; jp[3] = j  ; jp[6] = j+1;
    jp[1] = j-1; jp[4] = j  ; jp[7] = j+1;
    jp[2] = j-1; jp[5] = j  ; jp[8] = j+1;

    kp[0] = k+1; kp[1] = k+1; kp[2] = k+1;
    kp[3] = k+1; kp[4] = k+1; kp[5] = k+1;
    kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

    for (nif=0; nif<9; nif++) {
      pc[nif].x = coor[kp[nif]][jp[nif]][ip[nif]].x;
      pc[nif].y = coor[kp[nif]][jp[nif]][ip[nif]].y;
      pc[nif].z = coor[kp[nif]][jp[nif]][ip[nif]].z;
    }
    ICP2(p, pc, nfx, nfy, nfz, ibminfo, ip, jp, kp);
    switch (ibminfo->iimode) {
    case(0): {
      ibminfo->i11=ip[0]; ibminfo->j11 = jp[0]; ibminfo->k11 = kp[0];
      ibminfo->i22=ip[1]; ibminfo->j22 = jp[1]; ibminfo->k22 = kp[1];
      ibminfo->i33=ip[3]; ibminfo->j33 = jp[3]; ibminfo->k33 = kp[3];
      break;
    }
    case (1): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[2]; ibminfo->j22 = jp[2]; ibminfo->k22 = kp[2];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (2): {
      ibminfo->i11=ip[1]; ibminfo->j11 = jp[1]; ibminfo->k11 = kp[1];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[5]; ibminfo->j33 = jp[5]; ibminfo->k33 = kp[5];
      break;
    }
    case (3): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[5]; ibminfo->j22 = jp[5]; ibminfo->k22 = kp[5];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (4): {
      ibminfo->i11=ip[5]; ibminfo->j11 = jp[5]; ibminfo->k11 = kp[5];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[8]; ibminfo->j33 = jp[8]; ibminfo->k33 = kp[8];
      break;
    }
    case (5): {
      ibminfo->i11=ip[4]; ibminfo->j11 = jp[4]; ibminfo->k11 = kp[4];
      ibminfo->i22=ip[3]; ibminfo->j22 = jp[3]; ibminfo->k22 = kp[3];
      ibminfo->i33=ip[7]; ibminfo->j33 = jp[7]; ibminfo->k33 = kp[7];
      break;
    }
    case (6): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[7]; ibminfo->j22 = jp[7]; ibminfo->k22 = kp[7];
      ibminfo->i33=ip[6]; ibminfo->j33 = jp[6]; ibminfo->k33 = kp[6];
      break;
    }
    case (7): {
      ibminfo->i11=ip[3]; ibminfo->j11 = jp[3]; ibminfo->k11 = kp[3];
      ibminfo->i22=ip[4]; ibminfo->j22 = jp[4]; ibminfo->k22 = kp[4];
      ibminfo->i33=ip[1]; ibminfo->j33 = jp[1]; ibminfo->k33 = kp[1];
      break;
    }
    }
    if (ibminfo->iimode >=0) {
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

  DMDAVecRestoreArray(fda, user->lCent, &coor);
  
  return(0);
}




PetscErrorCode ibm_interpolation_advanced(UserCtx *user, 
					  IBMNodes *ibm, 
					  //  FSInfo *fsi,
					  PetscInt ibi,
					  PetscInt Add_dUndt)
{
  DM		da = user->da, fda = user->fda;
 
  Cmpnts	***ucont;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt nbn;
  PetscInt nbnumber = user->ibmnumber;
  PetscInt i,j,k;
  PetscReal sb, sc, lhs[3][3], rhs_l[3][3];
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscReal sk1, sk2, sk3, dtm, cv1, cv2, cv3, phia, phic;
  Cmpnts	***ucat;
  PetscReal	***nvert, ***nvert_o, ***p;
  PetscReal cs1, cs2, cs3;
  PetscInt	ni;
  PetscReal	ucx, ucy, ucz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  IBMInfo *ibminfo;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  PetscInt tmp;
  PetscReal nfx, nfy, nfz;
  PetscReal un, unold, unrm1, un_b, un_c;
  IBMListNode *current;

  for (tmp=0; tmp<9; tmp++) {
    
    DMDAVecGetArray(fda, user->Ucat, &ucat);
    DMDAVecGetArray(da, user->P, &p);
    
    PetscPrintf(PETSC_COMM_WORLD, "intp 011\n");
    
    current = user->ibmlist[ibi].head;
    while (current) {
      ibminfo = &current->ibm_intp;
      current = current->next;
      i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
      sb = ibminfo->d_s; sc = sb + ibminfo->d_i;

      lhs[0][0] = 1.;
      lhs[0][1] = - sb*sb;
      lhs[0][2] = -sb;
      lhs[1][0] = 0.;
      lhs[1][1] = sc*sc; 
      lhs[1][2] = sc;
      lhs[2][0] = -(sc-sb)/sc/sb+(1.-(sc-sb)/sc)/(sc-sb);
      lhs[2][1] = -2.*sb;
      lhs[2][2] = 1.;

      ip1 = ibminfo->i1; jp1 = ibminfo->j1; kp1 = ibminfo->k1;
      ip2 = ibminfo->i2; jp2 = ibminfo->j2; kp2 = ibminfo->k2;
      ip3 = ibminfo->i3; jp3 = ibminfo->j3; kp3 = ibminfo->k3;

      sk1  = ibminfo->cr1; sk2 = ibminfo->cr2; sk3 = ibminfo->cr3;
      cs1 = ibminfo->cs1; cs2 = ibminfo->cs2; cs3 = ibminfo->cs3;
      dtm = detmnt(lhs);
      rhs_l[0][1] = lhs[0][1];
      rhs_l[0][2] = lhs[0][2];
      rhs_l[1][1] = lhs[1][1];
      rhs_l[1][2] = lhs[1][2];
      rhs_l[2][1] = lhs[2][1];
      rhs_l[2][2] = lhs[2][2];

      ni = ibminfo->cell;

      nfx = ibm->nf_x[ni]; nfy = ibm->nf_y[ni]; nfz = ibm->nf_z[ni];

      PetscPrintf(PETSC_COMM_WORLD, "intp 012\n");
      if (ni>=0) {
	cv1 = ibm->u[ibm->nv1[ni]].x;
	cv2 = ibm->u[ibm->nv2[ni]].x;
	cv3 = ibm->u[ibm->nv3[ni]].x;
      }
      else {
	cv1 = 0;
	cv2 = 0;
	cv3=0;
      }
      /*       cv1 = 0; */
      /*       cv2 = 0; */
      /*       cv3=0; */

      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);
      un = phia * nfx;
      if (ni>=0)
	unold = (ibm->uold[ibm->nv1[ni]].x * cs1 +
		 ibm->uold[ibm->nv2[ni]].x * cs2 +
		 ibm->uold[ibm->nv3[ni]].x * cs3) * nfx;
      else
	unold =0;

      if (ni>=0)
	unrm1 = (ibm->urm1[ibm->nv1[ni]].x * cs1 +
		 ibm->urm1[ibm->nv2[ni]].x * cs2 +
		 ibm->urm1[ibm->nv3[ni]].x * cs3) * nfx;
      else
	unrm1 =0;

      cv1 = ucat[kp1][jp1][ip1].x;
      cv2 = ucat[kp2][jp2][ip2].x;
      cv3 = ucat[kp3][jp3][ip3].x;
      phic = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);

      if (invicid) {
	ucat[k][j][i].x = phic;	
      } else {
	rhs_l[0][0] = phia;
	rhs_l[1][0] = phic-phia;

	rhs_l[2][0] = -(sc-sb)/sc/sb*phia + (1.-(sc-sb)/sc)/(sc-sb)*phic;
	ucat[k][j][i].x = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc;//detmnt(rhs_l)/dtm; 
      }

      if (ni>=0) {
	cv1 = ibm->u[ibm->nv1[ni]].y;
	cv2 = ibm->u[ibm->nv2[ni]].y;
	cv3 = ibm->u[ibm->nv3[ni]].y;
      }
      else {
	cv1 = 0;
	cv2 = 0;
	cv3=0;
      }
   
      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);

      un += phia * nfy;
      if (ni>=0)
	unold += (ibm->uold[ibm->nv1[ni]].y * cs1 +
		  ibm->uold[ibm->nv2[ni]].y * cs2 +
		  ibm->uold[ibm->nv3[ni]].y * cs3) * nfy;
      else
	unold =0;
      if (ni>=0)
	unrm1 += (ibm->urm1[ibm->nv1[ni]].y * cs1 +
		  ibm->urm1[ibm->nv2[ni]].y * cs2 +
		  ibm->urm1[ibm->nv3[ni]].y * cs3) * nfy;
      else
	unrm1 =0;

      cv1 = ucat[kp1][jp1][ip1].y;
      cv2 = ucat[kp2][jp2][ip2].y;
      cv3 = ucat[kp3][jp3][ip3].y;
      phic = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);

      if (invicid){ 
	ucat[k][j][i].y = phic;	
      } else {
	rhs_l[0][0] = phia;
	rhs_l[1][0] = phic-phia;

	rhs_l[2][0] = -(sc-sb)/sc/sb*phia + (1.-(sc-sb)/sc)/(sc-sb)*phic;
	ucat[k][j][i].y = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc; //detmnt(rhs_l)/dtm;
      }

      if (ni>=0) {
	cv1 = ibm->u[ibm->nv1[ni]].z;
	cv2 = ibm->u[ibm->nv2[ni]].z;
	cv3 = ibm->u[ibm->nv3[ni]].z;
      }
      else {
	cv1 = 0;
	cv2 = 0;
	cv3 = 0;
      }
     
      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);

      un += phia * nfz;
      if (ni>=0)
	unold += (ibm->uold[ibm->nv1[ni]].z * cs1 +
		  ibm->uold[ibm->nv2[ni]].z * cs2 +
		  ibm->uold[ibm->nv3[ni]].z * cs3) * nfz;
      else
	unold =0;
      if (ni>=0)
	unrm1 += (ibm->urm1[ibm->nv1[ni]].z * cs1 +
		  ibm->urm1[ibm->nv2[ni]].z * cs2 +
		  ibm->urm1[ibm->nv3[ni]].z * cs3) * nfz;
      else
	unrm1 =0;

      cv1 = ucat[kp1][jp1][ip1].z;
      cv2 = ucat[kp2][jp2][ip2].z;
      cv3 = ucat[kp3][jp3][ip3].z;
      phic = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);

      if (invicid) {
	ucat[k][j][i].z = phic;		
      } else {
	rhs_l[0][0] = phia;
	rhs_l[1][0] = phic-phia;

	rhs_l[2][0] = -(sc-sb)/sc/sb*phia + (1.-(sc-sb)/sc)/(sc-sb)*phic;
	ucat[k][j][i].z = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc; //detmnt(rhs_l)/dtm;
      }

      if (invicid) {
	un_c= ucat[k][j][i].x*nfx +
	      ucat[k][j][i].y*nfy+
  	      ucat[k][j][i].z*nfz;
	un_b=un + (un_c-un)*sb/sc;
	ucat[k][j][i].x += un_b*nfx - un_c*nfx ;	
	ucat[k][j][i].y += un_b*nfy - un_c*nfy ;	
	ucat[k][j][i].z += un_b*nfz - un_c*nfz ;	
      }

      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];
      p[k][j][i] = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);
      
      p[k][j][i] += (un - unold) / user->dt * (sc-sb);

    }
    
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
  }

  return(0);  
}


PetscErrorCode ibm_prs_span(UserCtx *user, 
			    IBMNodes *ibm,
			    PetscInt ibi)
{
  DM		da = user->da, fda = user->fda;
 
  Cmpnts	***ucont;

  PetscInt nbn;
  PetscInt nbnumber = user->ibmnumber;
  PetscInt i,j,k;
  PetscReal sb, sc, lhs[3][3], rhs_l[3][3];
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscReal sk1, sk2, sk3, dtm, cv1, cv2, cv3, phia, phic;
  Cmpnts	***ucat;
  PetscReal	***nvert, ***nvert_o, ***p,p_inf;
  PetscReal cs1, cs2, cs3;
  PetscInt	ni,nii;

  Cmpnts        pmin, span_b[5], span_e[5], span[5];//, p_max[5], p_min[5];
  PetscReal     spanlength[5], spanloc, p_a, p_b, p_c, angvel=-0.361;

  Cmpnts ***coor;

  IBMInfo *ibminfo;

  
  PetscInt tmp;
  PetscReal nfx, nfy, nfz;
 
  IBMListNode *current;
  PetscInt    J_p[5];

  PetscInt rank;
  FILE            *f;
  char            filen[80];
  
  PetscPrintf(PETSC_COMM_WORLD, "intp 001\n");
  for (tmp=0; tmp<1; tmp++) {

    DMDAVecGetArray(da, user->P, &p);
    p_inf=p[2][71][71];
    DMDAVecRestoreArray(da, user->P, &p);
    VecShift(user->P, p_inf);

    DMDAVecGetArray(fda, user->Cent, &coor);

    DMDAVecGetArray(da, user->P, &p);    
    PetscPrintf(PETSC_COMM_WORLD, "intp 011\n");

    current = user->ibmlist[ibi].head;

    J_p[0]=140;//21;//254;
    J_p[1]=131;//30;//236;
    J_p[2]=120;//41;//217;
    J_p[3]=110;//51;//199;
    J_p[4]=100;//61;//179;
    span_b[0].z=0.   ;span_b[0].y=15.57;span_b[0].x=0.379;span_e[0].x=-0.803;span_e[0].y=14.25;span_e[0].z=-0.031;
    span_b[1].z=0.   ;span_b[1].y=13.11;span_b[1].x=0.449;span_e[1].x=-1.049;span_e[1].y=12.00;span_e[1].z=-0.00996;
    span_b[2].z=0.011;span_b[2].y=10.33;span_b[2].x=0.534;span_e[2].x=-1.246;span_e[2].y= 9.45;span_e[2].z= 0.0464;
    span_b[3].z=0.051;span_b[3].y= 7.7 ;span_b[3].x=0.615;span_e[3].x=-1.433;span_e[3].y= 7.05;span_e[3].z= 0.2196;
    span_b[4].z=0.173;span_b[4].y= 4.92;span_b[4].x=0.677;span_e[4].x=-1.582;span_e[4].y= 4.5 ;span_e[4].z= 0.748;

    for (ni=0; ni<5; ni++) {

      //scale span
      span_b[ni].x /=1.0928;
      span_b[ni].y /=1.0928;
      span_b[ni].z /=1.0928;
      span_e[ni].x /=1.0928;
      span_e[ni].y /=1.0928;
      span_e[ni].z /=1.0928;

      span[ni].x=span_e[ni].x - span_b[ni].x;
      span[ni].y=span_e[ni].y - span_b[ni].y;
      span[ni].z=span_e[ni].z - span_b[ni].z;
      spanlength[ni]=(span[ni].x * span[ni].x +
		      span[ni].z * span[ni].z );
    }

    while (current) {
      ibminfo = &current->ibm_intp;
      current = current->next;
      i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
      sb = ibminfo->d_s; sc = sb + ibminfo->d_i;
    
      ip1 = ibminfo->i1; jp1 = ibminfo->j1; kp1 = ibminfo->k1;
      ip2 = ibminfo->i2; jp2 = ibminfo->j2; kp2 = ibminfo->k2;
      ip3 = ibminfo->i3; jp3 = ibminfo->j3; kp3 = ibminfo->k3;

      sk1  = ibminfo->cr1; sk2 = ibminfo->cr2; sk3 = ibminfo->cr3;
         
      pmin = ibminfo->pmin;

      for (ni=0; ni<5; ni++) {
	if (J_p[ni]==j) {
	  // calculate x/c 
	  PetscPrintf(PETSC_COMM_WORLD, "location %le %le \n", pmin.y, pmin.y/15*100); 
	  spanloc= ((pmin.x - span_b[ni].x)*span[ni].x + 
		    (pmin.z - span_b[ni].z)*span[ni].z)/spanlength[ni];

	  cv1 = p[kp1][jp1][ip1];
	  cv2 = p[kp2][jp2][ip2];
	  cv3 = p[kp3][jp3][ip3];

	  p_c = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);
	  p_b=p[k][j][i];
	  p_a=p_c-(p_c-p_b)*sc/(sc-sb);//+(span_b[ni].y*angvel*angvel*
	  //(-pmin.y*nfy-pmin.x*nfx)*sb);   //pressure on surf elmt ni   
	  
	  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	  
	  if (!rank) {
	    sprintf(filen, "PressSpan%3.3d.dat",ni);
	    f = fopen(filen, "a");
	    
	    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le %le %d %d %d\n", p_a, spanloc, pmin.x, pmin.y, pmin.z, i,J_p[ni],k);
	    fclose(f);	
	  }
	}
      }

    }

    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(fda, user->Cent, &coor);
  
  }
 
  return(0);  
}

PetscErrorCode ibm_power(IBMNodes *ibm,
			 PetscInt ti)
{

  PetscInt   elmt, n_elmt=ibm->n_elmt;
  PetscInt   nv1,nv2,nv3;
  PetscReal  M_z, p, dA;
  PetscReal  n_x,n_y,n_z, F_x, F_y,F_z,Pw,Pw_y,r_x, r_y;
  PetscReal  u_x,u_y,u_z;
  PetscReal  F_xSum,F_ySum,F_zSum,A; //Surface Force

  M_z= 0.; Pw=0;  Pw_y=0.;
  F_xSum=0.; F_ySum=0.; F_zSum=0.;

  for (elmt=0; elmt<n_elmt; elmt++) {
    nv1=ibm->nv1[elmt];
    nv2=ibm->nv2[elmt];
    nv3=ibm->nv3[elmt];
    
    r_x= ibm->cent_x[elmt];
    r_y= ibm->cent_y[elmt];
    
    dA = ibm->dA[elmt];

    n_x= ibm->nf_x[elmt];
    n_y= ibm->nf_y[elmt];
    n_z= ibm->nf_z[elmt];

    p = (ibm->nt_y[elmt]);

    F_x = -p*dA*n_x;
    F_y = -p*dA*n_y;
    F_z = -p*dA*n_z;

    F_xSum += F_x;
    F_ySum += F_y;
    F_zSum += F_z;

    u_x = (ibm->u[nv1].x + ibm->u[nv2].x + ibm->u[nv3].x)/3.;
    u_y = (ibm->u[nv1].y + ibm->u[nv2].y + ibm->u[nv3].y)/3.;
    u_z = (ibm->u[nv1].z + ibm->u[nv2].z + ibm->u[nv3].z)/3.;

    Pw += F_x*u_x + F_y*u_y + F_z*u_z;    
    Pw_y += F_y*u_y;    
  }

  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "power_projected.dat");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le\n", ti,Pw, Pw_y,F_xSum, F_ySum, F_zSum);
    fclose(f);
  }
  PetscPrintf(PETSC_COMM_WORLD, "The Power %le forces %le %le %le\n", Pw, F_xSum, F_ySum, F_zSum); 

  return(0); 
}

PetscErrorCode ibm_prs_projection(UserCtx *user, 
				  IBMNodes *ibm,
				  PetscInt ibi)
{
  DM		da = user->da, fda = user->fda;
 
  Cmpnts	***ucont;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt nbn;
  PetscInt nbnumber = user->ibmnumber;
  PetscInt i,j,k;
  PetscReal sb, sc, lhs[3][3], rhs_l[3][3];
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscReal sk1, sk2, sk3, dtm, cv1, cv2, cv3, phia, phic;
  Cmpnts	***ucat;
  PetscReal	***nvert, ***nvert_o, ***p,p_inf;
  PetscReal cs1, cs2, cs3;
  PetscInt	ni,nii;
  PetscReal	ucx, ucy, ucz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  IBMInfo *ibminfo;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  PetscInt tmp;
  PetscReal nfx, nfy, nfz;
  PetscReal un, unold, unrm1, un_b, un_c, p_min[my], p_max[my];
  PetscReal c_min[my], c_max[my];
  IBMListNode *current;
  PetscInt   elmt;
  Cmpnts     ***cent;
  FILE            *f;
  char            filen[80];
  
  DMDAVecGetArray(user->fda, user->Cent, &cent)  ;
  
  for (j=ys; j<ye; j++) {
    p_min[j]=1.e20;
    p_max[j]=-1.e20;
    c_min[j]=1.e20;
    c_max[j]=-1.e20;
  }
  for (elmt=0; elmt<ibm->n_elmt; elmt++) {
    ibm->nt_x[elmt]=0.;
    ibm->nt_y[elmt]=0.;
    ibm->nt_z[elmt]=0.;
    ibm->ns_x[elmt]=0.;
  }

 
  for (tmp=0; tmp<1; tmp++) {

    DMDAVecGetArray(da, user->P, &p);    
   
    current = user->ibmlist[ibi].head;
    while (current) {
      ibminfo = &current->ibm_intp;
      current = current->next;
      i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
      sb = ibminfo->d_s; sc = sb + ibminfo->d_i;

      ip1 = ibminfo->i1; jp1 = ibminfo->j1; kp1 = ibminfo->k1;
      ip2 = ibminfo->i2; jp2 = ibminfo->j2; kp2 = ibminfo->k2;
      ip3 = ibminfo->i3; jp3 = ibminfo->j3; kp3 = ibminfo->k3;

      sk1  = ibminfo->cr1; sk2 = ibminfo->cr2; sk3 = ibminfo->cr3;
      cs1  = ibminfo->cs1; cs2 = ibminfo->cs2; cs3 = ibminfo->cs3;

      ni = ibminfo->cell;

      nfx = ibm->nf_x[ni]; nfy = ibm->nf_y[ni]; nfz = ibm->nf_z[ni];

      if (p[k][j][i]<p_min[j]) p_min[j]=p[k][j][i];
      if (p[k][j][i]>p_max[j]) p_max[j]=p[k][j][i];

      if (cent[k][j][i].x<c_min[j]) c_min[j]=cent[k][j][i].x;
      if (cent[k][j][i].x>c_max[j]) c_max[j]=cent[k][j][i].x;

      ibm->nt_x[ni]+=1.; //number of ibm nodes close to surf elmt ni
      ibm->nt_y[ni]= (p[k][j][i]+ibm->nt_y[ni]*(ibm->nt_x[ni]-1.))/ibm->nt_x[ni];   //pressure on surf elmt ni   
    }


    DMDAVecRestoreArray(da, user->P, &p);
  }
  
  for (ni=0; ni<ibm->n_elmt; ni++){
    if (fabs(ibm->nt_x[ni])<0.5) {
     
      cs1=0;
      for (nii=0; nii<ibm->n_elmt; nii++){
	// if different elements & on the same side n1.n2>0
	if (nii!=ni && ((ibm->nf_x[ni]*ibm->nf_x[nii] +
			 ibm->nf_y[ni]*ibm->nf_y[nii] +
			 ibm->nf_z[ni]*ibm->nf_z[nii])>0.)) { 
	  if (ibm->nv1[ni]==ibm->nv1[nii] || ibm->nv1[ni]==ibm->nv2[nii] || ibm->nv1[ni]==ibm->nv3[nii] ||
	      ibm->nv2[ni]==ibm->nv1[nii] || ibm->nv2[ni]==ibm->nv2[nii] || ibm->nv2[ni]==ibm->nv3[nii] ||
	      ibm->nv3[ni]==ibm->nv1[nii] || ibm->nv3[ni]==ibm->nv2[nii] || ibm->nv3[ni]==ibm->nv3[nii]) {
	    cs1=cs1+1.;
	    
	    ibm->nt_y[ni]= (ibm->nt_y[nii]+ibm->nt_y[ni]*(cs1-1.))/cs1;
	  }
	}
      }
    }
  }

  // smoothing
  for(i=0;i<2;i++){
    for (ni=0; ni<ibm->n_elmt; ni++){
      if (fabs(ibm->nt_x[ni])<0.5) {
	
	cs1=0;
	for (nii=0; nii<ibm->n_elmt; nii++){
	  if (nii!=ni) {
	    if (ibm->nv1[ni]==ibm->nv1[nii] || ibm->nv1[ni]==ibm->nv2[nii] || ibm->nv1[ni]==ibm->nv3[nii] ||
		ibm->nv2[ni]==ibm->nv1[nii] || ibm->nv2[ni]==ibm->nv2[nii] || ibm->nv2[ni]==ibm->nv3[nii] ||
		ibm->nv3[ni]==ibm->nv1[nii] || ibm->nv3[ni]==ibm->nv2[nii] || ibm->nv3[ni]==ibm->nv3[nii]) {
	      cs1=cs1+1.;
	      ibm->nt_y[ni]= (ibm->nt_y[nii]+ibm->nt_y[ni]*(cs1-1.))/cs1;
	    }
	  }
	}
      }
    }
  }
  
  return(0);
}

PetscErrorCode ibm_interpolation_advanced2(UserCtx *user, IBMNodes *ibm)
{
  DM		da = user->da, fda = user->fda;
 
  Cmpnts	***ucont;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***icsi, ***jeta, ***kzet;

  PetscInt nbn;
  PetscInt nbnumber = user->ibmnumber;
  PetscInt i,j,k;
  PetscReal sb, sc, lhs[3][3], rhs_l[3][3];
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscInt	ip11, ip22, ip33, jp11, jp22, jp33, kp11, kp22, kp33;
  PetscReal sk1, sk2, sk3, dtm, cv1, cv2, cv3, phia, phic;
  PetscReal sk11,sk22,sk33, cv11,cv22,cv33;
  Cmpnts	***ucat;
  PetscReal	***nvert, ***nvert_o, ***p;
  PetscReal cs1, cs2, cs3;
  PetscInt	ni;
  PetscReal	ucx, ucy, ucz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  IBMInfo *ibminfo;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  PetscInt tmp;
  PetscReal nfx, nfy, nfz;
  PetscReal un, unold;
  IBMListNode *current;
  for (tmp=0; tmp<9; tmp++) {

    DMDAVecGetArray(fda, user->lUcat, &ucat);
    DMDAVecGetArray(da, user->lP, &p);

    DMDAVecGetArray(fda, user->Ucont, &ucont);
    DMDAVecGetArray(fda, user->lICsi, &icsi);
    DMDAVecGetArray(fda, user->lJEta, &jeta);
    DMDAVecGetArray(fda, user->lKZet, &kzet);
    DMDAVecGetArray(da, user->lNvert, &nvert);
    DMDAVecGetArray(da, user->lNvert_o, &nvert_o);

    //current = user->ibmlist.head;
    while (current) {
      ibminfo = &current->ibm_intp;
      current = current->next;
      i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
      sb = ibminfo->d_s; sc = sb + ibminfo->d_i;
      lhs[0][0] = 1.;
      lhs[0][1] = - sb*sb;
      lhs[0][2] = -sb;
      lhs[1][0] = 0.;
      lhs[1][1] = sc*sc;
      lhs[1][2] = -sc;
      lhs[2][0] = ((sc-sb)*(sc-sb) - sb*sb)/(sc-sb);
      lhs[2][1] = -2 * sc * sb*sb;
      lhs[2][2] = -sb*sc;

          
      ip1 = ibminfo->i1; jp1 = ibminfo->j1; kp1 = ibminfo->k1;
      ip2 = ibminfo->i2; jp2 = ibminfo->j2; kp2 = ibminfo->k2;
      ip3 = ibminfo->i3; jp3 = ibminfo->j3; kp3 = ibminfo->k3;

      ip11 = ibminfo->i11; jp11 = ibminfo->j11; kp11 = ibminfo->k11;
      ip22 = ibminfo->i22; jp22 = ibminfo->j22; kp22 = ibminfo->k22;
      ip33 = ibminfo->i33; jp33 = ibminfo->j33; kp33 = ibminfo->k33;

      sk1  = ibminfo->cr1 ; sk2 = ibminfo->cr2 ; sk3 = ibminfo->cr3 ;
      sk11 = ibminfo->cr11; sk22= ibminfo->cr22; sk33= ibminfo->cr33;
      cs1 = ibminfo->cs1; cs2 = ibminfo->cs2; cs3 = ibminfo->cs3;
      dtm = detmnt(lhs);
      rhs_l[0][1] = lhs[0][1];
      rhs_l[0][2] = lhs[0][2];
      rhs_l[1][1] = lhs[1][1];
      rhs_l[1][2] = lhs[1][2];
      rhs_l[2][1] = lhs[2][1];
      rhs_l[2][2] = lhs[2][2];

      ni = ibminfo->cell;

      nfx = ibm->nf_x[ni]; nfy = ibm->nf_y[ni]; nfz = ibm->nf_z[ni];

      if (ni>=0) {
	cv1 = ibm->u[ibm->nv1[ni]].x;
	cv2 = ibm->u[ibm->nv2[ni]].x;
	cv3 = ibm->u[ibm->nv3[ni]].x;
      }
      else {
	cv1 = 0;
	cv2 = 0;
	cv3=0;
      }
     
      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);
      un = phia * nfx;
      
      cv1 = ucat[kp1][jp1][ip1].x;
      cv2 = ucat[kp2][jp2][ip2].x;
      cv3 = ucat[kp3][jp3][ip3].x;

      cv11 = ucat[kp11][jp11][ip11].x;
      cv22 = ucat[kp22][jp22][ip22].x;
      cv33 = ucat[kp33][jp33][ip33].x;

      phic = 0.5*(cv1 * sk1 + cv2 * sk2 + cv3 * sk3 + 
		  cv11*sk11 + cv22*sk22 + cv33*sk33) ;

      rhs_l[0][0] = phia;
      rhs_l[1][0] = phic-phia;
      rhs_l[2][0] = (sc-sb) * phia - (sb*sb/(sc-sb)) * phic;
      ucat[k][j][i].x = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc; 


      if (ni>=0) {
	cv1 = ibm->u[ibm->nv1[ni]].y;
	cv2 = ibm->u[ibm->nv2[ni]].y;
	cv3 = ibm->u[ibm->nv3[ni]].y;
      }
      else {
	cv1 = 0;
	cv2 = 0;
	cv3=0;
      }
     
      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);

      un += phia * nfy;

      cv1 = ucat[kp1][jp1][ip1].y;
      cv2 = ucat[kp2][jp2][ip2].y;
      cv3 = ucat[kp3][jp3][ip3].y;

      cv11 = ucat[kp11][jp11][ip11].y;
      cv22 = ucat[kp22][jp22][ip22].y;
      cv33 = ucat[kp33][jp33][ip33].y;

      phic = 0.5*(cv1 * sk1 + cv2 * sk2 + cv3 * sk3 + 
		  cv11*sk11 + cv22*sk22 + cv33*sk33) ;

      rhs_l[0][0] = phia;
      rhs_l[1][0] = phic-phia;
      rhs_l[2][0] = (sc-sb) * phia - (sb*sb/(sc-sb)) * phic;
      ucat[k][j][i].y = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc; 


      if (ni>=0) {
	cv1 = ibm->u[ibm->nv1[ni]].z;
	cv2 = ibm->u[ibm->nv2[ni]].z;
	cv3 = ibm->u[ibm->nv3[ni]].z;
      }
      else {
	cv1 = 0;
	cv2 = 0;
	cv3=0;
      }
     
      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);

      un += phia * nfz;

      cv1 = ucat[kp1][jp1][ip1].z;
      cv2 = ucat[kp2][jp2][ip2].z;
      cv3 = ucat[kp3][jp3][ip3].z;

      cv11 = ucat[kp11][jp11][ip11].z;
      cv22 = ucat[kp22][jp22][ip22].z;
      cv33 = ucat[kp33][jp33][ip33].z;

      phic = 0.5*(cv1 * sk1 + cv2 * sk2 + cv3 * sk3 + 
		  cv11*sk11 + cv22*sk22 + cv33*sk33) ;

      rhs_l[0][0] = phia;
      rhs_l[1][0] = phic-phia;
      rhs_l[2][0] = (sc-sb) * phia - (sb*sb/(sc-sb)) * phic;
      ucat[k][j][i].z = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc; 

      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];

      cv11 = p[kp11][jp11][ip11];
      cv22 = p[kp22][jp22][ip22];
      cv33 = p[kp33][jp33][ip33];

      p[k][j][i] = 0.5*(cv1 * sk1 + cv2 * sk2 + cv3 * sk3 + 
			cv11*sk11 + cv22*sk22 + cv33*sk33) ;

    }

    DMDAVecRestoreArray(fda, user->lUcat, &ucat);
    DMDAVecRestoreArray(da, user->lP, &p);

    DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
    DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

    DMDALocalToLocalBegin(da, user->lP, INSERT_VALUES, user->lP);
    DMDALocalToLocalEnd(da, user->lP, INSERT_VALUES, user->lP);

    PetscBarrier(PETSC_NULL);
    DMDAVecGetArray(fda, user->lUcat, &ucat);
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if ((int)(nvert[k][j][i]+0.5) ==1) {// || (int)(nvert_o[k][j][i]+0.5)==1) {
	    ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
	    ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
	    ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
	    ucont[k][j][i].x = ucx * icsi[k][j][i].x +
	      ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z;
	  
	    ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
	    ucont[k][j][i].y = ucx * jeta[k][j][i].x +
	      ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z;
	  
	    ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
	    ucont[k][j][i].z = ucx * kzet[k][j][i].x +
	      ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z;
	  }

	  if ((int)(nvert[k][j][i+1]+0.5)==1) {
	    ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
	    ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
	    ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
	    ucont[k][j][i].x = ucx * icsi[k][j][i].x +
	      ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z;
	  }
	
	  if ((int)(nvert[k][j+1][i]+0.5)==1) {
	    ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
	    ucont[k][j][i].y = ucx * jeta[k][j][i].x +
	      ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z;
	  }

	  if ((int)(nvert[k+1][j][i]+0.5)==1) {
	    ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
	    ucont[k][j][i].z = ucx * kzet[k][j][i].x +
	      ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z;
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->lUcat, &ucat);
    DMDAVecRestoreArray(fda, user->Ucont, &ucont);
    DMDAVecRestoreArray(fda, user->lICsi, &icsi);
    DMDAVecRestoreArray(fda, user->lJEta, &jeta);
    DMDAVecRestoreArray(fda, user->lKZet, &kzet);
    DMDAVecRestoreArray(da, user->lNvert, &nvert);
    DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);

    PetscBarrier(PETSC_NULL);

    DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

    Contra2Cart(user, user->lUcont, user->Ucat);
  }

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

}

PetscErrorCode ibm_interpolation_advanced_mg(UserCtx *user, IBMNodes *ibm, Vec X)
{
  DM		da = user->da, fda = user->fda;
 
  Cmpnts	***ucont;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***icsi, ***jeta, ***kzet;

  PetscInt nbn;
  PetscInt nbnumber = user->ibmnumber;
  PetscInt i,j,k;
  PetscReal sb, sc, lhs[3][3], rhs_l[3][3];
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscReal sk1, sk2, sk3, dtm, cv1, cv2, cv3, phia, phic;
  Cmpnts	***ucat;
  PetscReal	***nvert, ***nvert_o, ***p;
  PetscReal cs1, cs2, cs3;
  PetscInt	ni;
  PetscReal	ucx, ucy, ucz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  IBMInfo *ibminfo;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


  PetscInt tmp;
  IBMListNode *current;

  Vec lX;

  DMCreateLocalVector(da, &lX);

  DMGlobalToLocalBegin(da, X, INSERT_VALUES, lX);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, lX);
  DMDAVecGetArray(da, lX, &p);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  for (tmp=0; tmp<5; tmp++) {

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  PetscReal wa;
	  PetscInt na;
	  wa = 0;
	  na = 0;
	  if (nvert[k][j][i] > 0.1) {

	    if (nvert[k][j-1][i] <0.1 && j>1) {
	      wa += p[k][j-1][i];
	      na ++;	     
	    }			       
	    if (nvert[k][j+1][i] < 0.1 && j<my-2) {
	      wa += p[k][j+1][i];
	      na ++;
	    }
	    if (nvert[k][j][i-1] <0.1 && i>1) {
	      wa += p[k][j][i-1];
	      na ++;	        
	    }		   	       
	    if (nvert[k][j][i+1] < 0.1 && i<mx-2) {
	      wa += p[k][j][i+1];
	      na ++;
	    }
	    if (na) {
	      wa /= na;
	    }
	    p[k][j][i] = wa;
	  }
	}
      }
    }
  }
  DMDAVecRestoreArray(da, lX, &p);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMLocalToGlobalBegin(da, lX, INSERT_VALUES, X);
  DMLocalToGlobalEnd(da, lX, INSERT_VALUES, X);

  VecDestroy(&lX);
  return 0;
  
}
