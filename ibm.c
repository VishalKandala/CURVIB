#include "variables.h"
#include "stdlib.h"
#include "time.h"
extern PetscInt thin, block_number,invicid, sediment;
extern PetscInt NumberOfBodies, moveframe, wallfunction;

PetscErrorCode BoundingSphere(IBMNodes *ibm);
PetscErrorCode nearestcell1(Cmpnts p, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode InterceptionPoint(Cmpnts p, PetscInt i, PetscInt j, PetscInt k,
				 IBMInfo *ibminfo, UserCtx *user);
PetscInt point_cell_thin(Cmpnts p,Cmpnts p1,Cmpnts p2,Cmpnts p3,Cmpnts p4,
			 PetscInt ip, PetscInt jp, PetscInt kp,
			 IBMNodes *ibm, PetscInt ncx, PetscInt ncy,
			 PetscInt ncz, PetscReal dcx, PetscReal dcy,
			 PetscReal xbp_min, PetscReal ybp_min,
			 PetscReal zbp_max, List *cell_trg,
			 PetscInt flg);

PetscInt point_cell_advanced(Cmpnts p, PetscInt ip, PetscInt jp, PetscInt kp,
			     IBMNodes *ibm, PetscInt ncx, PetscInt ncy,
			     PetscInt ncz, PetscReal dcx, PetscReal dcy,
			     PetscReal xbp_min, PetscReal ybp_min,
			     PetscReal zbp_max, List *cell_trg,
			     PetscInt flg);

void wall_function_loglaw (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			   double nx, double ny, double nz);

void wall_function_Cabot (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			   double nx, double ny, double nz, double dpdx, double dpdy, double dpdz,int counter);

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

PetscErrorCode ibm_search_advanced(UserCtx *user, IBMNodes *ibm, 
				   PetscInt ibi)

/*      Note : Always go from ibi (immersed body number) 0 -> NumberOfBodies  */
/*             Nvert should be set to zero before any new search and this */
/*             happens if ibi==0--not anymore! set nvert=0 manually!*/
{
  DM	da = user->da, fda = user->fda;
  // mx,my,mz -- global number of grid points in each direction.
  // xm,ym,zm -- number of grid  points on this processor,excluding ghosts.
  // xe,ye,ze -- end point of  this processor in each direction,excluding ghost nodes.
  // xs,ys,zs -- start point of this processor in each direction,excluding ghost nodes.
  // gxs,gys.gzs -- start point of this processor including ghosts.
  // gxm,gym,gzm --  number of grid points on this processor including ghosts.
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

  PetscInt      gxs, gxe, gys, gye, gzs, gze, gxm,gym,gzm; 
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  gxm = info.gxm; 
  gym = info.gym; 
  gzm = info.gzm; 

  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;

  for(i=0; i<n_v; i++) {
//  The minimum and maximum  co-ordinates of the immersed boundary are determined.(Bounding box is determined)
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);

    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);

    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }
 
// Tolerances to bounding box are added.
  xbp_min -= 0.05; xbp_max += 0.05;
  ybp_min -= 0.05; ybp_max += 0.05;
  zbp_min -= 0.05; zbp_max += 0.05;
 
// Control cell size  in each direction is determined (ncx,ncy,ncz are the number of control cells in each direction).
  dcx = (xbp_max - xbp_min) / (ncx - 1.);
  dcy = (ybp_max - ybp_min) / (ncy - 1.);
  dcz = (zbp_max - zbp_min) / (ncz - 1.);

/*   PetscPrintf(PETSC_COMM_WORLD, "zbp min max %le %le\n",zbp_min,zbp_max); */
  PetscMalloc(ncz * ncy * ncx * sizeof(List), &cell_trg);   // A list object  called  cell_trg with each control cell as an element is created,
  
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);  //  The cell_trg list object is initialized.
      }
    }
  }

  for (ln_v=0; ln_v < ibm->n_elmt; ln_v++) {
// Going through each element of the  immersed boundary, determine what control cell each element overlaps with.
   // For each element, the min,max co-ordinate is identified in each direction.
    n1e = ibm->nv1[ln_v]; n2e = ibm->nv2[ln_v]; n3e = ibm->nv3[ln_v];

    xv_min = PetscMin(PetscMin(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);
    xv_max = PetscMax(PetscMax(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);

    yv_min = PetscMin(PetscMin(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    yv_max = PetscMax(PetscMax(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);

    zv_min = PetscMin(PetscMin(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    zv_max = PetscMax(PetscMax(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    
    // The control cell that overlaps with min x,y,z for each element in immersed boundary as well as that which overlaps with max x.y,z are identified.
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

   // All the control cells between the min and max in all directions are noted to be overlapping with the given element of the immersed boundary.
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
 // Going through all the cell centers(excluding boundary points which are at surface centers) 
 for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (ibi==0 && nvert[k][j][i]<5.) nvert[k][j][i] = 0; //reset nvert if new search

	if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&     
	    coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	    coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {            // Checking whether the grid point(fluid grid) is inside the bounding box.

	  // identifying  the control cell that the fluid grid point coincides with.
          ic = floor((coor[k][j][i].x - xbp_min )/ dcx);   
	  jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	  kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	  
	  nvert[k][j][i] =  PetscMax(nvert[k][j][i],
				     point_cell_advanced(coor[k][j][i], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg)); // Setting nvert.
	  if (nvert[k][j][i] < 0) nvert[k][j][i] = 0;
	}
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
	    
	    cutthrough = point_cell_thin(coor[k][j][i],coor[k][j][i+1],
					 coor[k][j+1][i],coor[k+1][j][i],
					 coor[k+1][j+1][i+1], ic, jc, kc, 
					 ibm, ncx, ncy, ncz, dcx, dcy, 
					 xbp_min, ybp_min, zbp_max, cell_trg, flg);

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
  //
  // Considering all the points,check whether it is classified as  a fluid node (i.e nvert = 0),  then  check if any of it's neighbours are classified as solid(nvert=4), then set this point to be a boundary point(nvert = 2).
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	//	if (ibi>0 && ibi<4 && i==39 && j==74 && k==124)   PetscPrintf(PETSC_COMM_SELF, "nvert [%d][%d][%d] %le \n", i, j, k,nvert[k][j][i]);
	if (nvert[k][j][i] <0) nvert[k][j][i] = 0;
	ip = (i<mx-1?(i+1):(i));  // neighbouring points i+1,i,i-1
	im = (i>0   ?(i-1):(i));

	jp = (j<my-1?(j+1):(j));  // neighbouring points  j+1,j,j-1 
	jm = (j>0   ?(j-1):(j));

	kp = (k<mz-1?(k+1):(k));  // neighbouring points k+1,k,k-1
	km = (k>0   ?(k-1):(k));

	if ((int)(nvert[k][j][i]+0.5) != 4) {  // check if i,j,k is  fluid.
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(nvert[kk][jj][ii] +0.5) == 4) {  // if any neighbour is solid.
		  nvert[k][j][i] = PetscMax(2, nvert[k][j][i]); // set i,j,k as boundary point.
		}
	      }
	    }
	  }
	}
      }
    }
  }

  PetscBarrier(PETSC_NULL);

// If there is a fluid node with boundaries around it in more than three directions (+x,-x,+y,-y,+z,-z), set it also to be a boundary point.
  PetscInt around, tmp;
  for (tmp=0;tmp<2; tmp++) {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] <0) nvert[k][j][i] = 0;
	  ip = (i<mx-1?(i+1):(i));
	  im = (i>0   ?(i-1):(i));

	  jp = (j<my-1?(j+1):(j));
	  jm = (j>0   ?(j-1):(j));
	  
	  kp = (k<mz-1?(k+1):(k));
	  km = (k>0   ?(k-1):(k));
	  
	  around=0;
	  if ((int)(nvert[k][j][i]+0.5) <0.5) {
	    for (kk=km; kk<kp+1; kk++)
	      if ((int)(nvert[kk][j][i] +0.5) >1 && kk!=k)
		around++;

	    for (jj=jm; jj<jp+1; jj++)
	      if ((int)(nvert[k][jj][i] +0.5) >1 && jj!=j)
		around++;

	    for (ii=im; ii<ip+1; ii++)
	      if ((int)(nvert[k][j][ii] +0.5) >1 && ii!=i)
		around++;
	  }

	  if (around>3) {
      //  PetscPrintf(PETSC_COMM_SELF, "Around 5! ijk %d %d %d nvert %e\n",i,j,k, nvert[k][j][i]);
	    nvert[k][j][i] = PetscMax(2, nvert[k][j][i]);
	  }
	}
      }
    }
  }  

  // If there is a massive change in nvert at a location, point out phase change.
  PetscReal	***nvert_o;
  DMDAVecGetArray(da, user->lNvert_o, &nvert_o);
  if (ibi==NumberOfBodies-1)
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert_o[k][j][i] >2.5 && nvert[k][j][i] < 0.5) {
	  PetscPrintf(PETSC_COMM_SELF, "Phase Change at %d, %d, %d! nvert_o %le nvert %le \n", i, j, k,nvert_o[k][j][i],nvert[k][j][i]);
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

  BoundingSphere(ibm);  // Find the radius and center of  bounding spheres for each element(qvec,radvec). 

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) {  // For each boundary point.
	  number ++;
          // Find the control cell associated with the grid point.
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
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) nvert[k][j][i]=1;
	if ((int)(nvert[k][j][i]+0.5) == 4) nvert[k][j][i]=3;
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
/*     PetscPrintf(PETSC_COMM_WORLD, "%i\n", i); */
/*     PetscPrintf(PETSC_COMM_WORLD, "%e\n", x_bp[i]); */
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
  //  zbp_min -= 0.09; zbp_max += 0.09;

  dcx = (xbp_max - xbp_min) / (ncx - 1.);
  dcy = (ybp_max - ybp_min) / (ncy - 1.);
  dcz = (zbp_max - zbp_min) / (ncz - 1.);

/*   PetscPrintf(PETSC_COMM_WORLD, "zbp min max %le %le\n",zbp_min,zbp_max); */
  PetscMalloc(ncz * ncy * ncx * sizeof(List), &cell_trg);
  //PetscPrintf(PETSC_COMM_SELF, "test00\n");
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);
      }
    }
  }
  //PetscPrintf(PETSC_COMM_SELF, "test0\n");
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

    /*    if (ln_v==25) {
      PetscPrintf(PETSC_COMM_SELF, "id25, %d %d %d %d %d %d\n", iv_min, iv_max, jv_min, jv_max, kv_min, kv_max);
      }*/
    // Insert IBM node information into a list
    for (k=kv_min; k<kv_max; k++) {
      for (j=jv_min; j<jv_max; j++) {
	for (i=iv_min; i<iv_max; i++) {
	  insertnode(&(cell_trg[k *ncx*ncy + j*ncx +i]), ln_v);
	}
      }
    }
  }

  //PetscPrintf(PETSC_COMM_SELF, "test001\n");
/*   List test; */
/*   insertnode(&test, 11); */
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
	//	if (ibi==0) nvert[k][j][i] = 0; //reset nvert if new search

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
  
  //PetscPrintf(PETSC_COMM_SELF, "test01 %d\n",rank);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);

/*   if (user->thislevel < user->mglevels-1) { */
/*     MyNvertRestriction(user->user_f, user); */
/*   } */

  //PetscPrintf(PETSC_COMM_SELF, "test010\n");

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
	    
/* 	  if (rank==3) PetscPrintf(PETSC_COMM_SELF, "test01 %d %d %d %d %d %d %d\n",i,j,k,ys,ye,zs,ze); */
/* 	  if (i==1 && j==132 && k==119) { */
/* 	    flg=1; */
/* 	    PetscPrintf(PETSC_COMM_SELF, "test01 %d %d %d %d %d %d %d\n",i,j,k,ys,ye,zs,ze); */
/* 	  } else */
/* 	    flg=0; */

	    cutthrough = point_cell_thin(coor[k][j][i],coor[k][j][i+1],
					 coor[k][j+1][i],coor[k+1][j][i],
					 coor[k+1][j+1][i+1], ic, jc, kc, 
					 ibm, ncx, ncy, ncz, dcx, dcy, 
					 xbp_min, ybp_min, zbp_max, cell_trg, flg);

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

  //PetscPrintf(PETSC_COMM_SELF, "test010\n");
  //PetscPrintf(PETSC_COMM_SELF, "test010\n");
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

/* 	  if (nvert[k][j][ip] == 3 || nvert[k][j][im]==3 || */
/* 	      nvert[k][jp][i] == 3 || nvert[k][jm][i]==3 || */
/* 	      nvert[kp][j][i] == 3 || nvert[km][j][i]==3) { */
/* 	    nvert[k][j][i] = 1; */
	}
      }
/* 	if ((int)(nvert[k][j][i]+0.5)!=0) { */
/* 	  PetscPrintf(PETSC_COMM_SELF, "%d %d %d %le\n", i, j, k, nvert[k][j][i]); */
/* 	} */
	
    }
  }
  PetscBarrier(PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "test11\n");
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
  //  PetscPrintf(PETSC_COMM_WORLD, "test21\n");

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

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test31\n"); */

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

/* 	  if (ibi==0 && i==39 && j==41 && k==63) { */
/* 	  if (i==1 && j==100 && k==119) { */
/* 	    PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d\n", i, j, k); */
/* 	    flagprint =1; */
/* 	  } */
	  // nearestcell(coor[k][j][i], ibm, &ibm_intp);
	  nearestcell1(coor[k][j][i], ibm, &ibm_intp);
	  //	  PetscPrintf(PETSC_COMM_WORLD, "nearest cell\n");
	  InterceptionPoint(coor[k][j][i], i, j, k, &ibm_intp, user);
	  //	  PetscPrintf(PETSC_COMM_WORLD, "inteception point\n");

	  //InterceptionPoint2(coor[k][j][i], i, j, k, &ibm_intp, user);

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

  //  PetscPrintf(PETSC_COMM_WORLD, "test22\n");
  //  PetscBarrier(PETSC_NULL);

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

/*   PetscPrintf(PETSC_COMM_WORLD, "test23\n"); */

/*   PetscBarrier(PETSC_NULL); */

/*   PetscPrintf(PETSC_COMM_WORLD, "test24\n"); */

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

  //PetscPrintf(PETSC_COMM_WORLD, "Interface pts blanked block!!!! %i,\n", block_number);

/*   if (block_number>1) */
/*     Blank_Interface(user); */

/*   PetscBarrier(PETSC_NULL); */
/*   PetscPrintf(PETSC_COMM_WORLD, "test25\n"); */

  //  PetscBarrier(PETSC_NULL);
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

    //    for (ln_v=0; ln_v<ibm->n_elmt; ln_v++) {x
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
	  /*
	  if (nvert_l>1) {
	    Singularity=PETSC_TRUE;
	    break;
	  }
	  */

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
		  //   || fabs(dirdotn)<eps_tangent) {
		//if (fabs(t-dt[temp]) < epsilon) {
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
/*   PetscInt	*nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3 = ibm->nv3; */
/*   PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z; */
/*   PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp; */

  PetscInt	i, j, k, ln_v;//, n1e, n2e, n3e, nintp;
  PetscInt    cut;
  PetscInt      ks, js,is;

  ks=PetscMax(kp-1,0);
  js=PetscMax(jp-1,0);
  is=PetscMax(ip-1,0);

  i = ip; j = jp; k = kp;

  node *current;
  PetscBool	*Element_Searched;
  PetscMalloc(ibm->n_elmt*sizeof(PetscBool), &Element_Searched);


  for (ln_v=0; ln_v<ibm->n_elmt; ln_v++) {
    Element_Searched[ln_v] = PETSC_FALSE;
  }

  if (flg)
    PetscPrintf(PETSC_COMM_SELF, "ip,jp,kp,nc_xyz %d %d %d %d %d %d\n",i,j,k,ncx,ncy,ncz);
  
  for (k=ks; k<kp+2 && k<ncz; k++) {
  for (j=js; j<jp+2 && j<ncy; j++) {
  for (i=is; i<ip+2 && i<ncx; i++) {
    current = cell_trg[k*ncx*ncy+j*ncx+i].head;
    while (current) {
      ln_v = current->Node;
      //  for (ln_v=0;ln_v<ibm->n_v;ln_v++) {
      
      if (flg) PetscPrintf(PETSC_COMM_SELF, "test010 ln_v %d \n",ln_v);
      if (!Element_Searched[ln_v]) {
	Element_Searched[ln_v] = PETSC_TRUE;
	cut = ISLineTriangleIntp(p,p1,ibm,ln_v);
	if (cut) {
	  PetscFree(Element_Searched);
	  return(2);
	}
	cut = ISLineTriangleIntp(p,p2,ibm,ln_v);
	if (cut) {
	  PetscFree(Element_Searched);
	  return(2);
	}
	cut = ISLineTriangleIntp(p,p3,ibm,ln_v);
	if (cut) {
	  PetscFree(Element_Searched);
	  return(2);
	}
	cut = ISLineTriangleIntp(p,p4,ibm,ln_v);
	if (cut) {
	  PetscFree(Element_Searched);
	  return(2);
	}
      } //if
      current = current->next;
    } //while
	//  }//for
  }
  }
  }

  PetscFree(Element_Searched);  
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

    // gama = (v^2 - u^2) / (2 d \dot (v - u));
    gama = -(Dist(pu, p0)*Dist(pu, p0) - Dist(pv, p0) * Dist(pv, p0));

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

/*   struct STRC{ */
/*     PetscReal x; */
/*     PetscReal y; */
/*     PetscReal z; */
/*   } pj,pmin,po; */
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
      //    if (flagprint==1) PetscPrintf(PETSC_COMM_WORLD, "ds %e %d %e %e %e\n", tf, ln_v, nfx, p.z, p.z-z_bp[n1e]);
      if (fabs(tf) < 1.e-10) tf = 1.e-10;
      if (tf>=0) { // Point p locates on the positive side of surface triangle
    //      PetscPrintf(PETSC_COMM_WORLD, "tf%e\n", tf);
    pj.x = p.x - tf * nfx;
    pj.y = p.y - tf * nfy;
    pj.z = p.z - tf * nfz;


    //      PetscPrintf(PETSC_COMM_WORLD, "Pj%e %d\n", tf, ln_v);
    if (ISPointInTriangle(pj, p1, p2, p3, nfx, nfy, nfz) == 1) { /* The projected point
                                    is inside the
                                    triangle */
      /*    if (flagprint)
        PetscPrintf(PETSC_COMM_WORLD, "%e dmin %e %d\n", d,dmin, ln_v);*/

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
  /*  if (flagprint==1)
      PetscPrintf(PETSC_COMM_WORLD, "L %d %e %e %e %d\n", cell_min, pmin.x, pmin.y, pmin.z, n1e); */
  
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
  return 0;
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
    //tf = tf/sqrt(nfx*nfx+nfy*nfy+nfz*nfz);
    //    if (flagprint==1) PetscPrintf(PETSC_COMM_WORLD, "ds %e %d %e %e %e\n", tf, ln_v, nfx, p.z, p.z-z_bp[n1e]);

    if (fabs(tf) < 1.e-10) tf = 1.e-10;
    if (tf>=0) { // Point p locates on the positive side of surface triangle
      //      PetscPrintf(PETSC_COMM_WORLD, "tf%e\n", tf);
      pj.x = p.x - tf * nfx;
      pj.y = p.y - tf * nfy;
      pj.z = p.z - tf * nfz;


      //      PetscPrintf(PETSC_COMM_WORLD, "Pj%e %d\n", tf, ln_v);
      if (ISPointInTriangle(pj, p1, p2, p3, nfx, nfy, nfz) == 1) { /* The projected point
							is inside the
							triangle */
	/*	if (flagprint)
		PetscPrintf(PETSC_COMM_WORLD, "%e dmin %e %d\n", d,dmin, ln_v);*/

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
  return 0;
}


PetscErrorCode ICP(Cmpnts p, Cmpnts pc[9], PetscReal nfx, PetscReal nfy,
	       PetscReal nfz, IBMInfo *ibminfo, PetscInt *ip,
	       PetscInt *jp  , PetscInt *kp)
{
  PetscInt 	triangles[3][8];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;

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

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/
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
  DM		fda=user->fda;
  PetscInt	ip[9], jp[9], kp[9];
  Cmpnts	pc[9];

  PetscInt	nif;

  PetscReal	nfx, nfy, nfz;
  PetscReal 	dr;
 
  Cmpnts	***coor;

/*   DAGetCoordinates(da, &Coor); */
/*   DMDAVecGetArray(fda, Coor, &coor); */
  DMDAVecGetArray(fda, user->lCent, &coor);

  nfx = p.x - ibminfo->pmin.x;
  nfy = p.y - ibminfo->pmin.y;
  nfz = p.z - ibminfo->pmin.z;

  dr = sqrt(nfx*nfx + nfy*nfy + nfz*nfz);
  nfx /= dr; nfy /= dr; nfz /=dr;

  flagprint = 0;
  //  if (number==18) flagprint = 1;
  //  if (user->nvert[lidx(i-1,j,k,user)] == 0) {
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
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }


    /*  }

    if (user->nvert[lidx(i+1,j,k,user)] == 0) {*/
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
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }


    /*  }

    if (user->nvert[lidx(i,j-1,k,user)] == 0) {*/
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
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    /*  }

    if (user->nvert[lidx(i,j+1,k,user)] == 0) {*/
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
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    /*  }

    if (user->nvert[lidx(i,j,k-1,user)] == 0) {*/
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
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }

    /*  }

    if (user->nvert[lidx(i,j,k+1,user)] == 0) {*/
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
        DMDAVecRestoreArray(fda, user->lCent, &coor);return(0);
    }
    //  }

  DMDAVecRestoreArray(fda, user->lCent, &coor);
  //  VecDestroy(Coor);
  //  PetscPrintf(PETSC_COMM_WORLD, "End %d %d %d\n", i, j, k);
  return(0);
}

PetscErrorCode ICP2(Cmpnts p, Cmpnts pc[9], PetscReal nfx, PetscReal nfy,
	       PetscReal nfz, IBMInfo *ibminfo, PetscInt *ip,
	       PetscInt *jp  , PetscInt *kp)
{
  PetscInt 	triangles[3][8];
  Cmpnts   	p1, p2, p3;

  PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
  PetscReal	rx1, ry1, rz1, rx2, ry2, rz2;

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

	  /*	  if (flagprint==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e \n", pint.x, pint.y, pint.z, d);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", nfx, nfy, nfz);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p2.x, p2.y, p2.z);
	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", p3.x, p3.y, p3.z);
	    }*/
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
  DM		fda=user->fda;
  PetscInt	ip[9], jp[9], kp[9];
  Cmpnts	pc[9];

  PetscInt	nif;

  PetscReal	nfx, nfy, nfz;
  PetscReal 	dr;
 
  Cmpnts	***coor;

/*   DAGetCoordinates(da, &Coor); */
/*   DMDAVecGetArray(fda, Coor, &coor); */
  DMDAVecGetArray(fda, user->lCent, &coor);

  nfx = p.x - ibminfo->pmin.x;
  nfy = p.y - ibminfo->pmin.y;
  nfz = p.z - ibminfo->pmin.z;

  dr = sqrt(nfx*nfx + nfy*nfy + nfz*nfz);
  nfx /= dr; nfy /= dr; nfz /=dr;

  flagprint = 0;
  //  if (number==18) flagprint = 1;
  //  if (user->nvert[lidx(i-1,j,k,user)] == 0) {
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


    /*  }

    if (user->nvert[lidx(i+1,j,k,user)] == 0) {*/
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

    /*  }

    if (user->nvert[lidx(i,j-1,k,user)] == 0) {*/
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

    /*  }

    if (user->nvert[lidx(i,j+1,k,user)] == 0) {*/
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

    /*  }

    if (user->nvert[lidx(i,j,k-1,user)] == 0) {*/
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

    /*  }

    if (user->nvert[lidx(i,j,k+1,user)] == 0) {*/
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


    //  }

  DMDAVecRestoreArray(fda, user->lCent, &coor);
  //  VecDestroy(Coor);
  //  PetscPrintf(PETSC_COMM_WORLD, "End %d %d %d\n", i, j, k);
  return(0);
}

PetscErrorCode ibm_interpolation_advanced(UserCtx *user,
					  IBMNodes *ibm, 
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

  Cmpnts	***icsi, ***jeta, ***kzet,***coor;
 
  PetscInt i,j,k;
 
  PetscReal  cv1, cv2, cv3;
  Cmpnts   ***ucat, ***lucat;
  PetscReal	***nvert, ***p, ***lp;

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

  IBMListNode *current;
  PetscReal us_max=-1e-10, us_min=1e10;	

  PetscReal ***ustar;
  PetscInt tmp, itr_tmp=6;

  for (tmp=0; tmp<itr_tmp; tmp++) {
	
  DMDAVecGetArray(fda, user->Ucat, &ucat);  
  DMDAVecGetArray(fda, user->lUcat, &lucat);
  if (wallfunction) DMDAVecGetArray(da, user->lUstar, &ustar);
  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(da, user->lP, &lp);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->lCsi, &icsi);
  DMDAVecGetArray(fda, user->lEta, &jeta);
  DMDAVecGetArray(fda, user->lZet, &kzet);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda,user->lCent,&coor);
  current = user->ibmlist[ibi].head;
		
    while (current) {
      ibminfo = &current->ibm_intp;
      current = current->next;
				
      int ni = ibminfo->cell;
      int ip1 = ibminfo->i1, jp1 = ibminfo->j1, kp1 = ibminfo->k1;
      int ip2 = ibminfo->i2, jp2 = ibminfo->j2, kp2 = ibminfo->k2;
      int ip3 = ibminfo->i3, jp3 = ibminfo->j3, kp3 = ibminfo->k3;
      i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
				
      double sb = ibminfo->d_s, sc = sb + ibminfo->d_i;
      double sk1  = ibminfo->cr1, sk2 = ibminfo->cr2, sk3 = ibminfo->cr3;
      double cs1 = ibminfo->cs1, cs2 = ibminfo->cs2, cs3 = ibminfo->cs3;
		
      double nfx,nfy,nfz;	
      Cmpnts  Ua, Uc;  
    
      PetscReal Ua_n, Ua_nold;
     
      nfx=coor[k][j][i].x-ibminfo->pmin.x;
      nfy=coor[k][j][i].y-ibminfo->pmin.y;
      nfz=coor[k][j][i].z-ibminfo->pmin.z;
      double dr=sqrt(nfx*nfx+nfy*nfy+nfz*nfz);	
      nfx=nfx/dr;
      nfy=nfy/dr;
      nfz=nfz/dr;
 
      if (ni>=0) {
	Ua.x = ibm->uold[ibm->nv1[ni]].x * cs1 + ibm->uold[ibm->nv2[ni]].x * cs2 + ibm->uold[ibm->nv3[ni]].x * cs3;
	Ua.y = ibm->uold[ibm->nv1[ni]].y * cs1 + ibm->uold[ibm->nv2[ni]].y * cs2 + ibm->uold[ibm->nv3[ni]].y * cs3;
	Ua.z = ibm->uold[ibm->nv1[ni]].z * cs1 + ibm->uold[ibm->nv2[ni]].z * cs2 + ibm->uold[ibm->nv3[ni]].z * cs3;
	
	Ua_nold= Ua.x*nfx + Ua.y*nfy + Ua.z*nfz;

	Ua.x = ibm->u[ibm->nv1[ni]].x * cs1 + ibm->u[ibm->nv2[ni]].x * cs2 + ibm->u[ibm->nv3[ni]].x * cs3;
	Ua.y = ibm->u[ibm->nv1[ni]].y * cs1 + ibm->u[ibm->nv2[ni]].y * cs2 + ibm->u[ibm->nv3[ni]].y * cs3;
	Ua.z = ibm->u[ibm->nv1[ni]].z * cs1 + ibm->u[ibm->nv2[ni]].z * cs2 + ibm->u[ibm->nv3[ni]].z * cs3;
	
	Ua_n= Ua.x*nfx + Ua.y*nfy + Ua.z*nfz;
      }
      else {
	Ua.x = Ua.y = Ua.z = 0;
	Ua_n = 0; Ua_nold = 0.;
      }
			
      Uc.x = (lucat[kp1][jp1][ip1].x * sk1 + lucat[kp2][jp2][ip2].x * sk2 +lucat[kp3][jp3][ip3].x * sk3);
      Uc.y = (lucat[kp1][jp1][ip1].y * sk1 + lucat[kp2][jp2][ip2].y * sk2 + lucat[kp3][jp3][ip3].y * sk3);
      Uc.z = (lucat[kp1][jp1][ip1].z * sk1 + lucat[kp2][jp2][ip2].z * sk2 + lucat[kp3][jp3][ip3].z * sk3);
				
      if(wallfunction==1) {
	wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], 
		       ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
      } else if(wallfunction==2) {
	wall_function_loglaw(user, 1.e-16, sc, sb, Ua, Uc, &ucat[k][j][i], 
			     &ustar[k][j][i], ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	if (ustar[k][j][i]>us_max) us_max=ustar[k][j][i];
	if (ustar[k][j][i]<us_min) us_min=ustar[k][j][i];
      } else  {
	if(invicid || wallfunction==3) 
	  freeslip (user, sc, sb, Ua, Uc, &ucat[k][j][i], 
		    ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
	else 
	  noslip (user, sc, sb, Ua, Uc, &ucat[k][j][i], 
		  ibm->nf_x[ni], ibm->nf_y[ni], ibm->nf_z[ni]);
      }
     /*  if (i==13 && j==22 && k==8)  { */
/* 	PetscPrintf(PETSC_COMM_WORLD, "sc  %le sb %le Ua.z %le Uc.z %le nf_x %le nf_y %le nf_z %le ucat.z %le \n",sc,sb,Ua.z,Uc.z,ibm->nf_x[ni],ibm->nf_y[ni],ibm->nf_z[ni],ucat[k][j][i].z); */
/*        	PetscPrintf(PETSC_COMM_WORLD, "kp1 %d kp2 %d kp3 %d  jp1 %d  jp2 %d jp3 %d ip1 %d ip2 %d ip3 %d \n",kp1,kp2,kp3,jp1,jp2,jp3,ip1,ip2,ip3); */
/*       } */
      cv1 = lp[kp1][jp1][ip1];
      cv2 = lp[kp2][jp2][ip2];
      cv3 = lp[kp3][jp3][ip3];
	
      p[k][j][i] = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);
      if (Add_dUndt)
	p[k][j][i] += (Ua_n - Ua_nold) / user->dt * (sc-sb);
    }

  DMDAVecRestoreArray(fda, user->Ucat, &ucat);	
  DMDAVecRestoreArray(fda, user->lUcat, &lucat);
  DMDAVecRestoreArray(fda,user->lCent,&coor);	
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);


  DMDAVecGetArray(fda, user->lUcat, &lucat);
  
  
  Cmpnts    uc;

  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	double f = 1.0;
/* 	//	if(immersed==3) f = 0; */
/* 	if (nvert[k][j][i]+nvert[k][j][i+1]> 2.1 && nvert[k][j][i]+nvert[k][j][i+1]< innerblank) */
/* 	  ucont[k][j][i].x=0.; */
/* 	if (nvert[k][j][i]+nvert[k][j+1][i]> 2.1 && nvert[k][j][i]+nvert[k][j+1][i]< innerblank) */
/* 	  ucont[k][j][i].y=0.; */
/* 	if (nvert[k][j][i]+nvert[k+1][j][i]> 2.1 && nvert[k][j][i]+nvert[k+1][j][i]< innerblank) */
/* 	  ucont[k][j][i].z=0.; */
	  
	if ((int)(nvert[k][j][i]+0.5) ==1) {// || (int)(nvert_o[k][j][i]+0.5)==1) {
	  uc.x = (lucat[k][j][i].x + lucat[k][j][i+1].x) * 0.5;
	  uc.y = (lucat[k][j][i].y + lucat[k][j][i+1].y) * 0.5;
	  uc.z = (lucat[k][j][i].z + lucat[k][j][i+1].z) * 0.5;
	  if (i<mx-2)
	  ucont[k][j][i].x = (uc.x * icsi[k][j][i].x + uc.y * icsi[k][j][i].y + uc.z * icsi[k][j][i].z) * f;
	    
	  uc.x = (lucat[k][j+1][i].x + lucat[k][j][i].x) * 0.5;
	  uc.y = (lucat[k][j+1][i].y + lucat[k][j][i].y) * 0.5;
	  uc.z = (lucat[k][j+1][i].z + lucat[k][j][i].z) * 0.5;
	  if (j<my-2)
	  ucont[k][j][i].y = (uc.x * jeta[k][j][i].x + uc.y * jeta[k][j][i].y + uc.z * jeta[k][j][i].z) * f;
		  
	  uc.x = (lucat[k+1][j][i].x + lucat[k][j][i].x) * 0.5;
	  uc.y = (lucat[k+1][j][i].y + lucat[k][j][i].y) * 0.5;
	  uc.z = (lucat[k+1][j][i].z + lucat[k][j][i].z) * 0.5;
	  if (k<mz-2)		
	  ucont[k][j][i].z = (uc.x * kzet[k][j][i].x + uc.y * kzet[k][j][i].y + uc.z * kzet[k][j][i].z) * f;
	}

	if ((int)(nvert[k][j][i+1]+0.5)==1) {
	  uc.x = (lucat[k][j][i].x + lucat[k][j][i+1].x) * 0.5;
	  uc.y = (lucat[k][j][i].y + lucat[k][j][i+1].y) * 0.5;
	  uc.z = (lucat[k][j][i].z + lucat[k][j][i+1].z) * 0.5;
	  if (i<mx-2)
	  ucont[k][j][i].x = (uc.x * icsi[k][j][i].x + uc.y * icsi[k][j][i].y + uc.z * icsi[k][j][i].z) * f;
	}
		
	if ((int)(nvert[k][j+1][i]+0.5)==1) {
	  uc.x = (lucat[k][j+1][i].x + lucat[k][j][i].x) * 0.5;
	  uc.y = (lucat[k][j+1][i].y + lucat[k][j][i].y) * 0.5;
	  uc.z = (lucat[k][j+1][i].z + lucat[k][j][i].z) * 0.5;
	  if (j<my-2)		
	  ucont[k][j][i].y = (uc.x * jeta[k][j][i].x + uc.y * jeta[k][j][i].y + uc.z * jeta[k][j][i].z) * f;
	}

	if ((int)(nvert[k+1][j][i]+0.5)==1) {
	  uc.x = (lucat[k+1][j][i].x + lucat[k][j][i].x) * 0.5;
	  uc.y = (lucat[k+1][j][i].y + lucat[k][j][i].y) * 0.5;
	  uc.z = (lucat[k+1][j][i].z + lucat[k][j][i].z) * 0.5;
	  if (k<mz-2)		
	  ucont[k][j][i].z = (uc.x * kzet[k][j][i].x + uc.y * kzet[k][j][i].y + uc.z * kzet[k][j][i].z )* f;
	}
      }
	

  DMDAVecRestoreArray(fda, user->lUcat, &lucat);
  if (wallfunction) {
    DMDAVecRestoreArray(da, user->lUstar, &ustar);
    PetscReal us_maxSum, us_minSum;
    MPI_Allreduce(&us_max, &us_maxSum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    MPI_Allreduce(&us_min, &us_minSum,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);
   /*  PetscGlobalMax(&us_max, &us_maxSum,PETSC_COMM_WORLD); */
/*     PetscGlobalMin(&us_min, &us_minSum,PETSC_COMM_WORLD); */
    PetscPrintf(PETSC_COMM_WORLD, "!!! Ustar Max Min %le %le\n", us_maxSum, us_minSum);
  }
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(da, user->lP, &lp);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->lCsi, &icsi);
  DMDAVecRestoreArray(fda, user->lEta, &jeta);
  DMDAVecRestoreArray(fda, user->lZet, &kzet);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

/*   PetscInt flg; */
/*   PetscReal ibm_Flux, ibm_Area; */
/*   if (sediment)  */
/*     flg=2; */
/*   else */
/*     flg=0; */
/*   VolumeFlux(user, &ibm_Flux, &ibm_Area, flg); */

  Contra2Cart(user);

  }

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  // PetscPrintf(PETSC_COMM_WORLD,"linear interpolation works \n");	 
  return(0);
}

PetscErrorCode ibm_interpolation_advanced_old(UserCtx *user, 
					      IBMNodes *ibm, 
					      //  FSInfo *fsi,
					      PetscInt ibi,
					      PetscInt Add_dUndt)
{
  DM		da = user->da, fda = user->fda;
  //  Vec		Ucont = user->Ucont;
  Cmpnts	***ucont;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***icsi, ***jeta, ***kzet;
  
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

/*   DMDAVecGetArray(da, user->lP, &p); */
/*   PetscReal Pconst; */
/*   if (xs==0 && ys==0 && zs==0) { */
/*     Pconst=p[2][2][2]; */
/*     MPI_Bcast(&Pconst, 1, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*   } else { */
/*     MPI_Bcast(&Pconst, 1, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*   }     */
/*   DMDAVecRestoreArray(da, user->lP, &p); */
  
/*   VecShift(user->lP, -Pconst); */

  for (tmp=0; tmp<9; tmp++) {

    DMDAVecGetArray(fda, user->lUcat, &ucat);
    DMDAVecGetArray(da, user->lP, &p);

    DMDAVecGetArray(fda, user->Ucont, &ucont);
    DMDAVecGetArray(fda, user->lICsi, &icsi);
    DMDAVecGetArray(fda, user->lJEta, &jeta);
    DMDAVecGetArray(fda, user->lKZet, &kzet);
    DMDAVecGetArray(da, user->lNvert, &nvert);
    DMDAVecGetArray(da, user->lNvert_o, &nvert_o);

    current = user->ibmlist[ibi].head;
    while (current) {
      ibminfo = &current->ibm_intp;
      current = current->next;
      i = ibminfo->ni; j= ibminfo->nj;    k = ibminfo->nk;
      sb = ibminfo->d_s; sc = sb + ibminfo->d_i;
/*       lhs[0][0] = 1.; */
/*       lhs[0][1] = - sb*sb; */
/*       lhs[0][2] = -sb; */
/*       lhs[1][0] = 0.; */
/*       lhs[1][1] =  sc*sc;// should be -sc*sc */
/*       lhs[1][2] =  -sc; */
/*       lhs[2][0] = ((sc-sb)*(sc-sb) - sb*sb)/(sc-sb); */
/*       lhs[2][1] = -2 * sc * sb*sb; */
/*       lhs[2][2] = -sb*sc;*/

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

/*       if ((int)(nvert[kp1][jp1][ip1]+0.5)>1) { */
/* 	ucat[kp1][jp1][ip1].z=0.; */
/* 	ucat[kp1][jp1][ip1].y=0.; */
/* 	ucat[kp1][jp1][ip1].x=0.; */
/*       } */
/*       if ((int)(nvert[kp2][jp2][ip2]+0.5)>1) { */
/* 	ucat[kp2][jp2][ip2].z=0.; */
/* 	ucat[kp2][jp2][ip2].y=0.; */
/* 	ucat[kp2][jp2][ip2].x=0.; */
/*       } */
/*       if ((int)(nvert[kp3][jp3][ip3]+0.5)>1) { */
/* 	ucat[kp3][jp3][ip3].z=0.; */
/* 	ucat[kp3][jp3][ip3].y=0.; */
/* 	ucat[kp3][jp3][ip3].x=0.; */
/*       } */

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
	/*       rhs_l[2][0] = (sc-sb) * phia - (sb*sb/(sc-sb)) * phic; */
	rhs_l[2][0] = -(sc-sb)/sc/sb*phia + (1.-(sc-sb)/sc)/(sc-sb)*phic;
	ucat[k][j][i].x = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc;//detmnt(rhs_l)/dtm; 
      }
/*       if (i==1 && j==75 && k==45) */
/*       PetscPrintf(PETSC_COMM_WORLD, "u_x %e %e %e %e %e %d nbn\n", sb, sc, sk3, ucat[k][j][i].x, phic ,ni); */

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
      /*       cv1 = 0; */
      /*       cv2 = 0; */
      /*       cv3=0; */

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
	/*       rhs_l[2][0] = (sc-sb) * phia - (sb*sb/(sc-sb)) * phic; */
	rhs_l[2][0] = -(sc-sb)/sc/sb*phia + (1.-(sc-sb)/sc)/(sc-sb)*phic;
	ucat[k][j][i].y = detmnt(rhs_l)/dtm;//phia + (phic-phia) * sb / sc; //detmnt(rhs_l)/dtm;
      }
/*       if (i==1 && j==75 && k==45) */
/*       PetscPrintf(PETSC_COMM_WORLD, "u_y %e %e %e %e %e %e %d nbn\n", sb, sc, sk3, ucat[k][j][i].y, phic ,phia + (phic-phia) * sb / sc,ni); */

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
      /*       cv1 = 0; */
      /*       cv2 = 0; */
      /*       cv3=0; */

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

/*       if (i==42 && j==199 && k==81) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "PHIA %e %e %e %e\n", phia, cv1, cv2, cv3); */
/*       } */
      cv1 = ucat[kp1][jp1][ip1].z;
      cv2 = ucat[kp2][jp2][ip2].z;
      cv3 = ucat[kp3][jp3][ip3].z;
      phic = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);

      if (invicid) {
	ucat[k][j][i].z = phic;		
      } else {
	rhs_l[0][0] = phia;
	rhs_l[1][0] = phic-phia;
	/*       rhs_l[2][0] = (sc-sb) * phia - (sb*sb/(sc-sb)) * phic; */
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

/*       if (i==1 && j==75 && k==45) */
/*       PetscPrintf(PETSC_COMM_WORLD, "u_z %e %e %e %e %e %d nbn\n", sb, sc, sk3, ucat[k][j][i].z, phic ,ni); */
/*       PetscPrintf(PETSC_COMM_WORLD, "u_z %e %e %e %e %e %e %d nbn\n", sb, sc, sk3, ucat[k][j][i].z, phic ,phia + (phic-phia) * sb / sc,ni); */

/*       if (ucat[k][j][i].x > 1 || ucat[k][j][i].y > 1 || ucat[k][j][i].z > 1) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "Larte %i %i %i %e %e %e %e %e %i\n", i, j, k,ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z, sb, sc, tmp); */
/*       } */
/*       if (ucat[k][j][i].x != ucat[k][j][i].x || */
/* 	  ucat[k][j][i].y != ucat[k][j][i].y || */
/* 	  ucat[k][j][i].z != ucat[k][j][i].z) { */
/*       if (k==72 && (j==114 || j==115) && xs ==135) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "SS %i %i %i %i %i %i\n", xs, xe, ys, ye, zs, ze); */
/* 	PetscPrintf(PETSC_COMM_SELF, "INTP %i %i %i %i %e %e %e %e %e %e\n", i, j, k, tmp, sk1, sk2, sk3, cs1, cs2, cs3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "KP1 %i %i %i\n", ip1, jp1, kp1); */
/* 	PetscPrintf(PETSC_COMM_SELF, "KP2 %i %i %i\n", ip2, jp2, kp2); */
/* 	PetscPrintf(PETSC_COMM_SELF, "KP3 %i %i %i\n", ip3, jp3, kp3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "FF %e %e %e %e %e %e\n", cv1, cv2, cv3, phia, phic, cv1*cs1+cv2*cs2+cv3*cs3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "FF %e %e %e\n", ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z); */
/*       } */

      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];
      p[k][j][i] = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);

/*       if (i==1 && j==100 && k==119) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "INTP IJK %d %d %d P %le %le %le %le\n", i, j, k, p[k][j][i],cv1,cv2,cv3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "INTP host1 %d %d %d\n", ip1,jp1,kp1); */
/* 	PetscPrintf(PETSC_COMM_SELF, "INTP host1 %d %d %d\n", ip2,jp2,kp2); */
/* 	PetscPrintf(PETSC_COMM_SELF, "INTP host1 %d %d %d\n", ip3,jp3,kp3); */
/* 	//	flagprint =1; */
/*       } */
      /*
      cv1 = 3.*fsi->S_new[1]-4.*fsi->S_real[1]+fsi->S_realm1[1];
      cv2 = 3.*fsi->S_new[3]-4.*fsi->S_real[3]+fsi->S_realm1[3];
      cv3 = 3.*fsi->S_new[5]-4.*fsi->S_real[5]+fsi->S_realm1[5];
      */
      // p(n+1)-p(n) = dn*-St(n.dU/dt) ==>
      // p(n) = p(n+1) + dn*St(n.dU/dt)
      //p[k][j][i] += sb*user->st*(cv1*nfx+cv2*nfy+cv3*nfz)*0.5/user->dt;

      // if (moveframe==0)
      p[k][j][i] += (un - unold) / user->dt * (sc-sb);
/*       if (Add_dUndt) */
/* 	p[k][j][i] += user->st*(3.*un - 4.*unold + unrm1)*0.5 / user->dt * (sc-sb); */

/*       if (i==1 && j==100 && k==119) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "INP IJK %d %d %d P %le Un %le %le %le\n", i, j, k, p[k][j][i],un,unold,unrm1); */
/* 	//	flagprint =1; */
/*       } */

    }

    DMDAVecRestoreArray(fda, user->lUcat, &ucat);
    DMDAVecRestoreArray(da, user->lP, &p);

    DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
    DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

/*     DALocalToLocalBegin(fda, user->lUcat, INSERT_VALUES, user->lUcat); */
/*     DALocalToLocalEnd(fda, user->lUcat, INSERT_VALUES, user->lUcat); */

    DMLocalToLocalBegin(da, user->lP, INSERT_VALUES, user->lP);
    DMLocalToLocalEnd(da, user->lP, INSERT_VALUES, user->lP);

/*     DALocalToGlobal(da, user->lP, INSERT_VALUES, user->P); */

/*     DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP); */
/*     DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP); */

    PetscBarrier(PETSC_NULL);
    DMDAVecGetArray(fda, user->lUcat, &ucat);
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if ((int)(nvert[k][j][i]+0.5) ==1) {// || (int)(nvert_o[k][j][i]+0.5)==1) {
	    ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
	    ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
	    ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
	    if (i<mx-2)
	    ucont[k][j][i].x = ucx * icsi[k][j][i].x +
	      ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z;
	  
	    ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
	    if (j<my-2)
	    ucont[k][j][i].y = ucx * jeta[k][j][i].x +
	      ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z;
	   
	    ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
	    if (k<mz-2)
	    ucont[k][j][i].z = ucx * kzet[k][j][i].x +
	      ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z;
	  }

	  if ((int)(nvert[k][j][i+1]+0.5)==1) {
	    ucx = (ucat[k][j][i].x + ucat[k][j][i+1].x) * 0.5;
	    ucy = (ucat[k][j][i].y + ucat[k][j][i+1].y) * 0.5;
	    ucz = (ucat[k][j][i].z + ucat[k][j][i+1].z) * 0.5;
	    if (i<mx-2)
	    ucont[k][j][i].x = ucx * icsi[k][j][i].x +
	      ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z;
	  }
	
	  if ((int)(nvert[k][j+1][i]+0.5)==1) {
	    ucx = (ucat[k][j+1][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k][j+1][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k][j+1][i].z + ucat[k][j][i].z) * 0.5;
	    if (j<my-2)
	    ucont[k][j][i].y = ucx * jeta[k][j][i].x +
	      ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z;
	  }

	  if ((int)(nvert[k+1][j][i]+0.5)==1) {
	    ucx = (ucat[k+1][j][i].x + ucat[k][j][i].x) * 0.5;
	    ucy = (ucat[k+1][j][i].y + ucat[k][j][i].y) * 0.5;
	    ucz = (ucat[k+1][j][i].z + ucat[k][j][i].z) * 0.5;
	    if (k<mz-2)
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

    /*   DALocalToGlobal(fda, user->lUcont, INSERT_VALUES, user->Ucont); */
  /*   DALocalToGlobal(fda, user->lUcat, INSERT_VALUES, user->Ucat); */


    PetscBarrier(PETSC_NULL);

/*     DALocalToGlobal(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
    
    DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
/*     DALocalToLocalBegin(fda, user->lUcont, INSERT_VALUES, user->lUcont); */
/*     DALocalToLocalEnd(fda, user->lUcont, INSERT_VALUES, user->lUcont); */

    Contra2Cart(user);
  }

/*   DALocalToGlobal(fda, user->lUcont, INSERT_VALUES, user->Ucont); */

/*   DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*   DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  return 0;
}

PetscErrorCode ibm_interpolation_advanced2(UserCtx *user, IBMNodes *ibm)
{
  DM		da = user->da, fda = user->fda;
  //  Vec		Ucont = user->Ucont;
  Cmpnts	***ucont;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***icsi, ***jeta, ***kzet;

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
  PetscReal un;
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

      //PetscPrintf(PETSC_COMM_WORLD, "%d %d mode\n", ibminfo->imode, ibminfo->iimode);

      
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
      /*       cv1 = 0; */
      /*       cv2 = 0; */
      /*       cv3=0; */

      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);
      un = phia * nfx;
/*       unold = (ibm->uold[ibm->nv1[ni]].x * cs1 + */
/* 	       ibm->uold[ibm->nv2[ni]].x * cs2 + */
/* 	       ibm->uold[ibm->nv3[ni]].x * cs3) * nfx; */
      
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


      //    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %e %e %d nbn\n", sb, sc, sk3, ucat[k][j][i].x, dtm ,nbn);
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
      /*       cv1 = 0; */
      /*       cv2 = 0; */
      /*       cv3=0; */

      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);

      un += phia * nfy;
/*       unold += (ibm->uold[ibm->nv1[ni]].x * cs1 + */
/* 	       ibm->uold[ibm->nv2[ni]].x * cs2 + */
/* 	       ibm->uold[ibm->nv3[ni]].x * cs3) * nfy; */

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
      /*       cv1 = 0; */
      /*       cv2 = 0; */
      /*       cv3=0; */

      phia = (cv1 * cs1 + cv2 * cs2 + cv3 * cs3);

      un += phia * nfz;
/*       unold += (ibm->uold[ibm->nv1[ni]].x * cs1 + */
/* 	       ibm->uold[ibm->nv2[ni]].x * cs2 + */
/* 	       ibm->uold[ibm->nv3[ni]].x * cs3) * nfz; */

/*       if (i==42 && j==199 && k==81) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "PHIA %e %e %e %e\n", phia, cv1, cv2, cv3); */
/*       } */
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

/*       if (ucat[k][j][i].x > 1 || ucat[k][j][i].y > 1 || ucat[k][j][i].z > 1) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "Larte %i %i %i %e %e %e %e %e %i\n", i, j, k,ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z, sb, sc, tmp); */
/*       } */
/*       if (ucat[k][j][i].x != ucat[k][j][i].x || */
/* 	  ucat[k][j][i].y != ucat[k][j][i].y || */
/* 	  ucat[k][j][i].z != ucat[k][j][i].z) { */
/*       if (k==72 && (j==114 || j==115) && xs ==135) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "SS %i %i %i %i %i %i\n", xs, xe, ys, ye, zs, ze); */
/* 	PetscPrintf(PETSC_COMM_SELF, "INTP %i %i %i %i %e %e %e %e %e %e\n", i, j, k, tmp, sk1, sk2, sk3, cs1, cs2, cs3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "KP1 %i %i %i\n", ip1, jp1, kp1); */
/* 	PetscPrintf(PETSC_COMM_SELF, "KP2 %i %i %i\n", ip2, jp2, kp2); */
/* 	PetscPrintf(PETSC_COMM_SELF, "KP3 %i %i %i\n", ip3, jp3, kp3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "FF %e %e %e %e %e %e\n", cv1, cv2, cv3, phia, phic, cv1*cs1+cv2*cs2+cv3*cs3); */
/* 	PetscPrintf(PETSC_COMM_SELF, "FF %e %e %e\n", ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z); */
/*       } */

      cv1 = p[kp1][jp1][ip1];
      cv2 = p[kp2][jp2][ip2];
      cv3 = p[kp3][jp3][ip3];

      cv11 = p[kp11][jp11][ip11];
      cv22 = p[kp22][jp22][ip22];
      cv33 = p[kp33][jp33][ip33];

      p[k][j][i] = 0.5*(cv1 * sk1 + cv2 * sk2 + cv3 * sk3 + 
			cv11*sk11 + cv22*sk22 + cv33*sk33) ;

/*       p[k][j][i] += (un - unold) / user->dt * sb; */
    }

    DMDAVecRestoreArray(fda, user->lUcat, &ucat);
    DMDAVecRestoreArray(da, user->lP, &p);

    DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
    DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

/*     DALocalToLocalBegin(fda, user->lUcat, INSERT_VALUES, user->lUcat); */
/*     DALocalToLocalEnd(fda, user->lUcat, INSERT_VALUES, user->lUcat); */

    DMLocalToLocalBegin(da, user->lP, INSERT_VALUES, user->lP);
    DMLocalToLocalEnd(da, user->lP, INSERT_VALUES, user->lP);

/*     DALocalToGlobal(da, user->lP, INSERT_VALUES, user->P); */

/*     DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP); */
/*     DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP); */


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

    /*   DALocalToGlobal(fda, user->lUcont, INSERT_VALUES, user->Ucont); */
  /*   DALocalToGlobal(fda, user->lUcat, INSERT_VALUES, user->Ucat); */


    PetscBarrier(PETSC_NULL);

/*     DALocalToGlobal(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
    
    DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
    DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
/*     DALocalToLocalBegin(fda, user->lUcont, INSERT_VALUES, user->lUcont); */
/*     DALocalToLocalEnd(fda, user->lUcont, INSERT_VALUES, user->lUcont); */

    Contra2Cart(user);
  }

/*   DALocalToGlobal(fda, user->lUcont, INSERT_VALUES, user->Ucont); */

/*   DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*   DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  return 0;
}

PetscErrorCode SetSolidNodeVelocityToZero(UserCtx *user) 
{
  DM		da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal	***nvert;
  Cmpnts	***ucat, ***ucont;
  PetscReal     solid=2.99, blanking=7.0;
  PetscInt      i,j,k;

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] > solid && nvert[k][j][i]<blanking)  {
	  ucat[k][j][i].x=0.;
	  ucat[k][j][i].y=0.;  
	  ucat[k][j][i].z=0.;
	}
      }
    }
  }

  for (k=zs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	if (nvert[k][j][i] > solid && nvert[k][j][i]<blanking)  {
	  ucat[k][j][i].x=0.;
	  ucat[k][j][i].y=0.;  
	  ucat[k][j][i].z=0.;

	  if (k>0 && j>0) ucont[k][j][i].x=0.;
	  if (k>0 && i>0) ucont[k][j][i].y=0.;  
	  if (j>0 && i>0) ucont[k][j][i].z=0.;
	}

	if (i>0) 
	  if (nvert[k][j][i]<solid && (nvert[k][j][i-1]>solid && 
				       nvert[k][j][i-1]<blanking)) 
	    ucont[k][j][i].x=0.;

	if (j>0) 
	  if (nvert[k][j][i]<solid && (nvert[k][j-1][i]>solid && 
				       nvert[k][j-1][i]<blanking)) 
	    ucont[k][j][i].y=0.;

	if (k>0) 
	  if (nvert[k][j][i]<solid && (nvert[k-1][j][i]>solid && 
				       nvert[k-1][j][i]<blanking)) 
	    ucont[k][j][i].z=0.;

	if (i==mx-2 && 
	    nvert[k][j][i-1]>solid && nvert[k][j][i-1]<blanking) 
	  ucont[k][j][i].x=0;

	if (j==my-2 && 
	    nvert[k][j-1][i]>solid && nvert[k][j-1][i]<blanking) 
	  ucont[k][j][i].y=0;

	if (k==mz-2 && 
	    nvert[k-1][j][i]>solid && nvert[k-1][j][i]<blanking) 
	  ucont[k][j][i].z=0;

      }
    }
  }
 
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  return(0);
}
