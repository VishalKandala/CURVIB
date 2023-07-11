static char help[] = "Testing programming!";

#include "variables.h"
#include "petscdmda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

/* #ifdef TECIO */
/* #include "TECIO.h" */
/* #endif */

PetscReal  FluxInSum, FluxOutSum;
PetscInt   ti, block_number, Flux_in;
PetscReal  angle;
PetscInt   invicid = 0;
PetscInt   thin=0;
PetscInt   wing=0, cop=0, rotateframe=0;
PetscReal  CMx_c=0., CMy_c=0., CMz_c=0.;
PetscInt   NumberOfBodies=1;
PetscReal  L_dim=1.;
PetscInt   dgf_z=0,dgf_y=1,dgf_x=0;
PetscInt   dgf_az=0,dgf_ay=0,dgf_ax=1 ;
PetscReal  St_exp=0.35,wavelength=0.65;
PetscInt   tiout = 10;
PetscInt   STRONG_COUPLING=0, wallfunction=0;
PetscInt   fish=0, eel=0,rheology=0,pizza=0,sediment=0;
PetscInt   radi=10,tistart,Ddata=0,qv=0,unorm=0,vx=0,nudata=1,averaging,max_angle;
PetscInt   visflg=0;
char       orient[] = "xx00";
/* typedef struct { */
/*   PetscReal t, f; */
/* } FlowWave; */

typedef struct {
  PassiveScalar u, v, w, p;
} PassiveField;

typedef struct {
  PetscScalar u, v, w;
} Field;

Cmpnts rotateT_y(Cmpnts v0, Cmpnts v_c, Cmpnts rot);
Cmpnts rotate_xyz(Cmpnts v, Cmpnts v_c, Cmpnts rot);

Cmpnts MINUS(Cmpnts v1,Cmpnts v2)
{
  Cmpnts v4;
  v4.x  = v1.x - v2.x ;
  v4.y  = v1.y - v2.y ;
  v4.z  = v1.z - v2.z ;
 
  return(v4);
} 

Cmpnts Cross(Cmpnts v1,Cmpnts v2)
{
  // output = v1 x v2
  Cmpnts v1cv2;
  v1cv2.x  = v1.y * v2.z - v2.y * v1.z;
  v1cv2.y  =-v1.x * v2.z + v2.x * v1.z;
  v1cv2.z  = v1.x * v2.y - v2.x * v1.y;
 
  return(v1cv2);
}

//----------------------------------------------
PetscErrorCode ReadCoordinates(UserCtx *user) {

  Cmpnts ***coor;

  Vec Coor;
  PetscInt bi, i, j, k, rank, IM, JM, KM;
  PetscReal *gc;
  FILE *fd;
  PetscReal	d0 = 1.;
  PetscInt    generate_grid=0, grid1d=0, nblk=block_number;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscOptionsGetInt(PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);

  PetscReal	cl = 1.;
  PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);

  PetscReal L_x,L_y,L_z;

  if (generate_grid) {
    PetscOptionsGetReal(PETSC_NULL, "-L_x", &L_x, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-L_y", &L_y, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-L_z", &L_z, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD,"ksi eta zeta length  are %le %le %le \n",L_x,L_y,L_z);
    block_number=1;
  } else {
    if (!rank) {
      fd = fopen("grid.dat", "r");
      fscanf(fd, "%i\n", &block_number);
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }
  PetscInt imm[block_number], kmm[block_number], jmm[block_number];
  if (generate_grid) {
    PetscOptionsGetIntArray(PETSC_NULL, "-im", imm, &nblk, PETSC_NULL);
    PetscOptionsGetIntArray(PETSC_NULL, "-jm", jmm, &nblk, PETSC_NULL);
    PetscOptionsGetIntArray(PETSC_NULL, "-km", kmm, &nblk, PETSC_NULL);
  }
  for (bi=0; bi<block_number; bi++) {

    //if (bi==0)  cl = 1.;
    //if (bi==1)  cl = 0.85;

    if (!rank) {
      if (!generate_grid)
	fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
      
      else {
	user[bi].IM=imm[bi];
	user[bi].JM=jmm[bi];
	user[bi].KM=kmm[bi];
      }	
      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
      
      MPI_Bcast(&(user[bi].IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&(user[bi].IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
    }
    
    DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
	       user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
		 PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,&(user[bi].da));
    
    DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DMGetCoordinateDM(user[bi].da, &(user[bi].fda));
    DMDAGetLocalInfo(user[bi].da, &(user[bi].info));
    
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  
    PetscOptionsGetInt(PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);
    if (grid1d) PetscMalloc((IM+JM+KM)*sizeof(PetscReal), &gc);
    else        PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);
    DMGetCoordinatesLocal(user[bi].da, &Coor);
    DMDAVecGetArray(user[bi].fda, Coor, &coor);

    if (!rank) {
      if (grid1d) {
	PetscReal xx;
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
		coor[k][j][i].x = *(gc + i)/cl*L_dim;
		coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
		coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
	      }
	    }
	  }
	}
	
      } else { // if 3d gridgen file
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3) = L_x/(IM-1.) * i;
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3+1) = L_y/(JM-1.) * j;
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
	    else
	      *(gc+(k*JM*IM + j*IM + i)*3+2) = L_z/(KM-1.) * k;
	    }
	  }
	}	
	
	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
      
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
      }  
    }
    else {
      if (grid1d) {
	MPI_Bcast(gc, (IM+JM+KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + i)/cl*L_dim;
		coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
		coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
	      }
	    }
	  }
	}

      } else { // if 3d gridgen file

	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
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
      }
    }
    PetscFree(gc);
    DMDAVecRestoreArray(user[bi].fda, Coor, &coor);
    
    Vec	gCoor;
    DMGetCoordinates(user[bi].da, &gCoor);
    DMLocalToGlobalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    DMLocalToGlobalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    
    DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);

  }
  if (!rank) {
    if(!generate_grid)
      fclose(fd);
  }
  
  for (bi=0; bi<block_number; bi++) {
    user[bi]._this = bi;
  }

  return(0);
}

//-------------------------------------------------------
PetscErrorCode SetCoordinates3d(UserCtx *user)
{
  /* Set 3d grid coordinates. */
  DM	cda = user->fda, da = user->da;
  Vec		Coor;
  Cmpnts	***coor;
  PetscInt	IM, JM, KM;
  PetscInt	rank;
  
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
  PetscReal	*gc;
  PetscReal	d0 = 1.;

  DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(cda, Coor, &coor);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) {
    FILE *fd;
    fd = fopen("grid.dat", "r");

    fscanf(fd, "%i\n", &block_number);
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    fclose(fd);
  }
  else {
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  }
  PetscMalloc(block_number*sizeof(UserCtx), &user);

  ReadCoordinates(user);

  Vec	gCoor;
  DMGetCoordinates(da, &gCoor);
  DMLocalToGlobalBegin(cda, Coor, INSERT_VALUES, gCoor);
  DMLocalToGlobalEnd(cda, Coor, INSERT_VALUES, gCoor);

  DMGlobalToLocalBegin(cda, gCoor, INSERT_VALUES, Coor);
  DMGlobalToLocalEnd(cda, gCoor, INSERT_VALUES, Coor);

  return(0);
}

PetscErrorCode Contra2Cart(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	mat[3][3], det, det0, det1, det2;

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  Vec		Aj = user->lAj;
  Vec		Ucont = user->lUcont, Ucat = user->Ucat;
  Vec		Ubcs = user->Bcs.Ubcs;
  Cmpnts	***csi, ***eta, ***zet;
  PetscReal	***aj;
  Cmpnts	***ucont, ***ucat, ***ubcs;

  PetscReal	***nvert;
  PetscReal	coef = 0.125;

  PetscReal	q[3]; //local working array
  PetscInt	i, j, k;
  PetscReal	c1, c2, c3, c4, g1, g2, g3, g4, ucon;
  Cmpnts	***gs;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(fda, Csi, &csi);
  DMDAVecGetArray(fda, Eta, &eta);
  DMDAVecGetArray(fda, Zet, &zet);
  DMDAVecGetArray(da,  Aj,  &aj);

  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(fda, Ubcs,  &ubcs);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5)==0) {
	  mat[0][0] = 0.5 * (csi[k][j][i-1].x + csi[k][j][i].x);
	  mat[0][1] = 0.5 * (csi[k][j][i-1].y + csi[k][j][i].y);
	  mat[0][2] = 0.5 * (csi[k][j][i-1].z + csi[k][j][i].z);
	         	      
	  mat[1][0] = 0.5 * (eta[k][j-1][i].x + eta[k][j][i].x);
	  mat[1][1] = 0.5 * (eta[k][j-1][i].y + eta[k][j][i].y);
	  mat[1][2] = 0.5 * (eta[k][j-1][i].z + eta[k][j][i].z);
	         	      
	  mat[2][0] = 0.5 * (zet[k-1][j][i].x + zet[k][j][i].x);
	  mat[2][1] = 0.5 * (zet[k-1][j][i].y + zet[k][j][i].y);
	  mat[2][2] = 0.5 * (zet[k-1][j][i].z + zet[k][j][i].z);

	  ucon = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	  if (ucon  > 0) {
	    if (i>1 && (int)(nvert[k][j][i-2]+0.5)==0) {

	      q[0] = coef * (-    ucont[k][j][i-2].x -
			     2. * ucont[k][j][i-1].x +
			     3. * ucont[k][j][i  ].x) + ucont[k][j][i-1].x;
	    }
	    else {
	      q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	    }
	  }
	  else {
	    if (i < mx-2 &&(int)(nvert[k][j][i+1]+0.5)==0) {

	      q[0] = coef * (-    ucont[k][j][i+1].x -
			     2. * ucont[k][j][i  ].x +
			     3. * ucont[k][j][i-1].x) + ucont[k][j][i].x;
	    }
	    else {
	      q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	    }
	  }


	  ucon = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	  if (ucon > 0) {
	    if (j>1&&(int)(nvert[k][j-2][i]+0.5)==0) {

	      q[1] = coef * (-    ucont[k][j-2][i].y -
			     2. * ucont[k][j-1][i].y +
			     3. * ucont[k][j  ][i].y) + ucont[k][j-1][i].y;
	    }
	    else {
	      q[1] = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	    }
	  }
	  else {
	    if (j < my-2 &&(int)(nvert[k][j+1][i]+0.5)==0) {

	      q[1] = coef * (-    ucont[k][j+1][i].y -
			     2. * ucont[k][j  ][i].y +
			     3. * ucont[k][j-1][i].y) + ucont[k][j][i].y;
	    }
	    else {
	      q[1] = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	    }
	  }

	  ucon = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	  if (ucon > 0) {
	    if (k>1&&(int)(nvert[k-2][j][i]+0.5)==0) {

	      q[2] = coef * (-    ucont[k-2][j][i].z -
			     2. * ucont[k-1][j][i].z +
			     3. * ucont[k  ][j][i].z) + ucont[k-1][j][i].z;
	    }
	    else {
	      q[2] = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	    }
	  }
	  else {
	    if (k < mz-2 &&(int)(nvert[k+1][j][i]+0.5)==0) {

	      q[2] = coef * (-    ucont[k+1][j][i].z -
			     2. * ucont[k  ][j][i].z +
			     3. * ucont[k-1][j][i].z) + ucont[k][j][i].z;
	    }
	    else {
	      q[2] = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	    }
	  }


	  det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

	  det0 = q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
	    q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

	  det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
	    q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

	  det2 = q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
	    q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
	    q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

	  ucat[k][j][i].x = det0 / det;
	  ucat[k][j][i].y = det1 / det;
	  ucat[k][j][i].z = det2 / det;

	}
      }
    }
  }

  DMDAVecRestoreArray(fda, Csi, &csi);
  DMDAVecRestoreArray(fda, Eta, &eta);
  DMDAVecRestoreArray(fda, Zet, &zet);
  DMDAVecRestoreArray(da,  Aj,  &aj);

  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMDAVecRestoreArray(fda, Ubcs,  &ubcs);

  VecAssemblyBegin(user->Ucat);
  VecAssemblyEnd(user->Ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
 
  return(0);
}

////////////////////////////////////////////////////////// binary
PetscErrorCode VTKOut(UserCtx *user) {
  
  const char *byte_order = "LittleEndian";
  const char precision[] = "Float64";
  MPI_Comm                 comm;
  FILE                     *fp;
  PetscMPIInt              rank,size,tag;
  PetscInt bi;
  PetscInt                 dim,mx,my,mz,cdim,bs,boffset,maxnnodes,maxbs,i,j,k,r,bytes,nnodes,dsize;
  PetscScalar              *array,*array2;
  PetscErrorCode           ierr;
  
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  MPI_Type_size(MPI_DOUBLE,&dsize);
  
  for (bi=0; bi<block_number; bi++) {
    DM da = user[bi].da, fda = user[bi].fda;
    DMDALocalInfo info = user[bi].info;
    
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    PetscInt	lxs, lys, lzs, lxe, lye, lze;
    PetscInt	i, j, k, n;
    PetscReal	***aj, x, y, z;
    Cmpnts	***ucat, ***coor, ***ucat_o, ***uxx,***uxy;
    PetscReal	***p, ***nvert,***nvert1,***nvert2,***nvert3, ***p_o, ***nvert_o;
    Vec		Coor, zCoor, nCoor;
    
    PetscReal theta_x=0.0;
    PetscOptionsGetReal(PETSC_NULL, "-theta", &theta_x, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD,"theta %le \n",theta_x);
    PetscInt Pres=0,airfoil=0;
    
    PetscOptionsGetInt(PETSC_NULL, "-Pres", &Pres, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-airfoil", &airfoil, PETSC_NULL);
    
    if (!Pres){
      nnodes = (mx-1)*(my-1)*(mz-1);
      boffset   = 0;
      
      cdim=3; //dimension of the grid
      //////////////////////
      char filen[80];
      sprintf(filen, "results/Result%2.2d_%5.5d.vts", bi,ti);
      
      PetscFOpen(MPI_COMM_WORLD,filen,"wb",&fp);
      PetscFPrintf(MPI_COMM_WORLD,fp,"<?xml version=\"1.0\"?>\n");
      PetscFPrintf(MPI_COMM_WORLD,fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",byte_order);
      PetscFPrintf(MPI_COMM_WORLD,fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",0,mx-2,0,my-2,0,mz-2);
      PetscFPrintf(MPI_COMM_WORLD,fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",0,mx-2,0,my-2,0,mz-2);
      ////////////////
      PetscFPrintf(MPI_COMM_WORLD,fp,"      <Points>\n");
      PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%D\" />\n",precision,boffset);
      boffset += 3*nnodes*sizeof(PetscScalar) + sizeof(int);
      PetscFPrintf(MPI_COMM_WORLD,fp,"      </Points>\n");
      ////////////////
      PetscFPrintf(MPI_COMM_WORLD,fp,"      <PointData Scalars=\"ScalarPointData\">\n");
      //////
      PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",precision,"nvert",1 ,boffset);
      boffset += nnodes*sizeof(PetscScalar) + sizeof(int);
      ////
      PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%D\" />\n",precision,"p",boffset);
      boffset += nnodes*sizeof(PetscScalar) + sizeof(int);
      ////////
      PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",precision,"u",3 ,boffset);
      boffset += 3*nnodes*sizeof(PetscScalar) + sizeof(int);
      ////////////////////////////////////
      PetscFPrintf(MPI_COMM_WORLD,fp,"      </PointData>\n");
      PetscFPrintf(MPI_COMM_WORLD,fp,"    </Piece>\n");
      ///
      PetscFPrintf(MPI_COMM_WORLD,fp,"  </StructuredGrid>\n");
      PetscFPrintf(MPI_COMM_WORLD,fp,"  <AppendedData encoding=\"raw\">\n");
      PetscFPrintf(MPI_COMM_WORLD,fp,"_");
      
      /* Now write the arrays. */
      PetscMalloc2(nnodes,&array,nnodes*3,&array2);
      
      /////////////////////////////////////////////////////// coordinate
      DMGetCoordinates(da, &Coor);
      DMDAVecGetArray(fda, Coor, &coor);
      
      for (k=0; k<mz-1; k++) {
	for (j=0; j<my-1; j++) {
	  for (i=0; i<mx-1; i++) {
	    PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
	    
	    x=coor[k][j][i].x;
	    y=coor[k][j][i].y;
	    
	    if (bi==2)
	      {
		x=coor[k][j][i].x*cos(theta_x)-coor[k][j][i].y*sin(theta_x);
		y=coor[k][j][i].y*cos(theta_x)+coor[k][j][i].x*sin(theta_x);
	      }
	    
	    array2[Iloc*3+0] = x;
	    array2[Iloc*3+1] = y;
	    array2[Iloc*3+2] = coor[k][j][i].z;
	  }
	}
      }
      
      bytes=3*nnodes*dsize;
      fwrite(&bytes,sizeof(int),1,fp);
      fwrite(array2,dsize,3*nnodes,fp);
      
      DMDAVecRestoreArray(fda, Coor, &coor);
      
      /////////////////////////////////////////////////////// nvert
      DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
      
      for (k=0; k<mz-1; k++) {
	for (j=0; j<my-1; j++) {
	  for (i=0; i<mx-1; i++) {
	    PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
//	    array[Iloc] = nvert[k][j][i];  // non-averaged nvert
	    array[Iloc] = 0.125 *
	   (nvert[k][j][i] + nvert[k][j][i+1] +
	    nvert[k][j+1][i] + nvert[k][j+1][i+1] +
	    nvert[k+1][j][i] + nvert[k+1][j][i+1] +
	    nvert[k+1][j+1][i] + nvert[k+1][j+1][i+1]);
	  }
	}
      }
      
      bytes=nnodes*dsize;
      fwrite(&bytes,sizeof(int),1,fp);
      fwrite(array,dsize,nnodes,fp);
      
      DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
      
      /////////////////////////////////////////////////////// p
      DMDAVecGetArray(user[bi].da, user[bi].P, &p);

      for (k=0; k<mz-1; k++) {
	for (j=0; j<my-1; j++) {
	  for (i=0; i<mx-1; i++) {
	    PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
	    array[Iloc] = 0.125 *
	   (p[k][j][i] + p[k][j][i+1] +
	    p[k][j+1][i] + p[k][j+1][i+1] +
	    p[k+1][j][i] + p[k+1][j][i+1] +
	    p[k+1][j+1][i] + p[k+1][j+1][i+1]);
	  }
	}
      }
      
      bytes=nnodes*dsize;
      fwrite(&bytes,sizeof(int),1,fp);
      fwrite(array,dsize,nnodes,fp);
      
      DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
      
      /////////////////////////////////////////////////////// velocity
      DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
      DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);
      
      for (k=0; k<mz-1; k++) {
	for (j=0; j<my-1; j++) {
	  for (i=0; i<mx-1; i++) {
	    ucat_o[k][j][i].x =  0.125 *
	      (ucat[k][j][i].x + ucat[k][j][i+1].x +
	       ucat[k][j+1][i].x + ucat[k][j+1][i+1].x +
	       ucat[k+1][j][i].x + ucat[k+1][j][i+1].x +
	       ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);
	    ucat_o[k][j][i].y = 0.125 *
	      (ucat[k][j][i].y + ucat[k][j][i+1].y +
	       ucat[k][j+1][i].y + ucat[k][j+1][i+1].y +
	       ucat[k+1][j][i].y + ucat[k+1][j][i+1].y +
	       ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y);
	    ucat_o[k][j][i].z = 0.125 *
	      (ucat[k][j][i].z + ucat[k][j][i+1].z +
	       ucat[k][j+1][i].z + ucat[k][j+1][i+1].z +
	       ucat[k+1][j][i].z + ucat[k+1][j][i+1].z +
	       ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z);
	    //////////
	    
	    x=ucat_o[k][j][i].x;
	    y=ucat_o[k][j][i].y;
	    
	    if (bi==2){
	    x=ucat_o[k][j][i].x*cos(theta_x)-ucat_o[k][j][i].y*sin(theta_x);
	    y=ucat_o[k][j][i].y*cos(theta_x)+ucat_o[k][j][i].x*sin(theta_x);
	    }
	    
	    PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
	    array2[Iloc*3+0] = x;
	    array2[Iloc*3+1] = y;
	    array2[Iloc*3+2] = ucat_o[k][j][i].z;
	  }
	}
      }
      
      bytes=3*nnodes*dsize;
      fwrite(&bytes,sizeof(int),1,fp);
      fwrite(array2,dsize,3*nnodes,fp);
      
      DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
      DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);
      
      PetscFPrintf(MPI_COMM_WORLD,fp,"\n </AppendedData>\n");
      PetscFPrintf(MPI_COMM_WORLD,fp,"</VTKFile>\n");
      PetscFClose(MPI_COMM_WORLD,fp);
      
      PetscFree2(array,array2);
    }else{
      
      PetscPrintf(PETSC_COMM_WORLD, " PL \n ");
      double pl,pu;
      /*	FILE *gp;
		char filet[80];
		gp = fopen(filet, "w");
		sprintf(filet, "PPPP.dat");
		PetscFPrintf(PETSC_COMM_WORLD, gp," PL \n ");
      */
      if (bi==2  || Pres || airfoil){
	
	DMDAVecGetArray(user[bi].da, user[bi].P, &p);
	DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
	DMGetCoordinates(da, &Coor);
	DMDAVecGetArray(fda, Coor, &coor);
	
	if ( bi==2)
	  {
	    i=20;
	    for (j=ys+1; j<ye-2; j++) {
	      for (k=zs+1; k<ze-2; k++) {
	 if (nvert[k][j][i]==0.0 && nvert[k+1][j][i]>0.0){
	   pl=p[k][j][i];
	   PetscPrintf(PETSC_COMM_WORLD,  "%f %f  %f  %f\n ",coor[k][j][i].x,coor[k][j][i].z,coor[k][j][i].y,pl);
	 }
	      }
	    }
	    PetscPrintf(PETSC_COMM_WORLD, " PU \n ");
	    
	    for (j=ys+1; j<ye-2; j++) {
	      for (k=zs+1; k<ze-2; k++) {
		
		if (nvert[k][j][i]==0.0 && nvert[k-1][j][i]>0.0){
		  pu=p[k][j][i];
		  PetscPrintf(PETSC_COMM_WORLD,  "%f  %f  %f   \n ",coor[k][j][i].z,coor[k][j][i].y,pu);
		}
	      }
	    }
	  }
	if (airfoil) {//for other cases
	  
	  i=8;
	  for (k=zs+1; k<ze-2; k++) {
	    for (j=ys+1; j<ye-2; j++) {
	      if (nvert[k][j][i]==0.0 && nvert[k][j+1][i]>0.0){
		pl=p[k][j][i];
		PetscPrintf(PETSC_COMM_WORLD,  "%f %f  %f  %f\n ",coor[k][j][i].x,coor[k][j][i].z,coor[k][j][i].y,pl);
	      }
	    }
	  }
	  PetscPrintf(PETSC_COMM_WORLD, " PU \n ");
	  
	  for (k=zs+1; k<ze-2; k++) {
	    for (j=ys+1; j<ye-2; j++) {
	 
	      if (nvert[k][j][i]==0.0 && nvert[k][j-1][i]>0.0){
		pu=p[k][j][i];
		PetscPrintf(PETSC_COMM_WORLD,  "%f  %f  %f   \n ",coor[k][j][i].z,coor[k][j][i].y,pu);
	      }
	    }
	  }
	  
	}

	DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
	DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
	DMDAVecRestoreArray(fda, Coor, &coor);
	
	//fclose(gp);
	
      }
    }
    
  }
 
  return 0;
}
//-----------------------------------------------------------------------------------------------
PetscErrorCode VTKOutShr(UserCtx *user)
{
  const char *byte_order = "LittleEndian";
  const char precision[] = "Float64";
  MPI_Comm                 comm;
  FILE                     *fp;
  PetscMPIInt              rank,size,tag;
  PetscInt bi;
  PetscInt                 dim,mx,my,mz,cdim,bs,boffset,maxnnodes,maxbs,i,j,k,r,bytes,nnodes,dsize;
  PetscScalar              *array,*array2;
  PetscErrorCode           ierr;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  
  MPI_Type_size(MPI_DOUBLE,&dsize);
  
  for (bi=0; bi<block_number; bi++) {
    DM da = user[bi].da, fda = user[bi].fda;
    DMDALocalInfo info = user[bi].info;
    
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    PetscInt	lxs, lys, lzs, lxe, lye, lze;
    PetscInt	i, j, k, n;
    PetscReal	***aj, x, y, z;
    Cmpnts	***ucat, ***coor, ***ucat_o, ***uxx,***uxy;
    PetscReal	***p, ***nvert,***nvert1,***nvert2,***nvert3, ***p_o, ***nvert_o, ***q;
    Vec		Coor, zCoor, nCoor;

  PetscReal theta_x=0.0;
  PetscOptionsGetReal(PETSC_NULL, "-theta", &theta_x, PETSC_NULL); 
  PetscPrintf(PETSC_COMM_WORLD,"theta %le \n",theta_x);
  PetscInt Pres=0,airfoil=0;

  PetscOptionsGetInt(PETSC_NULL, "-Pres", &Pres, PETSC_NULL); 
  PetscOptionsGetInt(PETSC_NULL, "-airfoil", &airfoil, PETSC_NULL); 
 

if (!Pres){   
    nnodes = (mx-1)*(my-1)*(mz-1);
    boffset   = 0;

    cdim=3; //dimension of the grid
    //////////////////////
    char filen[80];
    sprintf(filen, "Shear%2.2d_%5.5d.vts", bi,ti);
    
    PetscFOpen(MPI_COMM_WORLD,filen,"wb",&fp);
    PetscFPrintf(MPI_COMM_WORLD,fp,"<?xml version=\"1.0\"?>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",byte_order);
    PetscFPrintf(MPI_COMM_WORLD,fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",0,mx-2,0,my-2,0,mz-2);
    PetscFPrintf(MPI_COMM_WORLD,fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",0,mx-2,0,my-2,0,mz-2);
    ////////////////
    PetscFPrintf(MPI_COMM_WORLD,fp,"      <Points>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%D\" />\n",precision,boffset);
    boffset += 3*nnodes*sizeof(PetscScalar) + sizeof(int);
    PetscFPrintf(MPI_COMM_WORLD,fp,"      </Points>\n");
    ////////////////
    PetscFPrintf(MPI_COMM_WORLD,fp,"      <PointData Scalars=\"ScalarPointData\">\n");
    //////
    PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",precision,"nvert",1 ,boffset);
    boffset += nnodes*sizeof(PetscScalar) + sizeof(int);
    ////
    PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%D\" />\n",precision,"shear",boffset);
    boffset += nnodes*sizeof(PetscScalar) + sizeof(int);
    ////////
    /* PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",precision,"u",3 ,boffset); */
    /* boffset += 3*nnodes*sizeof(PetscScalar) + sizeof(int); */
    ////////////////////////////////////
    PetscFPrintf(MPI_COMM_WORLD,fp,"      </PointData>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"    </Piece>\n");
    ///
    PetscFPrintf(MPI_COMM_WORLD,fp,"  </StructuredGrid>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"  <AppendedData encoding=\"raw\">\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"_");


    /* Now write the arrays. */
    PetscMalloc2(nnodes,&array,nnodes*3,&array2);
    
    /////////////////////////////////////////////////////// coordinate
    DMGetCoordinates(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);
    
    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
	for (i=0; i<mx-1; i++) {
	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
	  
	  x=coor[k][j][i].x;
	  y=coor[k][j][i].y;
	  
	  if (bi==2)
	    {
	      x=coor[k][j][i].x*cos(theta_x)-coor[k][j][i].y*sin(theta_x);
	      y=coor[k][j][i].y*cos(theta_x)+coor[k][j][i].x*sin(theta_x);
	    }
	  
	  array2[Iloc*3+0] = x;
	  array2[Iloc*3+1] = y;
	  array2[Iloc*3+2] = coor[k][j][i].z;
	}
      }
    }
    
    bytes=3*nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array2,dsize,3*nnodes,fp);
    
    DMDAVecRestoreArray(fda, Coor, &coor);
    
    /////////////////////////////////////////////////////// nvert
    DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);

    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
    	for (i=0; i<mx-1; i++) {
    	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
    	  array[Iloc] = nvert[k][j][i];
    	}
      }
    }

    bytes=nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array,dsize,nnodes,fp);
    
    DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
    
    /////////////////////////////////////////////////////// p
    DMDAVecGetArray(user[bi].da, user[bi].shrr, &q);

    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
    	for (i=0; i<mx-1; i++) {
    	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
    	  array[Iloc] = q[k][j][i];
    	}
      }
    }

    bytes=nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array,dsize,nnodes,fp);
    
    DMDAVecRestoreArray(user[bi].da, user[bi].shrr, &q);
        
 
      /////////////////////////////////////////////////////// velocity
    /* DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
    /* DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */

    /* for (k=0; k<mz-1; k++) { */
    /*   for (j=0; j<my-1; j++) { */
    /* 	for (i=0; i<mx-1; i++) { */
    /* 	  ucat_o[k][j][i].x = 0.125 * */
    /* 	    (ucat[k][j][i].x + ucat[k][j][i+1].x + */
    /* 	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
    /* 	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
    /* 	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x); */
    /* 	  ucat_o[k][j][i].y = 0.125 * */
    /* 	    (ucat[k][j][i].y + ucat[k][j][i+1].y + */
    /* 	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
    /* 	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
    /* 	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
    /* 	  ucat_o[k][j][i].z = 0.125 * */
    /* 	    (ucat[k][j][i].z + ucat[k][j][i+1].z + */
    /* 	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
    /* 	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
    /* 	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */
    /* 	  ////////// */

    /* 	   x=ucat_o[k][j][i].x; */
    /* 	  y=ucat_o[k][j][i].y; */

    /* 	  if (bi==2){ */
    /* 	    x=ucat_o[k][j][i].x*cos(theta_x)-ucat_o[k][j][i].y*sin(theta_x); */
    /* 	    y=ucat_o[k][j][i].y*cos(theta_x)+ucat_o[k][j][i].x*sin(theta_x); */
    /* 	  } */
	  
    /* 	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k); */
    /* 	  array2[Iloc*3+0] = x; */
    /* 	  array2[Iloc*3+1] = y; */
    /* 	  array2[Iloc*3+2] = ucat_o[k][j][i].z; */
    /* 	} */
    /*   } */
    /* } */

    /* bytes=3*nnodes*dsize; */
    /* fwrite(&bytes,sizeof(int),1,fp); */
    /* fwrite(array2,dsize,3*nnodes,fp); */
    
    /* DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat); */
    /* DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */

    
    
    PetscFPrintf(MPI_COMM_WORLD,fp,"\n </AppendedData>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"</VTKFile>\n");
    PetscFClose(MPI_COMM_WORLD,fp);
    
    PetscFree2(array,array2);
 }else{

 
    PetscPrintf(PETSC_COMM_WORLD, " PL \n ");
   	double pl,pu;
     /*	FILE *gp;
	char filet[80]; 
	gp = fopen(filet, "w"); 
	sprintf(filet, "PPPP.dat");
	PetscFPrintf(PETSC_COMM_WORLD, gp," PL \n ");
   */
   if (bi==2  || Pres || airfoil){


   DMDAVecGetArray(user[bi].da, user[bi].P, &p);
   DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
   DMGetCoordinates(da, &Coor);
   DMDAVecGetArray(fda, Coor, &coor);
 
 if ( bi==2)   
  { 
   i=20;
     for (j=ys+1; j<ye-2; j++) {   
       for (k=zs+1; k<ze-2; k++) {	 
	 if (nvert[k][j][i]==0.0 && nvert[k+1][j][i]>0.0){
	   pl=p[k][j][i];
	   PetscPrintf(PETSC_COMM_WORLD,  "%f %f  %f  %f\n ",coor[k][j][i].x,coor[k][j][i].z,coor[k][j][i].y,pl);
	 }
       }
     }
      PetscPrintf(PETSC_COMM_WORLD, " PU \n ");

     for (j=ys+1; j<ye-2; j++) {   
       for (k=zs+1; k<ze-2; k++) {
	 
	 if (nvert[k][j][i]==0.0 && nvert[k-1][j][i]>0.0){
	   pu=p[k][j][i];
	   PetscPrintf(PETSC_COMM_WORLD,  "%f  %f  %f   \n ",coor[k][j][i].z,coor[k][j][i].y,pu);
	 }
       }
     }
  }
 if (airfoil) {//for other cases
 
  i=8; 
   for (k=zs+1; k<ze-2; k++) {	 	
     for (j=ys+1; j<ye-2; j++) {   
      if (nvert[k][j][i]==0.0 && nvert[k][j+1][i]>0.0){
	   pl=p[k][j][i];
	   PetscPrintf(PETSC_COMM_WORLD,  "%f %f  %f  %f\n ",coor[k][j][i].x,coor[k][j][i].z,coor[k][j][i].y,pl);
	 }
       }
     }
      PetscPrintf(PETSC_COMM_WORLD, " PU \n ");
  
       for (k=zs+1; k<ze-2; k++) {
        for (j=ys+1; j<ye-2; j++) { 
	 
	 if (nvert[k][j][i]==0.0 && nvert[k][j-1][i]>0.0){
	   pu=p[k][j][i];
	   PetscPrintf(PETSC_COMM_WORLD,  "%f  %f  %f   \n ",coor[k][j][i].z,coor[k][j][i].y,pu);
	 }
       }
     }

   }

   DMDAVecRestoreArray(user[bi].da, user[bi].P, &p);
   DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
   DMDAVecRestoreArray(fda, Coor, &coor);
   
   //fclose(gp); 

   }
 }
 
}
 



 
  return 0;
}

////////////////////////////////////////////////////////// ASCII output

/* PetscErrorCode VTKOut(UserCtx *user) */
/* { */
/*   PetscInt bi; */
/*   FILE *f; */

/*   PetscInt	rans = 0; */
/*   PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL); */
/*   PetscOptionsGetInt(PETSC_NULL, "-vx", &vx, PETSC_NULL); */
/*   PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */
    
/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k, n; */
/*     PetscReal	***aj, x, y, z; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o, ***uxx,***uxy; */
/*     PetscReal	***p, ***nvert,***nvert1,***nvert2,***nvert3, ***p_o, ***nvert_o; */
/*     Vec		Coor, zCoor, nCoor; */
/*     VecScatter	ctx; */

/*     DMGetCoordinates(da, &Coor); */
/*     DMDAVecGetArray(fda, Coor, &coor); */
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
/*     DMDAVecGetArray(user[bi].da, user[bi].P, &p); */
/*     DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */

    
/*     Cmpnts2 ***komega; */
/*     if (rans) */
/*     DMDAVecGetArray(user[bi].fda2, user[bi].K_Omega, &komega); */

/*     if (rans) */
/*       n=10; */
/*     else */
/*       n=8; */
    
/*     PetscInt	LOC[n] ;//  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0}; */
/*     /\* 1 is cell-centered */
/*        0 is node centered *\/ */
/*     for (i=0;i<n;i++){ */
/*       if (i<6) */
/* 	LOC[i]=1; */
/*       else */
/* 	LOC[i]=0; */
/*     } */

/*     PetscInt	vc = 0; */
/*     PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL); */

/*     PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL); */
/*     if (vc) { */
/*       LOC[3]=0; LOC[4]=0; LOC[5]=0; */
/*     } */


/*     // vtk file name */
/*     char filen[80]; */
/*     sprintf(filen, "Result%2.2d_%5.5d.vtk", bi,ti); */

/*     f = fopen(filen, "w"); // open file */
/*     /\* */
/*     You can use PetscFPrintf function to wirte to file *f similar to fprintf in C */
/*     Exp.this writes ti to file f: */
   
/*        coor[k][j][i].x give the x of node i,j,k */
/*       .y and .z components give y and z of the that node */
/*       Exp. */
/*     *\/ */
/*     //amir for saving just few points not all points in the results */
/*     int ii; */
/*     PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE); */

        
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "Structured Grid\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET STRUCTURED_GRID\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "DIMENSIONS %d %d %d\n",mx-1,my-1,mz-1); */
    
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d double\n",(mx-1)*(my-1)*(mz-1)); */
/*     for (k=0; k<mz-1; k++) { */
/*       for (j=0; j<my-1; j++) { */
/* 	for (i=0; i<mx-1; i++) { */
/* 	  x = coor[k][j][i].x; */
/* 	  y = coor[k][j][i].y; */
/* 	  z = coor[k][j][i].z; */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n ",x,y,z); */
/* 	} */
/*       } */
/*     } */
   
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n",(mx-1)*(my-1)*(mz-1)); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u double\n"); */
		 
/*     for (k=0; k<mz-1; k++) { */
/*       for (j=0; j<my-1; j++) { */
/* 	for (i=0; i<mx-1; i++) { */
/* 	  ucat_o[k][j][i].x = 0.125 * */
/* 	    (ucat[k][j][i].x + ucat[k][j][i+1].x + */
/* 	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
/* 	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
/* 	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x); */
/* 	  ucat_o[k][j][i].y = 0.125 * */
/* 	    (ucat[k][j][i].y + ucat[k][j][i+1].y + */
/* 	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
/* 	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
/* 	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */
/* 	  x = ucat_o[k][j][i].x; */
/* 	  y = ucat_o[k][j][i].y; */
/* 	  z = ucat_o[k][j][i].z; */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n ",x,y,z); */
/* 	} */
/*       } */
/*     } */
    
/*    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS p double\n"); */
/*    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*    for (k=zs; k<ze-1; k++) { */
/*      for (j=ys; j<ye-1; j++) { */
/*        for (i=xs; i<xe-1; i++) { */
/* 	 x = 0.125 * */
/* 	   (p[k][j][i] + p[k][j][i+1] + */
/* 	    p[k][j+1][i] + p[k][j+1][i+1] + */
/* 	    p[k+1][j][i] + p[k+1][j][i+1] + */
/* 	    p[k+1][j+1][i] + p[k+1][j+1][i+1]); */
/* 	 // x= p[k+1][j+1][i+1]; */
/* 	 PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n ",x); */
/*        } */
/*      } */
/*    } */

/*    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS nvert double\n"); */
/*    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*    for (k=zs; k<ze-1; k++) { */
/*      for (j=ys; j<ye-1; j++) { */
/*        for (i=xs; i<xe-1; i++) { */
/* 	 /\* x = 0.125 * *\/ */
/* 	 /\*   (nvert[k][j][i] + nvert[k][j][i+1] + *\/ */
/* 	 /\*    nvert[k][j+1][i] + nvert[k][j+1][i+1] + *\/ */
/* 	 /\*    nvert[k+1][j][i] + nvert[k+1][j][i+1] + *\/ */
/* 	 /\*    nvert[k+1][j+1][i] + nvert[k+1][j+1][i+1]); *\/ */
/* 	  x = nvert[k+1][j+1][i+1]; */
/* 	 PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n ",x); */
/*        } */
/*      } */
/*    } */
    
      
  
/*     fclose(f); // close the file */

/*     DMDAVecRestoreArray(fda, Coor, &coor); */
/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].P, &p); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */

/*   } */

 
/*   return 0; */
/* } */

PetscErrorCode VTKOutQ(UserCtx *user)
{
  const char *byte_order = "LittleEndian";
  const char precision[] = "Float64";
  MPI_Comm                 comm;
  FILE                     *fp;
  PetscMPIInt              rank,size,tag;
  PetscInt bi;
  PetscInt                 dim,mx,my,mz,cdim,bs,boffset,maxnnodes,maxbs,i,j,k,r,bytes,nnodes,dsize;
  PetscScalar              *array,*array2;
  PetscErrorCode           ierr;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  
  MPI_Type_size(MPI_DOUBLE,&dsize);
  
  for (bi=0; bi<block_number; bi++) {
    DM da = user[bi].da, fda = user[bi].fda;
    DMDALocalInfo info = user[bi].info;
    
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    PetscInt	lxs, lys, lzs, lxe, lye, lze;
    PetscInt	i, j, k, n;
    PetscReal	***aj, x, y, z;
    Cmpnts	***ucat, ***coor, ***ucat_o, ***uxx,***uxy;
    PetscReal	***q, ***nvert,***nvert1,***nvert2,***nvert3, ***p_o, ***nvert_o;
    Vec		Coor, zCoor, nCoor;

    nnodes = (mx-1)*(my-1)*(mz-1);
    boffset   = 0;

    cdim=3; //dimension of the grid
    //////////////////////
    char filen[80];
  PetscReal theta_x=0.0;
  PetscOptionsGetReal(PETSC_NULL, "-theta", &theta_x, PETSC_NULL); 
  PetscPrintf(PETSC_COMM_WORLD,"theta %le \n",theta_x);
 
 
    sprintf(filen, "QCriteria%2.2d_%5.5d.vts", bi,ti);
       
    PetscFOpen(MPI_COMM_WORLD,filen,"wb",&fp);
    PetscFPrintf(MPI_COMM_WORLD,fp,"<?xml version=\"1.0\"?>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",byte_order);
    PetscFPrintf(MPI_COMM_WORLD,fp,"  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",0,mx-2,0,my-2,0,mz-2);
    PetscFPrintf(MPI_COMM_WORLD,fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",0,mx-2,0,my-2,0,mz-2);
    ////////////////
    PetscFPrintf(MPI_COMM_WORLD,fp,"      <Points>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%D\" />\n",precision,boffset);
    boffset += 3*nnodes*sizeof(PetscScalar) + sizeof(int);
    PetscFPrintf(MPI_COMM_WORLD,fp,"      </Points>\n");
    ////////////////
    PetscFPrintf(MPI_COMM_WORLD,fp,"      <PointData Scalars=\"ScalarPointData\">\n");
    //////
    PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",precision,"nvert",1 ,boffset);
    boffset += nnodes*sizeof(PetscScalar) + sizeof(int);
    ////
    PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%D\" />\n",precision,"q",boffset);
    boffset += nnodes*sizeof(PetscScalar) + sizeof(int);
    ////////
    //   PetscFPrintf(MPI_COMM_WORLD,fp,"        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",precision,"u",3 ,boffset);
    //boffset += 3*nnodes*sizeof(PetscScalar) + sizeof(int);
    ////////////////////////////////////
    PetscFPrintf(MPI_COMM_WORLD,fp,"      </PointData>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"    </Piece>\n");
    ///
    PetscFPrintf(MPI_COMM_WORLD,fp,"  </StructuredGrid>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"  <AppendedData encoding=\"raw\">\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"_");


    /* Now write the arrays. */
    PetscMalloc2(nnodes,&array,nnodes*3,&array2);
    
    /////////////////////////////////////////////////////// coordinate
    DMGetCoordinates(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);
   
    
    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
	for (i=0; i<mx-1; i++) {


	x=coor[k][j][i].x;
        y=coor[k][j][i].y;
	  
	if (bi==2)
	    {
	      x=coor[k][j][i].x*cos(theta_x)-coor[k][j][i].y*sin(theta_x);
	      y=coor[k][j][i].y*cos(theta_x)+coor[k][j][i].x*sin(theta_x);
	    }
	

	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
	  array2[Iloc*3+0] = x;
	  array2[Iloc*3+1] = y;
	  array2[Iloc*3+2] = coor[k][j][i].z;
	}
      }
    }
    
    bytes=3*nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array2,dsize,3*nnodes,fp);
    
    DMDAVecRestoreArray(fda, Coor, &coor);
    
    /////////////////////////////////////////////////////// nvert
    DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);

    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
    	for (i=0; i<mx-1; i++) {
    	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
    	  array[Iloc] = nvert[k][j][i];
    	}
      }
    }

    bytes=nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array,dsize,nnodes,fp);
    
    DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
    
    /////////////////////////////////////////////////////// p
    DMDAVecGetArray(user[bi].da, user[bi].P, &q);

    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
    	for (i=0; i<mx-1; i++) {
    	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
    	  array[Iloc] = q[k][j][i];
    	}
      }
    }

    bytes=nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array,dsize,nnodes,fp);
    
    DMDAVecRestoreArray(user[bi].da, user[bi].P, &q);
    
    /////////////////////////////////////////////////////// velocity
    /* DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);

    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
    	for (i=0; i<mx-1; i++) {
    	  ucat_o[k][j][i].x = 0.125 *
    	    (ucat[k][j][i].x + ucat[k][j][i+1].x +
    	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x +
    	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x +
    	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);
    	  ucat_o[k][j][i].y = 0.125 *
    	    (ucat[k][j][i].y + ucat[k][j][i+1].y +
    	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y +
    	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y +
    	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y);
    	  ucat_o[k][j][i].z = 0.125 *
    	    (ucat[k][j][i].z + ucat[k][j][i+1].z +
    	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z +
    	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z +
    	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z);
    	  //////////
    	  PetscInt Iloc = i+(mx-1)*(j+(my-1)*k);
    	  array2[Iloc*3+0] = ucat_o[k][j][i].x;
    	  array2[Iloc*3+1] = ucat_o[k][j][i].y;
    	  array2[Iloc*3+2] = ucat_o[k][j][i].z;
    	}
      }
    }

    bytes=3*nnodes*dsize;
    fwrite(&bytes,sizeof(int),1,fp);
    fwrite(array2,dsize,3*nnodes,fp);
    
    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);

    */
    
    PetscFPrintf(MPI_COMM_WORLD,fp,"\n </AppendedData>\n");
    PetscFPrintf(MPI_COMM_WORLD,fp,"</VTKFile>\n");
    PetscFClose(MPI_COMM_WORLD,fp);
    
    PetscFree2(array,array2);
  }
  
  return 0;
}





/* PetscErrorCode VTKOutQ(UserCtx *user) */
/* { */
/*   PetscInt bi; */
/*   FILE *f; */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */
    
/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k, n; */
/*     PetscReal	***aj, x, y, z; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o; */
/*     PetscReal	***q, ***nvert, ***l2; */
/*     Vec		Coor, zCoor, nCoor; */
/*     VecScatter	ctx; */

/*     PetscInt	vc = 0; */
/*     PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL); */

/*     lxs = xs; lxe = xe; */
/*     lys = ys; lye = ye; */
/*     lzs = zs; lze = ze; */
    
/*     if (xs==0 && vc) lxs = xs+1; */
/*     if (ys==0 && vc) lys = ys+1; */
/*     if (zs==0 && vc) lzs = zs+1; */
    

/*     // vtk file name */
/*     char filen[80]; */
/*     sprintf(filen, "QCriteria%2.2d_%5.5d.vtk", bi,ti); */

/*     f = fopen(filen, "w"); // open file */
/*     /\* */
/*     You can use PetscFPrintf function to wirte to file *f similar to fprintf in C  */
/*     Exp.this writes ti to file f: */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d ", ti);  */
/*     *\/ */
/*     // PetscFPrintf(PETSC_COMM_WORLD, f, "%d ", ti);  */

/*     /\*  */
/*        coor[k][j][i].x give the x of node i,j,k */
/*       .y and .z components give y and z of the that node */
/*       Exp. */
/*     *\/ */
/*     /\*for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x = coor[k][j][i].x; */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", x) */
/* 	    } */
/*       } */
/*     }*\/ */
/*     PetscInt ii; */
/*     PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE); */
/*     PetscOptionsGetInt(PETSC_NULL, "-Ddata", &Ddata, PETSC_NULL); */
/*     PetscOptionsGetInt(PETSC_NULL, "-i_data", &ii, PETSC_NULL);   */
/*     PetscOptionsGetInt(PETSC_NULL, "-qv", &qv, PETSC_NULL); */
/*     PetscOptionsGetInt(PETSC_NULL, "-vx", &vx, PETSC_NULL); */
/*     if (Ddata==1) { */
/*       mx=2; */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "Structured Grid\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET STRUCTURED_GRID\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "DIMENSIONS %d %d %d\n",mx-1,my-1,mz-1); */
      
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d double\n",(mx-1)*(my-1)*(mz-1)); */

/*       DMGetCoordinates(da, &Coor); */
/*       DMDAVecGetArray(fda, Coor, &coor); */
/*       for (k=0; k<mz-1; k++) { */
/* 	for (j=0; j<my-1; j++) { */
/* 	  for (i=ii-1; i<ii+mx-2; i++) { */
/* 	    x = coor[k][j][i].x; */
/* 	    y = coor[k][j][i].y; */
/* 	    z = coor[k][j][i].z; */
/* 	    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",x,y,z); */
/* 	  } */
/* 	} */
/*       } */
/*       DMDAVecRestoreArray(fda, Coor, &coor); */
/*             qv=1; */
/*       if (qv==1) */
/* 	{    */
/* 	  DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
/* 	  DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);  */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n",(mx-1)*(my-1)*(mz-1)); */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u double\n"); */
/* 	  for (k=0; k<mz-1; k++) { */
/* 	    for (j=0; j<my-1; j++) { */
/* 	      for (i=ii-1; i<ii+mx-2; i++) { */
/* 		ucat_o[k][j][i].x = 0.125 * */
/* 		  (ucat[k][j][i].x + ucat[k][j][i+1].x + */
/* 		   ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
/* 		   ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
/* 		   ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);	   */
/* 		ucat_o[k][j][i].y = 0.125 * */
/* 		  (ucat[k][j][i].y + ucat[k][j][i+1].y + */
/* 		   ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
/* 		   ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
/* 		   ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
/* 		ucat_o[k][j][i].z = 0.125 * */
/* 		  (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 		   ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 		   ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 		   ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */
/* 		x = ucat_o[k][j][i].x;  */
/* 			if (vx==1) x=0.0;  */
/* 		y = ucat_o[k][j][i].y; */
/* 		z = ucat_o[k][j][i].z; */
/* 		PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n ",x,y,z); */
/* 	      } */
/* 	    } */
/* 	  } */
	  
	  
/* 	} */
/*       //if (vc) PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n",(mx-2)*(my-2)*(mz-2));   */
/*       if (qv!=1) PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n",(mx-1)*(my-1)*(mz-1));   */
/*       // write Q */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS q double\n");     */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*       DMDAVecGetArray(user[bi].da, user[bi].P, &q); */
/*       for (k=lzs; k<ze-1; k++) { */
/* 	for (j=lys; j<ye-1; j++) { */
/* 	  for (i=ii-1; i<ii+mx-2; i++) { */
/* 	  if (vc) x= q[k][j][i];   */
/* 	  else x = 0.125 * */
/* 		 (q[k][j][i] + q[k][j][i+1] + */
/* 		  q[k][j+1][i] + q[k][j+1][i+1] + */
/* 		  q[k+1][j][i] + q[k+1][j][i+1] + */
/* 		  q[k+1][j+1][i] + q[k+1][j+1][i+1]); */
	  
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	  } */
/* 	} */
/*       } */
/*       DMDAVecRestoreArray(user[bi].da, user[bi].P, &q); */
      
/*       // write l2 */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS l2 double\n");     */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*       DMDAVecGetArray(user[bi].da, user[bi].Phi, &q); */
/*       for (k=lzs; k<ze-1; k++) { */
/* 	for (j=lys; j<ye-1; j++) { */
/* 	  for (i=ii-1; i<ii+mx-2; i++) { */
/* 	    if (vc) x= q[k][j][i];   */
/* 	    else x = 0.125 * */
/* 		   (q[k][j][i] + q[k][j][i+1] + */
/* 		    q[k][j+1][i] + q[k][j+1][i+1] + */
/* 		    q[k+1][j][i] + q[k+1][j][i+1] + */
/* 		    q[k+1][j+1][i] + q[k+1][j+1][i+1]);	    */
/* 	    PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	  } */
/* 	} */
/*       } */
      
/*       DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &q); */
      
/*     // write nvert */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS nvert double\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*       DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/*       for (k=lzs; k<ze-1; k++) { */
/* 	for (j=lys; j<ye-1; j++) { */
/* 	  for (i=ii-1; i<ii+mx-2; i++) { */
/* 	    if (vc) x= nvert[k][j][i];   */
/* 	    else x = 0.125 * */
/* 		   (nvert[k][j][i] + nvert[k][j][i+1] + */
/* 		    nvert[k][j+1][i] + nvert[k][j+1][i+1] + */
/* 		    nvert[k+1][j][i] + nvert[k+1][j][i+1] + */
/* 		    nvert[k+1][j+1][i] + nvert[k+1][j+1][i+1]);	   */
/* 	    PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	  } */
/* 	} */
/*       } */
     
/*       PetscPrintf(PETSC_COMM_WORLD,  " %f uz\n",ucat[10][10][0].z); */


/*     }else */
/*       { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n"); */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "Structured Grid\n"); */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n"); */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET STRUCTURED_GRID\n"); */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "DIMENSIONS %d %d %d\n",mx-1,my-1,mz-1); */
	
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d double\n",(mx-1)*(my-1)*(mz-1)); */
	
/* 	DMGetCoordinates(da, &Coor); */
/* 	DMDAVecGetArray(fda, Coor, &coor); */
/* 	for (k=0; k<mz-1; k++) { */
/* 	  for (j=0; j<my-1; j++) { */
/* 	    for (i=0; i<mx-1; i++) { */
/* 	      x = coor[k][j][i].x; */
/* 	      y = coor[k][j][i].y; */
/* 	      z = coor[k][j][i].z; */
/* 	      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",x,y,z); */
/* 	    } */
/* 	  } */
/* 	} */
/* 	DMDAVecRestoreArray(fda, Coor, &coor); */

/* 	if (vc) PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n",(mx-2)*(my-2)*(mz-2));   */
/* 	else PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n",(mx-1)*(my-1)*(mz-1));   */
/* 	// write Q  */
/* 	if (qv==1) */
/* 	{    */
/* 	  DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
/* 	  DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o);  */
/* 	  //PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n",(mx-1)*(my-1)*(mz-1)); */
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u double\n"); */
/* 	  for (k=0; k<mz-1; k++) { */
/* 	    for (j=0; j<my-1; j++) { */
/* 	      for (i=0; i<mx-1; i++) { */
/* 		ucat_o[k][j][i].x = 0.125 * */
/* 		  (ucat[k][j][i].x + ucat[k][j][i+1].x + */
/* 		   ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
/* 		   ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
/* 		   ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);	   */
/* 		ucat_o[k][j][i].y = 0.125 * */
/* 		  (ucat[k][j][i].y + ucat[k][j][i+1].y + */
/* 		   ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
/* 		   ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
/* 		   ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
/* 		ucat_o[k][j][i].z = 0.125 * */
/* 		  (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 		   ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 		   ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 		   ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */
/* 		x = ucat_o[k][j][i].x; */
/* 		y = ucat_o[k][j][i].y; */
/* 		z = ucat_o[k][j][i].z; */
/* 		PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n ",x,y,z); */
/* 	      } */
/* 	    } */
/* 	  } */

/* 	} */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS q double\n");     */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/* 	DMDAVecGetArray(user[bi].da, user[bi].P, &q); */
/* 	for (k=lzs; k<ze-1; k++) { */
/* 	  for (j=lys; j<ye-1; j++) { */
/* 	    for (i=lxs; i<xe-1; i++) { */
/* 	      if (vc) x= q[k][j][i];   */
/* 	      else x = 0.125 * */
/* 		 (q[k][j][i] + q[k][j][i+1] + */
/* 		  q[k][j+1][i] + q[k][j+1][i+1] + */
/* 		  q[k+1][j][i] + q[k+1][j][i+1] + */
/* 		  q[k+1][j+1][i] + q[k+1][j+1][i+1]); */
	      
/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	    } */
/* 	  } */
/* 	} */
/* 	DMDAVecRestoreArray(user[bi].da, user[bi].P, &q); */
	
/* 	// write l2 */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS l2 double\n");     */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/* 	DMDAVecGetArray(user[bi].da, user[bi].Phi, &q); */
/* 	for (k=lzs; k<ze-1; k++) { */
/* 	  for (j=lys; j<ye-1; j++) { */
/* 	    for (i=lxs; i<xe-1; i++) { */
/* 	   if (vc) x= q[k][j][i];   */
/* 	   else x = 0.125 * */
/* 	    (q[k][j][i] + q[k][j][i+1] + */
/* 	     q[k][j+1][i] + q[k][j+1][i+1] + */
/* 	     q[k+1][j][i] + q[k+1][j][i+1] + */
/* 	     q[k+1][j+1][i] + q[k+1][j+1][i+1]);	    */
/* 	   PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	} */
/*       } */
/*     } */

/* 	DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &q); */
	
/* 	// write nvert */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS nvert double\n"); */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/* 	DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/* 	for (k=lzs; k<ze-1; k++) { */
/* 	  for (j=lys; j<ye-1; j++) { */
/* 	    for (i=lxs; i<xe-1; i++) { */
/* 	      if (vc) x= nvert[k][j][i];   */
/* 	      else x = 0.125 * */
/* 		     (nvert[k][j][i] + nvert[k][j][i+1] + */
/* 		      nvert[k][j+1][i] + nvert[k][j+1][i+1] + */
/* 		      nvert[k+1][j][i] + nvert[k+1][j][i+1] + */
/* 		      nvert[k+1][j+1][i] + nvert[k+1][j+1][i+1]);	   */
/* 	      PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	    } */
/* 	  } */
/* 	} */
/*       } //else   */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */
    
/*     fclose(f); // close the file     */
    
/*   } */

/*   return 0; */
/* } */

PetscErrorCode VTKOutQ_Binary(UserCtx *user)
{
  PetscInt bi;
  FILE *f;

  for (bi=0; bi<block_number; bi++) {
    DM da = user[bi].da, fda = user[bi].fda;
    DMDALocalInfo info = user[bi].info;
    
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    PetscInt	lxs, lys, lzs, lxe, lye, lze;
    PetscInt	i, j, k, n;
    PetscReal	***aj, x, y, z;
    Cmpnts	***ucat, ***coor, ***ucat_o;
    PetscReal	***q, ***nvert, ***l2;
    Vec		Coor, zCoor, nCoor;
    VecScatter	ctx;
  
    PetscInt	vc = 0;
    PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL);

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0 && vc) lxs = xs+1;
    if (ys==0 && vc) lys = ys+1;
    if (zs==0 && vc) lzs = zs+1;
    

    // vtk file name
    char filen[80];
    sprintf(filen, "results/QCriteria%2.2d_%5.5d.vtk", bi,ti);

    f = fopen(filen, "wb"); // open file
     
   /*  PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 4.0\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "Structured Grid\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "BINARY\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET STRUCTURED_GRID\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "DIMENSIONS %d %d %d\n",mx-1,my-1,mz-1); */
  
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d double\n",(mx-1)*(my-1)*(mz-1)); */
    //  char xx[10]="ABCDEFGHIJ";
    // int xx=1;
   
    DMGetCoordinates(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);
    for (k=0; k<mz-1; k++) {
      for (j=0; j<my-1; j++) {
	for (i=0; i<mx-1; i++) {
	  x = coor[k][j][i].x;
	  y = coor[k][j][i].y;
	  z = coor[k][j][i].z;
	//  PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",x,y,z);
	  fwrite(&x,sizeof(double),1,f);
	  fwrite(&y,sizeof(double),1,f);
	  fwrite(&z,sizeof(double),1,f);
	}
      }
    }
    DMDAVecRestoreArray(fda, Coor, &coor);
 
   /*  if (vc) PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n",(mx-2)*(my-2)*(mz-2)); */
/*     else PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n",(mx-1)*(my-1)*(mz-1)); */
    // write Q
    //  double t;
  /*   PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS q double\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
    DMDAVecGetArray(user[bi].da, user[bi].P, &q);
    for (k=lzs; k<ze-1; k++) {
      for (j=lys; j<ye-1; j++) {
	for (i=lxs; i<xe-1; i++) {
	  if (vc) x= q[k][j][i];
	  else x = 0.125 *
	    (q[k][j][i] + q[k][j][i+1] +
	     q[k][j+1][i] + q[k][j+1][i+1] +
	     q[k+1][j][i] + q[k+1][j][i+1] +
	     q[k+1][j+1][i] + q[k+1][j+1][i+1]);
	   fwrite(&x,sizeof(double),1,f);
	  // PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x);
	}
      }
    }
    DMDAVecRestoreArray(user[bi].da, user[bi].P, &q);
 
   /*  // write l2 */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS l2 double\n");     */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*     DMDAVecGetArray(user[bi].da, user[bi].Phi, &q); */
/*     for (k=lzs; k<ze-1; k++) { */
/*       for (j=lys; j<ye-1; j++) { */
/* 	for (i=lxs; i<xe-1; i++) { */
/* 	   if (vc) x= q[k][j][i];   */
/* 	   else x = 0.125 * */
/* 	    (q[k][j][i] + q[k][j][i+1] + */
/* 	     q[k][j+1][i] + q[k][j+1][i+1] + */
/* 	     q[k+1][j][i] + q[k+1][j][i+1] + */
/* 	     q[k+1][j+1][i] + q[k+1][j+1][i+1]);	    */
/* 	   PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	} */
/*       } */
/*     } */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &q); */

  /*   // write nvert */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS nvert double\n"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n"); */
/*     DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/*     for (k=lzs; k<ze-1; k++) { */
/*       for (j=lys; j<ye-1; j++) { */
/* 	for (i=lxs; i<xe-1; i++) { */
/* 	   if (vc) x= nvert[k][j][i];   */
/* 	   else x = 0.125 * */
/* 	    (nvert[k][j][i] + nvert[k][j][i+1] + */
/* 	     nvert[k][j+1][i] + nvert[k][j+1][i+1] + */
/* 	     nvert[k+1][j][i] + nvert[k+1][j][i+1] + */
/* 	     nvert[k+1][j+1][i] + nvert[k+1][j+1][i+1]);	   */
/* 	   PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",x); */
/* 	} */
/*       } */
/*     } */
    
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */
    
    fclose(f); // close the file    
    
  }

  return 0;
}
/* PetscErrorCode TECIOOut(UserCtx *user) */
/* { */
/*   PetscInt bi; */

/*   char filen[80]; */
/*   sprintf(filen, "Result%3.3d.plt", ti); */

/*   INTEGER4  II, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax; */
/*   VIsDouble = 0; */
/*   DIsDouble = 0; */
/*   Debug = 0; */
  
/*   PetscInt	rans = 0; */
/*   PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL); */
/*   if (rans) //rans */
/*   II = TECINI100("Heart Valve", */
/* 	     "X Y Z U V W P Nv K O", */
/* 	     filen, */
/* 	     ".", */
/* 	     &Debug, */
/* 	     &VIsDouble); */
/*   else     */
/*   II = TECINI100("Heart Valve", */
/* 	     "X Y Z U V W P Nv", */
/* 	     filen, */
/* 	     ".", */
/* 	     &Debug, */
/* 	     &VIsDouble); */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */

/*     DM da1, da2; */

/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k, n; */
/*     PetscReal	***aj; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o; */
/*     PetscReal	***p, ***nvert; */
/*     Vec		Coor, zCoor, nCoor; */
/*     VecScatter	ctx; */

/*     Vec X, Y, Z, U, V, W; */

/*     INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0; */
/*     INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0; */
/*     INTEGER4    ShareConnectivityFromZone=0; */
/*     if (rans)  */
/*       n=10; */
/*     else */
/*       n=8; */
    
/*     INTEGER4	LOC[n] ;//  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};  */
/*     /\* 1 is cell-centered */
/*        0 is node centered *\/ */
/*     for (i=0;i<n;i++){ */
/*       if (i<6) */
/* 	LOC[i]=1; */
/*       else */
/* 	LOC[i]=0; */
/*     } */

/*     PetscInt	vc = 0; */
/*     PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL); */

/*     if (vc) { */
/*       LOC[3]=0; LOC[4]=0; LOC[5]=0; */
/*     } */
/*     float	*x; */
/*     PetscMalloc(mx*my*mz*sizeof(float), &x); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da1); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da2); */

/*     IMax = mx-1; JMax = my-1; KMax = mz-1; */

/*     II = TECZNE100("Block 1", */
/* 		  &ZoneType, 	/\* Ordered zone *\/ */
/* 		  &IMax, */
/* 		  &JMax, */
/* 		  &KMax, */
/* 		  &ICellMax, */
/* 		  &JCellMax, */
/* 		  &KCellMax, */
/* 		  &IsBlock,	/\* ISBLOCK  BLOCK format *\/ */
/* 		  &NumFaceConnections, */
/* 		  &FaceNeighborMode, */
/* 		  LOC, */
/* 		  NULL, */
/* 		  &ShareConnectivityFromZone); /\* No connectivity sharing *\/ */

/*     III = (mx-1) * (my-1) * (mz-1); */
/*     DMCreateGlobalVector(da1, &X); */

/*     DMGetCoordinates(da, &Coor); */
/*     DMDAVecGetArray(fda, Coor, &coor); */


/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(fda, Coor, &coor); */
    
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
/*     for (k=0; k<mz-1; k++) { */
/*       for (j=0; j<my-1; j++) { */
/* 	for (i=0; i<mx-1; i++) { */
/* 	  ucat_o[k][j][i].x = 0.125 * */
/* 	    (ucat[k][j][i].x + ucat[k][j][i+1].x + */
/* 	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
/* 	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
/* 	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x); */
/* 	  ucat_o[k][j][i].y = 0.125 * */
/* 	    (ucat[k][j][i].y + ucat[k][j][i+1].y + */
/* 	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
/* 	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
/* 	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */
/* 	} */
/*       } */
/*     } */

/*     if (!vc) { */
/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */



/*     } */
/*     else { */



/*       III = (mx-2) * (my-2) * (mz-2); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*     } */

/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */

/*     VecDestroy(&X); */


/*     III = (mx-2) * (my-2) * (mz-2); */
/*     DMCreateGlobalVector(da2, &X); */

/*     DMDAVecGetArray(user[bi].da, user[bi].P, &p); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1]; */

/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].P, &p); */

/*     DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1]; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */

/*     if (rans) { */
/*     Cmpnts2 ***komega; */
/*     DMDAVecGetArray(user[bi].fda2, user[bi].K_Omega, &komega); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = komega[k+1][j+1][i+1].x; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = komega[k+1][j+1][i+1].y; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     DMDAVecRestoreArray(user[bi].fda2, user[bi].K_Omega, &komega); */
/*     } */

/*     VecDestroy(&X); */

      

/*     PetscFree(x); */
/*     DMDestroy(&da1); */
/*     DMDestroy(&da2); */
/*   } */
/*   II = TECEND100(); */
/*   return 0; */
/* } */

/* PetscErrorCode TECIOOut_AVE(UserCtx *user) */
/* { */
/*   PetscInt bi; */

/*   char filen[80]; */
/*   sprintf(filen, "Result_AVE%5.5d.plt", ti); */

/*   INTEGER4 II, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax; */
/*   VIsDouble = 0; */
/*   DIsDouble = 0; */
/*   Debug = 0; */
  
/*   PetscInt	rans = 0; */
/*   PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL); */
  
/*   II = TECINI100("Heart Valve", */
/* 	     "X Y Z U V W UU VV WW UV VW WU P Nv", */
/* 	     filen, */
/* 	     ".", */
/* 	     &Debug, */
/* 	     &VIsDouble); */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */

/*     DM da1, da2; */

/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k, n; */
/*     PetscReal	***aj; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o, ***ucat_sq, ***ucat_cross; */
/*     PetscReal	***p, ***nvert; */
/*     Vec		Coor, zCoor, nCoor; */
/*     VecScatter	ctx; */

/*     Vec X, Y, Z, U, V, W; */

/*     INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0; */
/*     INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0; */
/*     INTEGER4    ShareConnectivityFromZone=0; */

/*     n=14; */
    
/*     INTEGER4	LOC[n] ;//  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};  */
/*     /\* 1 is cell-centered */
/*        0 is node centered *\/ */
/*     for (i=0;i<n;i++){ */
/*       if (i<n-2) */
/* 	LOC[i]=1; */
/*       else */
/* 	LOC[i]=0; */
/*     } */

/*     PetscInt	vc = 0; */
/*     PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL); */

/*     if (vc) { */
/*       for (i=3;i<13;i++) LOC[i]=0; */
/*     } */
/*     float	*x; */
/*     PetscMalloc(mx*my*mz*sizeof(float), &x); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da1); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da2); */

/*     IMax = mx-1; JMax = my-1; KMax = mz-1; */

/*     II = TECZNE100("Block 1", */
/* 		  &ZoneType, 	/\* Ordered zone *\/ */
/* 		  &IMax, */
/* 		  &JMax, */
/* 		  &KMax, */
/* 		  &ICellMax, */
/* 		  &JCellMax, */
/* 		  &KCellMax, */
/* 		  &IsBlock,	/\* ISBLOCK  BLOCK format *\/ */
/* 		  &NumFaceConnections, */
/* 		  &FaceNeighborMode, */
/* 		  LOC, */
/* 		  NULL, */
/* 		  &ShareConnectivityFromZone); /\* No connectivity sharing *\/ */

/*     III = (mx-1) * (my-1) * (mz-1); */
/*     DMCreateGlobalVector(da1, &X); */

/*     DMGetCoordinates(da, &Coor); */
/*     DMDAVecGetArray(fda, Coor, &coor); */


/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(fda, Coor, &coor); */
    
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
/*     for (k=0; k<mz-1; k++) { */
/*       for (j=0; j<my-1; j++) { */
/* 	for (i=0; i<mx-1; i++) { */
/* 	  ucat_o[k][j][i].x = 0.125 * */
/* 	    (ucat[k][j][i].x + ucat[k][j][i+1].x + */
/* 	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
/* 	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
/* 	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x); */
/* 	  ucat_o[k][j][i].y = 0.125 * */
/* 	    (ucat[k][j][i].y + ucat[k][j][i+1].y + */
/* 	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
/* 	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
/* 	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */
/* 	} */
/*       } */
/*     } */

/*     if (!vc) { */
/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */
/*     } */
/*     else { */
/*       III = (mx-2) * (my-2) * (mz-2); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*     } */

/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat); */

/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat_square_sum, &ucat_sq); */

/*     // UU VV WW */
/*     for (k=0; k<mz-1; k++) { */
/*       for (j=0; j<my-1; j++) { */
/* 	for (i=0; i<mx-1; i++) { */
/* 	  ucat_o[k][j][i].x = 0.125 * */
/* 	    (ucat_sq[k][j][i].x + ucat_sq[k][j][i+1].x + */
/* 	     ucat_sq[k][j+1][i].x + ucat_sq[k][j+1][i+1].x + */
/* 	     ucat_sq[k+1][j][i].x + ucat_sq[k+1][j][i+1].x + */
/* 	     ucat_sq[k+1][j+1][i].x + ucat_sq[k+1][j+1][i+1].x); */
/* 	  ucat_o[k][j][i].y = 0.125 * */
/* 	    (ucat_sq[k][j][i].y + ucat_sq[k][j][i+1].y + */
/* 	     ucat_sq[k][j+1][i].y + ucat_sq[k][j+1][i+1].y + */
/* 	     ucat_sq[k+1][j][i].y + ucat_sq[k+1][j][i+1].y + */
/* 	     ucat_sq[k+1][j+1][i].y + ucat_sq[k+1][j+1][i+1].y); */
/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat_sq[k][j][i].z + ucat_sq[k][j][i+1].z + */
/* 	     ucat_sq[k][j+1][i].z + ucat_sq[k][j+1][i+1].z + */
/* 	     ucat_sq[k+1][j][i].z + ucat_sq[k+1][j][i+1].z + */
/* 	     ucat_sq[k+1][j+1][i].z + ucat_sq[k+1][j+1][i+1].z); */
/* 	} */
/*       } */
/*     } */

/*     if (!vc) { */
/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */
/*     } */
/*     else { */
/*       III = (mx-2) * (my-2) * (mz-2); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat_sq[k+1][j+1][i+1].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat_sq[k+1][j+1][i+1].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat_sq[k+1][j+1][i+1].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*     } */

/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_square_sum, &ucat_sq); */

/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat_cross_sum, &ucat_cross); */

/*     // UV VW WU */
/*     for (k=0; k<mz-1; k++) { */
/*       for (j=0; j<my-1; j++) { */
/* 	for (i=0; i<mx-1; i++) { */
/* 	  ucat_o[k][j][i].x = 0.125 * */
/* 	    (ucat_cross[k][j][i].x + ucat_cross[k][j][i+1].x + */
/* 	     ucat_cross[k][j+1][i].x + ucat_cross[k][j+1][i+1].x + */
/* 	     ucat_cross[k+1][j][i].x + ucat_cross[k+1][j][i+1].x + */
/* 	     ucat_cross[k+1][j+1][i].x + ucat_cross[k+1][j+1][i+1].x); */
/* 	  ucat_o[k][j][i].y = 0.125 * */
/* 	    (ucat_cross[k][j][i].y + ucat_cross[k][j][i+1].y + */
/* 	     ucat_cross[k][j+1][i].y + ucat_cross[k][j+1][i+1].y + */
/* 	     ucat_cross[k+1][j][i].y + ucat_cross[k+1][j][i+1].y + */
/* 	     ucat_cross[k+1][j+1][i].y + ucat_cross[k+1][j+1][i+1].y); */
/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat_cross[k][j][i].z + ucat_cross[k][j][i+1].z + */
/* 	     ucat_cross[k][j+1][i].z + ucat_cross[k][j+1][i+1].z + */
/* 	     ucat_cross[k+1][j][i].z + ucat_cross[k+1][j][i+1].z + */
/* 	     ucat_cross[k+1][j+1][i].z + ucat_cross[k+1][j+1][i+1].z); */
/* 	} */
/*       } */
/*     } */

/*     if (!vc) { */
/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-1; k++) { */
/* 	for (j=ys; j<ye-1; j++) { */
/* 	  for (i=xs; i<xe-1; i++) { */
/* 	    x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */



/*     } */
/*     else { */



/*       III = (mx-2) * (my-2) * (mz-2); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat_cross[k+1][j+1][i+1].x; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat_cross[k+1][j+1][i+1].y; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*       for (k=zs; k<ze-2; k++) { */
/* 	for (j=ys; j<ye-2; j++) { */
/* 	  for (i=xs; i<xe-2; i++) { */
/* 	    x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat_cross[k+1][j+1][i+1].z; */
/* 	  } */
/* 	} */
/*       } */
/*       II = TECDAT100(&III, x, &DIsDouble); */

/*     } */


/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_cross_sum, &ucat_cross); */
  
/*     VecDestroy(&X); */


/*     III = (mx-2) * (my-2) * (mz-2); */
/*     DMCreateGlobalVector(da2, &X); */

/*     DMDAVecGetArray(user[bi].da, user[bi].P, &p); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1]; */

/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].P, &p); */

/*     DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1]; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */


/*     VecDestroy(&X); */

      

/*     PetscFree(x); */
/*     DMDestroy(&da1); */
/*     DMDestroy(&da2); */
/*   } */
/*   II = TECEND100(); */
/*   return 0; */
/* } */

/* PetscErrorCode ExtractIO(UserCtx *user) */
/* { */
/*   PetscInt bi; */
  
/*   FILE *f; */
/*   char filen[80]; */
/*   sprintf(filen, "ExtractIO_u.dat"); */

/*   f = fopen(filen, "a"); */

/*   PetscFPrintf(PETSC_COMM_WORLD, f, "%d ", ti); */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */

/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k; */

/*     PetscInt	i_min=1, j_min=120, k_min=110; */
/*     PetscInt	i_max=2, j_max=121, k_max=129; */


/*     PetscReal	***aj; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o; */
/*     PetscReal	***p, ***nvert; */

/*     INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0; */
/*     INTEGER4	IsBlock=1, NusnmFaceConnections=0, FaceNeighborMode=0; */
/*     INTEGER4    ShareConnectivityFromZone=0; */
/*     INTEGER4	LOC[8] = {1, 1, 1, 1, 1, 1, 0, 0}; /\* 1 is cell-centered */
/* 						      0 is node centered *\/ */

/*     PetscInt	vc = 0; */
/*     PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL); */
    
/*     if (vc) { */
/*       LOC[3]=0; LOC[4]=0; LOC[5]=0; */
/*     } */
    
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecGetArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
/*     for (k=k_min; k<k_max; k++) { */
/*       for (j=j_min; j<j_max; j++) { */
/* 	for (i=i_min; i<i_max; i++) { */

/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */

/* 	  PetscFPrintf(PETSC_COMM_WORLD, f, "%le ", ucat_o[k][j][i].z); */

/* 	} */
/*       } */
/*     } */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat); */
/*     DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat_o, &ucat_o); */
    
/*   } */

/*   fclose(f); */
/*   return 0;   */
/* } */


/* PetscErrorCode TECIOOutQ(UserCtx *user) */
/* { */
/*   PetscInt bi; */

/*   char filen[80]; */
/*   sprintf(filen, "QCriteria%3.3d.plt", ti); */

/*   INTEGER4 II, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax; */
/*   VIsDouble = 0; */
/*   DIsDouble = 0; */
/*   Debug = 0; */

/*   II = TECINI100("Heart Valve", */
/* 	     "X Y Z P l2 Nv", */
/* 	     filen, */
/* 	     ".", */
/* 	     &Debug, */
/* 	     &VIsDouble); */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */

/*     DM da1, da2; */

/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k; */
/*     PetscReal	***aj; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o; */
/*     PetscReal	***p, ***nvert, ***l2; */
/*     Vec		Coor, zCoor, nCoor; */
/*     VecScatter	ctx; */

/*     Vec X, Y, Z, U, V, W; */

/*     INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0; */
/*     INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0; */
/*     INTEGER4    ShareConnectivityFromZone=0; */
/*     INTEGER4	LOC[8] = {1, 1, 1, 0, 0}; /\* 1 is cell-centered */
/* 						      0 is node centered *\/ */


/*     float	*x; */
/*     PetscMalloc(mx*my*mz*sizeof(float), &x); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da1); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da2); */

/*     IMax = mx-1; JMax = my-1; KMax = mz-1; */

/*     II = TECZNE100("Block 1", */
/* 		  &ZoneType, 	/\* Ordered zone *\/ */
/* 		  &IMax, */
/* 		  &JMax, */
/* 		  &KMax, */
/* 		  &ICellMax, */
/* 		  &JCellMax, */
/* 		  &KCellMax, */
/* 		  &IsBlock,	/\* ISBLOCK  BLOCK format *\/ */
/* 		  &NumFaceConnections, */
/* 		  &FaceNeighborMode, */
/* 		  LOC, */
/* 		  NULL, */
/* 		  &ShareConnectivityFromZone); /\* No connectivity sharing *\/ */

/*     III = (mx-1) * (my-1) * (mz-1); */
/*     DMCreateGlobalVector(da1, &X); */

/*     DMGetCoordinates(da, &Coor); */
/*     DMDAVecGetArray(fda, Coor, &coor); */


/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(fda, Coor, &coor); */
    

/*     VecDestroy(&X); */


/*     III = (mx-2) * (my-2) * (mz-2); */
/*     DMCreateGlobalVector(da2, &X); */

/*     DMDAVecGetArray(user[bi].da, user[bi].P, &p); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1]; */

/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].P, &p); */

/*     DMDAVecGetArray(user[bi].da, user[bi].Phi, &l2); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = l2[k+1][j+1][i+1]; */

/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &l2); */

/*     DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1]; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */

/*     VecDestroy(&X); */


/*     PetscFree(x); */
/*     DMDestroy(&da1); */
/*     DMDestroy(&da2); */
/*   }  //bi */
/*   II = TECEND100(); */
/*   return 0; */
/* } */

/* PetscErrorCode TECIOOutS(UserCtx *user) */
/* { */
/*   PetscInt bi; */

/*   char filen[80]; */
/*   sprintf(filen, "MaxStress%3.3d.plt", ti); */

/*   INTEGER4 II, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax; */
/*   VIsDouble = 0; */
/*   DIsDouble = 0; */
/*   Debug = 0; */

/*   II = TECINI100("Heart Valve", */
/* 	     "X Y Z P Nv", */
/* 	     filen, */
/* 	     ".", */
/* 	     &Debug, */
/* 	     &VIsDouble); */

/*   for (bi=0; bi<block_number; bi++) { */
/*     DM da = user[bi].da, fda = user[bi].fda; */
/*     DMDALocalInfo info = user[bi].info; */

/*     DM da1, da2; */

/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*     PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
    
/*     PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*     PetscInt	i, j, k; */
/*     PetscReal	***aj; */
/*     Cmpnts	***ucat, ***coor, ***ucat_o; */
/*     PetscReal	***p, ***nvert, ***l2; */
/*     Vec		Coor, zCoor, nCoor; */
/*     VecScatter	ctx; */

/*     Vec X, Y, Z, U, V, W; */

/*     INTEGER4	ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0; */
/*     INTEGER4	IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0; */
/*     INTEGER4    ShareConnectivityFromZone=0; */
/*     INTEGER4	LOC[8] = {1, 1, 1, 0, 0}; /\* 1 is cell-centered */
/* 						      0 is node centered *\/ */


/*     float	*x; */
/*     PetscMalloc(mx*my*mz*sizeof(float), &x); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da1); */
/*     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, */
/* 	       mx-1, my-1, mz-1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
/* 	       1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da2); */

/*     IMax = mx-1; JMax = my-1; KMax = mz-1; */

/*     II = TECZNE100("Block 1", */
/* 		  &ZoneType, 	/\* Ordered zone *\/ */
/* 		  &IMax, */
/* 		  &JMax, */
/* 		  &KMax, */
/* 		  &ICellMax, */
/* 		  &JCellMax, */
/* 		  &KCellMax, */
/* 		  &IsBlock,	/\* ISBLOCK  BLOCK format *\/ */
/* 		  &NumFaceConnections, */
/* 		  &FaceNeighborMode, */
/* 		  LOC, */
/* 		  NULL, */
/* 		  &ShareConnectivityFromZone); /\* No connectivity sharing *\/ */

/*     III = (mx-1) * (my-1) * (mz-1); */
/*     DMCreateGlobalVector(da1, &X); */

/*     DMGetCoordinates(da, &Coor); */
/*     DMDAVecGetArray(fda, Coor, &coor); */


/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */

/*     for (k=zs; k<ze-1; k++) { */
/*       for (j=ys; j<ye-1; j++) { */
/* 	for (i=xs; i<xe-1; i++) { */
/* 	  x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(fda, Coor, &coor); */
    

/*     VecDestroy(&X); */


/*     III = (mx-2) * (my-2) * (mz-2); */
/*     DMCreateGlobalVector(da2, &X); */

/*     DMDAVecGetArray(user[bi].da, user[bi].P, &p); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1]; */

/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].P, &p); */

/*     DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert); */
/*     for (k=zs; k<ze-2; k++) { */
/*       for (j=ys; j<ye-2; j++) { */
/* 	for (i=xs; i<xe-2; i++) { */
/* 	  x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1]; */
/* 	} */
/*       } */
/*     } */
/*     II = TECDAT100(&III, x, &DIsDouble); */
/*     DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert); */

/*     VecDestroy(&X); */


/*     PetscFree(x); */
/*     DMDestroy(&da1); */
/*     DMDestroy(&da2); */
/*   } */
/*   II = TECEND100(); */
/*   return 0; */
/* } */

/* PetscErrorCode TecOut(UserCtx *user) */
/* { */
/*   DM		da = user->da, fda = user->fda; */
/*   DMDALocalInfo	info = user->info; */

/*   PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*   PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*   PetscInt	zs = info.zs, ze = info.zs + info.zm; */
/*   PetscInt	mx = info.mx, my = info.my, mz = info.mz; */

/*   PetscInt	lxs, lys, lzs, lxe, lye, lze; */
/*   PetscInt	i, j, k; */

/*   PetscReal	***aj; */
/*   Cmpnts	***ucat, ***coor, ***ucat_o; */
/*   PetscReal	***p, ***nvert; */
/*   Vec		Coor, zCoor, nCoor; */
/*   VecScatter	ctx; */

/*   FILE	*f; */
/*   char filen[80]; */

/*   PetscInt	tt=0; */


/*   DMGetCoordinates(da, &Coor); */

/*   PetscInt rank; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */

/*   if (!rank) { */
/*     DMDAVecGetArray(user->fda, Coor, &coor); */
/*     DMDAVecGetArray(user->fda, user->Ucat, &ucat); */
/*     DMDAVecGetArray(user->fda, user->Ucat_o, &ucat_o); */
/*     DMDAVecGetArray(user->da, user->Nvert, &nvert); */
/*     DMDAVecGetArray(user->da, user->P, &p); */

/*   for (k=1; k<mz-1; k++) { */
/*     for (j=0; j<my-1; j++) { */
/*       for (i=0; i<mx-1; i++) { */

/* 	  ucat_o[k][j][i].x = 0.125 * */
/* 	    (ucat[k][j][i].x + ucat[k][j][i+1].x + */
/* 	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x + */
/* 	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x + */
/* 	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x); */
/* 	  ucat_o[k][j][i].y = 0.125 * */
/* 	    (ucat[k][j][i].y + ucat[k][j][i+1].y + */
/* 	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y + */
/* 	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y + */
/* 	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y); */
/* 	  ucat_o[k][j][i].z = 0.125 * */
/* 	    (ucat[k][j][i].z + ucat[k][j][i+1].z + */
/* 	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z + */
/* 	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z + */
/* 	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z); */

/*       } */
/*     } */
/*   } */

/*   sprintf(filen, "Result%3.3d.dat", ti); */
/*   f = fopen(filen, "a"); */

  
/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     coor[k-1][j-1][i-1].x); */
/*       } */
/*     } */
/*   } */
/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     coor[k-1][j-1][i-1].y); */
/*       } */
/*     } */
/*   } */

/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     coor[k-1][j-1][i-1].z); */
/*       } */
/*     } */
/*   } */
/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     ucat_o[k-1][j-1][i-1].x); */
/*       } */
/*     } */
/*   } */
/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     ucat_o[k-1][j-1][i-1].y); */
/*       } */
/*     } */
/*   } */
/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n");  */
/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     ucat_o[k-1][j-1][i-1].z); */
/*       } */
/*     } */
/*   } */

/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*   for (k=1; k<mz-1; k++) { */
/*     for (j=1; j<my-1; j++) { */
/*       for (i=1; i<mx-1; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", */
/* 		     p[k][j][i]); */
/*       } */
/*     } */
/*   } */

/*   PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*   for (k=1; k<mz-1; k++) { */
/*     for (j=1; j<my-1; j++) { */
/*       for (i=1; i<mx-1; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", */
/* 		     (int)(nvert[k][j][i]+0.5)); */
/*       } */
/*     } */
/*   } */


/*   fclose(f); */
  
/*   DMDAVecRestoreArray(user->fda, Coor, &coor); */
/*   DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); */
/*   DMDAVecRestoreArray(user->fda, user->Ucat_o, &ucat_o); */
/*   DMDAVecRestoreArray(user->da, user->Nvert, &nvert); */
/*   DMDAVecRestoreArray(user->da, user->P, &p); */
  
/*   } */
  
/*   return(0); */
/* } */

PetscErrorCode FormMetrics(UserCtx *user)
{
  DM		cda;
  Cmpnts	***csi, ***eta, ***zet;
  PetscScalar	***aj;
  Vec		coords;
  Cmpnts	***coor;

  DM		da = user->da, fda = user->fda;
  Vec		Csi = user->Csi, Eta = user->Eta, Zet = user->Zet;
  Vec		Aj = user->Aj;
  Vec		ICsi = user->ICsi, IEta = user->IEta, IZet = user->IZet;
  Vec		JCsi = user->JCsi, JEta = user->JEta, JZet = user->JZet;
  Vec		KCsi = user->KCsi, KEta = user->KEta, KZet = user->KZet;
  Vec		IAj = user->IAj, JAj = user->JAj, KAj = user->KAj;

  
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  Cmpnts	***gs;
  PetscReal	***iaj, ***jaj, ***kaj;

  Vec		Cent = user->Cent; //local working array for storing cell center geometry

  Vec		Centx, Centy, Centz, lCoor;
  Cmpnts	***cent, ***centx, ***centy, ***centz;

  PetscInt	xs, ys, zs, xe, ye, ze;
  DMDALocalInfo	info;

  PetscInt	mx, my, mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
  PetscInt	i, j, k, ia, ja, ka, ib, jb, kb;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscErrorCode	ierr;

  PetscReal	xcp, ycp, zcp, xcm, ycm, zcm;
  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DMGetCoordinateDM(da, &cda);
  DMDAVecGetArray(cda, Csi, &csi);
  DMDAVecGetArray(cda, Eta, &eta);
  DMDAVecGetArray(cda, Zet, &zet);
  ierr = DMDAVecGetArray(da, Aj,  &aj); CHKERRQ(ierr);

  DMGetCoordinatesLocal(da, &coords);
  DMDAVecGetArray(fda, coords, &coor);


  //  VecDuplicate(coords, &Cent);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  /* Calculating transformation metrics in i direction */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=xs; i<lxe; i++){
	/* csi = de X dz */
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j-1][i  ].x - coor[k-1][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j-1][i  ].y - coor[k-1][j-1][i  ].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j-1][i  ].z - coor[k-1][j-1][i  ].z);
				       		    	    	 
	dxdz = 0.5 * (coor[k  ][j-1][i  ].x + coor[k  ][j  ][i  ].x -
		      coor[k-1][j-1][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k  ][j-1][i  ].y + coor[k  ][j  ][i  ].y -
		      coor[k-1][j-1][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k  ][j-1][i  ].z + coor[k  ][j  ][i  ].z -
		      coor[k-1][j-1][i  ].z - coor[k-1][j  ][i  ].z);
	  
	csi[k][j][i].x = dyde * dzdz - dzde * dydz;
	csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	csi[k][j][i].z = dxde * dydz - dyde * dxdz;

	
      }
    }
  }

  // Need more work -- lg65
  /* calculating j direction metrics */
  for (k=lzs; k<lze; k++){
    for (j=ys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	/* eta = dz X de */
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z);
			    		         	 		   	 
	dxdz = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x);
	dydz = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y);
	dzdz = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z);
	  
	eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

      }
    }
  }

  /* calculating k direction metrics */
  for (k=zs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z);
			    		    	     	     	 
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z);
	  
	zet[k][j][i].x = dydc * dzde - dzdc * dyde;
	zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	zet[k][j][i].z = dxdc * dyde - dydc * dxde;
      }
    }
  }

  /* calculating Jacobian of the transformation */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	if (i>0 && j>0 && k>0) {
	  dxdc = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
			 coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydc = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
			 coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdc = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
			 coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

	  dxde = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x - 
			 coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x);
	  dyde = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y - 
			 coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y);
	  dzde = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z - 
			 coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z);

	  dxdz = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydz = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdz = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i  ].z - coor[k-1][j-1][i  ].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);
	  
	  aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	    dydc * (dxde * dzdz - dzde * dxdz) +
	    dzdc * (dxde * dydz - dyde * dxdz);
	  aj[k][j][i] = 1./aj[k][j][i];
	}
      }
    }
  }

  // mirror grid outside the boundary
  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	eta[k][j][i].x = eta[k][j][i+1].x;
	eta[k][j][i].y = eta[k][j][i+1].y;
	eta[k][j][i].z = eta[k][j][i+1].z;
	  
	zet[k][j][i].x = zet[k][j][i+1].x;
	zet[k][j][i].y = zet[k][j][i+1].y;
	zet[k][j][i].z = zet[k][j][i+1].z;

	aj[k][j][i] = aj[k][j][i+1];
      }
    }
  }

  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	eta[k][j][i].x = eta[k][j][i-1].x;
	eta[k][j][i].y = eta[k][j][i-1].y;
	eta[k][j][i].z = eta[k][j][i-1].z;
	  
	zet[k][j][i].x = zet[k][j][i-1].x;
	zet[k][j][i].y = zet[k][j][i-1].y;
	zet[k][j][i].z = zet[k][j][i-1].z;

	aj[k][j][i] = aj[k][j][i-1];
      }
    }
  }

  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	csi[k][j][i].x = csi[k][j+1][i].x;
	csi[k][j][i].y = csi[k][j+1][i].y;
	csi[k][j][i].z = csi[k][j+1][i].z;
	  
	zet[k][j][i].x = zet[k][j+1][i].x;
	zet[k][j][i].y = zet[k][j+1][i].y;
	zet[k][j][i].z = zet[k][j+1][i].z;

	aj[k][j][i] = aj[k][j+1][i];
      }
    }
  }

  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	csi[k][j][i].x = csi[k][j-1][i].x;
	csi[k][j][i].y = csi[k][j-1][i].y;
	csi[k][j][i].z = csi[k][j-1][i].z;
	  
	zet[k][j][i].x = zet[k][j-1][i].x;
	zet[k][j][i].y = zet[k][j-1][i].y;
	zet[k][j][i].z = zet[k][j-1][i].z;

	aj[k][j][i] = aj[k][j-1][i];
      }
    }
  }

  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	eta[k][j][i].x = eta[k+1][j][i].x;
	eta[k][j][i].y = eta[k+1][j][i].y;
	eta[k][j][i].z = eta[k+1][j][i].z;
	  		      			   
	csi[k][j][i].x = csi[k+1][j][i].x;
	csi[k][j][i].y = csi[k+1][j][i].y;
	csi[k][j][i].z = csi[k+1][j][i].z;

	aj[k][j][i] = aj[k+1][j][i];
      }
    }
  }

  if (ze==mz){
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	eta[k][j][i].x = eta[k-1][j][i].x;
	eta[k][j][i].y = eta[k-1][j][i].y;
	eta[k][j][i].z = eta[k-1][j][i].z;
	  		      			   
	csi[k][j][i].x = csi[k-1][j][i].x;
	csi[k][j][i].y = csi[k-1][j][i].y;
	csi[k][j][i].z = csi[k-1][j][i].z;

	aj[k][j][i] = aj[k-1][j][i];
      }
    }
  }

  DMDAVecGetArray(cda, user->Cent, &cent);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	cent[k][j][i].x = 0.125 *
	  (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
	   coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
	   coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
	   coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
	cent[k][j][i].y = 0.125 *
	  (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
	   coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
	   coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
	   coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
	cent[k][j][i].z = 0.125 *
	  (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
	   coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
	   coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
	   coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);
      }
    }
  }
  DMDAVecRestoreArray(cda, user->Cent, &cent);
 



  DMDAVecRestoreArray(cda, Csi, &csi);
  DMDAVecRestoreArray(cda, Eta, &eta);
  DMDAVecRestoreArray(cda, Zet, &zet);
  DMDAVecRestoreArray(da, Aj,  &aj);


  DMDAVecRestoreArray(cda, coords, &coor);


  VecAssemblyBegin(Csi);
  VecAssemblyEnd(Csi);
  VecAssemblyBegin(Eta);
  VecAssemblyEnd(Eta);
  VecAssemblyBegin(Zet);
  VecAssemblyEnd(Zet);
  VecAssemblyBegin(Aj);
  VecAssemblyEnd(Aj);

  PetscBarrier(PETSC_NULL);
  return 0;
}

PetscErrorCode Ucont_P_Binary_Input_Ave(UserCtx *user, PetscInt Ave)
{
  PetscViewer	viewer;
  char filen[90];
  PetscInt tie, mm;
  PetscBool flag,flag2;
  Vec Usub,Psub;
  PetscReal scale, min,phscale;
  int phave=0;

  
  PetscOptionsGetInt(PETSC_NULL, "-tim", &tie, &flag);
  PetscOptionsGetReal(PETSC_NULL, "-scale", &scale, &flag2);
  PetscOptionsGetInt(PETSC_NULL, "-unorm", &unorm, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL); 
  PetscOptionsGetInt(PETSC_NULL, "-phave", &phave, PETSC_NULL); 
 
 if (flag) 
    VecDuplicate(user->Ucat, &Usub);
    VecDuplicate(user->P, &Psub);
if (!phave){
  //  if (Ave==1) {
    sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->Ucat,viewer);
    PetscViewerDestroy(&viewer);

    if (flag) {
    sprintf(filen, "su0_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat, -1.,Usub);
    VecSet(Usub, 0.);
    }
    if (flag2) VecScale(user->Ucat,scale);
    VecMin(user->Ucat,&mm, &min) ;
    PetscPrintf(PETSC_COMM_WORLD,"Min averaging UVW... %le \n", min);

    //  } else if (Ave==2) {
    sprintf(filen, "su1_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->Ucat_cross_sum,viewer);
    PetscViewerDestroy(&viewer);

    if (flag) {
    sprintf(filen, "su1_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat_cross_sum, -1.,Usub);
    VecSet(Usub, 0.);
    }
    if (flag2) VecScale(user->Ucat_cross_sum,scale);
    VecMin(user->Ucat_cross_sum,&mm, &min) ;
    PetscPrintf(PETSC_COMM_WORLD,"Min averaging UVW cross... %le \n", min);
    //  } else if (Ave==3) {
    sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad((user->Ucat_square_sum),viewer);
    PetscViewerDestroy(&viewer);
    
    sprintf(filen, "sustar_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad((user->lUstar),viewer);
    PetscViewerDestroy(&viewer);

    if (flag) {
    sprintf(filen, "su2_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat_square_sum, -1.,Usub);
    VecSet(Usub, 0.);
    }
    //  }
    if (flag2) VecScale(user->Ucat_square_sum,scale);
    VecMin(user->Ucat_square_sum,&mm, &min) ;
    PetscPrintf(PETSC_COMM_WORLD,"Min averaging UVW... %le \n", min);

    if (flag) VecDestroy(&Usub);
  sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad(user->P,viewer);
  PetscViewerDestroy(&viewer);

   if (flag) {
    sprintf(filen, "sp_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Psub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->P, -1.,Psub);
    VecSet(Psub, 0.);
    }
 
  if (flag2) VecScale(user->P,scale);

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert),viewer);
  PetscViewerDestroy(&viewer);
  if (fish) {
   /* {sprintf(filen, "NVAVE.dat"); */
   /* PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer); */
   /* VecLoad((user->Nvert1),viewer); */
   /* PetscViewerDestroy(&viewer); */
   }

  PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... %d ***\n\n", filen, Ave);


 }//!phave
 if (phave) { 
   PetscInt tis, tsteps=5,phase=0;
   scale=1.0;flag=0;

   PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, PETSC_NULL);
   PetscOptionsGetInt(PETSC_NULL, "-phase", &phase, PETSC_NULL);
   PetscOptionsGetInt(PETSC_NULL, "-tim", &tie, &flag); 
   PetscOptionsGetInt(PETSC_NULL, "-tsph", &tsteps, PETSC_NULL);
   PetscOptionsGetReal(PETSC_NULL, "-scale", &scale, PETSC_NULL);
   VecDuplicate(user->Ucat, &Usub);
   VecDuplicate(user->P, &Psub);
  
if (phase==0)
   {
    sprintf(filen, "suph0_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->Ucat,viewer);
    PetscViewerDestroy(&viewer);

   if (flag){ 
    sprintf(filen, "suph0_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat, -1.,Usub);
    VecSet(Usub, 0.);
    }
    VecScale(user->Ucat,scale);
   
  sprintf(filen, "spph0_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad(user->P,viewer);
  PetscViewerDestroy(&viewer);

  if (flag){ 
    sprintf(filen, "spph0_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Psub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->P, -1.,Psub);
    VecSet(Psub, 0.);
    }
    VecScale(user->P,scale);
  sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert),viewer);
  PetscViewerDestroy(&viewer);
 
   }
if (phase==1)
   {
    sprintf(filen, "suph1_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->Ucat,viewer);
    PetscViewerDestroy(&viewer);

  if (flag){  
    sprintf(filen, "suph1_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat, -1.,Usub);
    VecSet(Usub, 0.);
    }
    VecScale(user->Ucat,scale);
  
  sprintf(filen, "spph1_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad(user->P,viewer);
  PetscViewerDestroy(&viewer);

   if (flag){ 
    sprintf(filen, "spph1_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Psub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->P, -1.,Psub);
    VecSet(Psub, 0.);
    }
    VecScale(user->P,scale);
 
   }

if (phase==2)
   {
    sprintf(filen, "suph2_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->Ucat,viewer);
    PetscViewerDestroy(&viewer);

     if (flag){
    sprintf(filen, "suph2_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat, -1.,Usub);
    VecSet(Usub, 0.);
    }
    VecScale(user->Ucat,scale);
    
  sprintf(filen, "spph2_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad(user->P,viewer);
  PetscViewerDestroy(&viewer);

    if (flag){
    sprintf(filen, "spph2_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Psub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->P, -1.,Psub);
    VecSet(Psub, 0.);
    }
    VecScale(user->P,scale);
 
   }


if (phase==3)
   {
    sprintf(filen, "suph3_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->Ucat,viewer);
    PetscViewerDestroy(&viewer);

    
    sprintf(filen, "suph3_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Usub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->Ucat, -1.,Usub);
    VecSet(Usub, 0.);
    
    VecScale(user->Ucat,scale);
   
  sprintf(filen, "spph3_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad(user->P,viewer);
  PetscViewerDestroy(&viewer);

   
    sprintf(filen, "spph3_%06d_%1d.dat", tie, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(Psub,viewer);
    PetscViewerDestroy(&viewer);
    VecAXPY(user->P, -1.,Psub);
    VecSet(Psub, 0.);
    
    VecScale(user->P,scale);
 
   }




}




/* 
   sprintf(filen, "ufield%5.5d_%1.1d.dat", tis, user->_this);
   PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
   VecLoad(user->Ucat,viewer);
   PetscViewerDestroy(&viewer);
   PetscPrintf(PETSC_COMM_WORLD,"first ufield...  \n" );


   int kk=1;
   for (ti=tis;ti<=tie;ti+=tsteps){ 
     kk++;
     sprintf(filen, "ufield%5.5d_%1.1d.dat", ti, user->_this);
     PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
     VecLoad(Usub,viewer);
     PetscViewerDestroy(&viewer);
     VecAXPY(user->Ucat, 1.,Usub);
     VecSet(Usub, 0.);
     PetscPrintf(PETSC_COMM_WORLD,"total ufield  \n" );


   }

    phscale =1./kk;

    VecScale(user->Ucat,phscale);

    VecDestroy(&Usub);
    sprintf(filen, "pfield%5.5d_%1.1d.dat", tis, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoad(user->P,viewer);
    PetscViewerDestroy(&viewer);

    for (ti=tis+tsteps;ti<=tie;ti+=tsteps){ 
      sprintf(filen, "pfield%5.5d_%1.1d.dat", tis, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(Psub,viewer);
      PetscViewerDestroy(&viewer);
      VecAXPY(user->P, 1.,Psub);
      VecSet(Psub, 0.);
       PetscPrintf(PETSC_COMM_WORLD,"pfield  \n" );


    }
 
   VecScale(user->P,phscale);
  sprintf(filen, "nvfield%5.5d_%1.1d.dat", tis, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert),viewer);
  PetscViewerDestroy(&viewer);
 
   ti=tie;

 }
*/

/*
   sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti-600, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert2),viewer);
  PetscViewerDestroy(&viewer);
 
   sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti-900, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert3),viewer);
  PetscViewerDestroy(&viewer);
 

*/

//amir 
//
//
//
/*




  VecAXPY(user->Nvert,1.0,user->Nvert1);
 
  VecAXPY(user->Nvert,1.0,user->Nvert2);
 
  VecAXPY(user->Nvert,1.0,user->Nvert3);
 
  sprintf(filen, "NVAVE.dat");

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);
 */


 
    return(0);
}

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{

  PetscViewer	viewer;



  char filen2[90];

  PetscOptionsClearValue("-vecload_block_size");
  sprintf(filen2, "results/pfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewer	pviewer;
  Vec temp;

  PetscInt rank;
  DMDACreateNaturalVector(user->da, &temp);


  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);


  VecLoad((user->P),pviewer);
  PetscReal norm;
  VecNorm(user->P, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);
  PetscViewerDestroy(&pviewer);

  VecDestroy(&temp);

  sprintf(filen2, "results/nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);

  VecLoad((user->Nvert),pviewer);
  PetscViewerDestroy(&pviewer);
  


}

PetscErrorCode Ucont_P_Binary_Input1(UserCtx *user, PetscInt flg)
{
  PetscViewer viewer;
  char filen[90];
  PetscInt bi=user->_this;
  
  sprintf(filen, "results/ufield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad((user->Ucat),viewer);
 
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);

}


/* ==================================================================================             */
PetscErrorCode rotate_force_data(FSInfo *fsi, PetscInt ibi, PetscInt bi)
{  
  Cmpnts    F, rot, a_c;
  PetscReal lift1, lift2, drag1, drag2, radial2, lift3, drag3;

  rot.x=fsi->S_ang_n[0];//-FSinfo->S_ang_o[0];
  rot.y=fsi->S_ang_n[2];//-FSinfo->S_ang_o[0];
  rot.z=fsi->S_ang_n[4];//-FSinfo->S_ang_o[0];

  F.x = fsi->F_x;
  F.y = fsi->F_y;
  F.z = fsi->F_z;

  lift1 =  F.z*cos(rot.y)-F.x*sin(rot.y);
  drag1 =  F.z*sin(rot.y)+F.x*cos(rot.y);

  lift3 =  F.z*cos(rot.y);
  drag3 =  F.z*sin(rot.y);

  a_c.x=a_c.y=a_c.z=0.;

  PetscPrintf(PETSC_COMM_WORLD, "Force Rotate: %le %le %le \n",rot.x,rot.y,rot.z);

  F = rotate_xyz(F, a_c, rot);
   
  lift2 =  F.z;
  drag2 =  F.x*cos(-rot.z)-F.y*sin(-rot.z);
  radial2= F.x*sin(-rot.z)+F.y*cos(-rot.z);

  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);                                                                                                                                                                                                              
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Inertial_data_%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le %le\n",ti, F.x,F.y,F.z, lift1, lift2,lift3, drag1, drag2, drag3, fsi->F_y, radial2);
    fclose(f);
  }

  return(0);
}


/* PetscErrorCode ibm_intp_pj(UserCtx *user, IBMNodes *ibm, FSInfo *fsi) */
/* { */
  
/*   DM	        da = user->da, fda = user->fda; */
/*   DMDALocalInfo	info = user->info; */
/*   PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*   PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*   PetscInt	zs = info.zs, ze = info.zs + info.zm; */
  
  
/*   PetscInt        i, j, k; */
/*   PetscInt	ncx = 20, ncy = 20, ncz = 20, bag=400; */
/*   List          *cell_trg; */
/*   node          *current; */
  
/*   PetscInt	n_elmt = ibm->n_elmt, n_v = ibm->n_v; */
/*   PetscInt	*nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3= ibm->nv3; */
/*   PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max; */
  
  
/*   PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp; */
/*   PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z; */
/*   PetscReal       *shear,*shear_x,*shear_y,*shear_z; */
/*   Vec             Coor; */
/*   Cmpnts          ***coor, ***ucat; */
/*   PetscReal       ***prs, ***nvert; */
/*   PetscReal       dcx, dcy, dcz; */
/*   PetscReal       dpj=0.1; */
/*   PetscReal       length_scale,specific_weight=1., Re,Vel=1.; */
/*   PetscInt        n1e, n2e, n3e; */
/*   PetscReal       xc, yc, zc; */
/*   Cmpnts          pj; */
/*   PetscInt        ivs, jvs, kvs, ive, jve, kve, nv; */
/*   PetscReal       *vx, *vy,*vz,*pressure; */
/*   PetscInt        rank; */
/*   FILE            *f; */
/*   char            filen[80]; */
/*   PetscBool      flag=PETSC_FALSE, flg; */
/*   PetscInt        k_Export_Min=zs+1, k_Export_Max=ze-1; */
  

/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
  
/*   PetscPrintf(PETSC_COMM_WORLD, "Intp dpj %d %d %d\n", rank,n_v,n_elmt); */


/*   DMDAVecGetArray(fda, user->Cent, &coor); */
/*   DMDAVecGetArray(fda, user->Ucat, &ucat); */
/*   DMDAVecGetArray(da, user->P, &prs); */
/*   PetscReal p_inf=prs[2][2][2]; */
/*   DMDAVecRestoreArray(da, user->P, &prs); */
/*   VecShift(user->P, p_inf); */
  
/*   DMDAVecGetArray(da, user->P, &prs); */
  
/*   DMDAVecGetArray(da, user->Nvert, &nvert); */

/*   //Allocate the memory for shear component */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&shear); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&shear_x); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&shear_y); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&shear_z); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&vx); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&vy); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&vz); */
/*   PetscMalloc(n_elmt*sizeof(PetscReal),&pressure); */
  
/*   // Import the k values */
  
/*   PetscOptionsGetInt(PETSC_NULL, "-k_min", &k_Export_Min, &flg); */
/*   PetscOptionsGetInt(PETSC_NULL, "-k_max", &k_Export_Max, &flg); */
/*   PetscOptionsGetReal(PETSC_NULL, "-dpj", &dpj, &flg); */
  
/*   PetscPrintf(PETSC_COMM_WORLD," ----\n ti %d ibm_intp_pj ncx, ncy, ncz %le %d %d %d-- \n",ti, dpj, ncx, ncy, ncz); */
  
/*   //Start finding the thing */
/*   xbp_min = 1.e20; xbp_max = -1.e20; */
/*   ybp_min = 1.e20; ybp_max = -1.e20; */
/*   zbp_min = 1.e20; zbp_max = -1.e20; */
  
/*   for (k=k_Export_Min; k<k_Export_Max; k++) { */
/*     for (j=ys+1; j<ye-1; j++) { */
/*       for (i=xs+1; i<xe-1; i++) { */
/* 	xbp_min = PetscMin(xbp_min, coor[k][j][i].x); */
/* 	xbp_max = PetscMax(xbp_max, coor[k][j][i].x); */

/* 	ybp_min = PetscMin(ybp_min, coor[k][j][i].y); */
/* 	ybp_max = PetscMax(ybp_max, coor[k][j][i].y); */

/* 	zbp_min = PetscMin(zbp_min, coor[k][j][i].z); */
/* 	zbp_max = PetscMax(zbp_max, coor[k][j][i].z); */
/*       } */
/*     } */
/*   } */

/*   //Adjust  */
/*   xbp_min -= 0.01; xbp_max += 0.01; */
/*   ybp_min -= 0.01; ybp_max += 0.01; */
/*   zbp_min -= 0.01; zbp_max += 0.01; */

/*   //Devide the block in to nc small boxes */
/*     dcx = (xbp_max - xbp_min) / (ncx - 1.); */
/*     dcy = (ybp_max - ybp_min) / (ncy - 1.); */
/*     dcz = (zbp_max - zbp_min) / (ncz - 1.); */

/*   PetscMalloc(ncz * ncy * ncx * sizeof(List), &cell_trg); */
/*   //PetscPrintf(PETSC_COMM_SELF, "test00\n"); */
/*   for (k=0; k<ncz; k++) { */
/*     for (j=0; j<ncy; j++) { */
/*       for (i=0; i<ncx; i++) { */
/* 	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]); */
/*       } */
/*     } */
/*   } */

/*   //Identify the box - control volume  that point is now in */
/*   PetscInt iv, jv, kv, ln; */
/*   for (k=k_Export_Min; k<k_Export_Max; k++) { */
/*     for (j=ys+1; j<ye-1; j++) { */
/*       for (i=xs+1; i<xe-1; i++) { */
/* 	iv = floor((coor[k][j][i].x - xbp_min) / dcx); */
/* 	jv = floor((coor[k][j][i].y - ybp_min) / dcy); */
/* 	kv = floor((coor[k][j][i].z - zbp_min) / dcz); */
/* 	//We can select only the IB and fluid nodes */
/* 	// if (nvert[k][j][i] <= 1) */
	  
/* 	ln = k*xe*ye + j*xe + i; */
/* 	insertnode(&(cell_trg[kv *ncx*ncy + jv*ncx +iv]), ln); */
	  
/*       } */
/*     } */
/*   } */

  
/*   for (nv=0; nv<n_elmt; nv++)  */
/*     {  */
/*       vx[nv] = 0;  */
/*       vy[nv] = 0;  */
/*       vz[nv] = 0; */
      
/*       pressure[nv] = 0; */
      
/*       //3 vertice of the triangle */
/*       n1e = nv1[nv];  */
/*       n2e = nv2[nv];  */
/*       n3e = nv3[nv]; */
      
/*       // C is the center point of the triangle */
/*       xc = (x_bp[n1e] + x_bp[n2e] + x_bp[n3e]) / 3.; */
/*       yc = (y_bp[n1e] + y_bp[n2e] + y_bp[n3e]) / 3.; */
/*       zc = (z_bp[n1e] + z_bp[n2e] + z_bp[n3e]) / 3.; */
      
/*       // recontruct a line from C normal to the cell */
/*       // For the aneurysm case we use - normal surface because the normal vector originally outward */
/*       // We want it inward */
/*       pj.x = xc + dpj * nf_x[nv]; */
/*       pj.y = yc + dpj * nf_y[nv]; */
/*       pj.z = zc + dpj * nf_z[nv]; */
      
/*       // pj is one point inside the fluid domain */
      
/*       // locate point P into control volume number v */
/*       iv = floor((pj.x - xbp_min) / dcx); */
/*       jv = floor((pj.y - ybp_min) / dcy); */
/*       kv = floor((pj.z - zbp_min) / dcz); */
      
/*       // Total 27 cubic around the current iv,jv, kv */
/*       ivs = PetscMax(iv-1, 0);  */
/*       jvs = PetscMax(jv-1, 0);  */
/*       kvs = PetscMax(kv-1, 0); */
      
/*       ive = PetscMin(iv+2, ncx);  */
/*       jve = PetscMin(jv+2, ncy);  */
/*       kve = PetscMin(kv+2, ncz); */

/*       // Now interpolate the velocity from box s- e */
/*  host: for (k=kvs; k<kve; k++) { */
/* 	for (j=jvs; j<jve; j++) { */
/* 	  for (i= ivs; i<ive; i++) { */
	    
/* 	    PetscInt        ivh, jvh, kvh, ic, ln_vh; */
/* 	    PetscReal       d[6]; */
/* 	    Cmpnts          p[8]; */
	    
/* 	    //If the point (i,j,k) has one cell */
/* 	    current = cell_trg[kv*ncx*ncy+jv*ncx+iv].head; */
/* 	    while (current) { */
/* 	      ln_vh = current->Node; */
      
/* 	      kvh=ln_vh/(xe*ye); */
/* 	      jvh=(ln_vh-kvh*xe*ye)/xe; */
/* 	      ivh=ln_vh-kvh*xe*ye-jvh*xe; */
	      
/* 	      if (ivh < xe-2 && jvh < ye-2 && kvh < ze-2) { */
	      
/* 	      //Set 8 points around this point */
/* 	      // p[i] is one point in the bucket */
/* 	      p[0] = coor[kvh][jvh][ivh]; */
/* 	      p[1] = coor[kvh][jvh][ivh+1]; */
/* 	      p[2] = coor[kvh][jvh+1][ivh+1]; */
/* 	      p[3] = coor[kvh][jvh+1][ivh]; */
	      
/* 	      p[4] = coor[kvh+1][jvh][ivh]; */
/* 	      p[5] = coor[kvh+1][jvh][ivh+1]; */
/* 	      p[6] = coor[kvh+1][jvh+1][ivh+1]; */
/* 	      p[7] = coor[kvh+1][jvh+1][ivh]; */
	      
/* 	      //if point pj - projection  */
	      
/* 	      if (ISInsideCell(pj, p, d)) { */
/* 		// do the interpolation and then break */
		
/* 		PetscReal x, y, z; */
/* 		x = d[0] / (d[0]+d[1]); */
/* 		y = d[2] / (d[2]+d[3]); */
/* 		z = d[4] / (d[4]+d[5]); */
		
/* 		// Tetra - quaradtic interpolation */
/* 		vx[nv] = (ucat[kvh  ][jvh  ][ivh  ].x * (1-x) * (1-y) * (1-z) + */
/* 			  ucat[kvh  ][jvh  ][ivh+1].x * x * (1-y) * (1-z) + */
/* 			  ucat[kvh  ][jvh+1][ivh  ].x * (1-x) * y * (1-z) + */
/* 			  ucat[kvh+1][jvh  ][ivh  ].x * (1-x) * (1-y) * z + */
/* 			  ucat[kvh+1][jvh  ][ivh+1].x * x * (1-y) * z + */
/* 			  ucat[kvh+1][jvh+1][ivh  ].x * (1-x) * y * z + */
/* 			  ucat[kvh  ][jvh+1][ivh+1].x * x * y * (1-z) + */
/* 			  ucat[kvh+1][jvh+1][ivh+1].x * x * y * z); */
		
/* 		vy[nv] = (ucat[kvh  ][jvh  ][ivh  ].y * (1-x) * (1-y) * (1-z) + */
/* 			  ucat[kvh  ][jvh  ][ivh+1].y * x * (1-y) * (1-z) + */
/* 			  ucat[kvh  ][jvh+1][ivh  ].y * (1-x) * y * (1-z) + */
/* 			  ucat[kvh+1][jvh  ][ivh  ].y * (1-x) * (1-y) * z + */
/* 			  ucat[kvh+1][jvh  ][ivh+1].y * x * (1-y) * z + */
/* 			  ucat[kvh+1][jvh+1][ivh  ].y * (1-x) * y * z + */
/* 			  ucat[kvh  ][jvh+1][ivh+1].y * x * y * (1-z) + */
/* 			  ucat[kvh+1][jvh+1][ivh+1].y * x * y * z); */
		
/* 		vz[nv] = (ucat[kvh  ][jvh  ][ivh  ].z * (1-x) * (1-y) * (1-z) + */
/* 			  ucat[kvh  ][jvh  ][ivh+1].z * x * (1-y) * (1-z) + */
/* 			  ucat[kvh  ][jvh+1][ivh  ].z * (1-x) * y * (1-z) + */
/* 			  ucat[kvh+1][jvh  ][ivh  ].z * (1-x) * (1-y) * z + */
/* 			  ucat[kvh+1][jvh  ][ivh+1].z * x * (1-y) * z + */
/* 			  ucat[kvh+1][jvh+1][ivh  ].z * (1-x) * y * z + */
/* 			  ucat[kvh  ][jvh+1][ivh+1].z * x * y * (1-z) + */
/* 			  ucat[kvh+1][jvh+1][ivh+1].z * x * y * z); */
		
		
/* 		pressure[nv] = (prs[kvh  ][jvh  ][ivh  ] * (1-x) * (1-y) * (1-z) + */
/* 				prs[kvh  ][jvh  ][ivh+1] * x * (1-y) * (1-z) + */
/* 				prs[kvh  ][jvh+1][ivh  ] * (1-x) * y * (1-z) + */
/* 				prs[kvh+1][jvh  ][ivh  ] * (1-x) * (1-y) * z + */
/* 				prs[kvh+1][jvh  ][ivh+1] * x * (1-y) * z + */
/* 				prs[kvh+1][jvh+1][ivh  ] * (1-x) * y * z + */
/* 				prs[kvh  ][jvh+1][ivh+1] * x * y * (1-z) + */
/* 				prs[kvh+1][jvh+1][ivh+1] * x * y * z); */
		
		
/* 		flag = PETSC_TRUE; */
		
/* 		goto next; */
/* 	      } //if inside		 */
/* 	      } // if ivh< */
/* 	    current = current->next;  */
	   
/* 	  } // while current */
/* 	} */
/*       }	 */
/*     } */
/*  next: vx[nv] = vx[nv]; */
 

/*  if (flag == PETSC_FALSE) */
/*    PetscPrintf(PETSC_COMM_WORLD," ----\n  Error in interpolation -- \n"); */
/*  } */



/* // Now write the interpolation results to a file-------------------------------------------------------------- */


/*   if (!rank) { */
        
/*       PetscBool dyn=PETSC_FALSE,Hg=PETSC_FALSE,OSI=PETSC_FALSE; */

/*       PetscOptionsGetReal(PETSC_NULL, "-real_chact_leng", &length_scale, PETSC_NULL); */
/*       PetscOptionsGetReal(PETSC_NULL, "-specific_weight", &specific_weight, PETSC_NULL); */
/*       PetscOptionsGetReal(PETSC_NULL, "-ren", &Re, PETSC_NULL); */
/*       PetscOptionsGetReal(PETSC_NULL, "-real_chact_vel", &Vel, PETSC_NULL); */
/*       PetscOptionsGetBool(PETSC_NULL,"-dynes",&dyn,PETSC_NULL); */
/*       PetscOptionsGetBool(PETSC_NULL,"-Hg",&Hg,PETSC_NULL); */
/*       PetscOptionsGetBool(PETSC_NULL,"-osi",&OSI,PETSC_NULL); */

/*       sprintf(filen, "stress_pj%3.3d.dat",ti); */
/*       f = fopen(filen, "w"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z, U, V, W, Shear, Pressure, Shear_X, Shear_Y, Shear_Z\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE VARLOCATION=([4-11]=CELLCENTERED)\n", n_v, n_elmt); */
/*       // x - component */
/*       for (i=0; i<n_v; i++) { */

/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e \n", ibm->x_bp[i]); */
/*       } */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
/*       // y - component */
/*       for (i=0; i<n_v; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e \n", ibm->y_bp[i]); */
/*       } */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
/*       // z - component */
/*       for (i=0; i<n_v; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e \n", ibm->z_bp[i]); */
/*       } */

/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*       // V - component */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", vx[i]); */
/*       } */
    
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
     
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", vy[i]); */
/*       } */

/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", vz[i]); */
/*       } */
    
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*       PetscReal temp,vn,v_x,v_y,v_z,shear_mag, u_x,u_y,u_z; */
/*       PetscReal un, un_old, u_xold,u_yold,u_zold; */

/*       //Write out the shear stress */
/*       for (i=0; i<n_elmt; i++)  */
/*       { */
/* 	u_x=(ibm->u[ibm->nv1[i]].x+ibm->u[ibm->nv2[i]].x+ibm->u[ibm->nv3[i]].x)/3.; */
/* 	u_y=(ibm->u[ibm->nv1[i]].y+ibm->u[ibm->nv2[i]].y+ibm->u[ibm->nv3[i]].y)/3.; */
/* 	u_z=(ibm->u[ibm->nv1[i]].z+ibm->u[ibm->nv2[i]].z+ibm->u[ibm->nv3[i]].z)/3.; */

/* 	u_xold=(ibm->uold[ibm->nv1[i]].x+ibm->uold[ibm->nv2[i]].x+ibm->uold[ibm->nv3[i]].x)/3.; */
/* 	u_yold=(ibm->uold[ibm->nv1[i]].y+ibm->uold[ibm->nv2[i]].y+ibm->uold[ibm->nv3[i]].y)/3.; */
/* 	u_zold=(ibm->uold[ibm->nv1[i]].z+ibm->uold[ibm->nv2[i]].z+ibm->uold[ibm->nv3[i]].z)/3.; */

/* 	//Normal velocity magnitude */
/* 	vn = abs((vx[i]-u_x)*nf_x[i] + (vy[i]-u_y)*nf_y[i]+(vz[i]-u_z)*nf_z[i]); */

/* 	un = u_x * nf_x[i] + u_y* nf_y[i] + u_z* nf_z[i]; */
/* 	un_old = u_xold * nf_x[i] + u_yold * nf_y[i] + u_zold * nf_z[i]; */

/* 	//Tangential velocity */
/* 	v_x = vx[i] - u_x - nf_x[i]*vn; */
/* 	v_y = vy[i] - u_y - nf_y[i]*vn; */
/* 	v_z = vz[i] - u_z - nf_z[i]*vn; */

/* 	//Velocity magnitude */
/* 	shear_mag = sqrt(v_x*v_x + v_y*v_y + v_z*v_z); */
/* 	shear_mag = shear_mag/dpj * specific_weight / Re * Vel * Vel; */

/* 	//Convert to dynes/cm2 */
/* 	if (dyn) */
/* 	  shear_mag= shear_mag*10; */

/* 	//Shear of the element */
/* 	shear[i] = shear_mag; */
  
/* 	//Velocity magnitude */
/* 	temp=sqrt(v_x*v_x + v_y*v_y +v_z*v_z); */

/* 	if (temp > 0) */
/* 	  { */
/* 	    //Direction of shear */
/* 	    shear_x[i] = v_x*shear_mag/temp; */
/* 	    shear_y[i] = v_y*shear_mag/temp; */
/* 	    shear_z[i] = v_z*shear_mag/temp; */
/* 	  } */
/* 	else */
/* 	  { */
/* 	    shear_x[i]=0; */
/* 	    shear_y[i]=0; */
/* 	    shear_z[i]=0; */
/* 	  } */

/* 	pressure[nv] += (un - un_old) / user->dt; */

/* 	//1000 is specific weight of blood = appox water */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", shear[i]); */

/*       } */

/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */


/*       //Write out the Pressure */
/*        for (i=0; i<n_elmt; i++)  */
/*          { */
/* 	   //1000 is specific weight of blood = appox water */
/* 	   temp = pressure[i] * specific_weight * Vel * Vel; */
/* 	  if (Hg) */
/* 	    temp =temp*760/101325; */

/* 	    PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", temp); */

/*           } */

/*       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */


/*        //Write out the shear of the surface - x direction */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", shear_x[i]); */
/*       } */

/*        PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*        //Write out the shear of the surface - y direction */
/*       for (i=0; i<n_elmt; i++)  */
/*       { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", shear_y[i]); */
/*       } */

/*        PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */

/*        //Write out the shear of the surface - z direction */
/*       for (i=0; i<n_elmt; i++)  */
/*       { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", shear_z[i]); */
/*       } */

/*        PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); */
/*       //Write out the link nodes */

/*       for (i=0; i<n_elmt; i++)  */
/*            { */

/*                	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*             }    */


/*       fclose(f); */

/*       //Write the OSI information out */

      
/*       //calcualte total force and write */
/*       PetscReal Fp_x=0., Fv_x=0., dA; */
/*       PetscReal Fp_y=0., Fv_y=0.; */
/*       PetscReal Fp_z=0., Fv_z=0.; */
/*       for (i=0;i<n_elmt;i++) { */
/* 	dA =  ibm->dA[i]; */
/* 	Fp_x -= pressure[i]*nf_x[i]*dA; */
/* 	Fp_y -= pressure[i]*nf_y[i]*dA; */
/* 	Fp_z -= pressure[i]*nf_z[i]*dA; */
   
/* 	Fv_x += shear_x[i]*dA; */
/* 	Fv_y += shear_y[i]*dA; */
/* 	Fv_z += shear_z[i]*dA; */
/*       } */

/*       fsi->F_x=Fp_x+Fv_x; */
/*       fsi->F_y=Fp_y+Fv_y; */
/*       fsi->F_z=Fp_z+Fv_z; */

/*       PetscInt   nv1,nv2,nv3, elmt; */
/*       PetscReal  M_z, p; */
/*       PetscReal  n_x,n_y, F_x, F_y, r_x, r_y; */
      
/*       M_z= 0.; */
/*       for (elmt=0; elmt<n_elmt; elmt++) { */
/* 	nv1=ibm->nv1[elmt]; */
/* 	nv2=ibm->nv2[elmt]; */
/* 	nv3=ibm->nv3[elmt]; */
	
/* 	r_x= ibm->cent_x[elmt]; */
/* 	r_y= ibm->cent_y[elmt]; */
	
/* 	dA = ibm->dA[elmt]; */
	
/* 	n_x= ibm->nf_x[elmt]; */
/* 	n_y= ibm->nf_y[elmt]; */
	
/* 	p = pressure[elmt]; */
	
/* 	F_x = -p*dA*n_x; */
/* 	F_y = -p*dA*n_y; */
	
/* 	M_z +=   r_x*F_y - r_y*F_x; */
/*       } */

/*       PetscPrintf(PETSC_COMM_WORLD, "The Moment %le\n", M_z);  */

/*       sprintf(filen, "Force_IBM%2.2d.dat", user->_this); */
/*       f = fopen(filen, "a"); */
      
/*       fprintf(f,"%d %le %le %le %le %le %le %le %le %le\n",ti,Fp_x, Fp_y,Fp_z, */
/* 	      Fv_x, Fv_y, Fv_z, Fp_x+Fv_x, Fp_y+Fv_y, Fp_z+Fv_z); */
/*       fclose(f);       */

/*       PetscPrintf(PETSC_COMM_WORLD, "IBM forces: %d %le %le %le %le %le %le %le %le %le\n",ti,Fp_x, Fp_y,Fp_z, */
/* 	      Fv_x, Fv_y, Fv_z, Fp_x+Fv_x, Fp_y+Fv_y, Fp_z+Fv_z); */
                
/*   } //End of root work */

/*   //Erase all memory */
/*   PetscFree(shear); */
/*   PetscFree(shear_x); */
/*   PetscFree(shear_y); */
/*   PetscFree(shear_z); */
/*   PetscFree(vx); */
/*   PetscFree(vy); */
/*   PetscFree(vz); */
/*   PetscFree(pressure); */

/*   for (k=0; k<ncz; k++) { */
/*     for (j=0; j<ncy; j++) { */
/*       for (i=0; i<ncx; i++) { */
/* 	destroy(&cell_trg[k*ncx*ncy+j*ncx+i]); */
/*       } */
/*     } */
/*   } */

/*   PetscFree(cell_trg); */

  
/*   DMDAVecRestoreArray(fda, user->Cent, &coor); */
/*   DMDAVecRestoreArray(fda, user->Ucat, &ucat); */
/*   DMDAVecRestoreArray(da, user->P, &prs); */
/*   DMDAVecRestoreArray(da, user->Nvert, &nvert); */

/*   return(0); */

/* } */

PetscErrorCode jacobi(PetscReal a[3][3],PetscInt n, 
		      PetscReal d[3],PetscReal v[3][3],
		      PetscInt *nrot) 
//!--------------------------------------------------------------------  
{
  // PetscPrintf(PETSC_COMM_WORLD, "jacobi \n");

  PetscInt  NMAX=500;  
  PetscInt  i, ip, iq, j;
  PetscReal c, g, h, s, sm,t,tau,theta,tresh,b[NMAX],z2[NMAX];
 
  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) {
      v[ip][iq]=0.;
      v[ip][ip]=1.;
    }
  }


  for (ip=0; ip<n; ip++) {
    b[ip]=a[ip][ip];
    d[ip]=b[ip];
    z2[ip]=0.;
  }

  *nrot=0;

  for (i=0; i<5; i++) {
    sm=0.;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {                     
	sm=sm+fabs(a[ip][iq]);
      }
    }

    if(sm==0.)return 0;
    
    if(i<4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.;

    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++){
	g=100.*fabs(a[ip][iq]);
	if((i>4) && (fabs(d[ip])+g==fabs(d[ip]))&&
	   (fabs(d[iq])+g==fabs(d[iq])))
	  a[ip][iq]=0.;
	else if(fabs(a[ip][iq])>tresh){
	  h=d[iq]-d[ip];
	  if(abs(h)+g==abs(h))
	    t=a[ip][iq]/h;
	  else { 
	    theta=0.5*h/a[ip][iq];
	    t=1./(abs(theta)+sqrt(1.+theta*theta));
	    if(theta<0.)t=-t;
	  }
	  c=1./sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.+c);
	  h=t*a[ip][iq];
	  z2[ip]=z2[ip]-h;
	  z2[iq]=z2[iq]+h;
	  d[ip]=d[ip]-h;
	  d[iq]=d[iq]+h;
	  a[ip][iq]=0.;
	  
	  for (j=0; j<ip-1; j++) {
	    g=a[j][ip];
	    h=a[j][iq];
	    a[j][ip]=g-s*(h+g*tau);
	    a[j][iq]=h+s*(g-h*tau);
	  }
	  
	  
	  for (j=ip+1; j<iq-1; j++) {                 
	    g=a[ip][j];
	    h=a[j][iq];
	    a[ip][j]=g-s*(h+g*tau);
	    a[j][iq]=h+s*(g-h*tau);
	  }
	  
	  for (j=iq+1; j<n; j++) {                  
	    g=a[ip][j];
	    h=a[iq][j];
	    a[ip][j]=g-s*(h+g*tau);
	    a[iq][j]=h+s*(g-h*tau);
	  }
	  
	  for (j=0; j<n; j++) {                 
	    g=v[j][ip];
	    h=v[j][iq];
	    v[j][ip]=g-s*(h+g*tau);
	    v[j][iq]=h+s*(g-h*tau);
	  }
	  
	  *nrot++;
	}	
      }        
    }


    for (ip=0; ip<n; ip++) { 
      b[ip]=b[ip]+z2[ip];
      d[ip]=b[ip];
      z2[ip]=0.;
    }
    
  }

  PetscPrintf(PETSC_COMM_WORLD, "too many iterations in jacobi! code aborted!");
  return 1;
}

PetscErrorCode eigsrt(PetscReal d[3], PetscReal v[3][3], 
		      PetscInt n, PetscInt npp) 
//!--------------------------------------------------------------------
{
  PetscInt  i,j,k;
  PetscReal  p;

  for (i=0; i<n-1; i++) {  

    k=i;
    p=d[i];

    for (j=i+1; j<n; j++) {       
      if(d[j]>p) {
	k=j;
	p=d[j];
      }
    }

    if(k!=i){
      d[k]=d[i];
      d[i]=p;
      for (j=0; j<n; j++) {	          
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }

  return 0;
}

PetscErrorCode QCriteria(UserCtx *user)
{

  PetscInt bi;

  for (bi=0; bi<block_number; bi++) {
  DMDALocalInfo	info = user[bi].info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert, ***l2;

  PetscInt   njac, nrot;
  PetscReal  a[3][3], v[3][3], d[3];

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
  DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
  DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
  DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);

  DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
  DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
  DMDAVecGetArray(user[bi].da, user[bi].P, &q);
  DMDAVecGetArray(user[bi].da, user[bi].Phi, &l2);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal ss11, ss12, ss13, ss21, ss22, ss23, ss31, ss32, ss33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  PetscReal ww11, ww12, ww13, ww21, ww22, ww23, ww31, ww32, ww33;
  PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (i==1 || nvert[k][j][i-1] > 2.) {
	  uc = (ucat[k][j][i+1].x - ucat[k][j][i].x);
	  vc = (ucat[k][j][i+1].y - ucat[k][j][i].y);
	  wc = (ucat[k][j][i+1].z - ucat[k][j][i].z);
	}
	else if (i==mx-2 || nvert[k][j][i+1] > 2.) {
	  uc = (ucat[k][j][i].x - ucat[k][j][i-1].x);
	  vc = (ucat[k][j][i].y - ucat[k][j][i-1].y);
	  wc = (ucat[k][j][i].z - ucat[k][j][i-1].z);
	}
	else {
	  uc = (ucat[k][j][i+1].x - ucat[k][j][i-1].x)*0.5;
	  vc = (ucat[k][j][i+1].y - ucat[k][j][i-1].y)*0.5;
	  wc = (ucat[k][j][i+1].z - ucat[k][j][i-1].z)*0.5;
	}

	if (j==1 || nvert[k][j-1][i] > 2) {
	  ue = (ucat[k][j+1][i].x - ucat[k][j][i].x);
	  ve = (ucat[k][j+1][i].y - ucat[k][j][i].y);
	  we = (ucat[k][j+1][i].z - ucat[k][j][i].z);
	}
	else if (j==my-2 || nvert[k][j+1][i] > 2) {
	  ue = (ucat[k][j][i].x - ucat[k][j-1][i].x);
	  ve = (ucat[k][j][i].y - ucat[k][j-1][i].y);
	  we = (ucat[k][j][i].z - ucat[k][j-1][i].z);
	}
	else {
	  ue = (ucat[k][j+1][i].x - ucat[k][j-1][i].x) * 0.5;
	  ve = (ucat[k][j+1][i].y - ucat[k][j-1][i].y) * 0.5;
	  we = (ucat[k][j+1][i].z - ucat[k][j-1][i].z) * 0.5;
	}

	if (k==1 || nvert[k-1][j][i] > 2) {
	  uz = (ucat[k+1][j][i].x - ucat[k][j][i].x);
	  vz = (ucat[k+1][j][i].y - ucat[k][j][i].y);
	  wz = (ucat[k+1][j][i].z - ucat[k][j][i].z);
	}
	else if (k==mz-2 || nvert[k+1][j][i] > 2) {
	  uz = (ucat[k][j][i].x - ucat[k-1][j][i].x);
	  vz = (ucat[k][j][i].y - ucat[k-1][j][i].y);
	  wz = (ucat[k][j][i].z - ucat[k-1][j][i].z);
	}
	else {
	  uz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x) * 0.5;
	  vz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y) * 0.5;
	  wz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z) * 0.5;
	}

	csi1 = (csi[k][j][i].x + csi[k][j][i-1].x) * 0.5 * aj[k][j][i];
	csi2 = (csi[k][j][i].y + csi[k][j][i-1].y) * 0.5 * aj[k][j][i];
	csi3 = (csi[k][j][i].z + csi[k][j][i-1].z) * 0.5 * aj[k][j][i];

	eta1 = (eta[k][j][i].x + eta[k][j-1][i].x) * 0.5 * aj[k][j][i];
	eta2 = (eta[k][j][i].y + eta[k][j-1][i].y) * 0.5 * aj[k][j][i];
	eta3 = (eta[k][j][i].z + eta[k][j-1][i].z) * 0.5 * aj[k][j][i];

	zet1 = (zet[k][j][i].x + zet[k-1][j][i].x) * 0.5 * aj[k][j][i];
	zet2 = (zet[k][j][i].y + zet[k-1][j][i].y) * 0.5 * aj[k][j][i];
	zet3 = (zet[k][j][i].z + zet[k-1][j][i].z) * 0.5 * aj[k][j][i];

	//  derivatives for vorticity
	d11 = uc * csi1 + ue * eta1 + uz * zet1;
	d12 = uc * csi2 + ue * eta2 + uz * zet2;
	d13 = uc * csi3 + ue * eta3 + uz * zet3;

	d21 = vc * csi1 + ve * eta1 + vz * zet1;
	d22 = vc * csi2 + ve * eta2 + vz * zet2;
	d23 = vc * csi3 + ve * eta3 + vz * zet3;

	d31 = wc * csi1 + we * eta1 + wz * zet1;
	d32 = wc * csi2 + we * eta2 + wz * zet2;
	d33 = wc * csi3 + we * eta3 + wz * zet3;

	//!...Sij, Wij tensors
	s11 = d11 + d11;
	s12 = d12 + d21;
	s13 = d13 + d31;

	s21 = s12;
	s22 = d22 + d22;
	s23 = d23 + d32;

	s31 = s13;
	s32 = s23;
	s33 = d33 + d33;

	so = (s11 * s11 +
	      s22 * s22 +
	      s33 * s33) / 4. +
	  (s12 * s12 + s13 * s13 + s23 * s23) / 2.;

	//!...Sij^2

	ss11=s11*s11+s21*s12+s31*s13;
	ss12=s12*s11+s22*s12*s32*s13;
	ss13=s13*s11+s23*s12*s33*s13;
	ss21=s11*s21+s21*s22+s31*s23;
	ss22=s12*s21+s22*s22*s32*s23;
	ss23=s13*s21+s23*s22*s33*s23;
	ss31=s11*s31+s21*s32+s31*s33;
	ss32=s12*s31+s22*s32*s32*s33;
	ss33=s13*s31+s23*s32*s33*s33;

	w11 = 0;
	w12 = d12 - d21;
	w13 = d13 - d31;
	w21 = -w12;
	w22 = 0.;
	w23 = d23 - d32;
	w31 = -w13;
	w32 = -w23;
	w33 = 0.;

	wo = (w12 * w12 + w13 * w13 + w23 * w23) / 2.;

	//!...Wij^2
	//!
	ww11=w11*w11+w21*w12+w31*w13;
	ww12=w12*w11+w22*w12*w32*w13;
	ww13=w13*w11+w23*w12*w33*w13;
	ww21=w11*w21+w21*w22+w31*w23;
	ww22=w12*w21+w22*w22*w32*w23;
	ww23=w13*w21+w23*w22*w33*w23;
	ww31=w11*w31+w21*w32+w31*w33;
	ww32=w12*w31+w22*w32*w32*w33;
	ww33=w13*w31+w23*w32*w33*w33;

	a[0][0]=ss11+ww11;
	a[0][1]=ss12+ww12;
	a[0][2]=ss13+ww13;
	a[1][0]=ss21+ww21;
	a[1][1]=ss22+ww22;
	a[1][2]=ss23+ww23;
	a[2][0]=ss31+ww31;
	a[2][1]=ss32+ww32;
	a[2][2]=ss33+ww33;

	njac = 3;

	d[0]=0.;
	d[1]=0.;
	d[2]=0.;

	jacobi(a, njac, d, v, &nrot);
	eigsrt(d, v, njac, njac);


	l2[k][j][i]= d[1];

	q[k][j][i] = (wo - so) / 2.;
      }
    }
  }

  DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
  DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
  DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
  DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);

  DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
  DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
  DMDAVecRestoreArray(user[bi].da, user[bi].P, &q);
  DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &l2);
  } //bi
  return 0;
}

PetscErrorCode MaxVelocity(UserCtx *user, PetscInt ti)
{
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k, ii,jj,kk;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert, ***l2;

  PetscInt   njac, nrot;
  PetscReal  a[3][3], v[3][3], d[3], II[4],x[3];

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user->fda, user->Ucat, &ucat);
  DMDAVecGetArray(user->fda, user->Csi, &csi);
  DMDAVecGetArray(user->fda, user->Eta, &eta);
  DMDAVecGetArray(user->fda, user->Zet, &zet);

  DMDAVecGetArray(user->da, user->Aj, &aj);
  DMDAVecGetArray(user->da, user->Nvert, &nvert);
  DMDAVecGetArray(user->da, user->P, &q);
  DMDAVecGetArray(user->da, user->Phi, &l2);

  PetscErrorCode imag;

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
 
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal smax=0., smax_global, smin;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] <.1) {
	  if (smax < fabs(ucat[k][j][i].x)) 
	    smax = fabs(ucat[k][j][i].x) ;	  
	}
      }
    }
  }
  MPI_Allreduce(&smax, &smax_global,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
 
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Max_U.dat");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le\n",ti, smax_global);
    fclose(f);
  }

  DMDAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(user->fda, user->Csi, &csi);
  DMDAVecRestoreArray(user->fda, user->Eta, &eta);
  DMDAVecRestoreArray(user->fda, user->Zet, &zet);

  DMDAVecRestoreArray(user->da, user->Aj, &aj);
  DMDAVecRestoreArray(user->da, user->Nvert, &nvert);
  DMDAVecRestoreArray(user->da, user->P, &q);
  DMDAVecRestoreArray(user->da, user->Phi, &l2);

  return 0;
}

//--------------------------------------------------------------------
PetscErrorCode Pressure_Difference(UserCtx *user) {

  PetscInt bi,ip=1,jp=1,kp=1;
  
  PetscOptionsGetInt(PETSC_NULL, "-ip", &ip, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-jp", &jp, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-kp", &kp, PETSC_NULL);
  
  PetscPrintf(PETSC_COMM_WORLD, "pressure reffrence point %d %d %d \n",ip,jp,kp);

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    PetscInt	lxs, lxe, lys, lye, lzs, lze;
    
    PetscInt i, j, k;
    PetscReal p_ref=0.0;
    lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
    
    
    PetscReal  ***q;
    
    if (lxs==0) lxs++;
    if (lxe==mx) lxe--;
    if (lys==0) lys++;
    if (lye==my) lye--;
    if (lzs==0) lzs++;
    if (lze==mz) lze--;
    
    
    DMDAVecGetArray(user[bi].da, user[bi].P, &q);
    p_ref=q[kp][jp][ip];
    PetscPrintf(PETSC_COMM_WORLD, "pressure at reffrence point %le \n",p_ref);
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  q[k][j][i]=q[k][j][i]-p_ref;
	}
      }
    }
    
    DMDAVecRestoreArray(user[bi].da, user[bi].P, &q);
    
  } //bi

  return 0;
}

//------------------------------------------------------------
PetscErrorCode pressurediff(UserCtx *user, PetscInt ti)
{
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k, ii,jj,kk;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  PetscReal pdiff,***p;

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user->da, user->P, &p);
  // Anatomic

  // Straight
  k=50 ; j=101; i=101;
  kk=83;jj=101;ii=101;

  pdiff=p[k][j][i]-p[kk][jj][ii];

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Press_diff.dat");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le\n",ti, pdiff);
    fclose(f);
  }


  DMDAVecRestoreArray(user->da, user->P, &p);

  return 0;
}

PetscErrorCode histogram(UserCtx *user, PetscInt ti)
{
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze, *h;

  PetscInt      i, j, k, ii,im=0,jm=0,km=0, N=200;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  PetscReal    ***p, ***nvert, MaxS=-1e10;
  PetscReal    Rmin=0., Rmax=0.02, dR;

  dR= (Rmax-Rmin)/ N;

  PetscMalloc(N*sizeof(PetscInt), &h);


  for (ii=0;ii<N;ii++) {
    h[ii]=0;
  }
    
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user->da, user->P, &p);
  DMDAVecGetArray(user->da, user->Nvert, &nvert);
  PetscInt count=0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] <.1) {
	  count++;
	  if (p[k][j][i]>MaxS) {
	    MaxS=p[k][j][i];
	    im=i;jm=j;km=k;
	  }
	  ii=(int)((p[k][j][i] - Rmin)/dR);
	  if (ii>N-1) ii=N-1;
	  if (ii<0) ii=0;
	  h[ii]=h[ii]+1;
	}
      }
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "dR %le %le %d %d %d\n",dR, MaxS, h[0], h[10], count);    

  DMDAVecRestoreArray(user->da, user->P, &p);
  DMDAVecRestoreArray(user->da, user->Nvert, &nvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "histogram%3.1d.dat",ti);
    f = fopen(filen, "a");
    
    for (ii=0;ii<N;ii++) {
     
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %d\n",(ii+0.5)*dR+Rmin,h[ii]);
    }
       PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); 
    fclose(f);

  }

  return 0;
}

//--------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv) {

  DM	da, fda;
  Vec	qn, qnm;
  Vec	c;
  UserCtx	*user;
  
  PetscErrorCode ierr;
  
  IBMNodes	*ibm, *ibm0, *ibm1;
  FSInfo        *fsi;
  
  PetscInitialize(&argc, &argv, (char *)0, help);
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  
  PetscInt rank, bi, ibi,i;
  
  //PetscMalloc(sizeof(IBMNodes), &ibm0);
   
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBool fish = PETSC_FALSE;
  PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);
  PetscBool QCR = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-qcr", &QCR, PETSC_NULL);
  PetscBool SHR = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-shr", &SHR, PETSC_NULL);
  PetscBool EXT = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-ext", &EXT, PETSC_NULL);
  PetscBool UUU = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-uuu", &UUU, PETSC_NULL);
  PetscBool PRSS = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-prs", &PRSS, PETSC_NULL);
  PetscBool PRS = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-prsd", &PRS, PETSC_NULL);
  PetscBool HST = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-hst", &HST, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL);
  PetscInt phave=0;
  PetscOptionsGetInt(PETSC_NULL, "-phave", &phave, PETSC_NULL); 
  
  PetscInt rans=0;
  PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);
  PetscInt vtk=0;
  PetscOptionsGetInt(PETSC_NULL, "-vtk", &vtk, PETSC_NULL);
  
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  
  PetscOptionsGetInt(PETSC_NULL, "-body", &NumberOfBodies, PETSC_NULL);
  
  PetscOptionsGetReal(PETSC_NULL, "-x_c", &(CMx_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-y_c", &(CMy_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-z_c", &(CMz_c), PETSC_NULL);
  
  PetscOptionsGetInt(PETSC_NULL, "-vis_flg", &visflg, PETSC_NULL);
  PetscOptionsGetString(PETSC_NULL, "-orient", orient,sizeof(orient),PETSC_NULL);
  
  PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
  
  PetscInt generate_grid=0;
  PetscOptionsGetInt(PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);
  
  if (generate_grid) {
    block_number=1;
    PetscOptionsGetInt(PETSC_NULL, "-block_number", &block_number, PETSC_NULL);
  } else {
    if (!rank) {
      FILE *fd;
      fd = fopen("grid.dat", "r");
      
      fscanf(fd, "%i\n", &block_number);
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      fclose(fd);
    }
    else {
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }
  
  PetscMalloc(block_number*sizeof(UserCtx), &user);
  
  ReadCoordinates(user);
  
  PetscPrintf(PETSC_COMM_WORLD, "read coord!\n");
  
  for (bi=0; bi<block_number; bi++) {
    ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat); CHKERRQ(ierr);
    VecDuplicate(user[bi].Ucat, &user[bi].Csi);
    VecDuplicate(user[bi].Ucat, &user[bi].Eta);
    VecDuplicate(user[bi].Ucat, &user[bi].Zet);
    VecDuplicate(user[bi].Ucat, &user[bi].Cent);
    VecDuplicate(user[bi].Ucat, &user[bi].Ucat_o);
    
    if (averaging) {
      VecDuplicate(user[bi].Ucat, &user[bi].Ucat_cross_sum);
      VecDuplicate(user[bi].Ucat, &user[bi].Ucat_square_sum);
    }

    ierr = DMCreateGlobalVector(user[bi].da, &user[bi].P); CHKERRQ(ierr);
    VecDuplicate(user[bi].P, &user[bi].shrr);

    VecDuplicate(user[bi].P, &user[bi].Nvert);
    //amir 

    /* VecDuplicate(user[bi].P, &user[bi].Nvert1); */
    /* VecDuplicate(user[bi].P, &user[bi].Nvert2); */
    
    /* VecDuplicate(user[bi].P, &user[bi].Nvert3); */
    
    VecDuplicate(user[bi].P, &user[bi].Aj);
    VecDuplicate(user[bi].P, &user[bi].Phi);
    VecDuplicate(user[bi].P, &user[bi].lNu_t);
    VecDuplicate(user[bi].P, &user[bi].lCs);
    if (averaging)
      VecDuplicate(user[bi].P, &user[bi].lUstar);
    
    ierr = FormMetrics(&(user[bi]));
    
    PetscOptionsGetReal(PETSC_NULL, "-ren",&user[bi].ren, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-dt", &user[bi].dt, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "Re %le \n",user[bi].ren);
  }

  PetscBool flag;

  PetscInt tis, tie, tsteps=5;
  PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");    
  }
  
  PetscOptionsGetInt(PETSC_NULL, "-tie", &tie, &flag);
  if (!flag) {
    tie = tis;
  }
  
  PetscOptionsGetInt(PETSC_NULL, "-ts", &tsteps, &flag);
  if (!flag) {
    tsteps = 10; /* Default increasement is 5 */
  }
  
  if (PRSS) {
    bi=0;
    PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
    
    PetscMalloc(NumberOfBodies*sizeof(IBMList), &(user[bi].ibmlist));
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      InitIBMList(&(user[bi].ibmlist[ibi]));
      
      if (fish) ibm_read_Icem(&ibm[ibi],ibi);
      else      ibm_read_Icem(&ibm[ibi], ibi);
      
      if (fish) fish_init(&(user[bi].dt));
      PetscPrintf(PETSC_COMM_WORLD, "fish initial %d\n", ibi);    
      
      /*   if (fish)  fish_swim(&ibm[ibi], (ti)*user[bi].dt, user[bi].dt); */
      /*      PetscPrintf(PETSC_COMM_WORLD, "fish swim %d\n",ti);     */

      ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
      PetscPrintf(PETSC_COMM_WORLD, "ibm search  %d\n",ibi); 
      
      PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
      FsiInitialize(0, &fsi[ibi], ibi);
      fsi[ibi].x_c=0.;
      fsi[ibi].y_c=0.;
      fsi[ibi].z_c=0.;
    }
  }
  
  PetscOptionsGetInt(PETSC_NULL, "-wing", &wing, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rframe", &rotateframe, PETSC_NULL);
  
  if (wing) {
    ibi=0;bi=0;
    ibm_read_Ansys(&ibm[ibi],ibi);
    
    PetscPrintf(PETSC_COMM_WORLD, "Ibm read wing!\n");
    PetscReal scale[7];
    ibm_scale(&ibm[ibi],scale);
    calc_ibm_normal(&ibm[ibi]);
    
    PetscBarrier(PETSC_NULL);
    PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
    FsiInitialize(0, &fsi[ibi], ibi);
    fsi[ibi].x_c=0.;
    fsi[ibi].y_c=0.;
    fsi[ibi].z_c=0.;
    PetscPrintf(PETSC_COMM_WORLD, "wing init done!! dt %le \n", user[0].dt);	            
  }
  
  /* Time STEPPING */
 
  FILE *f;
  char filen[80];
  for (ti=tis; ti<=tie; ti+=tsteps) {
     
    /* #ifndef TECIO */
    /*     sprintf(filen, "Result%3.3d.dat", ti); */
    /*     f = fopen(filen, "w"); */
    
    /*     if (rans)  */
    /*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables = X, Y, Z, U, V, W, P, Nv, K, O\n"); */
    /*     else */
    /*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables = X, Y, Z, U, V, W, P, Nv\n"); */
    /*     fclose(f); */
    /* #endif */
    for (bi=0; bi<block_number; bi++) {
	if (averaging || phave) 
	  Ucont_P_Binary_Input_Ave(&user[bi], averaging);
	else {
	  PetscInt flgI=0;
	  if (rans) flgI=1;
	  Ucont_P_Binary_Input(&user[bi]);
	  Ucont_P_Binary_Input1(&user[bi], flgI);
	}
    }
    
    //#ifdef TECIO
    
    if (wing) {
      Cmpnts alfa_c, a_c, omega_c, rr;
      bi=0;
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	wing_motion(&ibm[ibi], &fsi[ibi], ti-1, &user[bi].dt);
	if (block_number>1)
	  for (bi=1; bi<block_number; bi++)  user[bi].dt=user[0].dt;
	
	if (rotateframe) {
	  
	  PetscPrintf(PETSC_COMM_WORLD, "RFRAME IBM WING Motion \n");
	  // omega_c is frame angular speed
	  omega_c.x=fsi[ibi].S_ang_n[1];
	  omega_c.y=fsi[ibi].S_ang_n[3];
	  omega_c.z=fsi[ibi].S_ang_n[5];	     
	    
	  // axis of frame rotation
	  a_c.x=0.;//fsi[ibi].x_c;
	  a_c.y=0.;//fsi[ibi].y_c;
	  a_c.z=0.;//fsi[ibi].z_c;	    
	  
	  // angle of rotation
	  alfa_c.x=fsi[ibi].S_ang_n[0];
	  alfa_c.y=fsi[ibi].S_ang_n[2];
	  alfa_c.z=fsi[ibi].S_ang_n[4];
	  
	  omega_c = rotateT_y(omega_c,a_c, alfa_c);
	  PetscPrintf(PETSC_COMM_WORLD, "frame omega %le %le %le wing %le %le %le\n", omega_c.x, omega_c.y, omega_c.z,fsi[ibi].S_ang_n[1],fsi[ibi].S_ang_n[3],fsi[ibi].S_ang_n[5] );
	  
	  // wing velocity
	  for (i=0; i<ibm[ibi].n_v; i++) {	  
	    rr.x = ibm[ibi].x_bp0[i];
	    rr.y = ibm[ibi].y_bp0[i];
	    rr.z = ibm[ibi].z_bp0[i];
	    //	  rr   = rotate_xyz(rr, a_c, alfa_c);
	    rr   = MINUS(rr, a_c);
	    
	    ibm[ibi].uold[i]= ibm[ibi].u[i];
	  }
	}
	
	wing_motion(&ibm[ibi], &fsi[ibi], ti, &user[bi].dt);
	
	if (rotateframe) {
	    
	  PetscPrintf(PETSC_COMM_WORLD, "RFRAME IBM WING Motion \n");
	  // omega_c is frame angular speed
	  omega_c.x=fsi[ibi].S_ang_n[1];
	  omega_c.y=fsi[ibi].S_ang_n[3];
	  omega_c.z=fsi[ibi].S_ang_n[5];	     
	  
	  // axis of frame rotation
	  a_c.x=0.;//fsi[ibi].x_c;
	  a_c.y=0.;//fsi[ibi].y_c;
	  a_c.z=0.;//fsi[ibi].z_c;	    
	  
	  // angle of rotation
	  alfa_c.x=fsi[ibi].S_ang_n[0];
	  alfa_c.y=fsi[ibi].S_ang_n[2];
	  alfa_c.z=fsi[ibi].S_ang_n[4];
	  
	  omega_c = rotateT_y(omega_c,a_c, alfa_c);
	  PetscPrintf(PETSC_COMM_WORLD, "frame omega %le %le %le wing %le %le %le\n", omega_c.x, omega_c.y, omega_c.z,fsi[ibi].S_ang_n[1],fsi[ibi].S_ang_n[3],fsi[ibi].S_ang_n[5] );
	    
	  // wing velocity
	  for (i=0; i<ibm[ibi].n_v; i++) {	  
	    rr.x = ibm[ibi].x_bp0[i];
	    rr.y = ibm[ibi].y_bp0[i];
	    rr.z = ibm[ibi].z_bp0[i];
	    //	  rr   = rotate_xyz(rr, a_c, alfa_c);
	    rr   = MINUS(rr, a_c);
	    
	    ibm[ibi].u[i] = Cross(omega_c, rr);	      
	  }
	  
	  for (bi=0; bi<block_number; bi++) {
	    // ibm_intp_pj(&user[bi], &ibm[ibi], &fsi[ibi]);
	      rotate_force_data(&fsi[ibi], ibi, bi);
	  }
	  
	} else {
	  rotate_ibm(&ibm[ibi],&fsi[ibi]);
	  calc_ibm_normal(&ibm[ibi]);	
	  if (ti>0) calc_ibm_velocity(&ibm[ibi], user[bi].dt);
	  
	  // ibm_intp_pj(&user[bi], &ibm[ibi], &fsi[ibi]);
	  PetscPrintf(PETSC_COMM_WORLD, "IBM WING Motion \n");
	}
      }
    }  // if wing
    else if (EXT) {
      //ExtractIO(user);
    }
    else if (QCR) {
      PetscPrintf(PETSC_COMM_WORLD, "Q Crit!\n");
      
      QCriteria(user);
      PetscPrintf(PETSC_COMM_WORLD, "Q Crit done!\n");
      if (vtk) 
	VTKOutQ(user);
      else 
	// TECIOOutQ(user);
	PetscPrintf(PETSC_COMM_WORLD, "Q Crit!\n");
    }
    else if (SHR) {
      PetscPrintf(PETSC_COMM_WORLD, "Shear max!\n");
      
      //	MaxStress(user);
      shear(user, ti);
      PetscPrintf(PETSC_COMM_WORLD, "Shear max done!\n");
      VTKOutShr(user);
      //	TECIOOutS(user);
      
    }
    else if (UUU) {
      PetscPrintf(PETSC_COMM_WORLD, "span U  max!\n");
      
      MaxVelocity(user,ti);
      PetscPrintf(PETSC_COMM_WORLD, "U max done!\n");
      
    } else if (PRSS) {
      Pressure_Difference(user);//Mohsen
      bi=0;
      for (ibi=0; ibi<NumberOfBodies; ibi++) {    
	if (fish) {
	  ibm_surface_VTKOut(&ibm[ibi], ibi, tis);  //Mohsen
	  PetscPrintf(PETSC_COMM_WORLD," user[bi].dt is %le \n",user[bi].dt);
	  fish_swim(&ibm[ibi], (ti-1)*user[bi].dt, user[bi].dt);
	  for (i=0; i<ibm[ibi].n_v; i++) {
	    ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
	    ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
	    ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];      
	  }
	  fish_swim(&ibm[ibi], (ti)*user[bi].dt, user[bi].dt);
	  PetscPrintf(PETSC_COMM_WORLD, "fish swim %d\n",ti); 
	  ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);	    
	}

	PetscPrintf(PETSC_COMM_WORLD, "pressure proj %d\n",ibi);    
	//if (ti==tis) Initialize_Projecting(&ibm[ibi]);
	PetscPrintf(PETSC_COMM_WORLD, "Initialization is done\n");
	Projecting(&user[bi],&ibm[ibi]);
	PetscPrintf(PETSC_COMM_WORLD, "Projecting is done\n");
	Projecting_VTKOut(&ibm[ibi],ibi,ti);
	if (fish)  Projecting_Power(&ibm[ibi],&fsi[ibi],ti);

	PetscPrintf(PETSC_COMM_WORLD, "Projecting VTKOut works \n");    
      }
    }
    else if (PRS) {
      PetscPrintf(PETSC_COMM_WORLD, "PRS!\n");
      
      pressurediff(user,ti);
      PetscPrintf(PETSC_COMM_WORLD, "PRS done!\n");
      
    }
    else if (HST) {
      PetscPrintf(PETSC_COMM_WORLD, "Histogram!\n");
      
      //	MaxStress(user);
      //	TECIOOutS(user);
      histogram(user,ti);
      PetscPrintf(PETSC_COMM_WORLD, "HST done!\n");
      
    } else { // if (!QCR & !SHR & !EXT & !UUU & !PRSS & !PRS & !HST) {
      if (vtk){
	Pressure_Difference(user);//Mohsen
	VTKOut(user);
      }else {
	/*   if (averaging)  */
	/* 	    TECIOOut_AVE(user); */
	/* 	  else */
	/* 	    TECIOOut(user); */
	  PetscPrintf(PETSC_COMM_WORLD, "HST done!\n");
      } 
    }
    /* #else */
    /*       for (bi=0; bi<block_number; bi++) { */
    /* 	f = fopen(filen, "a"); */
    
    /* 	PetscFPrintf(PETSC_COMM_WORLD, f, "Zone I=%d, J=%d, K=%d F=BLOCK VARLOCATION=([7-8]=CELLCENTERED)\n", user[bi].info.mx-1, user[bi].info.my-1, user[bi].info.mz-1); */
    /* 	fclose(f); */
    /* 	TecOut(&(user[bi])); */
    /*       } */
    /* #endif */
  }
  
  PetscFinalize();
}

//!--------------------------------------------------------------------

PetscErrorCode SolveCubic(PetscReal II[4],PetscReal x[3]) 
//!--------------------------------------------------------------------  
/*
Solution: To find roots to the following cubic equation where a, b, c, and d are real.

   a*x3 + b*x2 + c*x + d = 0

Formula:

  Step 1: Calculate p and q
          p = ( 3*c/a - (b/a)^2 ) / 3
          q = ( 2*(b/a)^3 - 9*b*c/a/a + 27*d/a ) / 27
  Step 2: Calculate discriminant D
          D = (p/3)^3 + (q/2)^2
  Step 3: Depending on the sign of D, you follow different strategy.
          If D<0, three distinct real roots.
          If D=0, three real roots of which at least two are equal.
          If D>0, one real and two complex roots.
  Step 3a: For D>0 and D=0
          Calculate u and v
          u = cubic_root(-q/2 + sqrt(D))
          v = cubic_root(-q/2 - sqrt(D))
          Find the three transformed roots
          y1 = u + v
          y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
          y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
  Step 3b: Alternately, for D<0, a trigonometric formulation is more convenient
          y1 =  2 * sqrt(|p|/3) * cos(phi/3)
          y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
          y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
          where phi = acos(-q/2/sqrt(|p|^3/27))
                pi  = 3.141592654...
  Step 4  Finally, find the three roots
          x = y - b/a/3

  Things to watch out for:
    1. Make sure you know what is integer, real, and complex.
    2. FORTRAN's SQRT (square root) function takes only non-negative arguments.
    3. FORTRAN's exponentiation is "**", not "^".
    4. There is no "cubic_root" function in FORTRAN; you do it by raising a
       number to a factional power, e.g., 0.333... for cubic root.
    5. Your can only raise a non-negative number to a fractional power.
       e.g., (-2.)**(1./3.) will not work; you need to issue -(2.)**(1./3.)
       And do not forget the decimal point.  For example, in computing 2.**(1/3),
       FORTRAN evaluates 1/3 first, which is 0, not 0.33... as you might expect.
    6  There is no "|" in FORTRAN.  The absolute value function is "ABS".
    7. FORTRAN is case-insensitive; do not simply use "d" (for coefficient)
       and "D" (for discriminant), as FORTRAN thinks both variables are equivalent.
*/

{
  PetscReal p,q,D,u,y1,y2,y3;
  PetscReal phi, pi  = 3.141592654;
  
  p = ( 3*II[2]/II[0] - pow((II[1]/II[0]), 2) ) / 3.;
  q = ( 2*pow((II[1]/II[0]),3) - 9*II[1]*II[2]/II[0]/II[0] + 27*II[3]/II[0] ) / 27.;

  D = pow(p/3,3) + pow(q/2,2);

  if (D<0.) {
    phi = acos(-q/2/sqrt(pow(fabs(p),3)/27));
    y1 =  2 * sqrt(fabs(p)/3) * cos(phi/3);
    y2 = -2 * sqrt(fabs(p)/3) * cos((phi+pi)/3);
    y3 = -2 * sqrt(fabs(p)/3) * cos((phi-pi)/3);
  } else if (fabs(D)<1e-6) {

    u = -pow(fabs(q)/2,1./3.);
    if (q<0) u=-u;

   
    y1 = 2*u ;
    y2 = -u;
    y3 = -u;
  } else {
    
    y1=0;
    y2=0.;
    y3=0.;
    PetscPrintf(PETSC_COMM_WORLD, "Imaginary Roots!\n");    
    return(1);
  }

 

  x[0] = y1 - II[1]/II[0]/3;

  x[1] = y2 - II[1]/II[0]/3;

  x[2] = y3 - II[1]/II[0]/3;

 
  return(0);
}

PetscErrorCode MaxStress(UserCtx *user)
{
  PetscInt bi;

  for (bi=0; bi<block_number; bi++) {

  DMDALocalInfo	info = user[bi].info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k, ii,jj,kk;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert, ***l2;

  PetscInt   njac, nrot;
  PetscReal  a[3][3], v[3][3], d[3], II[4],x[3];

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
  DMDAVecGetArray(user[bi].fda, user[bi].Csi, &csi);
  DMDAVecGetArray(user[bi].fda, user[bi].Eta, &eta);
  DMDAVecGetArray(user[bi].fda, user[bi].Zet, &zet);

  DMDAVecGetArray(user[bi].da, user[bi].Aj, &aj);
  DMDAVecGetArray(user[bi].da, user[bi].Nvert, &nvert);
  DMDAVecGetArray(user[bi].da, user[bi].P, &q);
  DMDAVecGetArray(user[bi].da, user[bi].Phi, &l2);

  PetscErrorCode imag;

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
 
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal smax, smin;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (i==1 || nvert[k][j][i-1] > 2.) {
	  uc = (ucat[k][j][i+1].x - ucat[k][j][i].x);
	  vc = (ucat[k][j][i+1].y - ucat[k][j][i].y);
	  wc = (ucat[k][j][i+1].z - ucat[k][j][i].z);
	}
	else if (i==mx-2 || nvert[k][j][i+1] > 2.) {
	  uc = (ucat[k][j][i].x - ucat[k][j][i-1].x);
	  vc = (ucat[k][j][i].y - ucat[k][j][i-1].y);
	  wc = (ucat[k][j][i].z - ucat[k][j][i-1].z);
	}
	else {
	  uc = (ucat[k][j][i+1].x - ucat[k][j][i-1].x)*0.5;
	  vc = (ucat[k][j][i+1].y - ucat[k][j][i-1].y)*0.5;
	  wc = (ucat[k][j][i+1].z - ucat[k][j][i-1].z)*0.5;
	}

	if (j==1 || nvert[k][j-1][i] > 2) {
	  ue = (ucat[k][j+1][i].x - ucat[k][j][i].x);
	  ve = (ucat[k][j+1][i].y - ucat[k][j][i].y);
	  we = (ucat[k][j+1][i].z - ucat[k][j][i].z);
	}
	else if (j==my-2 || nvert[k][j+1][i] > 2) {
	  ue = (ucat[k][j][i].x - ucat[k][j-1][i].x);
	  ve = (ucat[k][j][i].y - ucat[k][j-1][i].y);
	  we = (ucat[k][j][i].z - ucat[k][j-1][i].z);
	}
	else {
	  ue = (ucat[k][j+1][i].x - ucat[k][j-1][i].x) * 0.5;
	  ve = (ucat[k][j+1][i].y - ucat[k][j-1][i].y) * 0.5;
	  we = (ucat[k][j+1][i].z - ucat[k][j-1][i].z) * 0.5;
	}

	if (k==1 || nvert[k-1][j][i] > 2) {
	  uz = (ucat[k+1][j][i].x - ucat[k][j][i].x);
	  vz = (ucat[k+1][j][i].y - ucat[k][j][i].y);
	  wz = (ucat[k+1][j][i].z - ucat[k][j][i].z);
	}
	else if (k==mz-2 || nvert[k+1][j][i] > 2) {
	  uz = (ucat[k][j][i].x - ucat[k-1][j][i].x);
	  vz = (ucat[k][j][i].y - ucat[k-1][j][i].y);
	  wz = (ucat[k][j][i].z - ucat[k-1][j][i].z);
	}
	else {
	  uz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x) * 0.5;
	  vz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y) * 0.5;
	  wz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z) * 0.5;
	}

	csi1 = (csi[k][j][i].x + csi[k][j][i-1].x) * 0.5 * aj[k][j][i];
	csi2 = (csi[k][j][i].y + csi[k][j][i-1].y) * 0.5 * aj[k][j][i];
	csi3 = (csi[k][j][i].z + csi[k][j][i-1].z) * 0.5 * aj[k][j][i];

	eta1 = (eta[k][j][i].x + eta[k][j-1][i].x) * 0.5 * aj[k][j][i];
	eta2 = (eta[k][j][i].y + eta[k][j-1][i].y) * 0.5 * aj[k][j][i];
	eta3 = (eta[k][j][i].z + eta[k][j-1][i].z) * 0.5 * aj[k][j][i];

	zet1 = (zet[k][j][i].x + zet[k-1][j][i].x) * 0.5 * aj[k][j][i];
	zet2 = (zet[k][j][i].y + zet[k-1][j][i].y) * 0.5 * aj[k][j][i];
	zet3 = (zet[k][j][i].z + zet[k-1][j][i].z) * 0.5 * aj[k][j][i];

	//  derivatives for vorticity
	d11 = uc * csi1 + ue * eta1 + uz * zet1;
	d12 = uc * csi2 + ue * eta2 + uz * zet2;
	d13 = uc * csi3 + ue * eta3 + uz * zet3;

	d21 = vc * csi1 + ve * eta1 + vz * zet1;
	d22 = vc * csi2 + ve * eta2 + vz * zet2;
	d23 = vc * csi3 + ve * eta3 + vz * zet3;

	d31 = wc * csi1 + we * eta1 + wz * zet1;
	d32 = wc * csi2 + we * eta2 + wz * zet2;
	d33 = wc * csi3 + we * eta3 + wz * zet3;

	//!...Sij, Wij tensors
	s11 = d11 + d11;
	s12 = d12 + d21;
	s13 = d13 + d31;

	s21 = s12;
	s22 = d22 + d22;
	s23 = d23 + d32;

	s31 = s13;
	s32 = s23;
	s33 = d33 + d33;

	// a=tow_ij=1/Re*Sij
	a[0][0]=1./user[bi].ren*s11;
	a[0][1]=1./user[bi].ren*s12;
	a[0][2]=1./user[bi].ren*s13;
	a[1][0]=1./user[bi].ren*s21;
	a[1][1]=1./user[bi].ren*s22;
	a[1][2]=1./user[bi].ren*s23;
	a[2][0]=1./user[bi].ren*s31;
	a[2][1]=1./user[bi].ren*s32;
	a[2][2]=1./user[bi].ren*s33;

	// invariants
	II[0]=1.;II[1]=0.;II[2]=0.;II[3]=0.;
	for (ii=0; ii<3; ii++) {
	  II[1]+=a[ii][ii];
	  for (jj=0; jj<3; jj++) {
	    II[2]+=a[ii][jj]*a[ii][jj];
	    for (kk=0; kk<3; kk++) {
	      II[3]+=a[ii][jj]*a[jj][kk]*a[jj][kk];
	    }
	  }
	}

	II[1]=-II[1];
	II[3]=-II[3];
	imag=SolveCubic(II,x);
	
	if (imag==1)
	  PetscPrintf(PETSC_COMM_WORLD, "Imaginary Roots! at %d %d %d\n",i,j,k);    
	
	smax=x[0];smin=x[0];
	if (smax<x[1]) smax=x[1];
	if (smax<x[2]) smax=x[2];
	
	if (smin>x[1]) smin=x[1];
	if (smin>x[2]) smin=x[2];

	q[k][j][i] = (smax - smin) / 2.;
      }
    }
  }

  DMDAVecRestoreArray(user[bi].fda, user[bi].Ucat, &ucat);
  DMDAVecRestoreArray(user[bi].fda, user[bi].Csi, &csi);
  DMDAVecRestoreArray(user[bi].fda, user[bi].Eta, &eta);
  DMDAVecRestoreArray(user[bi].fda, user[bi].Zet, &zet);

  DMDAVecRestoreArray(user[bi].da, user[bi].Aj, &aj);
  DMDAVecRestoreArray(user[bi].da, user[bi].Nvert, &nvert);
  DMDAVecRestoreArray(user[bi].da, user[bi].P, &q);
  DMDAVecRestoreArray(user[bi].da, user[bi].Phi, &l2);
  } //bi
  return 0;
}

//!--------------------------------------------------------------------



PetscErrorCode Combine_Elmt(IBMNodes *ibm, IBMNodes *ibm0, IBMNodes *ibm1)
{

  PetscInt i;

  ibm->n_v = ibm0->n_v + ibm1->n_v;
  ibm->n_elmt = ibm0->n_elmt + ibm1->n_elmt;

  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  for (i=0; i<ibm0->n_v; i++) {
    ibm->x_bp[i] = ibm0->x_bp[i];
    ibm->y_bp[i] = ibm0->y_bp[i];
    ibm->z_bp[i] = ibm0->z_bp[i];

    ibm->u[i] = ibm0->u[i];
    ibm->uold[i] = ibm0->uold[i];
   
  }
  for (i=0; i<ibm0->n_elmt; i++) {
    ibm->nv1[i] = ibm0->nv1[i];
    ibm->nv2[i] = ibm0->nv2[i];
    ibm->nv3[i] = ibm0->nv3[i];

    ibm->nf_x[i] = ibm0->nf_x[i];
    ibm->nf_y[i] = ibm0->nf_y[i];
    ibm->nf_z[i] = ibm0->nf_z[i];
  }

  for (i=ibm0->n_v; i<n_v; i++) {
    ibm->x_bp[i] = ibm1->x_bp[i-ibm0->n_v];
    ibm->y_bp[i] = ibm1->y_bp[i-ibm0->n_v];
    ibm->z_bp[i] = ibm1->z_bp[i-ibm0->n_v];
    ibm->u[i].x = 0.;
    ibm->u[i].y = 0.;
    ibm->u[i].z = 0.;
  }

  for (i=ibm0->n_elmt; i<n_elmt; i++) {
    ibm->nv1[i] = ibm1->nv1[i-ibm0->n_elmt] + ibm0->n_v;
    ibm->nv2[i] = ibm1->nv2[i-ibm0->n_elmt] + ibm0->n_v;
    ibm->nv3[i] = ibm1->nv3[i-ibm0->n_elmt] + ibm0->n_v;

    ibm->nf_x[i] = ibm1->nf_x[i-ibm0->n_elmt];
    ibm->nf_y[i] = ibm1->nf_y[i-ibm0->n_elmt];
    ibm->nf_z[i] = ibm1->nf_z[i-ibm0->n_elmt];
  }

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == ti) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, 1670-96);
      for (i=0; i<n_v; i++) {

	/*    ibm->x_bp[i] = ibm->x_bp0[i];
	      ibm->y_bp[i] = ibm->y_bp0[i];
	      ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=96; i<n_elmt; i++) {
	if (fabs(ibm->nf_z[i]) > 0.5 ||
	    (fabs(ibm->nf_z[i]) < 0.5 &&
	     (ibm->x_bp[ibm->nv1[i]] * ibm->x_bp[ibm->nv1[i]] +
	      ibm->y_bp[ibm->nv1[i]] * ibm->y_bp[ibm->nv1[i]]) < 0.44*0.44)) {
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
	}
      }
      fclose(f);

      sprintf(filen, "leaflet%3.3d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, 96);
      for (i=0; i<n_v; i++) {

	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<96; i++) {
	if (fabs(ibm->nf_z[i]) > 0.5 ||
	    (fabs(ibm->nf_z[i]) < 0.5 &&
	     (ibm->x_bp[ibm->nv1[i]] * ibm->x_bp[ibm->nv1[i]] +
	      ibm->y_bp[ibm->nv1[i]] * ibm->y_bp[ibm->nv1[i]]) < 0.44*0.44)) {
	  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
	}
      }
      fclose(f);

    }
  }

  return 0;
}

PetscErrorCode Elmt_Move(IBMNodes *ibm, UserCtx *user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.47;
  PetscReal dz;
  dz = -0.031;
  rcz = rcz-dz;
  rcx = rcx - dz * sin(10./180.*3.1415926);
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  angle =-angle * 3.1415926/180.;
  //angle = 0;
  for (i=0; i<n_v/2; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;

  }
  rcx = -rcx;
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }

  /* Rotate 90 degree */
  for (i=0; i<n_v; i++) {
    temp = ibm->y_bp[i];
    ibm->y_bp[i] = ibm->x_bp[i];
    ibm->x_bp[i] = temp;
  }
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

  if (ti>0) {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt;
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }
  return 0;
}

PetscErrorCode Elmt_Move1(IBMNodes *ibm, UserCtx *user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.47;
  PetscReal dz;
  dz = -0.031;
  rcz = rcz-dz;
  rcx = rcx - dz * sin(10./180.*3.1415926);
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  angle =-angle * 3.1415926/180.;
  //angle = 0;
  for (i=0; i<n_v/2; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;

  }
  rcx = -rcx;
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }

  /* Rotate 90 degree */
  for (i=0; i<n_v; i++) {
    temp = ibm->y_bp[i];
    ibm->y_bp[i] = ibm->x_bp[i];
    ibm->x_bp[i] = temp;
  }
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
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      ibm->uold[i] = ibm->u[i];

      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt;
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }
  return 0;
}

PetscErrorCode fluxin1(UserCtx *user)
{
  PetscReal rg=54;
  PetscInt  iRotate;

  PetscInt ts_p_cycle, systole_steps;
  PetscInt opening, closing;
  PetscInt open_steps, close_steps;

  ts_p_cycle = 1250;
  opening = 10;
  open_steps = 100;

  closing = 580;
  close_steps = 80;

  iRotate = ti - ((ti / ts_p_cycle) * ts_p_cycle);

  PetscReal t_rel;
  t_rel = iRotate * (1. / ts_p_cycle);

  PetscInt i;
  PetscBool interpolated = PETSC_FALSE;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (iRotate >= opening && iRotate<=opening + open_steps) {
    angle = -rg + rg/(PetscReal)open_steps * (iRotate - opening);
  }
  else if (iRotate>(closing - close_steps) && iRotate <=closing) {
    angle = -rg/(PetscReal)close_steps  * (iRotate-(closing-close_steps));
  }
  else if (iRotate>closing) {
    angle = -rg;
  }
  else {
    angle = 0;
  }
/*   angle = 0.; */
  Flux_in = 1.;
  PetscPrintf(PETSC_COMM_WORLD, "Angle %i %le\n", iRotate, angle);
  return 0;
}

//------------------------------------------------------------------------------------
PetscErrorCode shear(UserCtx *user, PetscInt ti){

  PetscInt  bi;

  for (bi=0; bi<block_number; bi++) {
    
    DM            da = user[bi].da, fda = user[bi].fda;
    DMDALocalInfo    info = user[bi].info;
    PetscInt    xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt    zs = info.zs, ze = info.zs + info.zm;
    PetscInt    mx = info.mx, my = info.my, mz = info.mz;

    PetscInt    lxs, lxe, lys, lye, lzs, lze;
    
    PetscInt i, j, k, ii,jj,kk;
    lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
    
    Cmpnts ***ucat, ***csi, ***eta, ***zet;
    PetscReal ***aj, ***q, ***nvert;//***shrate,***sholdarr;//, ***l2;
    //Vec shold;
    
    PetscInt   njac, nrot;
    PetscReal  a[3][3], v[3][3], d[3],x[3], shearmax=0.;
    
    if (lxs==0) lxs++;
    if (lxe==mx) lxe--;
    if (lys==0) lys++;
    if (lye==my) lye--;
    if (lzs==0) lzs++;
    if (lze==mz) lze--;
    
    
    //VecDuplicate(user->shrr, &shold);
    //VecCopy(user->shrr, shold);
    
    
    DMDAVecGetArray(fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(fda, user[bi].Csi, &csi);
    DMDAVecGetArray(fda, user[bi].Eta, &eta);
    DMDAVecGetArray(fda, user[bi].Zet, &zet);
    
    DMDAVecGetArray(da, user[bi].Aj, &aj);
    DMDAVecGetArray(da, user[bi].Nvert, &nvert);
    DMDAVecGetArray(da, user[bi].shrr, &q);
    //DMDAVecGetArray(da, user->rate, &shrate);
    // DMDAVecGetArray(da, shold, &sholdarr);
    /*   DMDAVecGetArray(user[bi].da, user[bi].Phi111, &l2); */
    
    PetscErrorCode imag;
    
    PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;
    
    PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
    
    PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

    PetscReal smax, smin;
    PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
    PetscReal test=0.0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (i==1 || nvert[k][j][i-1] > 2.) {
	    uc = (ucat[k][j][i+1].x - ucat[k][j][i].x);
	    vc = (ucat[k][j][i+1].y - ucat[k][j][i].y);
	    wc = (ucat[k][j][i+1].z - ucat[k][j][i].z);
	  }
	  else if (i==mx-2 || nvert[k][j][i+1] > 2.) {
	    uc = (ucat[k][j][i].x - ucat[k][j][i-1].x);
	    vc = (ucat[k][j][i].y - ucat[k][j][i-1].y);
	    wc = (ucat[k][j][i].z - ucat[k][j][i-1].z);
	  }
	  else {
	    uc = (ucat[k][j][i+1].x - ucat[k][j][i-1].x)*0.5;
	    vc = (ucat[k][j][i+1].y - ucat[k][j][i-1].y)*0.5;
	    wc = (ucat[k][j][i+1].z - ucat[k][j][i-1].z)*0.5;
	  }
	  
	  if (j==1 || nvert[k][j-1][i] > 2) {
	    ue = (ucat[k][j+1][i].x - ucat[k][j][i].x);
	    ve = (ucat[k][j+1][i].y - ucat[k][j][i].y);
	    we = (ucat[k][j+1][i].z - ucat[k][j][i].z);
	  }
	  else if (j==my-2 || nvert[k][j+1][i] > 2) {
	    ue = (ucat[k][j][i].x - ucat[k][j-1][i].x);
	    ve = (ucat[k][j][i].y - ucat[k][j-1][i].y);
	    we = (ucat[k][j][i].z - ucat[k][j-1][i].z);
	  }
	  else {
	    ue = (ucat[k][j+1][i].x - ucat[k][j-1][i].x) * 0.5;
	    ve = (ucat[k][j+1][i].y - ucat[k][j-1][i].y) * 0.5;
	    we = (ucat[k][j+1][i].z - ucat[k][j-1][i].z) * 0.5;
	  }
	  
	  if (k==1 || nvert[k-1][j][i] > 2) {
	    uz = (ucat[k+1][j][i].x - ucat[k][j][i].x);
	    vz = (ucat[k+1][j][i].y - ucat[k][j][i].y);
	    wz = (ucat[k+1][j][i].z - ucat[k][j][i].z);
	  }
	  else if (k==mz-2 || nvert[k+1][j][i] > 2) {
	    uz = (ucat[k][j][i].x - ucat[k-1][j][i].x);
	    vz = (ucat[k][j][i].y - ucat[k-1][j][i].y);
	    wz = (ucat[k][j][i].z - ucat[k-1][j][i].z);
	  }
	  else {
	    uz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x) * 0.5;
	    vz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y) * 0.5;
	    wz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z) * 0.5;
	  }
	  
	  csi1 = (csi[k][j][i].x + csi[k][j][i-1].x) * 0.5 * aj[k][j][i];
	  csi2 = (csi[k][j][i].y + csi[k][j][i-1].y) * 0.5 * aj[k][j][i];
	  csi3 = (csi[k][j][i].z + csi[k][j][i-1].z) * 0.5 * aj[k][j][i];
	  
	  eta1 = (eta[k][j][i].x + eta[k][j-1][i].x) * 0.5 * aj[k][j][i];
	  eta2 = (eta[k][j][i].y + eta[k][j-1][i].y) * 0.5 * aj[k][j][i];
	  eta3 = (eta[k][j][i].z + eta[k][j-1][i].z) * 0.5 * aj[k][j][i];
	  
	  zet1 = (zet[k][j][i].x + zet[k-1][j][i].x) * 0.5 * aj[k][j][i];
	  zet2 = (zet[k][j][i].y + zet[k-1][j][i].y) * 0.5 * aj[k][j][i];
	  zet3 = (zet[k][j][i].z + zet[k-1][j][i].z) * 0.5 * aj[k][j][i];

	  //  derivatives for vorticity
	  d11 = uc * csi1 + ue * eta1 + uz * zet1;
	  d12 = uc * csi2 + ue * eta2 + uz * zet2;
	  d13 = uc * csi3 + ue * eta3 + uz * zet3;
	  
	  d21 = vc * csi1 + ve * eta1 + vz * zet1;
	  d22 = vc * csi2 + ve * eta2 + vz * zet2;
	  d23 = vc * csi3 + ve * eta3 + vz * zet3;
	  
	  d31 = wc * csi1 + we * eta1 + wz * zet1;
	  d32 = wc * csi2 + we * eta2 + wz * zet2;
	  d33 = wc * csi3 + we * eta3 + wz * zet3;
	  
	  //!...Sij, Wij tensors
	  s11 = d11 + d11;
	  s12 = d12 + d21;
	  s13 = d13 + d31;
	  
	  s21 = s12;
	  s22 = d22 + d22;
	  s23 = d23 + d32;
	  
	  s31 = s13;
	  s32 = s23;
	  s33 = d33 + d33;
	  
	  // a=tow_ij=1/Re*Sij
	  a[0][0]=1./user[bi].ren*s11;
	  a[0][1]=1./user[bi].ren*s12;
	  a[0][2]=1./user[bi].ren*s13;
	  a[1][0]=1./user[bi].ren*s21;
	  a[1][1]=1./user[bi].ren*s22;
	  a[1][2]=1./user[bi].ren*s23;
	  a[2][0]=1./user[bi].ren*s31;
	  a[2][1]=1./user[bi].ren*s32;
	  a[2][2]=1./user[bi].ren*s33;
	  
	  
	  smax=sqrt(1.0/3.0)*sqrt(a[0][0]*a[0][0]+a[1][1]*a[1][1]+a[2][2]*a[2][2]-a[0][0]*a[1][1]-a[1][1]*a[2][2]-a[0][0]*a[2][2]+3.0*(a[0][1]*a[0][1]+a[1][2]*a[1][2]+a[0][2]*a[0][2]));

	  if (smax>shearmax)  shearmax = smax;
	  
	  q[k][j][i] = (smax);
	  
	  if(test<fabs(ucat[k][j][i].z)) test=fabs(ucat[k][j][i].z);
	}
      }
    }


    
    /* for (k=lzs; k<lze; k++) { */
    /*     for (j=lys; j<lye; j++) { */
    /*       for (i=lxs; i<lxe; i++) { */
    /* shrate[k][j][i]=fabs((q[k][j][i]-sholdarr[k][j][i])/user->dt); */
    /*       } */
    /*     } */
    /*   } */
    
    PetscPrintf(PETSC_COMM_WORLD,"MAX ucat:%f MAX shear:%f \n",test,shearmax);

    DMDAVecRestoreArray(fda, user[bi].Ucat, &ucat);
    DMDAVecRestoreArray(fda, user[bi].Csi, &csi);
    DMDAVecRestoreArray(fda, user[bi].Eta, &eta);
    DMDAVecRestoreArray(fda, user[bi].Zet, &zet);
    
    DMDAVecRestoreArray(da, user[bi].Aj, &aj);
    DMDAVecRestoreArray(da, user[bi].Nvert, &nvert);
    DMDAVecRestoreArray(da, user[bi].shrr, &q);
    //DMDAVecRestoreArray(da, user->rate, &shrate);
    // DMDAVecRestoreArray(da, user->Phi, &l2);
    //DMDAVecRestoreArray(da, shold, &sholdarr);
    //VecDestroy(&shold);
  }

  return 0;
}

