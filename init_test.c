#include "variables.h"
extern PetscInt block_number, inletprofile, blank;
extern PetscReal L_dim;
extern PetscInt  moveframe,rotateframe;
extern PetscInt  les, poisson;
extern int averaging;
extern int cpu_size;

#include <petsctime.h>

PetscErrorCode InflowWaveFormRead(UserCtx *user);
PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user);
PetscErrorCode FormMetrics(UserCtx *user);
PetscErrorCode MetricsDivergence(UserCtx *user);

PetscErrorCode FormInitialize(UserCtx *user)
{
  user->assignedA = PETSC_FALSE;
  user->multinullspace=PETSC_FALSE;

  VecSet(user->Nvert,0.);
  VecSet(user->lNvert, 0.);
  VecSet(user->Phi, 0.);

  if (user->thislevel==user->mglevels-1) {
  VecSet(user->P,0.);
  VecSet(user->lP,0.);
  VecSet(user->Nvert_o,0.);
  VecSet(user->lNvert_o, 0.);

  VecSet(user->Ucont,0.);
  VecSet(user->lUcont,0.);
  VecSet(user->lNFace, 0.);
  VecSet(user->Ucat,0.);
  VecSet(user->Bcs.Ubcs,0.);
  VecSet(user->DUold, 0.);
  if (moveframe || rotateframe) {
    VecSet(user->Vcont,0.);
    VecSet(user->lVcont,0.);
  }	  
  if (rotateframe) {
    VecSet(user->Wcat,0.);
    VecSet(user->lWcat,0.);
  }
  } 

  user->FluxIntpSum=0.;
  user->FluxInSum=0.;
  user->FluxIntfcSum=0.;
  user->FluxOutSum=0.;

  user->ren = 3000.;
  user->dt = 0.01;
  user->cfl=0.2;
  user->vnn=0.2;

  PetscOptionsGetReal(PETSC_NULL, "-ren", &user->ren, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-dt", &user->dt, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-cfl", &user->cfl, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-vnn", &user->vnn, PETSC_NULL);

  user->st = 1.;//0.038406145;
  
  PetscPrintf(PETSC_COMM_WORLD, "Re %le St %le dt %le \n",user->ren,user->st,user->dt);

  return(0);
}

PetscErrorCode MGDACreate(UserMG *usermg, PetscInt bi)
{
  DMBoundaryType  xperiod,yperiod,zperiod;

  MGCtx *mgctx = usermg->mgctx;

  UserCtx *user, *user_high;

  PetscErrorCode ierr;

  PetscInt l;
  

  for (l=usermg->mglevels-2; l>=0; l--) {
    user = mgctx[l].user;
    user_high = mgctx[l+1].user;
         
    if (user[bi].isc) {	 
      user[bi].IM = user_high[bi].IM;
    }
    else {
      user[bi].IM = (user_high[bi].IM + 1) / 2;
    }
    if (user[bi].jsc) {
      user[bi].JM = user_high[bi].JM;
    }
    else {
      user[bi].JM = (user_high[bi].JM + 1) / 2;
    }

    if (user[bi].ksc) {
      user[bi].KM = user_high[bi].KM;
    }
    else {
      user[bi].KM = (user_high[bi].KM + 1) / 2;
    }
    PetscPrintf(PETSC_COMM_WORLD, "Grid dimension in i,j and k direction for level %d are %d %d %d :\n",l,user[bi].IM, user[bi].JM, user[bi].KM);
    if (user[bi].IM*(2- (user[bi].isc))-(user_high[bi].IM+1-(user[bi].isc)) +
	user[bi].JM*(2- (user[bi].jsc))-(user_high[bi].JM+1-(user[bi].jsc)) +
	user[bi].KM*(2- (user[bi].ksc))-(user_high[bi].KM+1-(user[bi].ksc))) {
      PetscPrintf(PETSC_COMM_WORLD, "Grid at level %d can't be further restricted!\n", l);
    }
  }

  l = 0;
  user = mgctx[l].user;

  PetscInt M,N,P,m,n,p,s,MM,NN,PP, *lx, *ly, *lz;
  PetscInt *lxSum, *lySum, *lzSum;
  //Mohsen Jan 12 
  
  s=3;
 
  if(user[bi].bctype[0]==7  && user[bi].bctype[2]==7  && user[bi].bctype[4]==7) {
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  
  else if(user[bi].bctype[0]==7  && user[bi].bctype[2]==7) {
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_NONE;
  }
  else if(user[bi].bctype[2]==7  && user[bi].bctype[4]==7) {
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  else if(user[bi].bctype[0]==7  && user[bi].bctype[4]==7) {
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_NONE;
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  else if(user[bi].bctype[0]==7){
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_NONE;
    zperiod = DM_BOUNDARY_NONE;
  }
  else if(user[bi].bctype[2]==7) {
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_NONE;
  }
  else if(user[bi].bctype[4]==7 ) {
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_NONE; 
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  else {
    s=2;
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_NONE;
    zperiod = DM_BOUNDARY_NONE;
  }
 
  M=user[bi].IM+1;
  N=user[bi].JM+1;
  P=user[bi].KM+1;

  ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,M,N,P,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,s,
		  PETSC_NULL, PETSC_NULL, PETSC_NULL,&(user[bi].da)); //Mohsen Jan 12  
  PetscPrintf(PETSC_COMM_WORLD, "DA is Created for coarsest level and global dimension in i,j and k direction are %d %d %d :\n",M,N,P); 
 
  if (ierr) 
    SETERRQ1(PETSC_COMM_SELF,1, "problem creating DA %d",ierr);
  user[bi].aotopetsc = PETSC_FALSE;
  DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  DMGetCoordinateDM(user[bi].da, &(user[bi].fda));
  DMDAGetLocalInfo(user[bi].da, &(user[bi].info));
  
  if (rans && usermg->mglevels-1==0) { // fda2 for k-omega at the finest level
    DMDAGetInfo(user[bi].da, PETSC_NULL, &MM,&NN,&PP,&m, &n, &p, PETSC_NULL, PETSC_NULL,PETSC_NULL, PETSC_NULL,PETSC_NULL,PETSC_NULL);
 
    ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,user[bi].IM+1, user[bi].JM+1, user[bi].KM+1,
		    m,n,p,2,s, PETSC_NULL, PETSC_NULL, PETSC_NULL,&(user[bi].fda2));
    if (ierr) SETERRQ1(PETSC_COMM_SELF,1, "problem creating FDA2 for RANS%d",ierr);
  }
  
  for (l=1; l<usermg->mglevels; l++) {
    user = mgctx[l].user;
    
    // Get info about the coarse DA
    DMDAGetInfo(mgctx[l-1].user[bi].da, PETSC_NULL, &MM,&NN,&PP,&m,&n,&p,PETSC_NULL,PETSC_NULL,PETSC_NULL, PETSC_NULL,PETSC_NULL,PETSC_NULL);
    if (l==1) PetscPrintf(PETSC_COMM_WORLD, "Corresponing numberof procs in i,j and k direction are %d %d %d :\n",m,n,p);
    // Find the distribution of the coarse DA
    PetscMalloc(m*sizeof(PetscInt), &lx);
    PetscMalloc(n*sizeof(PetscInt), &ly);
    PetscMalloc(p*sizeof(PetscInt), &lz);
    PetscMalloc(m*sizeof(PetscInt), &lxSum);
    PetscMalloc(n*sizeof(PetscInt), &lySum);
    PetscMalloc(p*sizeof(PetscInt), &lzSum);

    PetscInt rank;
    for (rank=0; rank<m; rank++) lx[rank]=0;
    for (rank=0; rank<n; rank++) ly[rank]=0;
    for (rank=0; rank<p; rank++) lz[rank]=0;

    DMDALocalInfo	info = mgctx[l-1].user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt    ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    // Calculationg the distribution of the fine DA 
    // such that the refined portion of the coarse DA is on the
    // same processor as the coarse DA
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if ( user[bi].isc) 
      lx[rank % m]= (xe-xs);
    else {
      if (m==1)
	lx[rank % m]= user[bi].IM+1;
      else if (xs==0) 
	lx[rank % m]= 2*xe-1;
      else if (xe==mx)
	lx[rank % m]= user[bi].IM+1-2*xs+1;
      else
	lx[rank % m]= (xe-xs)*2;
    }
    if ( user[bi].jsc) 
      ly[(rank % (m*n))/m]= (ye-ys);
    else {
      if (n==1)
	ly[(rank % (m*n))/m]= user[bi].JM+1;
      else if (ys==0)
	ly[(rank % (m*n))/m]= 2*ye-1;
      else if (ye==my) 
	ly[(rank % (m*n))/m]= user[bi].JM+1-2*ys+1;
      else
	ly[(rank % (m*n))/m]= (ye-ys)*2;
    }
    if ( user[bi].ksc) 
      lz[rank/(m*n)]=ze-zs;
    else{
      if (p==1)
	lz[rank/(m*n)]= user[bi].KM+1;
      else if (zs==0)
	lz[rank/(m*n)]= 2*ze-1;
      else if (ze==mz)
	lz[rank/(m*n)]= user[bi].KM+1-2*zs+1;
      else
	lz[rank/(m*n)]= (ze-zs)*2;
    }
    
    // All_reduce to all processes
    MM=0;
    for (rank=0; rank<m; rank++) {
      MPI_Allreduce(&lx[rank],&lxSum[rank],1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);
      MM+=lxSum[rank];
    }
    NN=0;
    for (rank=0; rank<n; rank++) {
      MPI_Allreduce(&ly[rank],&lySum[rank],1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);      
      NN+=lySum[rank];
    }
    PP=0;
    for (rank=0; rank<p; rank++) {
      MPI_Allreduce(&lz[rank],&lzSum[rank],1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);      
      PP+=lzSum[rank];
    }
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscBarrier(PETSC_NULL);
    
    // Create the refined DA based on the new distribution
   
    ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,user[bi].IM+1,user[bi].JM+1,user[bi].KM+1,m,n,p,1,s,lxSum,lySum,lzSum,&(user[bi].da));
    PetscPrintf(PETSC_COMM_WORLD, "Global dimension of DA in i,j and k direction for level %d are %d %d %d :\n",l,user[bi].IM+1,user[bi].JM+1,user[bi].KM+1);
   
    if (ierr) SETERRQ1(PETSC_COMM_SELF,1, "problem creating DA %d",ierr);
    
    if (rans &&  l == usermg->mglevels-1) { // fda2 for k-omega at the finest level
      ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,
		      user[bi].IM+1, user[bi].JM+1, user[bi].KM+1,
		      m, n, p,
		      2, s, lxSum, lySum, lzSum,
		      &(user[bi].fda2));
      if (ierr) SETERRQ1(PETSC_COMM_SELF,1, "problem creating FDA2 for RANS%d",ierr);
    }
    user[bi].aotopetsc = PETSC_FALSE;
    DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DMGetCoordinateDM(user[bi].da, &(user[bi].fda));
    DMDAGetLocalInfo(user[bi].da, &(user[bi].info));
    
    PetscFree(lx);PetscFree(ly);PetscFree(lz);  
    PetscFree(lxSum);PetscFree(lySum);PetscFree(lzSum);  
  }
  return 0;
}

PetscErrorCode MG_Initial(UserMG *usermg, IBMNodes *ibm)
{
  MGCtx *mgctx;

  PetscInt level;

  PetscInt rank;
  PetscInt IM, JM, KM;

  PetscInt i, j, k;
  PetscInt ih, jh, kh;
  PetscInt bi,sb;

  PetscErrorCode ierr; 

  Vec Coor, gCoor, Coor_high;
  Cmpnts ***coor, ***coor_high;
  PetscReal *gc;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  PetscReal cl = 1.;
  
  PetscOptionsGetReal(PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);

  /* How many MG levels, the default is 3 */
  if (poisson) 
    usermg->mglevels = 1;
  else
  usermg->mglevels = 3;
  PetscOptionsGetInt(PETSC_NULL, "-mg_level", &usermg->mglevels, PETSC_NULL);
  
  usermg->ksc = PETSC_FALSE;
  usermg->jsc = PETSC_FALSE;
  usermg->isc = PETSC_FALSE;
  
/*   PetscOptionsGetTruth(PETSC_NULL, "-mg_k_semi", &usermg->ksc, PETSC_NULL); */
/*   PetscOptionsGetTruth(PETSC_NULL, "-mg_j_semi", &usermg->jsc, PETSC_NULL); */
/*   PetscOptionsGetTruth(PETSC_NULL, "-mg_i_semi", &usermg->isc, PETSC_NULL); */
  

  PetscMalloc(usermg->mglevels*sizeof(MGCtx), &(usermg->mgctx));
  mgctx = usermg->mgctx;
  FILE *fd, *fd1;
  
  /* Read in number of blocks and allocate memory for UserCtx */
  PetscInt generate_grid=0, grid1d=0;
  PetscOptionsGetInt(PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);
  if (generate_grid) {
    block_number=1;
    PetscOptionsGetInt(PETSC_NULL, "-block_number", &block_number, PETSC_NULL);
  } else {
    if (!rank) {
    fd = fopen("grid.dat", "r");
    
    fscanf(fd, "%i\n", &block_number);
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
/*     fclose(fd); */
    }
    else {
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }
  
  PetscInt ksc[block_number], jsc[block_number], isc[block_number], nblk=block_number;
  PetscOptionsGetIntArray(PETSC_NULL, "-mg_k_semi", ksc, &nblk, PETSC_NULL);
  PetscOptionsGetIntArray(PETSC_NULL, "-mg_j_semi", jsc, &nblk, PETSC_NULL);
  PetscOptionsGetIntArray(PETSC_NULL, "-mg_i_semi", isc, &nblk, PETSC_NULL);
 
  PetscInt cgrid[block_number];
  
  PetscOptionsGetIntArray(PETSC_NULL, "-cgrid",cgrid, &nblk, PETSC_NULL);
 
  for (level=0; level<usermg->mglevels; level++) {
    PetscMalloc(block_number*sizeof(UserCtx), &mgctx[level].user);
    for (bi=0; bi<block_number; bi++) {
      mgctx[level].user[bi].ibm = ibm;
      mgctx[level].user[bi].isc = isc[bi];//&usermg->isc;
      mgctx[level].user[bi].jsc = jsc[bi];//&usermg->jsc;
      mgctx[level].user[bi].ksc = ksc[bi];//&usermg->ksc;
      mgctx[level].user[bi].cgrid = cgrid[bi];
      if (level==0) PetscPrintf(PETSC_COMM_WORLD, "SEMI_COARSENING for block %d in i=%d j=%d k=%d\n", bi, isc[bi], jsc[bi], ksc[bi]);
      if (level==0) PetscPrintf(PETSC_COMM_WORLD, "C-grid for block %d is %d\n", bi,cgrid[bi]);
    }
  }
  
  /* Read BC type and interface control file at the 
     finest grid level*/
  level = usermg->mglevels-1;
  UserCtx *user;

  user = mgctx[level].user;
  if (!rank) {
    
    if (block_number>1) {
      fd1 = fopen("interface.dat", "r"); if (!fd1) SETERRQ(PETSC_COMM_SELF,1, "Cannot open interface.dat to write!!");
      for (bi=0; bi<block_number; bi++) {
	if (blank && bi==0) {    
	  for (sb=1; sb<block_number; sb++) { 
	    fscanf(fd1, "%i %i %i %i %i %i %i\n", &user[bi].ip[sb],&user[bi].jp[sb],&user[bi].kp[sb],
		   &user[bi].dispx[sb],&user[bi].dispy[sb],&user[bi].dispz[sb],&user[bi].dispnn[sb]);
	    PetscPrintf(PETSC_COMM_WORLD,"BLANK ip,jp,kp,dispy,dispz,dispnn - %i,%i,%i,%i,%i,%i\n", user[bi].ip[sb], user[bi].jp[sb], 
			user[bi].kp[sb], user[bi].dispy[sb], user[bi].dispz[sb], user[bi].dispnn[sb]);
	  }
	  MPI_Bcast(&user[bi].ip, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].jp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].kp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispx, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispy, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispz, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispnn, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  
	}   
	fscanf(fd1, "%i\n", &(user[bi].itfcptsnumber));
	MPI_Bcast(&(user[bi].itfcptsnumber), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "itfcp number # %d\n",user[bi].itfcptsnumber);

	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcK));
	
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostK));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostB));
	
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostx));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchosty));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostz));
      
	for (i=0; i<user[bi].itfcptsnumber; i++) {
	  fscanf(fd1, "%i %i %i\n", &(user[bi].itfcI[i]), &(user[bi].itfcJ[i]),
		 &(user[bi].itfcK[i]));
	  fscanf(fd1, "%i %i %i %i\n", &(user[bi].itfchostI[i]),
		 &(user[bi].itfchostJ[i]), &(user[bi].itfchostK[i]),
		 &(user[bi].itfchostB[i]));
	  fscanf(fd1, "%le %le %le\n", &(user[bi].itfchostx[i]),
		 &(user[bi].itfchosty[i]), &(user[bi].itfchostz[i]));
	}

	MPI_Bcast(user[bi].itfcI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostB, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostx, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchosty, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostz, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
      }
      fclose(fd1);
      
      /* Read the interface control file */
      fd1 = fopen("intcontrol.dat", "r"); if (!fd1) SETERRQ(PETSC_COMM_SELF,1, "Cannot open intcontrol.dat to write!!");
      for (bi=0; bi<block_number; bi++) {
	fscanf(fd1, "%i %i %i %i %i %i %i\n", &(user[bi].inttype[0]),
	       &(user[bi].inttype[1]), &(user[bi].inttype[2]),
	       &(user[bi].inttype[3]), &(user[bi].inttype[4]),
	       &(user[bi].inttype[5]), &(user[bi].inttype[6]));
	MPI_Bcast(&(user[bi].inttype[0]), 7, MPI_INT, 0, PETSC_COMM_WORLD);
      }
      fclose(fd1);
    } // if block_number>1
    
    /* Read in bcs.dat for boundary conditions at 6 boundary surfaces 
     First put the data onto the finest level and restrict to the coarser
     levels */
    
    fd1 = fopen("bcs.dat", "r"); if (!fd1) SETERRQ(PETSC_COMM_SELF,1, "Cannot open bcs.dat to write!!");

    for (bi=0; bi<block_number; bi++) {
      fscanf(fd1, "%i %i %i %i %i %i\n", &(user[bi].bctype[0]),
	     &(user[bi].bctype[1]), &(user[bi].bctype[2]),
	     &(user[bi].bctype[3]), &(user[bi].bctype[4]),
	     &(user[bi].bctype[5]));
      if (user[bi].bctype[0]==7 || user[bi].bctype[1]==7) i_periodic=1;
      if (user[bi].bctype[2]==7 || user[bi].bctype[3]==7) j_periodic=1;
      if (user[bi].bctype[4]==7 || user[bi].bctype[5]==7) k_periodic=1;

      MPI_Bcast(&(user[bi].bctype[0]), 6, MPI_INT, 0, PETSC_COMM_WORLD);
     
    }
    fclose(fd1);

  }  else { // if !rank
    if (block_number>1) {  
      for (bi=0; bi<block_number; bi++) {
  	if (blank && bi==0) {     
	  MPI_Bcast(&user[bi].ip, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].jp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].kp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispx, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispy, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispz, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispnn, 10, MPI_INT, 0, PETSC_COMM_WORLD);	 
	}    

	MPI_Bcast(&(user[bi].itfcptsnumber), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcK));
      
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostK));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostB));
      
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostx));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchosty));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostz));

      
	MPI_Bcast(user[bi].itfcI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostB, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostx, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchosty, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostz, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
      }

      for (bi=0; bi<block_number; bi++) {
	MPI_Bcast(&(user[bi].inttype[0]), 7, MPI_INT, 0, PETSC_COMM_WORLD);
      }
    }

    for (bi=0; bi<block_number; bi++) {
      MPI_Bcast(&(user[bi].bctype[0]), 6, MPI_INT, 0, PETSC_COMM_WORLD);
      if (user[bi].bctype[0]==7 || user[bi].bctype[1]==7) i_periodic=1;
      if (user[bi].bctype[2]==7 || user[bi].bctype[3]==7) j_periodic=1;
      if (user[bi].bctype[4]==7 || user[bi].bctype[5]==7) k_periodic=1;
    }
  }
  
  PetscPrintf(PETSC_COMM_WORLD, "i_periodic is %d  , j_periodic is %d and k_periodic is %d \n", i_periodic, j_periodic, k_periodic);
 
  /* The bcs.dat for boundary conditions at 6 boundary surfaces 
     on the finest level is restricted to the coarser
     levels */

  for (level = usermg->mglevels-2; level>=0; level--) {
      for (bi=0; bi<block_number; bi++) {
      mgctx[level].user[bi].bctype[0] = mgctx[level+1].user[bi].bctype[0];
      mgctx[level].user[bi].bctype[1] = mgctx[level+1].user[bi].bctype[1];
      mgctx[level].user[bi].bctype[2] = mgctx[level+1].user[bi].bctype[2];
      mgctx[level].user[bi].bctype[3] = mgctx[level+1].user[bi].bctype[3];
      mgctx[level].user[bi].bctype[4] = mgctx[level+1].user[bi].bctype[4];
      mgctx[level].user[bi].bctype[5] = mgctx[level+1].user[bi].bctype[5];
    }
  }

  /* Read from grid.dat the number of grid points along I, J, K directions
     and the detail coordinates of grid nodes */
  level = usermg->mglevels-1;
  user = mgctx[level].user;

  if (inletprofile==3) {
  PetscPrintf(PETSC_COMM_WORLD, "READ INFLOW WAVE FORM!!!! mg_init\n");
  for (bi=0; bi<block_number; bi++) {
    InflowWaveFormRead(&user[bi]);
  }
  PetscPrintf(PETSC_COMM_WORLD, "READ INFLOW WAVE FORM!!!! mg_init\n");
  }

  PetscInt imm[block_number], kmm[block_number], jmm[block_number];
  if (generate_grid) {
    PetscOptionsGetIntArray(PETSC_NULL, "-im", imm, &nblk, PETSC_NULL);
    PetscOptionsGetIntArray(PETSC_NULL, "-jm", jmm, &nblk, PETSC_NULL);
    PetscOptionsGetIntArray(PETSC_NULL, "-km", kmm, &nblk, PETSC_NULL);
  }
  
  for (bi=0; bi<block_number; bi++) {

    if (!rank) {
      if (!generate_grid)
	fscanf(fd, "%i %i %i\n", &(user[bi].IM),
	     &(user[bi].JM), &(user[bi].KM));
      else {
	user[bi].IM=imm[bi];
	user[bi].JM=jmm[bi];
	user[bi].KM=kmm[bi];
      }	
      IM = user[bi].IM;
      JM = user[bi].JM;
      KM = user[bi].KM;
      
      MPI_Bcast(&(IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&(IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      user[bi].IM = IM;
      user[bi].JM = JM;
      user[bi].KM = KM;
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "READ GRID.dat %le\n", L_dim);

    MGDACreate(usermg, bi);

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    //    PetscInt	mx = info.mx, my = info.my, mz = info.mz;

    if (grid1d) PetscMalloc((IM+JM+KM)*sizeof(PetscReal), &gc);
    else        PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);
  
    //DMGetGhostedCoordinates(user[bi].da, &Coor);
    DMGetCoordinatesLocal(user[bi].da, &Coor);
    VecSet(Coor,0.0);//Mohsen Feb 12//
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
	PetscReal Lx=0.1,Ly=10,Lz=29.5,X0=0.1,Y0=-0.5,Z0=0.5;	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3) = (Lx)/(IM-1.) * i+X0;	
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3+1) = (Ly)/(JM-1.) * j+Y0;
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3+2) = (Lz)/(KM-1.) * k+Z0 ;
	    }
	  }
	}
	
	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl*L_dim;
		coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl*L_dim;
		coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl*L_dim;
	      }
	    }
	  }
	}
	
      } // grid-1d
    } // rank
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
		coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl*L_dim;
		coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl*L_dim;
		coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl*L_dim;
	      }
	    }
	  }
	}
      } // grid-1d
    }
    PetscFree(gc);
    DMDAVecRestoreArray(user[bi].fda, Coor, &coor);
    

    DMGetCoordinates(user[bi].da, &gCoor);
    DMLocalToGlobalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    DMLocalToGlobalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    
    DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    
    // VecDestroy(&Coor);
    // VecDestroy(&gCoor);
  } //bi

  UserCtx *user_high;
  for (level=usermg->mglevels-2; level>-1; level--) {
    user = mgctx[level].user;
    user_high = mgctx[level+1].user;
    for (bi = 0; bi<block_number; bi++) {
      
      DMDALocalInfo	info = user[bi].info;
      PetscInt	xs = info.xs, xe = info.xs + info.xm;
      PetscInt  ys = info.ys, ye = info.ys + info.ym;
      PetscInt	zs = info.zs, ze = info.zs + info.zm;
      PetscInt	mx = info.mx, my = info.my, mz = info.mz;
      
      //DMGetGhostedCoordinates(user_high[bi].da, &Coor_high);
      DMGetCoordinatesLocal(user_high[bi].da, &Coor_high);
      DMGetCoordinates(user[bi].da, &Coor);

      DMDAVecGetArray(user_high[bi].fda,Coor_high, &coor_high);
      DMDAVecGetArray(user[bi].fda, Coor, &coor);

      
      if (xe==mx) xe--;
      if (ye==my) ye--;
      if (ze==mz) ze--;
     
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {

	    GridRestriction(i, j, k, &ih, &jh, &kh, &user[bi]);
	    coor[k][j][i].x = coor_high[kh][jh][ih].x;
	    coor[k][j][i].y = coor_high[kh][jh][ih].y;
	    coor[k][j][i].z = coor_high[kh][jh][ih].z;
	  }
	}
      }
      
      DMDAVecRestoreArray(user[bi].fda, Coor, &coor);
      DMDAVecRestoreArray(user_high[bi].fda, Coor_high, &coor_high);

      // DMGetGhostedCoordinates(user[bi].da, &gCoor);
      DMGetCoordinatesLocal(user[bi].da, &gCoor);

      DMGlobalToLocalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
      DMGlobalToLocalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);

      // VecDestroy(&Coor);
      // VecDestroy(&Coor_high);
      //  VecDestroy(&gCoor);
      
    }
  }

  if (!rank && !generate_grid) fclose(fd);


  /* Create the global and local vectors at all grid levels */

  for (level=usermg->mglevels-1; level>=0; level--) {
    user = mgctx[level].user;
    
    for (bi=0; bi<block_number; bi++) {

      user[bi].thislevel = level;
      user[bi]._this = bi;
      user[bi].mglevels = usermg->mglevels;
      if (level > 0) {
	user[bi].da_c = &mgctx[level-1].user[bi].da;
	user[bi].lNvert_c = &mgctx[level-1].user[bi].lNvert;
	user[bi].user_c = &mgctx[level-1].user[bi];
      }
      if (level < usermg->mglevels-1) {
	user[bi].da_f = &mgctx[level+1].user[bi].da;
	user[bi].user_f = &mgctx[level+1].user[bi];
      }
     
      //Destroyed in metrics.c from here !!!1
      ierr = DMCreateGlobalVector(user[bi].fda, &(user[bi].Csi));
      ierr = VecDuplicate(user[bi].Csi, &(user[bi].Eta));
      ierr = VecDuplicate(user[bi].Csi, &(user[bi].Zet));

      VecDuplicate(user[bi].Csi, &(user[bi].ICsi));
      VecDuplicate(user[bi].Csi, &(user[bi].IEta));
      VecDuplicate(user[bi].Csi, &(user[bi].IZet));
      VecDuplicate(user[bi].Csi, &(user[bi].JCsi));
      VecDuplicate(user[bi].Csi, &(user[bi].JEta));
      VecDuplicate(user[bi].Csi, &(user[bi].JZet));
      VecDuplicate(user[bi].Csi, &(user[bi].KCsi));
      VecDuplicate(user[bi].Csi, &(user[bi].KEta));
      VecDuplicate(user[bi].Csi, &(user[bi].KZet));
      //Destroyed in metrics.c to here!!!

      VecDuplicate(user[bi].Csi, &(user[bi].Cent));
      
      VecDuplicate(user[bi].Csi, &(user[bi].GridSpace));

      if(level==usermg->mglevels-1) {		     				     
	VecDuplicate(user[bi].Csi, &(user[bi].Ucont));
	VecDuplicate(user[bi].Csi, &(user[bi].Ucont_o));
	VecDuplicate(user[bi].Csi, &(user[bi].Ucont_rm1));
	VecDuplicate(user[bi].Csi, &(user[bi].Ucat));
	VecDuplicate(user[bi].Csi, &(user[bi].Ucat_o));
	VecDuplicate(user[bi].Csi, &(user[bi].DUold));
	VecDuplicate(user[bi].Csi, &(user[bi].Bcs.Ubcs));
	VecDuplicate(user[bi].Csi, &(user[bi].psuedot));
	
	VecDuplicate(user[bi].Csi, &(user[bi].Vcont));
	if (rotateframe)
	  VecDuplicate(user[bi].Csi, &(user[bi].Wcat));

	VecDuplicate(user[bi].Csi, &(user[bi].Itfc));
	///////////////////////////////////////////
	VecDuplicate(user[bi].Csi, &(user[bi].ItfcQ));
	///////////////////////////////////////////
      }
      //      VecDuplicate(user[bi].Csi, &(user[bi].Rhs));
      if (level < usermg->mglevels-1) {
	VecDuplicate(user[bi].Csi, &(user[bi].Forcing));
	VecDuplicate(user[bi].Csi, &(user[bi].Ucont_MG));
      }

      ierr = DMCreateGlobalVector(user[bi].da, &(user[bi].Aj)); CHKERRQ(ierr); //Destroyed in metrics.c
      VecDuplicate(user[bi].Aj, &(user[bi].IAj));//Destroyed in metrics.c
      VecDuplicate(user[bi].Aj, &(user[bi].JAj));//Destroyed in metrics.c
      VecDuplicate(user[bi].Aj, &(user[bi].KAj));//Destroyed in metrics.c

      VecDuplicate(user[bi].Aj, &(user[bi].Nvert));
      VecDuplicate(user[bi].Aj, &(user[bi].Phi));

      if(level==usermg->mglevels-1) {
      VecDuplicate(user[bi].Aj, &(user[bi].P));
      //////////////////////////////////////////
      VecDuplicate(user[bi].P, &(user[bi].Qnew));
      VecDuplicate(user[bi].P, &(user[bi].shrr));
      //////////////////////////////////////////
      VecDuplicate(user[bi].Aj, &(user[bi].P_o));
      VecDuplicate(user[bi].Aj, &(user[bi].Nvert_o));
      }

      if(averaging && level==usermg->mglevels-1) {// New for Averaging
	VecDuplicate(user[bi].Csi, &user[bi].Ucat_sum);
	VecDuplicate(user[bi].Csi, &user[bi].Ucat_cross_sum);
	VecDuplicate(user[bi].Csi, &user[bi].Ucat_square_sum);
	VecSet(user[bi].Ucat_sum,0);
	VecSet(user[bi].Ucat_cross_sum,0);
	VecSet(user[bi].Ucat_square_sum,0);
	  
	VecDuplicate(user[bi].Aj, &user[bi].P_sum);
	VecSet(user[bi].P_sum,0);
      }
    
      //      VecDuplicate(user[bi].Aj, &(user[bi].Volume));

      DMCreateLocalVector(user[bi].fda, &(user[bi].lCsi));

      VecDuplicate(user[bi].lCsi, &(user[bi].lEta));
      VecDuplicate(user[bi].lCsi, &(user[bi].lZet));
      VecDuplicate(user[bi].lCsi, &(user[bi].lICsi));
      VecDuplicate(user[bi].lCsi, &(user[bi].lIEta));
      VecDuplicate(user[bi].lCsi, &(user[bi].lIZet));
      VecDuplicate(user[bi].lCsi, &(user[bi].lJCsi));
      VecDuplicate(user[bi].lCsi, &(user[bi].lJEta));
      VecDuplicate(user[bi].lCsi, &(user[bi].lJZet));
      VecDuplicate(user[bi].lCsi, &(user[bi].lKCsi));
      VecDuplicate(user[bi].lCsi, &(user[bi].lKEta));
      VecDuplicate(user[bi].lCsi, &(user[bi].lKZet));
      VecDuplicate(user[bi].lCsi, &(user[bi].lGridSpace));
      VecDuplicate(user[bi].lCsi, &(user[bi].lCent));
      VecDuplicate(user[bi].lCsi, &(user[bi].Centx));
      VecDuplicate(user[bi].lCsi, &(user[bi].Centy));
      VecDuplicate(user[bi].lCsi, &(user[bi].Centz));	
      if(level==usermg->mglevels-1) {			 
      VecDuplicate(user[bi].lCsi, &(user[bi].lUcont));
      VecDuplicate(user[bi].lUcont, &(user[bi].lNFace));
      VecDuplicate(user[bi].lCsi, &(user[bi].lUcat));
      VecDuplicate(user[bi].lCsi, &(user[bi].lItfc));
      VecDuplicate(user[bi].lCsi, &(user[bi].lVcont));
      if (rotateframe) {
      VecDuplicate(user[bi].lCsi, &(user[bi].lWcat));
      VecDuplicate(user[bi].lCsi, &(user[bi].lMAreaCsi));
      VecDuplicate(user[bi].lCsi, &(user[bi].lMAreaEta));
      VecDuplicate(user[bi].lCsi, &(user[bi].lMAreaZet));
      //   FormAreaMoment(&user[bi]);
      }
      VecDuplicate(user[bi].lCsi, &(user[bi].lUcont_o));
      VecDuplicate(user[bi].lCsi, &(user[bi].lUcont_rm1));
      //      VecDuplicate(user[bi].lCsi, &(user[bi].lArea));
      }
      DMCreateLocalVector(user[bi].da, &(user[bi].lAj));
      VecDuplicate(user[bi].lAj, &(user[bi].lIAj));
      VecDuplicate(user[bi].lAj, &(user[bi].lJAj));
      VecDuplicate(user[bi].lAj, &(user[bi].lKAj));
      //    VecDuplicate(user[bi].lAj, &(user[bi].lItfcP));

      VecDuplicate(user[bi].lAj, &(user[bi].lNvert));
        VecDuplicate(user[bi].lAj, &(user[bi].lPhi));
 
      if(level==usermg->mglevels-1) {
      VecDuplicate(user[bi].lAj, &(user[bi].lP));
      //////////////////////////////////////////
      VecDuplicate(user[bi].lP, &(user[bi].Ql));
      VecDuplicate(user[bi].lP, &(user[bi].lItfcQ));
      //////////////////////////////////////////
      VecDuplicate(user[bi].lAj, &(user[bi].lNvert_o));
      }
      //      VecDuplicate(user[bi].lAj, &(user[bi].lVolume));

      if(rans  && level==usermg->mglevels-1) {
	DMCreateLocalVector(user[bi].fda2, &user[bi].lK_Omega);	
	VecSet(user[bi].lK_Omega, 0);	// seokkoo
	DMCreateLocalVector(user[bi].fda2, &user[bi].lK_Omega_o);	
	VecSet(user[bi].lK_Omega_o, 0);	// seokkoo	
	VecDuplicate(user[bi].lP, &user[bi].Distance);
	VecSet(user[bi].Distance, 0);

	DMCreateGlobalVector(user[bi].fda2, &user[bi].K_Omega);	
	VecSet(user[bi].K_Omega, 0);	// seokkoo
	VecDuplicate(user[bi].K_Omega, &(user[bi].K_Omega_o));	
	VecSet(user[bi].K_Omega_o, 0);// seokkoo
	//VecDuplicate(user[bi].K_Omega, &(user[bi].K_Omega_rm1));
	//VecDuplicate(user[bi].lP, &(user[bi].lSrans));		VecSet(user[bi].lSrans, 0);// seokkoo
	VecDuplicate(user[bi].lP, &(user[bi].lNu_t));		
	VecSet(user[bi].lNu_t, 0);// seokkoo

	if(rans==3) {
	  VecDuplicate(user[bi].lP, &(user[bi].lF1));
	  VecSet(user[bi].lF1, 0);
	}
      }
      
      // New for LES
      if(les && level==usermg->mglevels-1) {
	VecDuplicate(user[bi].lAj, &user[bi].lCs);	
	VecDuplicate(user[bi].lAj, &(user[bi].lNu_t));		
	VecSet(user[bi].lNu_t, 0);
      }
      if ((wallfunction || rans) && level==usermg->mglevels-1) 
	VecDuplicate(user[bi].lAj, &user[bi].lUstar);

      PetscPrintf(PETSC_COMM_WORLD,"global and local vectors at all level %d  are created\n",level);

      FormInitialize(&(user[bi]));

      PetscPrintf(PETSC_COMM_WORLD,"Initialization for level %d is done \n",level);

      // read Nvert for blanking in overset grids
      if (block_number>1 && bi==0 && blank && level==usermg->mglevels-1) {
	PetscViewer	viewer;
	char filen[80];
	sprintf(filen, "nvfield_blank.dat");
	
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad((user[bi].Nvert_o),viewer);
	PetscViewerDestroy(&viewer);
	
	DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
	DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
	
	VecCopy(user[bi].Nvert_o, user[bi].Nvert);
	VecCopy(user[bi].lNvert_o, user[bi].lNvert);    
      }
	
      FormMetrics(&(user[bi]));  
      PetscPrintf(PETSC_COMM_WORLD,"Metrics formation for level %d is done \n",level);

      if (level==usermg->mglevels-1) MetricsDivergence(&(user[bi]));
    }

  }

  return(0);
}


PetscErrorCode MG_Finalize(UserMG *usermg)
{

  MGCtx *mgctx;

  PetscInt level, bi;

  UserCtx *user;  

  mgctx = usermg->mgctx;
  for (level=usermg->mglevels-1; level>=0; level--) {
    user=mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {

      if(level==usermg->mglevels-1 && averaging) {
	VecDestroy(&user[bi].Ucat_square_sum);
	VecDestroy(&user[bi].Ucat_cross_sum);
	VecDestroy(&user[bi].Ucat_sum);
      }

      if ((wallfunction || rans) && level==usermg->mglevels-1) 
	VecDestroy(&user[bi].lUstar);

      VecDestroy(&user[bi].Cent);
      VecDestroy(&user[bi].Centx);
      VecDestroy(&user[bi].Centy);
      VecDestroy(&user[bi].Centz);
      VecDestroy(&user[bi].GridSpace);

      if(level==usermg->mglevels-1) {
      VecDestroy(&user[bi].Ucont);
      VecDestroy(&user[bi].Ucont_o);
      VecDestroy(&user[bi].Ucont_rm1);
      VecDestroy(&user[bi].Ucat);
      VecDestroy(&user[bi].Ucat_o);
      VecDestroy(&user[bi].DUold);
      VecDestroy(&user[bi].Bcs.Ubcs);
      VecDestroy(&user[bi].psuedot);

      VecDestroy(&user[bi].Vcont);
      if (rotateframe)
      VecDestroy(&user[bi].Wcat);

      VecDestroy(&user[bi].Itfc);
      //////////////////////////
      VecDestroy(&user[bi].ItfcQ);
      //////////////////////////
      }
     
      if (level < usermg->mglevels-1) {
	VecDestroy(&user[bi].Forcing);
	VecDestroy(&user[bi].Ucont_MG);
      }

      VecDestroy(&user[bi].Nvert);
      VecDestroy(&user[bi].Phi);

      if(level==usermg->mglevels-1) {
      VecDestroy(&user[bi].P);
      ///////////////////////////////
      VecDestroy(&user[bi].Qnew);
      VecDestroy(&user[bi].Ql);
      VecDestroy(&user[bi].shrr);
      ///////////////////////////////
      VecDestroy(&user[bi].P_o);
      VecDestroy(&user[bi].Nvert_o);
      }

      VecDestroy(&user[bi].lCsi);
      VecDestroy(&user[bi].lEta);
      VecDestroy(&user[bi].lZet);
      VecDestroy(&user[bi].lICsi);
      VecDestroy(&user[bi].lIEta);
      VecDestroy(&user[bi].lIZet);
      VecDestroy(&user[bi].lJCsi);
      VecDestroy(&user[bi].lJEta);
      VecDestroy(&user[bi].lJZet);
      VecDestroy(&user[bi].lKCsi);
      VecDestroy(&user[bi].lKEta);
      VecDestroy(&user[bi].lKZet);
      VecDestroy(&user[bi].lGridSpace);
      VecDestroy(&user[bi].lCent);

      if(level==usermg->mglevels-1) {
      VecDestroy(&user[bi].lUcont);
      VecDestroy(&user[bi].lUcat);
      VecDestroy(&user[bi].lNFace);
      VecDestroy(&user[bi].lUcont_o);
      VecDestroy(&user[bi].lUcont_rm1);

      VecDestroy(&user[bi].lVcont);
      if (rotateframe)
      VecDestroy(&user[bi].lWcat);

      VecDestroy(&user[bi].lItfc); 
      }
      VecDestroy(&user[bi].lAj);
      VecDestroy(&user[bi].lIAj);
      VecDestroy(&user[bi].lJAj);
      VecDestroy(&user[bi].lKAj);
      VecDestroy(&user[bi].lNvert);
      VecDestroy(&user[bi].lPhi);

      if(level==usermg->mglevels-1) {
      VecDestroy(&user[bi].lP);
      /////////////////////////////
      VecDestroy(&user[bi].lItfcQ); 
      /////////////////////////////
      VecDestroy(&user[bi].lNvert_o);
      }     
     
      if(les && level==usermg->mglevels-1) {
	VecDestroy(&user[bi].lCs);	
	VecDestroy(&user[bi].lNu_t);		
      }
      if(rans  && level==usermg->mglevels-1) {
	VecDestroy(&user[bi].lK_Omega);	
	VecDestroy(&user[bi].lK_Omega_o);	
	VecDestroy(&user[bi].Distance);

	VecDestroy(&user[bi].K_Omega);	
	VecDestroy(&user[bi].K_Omega_o);	
	VecDestroy(&user[bi].lNu_t);		
	
	if(rans==3) {
	  VecDestroy(&user[bi].lF1);
	}
	//	DMDestroy(&user[bi].fda2);
      }

      //DMDestroy(&user[bi].fda);
      DMDestroy(&user[bi].da);
    }//for bi

    ////////////////////////////
    VecDestroy(&Umult_int);
    VecDestroy(&U_int);
    DMDestroy(&int_packer);
    MatDestroy(&Int_matrix);
    ////////////////////////////
    
    PetscInt implicit;
    PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
    /* if(level==usermg->mglevels-1 && implicit==5) { */
/*       DMCompositeDestroy(usermg->mgctx[level].packer); */
/*       DMMGDestroy(usermg->mgctx[level].dmmg); */
/*     } */
    PetscFree(user);
  }//for level
  PetscFree(usermg->mgctx);
  //  PetscFree(usermg);
  return 0;
}





PetscErrorCode Interpolation_matrix(UserMG *usermg) // for two block
{
  
  PetscErrorCode       ierr;
  PetscInt             its;
  PetscInt             i,j,k;
  PetscInt             blk;
  FILE                 *fd1;
  char                 filen[80];
  PetscInt             rank,size,bi,level;
  //DM                   packer;
  DMDALocalInfo        info1,info2;
  Vec                  Ub[block_number];
  PetscViewer          viewer;
  

  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  level = usermg->mglevels-1;
  UserCtx *user = usermg->mgctx[level].user;

  
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&int_packer);CHKERRQ(ierr);
  
  for(bi=0;bi<block_number;bi++){
    ierr = DMCompositeAddDM(int_packer,user[bi].fda);CHKERRQ(ierr);
  }
  
  DMSetUp(int_packer);
  
  ierr = DMCreateGlobalVector(int_packer,&U_int);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(int_packer,&Umult_int);CHKERRQ(ierr);
  
  ////////////////////
  DMDAGetLocalInfo(user[0].da,&info1);
  DMDAGetLocalInfo(user[1].da,&info2);
  
  PetscInt	xs1 = info1.xs, xe1 = info1.xs + info1.xm;
  PetscInt  	ys1 = info1.ys, ye1 = info1.ys + info1.ym;
  PetscInt	zs1 = info1.zs, ze1 = info1.zs + info1.zm;
  PetscInt	mx1 = info1.mx, my1 = info1.my, mz1 = info1.mz;
  
  PetscInt	xs2 = info2.xs, xe2 = info2.xs + info2.xm;
  PetscInt  	ys2 = info2.ys, ye2 = info2.ys + info2.ym;
  PetscInt	zs2 = info2.zs, ze2 = info2.zs + info2.zm;
  PetscInt	mx2 = info2.mx, my2 = info2.my, mz2 = info2.mz;

  // PetscPrintf(PETSC_COMM_SELF, "rank:%d---xs1=%d ys1=%d zs1=%d xe1=%d ye1=%d ze1=%d xm1=%d ym1=%d zm1=%d !\n",rank,xs1,ys1,zs1,xe1,ye1,ze1,info1.xm,info1.ym,info1.zm);
  //PetscPrintf(PETSC_COMM_SELF, "rank:%d---xs2=%d ys2=%d zs2=%d xe2=%d ye2=%d ze2=%d xm2=%d ym2=%d zm2=%d \n",rank,xs2,ys2,zs2,xe2,ye2,ze2,info2.xm,info2.ym,info2.zm);
  
  /////////////////////////////////////////////////////////////////// allocate local vectors
  //PetscPrintf(PETSC_COMM_SELF, "%d %d %d !\n",info1.xm,info1.ym,info1.zm);
  //PetscPrintf(PETSC_COMM_SELF, "%d %d %d !\n",info2.xm,info2.ym,info2.zm);


  DMCompositeGetAccess(int_packer,U_int,&Ub[0],&Ub[1]);
  
  for (bi=0; bi<block_number; bi++) {
    VecCopy(user[bi].Ucat, Ub[bi]);
  }
  
  DMCompositeRestoreAccess(int_packer,U_int,&Ub[0],&Ub[1]);


  PetscInt N11=0,N2=0,N1_local=0,N2_local=0,NN=0,NN_local=0;
  //Mat Int_matrix;
  
  VecGetSize(U_int,&N11);
  NN+=N11;
  VecGetLocalSize(U_int,&N1_local);
  NN_local+=N1_local;

  //PetscPrintf(PETSC_COMM_SELF, "Size %d !\n",NN);
  //PetscPrintf(PETSC_COMM_SELF, "rank:%d ,Size local %d !\n",rank, N1_local);

/////////////////////////////////////////////////////////////////// create and allocate a matrix
  ISLocalToGlobalMapping *ltogs,ltog;
  DMCompositeGetISLocalToGlobalMappings(int_packer,&ltogs);
  ISLocalToGlobalMappingConcatenate(PETSC_COMM_SELF,block_number,ltogs,&ltog);
  
  MatCreateAIJ(PETSC_COMM_WORLD, NN_local, NN_local, NN, NN,8, PETSC_NULL,8, PETSC_NULL,&Int_matrix);
  MatSetOption(Int_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  MatSetLocalToGlobalMapping(Int_matrix,ltog,ltog);
  for (bi=0; bi<block_number; bi++) {
    ISLocalToGlobalMappingDestroy(&ltogs[bi]);
  }
  ISLocalToGlobalMappingDestroy(&ltog);

  
  /////////////////////////////////////////////////////////////// obtaining global indeces(starting index for each cpu)

  PetscInt *dg_start,*dgl_start,*dl_start_blk0,*dl_start_blk1;
  
  PetscInt dl_start_blk00,dl_start_blk11;
  PetscMalloc(size*sizeof(PetscInt), &(dgl_start));
  PetscMalloc(size*sizeof(PetscInt), &(dg_start));
  PetscMalloc(size*sizeof(PetscInt), &(dl_start_blk0));
  PetscMalloc(size*sizeof(PetscInt), &(dl_start_blk1));
  
  PetscInt cpu=0;
  
  dl_start_blk00=3*(info1.xm*info1.ym*info1.zm);
  dl_start_blk11=3*(info2.xm*info2.ym*info2.zm);

  
  for(cpu=0;cpu<size;cpu++){
    dg_start[cpu]=0;
    dgl_start[cpu]=0;
    dl_start_blk0[cpu]=0;
    dl_start_blk1[cpu]=0;
  }
  

  MPI_Gather(&dl_start_blk00, 1, MPIU_INT,dl_start_blk0, 1, MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&dl_start_blk11, 1, MPIU_INT,dl_start_blk1, 1, MPIU_INT, 0, MPI_COMM_WORLD);
  PetscBarrier(NULL);
  
  MPI_Bcast(dl_start_blk0,size , MPIU_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dl_start_blk1,size , MPIU_INT, 0, MPI_COMM_WORLD);
  PetscBarrier(NULL);

  cpu=0;
  while(cpu<(rank)){
    dgl_start[rank]+=(dl_start_blk0[cpu]+dl_start_blk1[cpu]);
    cpu++;
   }

  for (i=0;i<size;i++){
    MPI_Allreduce(&dgl_start[i],&dg_start[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
  }

  PetscBarrier(NULL);
  //PetscPrintf(PETSC_COMM_SELF, "rank:%d, dg: %d %d %d %d!\n",rank,dg_start[0],dg_start[1],dg_start[2],dg_start[3]);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////// sending distribution info to all cpus
  PetscInt *cpu_lxs11,*cpu_lys11,*cpu_lzs11,*cpu_lxe11,*cpu_lye11,*cpu_lze11,*cpu_lxm11,*cpu_lym11,*cpu_lzm11;
  PetscInt *cpu_xs11,*cpu_ys11,*cpu_zs11,*cpu_xe11,*cpu_ye11,*cpu_ze11,*cpu_xm11,*cpu_ym11,*cpu_zm11;
  PetscInt *cpu_lxs22,*cpu_lys22,*cpu_lzs22,*cpu_lxe22,*cpu_lye22,*cpu_lze22,*cpu_lxm22,*cpu_lym22,*cpu_lzm22;
  PetscInt *cpu_xs22,*cpu_ys22,*cpu_zs22,*cpu_xe22,*cpu_ye22,*cpu_ze22,*cpu_xm22,*cpu_ym22,*cpu_zm22;
  
  PetscMalloc(size*sizeof(PetscInt), &(cpu_lxs11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lys11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzs11));
  PetscMalloc(size*sizeof(PetscInt), &(cpu_lxm11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lym11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzm11));
  PetscMalloc(size*sizeof(PetscInt), &(cpu_xs11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ys11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zs11));
  PetscMalloc(size*sizeof(PetscInt), &(cpu_xm11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ym11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zm11));
  /////////////////////////////blk2
  PetscMalloc(size*sizeof(PetscInt), &(cpu_lxs22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lys22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzs22));
  PetscMalloc(size*sizeof(PetscInt), &(cpu_lxm22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lym22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzm22));
  PetscMalloc(size*sizeof(PetscInt), &(cpu_xs22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ys22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zs22));
  PetscMalloc(size*sizeof(PetscInt), &(cpu_xm22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ym22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zm22));
  
  for(cpu=0;cpu<size;cpu++){
    cpu_lxs11[cpu]=0;cpu_lys11[cpu]=0;cpu_lzs11[cpu]=0;cpu_xs11[cpu]=0;cpu_ys11[cpu]=0;cpu_zs11[cpu]=0;
    cpu_lxm11[cpu]=0;cpu_lym11[cpu]=0;cpu_lzm11[cpu]=0;cpu_xm11[cpu]=0;cpu_ym11[cpu]=0;cpu_zm11[cpu]=0;
    cpu_lxs22[cpu]=0;cpu_lys22[cpu]=0;cpu_lzs22[cpu]=0;cpu_xs22[cpu]=0;cpu_ys22[cpu]=0;cpu_zs22[cpu]=0;
    cpu_lxm22[cpu]=0;cpu_lym22[cpu]=0;cpu_lzm22[cpu]=0;cpu_xm22[cpu]=0;cpu_ym22[cpu]=0;cpu_zm22[cpu]=0;
  }

  ///------------------------///
  cpu_lxs11[rank]=info1.xs;cpu_lys11[rank]=info1.ys;cpu_lzs11[rank]=info1.zs;
  cpu_lxm11[rank]=info1.xm;cpu_lym11[rank]=info1.ym;cpu_lzm11[rank]=info1.zm;
  
  cpu_lxs22[rank]=info2.xs;cpu_lys22[rank]=info2.ys;cpu_lzs22[rank]=info2.zs;
  cpu_lxm22[rank]=info2.xm;cpu_lym22[rank]=info2.ym;cpu_lzm22[rank]=info2.zm;

  PetscBarrier(NULL);
  
  for (i=0;i<size;i++){
    MPI_Allreduce(&cpu_lxs11[i],&cpu_xs11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lys11[i],&cpu_ys11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lzs11[i],&cpu_zs11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lxm11[i],&cpu_xm11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lym11[i],&cpu_ym11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lzm11[i],&cpu_zm11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

    MPI_Allreduce(&cpu_lxs22[i],&cpu_xs22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lys22[i],&cpu_ys22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lzs22[i],&cpu_zs22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lxm22[i],&cpu_xm22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lym22[i],&cpu_ym22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&cpu_lzm22[i],&cpu_zm22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
  }
  
  PetscBarrier(NULL);
  //PetscPrintf(PETSC_COMM_SELF, "rank:%d, lx: %d %d!\n",rank,cpu_xm11[3],cpu_xm22[3]);
  
  //////----------------------------------------------------///
  PetscInt cpu_xs[block_number][size];PetscInt cpu_xm[block_number][size];
  PetscInt cpu_ys[block_number][size];PetscInt cpu_ym[block_number][size];
  PetscInt cpu_zs[block_number][size];PetscInt cpu_zm[block_number][size];

  for (bi=0; bi<block_number; bi++) {
    
    for (i=0;i<size;i++){
      if(bi==0){
	cpu_xs[bi][i]=cpu_xs11[i];
	cpu_xm[bi][i]=cpu_xm11[i];
	cpu_ys[bi][i]=cpu_ys11[i];
	cpu_ym[bi][i]=cpu_ym11[i];
	cpu_zs[bi][i]=cpu_zs11[i];
	cpu_zm[bi][i]=cpu_zm11[i];
      }else if(bi==1){
	cpu_xs[bi][i]=cpu_xs22[i];
	cpu_xm[bi][i]=cpu_xm22[i];
	cpu_ys[bi][i]=cpu_ys22[i];
	cpu_ym[bi][i]=cpu_ym22[i];
	cpu_zs[bi][i]=cpu_zs22[i];
	cpu_zm[bi][i]=cpu_zm22[i];
      }
      //PetscPrintf(PETSC_COMM_SELF, "rank:%d  bi: %d, lx: %d %d %d %d %d %d !\n",i,bi,cpu_xs[bi][i],cpu_xm[bi][i],cpu_ys[bi][i],cpu_ym[bi][i],cpu_zs[bi][i],cpu_zm[bi][i]);
    }
  }
  PetscBarrier(NULL);

  /////////////////////////////////////----------------------------///////////////////////////////////////////////// reading the coeficients and host
  PetscLogDouble v1,v2,elapsed_time;
  PetscTime(&v1);
  
  PetscInt  hb=0;
  
  for (bi=0; bi<block_number; bi++) {
    
    DMDALocalInfo	info = user[bi].info;
    
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt    lmx, lmy, lmz;



    PetscInt d,dg,cll=0;
    PetscInt sum=0;
    PetscInt itfc_number=0;
    PetscReal lval[8];
    PetscInt	row=0,lcol[8];
    PetscInt col[8];
    PetscReal val[8];
    PetscReal x, y, z;
    PetscInt  itfnumber=0,itfnumber1=0;
    PetscInt d_freedom=0;
    
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes
       from the host nodes.
       itfc is the velocity at the cell corners
       hostU is the velocity at the cell centers of the host block */
    /**************************************************************************************************************************/
    
        
    for (itfc_number=0; itfc_number<user[bi].itfcptsnumber; itfc_number++) {
      
      x = user[bi].itfchostx[itfc_number];
      y = user[bi].itfchosty[itfc_number];
      z = user[bi].itfchostz[itfc_number];
      
      
      dg=0;d=0;cll=0;
      
      for (i=0;i<8;i++){
	lcol[i]=0;
	col[i]=0;
	lval[i]=0.0;
	val[i]=0.0;
      }
      
      cpu=0;
      
      hb = user[bi].itfchostB[itfc_number];

      if(user[bi].itfchostx[itfc_number]>-1.25){

	if((user[bi].itfcK[itfc_number])>=zs&&(user[bi].itfcK[itfc_number])<ze)
	  if((user[bi].itfcJ[itfc_number])>=ys&&(user[bi].itfcJ[itfc_number])<ye)
	    if((user[bi].itfcI[itfc_number])>=xs&&(user[bi].itfcI[itfc_number])<xe)
	      {

		for(cpu=0;cpu<size;cpu++){
		  
		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number];
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[0]=dg;
			  val[0]=(1-x) * (1-y) * (1-z);
			  itfnumber++;
			}
		  
		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number]+1;
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[1]=dg;
			  val[1]=x * (1-y) * (1-z);
			  itfnumber++;
			}
		  
		  
		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number];
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[2]=dg;
			  val[2]=(1-x) * y * (1-z);
			  itfnumber++;
			}
		  
		  
		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number];
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[3]=dg;
			  val[3]=(1-x) * (1-y) * z;
			  itfnumber++;
			}
		  
		  
		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number]+1;
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[4]=dg;
			  val[4]= x * y * (1-z);
			  itfnumber++;
			}
		  
		  
		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number]+1;
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[5]=dg;
			  val[5]=x * (1-y) * z;
			  itfnumber++;
			}
		  
		  
		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number];
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[6]=dg;
			  val[6]=(1-x) * y * z;
			  itfnumber++;
			}
		  
		  
		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu]))
		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu]))
		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu]))
			{
			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number]+1;
			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
			  dg=dg_start[cpu]+hb*dl_start_blk0[cpu]+3*d;
			  col[7]=dg;
			  val[7]=x * y * z;
			  itfnumber++;
			}
		  
		}
		
		
		//-------------------------------------------------------------------------//
		k=user[bi].itfcK[itfc_number];j=user[bi].itfcJ[itfc_number];i=user[bi].itfcI[itfc_number];
		d=lidxLocal1_matrix(i, j, k, &user[bi],bi);
		dg=dg_start[rank]+bi*dl_start_blk0[rank]+3*d;
		
		for(d_freedom=0;d_freedom<3;d_freedom++){
		  row=dg+d_freedom;//lidxLocal(8, 8, 8, &user,0);
		  itfnumber1++;
		  
		  for(i=0;i<8;i++){
		    col[i]= col[i]+d_freedom;
		  }
		  
		  MatSetValues(Int_matrix,1,&row,8,col,val,INSERT_VALUES);
		}
		
	      }
      }
      
    }
    
    //PetscPrintf(PETSC_COMM_SELF, "bi:%d number: %d %d!\n",bi,itfnumber,itfnumber1);
  }


  MatAssemblyBegin(Int_matrix,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Int_matrix,MAT_FINAL_ASSEMBLY);


  PetscTime(&v2);

  elapsed_time = v2 - v1;

  PetscPrintf(PETSC_COMM_WORLD,"elapsed_time:%le \n", elapsed_time);

/*   sprintf(filen, "Matfield.dat"); */
/*   PetscViewerASCIIOpen(PETSC_COMM_WORLD, filen, &viewer); */
/*   MatView(Int_matrix, viewer); */



/*   /////////////////////////////////////////////////////////////////// */
/*  ierr = DMCompositeGetAccess(packer,Umult,&Ub[0],&Ub[1]);CHKERRQ(ierr); */

/*  //sprintf(filen, "umult%5.5d_%1.1d.dat", 0, 1); */
/*   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"umult1" , &viewer); */
/*   VecView(user[1].Ucat, viewer); */
/*   PetscViewerDestroy(&viewer); */

/*   //sprintf(filen, "umult%5.5d_%1.1d.dat", 0, 0); */
/*   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"umult0", &viewer); */
/*   VecView(user[0].Ucat, viewer); */
/*   PetscViewerDestroy(&viewer); */

/*   ierr = DMCompositeRestoreAccess(packer,Umult,&Ub[0],&Ub[1]);CHKERRQ(ierr); */
  
/*   ////////////////////////////////////////////////////////////////////////////////////// */

  for (bi=0; bi<block_number; bi++) {
    ierr = VecDestroy(&Ub[bi]);CHKERRQ(ierr);
  }
  
  //PetscViewerDestroy(&viewer)

  //////////////////////////////////////////////////////////////////////
  PetscFree(dgl_start);PetscFree(dg_start);PetscFree(dl_start_blk0); PetscFree(dl_start_blk1);
  PetscFree(cpu_lxs11);PetscFree(cpu_lys11);PetscFree(cpu_lzs11);PetscFree(cpu_lxm11);PetscFree(cpu_lym11);PetscFree(cpu_lzm11);
  PetscFree(cpu_lxs11);PetscFree(cpu_ys11);PetscFree(cpu_zs11);PetscFree(cpu_xm11);PetscFree(cpu_ym11);PetscFree(cpu_zm11);
  PetscFree(cpu_lxs22);PetscFree(cpu_lys22);PetscFree(cpu_lzs22);PetscFree(cpu_lxm22);PetscFree(cpu_lym22);PetscFree(cpu_lzm22);
  PetscFree(cpu_lxs22);PetscFree(cpu_ys22);PetscFree(cpu_zs22);PetscFree(cpu_xm22);PetscFree(cpu_ym22);PetscFree(cpu_zm22);
  //////////////////////////////////////////////////////////////////////


  PetscPrintf(PETSC_COMM_WORLD, "Done!\n");

  return 0;
}




PetscInt lidxLocal_matrix(PetscInt i, PetscInt j, PetscInt k,PetscInt blk,PetscInt cpu,PetscInt xs[block_number][cpu_size],PetscInt xm[block_number][cpu_size],PetscInt ys[block_number][cpu_size],PetscInt ym[block_number][cpu_size],PetscInt zs[block_number][cpu_size],PetscInt zm[block_number][cpu_size])
{
  
  return ((k-zs[blk][cpu]) * (xm[blk][cpu]*ym[blk][cpu]) + (j-ys[blk][cpu])*(xm[blk][cpu]) + (i-xs[blk][cpu]));
}



PetscInt lidxLocal1_matrix(PetscInt i, PetscInt j, PetscInt k, UserCtx *user,PetscInt blk)
{
  DMDALocalInfo	info;

  DMDAGetLocalInfo(user->da,&info); 
  
  PetscInt	xs, xe, ys, ye, zs, ze;
  
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  return ((k-zs) * (info.xm*info.ym) + (j-ys)*(info.xm) + (i-xs));
}
