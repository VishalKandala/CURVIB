#include "variables.h"
#include "petscvec.h" 
#include <petsctime.h>
#include "petscdmcomposite.h"
#include "petscdmda.h"
#include "petsctime.h"

extern PetscInt ti, moveframe, blank;
extern PetscReal FluxInSum,FluxOutSum,FluxOutSumcont;
extern PetscReal Flux_in, angle,CMy_c, CMx_c,CMz_c;
extern PetscInt block_number;
extern PetscInt inletprofile,visflg;
extern PetscInt inletface,outletface,catcorr;
PetscReal FluxInSumB[5],FluxOutSumB[5],AreaOutB[5],OFC[5];  // Maximum block number assumed 5,OFC:Outlet Flux Coefficient
PetscInt ts_p_cycle;
PetscReal *Flux_waveform;
PetscErrorCode Contra2Cart(UserCtx *user);
PetscErrorCode VTKOut(UserCtx *user);
void wall_function (UserCtx *user, double sc, double sb, 
		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, 
		    PetscReal *ustar, double nx, double ny, double nz);
void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], 
		      double nj[3], double nk[3]);
PetscErrorCode GhostNodeVelocity(UserCtx *user);

//////////////////////////////
Mat Int_matrix;
DM  int_packer;
Vec U_int,Umult_int,Flux_Waveform;
extern int cpu_size;

PetscInt lidxLocal1_matrix(PetscInt i, PetscInt j, PetscInt k, UserCtx *user,PetscInt blk);
PetscInt lidxLocal_matrix(PetscInt i, PetscInt j, PetscInt k,PetscInt blk,PetscInt cpu,PetscInt xs[block_number][cpu_size],PetscInt xm[block_number][cpu_size],PetscInt ys[block_number][cpu_size],PetscInt ym[block_number][cpu_size],PetscInt zs[block_number][cpu_size],PetscInt zm[block_number][cpu_size]);
/////////////////////////////



#define INLET 5
#define OUTLET 4
#define SOLIDWALL 1
#define SYMMETRIC 3
#define FARFIELD 6

#define FLUX_THRESHOLD 1e-16

PetscErrorCode Flux_Waveform_Read(UserCtx *user)
{
  PetscInt i,rank;
  PetscOptionsGetInt(PETSC_NULL,"-ts_p_cycle",&ts_p_cycle,PETSC_NULL); 
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscMalloc(ts_p_cycle*sizeof(PetscReal),&Flux_waveform);
  PetscPrintf(PETSC_COMM_WORLD,"Reading flux waveform from inlet.dat \n");
  if (!rank) {
  FILE *fd;
  fd = fopen("inlet.dat", "r");
  for(i=1;i<=ts_p_cycle;i++){
  fscanf(fd, "%le", &(Flux_waveform[i])); 
   }
   fclose(fd);
   MPI_Bcast(Flux_waveform, ts_p_cycle, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else { 
   MPI_Bcast(Flux_waveform, ts_p_cycle, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  PetscPrintf(PETSC_COMM_WORLD, "Waveform Read for %d timesteps per cycle \n",ts_p_cycle);
 return(0);
}

PetscReal Pulsatile_Plug_Inlet_Flux(UserCtx *user, PetscReal area)
{
  PetscReal uin;
  
  PetscInt tstep,cycle,flux_vel_switch;
    
  cycle = ((PetscInt)(ti / ts_p_cycle)); 
  tstep = ti - cycle*ts_p_cycle;
  cycle=cycle+1;
  FILE *fd;
  fd = fopen("inlet.dat","r");
  if(fd) fscanf(fd,"%i",&flux_vel_switch);  // Read the first line of inlet.dat, to get whether inflow is flux or velocity.
  fclose(fd);
  PetscPrintf(PETSC_COMM_WORLD,"inlet.dat read, vel(0)/flux(1): %d \n",flux_vel_switch);
  PetscPrintf(PETSC_COMM_WORLD, "Cycle: %d step: %d \n",cycle,tstep);
  if(flux_vel_switch) uin = Flux_waveform[tstep]/area;
  else uin = Flux_waveform[tstep];
  return (uin);
} 

Cmpnts InflowCenter(UserCtx *user)
{
  PetscInt i,j,k,rank; 
  Cmpnts center,***cent,***coor;
  Vec Coor;
  DM       da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert; //local working array
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  VecDuplicate(user->lNvert,&user->RFC);
//  PetscPrintf(PETSC_COMM_SELF,"xs,xe,mx,rank: %d,%d,%d,%d \n",xs,xe,mx,rank);
  // PetscMalloc(block_number*sizeof(PetscReal), &FluxInSumB);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  //DMDAGetGhostedCoordinates(da, &Coor);
  DMGetCoordinatesLocal(da, &Coor);
//  DMDAVecGetArray(da, user->RFC, &RR);
  DMDAVecGetArray(fda, user->Cent, &cent);
  
}
PetscErrorCode InflowFlux(UserCtx *user) 
{
  PetscInt     i, j, k,rank;
  PetscReal    FluxIn,r, uin0, uin, xc, yc,zc,***RR ;
  PetscReal    lAreaIn,AreaSumIn,Re;
  Vec          Coor;
  Cmpnts       ***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet, ***cent;  
  
  
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert; //local working array
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  VecDuplicate(user->lNvert,&user->RFC);
//  PetscPrintf(PETSC_COMM_SELF,"xs,xe,mx,rank: %d,%d,%d,%d \n",xs,xe,mx,rank);
  // PetscMalloc(block_number*sizeof(PetscReal), &FluxInSumB);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  //DMDAGetGhostedCoordinates(da, &Coor);
  DMGetCoordinatesLocal(da, &Coor);

  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->RFC, &RR);
  DMDAVecGetArray(fda, user->Cent, &cent);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

//  if (user->bctype[4] == FARFIELD) FluxIn=0.;
//  if (user->bctype[4] == SOLIDWALL) FluxIn=0.;
//  if (user->bctype[4] == SYMMETRIC) FluxIn=0.;

  
    // S-Calc uin
//  PetscPrintf(PETSC_COMM_WORLD,"Inletprofile : %d \n",inletprofile);  
  if (inletprofile == 1) {
    PetscOptionsGetReal(PETSC_NULL, "-uin", &uin, PETSC_NULL);
//      PetscPrintf(PETSC_COMM_WORLD,"Inlet Velocity : %f \n",uin);  
//  }else if(inletprofile=2) // Fully Developed flow
        
 //------------------------------------------------------------------------------
  }else if(inletprofile==3){ // LVAD Pulsatile plug flow inlet
    Flux_Waveform_Read(user);
    uin=0.;
  }
  else {
    PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
    uin = 0.;
  }
  // E-Calc uin
  //fn is face number which is 0 to 5
  PetscInt fn;PetscReal d;
  FluxIn=0.0;
  lAreaIn=0.0;
  for (fn=0; fn<6; fn++) {
   if (user->bctype[fn] == INLET) {
      inletface = fn;
   if(visflg) PetscPrintf(PETSC_COMM_WORLD,"Inlet detected at face: %d \n",inletface);
    switch(fn){
      // face 0
    case 0:
       if (xs==0) {
//	PetscPrintf(PETSC_COMM_WORLD,"Case 0 entered \n");
       	i = xs;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
              if(nvert[k][j][i+1]<0.1){
                 d=sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
                 lAreaIn+=d;
              }
          }
        }
        if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,lAreaIn);
        if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"Inlet Velocity : %f \n",uin);   
        for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    //S-Calc covarient velocity componenet if it is inside
	    if (nvert[k][j][i+1]<0.1) {
	      d=sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
//              lAreaIn+=d;
//              if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,d); 
              ucont[k][j][i].x = uin*d;
	      ubcs[k][j][i].x = uin*csi[k][j][i].x/d;
	      ubcs[k][j][i].y = uin*csi[k][j][i].y/d;
	      ubcs[k][j][i].z = uin*csi[k][j][i].z/d;
	      ucat[k][j][i+1].x = ubcs[k][j][i].x;
	      ucat[k][j][i+1].y = ubcs[k][j][i].y;
	      ucat[k][j][i+1].z = ubcs[k][j][i].z;
	      FluxIn += ucont[k][j][i].x;
            
	    } // nvert  
	    //E-Calc covarient velocity componenet if it is inside
	  }// location for loop
	}// location for loop

      } // if statement 
      break;
      // face 1
    case 1:
      if (xe==mx) {
	i = mx-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
              if(nvert[k][j][i]<0.1){
                 d=sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
                 lAreaIn+=d;
              }
          }
        }
        if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,lAreaIn);
        if(visflg){
          PetscPrintf(PETSC_COMM_SELF,"Inlet Velocity : %f \n",uin);
        }
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    //S-Calc covarient velocity componenet if it is inside
	    if (nvert[k][j][i]<0.1) {
              d=sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
//	      lAreaIn+=d;
              ucont[k][j][i].x = -uin*d;
	      ubcs[k][j][i+1].x = -uin*csi[k][j][i].x/d;
	      ubcs[k][j][i+1].y = -uin*csi[k][j][i].y/d;
	      ubcs[k][j][i+1].z = -uin*csi[k][j][i].z/d;
//------------------------------------------------------
//             if(k<=(PetscInt)(ceil((lze-lzs)/2)+5) &&  k>=(PetscInt)(ceil((lze-lzs)/2)-5) && j >= (PetscInt)(ceil((lye-lys)/2)-5) && j <= (PetscInt)(ceil((lye-lys)/2)+5)){
//              PetscPrintf(PETSC_COMM_SELF,"ubcs.x=%le,csi.x =%le \n",ubcs[k][j][i+1].x,csi[k][j][i].x);
//                }
//------------------------------------------------------
	      ucat[k][j][i].x = ubcs[k][j][i+1].x;
	      ucat[k][j][i].y = ubcs[k][j][i+1].y;
	      ucat[k][j][i].z = ubcs[k][j][i+1].z;

// ---------- diagnostics ---------- Vishal Kandala 
//             if(k<=(PetscInt)(ceil((lze-lzs)/2)+5) &&  k>=(PetscInt)(ceil((lze-lzs)/2)-5) && j >= (PetscInt)(ceil((lye-lys)/2)-5) && j <= (PetscInt)(ceil((lye-lys)/2)+5)){
//              PetscPrintf(PETSC_COMM_SELF,"ucat.x=%le \n",ucat[k][j][i].x);
//                }
// -----------------------------------------------
	      FluxIn += ucont[k][j][i].x;
	    }
	    //E-Calc covarient velocity componenet if it is inside
	  }// location for loop
	}// location for loop
      }
      break;
      // face 2
    case 2:
      if (ys==0) {
	j = 0;
        for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
              if(nvert[k][j+1][i]<0.1){
                 d=sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);
                 lAreaIn+=d;
              }
          }
        }
        if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,lAreaIn);
        if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"Inlet Velocity : %f \n",uin);
        
        for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    //S-Calc covarient velocity componenet if it is inside
	    if (nvert[k][j+1][i]<0.1) {
	      d=sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);
//              lAreaIn+=d;
              ucont[k][j][i].y = uin*d;
	      ubcs[k][j][i].x = uin*eta[k][j][i].x/d;
	      ubcs[k][j][i].y = uin*eta[k][j][i].y/d;
	      ubcs[k][j][i].z = uin*eta[k][j][i].z/d;
	      ucat[k][j+1][i].x = ubcs[k][j][i].x;
	      ucat[k][j+1][i].y = ubcs[k][j][i].y;
	      ucat[k][j+1][i].z = ubcs[k][j][i].z;
	      FluxIn += ucont[k][j][i].y;
	    }
	    //E-Calc covarient velocity componenet if it is inside
	  }// location for loop
	}// location for loop

      }
      break;
      // face 3
    case 3:
      if (ye==my) {

	j = my-2;
        for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
              if(nvert[k][j][i]<0.1){
                 d=sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);
                 lAreaIn+=d;
              }
          }
        }
        if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,lAreaIn);
        if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"Inlet Velocity : %f \n",uin);
        
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    //S-Calc covarient velocity componenet if it is inside
	    if (nvert[k][j][i]<0.1) {
 	      d=sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);   
//	      lAreaIn+=d;
              ucont[k][j][i].y = -1.0*uin*d;
	      ubcs[k][j+1][i].x = -uin*eta[k][j][i].x/d;
	      ubcs[k][j+1][i].y = -uin*eta[k][j][i].y/d;
	      ubcs[k][j+1][i].z = -uin*eta[k][j][i].z/d;
	      ucat[k][j][i].x = ubcs[k][j+1][i].x;
	      ucat[k][j][i].y = ubcs[k][j+1][i].y;
	      ucat[k][j][i].z = ubcs[k][j+1][i].z;
	      FluxIn += ucont[k][j][i].y;
	    }
	    //E-Calc covarient velocity componenet if it is inside
	  }// location for loop
	}// location for loop

      }
      break;
      // face 4
    case 4:
      if (zs==0) {
	
	k = 0;
        for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
              if(nvert[k+1][j][i]<0.1){
                 d=sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
                 lAreaIn+=d;
              }
          }
        }
        if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,lAreaIn);
        if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"Inlet Velocity : %f \n",uin);
        
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    //S-Calc covarient velocity componenet if it is inside
	    if (nvert[k+1][j][i]<0.1) {
	      d=sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
              ucont[k][j][i].z = uin*d;
	      ubcs[k][j][i].x = uin*zet[k][j][i].x/d;
	      ubcs[k][j][i].y = uin*zet[k][j][i].y/d;
	      ubcs[k][j][i].z = uin*zet[k][j][i].z/d;
	      ucat[k+1][j][i].x = ubcs[k][j][i].x;
	      ucat[k+1][j][i].y = ubcs[k][j][i].y;
	      ucat[k+1][j][i].z = ubcs[k][j][i].z;
	      FluxIn += ucont[k][j][i].z;
	    }
	    //E-Calc covarient velocity componenet if it is inside
	  }// location for loop
	}// location for loop
  	
      }
      break;
      // face 5
    case 5:
      if (ze==mz) {	
	k = mz-2;
        for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
              if(nvert[k][j][i]<0.1){
                 d=sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
                 lAreaIn+=d;
              }
          }
        }
        if(inletprofile==3) uin = Pulsatile_Plug_Inlet_Flux(user,lAreaIn);
        if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"Inlet Velocity : %f \n",uin);
        
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    //S-Calc covarient velocity componenet if it is inside
	    if (nvert[k][j][i]<0.1) {
	      d=sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
              ucont[k][j][i].z = -uin*d;
	      ubcs[k+1][j][i].x = -uin*zet[k][j][i].x/d;
	      ubcs[k+1][j][i].y = -uin*zet[k][j][i].y/d;
	      ubcs[k+1][j][i].z = -uin*zet[k][j][i].z/d;
	      ucat[k][j][i].x = ubcs[k+1][j][i].x;
	      ucat[k][j][i].y = ubcs[k+1][j][i].y;
	      ucat[k][j][i].z = ubcs[k+1][j][i].z;
	      FluxIn += ucont[k][j][i].z;
	    }
	    //E-Calc covarient velocity componenet if it is inside
	  }// location for loop
	}// location for loop
   
      }
     //  else{
     //        FluxIn=0;
     //        lAreaIn=0;
//	    } 
      break;
    }//end switch

   }// end inlet check
   else if(user->bctype[fn]==SOLIDWALL){
     if(visflg) PetscPrintf(PETSC_COMM_WORLD,"Solid Wall detected at face: %d \n",fn);
  }
  else if(user->bctype[fn]==OUTLET){
      outletface=fn;
   if(visflg) PetscPrintf(PETSC_COMM_WORLD,"Outlet  detected at face: %d \n",outletface);
 }
          
//    PetscPrintf(PETSC_COMM_WORLD,"face: %di \n",fn);
  }// end face counter 

//Parallel
 MPI_Allreduce(&FluxIn,&FluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 MPI_Allreduce(&lAreaIn,&AreaSumIn,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 PetscBarrier(PETSC_NULL);
 user->FluxInSum = FluxInSum;
 if(visflg) PetscPrintf(PETSC_COMM_WORLD,"Inflow Flux - Area:  %le - %le \n",FluxInSum,AreaSumIn);    
 FluxInSumB[user->_this]=FluxInSum;

  
  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);

  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->RFC, &RR);

  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  VecDestroy(&user->RFC);

  return 0; 
}


PetscErrorCode OutflowFlux(UserCtx *user) {


  PetscInt      i, j, k;
  PetscReal     FluxOut,FluxOutcont;
  Cmpnts	***ucont, ***ucat, ***csi, ***eta, ***zet;
  PetscReal     ***nvert;
  DM            da = user->da, fda = user->fda;
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
  

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);
  
  //fn is face number which is 0 to 5
  PetscInt fn;
  FluxOut = 0.0, FluxOutcont = 0.0;
  PetscReal lArea=0.0,lAreaSum=0.0;
//  for (fn=0; fn<6; fn++) {
//    if (user->bctype[fn] == 4) {
//      outletface=fn;
//      PetscPrintf(PETSC_COMM_WORLD,"Outlet detected at face: %d \n",outletface); 
      switch(outletface){
	// face 0
      case 0:
	if (xs==0) {	  
	  i = xs;
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
              if(nvert[k][j][i+1]<0.1){
	        FluxOutcont +=(ucat[k][j][i+1].x*(csi[k][j][i].x) +
                          ucat[k][j][i+1].y*(csi[k][j][i].y) + 
                          ucat[k][j][i+1].z*(csi[k][j][i].z));         
                
                FluxOut += ucont[k][j][i].x;
	        lArea += sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			     (csi[k][j][i].y) * (csi[k][j][i].y) +
			     (csi[k][j][i].z) * (csi[k][j][i].z));
	       } // nvert if condition.
             }// location for loop
	  }// location for loop
	  
	}
	break;
	// face 1
      case 1:
	if (xe==mx) {
	  
	  i = mx-2;
	  for (k=lzs; k<lze; k++) {
	    for (j=lys; j<lye; j++) {
	      if(nvert[k][j][i]<0.1){
                FluxOutcont +=(ucat[k][j][i].x*(csi[k][j][i].x) +
                           ucat[k][j][i].y*(csi[k][j][i].y) +
                           ucat[k][j][i].z*(csi[k][j][i].z)); 
 
                FluxOut += ucont[k][j][i].x;
	        lArea += sqrt((csi[k][j][i].x) * (csi[k][j][i].x) +
			     (csi[k][j][i].y) * (csi[k][j][i].y) +
			     (csi[k][j][i].z) * (csi[k][j][i].z));
	      } // nvert
            }// location for loop
	  }// location for loop
	  
	}
	break;
	// face 2
      case 2:
	if (ys==0) { 
	  j = ys;
	  for (k=lzs; k<lze; k++) {
	    for (i=lxs; i<lxe; i++) {
              if(nvert[k][j+1][i]<0.1){
	      FluxOutcont +=(ucat[k][j+1][i].x*(eta[k][j][i].x) + 
                         ucat[k][j+1][i].y*(eta[k][j][i].y) + 
                         ucat[k][j+1][i].z*(eta[k][j][i].z)); 
              
              FluxOut += ucont[k][j][i].y;
	      lArea += sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			     (eta[k][j][i].y) * (eta[k][j][i].y) +
			     (eta[k][j][i].z) * (eta[k][j][i].z));
	      } // nvert
            }// location for loop
	  }// location for loop
	  
	}
	break;
	// face 3
      case 3:
	if (ye==my) { 
	  j = my-2;
	  for (k=lzs; k<lze; k++) {
	    for (i=lxs; i<lxe; i++) {
	      if(nvert[k][j][i]<0.1){
                FluxOutcont +=(ucat[k][j][i].x*(eta[k][j][i].x) +
                           ucat[k][j][i].y*(eta[k][j][i].y) +
                           ucat[k][j][i].z*(eta[k][j][i].z)); 
 
                FluxOut += ucont[k][j][i].y;
               // PetscPrintf(PETSC_COMM_SELF,"ucont.x,y,z=%le,%le,%le \n",ucont[k][j][i].x,ucont[k][j][i].y,ucont[k][j][i].z);
                lArea += sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			     (eta[k][j][i].y) * (eta[k][j][i].y) +
			     (eta[k][j][i].z) * (eta[k][j][i].z));
	      } // nvert
            }// location for loop
	  }// location for loop
	  
	}
	break;
	// face 4
      case 4:
	if (zs==0) { 
	  k = zs;
	  for (j=lys; j<lye; j++) {
	    for (i=lxs; i<lxe; i++) {
	     if(nvert[k+1][j][i]<0.1){
                FluxOutcont +=(ucat[k+1][j][i].x*(zet[k][j][i].x) + 
                          ucat[k+1][j][i].y*(zet[k][j][i].y) + 
                          ucat[k+1][j][i].z*(zet[k][j][i].z)); 
                
                FluxOut += ucont[k][j][i].z;
                lArea += sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			     (zet[k][j][i].y) * (zet[k][j][i].y) +
			     (zet[k][j][i].z) * (zet[k][j][i].z));
	      } //nvert
            }// location for loop
	  }// location for loop
	  
	}
	break;
	// face 5
      case 5:
	if (ze==mz) { 
	  k = mz-2;
	  for (j=lys; j<lye; j++) {
	    for (i=lxs; i<lxe; i++) {
	      if(nvert[k][j][i]<0.1){
                FluxOutcont +=(ucat[k][j][i].x*(zet[k][j][i].x) +
                           ucat[k][j][i].y*(zet[k][j][i].y) +
                           ucat[k][j][i].z*(zet[k][j][i].z)); 
                 FluxOut += ucont[k][j][i].z;
                 lArea += sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			     (zet[k][j][i].y) * (zet[k][j][i].y) +
			     (zet[k][j][i].z) * (zet[k][j][i].z));
	      } // nvert
            }// location for loop
	  }// location for loop
	}
	break;
      }//end switch
      
 //   }// end outlet check
 // }// end face counter 
  
  MPI_Allreduce(&lArea,&lAreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FluxOutcont,&FluxOutSumcont,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"Outflow Cont Flux - Cart Flux -  Area:  %le - %le - %le \n",FluxOutSum,FluxOutSumcont,lAreaSum);    
  user->FluxOutSum = FluxOutSum;  // Contravariant sum
  FluxOutSumB[user->_this]=FluxOutSum;
  AreaOutB[user->_this]=lAreaSum;

  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);
 
  return 0;
}


/* 
 INTERFACE TYPE inttype
cont   Normal       0
   Mixed Inlet  1
   Mixed Outlet 2
   Side no flow 3

inttype_ratio is set by inttype[6]
   0      constant velocity correction
   1      proportional (needs FluxOutSum>0)
   2      proportional to flux
   3      proportional to normal velocity (flux/area)
*/

PetscErrorCode Calc_Correction(PetscReal lFluxOutSum[6], PetscReal lAreaSum[6], 
			       PetscReal lFluxOutSum_abs[6],
			       PetscReal *ratio,PetscReal *FluxIn, PetscReal *FluxOutSum,
			       PetscReal *AreaSum, PetscInt inttype, 
			       PetscInt inttype_ratio, PetscInt surf, UserCtx *user)
{
  PetscInt      bi,i;

  PetscReal     FluxOutSum_abs=0.;
  PetscScalar   epsilon=1.e-10, sign_surf;

   if (surf%2) 
    sign_surf=-1; 
   else 
    sign_surf=1.;
  bi=user->_this;

  *FluxIn = FluxInSum;
  
  if (inttype==0) {    
    *FluxIn = sign_surf*FluxInSum;
    *FluxOutSum = lFluxOutSum[surf];
    *AreaSum = lAreaSum[surf];
    FluxOutSum_abs=lFluxOutSum_abs[surf];
  } else  if (inttype==1) {
    *FluxIn = -FluxInSum;
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  } else if (inttype==2) {
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  } else if (inttype==3) {
    *FluxIn = 0.;
    
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  } else if (inttype==4) {
    *FluxIn = user[bi].FluxInSum;
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  }
    
  if (*AreaSum<epsilon) *AreaSum=1.;

  if (inttype_ratio>1) {
    if (fabs(*FluxOutSum)< 1.e-8)
      *ratio=0.;
    else
      *ratio = (*FluxIn - *FluxOutSum) / FluxOutSum_abs;
  } else if (inttype_ratio)
    *ratio = *FluxIn/ *FluxOutSum;
  else
    *ratio = (*FluxIn - *FluxOutSum) / *AreaSum;

  
 
  

  return(0);
}


/* INTERFACE TYPE inttype
   Normal       0
   Mixed Inlet  1
   Mixed Outlet 2
   Side no flow 3

   flg : assigns the correction type
   flg == 0 normal corection proportional to area (const vel)
   flg == 2 correction porportional to flux
   flg == 3 correction porportional to normal velocity
*/

PetscErrorCode Block_Blank_Correction_adv(UserCtx *user, //Vec lUcor, 					  
					  PetscInt flg) 
{
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      i, j, k, sb;
  PetscInt      lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscReal epsilon=-1.e-8;
  PetscReal ***nvert, ibmval,ibm_Flux,ibm_Area;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area, libm_Flux_abs, ibm_Flux_abs;

  for (sb=1; sb<block_number; sb++) {
  ibmval=sb*10.-1.;

  libm_Flux = 0;
  libm_area = 0;
  libm_Flux_abs=0.;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.5 && nvert[k][j][i+1] < ibmval+0.5 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux += ucor[k][j][i].x;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							  csi[k][j][i].y * csi[k][j][i].y +
							  csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
	    } else 
	      ucor[k][j][i].x=0.;
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.5 && nvert[k][j+1][i] < ibmval+0.5 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux += ucor[k][j][i].y;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/sqrt(eta[k][j][i].x * eta[k][j][i].x +
							  eta[k][j][i].y * eta[k][j][i].y +
							  eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    } else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] > ibmval-0.5 && nvert[k+1][j][i] < ibmval+0.5 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux += ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux -= ucor[k][j][i].x;
	    if (flg==3)
	    libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							csi[k][j][i].y * csi[k][j][i].y +
							csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	    }else 
	      ucor[k][j][i].x=0.;
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux -= ucor[k][j][i].y;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/ sqrt(eta[k][j][i].x * eta[k][j][i].x +
							   eta[k][j][i].y * eta[k][j][i].y +
							   eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    }else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux -= ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux,&ibm_Flux,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_Flux_abs,&ibm_Flux_abs,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area,&ibm_Area,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 /*  PetscGlobalSum(PETSC_COMM_WORLD,&libm_Flux, &ibm_Flux); */
/*   PetscGlobalSum(PETSC_COMM_WORLD,&libm_Flux_abs, &ibm_Flux_abs); */
/*   PetscGlobalSum(PETSC_COMM_WORLD,&libm_area, &ibm_Area); */

  PetscPrintf(PETSC_COMM_WORLD, "BLANKFlux %le %le\n", ibm_Flux, ibm_Area);

  PetscReal correction;

  if (ibm_Area > 1.e-15) {
    if (flg>1) 
      correction = ibm_Flux / ibm_Flux_abs;
    else if (flg)
      correction = (ibm_Flux + FluxInSum) / ibm_Area;
    else
      correction = ibm_Flux / ibm_Area;
  }
  else {
    correction = 0;
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.5 && nvert[k][j][i+1] < ibmval+0.5 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon){
	    if (flg==3) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.5 && nvert[k][j+1][i] < ibmval+0.5 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x + 
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] > ibmval-0.5 && nvert[k+1][j][i] < ibmval+0.5 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x +
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

      }
    }
  }
  
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.5 && nvert[k][j][i+1] < ibmval+0.5 && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.5 && nvert[k][j+1][i] < ibmval+0.5 && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > ibmval-0.5 && nvert[k+1][j][i] < ibmval+0.5 && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux,&ibm_Flux,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area,&ibm_Area,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 /*  PetscGlobalSum(PETSC_COMM_WORLD,&libm_Flux, &ibm_Flux); */
/*   PetscGlobalSum(PETSC_COMM_WORLD,&libm_area, &ibm_Area); */
  PetscPrintf(PETSC_COMM_WORLD, "BLANK Flux22 %le %le\n", ibm_Flux, ibm_Area);

  //  FormBCS(user);        
  } //sb

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  //DALocalToGlobal(user->fda, user->lUcont,INSERT_VALUES,user->Ucont);
 
  
 return 0;
}

/*
 intcontrol
 0 single inlet or outlet
 1 multiple inlet
 2 multiple outlet
 3 summation zero
*/
PetscErrorCode Block_Blank_Correction(UserCtx *user) {
  PetscInt      bi,sb;
  DM            da, fda;
  PetscInt      i, j, k;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts	***ucont,***csi, ***eta, ***zet;
  PetscScalar	FluxIn,ratio;
  PetscScalar	lFluxOut[6];//FluxOutKmax,FluxOutKmin,FluxOutJmax,FluxOutJmin,FluxOutImax,FluxOutImin;
  PetscScalar	lArea[6];//lAreaKmax,lAreaKmin,lAreaJmax,lAreaJmin,lAreaImax,lAreaImin;
  PetscScalar	lFluxOutSum[6]; //FluxOutKmaxSum=0.,FluxOutKminSum=0.,FluxOutJmaxSum=0.,FluxOutJminSum=0.,FluxOutImaxSum=0.,FluxOutIminSum=0.;
  PetscScalar	lAreaSum[6];//AreaKmaxSum=0.,AreaKminSum=0.,AreaJmaxSum=0.,AreaJminSum=0.,AreaImaxSum=0.,AreaIminSum=0.;
  PetscScalar   Area, AreaSum;
  PetscScalar   epsilon=-1.e-8;

  PetscInt ip,jp,kp;
  PetscInt dispx, dispy, dispz, dispnn;
  /**************************************************************************************************************************/
  /* Calculate the fluxIn*/
  /**************************************************************************************************************************/
  //  for (bi=0; bi<block_number; bi++) { 
  bi=user->_this;
  InflowFlux(&(user[bi]));
  //  }

  /**************************************************************************************************************************/
  /* Calculate fluxes and Correct fluxes at Interface
   ----currently only at k direction. easily extendable to other
   ----directions */
  /**************************************************************************************************************************/
 
    //  for (bi=0; bi<block_number; bi++) { 

  da = user[bi].da;
  fda = user[bi].fda;
  info = user[bi].info;
  
  xs = info.xs; xe = info.xs + info.xm;
  ys = info.ys; ye = info.ys + info.ym;
  zs = info.zs; ze = info.zs + info.zm;
  mx = info.mx; my = info.my; mz = info.mz;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  
  DMDAVecGetArray(fda, user[bi].lCsi,  &csi);
  DMDAVecGetArray(fda, user[bi].lEta,  &eta);
  DMDAVecGetArray(fda, user[bi].lZet,  &zet);
  
  //  DMDAVecGetArray(da, user[bi].Nvert,  &nvert);
  
  for (sb=1; sb<block_number; sb++) { 
    ip=user[bi].ip[sb]; jp=user[bi].jp[sb]; kp=user[bi].kp[sb];
    dispx=user[bi].dispx[sb]; dispy=user[bi].dispy[sb]; dispz=user[bi].dispz[sb];
    dispnn=user[bi].dispnn[sb]; 
  /**************************************************************************************************************************/
  /* Calculate the fluxOut and Area on all interfaces*/
  /**************************************************************************************************************************/
    for (i=0; i<6; i++) {
      lFluxOutSum[i]=0.;
      lAreaSum[i]=0.;
    }

    /**************************************************************************************************************************/
    /* Interface at Kmax */
    /**************************************************************************************************************************/

    DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
    //    DMDAVecGetArray(fda, user[bi].Bcs.Ubcs, &ubcs);
      
    lArea[5]=0.;
    lFluxOut[5] = 0.;
    k=kp+dispz;
    if (k>=lzs && k<lze){ 

      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){ 

	  if (i>=lxs && i<lxe && j>=lys && j<lye){
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon){
	      lFluxOut[5] -= ucont[k][j][i].z;
	      lArea[5] += Area;
	    }
	  }

	}
      }
    }
    
    MPI_Allreduce(&lFluxOut[5],&lFluxOutSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[5],&lAreaSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    //  PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[5], &lFluxOutSum[5]);
    //  PetscGlobalSum(PETSC_COMM_WORLD,&lArea[5], &lAreaSum[5]);


    /**************************************************************************************************************************/
    /* Interface at Kmin */
    /**************************************************************************************************************************/
    lArea[4]=0.;
    lFluxOut[4] = 0.;
    k=kp-dispz;
    if (k>=lzs && k<lze){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  if (i>=lxs && i<lxe && j>=lys && j<lye){	
            Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
                         (zet[k][j][i].y) * (zet[k][j][i].y) +
                         (zet[k][j][i].z) * (zet[k][j][i].z));
            if (fabs(ucont[k][j][i].z/Area)>epsilon){
              lFluxOut[4] += ucont[k][j][i].z; 
              lArea[4] += Area; 
	    }  
	  }      
	} 
      } 
    } 

    MPI_Allreduce(&lFluxOut[4],&lFluxOutSum[4],1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[4],&lAreaSum[4],1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);
 
    // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[4], &lFluxOutSum[4]);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[4], &lAreaSum[4]);
    
    /**************************************************************************************************************************/
    /* Interface at Jmax */
    /**************************************************************************************************************************/
    lArea[3]=0.;
    lFluxOut[3] = 0.;
    j=jp+dispy;
    if (j>=lys && j<lye){
      //for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){	  
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[3]==5)
	      ucont[k][j][i].y = 0.;
            if (fabs(ucont[k][j][i].y/Area)>epsilon){
              lFluxOut[3] -= ucont[k][j][i].y;
              lArea[3] += Area;
	    }	
	  }
	}
      }
    }

    MPI_Allreduce(&lFluxOut[3],&lFluxOutSum[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[3],&lAreaSum[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[3], &lFluxOutSum[3]);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[3], &lAreaSum[3]);


    /**************************************************************************************************************************/
    /* Interface at Jmin */
    /**************************************************************************************************************************/
    lArea[2]=0.;
    lFluxOut[2] = 0.;
    j=jp-dispy;
    if (j>=lys && j<lye){ 
      //for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[2]==5)
	      ucont[k][j][i].y = 0.;
            if (fabs(ucont[k][j][i].y/Area)>epsilon){
              lFluxOut[2] += ucont[k][j][i].y;
              lArea[2] += Area;
	    }
	  }
	} 
      }
    }

    MPI_Allreduce(&lFluxOut[2],&lFluxOutSum[2],1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[2],&lAreaSum[2],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    //  PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[2], &lFluxOutSum[2]);
    //  PetscGlobalSum(PETSC_COMM_WORLD,&lArea[2], &lAreaSum[2]);
    
    /**************************************************************************************************************************/
    /* Interface at Imax */
    /**************************************************************************************************************************/
    lArea[1]=0.;
    lFluxOut[1] = 0.;
    i=ip+dispx;
    if (i>=lxs && i<lxe){
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
            Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
                         (csi[k][j][i].y) * (csi[k][j][i].y) +
                         (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[1]==5)
	      ucont[k][j][i].x = 0.;
            if (fabs(ucont[k][j][i].x/Area)>epsilon){
              lFluxOut[1] -= ucont[k][j][i].x;
              lArea[1] += Area;
	    }
	  }
	}
      }
    }

    MPI_Allreduce(&lFluxOut[1],&lFluxOutSum[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[1],&lAreaSum[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    //   PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[1], &lFluxOutSum[1]);
    //   PetscGlobalSum(PETSC_COMM_WORLD,&lArea[1], &lAreaSum[1]);

    /**************************************************************************************************************************/
    /* Interface at Imin */
    /**************************************************************************************************************************/
    lArea[0]=0.;
    lFluxOut[0] = 0.;
    i=ip-dispx;
    if (i>=lxs && i<lxe){
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[0]==5)
	      ucont[k][j][i].x = 0.;
	    if (fabs(ucont[k][j][i].x/Area)>epsilon){
	      lFluxOut[0] += ucont[k][j][i].x;
	      lArea[0] += Area;
	    }
	  }
	}
      }
    }

    MPI_Allreduce(&lFluxOut[0],&lFluxOutSum[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[0],&lAreaSum[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[0], &lFluxOutSum[0]);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[0], &lAreaSum[0]);

  /**************************************************************************************************************************/
  /* Correct fluxes at Interface
   ----Extended to all directions from the previous version
   ---- */
  /**************************************************************************************************************************/
      
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[5]==0) {
      FluxOutSum = lFluxOutSum[5];
      AreaSum = lAreaSum[5];
    } else  if (user[bi].inttype[5]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[5]==2) {
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[5]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[5]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }
    
    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Kmax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[5], lAreaSum[5]);

    k=kp+dispz;
    if (k>=lzs && k<lze){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  if (i>=lxs && i<lxe && j>=lys && j<lye){
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon){
	      ucont[k][j][i].z -= ratio* Area;
	    } else
	      ucont[k][j][i].z = 0.;
	  }
	}
      }
    }
    //////////// Kmax End
    
    //////////// Kmin begin    
    FluxIn = FluxInSum;

    if (user[bi].inttype[4]==0) {
      FluxIn = -FluxInSum;
      FluxOutSum = lFluxOutSum[4];
      AreaSum = lAreaSum[4];
    } else  if (user[bi].inttype[4]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[4]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[4]==3) {
      FluxIn = 0.;

      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[4]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }


    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else      
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Kmin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[4], lAreaSum[4]);
    k=kp-dispz;
    if (k>=lzs && k<lze){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  if (i>=lxs && i<lxe && j>=lys && j<lye){
            Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
                         (zet[k][j][i].y) * (zet[k][j][i].y) +
                         (zet[k][j][i].z) * (zet[k][j][i].z));
            if (fabs(ucont[k][j][i].z/Area)>epsilon){
              ucont[k][j][i].z += ratio* Area;
            } else
              ucont[k][j][i].z = 0.;
	  }
	}
      } 
    }
      
    // kmin 

    // jmax
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[3]==0) {
      FluxOutSum = lFluxOutSum[3];
      AreaSum = lAreaSum[3];
    } else  if (user[bi].inttype[3]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[3]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[3]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[3]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } 


    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Jmax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[3], lAreaSum[3]);
    j=jp+dispy;
    if (j>=lys && j<lye){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[3]==5){
	      ucont[k][j][i].y = 0.;
	    } else if (fabs(ucont[k][j][i].y/Area)>epsilon){
              ucont[k][j][i].y -= ratio * Area;
            } else
              ucont[k][j][i].y = 0.;    
	  }
	}
      }
    }
    // jmax

    // jmin
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[2]==0) {
      FluxIn = -FluxInSum;
      FluxOutSum = lFluxOutSum[2];
      AreaSum = lAreaSum[2];
    } else  if (user[bi].inttype[2]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[2]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[2]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[2]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }
    
    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Jmin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[2], lAreaSum[2]);
    j=jp-dispy;
    if (j>=lys && j<lye){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[2]==5){
	      ucont[k][j][i].y = 0.;
            } else if (fabs(ucont[k][j][i].y/Area)>epsilon){
              ucont[k][j][i].y += ratio * Area;
            } else
              ucont[k][j][i].y = 0.;
	  }
	}
      }
    }    
    // jmin

    // imax
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[1]==0) {
      FluxOutSum = lFluxOutSum[1];
      AreaSum = lAreaSum[1];
    } else  if (user[bi].inttype[1]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[1]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[1]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[1]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }
    
    
    if (AreaSum<epsilon) AreaSum=1.;

    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
 
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Imax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[1], lAreaSum[1]);
    i=ip+dispx;
    if (i>=lxs && i<lxe){ 
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
            Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
                         (csi[k][j][i].y) * (csi[k][j][i].y) +
                         (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[1]==5){
	      ucont[k][j][i].x = 0.;
/* 	    }else if (ratio*lFluxOutSum[1]< 0 && */
/* 		      fabs(ratio*lAreaSum[1])>fabs(lFluxOutSum[1])){ */
/* 	      ucont[k][j][i].x = 0.; */
	    }
            else if (fabs(ucont[k][j][i].x/Area)>epsilon){
              ucont[k][j][i].x -= ratio * Area;
	    } else
	      ucont[k][j][i].x = 0.;
	  }
	}
      }
    }    
    // imax 

    // imin
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[0]==0) {
      FluxIn = -FluxInSum;
      FluxOutSum = lFluxOutSum[0];
      AreaSum = lAreaSum[0];
    } else  if (user[bi].inttype[0]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[0]==2) {
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[0]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[0]==4) {
      FluxIn = user[bi].FluxInSum;
	FluxOutSum=0.;AreaSum=0.;
	for (i=0; i<6; i++) {
	  FluxOutSum += lFluxOutSum[i];
	  AreaSum += lAreaSum[i];
	}
    }
    
    
    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Imin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[0], lAreaSum[0]);
    
    i=ip-dispx;
    if (i>=lxs && i<lxe){
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
            Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
                         (csi[k][j][i].y) * (csi[k][j][i].y) +
                         (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[0]==5){
	      ucont[k][j][i].x = 0.;
/* 	    } else if (ratio*lFluxOutSum[0]< 0 && */
/* 		      fabs(ratio*lAreaSum[0])>fabs(lFluxOutSum[0])){ */
/* 	      ucont[k][j][i].x = 0.; */
	    }
	    else if (fabs(ucont[k][j][i].x/Area)>epsilon){
              ucont[k][j][i].x += ratio * Area;
	    } else
	      ucont[k][j][i].x = 0.;
	  }
	}
        }
    }
    
    // imin 

    /**************************************************************************************************************************/
    /* Interface at  */
    /**************************************************************************************************************************/

    DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
    //    DMDAVecRestoreArray(fda, user[bi].Bcs.Ubcs, &ubcs);
/*     DALocalToGlobalBegin(fda, user[bi].lUcont,  user[bi].Ucont); */
/*     DALocalToGlobalEnd(fda, user[bi].lUcont,  user[bi].Ucont); */
    DMGlobalToLocalBegin(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    
    } //for sb

    DMDAVecRestoreArray(fda, user[bi].lCsi,  &csi);
    DMDAVecRestoreArray(fda, user[bi].lEta,  &eta);
    DMDAVecRestoreArray(fda, user[bi].lZet,  &zet);

    //    DMDAVecRestoreArray(da, user[bi].Nvert,  &nvert);

    Contra2Cart(&(user[bi]));

    //    GhostNodeVelocity(&user[bi]);    

    //    FormBCS(&(user[bi]));        
    // } //bi
 
  return(0);
}

PetscErrorCode Block_Interface_Correction(UserCtx *user) {
  PetscInt      bi;
  DM            da, fda;
  PetscInt      i, j, k;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***ucont,***csi, ***eta, ***zet,***ucat,***ubcs;
  PetscReal     ***nvert;
  PetscScalar	FluxIn, FluxOut,ratio;
  PetscScalar	Area=0.0,lArea[6],AreaSum[6];
  PetscScalar   lFlux[6],FluxSum[6];
  PetscScalar   epsilon=1.e-8;
  PetscInt      ttemp;
  /**************************************************************************************************************************/
  /* Calculate the fluxIn*/
  /**************************************************************************************************************************/
  //  for (bi=0; bi<block_number; bi++) {
  //    InflowFlux(&(user[bi]));
  // }

  /**************************************************************************************************************************/
  /* Calculate fluxes and Correct fluxes at Interface
   ----currently only at k direction. easily extendable to other
   ----directions */
  /**************************************************************************************************************************/

   if(ti<2) epsilon=(-epsilon);
   
/*    if(ti>200&&ti<750){  */
/*      epsilon1=.05; */
/*    }else if(ti>=750&&ti<930){ */
/*      epsilon1=.01; */
/*    }else if(ti>=930){ */
/*      epsilon1=.5; */
/*    } */



  for (bi=0; bi<block_number; bi++) {

    InflowFlux(&(user[bi]));  //InflowFlux should be specified at each block ( Hafez )
    OutflowFlux(&(user[bi]));

    da = user[bi].da;
    fda = user[bi].fda;
    info = user[bi].info;

    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;
    
    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

 
    DMDAVecGetArray(fda, user[bi].lCsi,  &csi);
    DMDAVecGetArray(fda, user[bi].lEta,  &eta);
    DMDAVecGetArray(fda, user[bi].lZet,  &zet);

    DMDAVecGetArray(da, user[bi].Nvert,  &nvert);
    //////////////////////////////////
     DMDAVecGetArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    /////////////////////////////////

  /**************************************************************************************************************************/
  /* Calculate the fluxOut and Area on all interfaces*/
  /**************************************************************************************************************************/
    for (i=0; i<6; i++) {
      FluxSum[i]=0.;
      AreaSum[i]=0.;
      lArea[i]=0.;
      lFlux[i]=0.;
    }

    /**************************************************************************************************************************/
    /* Interface at Kmax */
    /**************************************************************************************************************************/
   

    DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
   
    if (user[bi].bctype[5]==0) {

      if (ze == mz) {
	k=ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>(epsilon) && nvert[k][j][i]<0.1){
		lFlux[5] += ucont[k][j][i].z;
		lArea[5] += Area;
	    }else{
	      
	    }
	  }
 
	}
      }
    }
   /**************************************************************************************************************************/
    /* Interface at Kmin */
    /**************************************************************************************************************************/

    if (user[bi].bctype[4]==0) {
      if (zs == 0) {
	k=0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>(epsilon) && nvert[k+1][j][i]<0.1){
	      lFlux[4] -= ucont[k][j][i].z;
	      lArea[4] += Area;
	    }
	  }
	  
	}
      }
    }
 
    /**************************************************************************************************************************/
    /* Interface at Jmax */
    /**************************************************************************************************************************/

    if (user[bi].bctype[3]==0) {
     
      if (ye == my) {
	j=ye-2;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon && nvert[k][j][i]<0.1){
	      lFlux[3] += ucont[k][j][i].y;
	      lArea[3] += Area;
	      }
	    }

	  }
	}
    }
 
    /**************************************************************************************************************************/
    /* Interface at Jmin */
    /**************************************************************************************************************************/

    if (user[bi].bctype[2]==0) {
     
      if (ys == 0) {
	j=0;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon && nvert[k][j+1][i]<0.1){
	      lFlux[2] -= ucont[k][j][i].y;
	       lArea[2] += Area;
	      }
	    }
	  }
	}
    }
   /**************************************************************************************************************************/
    /* Interface at Imax */
    /**************************************************************************************************************************/

    if (user[bi].bctype[1]==0) {
     
      if (xe == mx) {
	i=xe-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon && nvert[k][j][i]<0.1){
	      lFlux[1] += ucont[k][j][i].x;
	      lArea[1] += Area;
	      }
	    }
 
	  }
	}
    }
     /**************************************************************************************************************************/
    /* Interface at Imin */
    /**************************************************************************************************************************/

    if (user[bi].bctype[0]==0) {
     
      if (xs == 0) {
	i=0;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon && nvert[k][j][i+1]<0.1){
	      lFlux[0] -= ucont[k][j][i].x;
	       lArea[0] += Area;
	      }
	  }

	}
      }
    }
    for (i=0;i<6;i++){
      if (user[bi].bctype[i]==0) {
	MPI_Allreduce(&lFlux[i],&FluxSum[i],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	MPI_Allreduce(&lArea[i],&AreaSum[i],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      }
    }
    
  /**************************************************************************************************************************/
  /* Correct fluxes at Interface
   ----Extended to all directions from the previous version
   ---- */
  /**************************************************************************************************************************/
    PetscScalar ZFlux=0.,ZA=0.;
    for (i=0; i<6; i++) {
      ZFlux+=FluxSum[i]; //Flux Summation on whole control volume of each block(aim is to make this amountplus fluxes on boundaries zero)
      ZA+=AreaSum[i]; // Summation of Area's which need to be corrected (here is interfaces since physical boundary's are corrected before)
    }
 
    PetscInt flag=0;ratio=0.0;
    
    PetscInt period;
    PetscOptionsGetInt(PETSC_NULL, "-period", &period, PETSC_NULL);
    PetscReal LVflux=0;
    PetscReal psys=0.273;

    //////////////////////////////////////////////////////////////////////////////////
    if (((ti-ti/period*period)<psys*period)){
      LVflux= user[bi].FluxIntpSum;
    }else if (((ti-ti/period*period)<(psys*period+557))&&((ti-ti/period*period)>=psys*period)){
      if((ti-ti/period*period)<1851){
	LVflux=1./2*user[bi].FluxIntpSum;
      }else{
	LVflux=1./2*user[bi].FluxIntpSum+((ti-ti/period*period)-1851)*.0023*(1-2.*bi);
	PetscPrintf(PETSC_COMM_SELF, "bi:%d, LVFLUX:%le \n",bi ,LVflux);
      }
    }else if ((ti - ti/period*period) >= 1922 && (ti - ti/period*period) < 4980){
      LVflux=-0.032*(1-2.*bi);
    }else{
      LVflux = (- 0.032 + (0.032/20)*((ti-ti/period*period) - 4980))*(1 - 2.*bi);
    }
    //////////////////////////////////////////////////////////////////////////////////

    ratio = -(ZFlux+FluxInSumB[bi]+user[bi].FluxOutSum-LVflux) / (ZA);
    
    flag=1;
    
  
    PetscScalar Fcheck=0.0,Fchecksum=0.0,Fcheckout=0.0,Fcheckint=0.0,Fcheckoutsum=0.0,Fcheckintsum=0.0;
    if (user[bi].bctype[5]==0 ) {
      

      if (ze==mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>(epsilon) && nvert[k][j][i]<0.1){
	      ucont[k][j][i].z += ratio* Area;
	      
	      //	if (((ti-ti/period*period)>=psys*period)){
	      //	  ucont[k][j][i].z = 0.0;
	      //	}
	      // ucont[k][j][i].z *= ratio;
	      if(j==13 && i==10){
		PetscPrintf(PETSC_COMM_SELF, "Interface: U[%d][%d][%d]={%le,%le,%le} \n",ze-2,13,10,ucont[ze-2][13][10].x,ucont[ze-2][13][10].y,ucont[ze-2][13][10].z);
	      }
	      Fcheck+=ucont[k][j][i].z;
	      if(user[bi].bctype[5]==4){
		Fcheckout+=ucont[k][j][i].z;
	      }else{
		Fcheckint+=ucont[k][j][i].z;
	      }


	    } else {
	      ucont[k][j][i].z=0.0;
	    }
	  }
	  
	}
      }
      
    }// kmax bctype[5]== INTERFACE
    
    if (user[bi].bctype[4]==0 ) {
         
      
      if (zs == 0) {
	k = 0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>(epsilon) && nvert[k+1][j][i]<0.1){
	        ucont[k][j][i].z -= ratio* Area;

		//	if (((ti-ti/period*period)>=psys*period)){
		//	  ucont[k][j][i].z = 0.0;
		//	}
	      // ucont[k][j][i].z *= ratio;
	      Fcheck-=ucont[k][j][i].z;
	      if(user[bi].bctype[4]==4){
		Fcheckout-=ucont[k][j][i].z;
	      }else{
		Fcheckint-=ucont[k][j][i].z;
	      }

	    } else {
	      ucont[k][j][i].z=0.0;
	    }
	  }
	  
	}
      }
    } // kmin bctype[4]== INTERFACE

    if (user[bi].bctype[3]==0 ) {
     
      if (ye==my) {
	j = ye-2;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon && nvert[k][j][i]<0.1){
	       ucont[k][j][i].y += ratio* Area;

	       //    if (((ti-ti/period*period)>=psys*period)){
	       //	 ucont[k][j][i].y = 0.0;
	       //   }
	      //  ucont[k][j][i].y *= ratio;
	      Fcheck+=ucont[k][j][i].y;
	      if(user[bi].bctype[3]==4){
		Fcheckout+=ucont[k][j][i].y;
	      }else{
		Fcheckint+=ucont[k][j][i].y;
	      }
	    } else {
	      ucont[k][j][i].y=0.0;
	    }
	  }
	   
	}
      }
     
    } // jmax bctype[3]== INTERFACE

    if (user[bi].bctype[2]==0 ) {
   
      if (ys == 0) {
	j = 0;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon && nvert[k][j+1][i]<0.1){
	      ucont[k][j][i].y -= ratio* Area;

	       //   if (((ti-ti/period*period)>=psys*period)){
	       //	 ucont[k][j][i].y = 0.0;
	       //   }
	       //   ucont[k][j][i].y *= ratio;
	      Fcheck-=ucont[k][j][i].y;
	      if(user[bi].bctype[2]==4){
		Fcheckout-=ucont[k][j][i].y;
	      }else{
		Fcheckint-=ucont[k][j][i].y;
	      }
	    } else {
	      ucont[k][j][i].y=0.0;
	         }
	  }
	   
	}
      }
      
    } // jmin bctype[2]== INTERFACE


    if (user[bi].bctype[1]==0 ) {
       
      
      if (xe == mx) {
	i=xe-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon && nvert[k][j][i]<0.1){
	       ucont[k][j][i].x += ratio* Area;

	      //  if (((ti-ti/period*period)>=psys*period)){
	      //	ucont[k][j][i].x = 0.0;
	      //  }
	      // ucont[k][j][i].x *= ratio;
	      Fcheck+=ucont[k][j][i].x;
	      if(user[bi].bctype[1]==4){
		Fcheckout+=ucont[k][j][i].x;
	      }else{
		Fcheckint+=ucont[k][j][i].x;
	      }
	    } else {
	      ucont[k][j][i].x=0.0;
	    }
	  }
	  
	}
      }
      
    } // imax bctype[1]== INTERFACE

    if (user[bi].bctype[0]==0 ) {
         

      if (xs == 0) {
	i=0;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon && nvert[k][j][i+1]<0.1){
	      ucont[k][j][i].x -= ratio* Area;
	       //   if (((ti-ti/period*period)>=psys*period)){
	       //	 ucont[k][j][i].x = 0.0;
	       //   }
	       //  ucont[k][j][i].x *= ratio;
	      Fcheck-=ucont[k][j][i].x;
              if(user[bi].bctype[0]==4){
		Fcheckout-=ucont[k][j][i].x;
	      }else{
		Fcheckint-=ucont[k][j][i].x;
	      }
	    } else {
	      ucont[k][j][i].x=0.0;
	       }
	  }
	  
	}
      }
      
    } // imin bctype[0]== INTERFACE

   
    MPI_Allreduce(&Fcheck,&Fchecksum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&Fcheckint,&Fcheckintsum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&Fcheckout,&Fcheckoutsum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   
    if (flag){
      //  PetscPrintf(PETSC_COMM_WORLD, "Block:%d: ratio=%le,FInPhysicalBoundary= %le FOutPhysicalBoundaryPlusInterfaces=%le , FluxSum=%le \n", bi, ratio,FluxInSumB[bi], Fchecksum,FluxInSumB[bi]+Fchecksum);
      PetscPrintf(PETSC_COMM_WORLD, "Block:%d:Correction=%le,ratio:%le  FIn=%le FOut= %le->%le FInterfaces=%le->%le \n", bi,ZA,ratio,(ZFlux+FluxInSumB[bi]+user[bi].FluxOutSum-LVflux),Fcheckoutsum,LVflux,Fcheckintsum);
      
 PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------------------------------------------------------- \n");
       }





    DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
    DMDAVecRestoreArray(fda, user[bi].lCsi,  &csi);
    DMDAVecRestoreArray(fda, user[bi].lEta,  &eta);
    DMDAVecRestoreArray(fda, user[bi].lZet,  &zet);
    /////////////////////////
   DMDAVecRestoreArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    //////////////////////////

    DMDAVecRestoreArray(da, user[bi].Nvert,  &nvert);

    DMGlobalToLocalBegin(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    Contra2Cart(&(user[bi]));

    GhostNodeVelocity(&user[bi]);

          
  }//bi



  return(0);
}


PetscErrorCode Block_Interface_U(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	hostU;
  Cmpnts ***itfc, ***ucat,***ubcs, ucent;
  PetscReal *hostu,***nvert;

  VecScatter tolocalall;
  PetscErrorCode ierr;
  PetscLogDouble size1=0., size2=0., size3=0.;
    
  /* ucat is at cell centers while litfc is now on the cell corners! */
  
  for (bi=0; bi<block_number; bi++) {
    /* hostU is a parallel PETSc vector that will hold vector values 
       in the natural numbering, rather than in the PETSc parallel 
       numbering associated with the DA */
    ierr = DMDACreateNaturalVector(user[bi].fda, &hostU);
    

    // put the center node velocities in the natural ordering in nhostU 
    DMDAGlobalToNaturalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);
    DMDAGlobalToNaturalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);

    /* the velocities at cell centers are saved in user[bi].nhostU
       which is a sequential vector*/
   
    VecScatterCreateToAll(hostU, &tolocalall, &(user[bi].nhostU));

    VecScatterBegin(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		    SCATTER_FORWARD);
    VecScatterEnd(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		  SCATTER_FORWARD);
   
    VecScatterDestroy(&tolocalall);
    VecDestroy(&hostU);
   
  }

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt    lmx, lmy, lmz;
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes
       from the host nodes.
       itfc is the velocity at the cell corners
       hostU is the velocity at the cell centers of the host block */
    /**************************************************************************************************************************/
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
   
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];

      if ((x+y+z<-1.) &&
	  ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {
	itfc[ck][cj][ci].x = 0.;
	itfc[ck][cj][ci].y = 0.;
	itfc[ck][cj][ci].z = 0.;

      } else
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {

	VecGetArray(user[hb].nhostU, &hostu);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].z = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * z); //i+1,j+1,k+1

	VecRestoreArray(user[hb].nhostU, &hostu);

       /*  if(hb==0&&hk>120){ */
       /* 	  itfc[ck][cj][ci].x=0.0; */
       /* 	  itfc[ck][cj][ci].y=0.0; */
       /* 	  itfc[ck][cj][ci].z=0.0; */
       /* } */

      }
    } // for itfcnumber
  
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDALocalToLocalBegin(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
			user[bi].lItfc);
    DMDALocalToLocalEnd(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
		      user[bi].lItfc);
  } // for bi
  //  PetscBarrier(PETSC_NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");

  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
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

    Cmpnts ***ucont, ***kzet, ***jeta,***icsi;
    PetscReal ucx, ucy, ucz;
    //    DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecGetArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecGetArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecGetArray(user[bi].fda, user[bi].lCsi, &icsi);
    DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert);

    /**************************************************************************************************************************/
    /* Create boundary condition for flux (ucont) and cell surface
       center velocity (ubcs) from the interpolated velocity (itfc).
       
       itfc is the velocity at the cell surface centers */
    /**************************************************************************************************************************/
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  
	  if(nvert[k][j][i+1]<0.1){
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	 
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j-1][i  ].x) ;
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i].x*
			      icsi[k][j][i].x +
			      ubcs[k][j][i].y*
			      icsi[k][j][i].y +
			      ubcs[k][j][i].z*
			      icsi[k][j][i].z);
	  }
	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  
	  if(nvert[k][j][i]<0.1){
/* 	  ubcs[k][j][i+1].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i+1].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i+1].z = itfc[k][j][i].z; */
	  ubcs[k][j][i+1].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j-1][i  ].x);
	  ubcs[k][j][i+1].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i+1].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i+1].x  *
			      icsi[k][j][i].x +
			      ubcs[k][j][i+1].y  *
			      icsi[k][j][i].y +
			      ubcs[k][j][i+1].z *
			      icsi[k][j][i].z);
	  }

	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  
	  if(nvert[k][j+1][i]<0.1){
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j  ][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j  ][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j][i].x*
			       jeta[k][j][i].x +
			       ubcs[k][j][i].y *
			       jeta[k][j][i].y +
			       ubcs[k][j][i].z *
			       jeta[k][j][i].z);
	  }
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {

	  if(nvert[k][j][i]<0.1){
/* 	  ubcs[k][j+1][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j+1][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j+1][i].z = itfc[k][j][i].z; */
	  ubcs[k][j+1][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j  ][i-1].x);
	  ubcs[k][j+1][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j  ][i-1].y);
	  ubcs[k][j+1][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j+1][i].x*
			       jeta[k][j  ][i].x +
			       ubcs[k][j+1][i].y*
			       jeta[k][j  ][i].y +
			       ubcs[k][j+1][i].z*
			       jeta[k][j  ][i].z);
	 /*  if (i==1 && k==50)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==100)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==150)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==200)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==50)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==100)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==150)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==200)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
	  }
	}
      }
    }


    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  
	  if(nvert[k+1][j][i]<0.1){
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j-1][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j-1][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k][j][i].x*
			       kzet[k][j][i].x +
			       ubcs[k][j][i].y*
			       kzet[k][j][i].y +
			       ubcs[k][j][i].z *
			       kzet[k][j][i].z);
	  }
	  
	}
      }
    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  
	  if(nvert[k][j][i]<0.1){
/* 	  ubcs[k+1][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k+1][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k+1][j][i].z = itfc[k][j][i].z; */
	  ubcs[k+1][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j-1][i-1].x);
	  ubcs[k+1][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j-1][i-1].y);
	  ubcs[k+1][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k+1][j][i].x*
			       kzet[k  ][j][i].x +
			       ubcs[k+1][j][i].y*
			       kzet[k  ][j][i].y +
			       ubcs[k+1][j][i].z *
			       kzet[k  ][j][i].z);
	  }

	}
      }
    }

    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    // This part is for blanking
    if(blank && bi==0){
           
      DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecGetArray(user[bi].fda, user[bi].lUcont, &ucont);
      for (sb=1; sb<block_number; sb++) {
	PetscInt ip, im, jp, jm, kp, km;
	PetscInt ii, jj, kk;
	  
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    for (i=lxs; i<lxe; i++) {
	      ip = (i<mx-2?(i+1):(i));
	      im = (i>1   ?(i-1):(i));
	      
	      jp = (j<my-2?(j+1):(j));
	      jm = (j>1   ?(j-1):(j));
	      
	      kp = (k<mz-2?(k+1):(k));
	      km = (k>1   ?(k-1):(k));
	      
	      if (((int)(nvert[k][j][i]+0.5) < sb*10) &&
		  ((int)(nvert[k][j][i]+0.5) > sb*10-3) ) {
		// flux in x direction
		kk=k;	jj=j;
		for (ii=im; ii<ip; ii++) {
		  ucent.x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
				    itfc[kk  ][jj-1][ii  ].x +
				    itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk-1][jj-1][ii  ].x) ;
		  ucent.y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
				    itfc[kk  ][jj-1][ii  ].y +
				    itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk-1][jj-1][ii  ].y);
		  ucent.z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
				    itfc[kk  ][jj-1][ii  ].z +
				    itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk-1][jj-1][ii  ].z);
		  ucont[kk][jj][ii].x = (ucent.x*
					 icsi[kk][jj][ii].x +
					 ucent.y*
					 icsi[kk][jj][ii].y +
					 ucent.z*
					 icsi[kk][jj][ii].z);
		}
		// flux in y direction
		kk=k;   ii=i;
		for (jj=jm; jj<jp; jj++) {
		  ucent.x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
				    itfc[kk  ][jj  ][ii-1].x +
				    itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk-1][jj  ][ii-1].x);
		  ucent.y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
				    itfc[kk  ][jj  ][ii-1].y +
				    itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk-1][jj  ][ii-1].y);
		  ucent.z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
				    itfc[kk  ][jj  ][ii-1].z +
				    itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk-1][jj  ][ii-1].z);
		  ucont[kk][jj][ii].y = ( ucent.x*
					  jeta[kk][jj][ii].x +
					  ucent.y *
					  jeta[kk][jj][ii].y +
					  ucent.z *
					  jeta[kk][jj][ii].z);
		}
		// flux in z direction
		jj=j;  ii=i;
		for (kk=km; kk<kp; kk++) {
		  ucent.x = 0.25 * (itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk  ][jj  ][ii-1].x +
				    itfc[kk  ][jj-1][ii  ].x +
				    itfc[kk  ][jj-1][ii-1].x);
		  ucent.y = 0.25 * (itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk  ][jj  ][ii-1].y +
				    itfc[kk  ][jj-1][ii  ].y +
				    itfc[kk  ][jj-1][ii-1].y);
		  ucent.z = 0.25 * (itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk  ][jj  ][ii-1].z +
				    itfc[kk  ][jj-1][ii  ].z +
				    itfc[kk  ][jj-1][ii-1].z);
		  ucont[kk][jj][ii].z = ( ucent.x*
					  kzet[kk][jj][ii].x +
					  ucent.y*
					  kzet[kk][jj][ii].y +
					  ucent.z *
					  kzet[kk][jj][ii].z);
		}
	      }// if (nvert)
	    }
	  }
	}
      } //for sb
      //DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecRestoreArray(user[bi].fda, user[bi].lUcont, &ucont);
      //      PetscPrintf(PETSC_COMM_WORLD, "Local to global lUcont _U ");

      DMLocalToGlobalBegin(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);
      DMLocalToGlobalEnd(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);

      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    } // if blank
    DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lCsi, &icsi);
 
  } // bi

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostU);
    Contra2Cart(&(user[bi]));
    GhostNodeVelocity(&user[bi]);
  }

  return(0);      
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscErrorCode Block_Interface_Q(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	hostQ;
  PetscReal ***itfcq, ***qnew;
  PetscReal *hostq,***nvert;

  VecScatter tolocalall;
  PetscErrorCode ierr;
  PetscLogDouble size1=0., size2=0., size3=0.;
    
  /* ucat is at cell centers while litfc is now on the cell corners! */
  
  for (bi=0; bi<block_number; bi++) {
    /* hostU is a parallel PETSc vector that will hold vector values 
       in the natural numbering, rather than in the PETSc parallel 
       numbering associated with the DA */
    ierr = DMDACreateNaturalVector(user[bi].da, &hostQ);
    

    // put the center node velocities in the natural ordering in nhostU 
    DMDAGlobalToNaturalBegin(user[bi].da, user[bi].Qnew, INSERT_VALUES, hostQ);
    DMDAGlobalToNaturalEnd(user[bi].da, user[bi].Qnew, INSERT_VALUES, hostQ);

    /* the velocities at cell centers are saved in user[bi].nhostU
       which is a sequential vector*/
   
    VecScatterCreateToAll(hostQ, &tolocalall, &(user[bi].nhostQ));

    VecScatterBegin(tolocalall, hostQ, user[bi].nhostQ, INSERT_VALUES,
		    SCATTER_FORWARD);
    VecScatterEnd(tolocalall, hostQ, user[bi].nhostQ, INSERT_VALUES,
		  SCATTER_FORWARD);
   
    VecScatterDestroy(&tolocalall);
    VecDestroy(&hostQ);
   
  }

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt    lmx, lmy, lmz;
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes
       from the host nodes.
       itfc is the velocity at the cell corners
       hostU is the velocity at the cell centers of the host block */
    /**************************************************************************************************************************/
    DMDAVecGetArray(user[bi].da, user[bi].lItfcQ, &itfcq);
   
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];

      if ((x+y+z<-1.) &&
	  ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {
	itfcq[ck][cj][ci] = 0.;

      } else
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {

	VecGetArray(user[hb].nhostQ, &hostq);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfcq[ck][cj][ci] = (hostq[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostq[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostq[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostq[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostq[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostq[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostq[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostq[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))]
			      * x     * y     * z); //i+1,j+1,k+1


	VecRestoreArray(user[hb].nhostQ, &hostq);
      }
    } // for itfcnumber
  
    DMDAVecRestoreArray(user[bi].da, user[bi].lItfcQ, &itfcq);
    DMDALocalToLocalBegin(user[bi].da, user[bi].lItfcQ, INSERT_VALUES,
			user[bi].lItfcQ);
    DMDALocalToLocalEnd(user[bi].da, user[bi].lItfcQ, INSERT_VALUES,
		      user[bi].lItfcQ);
  } // for bi
  //  PetscBarrier(PETSC_NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");

  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
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

    PetscReal ucx, ucy, ucz;
    //    DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].da, user[bi].lItfcQ, &itfcq);
    DMDAVecGetArray(user[bi].da, user[bi].Qnew, &qnew);
    /**************************************************************************************************************************/
    /* Create boundary condition for flux (ucont) and cell surface
       center velocity (ubcs) from the interpolated velocity (itfc).
       
       itfc is the velocity at the cell surface centers */
    /**************************************************************************************************************************/
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  qnew[k][j][i] = 0.25 * (itfcq[k-1][j  ][i  ] +
				  itfcq[k  ][j-1][i  ] +
				  itfcq[k  ][j  ][i  ] +
				   itfcq[k-1][j-1][i  ]) ;

	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  qnew[k][j][i+1] = 0.25 * (itfcq[k-1][j  ][i  ]+
				     itfcq[k  ][j-1][i  ] +
				     itfcq[k  ][j  ][i  ] +
				     itfcq[k-1][j-1][i  ]);
	  
	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  qnew[k][j][i] = 0.25 * (itfcq[k-1][j  ][i  ] +
				   itfcq[k  ][j  ][i-1] +
				   itfcq[k  ][j  ][i  ] +
				   itfcq[k-1][j  ][i-1]);
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  qnew[k][j+1][i] = 0.25 * (itfcq[k-1][j  ][i  ] +
				     itfcq[k  ][j  ][i-1] +
				     itfcq[k  ][j  ][i  ] +
				     itfcq[k-1][j  ][i-1]);
	}
      }
    }


    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  qnew[k][j][i] = 0.25 * (itfcq[k  ][j  ][i  ] +
				   itfcq[k  ][j  ][i-1] +
				   itfcq[k  ][j-1][i  ] +
				   itfcq[k  ][j-1][i-1]);
	}
      }
    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  qnew[k+1][j][i] = 0.25 * (itfcq[k  ][j  ][i  ] +
				     itfcq[k  ][j  ][i-1] +
				     itfcq[k  ][j-1][i  ] +
				     itfcq[k  ][j-1][i-1]);

	}
      }
    }


    /////////blanking removed///////

   
    DMDAVecRestoreArray(user[bi].da, user[bi].lItfcQ, &itfcq);
    DMDAVecRestoreArray(user[bi].da, user[bi].Qnew, &qnew);
  } // bi

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostQ);
  }

  return(0);      
}


      
PetscErrorCode GhostNodeVelocity(UserCtx *user)
{

  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  //  Vec		Coor;
  Cmpnts	***ubcs, ***ucat, ***csi, ***eta, ***zet,***cent;

  PetscReal Un, nx,ny,nz,A;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);
 
  
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */
  extern PetscInt wallfunction;
  if (wallfunction) {
    PetscReal ***ustar, ***aj, ***nvert;
  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(da, user->lUstar, &ustar);
  DMDAVecGetArray(da, user->lAj,  &aj);
  DMDAVecGetArray(da, user->lNvert,  &nvert);

  // wall function for boundary
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if( nvert[k][j][i]<1.1 && ( (user->bctype[0]==-1 && i==1) || (user->bctype[1]==-1 &&  i==mx-2) ) ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double sb, sc; 
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(i==1) {
	    sc = 2*sb + 0.5/aj[k][j][i+1]/area;
	    Uc = ucat[k][j][i+1];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j][i-1]/area;
	    Uc = ucat[k][j][i-1];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
	  
	  //if(i==1) printf("%f %f, %f %f %f, %e %e %e, %e\n", sc, sb, Ua.x, Ua.y, Ua.z, Uc.x, Uc.y, Uc.z, ucat[k][j][i+2].z);
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];	  
	}
	
	if( nvert[k][j][i]<1.1 && ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 &&  j==my-2) ) ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc; 
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(j==1) {
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j-1][i]/area;
	    Uc = ucat[k][j-1][i];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	  
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];
	}
	
      }
  DMDAVecRestoreArray(da, user->lAj,  &aj);
  DMDAVecRestoreArray(da, user->lNvert,  &nvert);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);   
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  } // if wallfunction

/* ==================================================================================             */
/*   Driven Cavity  */
/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
  if (user->bctype[5]==2) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Cylinder O-grid */
  /*  ==================================================================================== */
  if (user->bctype[3]==12) {
  /* Designed to test O-grid for flow over a cylinder at jmax velocity is 1 (horizontal) 
   u_x = 1 at k==kmax */
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 1.;
	}
      }
    }
  }
/*  ==================================================================================== */
/*     Annulus */
/*  ==================================================================================== */
/* designed to test periodic boundary condition for c grid j=0 rotates */
  DMDAVecGetArray(fda, user->Cent, &cent);
  if (user->bctype[2]==11) {
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);;
	  // ubcs[k][j][i].x=0.0;
	  ubcs[k][j][i].y =-cent[k][j+1][i].x/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  // ubcs[k][j][i].y = -cent[k][j+1][i].z/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  ubcs[k][j][i].z = 0.0;
	  // ubcs[k][j][i].z =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y);
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Cent,  &cent);
/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i].z*csi[k][j][i].z +
	       csi[k][j][i].y*csi[k][j][i].y +
	       csi[k][j][i].x*csi[k][j][i].x);
	nx=csi[k][j][i].x/A;
	ny=csi[k][j][i].y/A;
	nz=csi[k][j][i].z/A;
	Un=ucat[k][j][i+1].x*nx+ucat[k][j][i+1].y*ny+ucat[k][j][i+1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i+1].x-Un*nx;//-V_frame.x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i+1].z-Un*nz;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i-1].z*csi[k][j][i-1].z +
	       csi[k][j][i-1].y*csi[k][j][i-1].y +
	       csi[k][j][i-1].x*csi[k][j][i-1].x);
	nx=csi[k][j][i-1].x/A;
	ny=csi[k][j][i-1].y/A;
	nz=csi[k][j][i-1].z/A;
	Un=ucat[k][j][i-1].x*nx+ucat[k][j][i-1].y*ny+ucat[k][j][i-1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i-1].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i-1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i-1].z-Un*nz;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
	       eta[k][j][i].y*eta[k][j][i].y +
	       eta[k][j][i].x*eta[k][j][i].x);
	nx=eta[k][j][i].x/A;
	ny=eta[k][j][i].y/A;
	nz=eta[k][j][i].z/A;
	Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j+1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j+1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j+1][i].z-Un*nz;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j-1][i].z*eta[k][j-1][i].z +
	       eta[k][j-1][i].y*eta[k][j-1][i].y +
	       eta[k][j-1][i].x*eta[k][j-1][i].x);
	nx=eta[k][j-1][i].x/A;
	ny=eta[k][j-1][i].y/A;
	nz=eta[k][j-1][i].z/A;
	Un=ucat[k][j-1][i].x*nx+ucat[k][j-1][i].y*ny+ucat[k][j-1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j-1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j-1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j-1][i].z-Un*nz;
      }
    }
    }
  }


  if (user->bctype[4]==3) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	A=sqrt(zet[k][j][i].z*zet[k][j][i].z +
	       zet[k][j][i].y*zet[k][j][i].y +
	       zet[k][j][i].x*zet[k][j][i].x);
	nx=zet[k][j][i].x/A;
	ny=zet[k][j][i].y/A;
	nz=zet[k][j][i].z/A;
	Un=ucat[k][j][i].x*nx+ucat[k][j][i].y*ny+ucat[k][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i].z-Un*nz;
	}
      }
    }
  }

  if (user->bctype[5]==3) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	A=sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z +
	       zet[k-1][j][i].y*zet[k-1][j][i].y +
	       zet[k-1][j][i].x*zet[k-1][j][i].x);
	nx=zet[k-1][j][i].x/A;
	ny=zet[k-1][j][i].y/A;
	nz=zet[k-1][j][i].z/A;
	Un=ucat[k-1][j][i].x*nx+ucat[k-1][j][i].y*ny+ucat[k-1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k-1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k-1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k-1][j][i].z-Un*nz;
	}
      }
    }
  }
/*  ==================================================================================== */
  /*     Rheology */
  /*  ==================================================================================== */
  if (user->bctype[2]==13){
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z =-1.;
	}
      }
    }
  }
  if (user->bctype[3]==13){
    if (ye==my){
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 1.;
	}
      }
    }
  }
/* ==================================================================================             */
  // boundary conditions on ghost nodes
/* ==================================================================================             */
//Mohsen Aug 2012
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }

  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }

  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,&ucat);
/* ==================================================================================             */
/*   Periodic BC Mohsen Aug 2012
/* ==================================================================================             */
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 
  DMDAVecGetArray(fda, user->lUcat, &ucat);
   
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  if(k>0 && k<user->KM && j>0 && j<user->JM){
	    ucat[k][j][i].x=ucat[k][j][i-2].x;
	    ucat[k][j][i].y=ucat[k][j][i-2].y;
	    ucat[k][j][i].z=ucat[k][j][i-2].z;
	  }
	}
      }
    }
  }
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  if(k>0 && k<user->KM && i>0 && i<user->IM){
	    ucat[k][j][i].x=ucat[k][j-2][i].x;
	    ucat[k][j][i].y=ucat[k][j-2][i].y;
	    ucat[k][j][i].z=ucat[k][j-2][i].z;
	  }
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  if(j>0 && j<user->JM && i>0 && i<user->IM){
	    ucat[k][j][i].x=ucat[k-2][j][i].x;
	    ucat[k][j][i].y=ucat[k-2][j][i].y;
	    ucat[k][j][i].z=ucat[k-2][j][i].z;
	  }
	}
      }
    }
  }
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xe==mx){
      i=mx-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  if(k>0 && k<user->KM && j>0 && j<user->JM){
	    ucat[k][j][i].x=ucat[k][j][i+2].x;
	    ucat[k][j][i].y=ucat[k][j][i+2].y;
	    ucat[k][j][i].z=ucat[k][j][i+2].z;
	  }
	}
      }
    }
  }
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ye==my){
      j=my-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  if(k>0 && k<user->KM && i>0 && i<user->IM){
	    ucat[k][j][i].x=ucat[k][j+2][i].x;
	    ucat[k][j][i].y=ucat[k][j+2][i].y;
	    ucat[k][j][i].z=ucat[k][j+2][i].z;
	  }
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (ze==mz){
      k=mz-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  if(j>0 && j<user->JM && i>0 && i<user->IM){
	    ucat[k][j][i].x=ucat[k+2][j][i].x;
	    ucat[k][j][i].y=ucat[k+2][j][i].y;
	    ucat[k][j][i].z=ucat[k+2][j][i].z;
	  }
	}
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  // velocity on the corner point
  if (zs==0 ) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
      }
    }
  }

  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i+1].z);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
      }
    }
    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
      }
    }
  }

  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i+1].z);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
      }
    }
  }

  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i+1].z);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

 //i-direction

  DMDAVecGetArray(fda, user->lUcat,  &ucat);
 
  if (user->bctype[0]==7){
    if (xs==0){
      i=xs;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i-2].x;
	  ucat[k][j][i].y=ucat[k][j][i-2].y;
	  ucat[k][j][i].z=ucat[k][j][i-2].z;
	}
      }
    }
  }
  if (user->bctype[1]==7){
    if (xe==mx){
      i=xe-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i+2].x;
	  ucat[k][j][i].y=ucat[k][j][i+2].y;
	  ucat[k][j][i].z=ucat[k][j][i+2].z;
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  //j-direction
  
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
 
  if (user->bctype[2]==7){
    if (ys==0){
      j=ys;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j-2][i].x;
	  ucat[k][j][i].y=ucat[k][j-2][i].y;
	  ucat[k][j][i].z=ucat[k][j-2][i].z;
	}
      }
    }
  }
  
  if (user->bctype[3]==7){
    if (ye==my){
      j=my-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j+2][i].x;
	  ucat[k][j][i].y=ucat[k][j+2][i].y;
	  ucat[k][j][i].z=ucat[k][j+2][i].z;
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  
  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  //k-direction
  
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
 
  if (user->bctype[4]==7){
    if (zs==0){
      k=zs;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k-2][j][i].x;
	  ucat[k][j][i].y=ucat[k-2][j][i].y;
	  ucat[k][j][i].z=ucat[k-2][j][i].z;
	}
      }
    }
  }
  if (user->bctype[5]==7){
    if (ze==mz){
      k=mz-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k+2][j][i].x;
	  ucat[k][j][i].y=ucat[k+2][j][i].y;
	  ucat[k][j][i].z=ucat[k+2][j][i].z;
	}
      }
    }
  }
   
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

 
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  return(0);
  
}

/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation/interface
   -1:  wallfunction
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   3:   slip wall/symmetry
   5:	Inlet
   4:	Outlet
   6:   farfield
   7:   periodic
   8:   Characteristic BC
   9:   Analytical Vortex
   10:  Oulet Junction
   11:  Annulus
   12:  Ogrid
   13:  Rheology
   14:  Outlet with Interface
  
*/

PetscErrorCode FormBCS(UserCtx *user)
{
DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  PetscReal	mat[3][3], det, det0, det1, det2;
  PetscReal	q[3],***nvert; //local working array

  PetscReal	***p;
   Cmpnts	***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet,***ucont_o;
   Cmpnts	***cent,***centx,***centy,***centz;
  PetscScalar	FluxIn, r, uin, ratio, xc, yc;
  PetscReal     FluxOut,FluxOutSumcont,FluxOutSumInside;
  PetscScalar   lArea, AreaSum;
  //  PetscReal     absn, n_x,n_y,n_z;
  PetscScalar   FarFluxIn=0., FarFluxOut=0., FarFluxInSum, FarFluxOutSum;
  PetscScalar   FarAreaIn=0., FarAreaOut=0., FarAreaInSum, FarAreaOutSum;
  PetscScalar   FluxDiff, VelDiffIn, VelDiffOut;
  Cmpnts        V_frame;
  //  PetscInt      moveframe=0;

  Vec		lUcont;
  PetscReal	***lucont;

  PetscReal Un, nx,ny,nz,A;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx ) lxe = xe-1;
  if (ye==my ) lye = ye-1;
  if (ze==mz ) lze = ze-1;


  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);

  DMDAVecGetArray(da, user->Nvert, &nvert);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  PetscInt ttemp;
  for (ttemp=0; ttemp<1; ttemp++) {
   
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
/* ==================================================================================             */
/*   FAR-FIELD BC */
/* ==================================================================================             */
 

  // reset FAR FLUXES
  FarFluxIn = 0.; FarFluxOut=0.;
  FarAreaIn = 0.; FarAreaOut=0.;

    V_frame.x=0.;
    V_frame.y=0.;
    V_frame.z=0.;

 

  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x;
	  ubcs[k][j][i].y = ucat[k][j][i+1].y;
	  ubcs[k][j][i].z = ucat[k][j][i+1].z;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;
	  FarFluxIn += ucont[k][j][i].x;
	  FarAreaIn += csi[k][j][i].x;
	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x;
	  ubcs[k][j][i].y = ucat[k][j][i-1].y;
	  ubcs[k][j][i].z = ucat[k][j][i-1].z;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	  FarFluxOut += ucont[k][j][i-1].x;
	  FarAreaOut += csi[k][j][i-1].x;
	}
      }
    }
  }

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j+1][i].x;
	  ubcs[k][j][i].y = ucat[k][j+1][i].y;
	  ubcs[k][j][i].z = ucat[k][j+1][i].z;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;
	  FarFluxIn += ucont[k][j][i].y;
	  FarAreaIn += eta[k][j][i].y;
	}
      }
    }
  }
  
  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j-1][i].x;
	  ubcs[k][j][i].y = ucat[k][j-1][i].y;
	  ubcs[k][j][i].z = ucat[k][j-1][i].z;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
	  FarFluxOut += ucont[k][j-1][i].y;
	  FarAreaOut += eta[k][j-1][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==6 || user->bctype[4]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;
	  FarFluxIn += ucont[k][j][i].z;
	  FarAreaIn += zet[k][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  FarFluxOut += ucont[k-1][j][i].z;
	  FarAreaOut += zet[k-1][j][i].z;
	}
      }
    }
  }
  
  MPI_Allreduce(&FarFluxIn,&FarFluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarFluxOut,&FarFluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    
  MPI_Allreduce(&FarAreaIn,&FarAreaInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarAreaOut,&FarAreaOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
  if (user->bctype[5]==6 || user->bctype[3]==6 || user->bctype[1]==6) {
    FluxDiff = 0.5*(FarFluxInSum - FarFluxOutSum) ;
    VelDiffIn  = FluxDiff / FarAreaInSum ;
    
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = FluxDiff / FarAreaOutSum ;
  
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    PetscPrintf(PETSC_COMM_WORLD, "Far Flux Diff %d %le %le %le %le %le %le %le\n", ti, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  if (user->bctype[5]==10) {
    FluxDiff = FluxInSum -( FarFluxOutSum -FarFluxInSum) ;
    VelDiffIn  = 1/3.*FluxDiff / (FarAreaInSum);// +  FarAreaOutSum);
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = 2./3.* FluxDiff / (FarAreaOutSum) ;//(FarAreaInSum +  FarAreaOutSum) ;
   
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    PetscPrintf(PETSC_COMM_WORLD, "Far Flux Diff %d %le %le %le %le %le %le %le\n", ti, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  
  // scale global mass conservation

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	  ubcs[k][j][i].z = ucat[k-1][j][i].z + VelDiffOut ;//+ V_frame.z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;


	}
      }
    }
  }

  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].y = ucat[k][j-1][i].y + VelDiffOut;// + V_frame.y;
/* 	  ubcs[k][j][i].y = ucat[k][j-1][i].y - VelDiffIn + V_frame.y; */
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;

	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x + VelDiffOut;// + V_frame.x;
/* 	  ubcs[k][j][i].x = ucat[k][j][i-1].x - VelDiffIn + V_frame.x; */
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;

	}
      }
    }
  }


  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x - VelDiffIn;// + V_frame.x;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;

	}
      }
    }
  }
  

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].y = ucat[k][j+1][i].y - VelDiffIn;// + V_frame.y;
/* 	  ubcs[k][j][i].y = ucat[k][j+1][i].y + VelDiffOut + V_frame.y; */
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;

	}
      }
    }
  }
  
  
  if (user->bctype[4]==6 || user->bctype[5]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	  ubcs[k][j][i].z = ucat[k+1][j][i].z - VelDiffIn;// + V_frame.z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;

	}
      }
    }
  }
 
/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i].z*csi[k][j][i].z +
	       csi[k][j][i].y*csi[k][j][i].y +
	       csi[k][j][i].x*csi[k][j][i].x);
	nx=csi[k][j][i].x/A;
	ny=csi[k][j][i].y/A;
	nz=csi[k][j][i].z/A;
	Un=ucat[k][j][i+1].x*nx+ucat[k][j][i+1].y*ny+ucat[k][j][i+1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i+1].x-Un*nx;//-V_frame.x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i+1].z-Un*nz;
	ucont[k][j][i].x = 0.;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i-1].z*csi[k][j][i-1].z +
	       csi[k][j][i-1].y*csi[k][j][i-1].y +
	       csi[k][j][i-1].x*csi[k][j][i-1].x);
	nx=csi[k][j][i-1].x/A;
	ny=csi[k][j][i-1].y/A;
	nz=csi[k][j][i-1].z/A;
	Un=ucat[k][j][i-1].x*nx+ucat[k][j][i-1].y*ny+ucat[k][j][i-1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i-1].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i-1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i-1].z-Un*nz;
	ucont[k][j][i-1].x = 0.;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
	       eta[k][j][i].y*eta[k][j][i].y +
	       eta[k][j][i].x*eta[k][j][i].x);
	nx=eta[k][j][i].x/A;
	ny=eta[k][j][i].y/A;
	nz=eta[k][j][i].z/A;
	Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j+1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j+1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j+1][i].z-Un*nz;
	ucont[k][j][i].y = 0.;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j-1][i].z*eta[k][j-1][i].z +
	       eta[k][j-1][i].y*eta[k][j-1][i].y +
	       eta[k][j-1][i].x*eta[k][j-1][i].x);
	nx=eta[k][j-1][i].x/A;
	ny=eta[k][j-1][i].y/A;
	nz=eta[k][j-1][i].z/A;
	Un=ucat[k][j-1][i].x*nx+ucat[k][j-1][i].y*ny+ucat[k][j-1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j-1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j-1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j-1][i].z-Un*nz;
	ucont[k][j-1][i].y = 0.;
      }
    }
    }
  }
  

  if (user->bctype[4]==3) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k][j][i].z*zet[k][j][i].z +
	       zet[k][j][i].y*zet[k][j][i].y +
	       zet[k][j][i].x*zet[k][j][i].x);
	nx=zet[k][j][i].x/A;
	ny=zet[k][j][i].y/A;
	nz=zet[k][j][i].z/A;
	Un=ucat[k+1][j][i].x*nx+ucat[k+1][j][i].y*ny+ucat[k+1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k+1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k+1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k+1][j][i].z-Un*nz;
	ucont[k][j][i].z = 0.;
	}
      }
    }
  }

  if (user->bctype[5]==3) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z +
	       zet[k-1][j][i].y*zet[k-1][j][i].y +
	       zet[k-1][j][i].x*zet[k-1][j][i].x);
	nx=zet[k-1][j][i].x/A;
	ny=zet[k-1][j][i].y/A;
	nz=zet[k-1][j][i].z/A;
	Un=ucat[k-1][j][i].x*nx+ucat[k-1][j][i].y*ny+ucat[k-1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k-1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k-1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k-1][j][i].z-Un*nz;
	ucont[k-1][j][i].z = 0.;
	}
      }
    }
  }
  /* for(k=zs;k<ze;k++){ */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if((i==1)&& j==1 && ( k==0 || k==3 ) && ttemp==1) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d ubcs[k][j][i].x is %le ubcs[k][j][i].y is %le ubcs[k][j][i].z is %le \n",i,j,k,ubcs[k][j][i].x,ubcs[k][j][i].y,ubcs[k][j][i].z); */
	
/*       } */
/*     } */
/*   } */
/* ==================================================================================             */
/*     CHARACTERISTIC OUTLET BC :8 */
/* ==================================================================================             */

  if (user->bctype[5]==8) {
    if (ze == mz) {
      k = ze-2;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    // PetscGlobalSum(PETSC_COMM_WORLD,&FluxOut, &FluxOutSum);

    //ratio = FluxInSum / FluxOutSum;
    ratio = FluxIn / FluxOutSum;
    if (fabs(FluxOutSum) < 1.e-6) ratio = 1.;
    //if (fabs(FluxInSum) <1.e-6) ratio = 0.;
    if (fabs(FluxIn) <1.e-6) ratio = 0.;
    PetscPrintf(PETSC_COMM_WORLD, "Char Ratio %d %le %le %le %le %d %d\n", ti, ratio, FluxIn, FluxOutSum, FarFluxInSum,zs, ze);

    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  if (ti==0 || ti==1)
	    if (inletprofile<0)
	      ubcs[k][j][i].z = -1.;
	    else if (user->bctype[4]==6)
	      ubcs[k][j][i].z = 0.;
	    else
	      ubcs[k][j][i].z = 1.;//ubcs[0][j][i].z;//-1.;//1.;
	  
	  else
	    ucont[k-1][j][i].z = ucont[k-1][j][i].z*ratio;
	  
	  ubcs[k][j][i].z = ucont[k-1][j][i].z / zet[k-1][j][i].z;
	}
      }
    }
  }

  
/* ==================================================================================             */
/*     OUTLET BC :4 */
/* ==================================================================================             */

  /*
// temp for plate !!
   for (k=lzs; k<lze; k++) { 
     for (j=lys; j<lye; j++) { 
       for (i=lxs; i<lxe; i++) { 
 	if (k==121 && j>96 && j<145)
 	  ucont[k][j][i].z = 0.;
       } 
     } 
   } 
*/
  PetscInt period=80;
  PetscReal psys=0.273;
  if(visflg) PetscOptionsGetInt(PETSC_NULL, "-period", &period, PETSC_NULL);
  PetscInt fn;
  PetscReal FluxOut2=0.0,cellflux=0.0,cellfluxcont=0.0;
  FluxOut=0.0,lArea=0.0,AreaSum=0.0;
//  for(fn=0;fn<6;fn++){
//    if(user->bctype[fn]==4){
      switch(outletface){
        case 0:  // outlet at x=0;
         if(xs==0){
	   i=xs;
           for(k=lzs;k<lze;k++){
             for(j=lys;j<lye;j++){
               if(nvert[k][j][i+1]<0.1){
               cellflux  = (ucat[k][j][i+1].x*(csi[k][j][i].x) + 
                          ucat[k][j][i+1].y*(csi[k][j][i].y) + 
                          ucat[k][j][i+1].z*(csi[k][j][i].z)); 

               cellfluxcont =ucont[k][j][i].x;
               
               if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
               if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;

               lArea += sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			     (csi[k][j][i].y) * (csi[k][j][i].y) +
			     (csi[k][j][i].z) * (csi[k][j][i].z));
               
               cellflux=0.0;
               cellfluxcont=0.0;

                }
              }
            }
          }
         break;
        case 1:
          if(xe==mx){
	   i=xe-2;
           for(k=lzs;k<lze;k++){
             for(j=lys;j<lye;j++){
	       if(nvert[k][j][i]<0.1){
               cellflux  = (ucat[k][j][i].x*(csi[k][j][i].x) + 
                          ucat[k][j][i].y*(csi[k][j][i].y) + 
                          ucat[k][j][i].z*(csi[k][j][i].z)); 

               cellfluxcont =ucont[k][j][i].x;
               
               if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
               if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;

               lArea += sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			     (csi[k][j][i].y) * (csi[k][j][i].y) +
			     (csi[k][j][i].z) * (csi[k][j][i].z));
               
               cellflux=0.0;
               cellfluxcont=0.0;


               }
              }
            }
          }       
         break;
        case 2:
//              PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
         if(ys==0){
	   j=ys;
           for(k=lzs;k<lze;k++){
             for(i=lxs;i<lxe;i++){
               if(nvert[k][j+1][i]<0.1){
               cellflux  = (ucat[k][j+1][i].x*(eta[k][j][i].x) + 
                          ucat[k][j+1][i].y*(eta[k][j][i].y) + 
                          ucat[k][j+1][i].z*(eta[k][j][i].z)); 

               cellfluxcont =ucont[k][j][i].y;
               
               if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
               if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;

               lArea += sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			     (eta[k][j][i].y) * (eta[k][j][i].y) +
			     (eta[k][j][i].z) * (eta[k][j][i].z));
               
               cellflux=0.0;
               cellfluxcont=0.0;


               }
              }
            }
          }
 
         break;
        case 3:
//             PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
          if(ye==my){
	   j=ye-2;
           for(k=lzs;k<lze;k++){
             for(i=lxs;i<lxe;i++){
	       if(nvert[k][j][i]<0.1){
               cellflux  = (ucat[k][j][i].x*(eta[k][j][i].x) + 
                          ucat[k][j][i].y*(eta[k][j][i].y) + 
                          ucat[k][j][i].z*(eta[k][j][i].z)); 

               cellfluxcont =ucont[k][j][i].y;
               
               if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
               if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;

               lArea += sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			     (eta[k][j][i].y) * (eta[k][j][i].y) +
			     (eta[k][j][i].z) * (eta[k][j][i].z));
               
               cellflux=0.0;
               cellfluxcont=0.0;



               }
              }
            }
          }
         break;
	case 4:
//              PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
         if(zs==0){
	   k=zs;
           for(j=lys;j<lye;j++){
             for(i=lxs;i<lxe;i++){
               if(nvert[k+1][j][i]<0.1){
                cellflux = (ucat[k+1][j][i].x*(zet[k][j][i].x) +
                           ucat[k+1][j][i].y*(zet[k][j][i].y) +
                           ucat[k+1][j][i].z*(zet[k][j][i].z)); 
  
                cellfluxcont = ucont[k][j][i].z;

                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                 
               lArea += sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			     (zet[k][j][i].y) * (zet[k][j][i].y) +
			     (zet[k][j][i].z) * (zet[k][j][i].z));
              
               cellflux = 0.0;
               cellfluxcont = 0.0; 
 
              } 
             }
            }
          }
	 break;
	case 5:   
//              PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
          if(ze==mz){
	   k=ze-2;
           for(j=lys;j<lye;j++){
             for(i=lxs;i<lxe;i++){
	       if(nvert[k][j][i]<0.1){
                cellflux = (ucat[k][j][i].x*(zet[k][j][i].x) +
                           ucat[k][j][i].y*(zet[k][j][i].y) +
                           ucat[k][j][i].z*(zet[k][j][i].z)); 
  
                cellfluxcont = ucont[k][j][i].z;

                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                 
               lArea += sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			     (zet[k][j][i].y) * (zet[k][j][i].y) +
			     (zet[k][j][i].z) * (zet[k][j][i].z));
              
               cellflux = 0.0;
               cellfluxcont = 0.0; 
               
               }
              }
	   }
          } 
	 break;
      } // switch
//    } // face check 
//  } //faces loop.
  
//  PetscPrintf(PETSC_COMM_WORLD,"FormBCS check 1 \n"); 
  MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FluxOut2,&FluxOutSumcont,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
//  PetscPrintf(PETSC_COMM_WORLD,"FormBCS check 2 \n"); 
  MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  if(catcorr) user->FluxOutSum = FluxOutSum;
  user->FluxOutSum = FluxOutSumcont; // Contravariant Outflow Sum.
//  user->AreaOutSum = AreaSum;
   if(visflg) PetscPrintf(PETSC_COMM_WORLD,"FormBCS Pre-correction - Outflow: Cartesian-Flux - %le, Contravariant-Flux - %le, Area - %le \n",FluxOutSum,FluxOutSumcont,AreaSum);
   // Correction 
  FluxIn = FluxInSum + FarFluxInSum + user->FluxIntpSum;
  // ratio conditions for generalized BCS --- Vishal Kandala
 //------------------------------------------------------------------------------------------------------- 
 if(inletface==outletface) PetscPrintf(PETSC_COMM_WORLD,"Inlet and Outlet cannot be at the same surface !!!!!!!!!!!!!!!!!!!!!!!!!");
  
  if(inletface%2!=0){ // Inlet is at +x,+y or +z
     if((inletface*outletface)%2!=0) {// if outlet is also at +x,+y or +z
        if(catcorr) ratio = (-FluxIn -FluxOutSum)/AreaSum;
        else ratio = (-FluxIn -FluxOutSumcont)/AreaSum;
    } else{ // if outlet is at -x,-y or -z
       if(catcorr) ratio = (FluxIn - FluxOutSum)/AreaSum;
       else   ratio = (FluxIn - FluxOutSumcont)/AreaSum;
    }
  }// Inlet is at +x,+y or +z
  else{ // if inlet is at -x,-y or -z
   if(inletface==0){ // if inlet at -x 
     if(outletface==1 || outletface==3 || outletface==5) {// outlet at +x,+y or +z
       if(catcorr) ratio = (FluxIn - FluxOutSum)/AreaSum;
       else ratio = (FluxIn - FluxOutSumcont)/AreaSum; 
     } else{ // if outlet at -y or -z
       if(catcorr)  ratio = (-FluxIn - FluxOutSum)/AreaSum;
       else   ratio = (-FluxIn - FluxOutSumcont)/AreaSum;  
     } // if outlet at -y or -z
    }else if(inletface==2){ // if inlet at -y
     if(outletface==0 || outletface==4){ // outlet at -x or -z.
       if(catcorr) ratio = (-FluxIn - FluxOutSum)/AreaSum;
       else ratio = (-FluxIn - FluxOutSumcont)/AreaSum;    
     }else{// outlet at +x,+y or +z.
       if(catcorr)  ratio = (FluxIn - FluxOutSum)/AreaSum;
       else  ratio = (FluxIn - FluxOutSumcont)/AreaSum;  
    } // outlet at +x,+y or +z.
   }else if(inletface==4){ // if inlet at -z
     if(outletface==0 || outletface==2){  // outlet at -x or -y
       if(catcorr) ratio = (-FluxIn - FluxOutSum)/AreaSum;
       else ratio = (-FluxIn - FluxOutSumcont)/AreaSum;   
     }else{// outlet at +x,+y,+z
       if(catcorr)  ratio = (FluxIn - FluxOutSum)/AreaSum;
       else  ratio = (FluxIn - FluxOutSumcont)/AreaSum;  
     }
    } // inlet at -z
   } // inlet at -x,-y or -z 
 //-------------------------------------------------------------------------------------------------- 
  if(visflg) PetscPrintf(PETSC_COMM_WORLD,"FormBCS Correction Ratio - %le \n",ratio);
  user->FluxOutSum=0.0;
  FluxOutSum=0.0,FluxOutSumcont=0.0;
  FluxOut=0.0,FluxOut2=0;
  PetscReal FluxOutInside = 0.0; 
//  for(fn=0;fn<6;fn++){
//    if(user->bctype[fn]==4){
      switch(outletface){
        case 0:  // outlet at x=0;
         if(xs==0){
	   i=xs;
           for(k=lzs;k<lze;k++){
             for(j=lys;j<lye;j++){
               if(nvert[k][j][i+1]<0.1){
                ubcs[k][j][i].x=ucat[k][j][i+1].x;
                ubcs[k][j][i].y=ucat[k][j][i+1].y;
                ubcs[k][j][i].z=ucat[k][j][i+1].z;
 
                if(catcorr) ucont[k][j][i].x = (ubcs[k][j][i].x*(csi[k][j][i].x) +
                                      ubcs[k][j][i].y*(csi[k][j][i].y) +
                                      ubcs[k][j][i].z*(csi[k][j][i].z)) + 
                                      ratio*sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			                          (csi[k][j][i].y) * (csi[k][j][i].y) +
			                          (csi[k][j][i].z) * (csi[k][j][i].z));
                
                else ucont[k][j][i].x = ucont[k][j][i].x + ratio*sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			                          (csi[k][j][i].y) * (csi[k][j][i].y) +
			                          (csi[k][j][i].z) * (csi[k][j][i].z));
                 cellfluxcont= ucont[k][j][i].x;

                cellflux=(ucat[k][j][i+1].x*(csi[k][j][i].x) +
                           ucat[k][j][i+1].y*(csi[k][j][i].y) +
                           ucat[k][j][i+1].z*(csi[k][j][i].z));
              
                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                FluxOutInside+= ucont[k][j][i+1].x;
                cellflux=0.0;
                cellfluxcont=0.0;
               }
              }
            }
          }
         break;
        case 1:
          if(xe==mx){
	   i=xe-2;
           for(k=lzs;k<lze;k++){
             for(j=lys;j<lye;j++){
	       if(nvert[k][j][i]<0.1){                
                ubcs[k][j][i].x=ucat[k][j][i].x;
                ubcs[k][j][i].y=ucat[k][j][i].y;
                ubcs[k][j][i].z=ucat[k][j][i].z;
                
                if(catcorr) ucont[k][j][i].x = (ubcs[k][j][i].x*(csi[k][j][i].x) +
                                      ubcs[k][j][i].y*(csi[k][j][i].y) +
                                      ubcs[k][j][i].z*(csi[k][j][i].z)) + 
                                      ratio*sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			                          (csi[k][j][i].y) * (csi[k][j][i].y) +
			                          (csi[k][j][i].z) * (csi[k][j][i].z));
                
                else ucont[k][j][i].x = ucont[k][j][i].x + ratio*sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			                          (csi[k][j][i].y) * (csi[k][j][i].y) +
			                          (csi[k][j][i].z) * (csi[k][j][i].z));
                 cellfluxcont= ucont[k][j][i].x;

                cellflux=(ucat[k][j][i].x*(csi[k][j][i].x) +
                           ucat[k][j][i].y*(csi[k][j][i].y) +
                           ucat[k][j][i].z*(csi[k][j][i].z));
              
                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                FluxOutInside+= ucont[k][j][i-1].x;
                cellflux=0.0;
                cellfluxcont=0.0;
 
               }
              }
            }
          }       
         break;
        case 2:
//              PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
         if(ys==0){
	   j=ys;
           for(k=lzs;k<lze;k++){
             for(i=lxs;i<lxe;i++){
               if(nvert[k][j+1][i]<0.1){
                ubcs[k][j][i].x=ucat[k][j+1][i].x;
                ubcs[k][j][i].y=ucat[k][j+1][i].y;
                ubcs[k][j][i].z=ucat[k][j+1][i].z;
                
                if(catcorr) ucont[k][j][i].y = (ubcs[k][j][i].x*(eta[k][j][i].x) +
                                      ubcs[k][j][i].y*(eta[k][j][i].y) +
                                      ubcs[k][j][i].z*(eta[k][j][i].z)) + 
                                      ratio*sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			                          (eta[k][j][i].y) * (eta[k][j][i].y) +
			                          (eta[k][j][i].z) * (eta[k][j][i].z));
                
                else ucont[k][j][i].y = ucont[k][j][i].y + ratio*sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			                          (eta[k][j][i].y) * (eta[k][j][i].y) +
			                          (eta[k][j][i].z) * (eta[k][j][i].z));
                 cellfluxcont= ucont[k][j][i].y;

                cellflux=(ucat[k][j+1][i].x*(eta[k][j][i].x) +
                           ucat[k][j+1][i].y*(eta[k][j][i].y) +
                           ucat[k][j+1][i].z*(eta[k][j][i].z));
              
                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                FluxOutInside+= ucont[k][j+1][i].y;
                cellflux=0.0;
                cellfluxcont=0.0;
 
               }
              }
            }
          } 
         break;
        case 3:
//             PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
          if(ye==my){
	   j=ye-2;
           for(k=lzs;k<lze;k++){
             for(i=lxs;i<lxe;i++){
	       if(nvert[k][j][i]<0.1){                
                ubcs[k][j][i].x=ucat[k][j][i].x;
                ubcs[k][j][i].y=ucat[k][j][i].y;
                ubcs[k][j][i].z=ucat[k][j][i].z;
                
                if(catcorr) ucont[k][j][i].y = (ubcs[k][j][i].x*(eta[k][j][i].x) +
                                      ubcs[k][j][i].y*(eta[k][j][i].y) +
                                      ubcs[k][j][i].z*(eta[k][j][i].z)) + 
                                      ratio*sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			                          (eta[k][j][i].y) * (eta[k][j][i].y) +
			                          (eta[k][j][i].z) * (eta[k][j][i].z));
                
                else ucont[k][j][i].y = ucont[k][j][i].y + ratio*sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			                          (eta[k][j][i].y) * (eta[k][j][i].y) +
			                          (eta[k][j][i].z) * (eta[k][j][i].z));
                 cellfluxcont= ucont[k][j][i].y;

                cellflux=(ucat[k][j][i].x*(eta[k][j][i].x) +
                           ucat[k][j][i].y*(eta[k][j][i].y) +
                           ucat[k][j][i].z*(eta[k][j][i].z));
              
                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                FluxOutInside+= ucont[k][j-1][i].y;
                cellflux=0.0;
                cellfluxcont=0.0;
             
              }
              }
            }
          } 
         break;
	case 4:
//              PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
         if(zs==0){
	   k=zs;
           for(j=lys;j<lye;j++){
             for(i=lxs;i<lxe;i++){
               if(nvert[k+1][j][i]<0.1){
                ubcs[k][j][i].x=ucat[k+1][j][i].x;
                ubcs[k][j][i].y=ucat[k+1][j][i].y;
                ubcs[k][j][i].z=ucat[k+1][j][i].z;
                
                if(catcorr) ucont[k][j][i].z = (ubcs[k][j][i].x*(zet[k][j][i].x) +
                                      ubcs[k][j][i].y*(zet[k][j][i].y) +
                                      ubcs[k][j][i].z*(zet[k][j][i].z)) + 
                                      ratio*sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			                          (zet[k][j][i].y) * (zet[k][j][i].y) +
			                          (zet[k][j][i].z) * (zet[k][j][i].z));
                
                else ucont[k][j][i].z = ucont[k][j][i].z + ratio*sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			                          (zet[k][j][i].y) * (zet[k][j][i].y) +
			                          (zet[k][j][i].z) * (zet[k][j][i].z));
                
                cellfluxcont= ucont[k][j][i].z;

                cellflux=(ucat[k+1][j][i].x*(zet[k][j][i].x) +
                           ucat[k+1][j][i].y*(zet[k][j][i].y) +
                           ucat[k+1][j][i].z*(zet[k][j][i].z));
              
                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                FluxOutInside+= ucont[k+1][j][i].z;
                cellflux=0.0;
                cellfluxcont=0.0;
               }
              }
            }
          }
 
	 break;
	case 5:
//              PetscPrintf(PETSC_COMM_WORLD,"Condition yet to be implemented \n");
          if(ze==mz){
	   i=ze-2;
           for(j=lys;j<lye;j++){
             for(i=lxs;i<lxe;i++){
	       if(nvert[k][j][i]<0.1){                
                ubcs[k][j][i].x=ucat[k][j][i].x;
                ubcs[k][j][i].y=ucat[k][j][i].y;
                ubcs[k][j][i].z=ucat[k][j][i].z;
                
                if(catcorr) ucont[k][j][i].z = (ubcs[k][j][i].x*(zet[k][j][i].x) +
                                      ubcs[k][j][i].y*(zet[k][j][i].y) +
                                      ubcs[k][j][i].z*(zet[k][j][i].z)) + 
                                      ratio*sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			                          (zet[k][j][i].y) * (zet[k][j][i].y) +
			                          (zet[k][j][i].z) * (zet[k][j][i].z));

                 else ucont[k][j][i].z = ucont[k][j][i].z + ratio*sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			                          (zet[k][j][i].y) * (zet[k][j][i].y) +
			                          (zet[k][j][i].z) * (zet[k][j][i].z));


               cellfluxcont= ucont[k][j][i].z;

               cellflux=(ucat[k][j][i].x*(zet[k][j][i].x) +
                           ucat[k][j][i].y*(zet[k][j][i].y) +
                           ucat[k][j][i].z*(zet[k][j][i].z)); 
                
                if(fabs(cellflux)>FLUX_THRESHOLD) FluxOut+=cellflux;
                if(fabs(cellfluxcont)>FLUX_THRESHOLD) FluxOut2+=cellfluxcont;
                FluxOutInside+= ucont[k-1][j][i].z;
                cellflux=0.0;
                cellfluxcont=0.0;
               }
              }
            }
          }       
	 break;
      } // switch  
//    }  // face check
//  } // faces loop.
  
  MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FluxOut2,&FluxOutSumcont,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FluxOutInside,&FluxOutSumInside,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  if(catcorr) user->FluxOutSum = FluxOutSum;
  else user->FluxOutSum = FluxOutSumcont;
  user->AreaOutSum = AreaSum;
  if(visflg) PetscPrintf(PETSC_COMM_WORLD,"FormBCS Post-correction - Inside Cont Flux - %le, Outflow: Cart Flux - %le,Cont Flux - %le, Area - %le \n",FluxOutSumInside,FluxOutSum,FluxOutSumcont,AreaSum);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
/*   DALocalToLocalBegin(fda, user->lUcont, INSERT_VALUES, user->lUcont); */
/*   DALocalToLocalEnd(fda, user->lUcont, INSERT_VALUES, user->lUcont); */
    
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);

  Contra2Cart(user); // it also does global to local for Ucat

  //  DAVecGetArray(fda, user->Ucont, &ucont);

 

 
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */
  extern PetscInt wallfunction;
  if (wallfunction) {
  PetscReal ***ustar, ***aj;
  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(da, user->lUstar, &ustar);
  DMDAVecGetArray(da, user->lAj,  &aj);

  // wall function for boundary
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if( nvert[k][j][i]<1.1 && ( (user->bctype[0]==-1 && i==1) ||
				    (user->bctype[1]==-1 && i==mx-2) ) ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double sb, sc;
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(i==1) {
	    sc = 2*sb + 0.5/aj[k][j][i+1]/area;
	    Uc = ucat[k][j][i+1];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j][i-1]/area;
	    Uc = ucat[k][j][i-1];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
	  
	  //if(i==1) printf("%f %f, %f %f %f, %e %e %e, %e\n", sc, sb, Ua.x, Ua.y, Ua.z, Uc.x, Uc.y, Uc.z, ucat[k][j][i+2].z);
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //PetscPrintf(PETSC_COMM_WORLD, "ustar %d %d %d %le !\n",i,j,k,ustar[k][j][i]);
	  //	  ucat[k][j][i] = lucat[k][j][i];
	}
	
	if( nvert[k][j][i]<1.1 && ( (user->bctype[2]==-1 && j==1) ||
				    (user->bctype[3]==-1 && j==my-2) ) ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc;
	  double ni[3], nj[3], nk[3];
	   Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(j==1) {
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j-1][i]/area;
	    Uc = ucat[k][j-1][i];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	  
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];
	}
	
      }
  DMDAVecRestoreArray(da, user->lAj,  &aj);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMDAVecGetArray(da, user->Nvert, &nvert);
  }

/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Ucat, &ucat);
 
/* ==================================================================================             */
 
  // boundary conditions on ghost nodes
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }
 
  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
 /* ==================================================================================             */
  /*   Periodic BC *///Mohsen
/* ==================================================================================             */
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  DMDAVecGetArray(da, user->lP, &p);
  DMDAVecGetArray(fda, user->lUcat, &ucat);
   
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  if(j>0 && k>0 && j<user->JM && k<user->KM){
	    ucat[k][j][i].x=ucat[k][j][i-2].x;
	    ucat[k][j][i].y=ucat[k][j][i-2].y;
	    ucat[k][j][i].z=ucat[k][j][i-2].z;
	    p[k][j][i]=p[k][j][i-2];
	  }
	}
      }
    }
  }
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && k>0 && i<user->IM && k<user->KM){
	    ucat[k][j][i].x=ucat[k][j-2][i].x;
	    ucat[k][j][i].y=ucat[k][j-2][i].y;
	    ucat[k][j][i].z=ucat[k][j-2][i].z;
	    p[k][j][i]=p[k][j-2][i];
	  }
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && j>0 && i<user->IM && j<user->JM){
	    ucat[k][j][i].x=ucat[k-2][j][i].x;
	    ucat[k][j][i].y=ucat[k-2][j][i].y;
	    ucat[k][j][i].z=ucat[k-2][j][i].z;
	    p[k][j][i]=p[k-2][j][i];
	  }
	}
      }
    }
  }
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xe==mx){
      i=mx-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  if(j>0 && k>0 && j<user->JM && k<user->KM){
	    ucat[k][j][i].x=ucat[k][j][i+2].x;
	    ucat[k][j][i].y=ucat[k][j][i+2].y;
	    ucat[k][j][i].z=ucat[k][j][i+2].z;
	    p[k][j][i]=p[k][j][i+2];
	  }
	}
      }
    }
  }
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ye==my){
      j=my-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && k>0 && i<user->IM && k<user->KM){
	    ucat[k][j][i].x=ucat[k][j+2][i].x;
	    ucat[k][j][i].y=ucat[k][j+2][i].y;
	    ucat[k][j][i].z=ucat[k][j+2][i].z;
	    p[k][j][i]=p[k][j+2][i];
	  }
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (ze==mz){
      k=mz-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  if(i>0 && j>0 && i<user->IM && j<user->JM){
	    ucat[k][j][i].x=ucat[k+2][j][i].x;
	    ucat[k][j][i].y=ucat[k+2][j][i].y;
	    ucat[k][j][i].z=ucat[k+2][j][i].z;
	    p[k][j][i]=p[k+2][j][i];
	  }
	}
      }
    }
  }
 
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(da, user->lP, &p);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

 // 0 velocity on the corner point

  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->P, &p);

  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k+1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k][j+1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j+1][i]+p[k][j][i-1]);
      }
    }
  }
 
  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i-1]);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da, user->P, &p);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  //Mohsen Nov 2012
  //Velocity and Presurre at corners for Periodic BC's
  //i-direction

  DMDAVecGetArray(fda, user->lUcat,  &ucat);
  DMDAVecGetArray(da, user->lP, &p);

  if (user->bctype[0]==7){
    if (xs==0){
      i=xs;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i-2].x;
	  ucat[k][j][i].y=ucat[k][j][i-2].y;
	  ucat[k][j][i].z=ucat[k][j][i-2].z;
	  p[k][j][i]=p[k][j][i-2];
	}
      }
    }
  }
  if (user->bctype[1]==7){
    if (xe==mx){
      i=xe-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i+2].x;
	  ucat[k][j][i].y=ucat[k][j][i+2].y;
	  ucat[k][j][i].z=ucat[k][j][i+2].z;
	  p[k][j][i]=p[k][j][i+2];
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(da, user->lP, &p);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);
  
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  //j-direction
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
  DMDAVecGetArray(da, user->lP, &p);

  if (user->bctype[2]==7){
    if (ys==0){
      j=ys;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j-2][i].x;
	  ucat[k][j][i].y=ucat[k][j-2][i].y;
	  ucat[k][j][i].z=ucat[k][j-2][i].z;
	  p[k][j][i]=p[k][j-2][i];
	}
      }
    }
  }
  
  if (user->bctype[3]==7){
    if (ye==my){
      j=my-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j+2][i].x;
	  ucat[k][j][i].y=ucat[k][j+2][i].y;
	  ucat[k][j][i].z=ucat[k][j+2][i].z;
	  p[k][j][i]=p[k][j+2][i];
	}
      }
    }
  }

  DMDAVecRestoreArray(da, user->lP, &p);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  
  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  //k-direction
  
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
  DMDAVecGetArray(da, user->lP, &p);
  
  if (user->bctype[4]==7){
    if (zs==0){
      k=zs;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k-2][j][i].x;
	  ucat[k][j][i].y=ucat[k-2][j][i].y;
	  ucat[k][j][i].z=ucat[k-2][j][i].z;
	  p[k][j][i]=p[k-2][j][i];
	}
      }
    }
  }
  if (user->bctype[5]==7){
    if (ze==mz){
      k=mz-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k+2][j][i].x;
	  ucat[k][j][i].y=ucat[k+2][j][i].y;
	  ucat[k][j][i].z=ucat[k+2][j][i].z;
	  p[k][j][i]=p[k+2][j][i];
	}
      }
    }
  }
   
  DMDAVecRestoreArray(da, user->lP, &p);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);

  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat);
  DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat);

  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P);
  DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
 
 /* ==================================================================================             */
/*   Analytical Vortex BC */
/* ==================================================================================             */
 
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Cent, &cent); 
  DMDAVecGetArray(fda, user->Centx, &centx);
  DMDAVecGetArray(fda, user->Centy, &centy);
  DMDAVecGetArray(fda, user->Centz, &centz);

  if (user->bctype[0]==9) {
    if (xs==0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i+1].x)*sin(cent[k][j][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i+1].x)*cos(cent[k][j][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	 
	  ucont[k][j][i].x =-(cos(centx[k][j][i].x)*sin(centx[k][j][i].y)*csi[k][j][i].x)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i+1].x)*sin(cent[k][j+1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i+1].x)*cos(cent[k][j+1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i+1].x)*sin(cent[k][j-1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i+1].x)*cos(cent[k][j-1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }
  if (user->bctype[1]==9) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i-1].x)*sin(cent[k][j][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i-1].x)*cos(cent[k][j][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j][i-1].x =(-cos(centx[k][j][i-1].x)*sin(centx[k][j][i-1].y)*csi[k][j][i-1].x)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i-1].x)*sin(cent[k][j+1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i-1].x)*cos(cent[k][j+1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i-1].x)*sin(cent[k][j-1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i-1].x)*cos(cent[k][j-1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }

  if (user->bctype[2]==9) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i].x)*sin(cent[k][j+1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=sin(cent[k][j+1][i].x)*cos(cent[k][j+1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;

	  ucont[k][j][i].y=(sin(centy[k][j][i].x)*cos(centy[k][j][i].y)*eta[k][j][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
    }
  }
 
 
  if (user->bctype[3]==9) {
    if (ye==my) {
      j= ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i].x)*sin(cent[k][j-1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=sin(cent[k][j-1][i].x)*cos(cent[k][j-1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j-1][i].y=(sin(centy[k][j-1][i].x)*cos(centy[k][j-1][i].y)*eta[k][j-1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
    }
  }
  if (user->bctype[4]==9) {
    if (zs==0) {
      k= zs;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k+1][j][i].x;
	  ucat[k][j][i].y=ucat[k+1][j][i].y;
	  ucat[k][j][i].z=ucat[k+1][j][i].z;

	  ucont[k][j][i].z=0.0;
	}
      }
    }
  }
  if (user->bctype[5]==9) {
    if (ze==mz) {
      k= ze-1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k-1][j][i].x;
	  ucat[k][j][i].y=ucat[k-1][j][i].y;
	  ucat[k][j][i].z=ucat[k-1][j][i].z;

	  ucont[k-1][j][i].z=0.0;
	}
      }
    }
  }
/*  for(k=zs;k<ze;k++){ */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if((i==1) && j==1 && (k==0 || k==1 ||k==2 || k==3) && ttemp==1) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d ucat[k][j][i].x is %le ucat[k][j][i].y is %le ucat[k][j][i].z is %le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z);   */
	
/*       } */
/*     } */
/*   } */
  /* for(k=zs;k<ze;k++){ */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if((i==0 || i==1) && j==1 && (k==0 || k==1 || k==2) && ttemp==2)  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d ucat[k][j][i].x  is %le ucat[k][j][i].y is %le ucat[k][j][i].z is %le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z); */
/*       } */
/*     } */
/*   } */

  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->Centx, &centx);
  DMDAVecRestoreArray(fda, user->Centy, &centy);
  DMDAVecRestoreArray(fda, user->Centz, &centz);

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(fda, user->Ucont,  &ucont);

  } // ttemp

 
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); 

  return(0);
}


PetscErrorCode ForceSolidBoundary(UserCtx *user)
{
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k,fn,flag;

  
  PetscReal	***nvert; 

  Cmpnts	***ucont;
 

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx ) lxe = xe-1;
  if (ye==my ) lye = ye-1;
  if (ze==mz ) lze = ze-1;


  DMDAVecGetArray(da, user->Nvert, &nvert);

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  flag=0;
 for (fn=0; fn<6; fn++) {
  if (user->bctype[fn] == 1) {
    flag=1;
    switch(fn){
     
    case 0:  // face 0
      if (xs==0) {

       	i = 0;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    
	    if (nvert[k][j][i+1]<0.1) {
	      ucont[k][j][i].x =0.0; 
	    }//end if 

	  }// location for loop
	}// location for loop

      }
      break;
    
    case 1:  // face 1
      if (xe==mx) {
	
	i = mx-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    
	    if (nvert[k][j][i]<0.1) {
	      ucont[k][j][i].x =0.0; ;  
	    }//end if 

	  }// location for loop
	}// location for loop

      }
      break;
     
    case 2: // face 2
      if (ys==0) {

	j = 0;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {

	    if (nvert[k][j+1][i]<0.1) {
	      ucont[k][j][i].y =0.0;  
	    }//end if 

	  }// location for loop
	}// location for loop

      }
      break;
     
    case 3: // face 3
      if (ye==my) {

	j = my-2;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {

	    if (nvert[k][j][i]<0.1) {
	      ucont[k][j][i].y =0.0; 
	    }//end if 

	  }// location for loop
	}// location for loop

      }
      break;
     
    case 4: // face 4
      if (zs==0) {
	
	k = 0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {

	    if (nvert[k+1][j][i]<0.1) {
	      ucont[k][j][i].z =0.0; 
	    }//end if 

	  }// location for loop
	}// location for loop
  	
      }
      break;
     
    case 5: // face 5
      if (ze==mz) {
	
	k = mz-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {

	    if (nvert[k][j][i]<0.1) {
	       ucont[k][j][i].z =0.0;
	    }//end if 

	  }// location for loop
	}// location for loop
      }
      break;
    }//end switch

  }// end outlet check
  }// end face counter 

  DMDAVecRestoreArray(fda, user->Ucont,  &ucont);
  if (flag)
    Contra2Cart(user); 
  DMDAVecRestoreArray(da, user->Nvert, &nvert);

 return(0);
}

PetscErrorCode fluxin(UserCtx *user)
{
  PetscInt  iRotate;

  PetscInt ts_p_cycle;
  PetscInt opening, closing;
  PetscInt open_steps, close_steps;

/*   ts_p_cycle = 500; */
/*   opening = 0; */
/*   open_steps = 50; */

/*   closing = 225; */
/*   close_steps = 50; */

  //ts_p_cycle = 1290;
  ts_p_cycle = 2500;//10000;//5000;

  opening = 10;
  open_steps = 100;

  closing = 580;
  close_steps = 80;

  PetscReal t_rel;

  iRotate = ti - ((ti / ts_p_cycle) * ts_p_cycle);

  // if (angle>.0 && iRotate>1058) iRotate-=angle;

  t_rel = iRotate * (1. / ts_p_cycle) * 860 + 6.8  ; //+7.15;
/*   t_rel = (iRotate-940) * (1. / ts_p_cycle) * 860/2 + */
/*                    940. * (1. / 2500.     ) * 860 + 6.8;     */

  PetscInt i;
  PetscBool interpolated = PETSC_FALSE;
  //PetscPrintf(PETSC_COMM_WORLD, "Inflow00 Rate %d %e %e %i\n",ti, Flux_in, t_rel, user->number_flowwave);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    for (i=0; i<user->number_flowwave-1; i++) {
      /*       PetscPrintf(PETSC_COMM_WORLD, "Inflow Rate %e %e\n", Flux_in, user->inflow[i].t); */
      if (t_rel >= user->inflow[i].t && t_rel <= user->inflow[i+1].t) {
	Flux_in = user->inflow[i].f + (user->inflow[i+1].f - user->inflow[i].f) /
	  (user->inflow[i+1].t - user->inflow[i].t) *
	  (t_rel - user->inflow[i].t);
	/* 	PetscPrintf(PETSC_COMM_SELF, "Inflow Rate %i %e %e %e %e %e %e\n", i, Flux_in, t_rel, user->inflow[i].f, user->inflow[i].t, user->inflow[i+1].f, user->inflow[i+1].t); */
	interpolated = PETSC_TRUE;
      }
      if (interpolated) break;
    }  
    //if (t_rel > 350) Flux_in = 0.;
    
    //Flux_in=Flux_in;
    
    MPI_Bcast(&Flux_in, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  }
  else {
    MPI_Bcast(&Flux_in, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  }
  
  if (user->bctype[4]>1) {
    if (Flux_in<0.) {
      user->bctype[5]= 5;
      user->bctype[4]= 4;
      PetscPrintf(PETSC_COMM_WORLD, "IINNFLOW change inlet!!!!!!!%e\n", Flux_in);
    } else {
      user->bctype[5]= 4;
      user->bctype[4]= 5;
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "Angle %d %le %le flux-in %le intp%d\n",ti, t_rel, angle, Flux_in, interpolated);

  return 0;
}

PetscErrorCode OutflowVelocity(UserCtx *user, Vec Ucont)
{
  DM            fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	lFluxOut = 0., ratio;
  Cmpnts        ***ucont, ***zet, ***ucat, ***ubcs;
  
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  if (user->bctype[5] == 4) {

    //    OutflowFlux(user);

    Contra2Cart(user);
    DMDAVecGetArray(fda, Ucont, &ucont);
    DMDAVecGetArray(fda, user->lZet, &zet);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
    DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
    /* Inflow flux at previous time step is 0, under this condition, it's assumed
       the flux difference at two time steps is uniformly distributed 
       to all outflow boundary nodes.*/
    if (ti==1) {//Flux_in_old < 1.e-6) {
      PetscReal lArea = 0., AreaSum;
      
      if (ze == mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lArea += zet[k][j][i].z;
	  }
	}
      }
      MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea, &AreaSum);

      PetscReal vd;
      vd = (user->FluxInSum) / AreaSum;
      PetscPrintf(PETSC_COMM_SELF, "FluxOOOO %e %e %e\n", vd, Flux_in, Flux_in);
      
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].z = vd;
	    ucont[k-1][j][i].z = vd * zet[k-1][j][i].z;
	  }
	}
      }
    }
    /* Scale the outflow flux to ensure global flux conservation */
    else {
      lFluxOut = 0.;
      if (ze==mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lFluxOut += (ucat[k][j][i].x * (zet[k][j][i].x + zet[k-1][j][i].x) +
			 ucat[k][j][i].y * (zet[k][j][i].y + zet[k-1][j][i].y) +
			 ucat[k][j][i].z * (zet[k][j][i].z + zet[k-1][j][i].z))
	      * 0.5;
	  }
	}
      }

      MPI_Allreduce(&lFluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      //  PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut, &FluxOutSum);
      PetscBarrier(PETSC_NULL);
      ratio = user->FluxInSum / FluxOutSum;

      PetscPrintf(PETSC_COMM_WORLD, "Ratio %e %e\n", ratio, FluxOutSum);
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z * ratio;
	    ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  }
	}
      }
      
    }
    DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(fda, Ucont, &ucont);
    DMDAVecRestoreArray(fda, user->lZet, &zet);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

/*     DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*     DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */

/*     Contra2Cart(user, user->lUcont, user->Ucat); */
  }

  return 0;
}





PetscErrorCode SetInitialGuessToOne(UserCtx *user)
{
  DM da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  Vec           Coor;
  PetscReal	***nvert; 

  Cmpnts ***ucont, ***zet,***coor, ***csi, ***eta, ***cent,***centx,***centy,***centz;
  
  PetscReal xc,yc,zc,r,uin,***p;
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

 
  //DMDAGetGhostedCoordinates(da, &Coor);
  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->lZet,  &zet);
  DMDAVecGetArray(fda, user->lCsi,  &csi);  
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  
  extern PetscInt InitialGuessOne;
  if(visflg)  PetscPrintf(PETSC_COMM_WORLD,"InitialGuess: %d \n",InitialGuessOne); 
  for (k=zs ; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	
        if (InitialGuessOne==1) {
	    ucont[k][j][i].x  = 0.;
	    ucont[k][j][i].y  = 0.;
            ucont[k][j][i].z  = 0.;
           
	    //if (nvert[k][j][i]+nvert[k+1][j][i]<0.1) 
	  //  ucont[k][j][i].z  	    
	  //    = uin* sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
	//		   (zet[k][j][i].y) * (zet[k][j][i].y) +
	//		   (zet[k][j][i].z) * (zet[k][j][i].z));
	  }
      }
    }
  }
  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);  
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  // VecDestroy(&Coor);

  DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);	
  DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);

  Contra2Cart(user); 
  return(0);

}

PetscErrorCode MomentumJet(UserCtx *user, PetscInt Kplane)
{
  DM            fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***ucont, ***ucat;   
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->lUcat , &ucat);
 
  k=Kplane;
  PetscReal lMomentum =0., Momentum=0.;
  if (k>zs && k<ze) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	 
	lMomentum += fabs(ucont[k][j][i].z) * 0.5 *
	  (ucat[k][j][i].z + ucat[k+1][j][i].z);
      }
    }
  }

  MPI_Allreduce(&lMomentum,&Momentum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  //  PetscGlobalSum(PETSC_COMM_WORLD,&lMomentum, &Momentum);
  
  PetscReal lMomentum_in =0., Momentum_in=0.;
  if (zs==0) {
    k=0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	 
	lMomentum_in += fabs(ucont[k][j][i].z) * 0.5 *
	  (ucat[k][j][i].z + ucat[k+1][j][i].z);
      }
    }
  }

  MPI_Allreduce(&lMomentum_in,&Momentum_in,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  //  PetscGlobalSum(PETSC_COMM_WORLD,&lMomentum_in, &Momentum_in);

  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->lUcat , &ucat);

  PetscPrintf(PETSC_COMM_WORLD, "Jet Momentum %d %d %le %le %le fluxin %le\n", ti, k, Momentum, Momentum_in, Flux_in*fabs(user->FluxInSum), Flux_in);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {	
    char filen[80];
    sprintf(filen, "JetMomentum_%2.2d.dat", user->_this);
    FILE *f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d\t%.7e %le %le %le\n", ti, Momentum, Momentum_in, Flux_in*fabs(user->FluxInSum), Flux_in);
    fclose(f);
  }

  return(0);
}









/* PetscErrorCode Interpolation_matrix(UserMG *usermg) // for two block */
/* { */
  
/*   PetscErrorCode       ierr; */
/*   PetscInt             i,j,k; */
/*   char                 filen[80]; */
/*   PetscInt             rank,size,bi,level; */
/*   //DM                   packer; */
/*   DMDALocalInfo        info1,info2,info3; */
/*   Vec                  Ub[block_number]; */
/*   PetscViewer          viewer; */
  
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   MPI_Comm_size(PETSC_COMM_WORLD, &size); */

/*   level = usermg->mglevels-1; */
/*   UserCtx *user = usermg->mgctx[level].user; */

  
/*   ierr = DMCompositeCreate(PETSC_COMM_WORLD,&int_packer);CHKERRQ(ierr); */
  
/*   for(bi=0;bi<block_number;bi++){ */
/*     ierr = DMCompositeAddDM(int_packer,user[bi].fda);CHKERRQ(ierr); */
/*   } */
  
/*   DMSetUp(int_packer); */
  
/*   ierr = DMCreateGlobalVector(int_packer,&U_int);CHKERRQ(ierr); */
/*   ierr = DMCreateGlobalVector(int_packer,&Umult_int);CHKERRQ(ierr); */
  
/*   //////////////////// */
/*   DMDAGetLocalInfo(user[0].da,&info1); */
/*   DMDAGetLocalInfo(user[1].da,&info2); */
/*   DMDAGetLocalInfo(user[2].da,&info3); */
  
/* /\*   PetscInt	xs1 = info1.xs, xe1 = info1.xs + info1.xm; *\/ */
/* /\*   PetscInt  	ys1 = info1.ys, ye1 = info1.ys + info1.ym; *\/ */
/* /\*   PetscInt	zs1 = info1.zs, ze1 = info1.zs + info1.zm; *\/ */
/* /\*   PetscInt	mx1 = info1.mx, my1 = info1.my, mz1 = info1.mz; *\/ */
  
/* /\*   PetscInt	xs2 = info2.xs, xe2 = info2.xs + info2.xm; *\/ */
/* /\*   PetscInt  	ys2 = info2.ys, ye2 = info2.ys + info2.ym; *\/ */
/* /\*   PetscInt	zs2 = info2.zs, ze2 = info2.zs + info2.zm; *\/ */
/* /\*   PetscInt	mx2 = info2.mx, my2 = info2.my, mz2 = info2.mz; *\/ */

/* /\*   PetscInt	xs3 = info3.xs, xe3 = info3.xs + info3.xm; *\/ */
/* /\*   PetscInt  	ys3 = info3.ys, ye3 = info3.ys + info3.ym; *\/ */
/* /\*   PetscInt	zs3 = info3.zs, ze3 = info3.zs + info3.zm; *\/ */
/* /\*   PetscInt	mx3 = info3.mx, my3 = info3.my, mz3 = info3.mz; *\/ */

/* /\*   PetscPrintf(PETSC_COMM_SELF, "rank:%d---xs1=%d ys1=%d zs1=%d xe1=%d ye1=%d ze1=%d xm1=%d ym1=%d zm1=%d !\n",rank,xs1,ys1,zs1,xe1,ye1,ze1,info1.xm,info1.ym,info1.zm); *\/ */
/* /\*   PetscPrintf(PETSC_COMM_SELF, "rank:%d---xs2=%d ys2=%d zs2=%d xe2=%d ye2=%d ze2=%d xm2=%d ym2=%d zm2=%d \n",rank,xs2,ys2,zs2,xe2,ye2,ze2,info2.xm,info2.ym,info2.zm); *\/ */
  
/*   /////////////////////////////////////////////////////////////////// allocate local vectors */
/*   //PetscPrintf(PETSC_COMM_SELF, "%d %d %d !\n",info1.xm,info1.ym,info1.zm); */
/*   //PetscPrintf(PETSC_COMM_SELF, "%d %d %d !\n",info2.xm,info2.ym,info2.zm); */


/*   DMCompositeGetAccess(int_packer,U_int,&Ub[0],&Ub[1],&Ub[2]); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(user[bi].Ucat, Ub[bi]); */
/*   } */
  
/*   DMCompositeRestoreAccess(int_packer,U_int,&Ub[0],&Ub[1],&Ub[2]); */


/*   PetscInt N11=0,N1_local=0,NN=0,NN_local=0; */
/*   //Mat Int_matrix; */
  
/*   VecGetSize(U_int,&N11); */
/*   NN+=N11; */
/*   VecGetLocalSize(U_int,&N1_local); */
/*   NN_local+=N1_local; */

/*   // PetscPrintf(PETSC_COMM_SELF, "Size %d !\n",NN); */
/*   //PetscPrintf(PETSC_COMM_SELF, "rank:%d ,Size local %d !\n",rank, N1_local); */

/* /////////////////////////////////////////////////////////////////// create and allocate a matrix */
/*   ISLocalToGlobalMapping *ltogs,ltog; */
/*   DMCompositeGetISLocalToGlobalMappings(int_packer,&ltogs); */
/*   ISLocalToGlobalMappingConcatenate(PETSC_COMM_SELF,block_number,ltogs,&ltog); */
  
/*   MatCreateAIJ(PETSC_COMM_WORLD, NN_local, NN_local, NN, NN,8, PETSC_NULL,8, PETSC_NULL,&Int_matrix); */
/*   MatSetOption(Int_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); */

/*   MatSetLocalToGlobalMapping(Int_matrix,ltog,ltog); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     ISLocalToGlobalMappingDestroy(&ltogs[bi]); */
/*   } */
/*   ISLocalToGlobalMappingDestroy(&ltog); */
/*   PetscFree(ltogs); */

  
/*   /////////////////////////////////////////////////////////////// obtaining global indeces(starting index for each cpu) */

/*   PetscInt *dg_start,*dgl_start,*dl_start_blk0,*dl_start_blk1,*dl_start_blk2; */
  
/*   PetscInt dl_start_blk00,dl_start_blk11,dl_start_blk22; */
/*   PetscMalloc(size*sizeof(PetscInt), &(dgl_start)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(dg_start)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(dl_start_blk0)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(dl_start_blk1)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(dl_start_blk2)); */
    
/*   PetscInt cpu=0; */
  
/*   dl_start_blk00=3*(info1.xm*info1.ym*info1.zm); */
/*   dl_start_blk11=3*(info2.xm*info2.ym*info2.zm); */
/*   dl_start_blk22=3*(info3.xm*info3.ym*info3.zm); */
  
/*   for(cpu=0;cpu<size;cpu++){ */
/*     dg_start[cpu]=0; */
/*     dgl_start[cpu]=0; */
/*     dl_start_blk0[cpu]=0; */
/*     dl_start_blk1[cpu]=0; */
/*     dl_start_blk2[cpu]=0; */
/*   } */
  

/*   MPI_Gather(&dl_start_blk00, 1, MPIU_INT,dl_start_blk0, 1, MPIU_INT, 0, MPI_COMM_WORLD); */
/*   MPI_Gather(&dl_start_blk11, 1, MPIU_INT,dl_start_blk1, 1, MPIU_INT, 0, MPI_COMM_WORLD); */
/*   MPI_Gather(&dl_start_blk22, 1, MPIU_INT,dl_start_blk2, 1, MPIU_INT, 0, MPI_COMM_WORLD); */
/*   PetscBarrier(NULL); */
  
/*   MPI_Bcast(dl_start_blk0,size , MPIU_INT, 0, MPI_COMM_WORLD); */
/*   MPI_Bcast(dl_start_blk1,size , MPIU_INT, 0, MPI_COMM_WORLD); */
/*   MPI_Bcast(dl_start_blk2,size , MPIU_INT, 0, MPI_COMM_WORLD); */
/*   PetscBarrier(NULL); */

/*   cpu=0; */
/*   while(cpu<(rank)){ */
/*     dgl_start[rank]+=(dl_start_blk0[cpu]+dl_start_blk1[cpu]+dl_start_blk2[cpu]); */
/*     cpu++; */
/*    } */

/*   for (i=0;i<size;i++){ */
/*     MPI_Allreduce(&dgl_start[i],&dg_start[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*   } */

/*   PetscBarrier(NULL); */
/*   //PetscPrintf(PETSC_COMM_SELF, "rank:%d, dg: %d %d %d %d!\n",rank,dg_start[0],dg_start[1],dg_start[2],dg_start[3]); */


/*   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////// sending distribution info to all cpus */
/*   PetscInt *cpu_lxs11,*cpu_lys11,*cpu_lzs11,*cpu_lxm11,*cpu_lym11,*cpu_lzm11; */
/*   PetscInt *cpu_xs11,*cpu_ys11,*cpu_zs11,*cpu_xm11,*cpu_ym11,*cpu_zm11; */
/*   PetscInt *cpu_lxs22,*cpu_lys22,*cpu_lzs22,*cpu_lxm22,*cpu_lym22,*cpu_lzm22; */
/*   PetscInt *cpu_xs22,*cpu_ys22,*cpu_zs22,*cpu_xm22,*cpu_ym22,*cpu_zm22; */
/*   PetscInt *cpu_lxs33,*cpu_lys33,*cpu_lzs33,*cpu_lxm33,*cpu_lym33,*cpu_lzm33; */
/*   PetscInt *cpu_xs33,*cpu_ys33,*cpu_zs33,*cpu_xm33,*cpu_ym33,*cpu_zm33; */
/*    /////////////////////////////blk1 */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_lxs11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lys11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzs11)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_lxm11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lym11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzm11)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_xs11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ys11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zs11)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_xm11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ym11)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zm11)); */
/*   /////////////////////////////blk2 */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_lxs22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lys22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzs22)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_lxm22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lym22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzm22)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_xs22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ys22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zs22)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_xm22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ym22)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zm22)); */
/*   /////////////////////////////blk3 */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_lxs33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lys33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzs33)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_lxm33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lym33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_lzm33)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_xs33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ys33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zs33)); */
/*   PetscMalloc(size*sizeof(PetscInt), &(cpu_xm33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_ym33)); PetscMalloc(size*sizeof(PetscInt), &(cpu_zm33)); */
  
/*   for(cpu=0;cpu<size;cpu++){ */
/*     cpu_lxs11[cpu]=0;cpu_lys11[cpu]=0;cpu_lzs11[cpu]=0;cpu_xs11[cpu]=0;cpu_ys11[cpu]=0;cpu_zs11[cpu]=0; */
/*     cpu_lxm11[cpu]=0;cpu_lym11[cpu]=0;cpu_lzm11[cpu]=0;cpu_xm11[cpu]=0;cpu_ym11[cpu]=0;cpu_zm11[cpu]=0; */
/*     cpu_lxs22[cpu]=0;cpu_lys22[cpu]=0;cpu_lzs22[cpu]=0;cpu_xs22[cpu]=0;cpu_ys22[cpu]=0;cpu_zs22[cpu]=0; */
/*     cpu_lxm22[cpu]=0;cpu_lym22[cpu]=0;cpu_lzm22[cpu]=0;cpu_xm22[cpu]=0;cpu_ym22[cpu]=0;cpu_zm22[cpu]=0; */
/*     cpu_lxs33[cpu]=0;cpu_lys33[cpu]=0;cpu_lzs33[cpu]=0;cpu_xs33[cpu]=0;cpu_ys33[cpu]=0;cpu_zs33[cpu]=0; */
/*     cpu_lxm33[cpu]=0;cpu_lym33[cpu]=0;cpu_lzm33[cpu]=0;cpu_xm33[cpu]=0;cpu_ym33[cpu]=0;cpu_zm33[cpu]=0; */
/*   } */

/*   ///------------------------/// */
/*   cpu_lxs11[rank]=info1.xs;cpu_lys11[rank]=info1.ys;cpu_lzs11[rank]=info1.zs; */
/*   cpu_lxm11[rank]=info1.xm;cpu_lym11[rank]=info1.ym;cpu_lzm11[rank]=info1.zm; */
  
/*   cpu_lxs22[rank]=info2.xs;cpu_lys22[rank]=info2.ys;cpu_lzs22[rank]=info2.zs; */
/*   cpu_lxm22[rank]=info2.xm;cpu_lym22[rank]=info2.ym;cpu_lzm22[rank]=info2.zm; */

/*   cpu_lxs33[rank]=info3.xs;cpu_lys33[rank]=info3.ys;cpu_lzs33[rank]=info3.zs; */
/*   cpu_lxm33[rank]=info3.xm;cpu_lym33[rank]=info3.ym;cpu_lzm33[rank]=info3.zm; */

/*   PetscBarrier(NULL); */
  
/*   for (i=0;i<size;i++){ */
/*     MPI_Allreduce(&cpu_lxs11[i],&cpu_xs11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lys11[i],&cpu_ys11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lzs11[i],&cpu_zs11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lxm11[i],&cpu_xm11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lym11[i],&cpu_ym11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lzm11[i],&cpu_zm11[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */

/*     MPI_Allreduce(&cpu_lxs22[i],&cpu_xs22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lys22[i],&cpu_ys22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lzs22[i],&cpu_zs22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lxm22[i],&cpu_xm22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lym22[i],&cpu_ym22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lzm22[i],&cpu_zm22[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */

/*     MPI_Allreduce(&cpu_lxs33[i],&cpu_xs33[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lys33[i],&cpu_ys33[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lzs33[i],&cpu_zs33[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lxm33[i],&cpu_xm33[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lym33[i],&cpu_ym33[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&cpu_lzm33[i],&cpu_zm33[i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); */
/*   } */
  
/*   PetscBarrier(NULL); */
/*   //PetscPrintf(PETSC_COMM_SELF, "rank:%d, lx: %d %d!\n",rank,cpu_xm11[3],cpu_xm22[3]); */
  
/*   //////----------------------------------------------------/// */
/*   PetscInt cpu_xs[block_number][size];PetscInt cpu_xm[block_number][size]; */
/*   PetscInt cpu_ys[block_number][size];PetscInt cpu_ym[block_number][size]; */
/*   PetscInt cpu_zs[block_number][size];PetscInt cpu_zm[block_number][size]; */

/*   for (bi=0; bi<block_number; bi++) { */
    
/*     for (i=0;i<size;i++){ */
/*       if(bi==0){ */
/* 	cpu_xs[bi][i]=cpu_xs11[i]; */
/* 	cpu_xm[bi][i]=cpu_xm11[i]; */
/* 	cpu_ys[bi][i]=cpu_ys11[i]; */
/* 	cpu_ym[bi][i]=cpu_ym11[i]; */
/* 	cpu_zs[bi][i]=cpu_zs11[i]; */
/* 	cpu_zm[bi][i]=cpu_zm11[i]; */
/*       }else if(bi==1){ */
/* 	cpu_xs[bi][i]=cpu_xs22[i]; */
/* 	cpu_xm[bi][i]=cpu_xm22[i]; */
/* 	cpu_ys[bi][i]=cpu_ys22[i]; */
/* 	cpu_ym[bi][i]=cpu_ym22[i]; */
/* 	cpu_zs[bi][i]=cpu_zs22[i]; */
/* 	cpu_zm[bi][i]=cpu_zm22[i]; */
/*       }else if(bi==2){ */
/* 	cpu_xs[bi][i]=cpu_xs33[i]; */
/* 	cpu_xm[bi][i]=cpu_xm33[i]; */
/* 	cpu_ys[bi][i]=cpu_ys33[i]; */
/* 	cpu_ym[bi][i]=cpu_ym33[i]; */
/* 	cpu_zs[bi][i]=cpu_zs33[i]; */
/* 	cpu_zm[bi][i]=cpu_zm33[i];	 */
/*       } */
/*       //PetscPrintf(PETSC_COMM_SELF, "rank:%d  bi: %d, lx: %d %d %d %d %d %d !\n",i,bi,cpu_xs[bi][i],cpu_xm[bi][i],cpu_ys[bi][i],cpu_ym[bi][i],cpu_zs[bi][i],cpu_zm[bi][i]); */
/*     } */
/*   } */
/*   PetscBarrier(NULL); */

/*   /////////////////////////////////////----------------------------///////////////////////////////////////////////// reading the coeficients and host */

/*   PetscLogDouble v1,v2,elapsed_time; */

/*   PetscTime(&v1); */
  
/*   PetscInt  hb=0; */
  
/*   for (bi=0; bi<block_number; bi++) { */
    
/*     DMDALocalInfo	info = user[bi].info; */
    
/*     PetscInt	xs = info.xs, xe = info.xs + info.xm; */
/*     PetscInt  	ys = info.ys, ye = info.ys + info.ym; */
/*     PetscInt	zs = info.zs, ze = info.zs + info.zm; */


/*     PetscInt d,dg,cll=0; */
/*     PetscInt itfc_number=0; */
/*     PetscReal lval[8]; */
/*     PetscInt	row=0,lcol[8]; */
/*     PetscInt col[8]; */
/*     PetscReal val[8]; */
/*     PetscReal x, y, z; */
/*     PetscInt  itfnumber=0,itfnumber1=0; */
/*     PetscInt d_freedom=0; */
/*     PetscInt coeff[block_number]; */
    
/*     /\**************************************************************************************************************************\/ */
/*     /\* Interpolate the Velocity on the interface nodes */
/*        from the host nodes. */
/*        itfc is the velocity at the cell corners */
/*        hostU is the velocity at the cell centers of the host block *\/ */
/*     /\**************************************************************************************************************************\/ */
    
        
/*     for (itfc_number=0; itfc_number<user[bi].itfcptsnumber; itfc_number++) { */
      
/*       x = user[bi].itfchostx[itfc_number]; */
/*       y = user[bi].itfchosty[itfc_number]; */
/*       z = user[bi].itfchostz[itfc_number]; */
      
      
/*       dg=0;d=0;cll=0; */
      
/*       for (i=0;i<8;i++){ */
/* 	lcol[i]=0; */
/* 	col[i]=0; */
/* 	lval[i]=0.0; */
/* 	val[i]=0.0; */
/*       } */
      
/*       cpu=0; */
/*       //itfnumber=0; */
         

/*       if(user[bi].itfchostx[itfc_number]>-1.25) */
/* 	if((user[bi].itfcK[itfc_number])>=zs&&(user[bi].itfcK[itfc_number])<ze) */
/* 	  if((user[bi].itfcJ[itfc_number])>=ys&&(user[bi].itfcJ[itfc_number])<ye) */
/* 	    if((user[bi].itfcI[itfc_number])>=xs&&(user[bi].itfcI[itfc_number])<xe) */
/* 	      { */
		
/* 		hb = user[bi].itfchostB[itfc_number]; */
/* 		/////////////////////////////////////// */
/* 		PetscInt find; */
/* 		for (find=0; find<block_number; find++) */
/* 		  { */
/* 		    if(hb>find) */
/* 		      coeff[find]=1; */
/* 		    else */
/* 		      coeff[find]=0; */
/* 		  } */
		  
/* 		  //////////////////////////////////// */
/* 		for(cpu=0;cpu<size;cpu++){ */
		  
/* 		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number]; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[0]=dg; */
/* 			  val[0]=(1-x) * (1-y) * (1-z); */
/* 			  itfnumber++; */
/* 			} */
		  
/* 		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number]+1; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[1]=dg; */
/* 			  val[1]=x * (1-y) * (1-z); */
/* 			  itfnumber++; */
/* 			} */
		  
		  
/* 		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number]; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[2]=dg; */
/* 			  val[2]=(1-x) * y * (1-z); */
/* 			  itfnumber++; */
/* 			} */
		  
		  
/* 		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number]; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[3]=dg; */
/* 			  val[3]=(1-x) * (1-y) * z; */
/* 			  itfnumber++; */
/* 			} */
	  
		  
/* 		  if((user[bi].itfchostK[itfc_number])>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number])<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number];j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number]+1; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[4]=dg; */
/* 			  val[4]= x * y * (1-z); */
/* 			  itfnumber++; */
/* 		} */
		  
		  
/* 		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number])>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number])<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number];i=user[bi].itfchostI[itfc_number]+1; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[5]=dg; */
/* 			  val[5]=x * (1-y) * z; */
/* 			  itfnumber++; */
/* 			} */
		  
		  
/* 		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number])>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number])<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number]; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[6]=dg; */
/* 			  val[6]=(1-x) * y * z; */
/* 			  itfnumber++; */
/* 			} */
		  
		  
/* 		  if((user[bi].itfchostK[itfc_number]+1)>=cpu_zs[hb][cpu]&&(user[bi].itfchostK[itfc_number]+1)<(cpu_zs[hb][cpu]+cpu_zm[hb][cpu])) */
/* 		    if((user[bi].itfchostJ[itfc_number]+1)>=cpu_ys[hb][cpu]&&(user[bi].itfchostJ[itfc_number]+1)<(cpu_ys[hb][cpu]+cpu_ym[hb][cpu])) */
/* 		      if((user[bi].itfchostI[itfc_number]+1)>=cpu_xs[hb][cpu]&&(user[bi].itfchostI[itfc_number]+1)<(cpu_xs[hb][cpu]+cpu_xm[hb][cpu])) */
/* 			{ */
/* 			  k=user[bi].itfchostK[itfc_number]+1;j=user[bi].itfchostJ[itfc_number]+1;i=user[bi].itfchostI[itfc_number]+1; */
/* 			  d=lidxLocal_matrix(i, j, k,hb,cpu,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm); */
/* 			  dg=dg_start[cpu]+coeff[0]*dl_start_blk0[cpu]+coeff[1]*dl_start_blk1[cpu]+3*d; */
/* 			  lcol[7]=dg; */
/* 			  val[7]=x * y * z; */
/* 			  itfnumber++; */
/* 			} */
		  
/* 		} */
		

/* 		/////////////////////////////////////// */
		
/* 		for (find=0; find<block_number; find++) */
/* 		  { */
/* 		    if(bi>find) */
/* 		      coeff[find]=1; */
/* 		    else */
/* 		      coeff[find]=0; */
/* 		  } */
		
/* 		///////////////////////////////////////////////////////////////////////////////////////// */
		
/* 		k=user[bi].itfcK[itfc_number];j=user[bi].itfcJ[itfc_number];i=user[bi].itfcI[itfc_number]; */
/* 		d=lidxLocal1_matrix(i, j, k, &user[bi],bi); */
/* 		dg=dg_start[rank]+coeff[0]*dl_start_blk0[rank]+coeff[1]*dl_start_blk1[rank]+3*d; */
		
/* 		for(d_freedom=0;d_freedom<3;d_freedom++){  */
/* 		  row=dg+d_freedom;//lidxLocal(8, 8, 8, &user,0); */
/* 		  itfnumber1++; */
		  
/* 		  for(i=0;i<8;i++){ */
/* 		    col[i]=lcol[i]+d_freedom; */
/* 		  } */
/* 		  // if(k==33&&j==39) */
/* 		  //PetscPrintf(PETSC_COMM_SELF, "row:%d col:%d %d %d %d %d %d %d %d!\n",row,col[0],col[1],col[2],col[3],col[4],col[5],col[6],col[7]); */
/* 		    //PetscPrintf(PETSC_COMM_SELF, "itfsearch:%d , col:%d %d %d %d %d %d %d %d -------row:%d!\n",itfc_number,col[0],col[1],col[2],col[3],col[4],col[5],col[6],col[7],row); */
/* 		  MatSetValues(Int_matrix,1,&row,8,col,val,INSERT_VALUES); */
		  
/* 		} */
		
/* 	      } */
/*     } */
    
/*     //PetscPrintf(PETSC_COMM_SELF, "bi:%d number: %d %d!\n",bi,itfnumber,itfnumber1); */
/*   } */


/*   MatAssemblyBegin(Int_matrix,MAT_FINAL_ASSEMBLY); */
/*   MatAssemblyEnd(Int_matrix,MAT_FINAL_ASSEMBLY); */


/*   PetscTime(&v2); */
/*   elapsed_time = v2 - v1; */
/*   //PetscPrintf(PETSC_COMM_WORLD,"elapsed_time:%le \n", elapsed_time); */

/* /\*   sprintf(filen, "Matfield.dat"); *\/ */
/* /\*   PetscViewerASCIIOpen(PETSC_COMM_WORLD, filen, &viewer); *\/ */
/* /\*   MatView(Int_matrix, viewer); *\/ */



/* /\*   /////////////////////////////////////////////////////////////////// *\/ */
/* /\*  ierr = DMCompositeGetAccess(packer,Umult,&Ub[0],&Ub[1]);CHKERRQ(ierr); *\/ */

/* /\*  //sprintf(filen, "umult%5.5d_%1.1d.dat", 0, 1); *\/ */
/* /\*   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"umult1" , &viewer); *\/ */
/* /\*   VecView(user[1].Ucat, viewer); *\/ */
/* /\*   PetscViewerDestroy(&viewer); *\/ */

/* /\*   //sprintf(filen, "umult%5.5d_%1.1d.dat", 0, 0); *\/ */
/* /\*   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"umult0", &viewer); *\/ */
/* /\*   VecView(user[0].Ucat, viewer); *\/ */
/* /\*   PetscViewerDestroy(&viewer); *\/ */

/* /\*   ierr = DMCompositeRestoreAccess(packer,Umult,&Ub[0],&Ub[1]);CHKERRQ(ierr); *\/ */
  
/* /\*   ////////////////////////////////////////////////////////////////////////////////////// *\/ */

/*   for (bi=0; bi<block_number; bi++) { */
/*     ierr = VecDestroy(&Ub[bi]);CHKERRQ(ierr); */
/*   } */
  
/*   //PetscViewerDestroy(&viewer) */

/*   ////////////////////////////////////////////////////////////////////// */
/*   PetscFree(dgl_start);PetscFree(dg_start);PetscFree(dl_start_blk0); PetscFree(dl_start_blk1);  PetscFree(dl_start_blk2); */
/*   PetscFree(cpu_lxs11);PetscFree(cpu_lys11);PetscFree(cpu_lzs11);PetscFree(cpu_lxm11);PetscFree(cpu_lym11);PetscFree(cpu_lzm11); */
/*   PetscFree(cpu_xs11);PetscFree(cpu_ys11);PetscFree(cpu_zs11);PetscFree(cpu_xm11);PetscFree(cpu_ym11);PetscFree(cpu_zm11); */
/*   PetscFree(cpu_lxs22);PetscFree(cpu_lys22);PetscFree(cpu_lzs22);PetscFree(cpu_lxm22);PetscFree(cpu_lym22);PetscFree(cpu_lzm22); */
/*   PetscFree(cpu_xs22);PetscFree(cpu_ys22);PetscFree(cpu_zs22);PetscFree(cpu_xm22);PetscFree(cpu_ym22);PetscFree(cpu_zm22); */
/*   PetscFree(cpu_lxs33);PetscFree(cpu_lys33);PetscFree(cpu_lzs33);PetscFree(cpu_lxm33);PetscFree(cpu_lym33);PetscFree(cpu_lzm33); */
/*   PetscFree(cpu_xs33);PetscFree(cpu_ys33);PetscFree(cpu_zs33);PetscFree(cpu_xm33);PetscFree(cpu_ym33);PetscFree(cpu_zm33); */

/*   ////////////////////////////////////////////////////////////////////// */


/*   PetscPrintf(PETSC_COMM_WORLD, "Done!\n"); */

/*   return 0; */
/* } */




/* PetscInt lidxLocal_matrix(PetscInt i, PetscInt j, PetscInt k,PetscInt blk,PetscInt cpu,PetscInt xs[block_number][cpu_size],PetscInt xm[block_number][cpu_size],PetscInt ys[block_number][cpu_size],PetscInt ym[block_number][cpu_size],PetscInt zs[block_number][cpu_size],PetscInt zm[block_number][cpu_size]) */
/* { */
  
/*   return ((k-zs[blk][cpu]) * (xm[blk][cpu]*ym[blk][cpu]) + (j-ys[blk][cpu])*(xm[blk][cpu]) + (i-xs[blk][cpu])); */
/* } */



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


PetscErrorCode Resistance(UserCtx *user)  //It is called in solver file
{

 PetscInt   i, j, k;
 PetscReal Q=0.0,pin=0.0,pout=0.0,Ql,plin,plout,R;
 PetscReal Arealin,Arealout,Areain=0.0,Areaout=0.0;
  Cmpnts    ***ucont,***csi, ***eta, ***zet;
 

  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert, ***pp; //local working array


  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->P, &pp);
  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  if (zs==0) { //Assumed this one is inlet
    k = 0;plin=0.0;Arealin=0.0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if (nvert[k+1][j][i]<0.1) {
	      Ql+= -ucont[k][j][i].z;
	      plin+=pp[k][j][i];
	      Arealin += sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			   (csi[k][j][i].y) * (csi[k][j][i].y) +
			   (csi[k][j][i].z) * (csi[k][j][i].z));
	    }
	  }
	}
  }
  MPI_Allreduce(&Ql,&Q,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&plin,&pin,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Arealin,&Areain,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  if (ze==mz) { //Assumed this one is outlet
    k = mz-2;plout=0.0;Arealout=0.0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if (nvert[k][j][i]<0.1) {
	     plout+=pp[k][j][i];
	     Arealout += sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			   (csi[k][j][i].y) * (csi[k][j][i].y) +
			   (csi[k][j][i].z) * (csi[k][j][i].z));
	    }
	  }
	}
  }
  MPI_Allreduce(&plout,&pout,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Arealout,&Areaout,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
 DMDAVecRestoreArray(fda, user->Ucont, &ucont);
 DMDAVecRestoreArray(da, user->lNvert, &nvert);
 DMDAVecRestoreArray(da, user->P, &pp);
 DMDAVecRestoreArray(fda, user->lCsi,  &csi);
 DMDAVecRestoreArray(fda, user->lEta,  &eta);
 DMDAVecRestoreArray(fda, user->lZet,  &zet);

 R=-(pin-pout)/Q;
 PetscPrintf(PETSC_COMM_WORLD, "Q=%le , pout=%le , pin=%le , R=%le\n",Q,pout,pin,R);
 PetscPrintf(PETSC_COMM_WORLD, "Ain=%le , Aout=%le\n",Areain,Areaout);
 return 0; 
}
