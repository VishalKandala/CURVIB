#include "variables.h"
extern PetscInt   block_number, ti,immersed,NumberOfBodies,STRONG_COUPLING,duplicate;

PetscReal     Sigma_yz=0.0, Sigma_xx=0.0, Sigma_yy=0.0, Sigma_zz=0.0;
PetscReal     S_yz=0.0,S_yzSum=0.0,S_xx,S_yy=0.0,S_zz=0.0,T_yz=0.0,S_xxSum=0.0,S_yySum=0.0,S_zzSum=0.0,T_yzSum=0.0; 
PetscReal     Rey_yz=0.0,Rey_xx=0.0,Rey_yy=0.0,Rey_zz=0.0;
PetscReal     AInertia_yz=0.0, Inertia_yz=0.0, Inertia_xx=0.0,  Inertia_yy=0.0, Inertia_zz=0.0;

PetscErrorCode Reynolds_Stress(UserCtx *user,FSInfo *FSinfo)
{

  DM             da = user->da, fda = user->fda;
  DMDALocalInfo  info = user->info;
  PetscInt	 xs = info.xs, xe = info.xs + info.xm;
  PetscInt       ys = info.ys, ye = info.ys + info.ym;
  PetscInt	 zs = info.zs, ze = info.zs + info.zm;
  PetscInt	 mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	 lxs, lxe, lys, lye, lzs, lze;
  PetscInt       i, j, k,ii;

  Cmpnts         ***ucat,***cent;
  PetscReal      ***nvert,***aj;
  PetscReal      U_inf=0.0,V_inf=0.0,W_inf=0.0;
 
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
      
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(fda, user->Cent, &cent);

  PetscReal   gammadot=2.0,U_bc;
  PetscOptionsGetReal(PETSC_NULL, "-gammadot", &gammadot, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-U_bc", &U_bc, PETSC_NULL);
  // PetscPrintf(PETSC_COMM_WORLD, "velocity of moving boundary is %le \n",U_bc);

  PetscReal  Rey_xx_f=0.0,Rey_xx_s=0.0;
  PetscReal  Rey_xx_fSum=0.0,Rey_xx_sSum=0.0;
  PetscReal  Rey_yy_f=0.0,Rey_yy_s=0.0;
  PetscReal  Rey_yy_fSum=0.0,Rey_yy_sSum=0.0;
  PetscReal  Rey_zz_f=0.0,Rey_zz_s=0.0;
  PetscReal  Rey_zz_fSum=0.0,Rey_zz_sSum=0.0;
  PetscReal  Rey_yz_f=0.0,Rey_yz_s=0.0;
  PetscReal  Rey_yz_fSum=0.0,Rey_yz_sSum=0.0;
 
  Cmpnts     omega_c,a_c,u_s;
  PetscReal  rx,ry,rz;
 

  //  PetscPrintf(PETSC_COMM_WORLD, "angular velosity of solid is %le %le %le \n",omega_c.x,omega_c.y,omega_c.z);
  

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	U_inf=0.0;
	V_inf=0.0;
	W_inf=(-U_bc+gammadot*cent[k][j][i].y);
	if (nvert[k][j][i]<0.9 ){
	  Rey_xx_f+=(ucat[k][j][i].x-U_inf)*(ucat[k][j][i].x-U_inf)/aj[k][j][i];
	  Rey_yy_f+=(ucat[k][j][i].y-V_inf)*(ucat[k][j][i].y-V_inf)/aj[k][j][i];
	  Rey_zz_f+=(ucat[k][j][i].z-W_inf)*(ucat[k][j][i].z-W_inf)/aj[k][j][i];
	  Rey_yz_f+=(ucat[k][j][i].y-V_inf)*(ucat[k][j][i].z-W_inf)/aj[k][j][i];
	}else {
	  if (nvert[k][j][i]<2.0) ii=(int)((nvert[k][j][i]-1.0)*1001);
	  else  ii=(int)((nvert[k][j][i]-3.0)*1001);
	

	  a_c.x=FSinfo[ii].x_c;
	  a_c.y=FSinfo[ii].y_c;
	  a_c.z=FSinfo[ii].z_c;

	  omega_c.x= FSinfo[ii].S_ang_n[1];
	  omega_c.y= FSinfo[ii].S_ang_n[3];
	  omega_c.z= FSinfo[ii].S_ang_n[5];
	  // PetscPrintf(PETSC_COMM_SELF, "center of mass for node i= %d ,j=%d ,k=%d  a_c.x %le a_c.y %le a_c.z %le \n",i,j,k,a_c.x,a_c.y,a_c.z);

	  rx = cent[k][j][i].x-a_c.x;
	  ry = cent[k][j][i].y-a_c.y;
	  rz = cent[k][j][i].z-a_c.z;
	  u_s.x =  rz*omega_c.y-omega_c.z*ry+FSinfo[ii].S_new[1];
	  u_s.y =  rx*omega_c.z-omega_c.x*rz+FSinfo[ii].S_new[3];
	  u_s.z =  ry*omega_c.x-omega_c.y*rx+FSinfo[ii].S_new[5];
	  Rey_xx_s+=(u_s.x-U_inf)*(u_s.x-U_inf)/aj[k][j][i];
	  Rey_yy_s+=(u_s.y-V_inf)*(u_s.y-V_inf)/aj[k][j][i];
	  Rey_zz_s+=(u_s.z-W_inf)*(u_s.z-W_inf)/aj[k][j][i];
	  Rey_yz_s+=(u_s.y-V_inf)*(u_s.z-W_inf)/aj[k][j][i];

	}
      }
    }
  }



  MPI_Allreduce(&Rey_xx_f,&Rey_xx_fSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Rey_yy_f,&Rey_yy_fSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Rey_zz_f,&Rey_zz_fSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Rey_yz_f,&Rey_yz_fSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_xx_f is %le  \n",Rey_xx_fSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_yy_f is %le  \n",Rey_yy_fSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_zz_f is %le  \n",Rey_zz_fSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_yz_f is %le  \n",Rey_yz_fSum); */
 
  MPI_Allreduce(&Rey_xx_s,&Rey_xx_sSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Rey_yy_s,&Rey_yy_sSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Rey_zz_s,&Rey_zz_sSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Rey_yz_s,&Rey_yz_sSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 /*  PetscPrintf(PETSC_COMM_WORLD, "Rey_xx_s is %le  \n",Rey_xx_sSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_yy_s is %le  \n",Rey_yy_sSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_zz_s is %le  \n",Rey_zz_sSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Rey_yz_s is %le  \n",Rey_yz_sSum); */

  Rey_xx=Rey_xx_fSum+Rey_xx_sSum;
  Rey_yy=Rey_yy_fSum+Rey_yy_sSum;
  Rey_zz=Rey_zz_fSum+Rey_zz_sSum;
  Rey_yz=Rey_yz_fSum+Rey_yz_sSum;

  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  return(0);
}

PetscErrorCode Inertia(UserCtx *user,FSInfo *FSinfo, PetscReal *Inertial_yz,PetscReal *Inertial_xx,PetscReal *Inertial_yy, PetscReal *Inertial_zz)

{

  DM             da = user->da, fda = user->fda;
  DMDALocalInfo  info = user->info;
  PetscInt	 xs = info.xs, xe = info.xs + info.xm;
  PetscInt       ys = info.ys, ye = info.ys + info.ym;
  PetscInt	 zs = info.zs, ze = info.zs + info.zm;
  PetscInt	 mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	 lxs, lxe, lys, lye, lzs, lze;
  PetscInt       i, j, k,ii;

  Cmpnts         ***cent,acc;
  PetscReal      ***nvert,***aj;
  PetscReal      r_z,r_y,r_x;
  PetscReal      X_c,Y_c,Z_c;
  PetscReal      w_x,w_y,w_z;
  PetscReal      a_x,a_y,a_z;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
      
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal *I_yz,*I_xx,*I_yy,*I_zz;


  I_yz=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
  I_xx=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
  I_yy=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
  I_zz=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(fda, user->Cent, &cent);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
  	if (nvert[k][j][i]>0.9 ){

	  if (nvert[k][j][i]<2.0) ii=(int)((nvert[k][j][i]-1.0)*1001);
	  else  ii=(int)((nvert[k][j][i]-3.0)*1001);
	  
	  
	  X_c=FSinfo[ii].x_c;
	  Y_c=FSinfo[ii].y_c;
	  Z_c=FSinfo[ii].z_c;
	  
	  w_x= FSinfo[ii].S_ang_n[1];
	  w_y= FSinfo[ii].S_ang_n[3];
	  w_z= FSinfo[ii].S_ang_n[5];

	  a_x= FSinfo[ii].alpha[0];
	  a_y= FSinfo[ii].alpha[1];
	  a_z= FSinfo[ii].alpha[2];

	  r_x=cent[k][j][i].x-X_c;
	  r_y=cent[k][j][i].y-Y_c;
	  r_z=cent[k][j][i].z-Z_c;
	
	
	  acc.x=r_z*a_y-r_y*a_z+r_y*w_x*w_y-r_x*w_y*w_y-r_x*w_z*w_z+r_z*w_x*w_z+FSinfo[ii].acc[0];
	  acc.y=r_x*a_z-r_z*a_x+r_x*w_x*w_y-r_y*w_x*w_x-r_y*w_z*w_z+r_z*w_y*w_z+FSinfo[ii].acc[1];
	  acc.z=r_y*a_x-r_x*a_y+r_x*w_x*w_z-r_z*w_x*w_x-r_z*w_y*w_y+r_y*w_y*w_z+FSinfo[ii].acc[2];
	  
	 /*  I_yz += 0.5*(acc.y*r_z+acc.z*r_y)/aj[k][j][i]; */
/* 	  J_yz += 0.5*(acc.y*r_z-acc.z*r_y)/aj[k][j][i]; */
/* 	  I_xx +=(acc.x*r_x)/aj[k][j][i]; */
/* 	  I_yy +=(acc.y*r_y)/aj[k][j][i]; */
/* 	  I_zz +=(acc.z*r_z)/aj[k][j][i]; */

	

	  I_yz[ii] += 0.5*(acc.y*r_z+acc.z*r_y)/aj[k][j][i];
	  I_xx[ii] +=(acc.x*r_x)/aj[k][j][i];
	  I_yy[ii] +=(acc.y*r_y)/aj[k][j][i];
	  I_zz[ii] +=(acc.z*r_z)/aj[k][j][i];
	}
      }
    }
  }

  MPI_Allreduce(I_xx,Inertial_xx,NumberOfBodies,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(I_yy,Inertial_yy,NumberOfBodies,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(I_zz,Inertial_zz,NumberOfBodies,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(I_yz,Inertial_yz,NumberOfBodies,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

/*   PetscPrintf(PETSC_COMM_WORLD, "Inertial_xx is %le  \n",Inertial_xx[0]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Inertial_yy is %le  \n",Inertial_yy[1]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Inertial_zz is %le  \n",Inertial_zz[2]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Inertial_yz is %le  \n",Inertial_yz[3]);  */

 /*  MPI_Allreduce(&I_xx,&Inertia_xx,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*   MPI_Allreduce(&I_yy,&Inertia_yy,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*   MPI_Allreduce(&I_zz,&Inertia_zz,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*   MPI_Allreduce(&I_yz,&Inertia_yz,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*   MPI_Allreduce(&J_yz,&AInertia_yz,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */


  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lAj, &aj);

  free(I_yz);
  free(I_xx);
  free(I_yy);
  free(I_zz);

  return 0;
}

PetscErrorCode Stresslet(UserCtx *user,FSInfo *FSinfo,IBMNodes *ibm, PetscInt ibi)
//PetscErrorCode Stresslet(UserCtx *user,FSInfo *FSinfo)

{
  DM             da = user->da, fda = user->fda;
  DMDALocalInfo  info = user->info;
  PetscInt	 xs = info.xs, xe = info.xs + info.xm;
  PetscInt       ys = info.ys, ye = info.ys + info.ym;
  PetscInt	 zs = info.zs, ze = info.zs + info.zm;
  PetscInt	 mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	 lxs, lxe, lys, lye, lzs, lze;
  PetscInt       i, j, k;
  
  Cmpnts         ***cent,***ucat;
  PetscReal      ***p;
  IBMInfo       *ibminfo;
  IBMListNode   *current;
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
      
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
 
  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(da, user->lP, &p);

  DMDAVecGetArray(fda, user->lCent, &cent);
  


      PetscReal     dwdz ,dvdz ,dudz ;
      PetscReal	    dwdy ,dvdy ,dudy ;
      PetscReal     dwdx ,dvdx ,dudx ;
     
      Cmpnts        ***csi,***eta,***zet;
      PetscReal     csi1,csi2,csi3;
      PetscReal     eta1,eta2,eta3;
      PetscReal     zet1,zet2,zet3;
      PetscReal     ***nvert,***iaj,***jaj,***kaj,***aj;
      PetscReal     Txx,Tyy,Tzz;
      PetscReal     Tzy,Tzx,Tyx;
      PetscReal     rei= 1./user->ren;
    
      PetscReal     r_z,z,r_y,y,r_x,x,X_c,Y_c,Z_c;
      PetscReal     A_x,A_y,A_z;

      PetscReal     Ap_x,Ap_y,Ap_z,Ap_xSum,Ap_ySum,Ap_zSum;
      PetscReal     An_x,An_y,An_z,An_xSum,An_ySum,An_zSum;

       
      Ap_x=0.;Ap_y=0.;Ap_z=0.;
      An_x=0.;An_y=0.;An_z=0.;
      Ap_xSum=0.;Ap_ySum=0.;Ap_zSum=0.;
      An_xSum=0.;An_ySum=0.;An_zSum=0.;

     
      S_xx=0.0,S_xxSum=0.0;
      S_yy=0.0,S_yySum=0.0;
      S_zz=0.0,S_zzSum=0.0;
      S_yz=0.0,S_yzSum=0.0;
      T_yz=0.0,T_yzSum=0.0;

      DMDAVecGetArray(da, user->lNvert, &nvert);
      DMDAVecGetArray(fda, user->lCsi, &csi);
      DMDAVecGetArray(fda, user->lEta, &eta);
      DMDAVecGetArray(fda, user->lZet, &zet);
   
      DMDAVecGetArray(da, user->lIAj, &iaj);
      DMDAVecGetArray(da, user->lJAj, &jaj);
      DMDAVecGetArray(da, user->lKAj, &kaj);
      DMDAVecGetArray(da, user->lAj, &aj);
   
   /*    for (k=lzs; k<lze; k++) { */
/* 	for (j=lys; j<lye; j++) { */
/* 	  for (i=lxs; i<lxe; i++) { */

      X_c=FSinfo->x_c; Y_c=FSinfo->y_c; Z_c=FSinfo->z_c;


      current = user->ibmlist[ibi].head;
      while (current) {
	ibminfo = &current->ibm_intp;
	current = current->next;
	i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk;
	if (i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze) {
	/*   if (nvert[k][j][i]<0.9 ) { */
	    
/* 	    V+=1.0/aj[k][j][i]; */
	    
/* 	    dvdz = 0.; */
/* 	    dwdy = 0.; */
	    
/* 	    dwdy =0.5*((ucat[k][j][i+1].z+ucat[k][j][i].z)*csi[k][j][i  ].y- */
/* 			 (ucat[k][j][i-1].z+ucat[k][j][i].z)*csi[k][j][i-1].y); */
/* 	    dvdz =0.5*((ucat[k][j][i+1].y+ucat[k][j][i].y)*csi[k][j][i  ].z- */
/* 			 (ucat[k][j][i-1].y+ucat[k][j][i].y)*csi[k][j][i-1].z); */
/* 	    dwdy+=0.5*((ucat[k][j+1][i].z+ucat[k][j][i].z)*eta[k][j  ][i].y- */
/* 		       (ucat[k][j-1][i].z+ucat[k][j][i].z)*eta[k][j-1][i].y); */
/* 	    dvdz+=0.5*((ucat[k][j+1][i].y+ucat[k][j][i].y)*eta[k][j  ][i].z- */
/* 		       (ucat[k][j-1][i].y+ucat[k][j][i].y)*eta[k][j-1][i].z); */
/* 	    dwdy+=0.5*((ucat[k+1][j][i].z+ucat[k][j][i].z)*zet[k][j  ][i].y- */
/* 		       (ucat[k-1][j][i].z+ucat[k][j][i].z)*zet[k-1][j][i].y); */
/* 	    dvdz+=0.5*((ucat[k+1][j][i].y+ucat[k][j][i].y)*zet[k  ][j][i].z- */
/* 		         (ucat[k-1][j][i].y+ucat[k][j][i].y)*zet[k-1][j][i].z); */
/* 	    dWdy_f +=dwdy; */
/* 	    dVdz_f +=dvdz; */
	    
/* 	  }else if (nvert[k][j][i]>0.9 && nvert[k][j][i]<1.1){ */
	    
/* 	    dvdz = 0.; */
/* 	    dwdy = 0.; */
	    
/* 	    if (nvert[k+1][j][i]<0.9){ */
/* 	      dwdy+=0.5*(ucat[k+1][j][i].z+ucat[k][j][i].z)*zet[k][j][i].y; */
/* 	      dvdz+=0.5*(ucat[k+1][j][i].y+ucat[k][j][i].y)*zet[k][j][i].z; */
/* 	    } */
/* 	    if ( nvert[k-1][j][i]<0.9){ */
	      
/* 	      dwdy-=0.5*(ucat[k-1][j][i].z+ucat[k][j][i].z)*zet[k-1][j][i].y; */
/* 	      dvdz-=0.5*(ucat[k-1][j][i].y+ucat[k][j][i].y)*zet[k-1][j][i].z; */
/* 	    } */
/* 	    if ( nvert[k][j+1][i]<0.9){ */
/* 	      dwdy+=0.5*(ucat[k][j+1][i].z+ucat[k][j][i].z)*eta[k][j][i].y; */
/* 	      dvdz+=0.5*(ucat[k][j+1][i].y+ucat[k][j][i].y)*eta[k][j][i].z; */
/* 	    } */
/* 	    if (nvert[k][j-1][i]<0.9){ */
/* 		dwdy-=0.5*(ucat[k][j-1][i].z+ucat[k][j][i].z)*eta[k][j-1][i].y; */
/* 		dvdz-=0.5*(ucat[k][j-1][i].y+ucat[k][j][i].y)*eta[k][j-1][i].z; */
/* 	    } */
/* 	    if ( nvert[k][j][i+1]<0.9){ */
/* 	      dwdy+=0.5*(ucat[k][j][i+1].z+ucat[k][j][i].z)*csi[k][j][i].y; */
/* 	      dvdz+=0.5*(ucat[k][j][i+1].y+ucat[k][j][i].y)*csi[k][j][i].z; */
/* 	    } */
/* 	    if ( nvert[k][j][i-1]<0.9){ */
/* 	      dwdy-=0.5*(ucat[k][j][i-1].z+ucat[k][j][i].z)*csi[k][j][i-1].y; */
/* 	      dvdz-=0.5*(ucat[k][j][i-1].y+ucat[k][j][i].y)*csi[k][j][i-1].z; */
/* 	    } */
	    
/* 	    dWdy_s +=dwdy; */
/* 	    dVdz_s +=dvdz; */
	    
	    /*   Min_dist=1.e10; */
	    /* 	      for (ibi=0;ibi<NumberOfBodies;ibi++){ */
	    /* 		r_x=(cent[k][j][i].x-FSinfo[ibi].x_c); */
/* 		r_y=(cent[k][j][i].y-FSinfo[ibi].y_c); */
/* 		r_z=(cent[k][j][i].z-FSinfo[ibi].z_c); */
/* 		distance=r_x*r_x+r_y*r_y+r_z*r_z; */
/* 		if (distance<Min_dist) { */
/* 		  ii=ibi; */
/* 		  Min_dist=distance; */
/* 		} */
/* 	      } */
	      
	  /*   X_c=FSinfo[ii].x_c; */
/* 	    Y_c=FSinfo[ii].y_c; */
/* 	    Z_c=FSinfo[ii].z_c; */
	    


	  dwdz = 0.;
	  dvdz = 0.;
	  dudz = 0.;
	  
	  dwdy = 0.;
	  dvdy = 0.;
	  dudy = 0.;
	  
	  dwdx = 0.;
	  dvdx = 0.;
	  dudx = 0.;
	  
	  if ((user->bctype[4]==7 && nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5) || (nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5 && k+1<mz-1 && k-1>0)) {
	// if (nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5 && k<zend && k>zstr) {
	    zet1 = 0.25*(zet[k][j][i].x+zet[k-1][j][i].x)*
	      (kaj[k][j][i]+kaj[k-1][j][i]);
	    zet2 = 0.25*(zet[k][j][i].y+zet[k-1][j][i].y)*
	      (kaj[k][j][i]+kaj[k-1][j][i]);
	    zet3 = 0.25*(zet[k][j][i].z+zet[k-1][j][i].z)*
	      (kaj[k][j][i]+kaj[k-1][j][i]);
	    dwdz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet3;
	    dvdz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet3;
	    dudz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet3;
	    
	    dwdy = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet2;
	    dvdy = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet2;
	    dudy = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet2;
	    
	    dwdx = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet1;
	    dvdx = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet1;
	    dudx = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet1;
	    
	  } else if (nvert[k+1][j][i]<2.5 && k+1<mz-1) {
	    // } else if (nvert[k+1][j][i]<2.5 && k<zend) {
	    zet1 = (zet[k][j][i].x)*kaj[k][j][i];
	    zet2 = (zet[k][j][i].y)*kaj[k][j][i];
	    zet3 = (zet[k][j][i].z)*kaj[k][j][i];
	    
	    dwdz = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet3;
	    dvdz = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet3;
	    dudz = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet3;
	    
	    dwdy = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet2;
	    dvdy = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet2;
	    dudy = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet2;
	    
	    dwdx = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet1;
	    dvdx = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet1;
	    dudx = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet1;
	    
	  } else if (nvert[k-1][j][i]<2.5 && k-1>0){
	    // } else if (nvert[k-1][j][i]<2.5 && k>zstr){
	    zet1 = (zet[k-1][j][i].x)*kaj[k-1][j][i];
	    zet2 = (zet[k-1][j][i].y)*kaj[k-1][j][i];
	    zet3 = (zet[k-1][j][i].z)*kaj[k-1][j][i];
	    
	    dwdz = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet3;
	    dvdz = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet3;
	    dudz = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet3;
	    
	    dwdy = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet2;
	    dvdy = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet2;
	    dudy = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet2;
	    
	    dwdx = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet1;
	    dvdx = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet1;
	    dudx = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet1;
	    
	  } else {
	    dwdz = 0.;
	    dvdz = 0.;
	    dudz = 0.;
	    
	    dwdy = 0.;
	    dvdy = 0.;
	    dudy = 0.;
	    
	    dwdx = 0.;
	    dvdx = 0.;
	    dudx = 0.;
	  }
	  
	  if ((user->bctype[2]==7 && nvert[k][j+1][i]<2.5 && nvert[k][j-1][i]<2.5)  || (nvert[k][j+1][i]<2.5 && j+1<my-1 && nvert[k][j-1][i]<2.5 && j-1>0)) {
	    // if (nvert[k][j+1][i]<2.5 && j<yend && nvert[k][j-1][i]<2.5 && j>ystr) {
	    eta1 = 0.25*(eta[k][j][i].x+eta[k][j-1][i].x)*
	      (jaj[k][j][i]+jaj[k][j-1][i]);
	    eta2 = 0.25*(eta[k][j][i].y+eta[k][j-1][i].y)*
	      (jaj[k][j][i]+jaj[k][j-1][i]);
	    eta3 = 0.25*(eta[k][j][i].z+eta[k][j-1][i].z)*
	      (jaj[k][j][i]+jaj[k][j-1][i]);
	    dwdz += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta3;
	    dvdz += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta3;
	    dudz += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta3;
	    
	    dwdy += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta2;
	    dvdy += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta2;
	    dudy += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta2;
	    
	    dwdx += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta1;
	    dvdx += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta1;
	    dudx += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta1;
	    
	  } else if (nvert[k][j+1][i]<2.5 && j+1<my-1) {
	    // } else if (nvert[k][j+1][i]<2.5 && j<yend) {
	    eta1 = eta[k][j][i].x*jaj[k][j][i];
	    eta2 = eta[k][j][i].y*jaj[k][j][i];
	    eta3 = eta[k][j][i].z*jaj[k][j][i];
	    
	    dwdz += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta3;
	    dvdz += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta3;
	    dudz += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta3;
	    
	    dwdy += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta2;
	    dvdy += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta2;
	    dudy += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta2;
	    
	    dwdx += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta1;
	    dvdx += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta1;
	    dudx += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta1;
	    
	  } else if (nvert[k][j-1][i]<2.5 && j-1>0){
	    // } else if (nvert[k][j-1][i]<2.5 && j>ystr){
	    eta1 = eta[k][j-1][i].x*jaj[k][j-1][i];
	    eta2 = eta[k][j-1][i].y*jaj[k][j-1][i];
	    eta3 = eta[k][j-1][i].z*jaj[k][j-1][i];
	    
	    dwdz += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta3;
	    dvdz += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta3;
	    dudz += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta3;
	    
	    dwdy += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta2;
	    dvdy += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta2;
	    dudy += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta2;
	    
	    dwdx += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta1;
	    dvdx += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta1;
	    dudx += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta1;
	  } 
	  if ((user->bctype[0]==7 && nvert[k][j][i+1]<2.5 && nvert[k][j][i-1]<2.5) || (nvert[k][j][i+1]<2.5 && i+1<mx-1 && nvert[k][j][i-1]<2.5 && i-1>0)) {
	//if (nvert[k][j][i+1]<2.5 && i<xend && nvert[k][j][i-1]<2.5 && i>xstr) {
	    csi1 = 0.25*(csi[k][j][i].x+csi[k][j][i-1].x)*
	      (iaj[k][j][i]+iaj[k][j][i-1]);
	    csi2 = 0.25*(csi[k][j][i].y+csi[k][j][i-1].y)*
	      (iaj[k][j][i]+iaj[k][j][i-1]);
	    csi3 = 0.25*(csi[k][j][i].z+csi[k][j][i-1].z)*
	      (iaj[k][j][i]+iaj[k][j][i-1]);
	    
	    dwdz += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi3;
	    dvdz += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi3;
	    dudz += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi3;
	    
	    dwdy += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi2;
	    dvdy += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi2;
	    dudy += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi2;
	    
	    dwdx += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi1;
	    dvdx += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi1;
	    dudx += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi1;
	    
	  } else if (nvert[k][j][i+1]<2.5 && i+1<mx-1) {
	    // } else if (nvert[k][j][i+1]<2.5 && i<xend) {
	    csi1 = csi[k][j][i].x*iaj[k][j][i];
	    csi2 = csi[k][j][i].y*iaj[k][j][i];
	    csi3 = csi[k][j][i].z*iaj[k][j][i];
	    
	    dwdz += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi3;
	    dvdz += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi3;
	    dudz += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi3;
	    
	    dwdy += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi2;
	    dvdy += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi2;
	    dudy += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi2;
	    
	    dwdx += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi1;
	    dvdx += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi1;
	    dudx += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi1;
	    
	  } else if (nvert[k][j][i-1]<2.5 && i-1>0){
	    // } else if (nvert[k][j][i-1]<2.5 && i>xstr){
	    csi1 = csi[k][j][i-1].x*iaj[k][j][i-1];
	    csi2 = csi[k][j][i-1].y*iaj[k][j][i-1];
	    csi3 = csi[k][j][i-1].z*iaj[k][j][i-1];
	    
	    dwdz += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi3;
	    dvdz += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi3;
	    dudz += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi3;
	    
	    dwdy += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi2;
	    dvdy += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi2;
	    dudy += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi2;
	    
	    dwdx += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi1;
	    dvdx += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi1;
	    dudx += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi1;
	  }
	  
	  Tzz = rei * (dwdz + dwdz);
	  Tyy = rei * (dvdy + dvdy);
	  Txx = rei * (dudx + dudx);
	  Tzy = rei * (dwdy + dvdz);
	  Tzx = rei * (dwdx + dudz);
	  Tyx = rei * (dvdx + dudy);
	  


	  r_x = 0.;
	  r_y = 0.;
	  r_z = 0.;

	  /*   K +  */

	  if (nvert[k+1][j][i]<0.9 ){
	      
	    z = 0.5*(cent[k][j][i].z+cent[k+1][j][i].z);
	    y = 0.5*(cent[k][j][i].y+cent[k+1][j][i].y);
	    x = 0.5*(cent[k][j][i].x+cent[k+1][j][i].x);

	    r_z = z-Z_c;
	    r_y = y-Y_c;
	    r_x = x-X_c;
	    
	    A_z=zet[k  ][j][i].z;
	    A_y=zet[k  ][j][i].y;
	    A_x=zet[k  ][j][i].x;
	    
	    Ap_z  += A_z;
	    Ap_y  += A_y;
	    Ap_x  += A_x;
	    
	 /*    S_xx += Tzx*A_z; */
/* 	    S_yy += Tzy*A_z; */
/* 	    S_zz +=(-p[k+1][j][i]+Tzz))*A_z; */
/* 	    S_yz +=0.5*(r_z*Tzy+r_y*(-p[k+1][j][i]+Tzz))*A_z; */

	    S_xx +=(r_x*Tzx)*A_z;
	    S_yy +=(r_y*Tzy)*A_z;
	    S_zz +=(r_z*(-p[k][j][i]+Tzz))*A_z;
	    S_yz +=0.5*(r_z*Tzy+r_y*(-p[k][j][i]+Tzz))*A_z;
	    //  T_yz +=0.5*(r_z*Tzy-r_y*(-p[k+1][j][i]+Tzz))*A_z;
	      //	S_yz +=(r_z*Tzy-rei*0.5*(ucat[k+1][j][i].y+ucat[k][j][i].y))*A_z;
	      // Sigma_yz +=r_z*Tyx*A_x+(r_z*Tyy- rei*0.5*(ucat[k+1][j][i].z+ucat[k][j][i].z))*A_y+(r_z*Tzy-rei*0.5*(ucat[k+1][j][i].y+ucat[k][j][i].y))*A_z;
	      
	  }
	    /* 	   K - */
	  if (nvert[k-1][j][i]<0.9) {
	    
	    z = 0.5*(cent[k][j][i].z+cent[k-1][j][i].z);
	    y = 0.5*(cent[k][j][i].y+cent[k-1][j][i].y);
	    x = 0.5*(cent[k][j][i].x+cent[k-1][j][i].x);
	    r_z = z-Z_c;
	    r_y = y-Y_c;
	    r_x = x-X_c;
	    
	    A_z=zet[k-1][j][i].z;
	    A_y=zet[k-1][j][i].y;
	    A_x=zet[k-1][j][i].x;
	    
	    An_z  += A_z;
	    An_y  += A_y;
	    An_x  += A_x;
	    
	   /*  S_xx -= Tzx*A_z; */
/* 	    S_yy -= Tzy*A_z; */
/* 	    S_zz -= (-p[k-1][j][i]+Tzz)*A_z; */
/* 	    S_yz -=0.5*(r_z*Tzy+r_y*(-p[k-1][j][i]+Tzz))*A_z; */

	    S_xx -=(r_x*Tzx)*A_z;
	    S_yy -=(r_y*Tzy)*A_z;
	    S_zz -=(r_z*(-p[k][j][i]+Tzz))*A_z;
	    S_yz -=0.5*(r_z*Tzy+r_y*(-p[k][j][i]+Tzz))*A_z;
	    //  T_yz -=0.5*(r_z*Tzy-r_y*(-p[k][j][i]+Tzz))*A_z;
	    //	S_yz -=(r_z*Tzy-rei*0.5*(ucat[k-1][j][i].y+ucat[k][j][i].y))*A_z;
	    // Sigma_yz -=r_z*Tyx*A_x+(r_z*Tyy- rei*0.5*(ucat[k-1][j][i].z+ucat[k][j][i].z))*A_y+(r_z*Tzy-rei*0.5*(ucat[k-1][j][i].y+ucat[k][j][i].y))*A_z;;
	  }
	  
	  /*       j +  */
	  if (nvert[k][j+1][i]<0.9 ){
	    
	    z = 0.5*(cent[k][j][i].z+cent[k][j+1][i].z);
	    y = 0.5*(cent[k][j][i].y+cent[k][j+1][i].y);
	    x = 0.5*(cent[k][j][i].x+cent[k][j+1][i].x);
	    r_z = z-Z_c;
	    r_y = y-Y_c;
	    r_x = x-X_c;
	    
	    A_z=eta[k][j  ][i].z;
	    A_y=eta[k][j  ][i].y;
	    A_x=eta[k][j  ][i].x;
	    
	    
	    Ap_z  += A_z;
	    Ap_y  += A_y;
	    Ap_x  += A_x;
	    
/* 	    S_xx += Tyx*A_y; */
/* 	    S_yy += (-p[k][j+1][i]+Tyy)*A_y; */
/* 	    S_zz += Tzy*A_y; */
/* 	    S_yz +=0.5*(r_z*(-p[k][j+1][i]+Tyy)+r_y*Tzy)*A_y; */

	    S_xx +=(r_x*Tyx)*A_y;
	    S_yy +=(r_y*(-p[k][j][i]+Tyy))*A_y;
	    S_zz +=(r_z*Tzy)*A_y;
	    S_yz +=0.5*(r_z*(-p[k][j][i]+Tyy)+r_y*Tzy)*A_y;
	    //  T_yz +=0.5*(r_z*(-p[k][j][i]+Tyy)-r_y*Tzy)*A_y;
	    //	S_yz +=(r_z*(Tyy-p[k][j+1][i])-rei*0.5*(ucat[k][j+1][i].z+ucat[k][j][i].z))*A_y;
	    // Sigma_yz +=r_z*Tyx*A_x+(r_z*Tyy- rei*0.5*(ucat[k][j+1][i].z+ucat[k][j][i].z))*A_y+(r_z*Tzy-rei*0.5*(ucat[k][j+1][i].y+ucat[k][j][i].y))*A_z;
	  }
	  /*       j -  */
	    
	  if (nvert[k][j-1][i]<0.9) {
	    
	    z = 0.5*(cent[k][j][i].z+cent[k][j-1][i].z);
	    y = 0.5*(cent[k][j][i].y+cent[k][j-1][i].y);
	    x = 0.5*(cent[k][j][i].x+cent[k][j-1][i].x);
	    r_z = z-Z_c;
	    r_y = y-Y_c;
	    r_x = x-X_c;
	    
	    
	    A_z=eta[k][j-1][i].z;
	    A_y=eta[k][j-1][i].y;
	    A_x=eta[k][j-1][i].x;
	    
	    An_z  += A_z;
	    An_y  += A_y;
	    An_x  += A_x;
	    
	  /*   S_xx -= Tyx*A_y; */
/* 	    S_yy -=(-p[k][j-1][i]+Tyy)*A_y; */
/* 	    S_zz -= Tzy*A_y; */
/* 	    S_yz -=0.5*(r_z*(-p[k][j-1][i]+Tyy)+r_y*Tzy)*A_y; */
	   
	    S_xx -=(r_x*Tyx)*A_y;
	    S_yy -=(r_y*(-p[k][j][i]+Tyy))*A_y;
	    S_zz -=(r_z*Tzy)*A_y;
	    S_yz -=0.5*(r_z*(-p[k][j][i]+Tyy)+r_y*Tzy)*A_y;
	    //  T_yz -=0.5*(r_z*(-p[k][j][i]+Tyy)-r_y*Tzy)*A_y;
	    //	S_yz -=(r_z*(Tyy-p[k][j-1][i])-rei*0.5*(ucat[k][j-1][i].z+ucat[k][j][i].z))*A_y;
	    // Sigma_yz -=r_z*Tyx*A_x+(r_z*Tyy- rei*0.5*(ucat[k][j-1][i].z+ucat[k][j][i].z))*A_y+(r_z*Tzy-rei*0.5*(ucat[k][j-1][i].y+ucat[k][j][i].y))*A_z;;;
	  }
	  /*       i +  */
	  
	  if (nvert[k][j][i+1]<0.9){
	    
	    z = 0.5*(cent[k][j][i].z+cent[k][j][i+1].z);
	    y = 0.5*(cent[k][j][i].y+cent[k][j][i+1].y);
	    x = 0.5*(cent[k][j][i].x+cent[k][j][i+1].x);

	    r_z = z-Z_c;
	    r_y = y-Y_c;
	    r_x = x-X_c;
	    
	    
	    A_z=csi[k][j][i].z;
	    A_y=csi[k][j][i].y;
	    A_x=csi[k][j][i].x;
	    
	    Ap_z  += A_z;
	    Ap_y  += A_y;
	    Ap_x  += A_x;
	    
	   /*  S_xx += (-p[k][j][i+1]+Txx)*A_x; */
/* 	    S_yy += Tyx*A_x; */
/* 	    S_zz += Tzx*A_x; */
/* 	    S_yz +=0.5*(r_z*Tyx+r_y*Tzx)*A_x; */
	    S_xx +=(r_x*(-p[k][j][i]+Txx))*A_x;
	    S_yy +=(r_y*Tyx)*A_x;
	    S_zz +=(r_z*Tzx)*A_x;
	    S_yz +=0.5*(r_z*Tyx+r_y*Tzx)*A_x;
	    // T_yz +=0.5*(r_z*Tyx-r_y*Tzx)*A_x;
	    //	S_yz +=r_z*Tyx*A_x;
	    //Sigma_yz +=r_z*Tyx*A_x+(r_z*Tyy- rei*0.5*(ucat[k][j][i+1].z+ucat[k][j][i].z))*A_y+(r_z*Tzy-rei*0.5*(ucat[k][j][i+1].y+ucat[k][j][i].y))*A_z;
	  }
	  /*       i -  */
	  if  (nvert[k][j][i-1]<0.9) {
	    
	    z = 0.5*(cent[k][j][i].z+cent[k][j][i-1].z);
	    y = 0.5*(cent[k][j][i].y+cent[k][j][i-1].y);
	    x = 0.5*(cent[k][j][i].x+cent[k][j][i-1].x);

	    r_z = z-Z_c;
	    r_y = y-Y_c;
	    r_x = x-X_c;
	    
	    A_z=csi[k][j][i-1].z;
	    A_y=csi[k][j][i-1].y;
	    A_x=csi[k][j][i-1].x;
	    
	    An_z  += A_z;
	    An_y  += A_y;
	    An_x  += A_x;
	    
	 /*    S_xx -= (-p[k][j][i-1]+Txx)*A_x; */
/* 	    S_yy -= Tyx*A_x; */
/* 	    S_zz -= Tzx*A_x; */
/* 	    S_yz -=0.5*(r_z*Tyx+r_y*Tzx)*A_x; */
	    S_xx -=(r_x*(-p[k][j][i]+Txx))*A_x;
	    S_yy -=(r_y*Tyx)*A_x;
	    S_zz -=(r_z*Tzx)*A_x;
	    S_yz -=0.5*(r_z*Tyx+r_y*Tzx)*A_x;
	    //   T_yz -=0.5*(r_z*Tyx-r_y*Tzx)*A_x;
	      //	S_yz -=r_z*Tyx*A_x;
	      // Sigma_yz -=r_z*Tyx*A_x+(r_z*Tyy- rei*0.5*(ucat[k][j][i-1].z+ucat[k][j][i].z))*A_y+(r_z*Tzy-rei*0.5*(ucat[k][j][i-1].y+ucat[k][j][i].y))*A_z;
	  }
	}
      }

/* 	  } */
/* 	} */
/*       } */

     /*  MPI_Allreduce(&V,&VSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       PetscPrintf(PETSC_COMM_WORLD, "V_f is %le  \n",VSum); */
     
   /*    MPI_Allreduce(&dWdy_f,&dWdySum_f,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&dVdz_f,&dVdzSum_f,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&dWdy_s,&dWdySum_s,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&dVdz_s,&dVdzSum_s,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */

/*       PetscPrintf(PETSC_COMM_WORLD, "dwdy of fluid is %le  \n",dWdySum_f); */
/*       PetscPrintf(PETSC_COMM_WORLD, "dwdy of solid is %le  \n",dWdySum_s); */

/*       PetscPrintf(PETSC_COMM_WORLD, "dvdz of fluid is %le  \n",dVdzSum_f); */
/*       PetscPrintf(PETSC_COMM_WORLD, "dvdz of solid is %le  \n",dVdzSum_s); */

   /*    MPI_Allreduce(&Ap_x, &Ap_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&Ap_y, &Ap_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&Ap_z, &Ap_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
      
/*       MPI_Allreduce(&An_x, &An_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&An_y, &An_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*       MPI_Allreduce(&An_z, &An_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */

/*       PetscPrintf(PETSC_COMM_WORLD, "Ap_x,An_x, %le ,%le ,Ap_y,An_y, %le ,%le ,Ap_z,An_z, %le, %le \n",Ap_xSum,An_xSum,Ap_ySum,An_ySum,Ap_zSum,An_zSum); */

      MPI_Allreduce(&S_xx,&S_xxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&S_yy,&S_yySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&S_zz,&S_zzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&S_yz,&S_yzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&T_yz,&T_yzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      
      /*  PetscPrintf(PETSC_COMM_WORLD, "S_xx is %le  \n",S_xxSum); */
      /*       PetscPrintf(PETSC_COMM_WORLD, "S_yy is %le  \n",S_yySum); */
      /*       PetscPrintf(PETSC_COMM_WORLD, "S_zz is %le  \n",S_zzSum); */
      /*       PetscPrintf(PETSC_COMM_WORLD, "S_yz is %le  \n",S_yzSum); */
      //  PetscPrintf(PETSC_COMM_WORLD, "T_yz is %le  \n",T_yzSum);
      
      
      DMDAVecRestoreArray(da, user->lNvert, &nvert);
      DMDAVecRestoreArray(fda, user->lCsi, &csi);
      DMDAVecRestoreArray(fda, user->lEta, &eta);
      DMDAVecRestoreArray(fda, user->lZet, &zet);
      DMDAVecRestoreArray(fda, user->lCent, &cent);
      DMDAVecRestoreArray(da, user->lIAj, &iaj);
      DMDAVecRestoreArray(da, user->lJAj, &jaj);
      DMDAVecRestoreArray(da, user->lKAj, &kaj);
      DMDAVecRestoreArray(da, user->lAj, &aj);
          
      DMDAVecRestoreArray(da, user->lP, &p);
      DMDAVecRestoreArray(fda, user->lUcat, &ucat);
   
      return(0);
}
PetscErrorCode Viscosity(UserCtx *user,FSInfo *FSinfo,IBMNodes *ibm){


  PetscInt ibi,P_No;
  PetscReal *Stresslet_yz,*Stresslet_xx,*Stresslet_yy,*Stresslet_zz;
  PetscReal *Inertial_yz,*Inertial_xx,*Inertial_yy,*Inertial_zz;
  
  Stresslet_yz= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));
  Stresslet_xx= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));
  Stresslet_yy= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));
  Stresslet_zz= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));

  Inertial_yz= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));
  Inertial_xx= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));
  Inertial_yy= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));
  Inertial_zz= (PetscReal *) calloc(NumberOfBodies,sizeof(PetscReal));

  for (ibi=0;ibi<NumberOfBodies;ibi++){
   /*  Stresslet_yz[ibi]=0.0; */
/*     Stresslet_xx[ibi]=0.0; */
/*     Stresslet_yy[ibi]=0.0; */
/*     Stresslet_zz[ibi]=0.0; */
    Stresslet(user,&FSinfo[ibi],&ibm[ibi],ibi);
    Stresslet_yz[ibi]=S_yzSum;
    Stresslet_xx[ibi]=S_xxSum;
    Stresslet_yy[ibi]=S_yySum;
    Stresslet_zz[ibi]=S_zzSum;
    PetscPrintf(PETSC_COMM_WORLD, "S_yz[%d] is %le  \n",ibi,Stresslet_yz[ibi]);
  }

  S_yzSum=0.0;
  S_xxSum=0.0;
  S_yySum=0.0;
  S_zzSum=0.0;

  Inertia_xx=0.0;
  Inertia_yy=0.0;
  Inertia_zz=0.0;
  Inertia_yz=0.0;

  Reynolds_Stress(user,FSinfo); 
  
  Inertia(user,FSinfo,Inertial_yz,Inertial_xx,Inertial_yy,Inertial_zz);

  for (ibi=0;ibi<NumberOfBodies;ibi++){
    
    S_yzSum +=Stresslet_yz[ibi];
    S_xxSum +=Stresslet_xx[ibi];
    S_yySum +=Stresslet_yy[ibi];
    S_zzSum +=Stresslet_zz[ibi];
    
    Inertia_xx +=Inertial_xx[ibi];
    Inertia_yy +=Inertial_yy[ibi];
    Inertia_zz +=Inertial_zz[ibi];
    Inertia_yz +=Inertial_yz[ibi];
    
    if (FSinfo[ibi].clone > 1.){
      P_No=(int)((FSinfo[ibi].clone-(int)FSinfo[ibi].clone)*1001)-1;
      
      Inertial_yz[P_No] +=Inertial_yz[ibi];
      Inertial_xx[P_No] +=Inertial_xx[ibi];
      Inertial_yy[P_No] +=Inertial_yy[ibi];
      Inertial_zz[P_No] +=Inertial_zz[ibi];

      Inertial_yz[ibi]  =Inertial_yz[P_No];
      Inertial_xx[ibi]  =Inertial_xx[P_No];
      Inertial_yy[ibi]  =Inertial_yy[P_No];
      Inertial_zz[ibi]  =Inertial_zz[P_No];

      Stresslet_yz[P_No] +=Stresslet_yz[ibi];
      Stresslet_xx[P_No] +=Stresslet_xx[ibi];
      Stresslet_yy[P_No] +=Stresslet_yy[ibi];
      Stresslet_zz[P_No] +=Stresslet_zz[ibi];
      
      Stresslet_yz[ibi]  =Stresslet_yz[P_No];
      Stresslet_xx[ibi]  =Stresslet_xx[P_No];
      Stresslet_yy[ibi]  =Stresslet_yy[P_No];
      Stresslet_zz[ibi]  =Stresslet_zz[P_No];
    }
     
  }

 

 

 /*  Sigma_xx=S_xxSum-Rey_xx-Inertia_xx; */
/*   Sigma_yy=S_yySum-Rey_yy-Inertia_yy; */
/*   Sigma_zz=S_zzSum-Rey_zz-Inertia_zz; */
/*   Sigma_yz=S_yzSum+T_yzSum-Rey_yz-Inertia_yz-AInertia_yz; */
 
/*   PetscPrintf(PETSC_COMM_WORLD, "Sigma_xx is %le  \n",S_xxSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Sigma_yy is %le  \n",S_yySum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Sigma_zz is %le  \n",S_zzSum); */
/*   PetscPrintf(PETSC_COMM_WORLD, "Sigma_yz is %le  \n",S_yzSum); */

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Sigma_rheo");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le %le %le \n",ti,S_yzSum,Inertia_yz,Rey_yz,S_xxSum,Inertia_xx,Rey_xx,S_yySum,Inertia_yy,Rey_yy,S_zzSum,Inertia_zz,Rey_zz);
    fclose(f); 
    for (ibi=0;ibi<NumberOfBodies;ibi++){
      sprintf(filen, "Sigma_rheo%2.2d",ibi);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD,f ,"%d %le %le %le %le %le %le %le %le %le %le %le %le \n",ti,Stresslet_yz[ibi],Inertial_yz[ibi],Rey_yz,Stresslet_xx[ibi],Inertial_xx[ibi],Rey_xx,Stresslet_yy[ibi],Inertial_yy[ibi],Rey_yy,Stresslet_zz[ibi],Inertial_zz[ibi],Rey_zz);
    }
    fclose(f);
  }

  free(Stresslet_yz);
  free(Stresslet_xx);
  free(Stresslet_yy);
  free(Stresslet_zz);

  free(Inertial_yz);
  free(Inertial_xx);
  free(Inertial_yy);
  free(Inertial_zz);

  return 0;
}

PetscErrorCode ibmv_cent_of_mass(IBMVNodes *ibmv, FSInfo *FSinfo,PetscInt ibi)

{

  PetscInt	n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ,*nv4;
  PetscInt	i,j,rank;
  PetscInt	n1e, n2e, n3e, n4e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13, dx14, dy14, dz14;
  PetscReal     dr,dV,V_tot=0.0;
  PetscReal     x_c=0.0,y_c=0.0,z_c=0.0; 

  PetscReal     det,J[3][3];
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  n_elmt= ibmv->n_elmt;
  nv1=ibmv->nv1; nv2=ibmv->nv2 ; nv3=ibmv->nv3 ; nv4=ibmv->nv4 ;
  x_bp=ibmv->x_bp;  y_bp=ibmv->y_bp ; z_bp = ibmv->z_bp ;
  
  
  PetscPrintf(PETSC_COMM_WORLD, "number of volume element  is %d \n", ibmv->n_elmt);
  i=0;
  PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d %d\n", nv1[i], nv2[i], nv3[i], nv4[i]);
  
    
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
    
    dV = fabs(dr/6.); 
    ibmv->dV0[i]=dV;
    V_tot+=dV;
    
    ibmv->cent_x[i]= (x_bp[n1e]+x_bp[n2e]+x_bp[n3e]+x_bp[n4e])/4.;
    ibmv->cent_y[i]= (y_bp[n1e]+y_bp[n2e]+y_bp[n3e]+y_bp[n4e])/4.;
    ibmv->cent_z[i]= (z_bp[n1e]+z_bp[n2e]+z_bp[n3e]+z_bp[n4e])/4.;
    
    x_c +=ibmv->cent_x[i]*ibmv->dV0[i];
    y_c +=ibmv->cent_y[i]*ibmv->dV0[i];
    z_c +=ibmv->cent_z[i]*ibmv->dV0[i];
    
  }
  
  x_c /=V_tot;
  y_c /=V_tot;
  z_c /=V_tot;

  FSinfo->x_c=x_c;
  FSinfo->y_c=y_c;
  FSinfo->z_c=z_c;

  FSinfo->a_c[0]=x_c;
  FSinfo->a_c[1]=y_c;
  FSinfo->a_c[2]=z_c;

  FSinfo->mu_s=V_tot;

  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      J[i][j]=0.0;
    }
  }
  
  for (i=0; i<n_elmt; i++) {
    J[0][0] +=((ibmv->cent_y[i]-y_c)*(ibmv->cent_y[i]-y_c)+(ibmv->cent_z[i]-z_c)*(ibmv->cent_z[i]-z_c))*ibmv->dV0[i];
    J[1][1] +=((ibmv->cent_x[i]-x_c)*(ibmv->cent_x[i]-x_c)+(ibmv->cent_z[i]-z_c)*(ibmv->cent_z[i]-z_c))*ibmv->dV0[i];
    J[2][2] +=((ibmv->cent_y[i]-y_c)*(ibmv->cent_y[i]-y_c)+(ibmv->cent_x[i]-x_c)*(ibmv->cent_x[i]-x_c))*ibmv->dV0[i];
    J[0][1] -=(ibmv->cent_x[i]-x_c)*(ibmv->cent_y[i]-y_c)*ibmv->dV0[i];
    J[0][2] -=(ibmv->cent_x[i]-x_c)*(ibmv->cent_z[i]-z_c)*ibmv->dV0[i];
    J[1][2] -=(ibmv->cent_z[i]-z_c)*(ibmv->cent_y[i]-y_c)*ibmv->dV0[i];
  }
  J[1][0]=J[0][1];
  J[2][0]=J[0][2];
  J[2][1]=J[1][2];
  
  ibmv->V=V_tot;

  ibmv->x_c=x_c;
  ibmv->y_c=y_c;
  ibmv->z_c=z_c;

  ibmv->J[0][0]=J[0][0];
  ibmv->J[0][1]=J[0][1];
  ibmv->J[0][2]=J[0][2];
  ibmv->J[1][0]=J[1][0];
  ibmv->J[1][1]=J[1][1];
  ibmv->J[1][2]=J[1][2];
  ibmv->J[2][0]=J[2][0];
  ibmv->J[2][1]=J[2][1];
  ibmv->J[2][2]=J[2][2];
 
  det=J[0][0]*J[1][1]*J[2][2]+2*J[0][1]*J[0][2]*J[1][2]-J[0][0]*J[1][2]*J[1][2]-J[1][1]*J[0][2]*J[0][2]-J[2][2]*J[0][1]*J[0][1];
 
  ibmv->I_inv[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[1][2])/det;
  ibmv->I_inv[0][1]=(J[0][2]*J[1][2]-J[0][1]*J[2][2])/det;
  ibmv->I_inv[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/det;
  ibmv->I_inv[1][0]=ibmv->I_inv[0][1];
  ibmv->I_inv[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[0][2])/det;
  ibmv->I_inv[1][2]=(J[0][1]*J[0][2]-J[0][0]*J[1][2])/det;
  ibmv->I_inv[2][0]=ibmv->I_inv[0][2];
  ibmv->I_inv[2][1]=ibmv->I_inv[1][2];
  ibmv->I_inv[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[0][1])/det;

  for (i=0;i<3;i++){
    for (j=0;j<3;j++){

      FSinfo->I_inv[i][j]=ibmv->I_inv[i][j]; 

    }
  }


  PetscPrintf(PETSC_COMM_WORLD, "Total Volume is  %le \n",ibmv->V);
  PetscPrintf(PETSC_COMM_WORLD, "Centre of mass: x_c  %le y_c %le z_c %le \n",ibmv->x_c,ibmv->y_c,ibmv->z_c);
  PetscPrintf(PETSC_COMM_WORLD, "Moment of Inertia: I_xx  %le I_xy %le I_xz %le I_yy %le I_yz %le I_zz %le \n",ibmv->J[0][0],ibmv->J[0][1],ibmv->J[0][2],ibmv->J[1][1],ibmv->J[1][2],ibmv->J[2][2]);
  PetscPrintf(PETSC_COMM_WORLD, "Moment of Inertia: ibmv->I_inv_xx  %le ibmv->I_inv_xy %le ibmv->I_inv_xz %le ibmv->I_inv_yy %le ibmv->I_inv_yz %le ibmv->I_inv_zz %le \n",ibmv->I_inv[0][0],ibmv->I_inv[0][1],ibmv->I_inv[0][2],ibmv->I_inv[1][1],ibmv->I_inv[1][2],ibmv->I_inv[2][2]);

  return 0;

}
PetscErrorCode calc_quarternion(FSInfo *FSinfo,PetscReal dt,PetscInt ibi)
{
  PetscInt     i,j,k;
 
  PetscReal    w_r[3],w_o[3],w_n[3];
  PetscReal    q_r[4],q_o[4],q_n[4];
  PetscReal    M_x,M_y,M_z; //Torque
  PetscReal    L_n[3],L_r[3],L_o[3]; //Angular momentum
  PetscReal    residue=0.0,w;
  PetscReal    dtaw=0.5*dt;
  PetscReal    rhs[4],R[3][3],RT[3][3],X[3][3],I_n[3][3];

  for (i=0;i<3;i++){
    w_r[i]=0.0;
    w_o[i]=0.0;
    w_n[i]=0.0;
    L_n[i]=0.0;
    L_o[i]=0.0;
    L_r[i]=0.0;
  }
  for (i=0;i<4;i++){
    q_r[i]=0.0;
    q_o[i]=0.0;
    q_n[i]=0.0;
  }

  w_r[0]=FSinfo->S_ang_r[1];
  w_r[1]=FSinfo->S_ang_r[3];
  w_r[2]=FSinfo->S_ang_r[5];

  w_o[0]=FSinfo->S_ang_r[1];
  w_o[1]=FSinfo->S_ang_r[3];
  w_o[2]=FSinfo->S_ang_r[5];

  for (i=0;i<4;i++){
    q_r[i]=FSinfo->q_r[i];
    q_o[i]=FSinfo->q_r[i];
  }
 // Trapezoidal rule
 /*  if (fabs(FSinfo->M_x_real)<1.0e-10) FSinfo->M_x_real= FSinfo->M_x; */
/*   if (fabs(FSinfo->M_y_real)<1.0e-10) FSinfo->M_y_real= FSinfo->M_y; */
/*   if (fabs(FSinfo->M_z_real)<1.0e-10) FSinfo->M_z_real= FSinfo->M_z; */


/*   PetscPrintf(PETSC_COMM_WORLD, "M_x %le M_y %le M_z %le \n",FSinfo->M_x,FSinfo->M_y,FSinfo->M_z); */
 
/*   PetscPrintf(PETSC_COMM_WORLD, "M_x_real %le M_y_real %le M_z_real %le \n",FSinfo->M_x_real,FSinfo->M_y_real,FSinfo->M_z_real); */

  M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_real);
  M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real);
  M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real);

/*   PetscPrintf(PETSC_COMM_WORLD, "q[0] %le q[1] %le q[2] %le q[3] %le \n",FSinfo->q[0],FSinfo->q[1],FSinfo->q[2],FSinfo->q[3]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "M_x %le M_y %le M_z %le \n",M_x,M_y,M_z); */
/*   PetscPrintf(PETSC_COMM_WORLD, "w_x_r %le w_y_r %le w_z_r %le \n",FSinfo->S_ang_r[1],FSinfo->S_ang_r[3],FSinfo->S_ang_r[5]); */

  for (i=0;i<3;i++){
    L_r[i] = FSinfo->L_r[i];
  }
  
  //  PetscPrintf(PETSC_COMM_WORLD, "L_r[0] %le L_r[1] %le L_r[2] %le \n",L_r[0],L_r[1],L_r[2]);
 
  
  L_n[0]=dt*M_x+L_r[0];
  L_n[1]=dt*M_y+L_r[1];
  L_n[2]=dt*M_z+L_r[2];

  // Relaxation
  if (STRONG_COUPLING) {
   /*  if (fabs(FSinfo->M_x_old)<1.0e-10) FSinfo->M_x_old= FSinfo->M_x; */
/*     if (fabs(FSinfo->M_y_old)<1.0e-10) FSinfo->M_y_old= FSinfo->M_y; */
/*     if (fabs(FSinfo->M_z_old)<1.0e-10) FSinfo->M_z_old= FSinfo->M_z; */

    M_x = 0.5*((0.5*FSinfo->M_x+0.5*FSinfo->M_x_old) + FSinfo->M_x_real);
    M_y = 0.5*((0.5*FSinfo->M_y+0.5*FSinfo->M_y_old) + FSinfo->M_y_real);
    M_z = 0.5*((0.5*FSinfo->M_z+0.5*FSinfo->M_z_old) + FSinfo->M_z_real);

    //  PetscPrintf(PETSC_COMM_WORLD, "Stroung_Coupling works for quaternion \n");

    for (i=0;i<3;i++){
      L_r[i] = FSinfo->L_r[i];
      L_o[i] = FSinfo->L_o[i];
    }

    FSinfo->atk=0.;
    for (i=1;i<3;i++){
      FSinfo->dS[i]=L_o[i]-L_n[i];
    
      if (fabs(FSinfo->dS[i]-FSinfo->dS_o[i])>1e-8  &&  FSinfo->atk_o!=0.9) {
	FSinfo->atk+=(FSinfo->dS[i])/(FSinfo->dS_o[i]-FSinfo->dS[i]);
      }
    }
    FSinfo->atk=FSinfo->atk_o+(FSinfo->atk_o-1)*FSinfo->atk;
    if (FSinfo->atk>.9) FSinfo->atk=.9;
    if (FSinfo->atk<0.1) FSinfo->atk=0.1;
    w=1.-FSinfo->atk;
    PetscOptionsGetReal(PETSC_NULL, "-w_str", &w, PETSC_NULL);
    //   w=0.1;
    PetscPrintf(PETSC_COMM_WORLD, "under relaxation coefficient %le \n",w);
 
    for (i=0;i<3;i++){
      L_n[i]=w*L_n[i]+(1.-w)*L_o[i];
    }
  
  }
  //

  //  PetscPrintf(PETSC_COMM_WORLD, "L_x %le L_y  %le L_z %le  \n",L_n[0],L_n[1],L_n[2]);

  for (k=0;k<50;k++){ /*  while (residue>0.001){ */
    rhs[0]=0.25*(-w_o[0]*q_o[1]-w_o[1]*q_o[2]-w_o[2]*q_o[3])+0.25*(-w_r[0]*q_r[1]-w_r[1]*q_r[2]-w_r[2]*q_r[3]);
    rhs[1]=0.25*( w_o[0]*q_o[0]+w_o[1]*q_o[3]-w_o[2]*q_o[2])+0.25*( w_r[0]*q_r[0]+w_r[1]*q_r[3]-w_r[2]*q_r[2]);
    rhs[2]=0.25*(-w_o[0]*q_o[3]+w_o[1]*q_o[0]+w_o[2]*q_o[1])+0.25*(-w_r[0]*q_r[3]+w_r[1]*q_r[0]+w_r[2]*q_r[1]);
    rhs[3]=0.25*( w_o[0]*q_o[2]-w_o[1]*q_o[1]+w_o[2]*q_o[0])+0.25*( w_r[0]*q_r[2]-w_r[1]*q_r[1]+w_r[2]*q_r[0]);
    
    for (i=0;i<4;i++){
      q_n[i]=q_o[i]+dtaw*(rhs[i]+(q_r[i]-q_o[i])/dt);
    }
    //  PetscPrintf(PETSC_COMM_WORLD, "q[0] %le q[1] %le q[2] %le q[3] %le \n",q_n[0],q_n[1],q_n[2],q_n[3]);
    R[0][0]=1.0-2.0*q_n[2]*q_n[2]-2.0*q_n[3]*q_n[3];
    R[0][1]=2.0*(q_n[1]*q_n[2]-q_n[0]*q_n[3]);
    R[0][2]=2.0*(q_n[1]*q_n[3]+q_n[0]*q_n[2]);
    R[1][0]=2.0*(q_n[1]*q_n[2]+q_n[0]*q_n[3]);
    R[1][1]=1.0-2.0*q_n[1]*q_n[1]-2.0*q_n[3]*q_n[3];
    R[1][2]=2.0*(q_n[2]*q_n[3]-q_n[0]*q_n[1]);
    R[2][0]=2.0*(q_n[1]*q_n[3]-q_n[0]*q_n[2]);
    R[2][1]=2.0*(q_n[2]*q_n[3]+q_n[0]*q_n[1]);
    R[2][2]=1.0-2.0*q_n[1]*q_n[1]-2.0*q_n[2]*q_n[2];
    //  PetscPrintf(PETSC_COMM_WORLD, "R: R[0][0] %le R[0][1] %le R[0][2] %le R[1][1] %le R[1][2] %le R[2][2] %le \n",R[0][0],R[0][1],R[0][2],R[1][1],R[1][2],R[2][2]);

    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	if (i==j) RT[i][j]=R[i][j];
	else      RT[i][j]=R[j][i];
      }
    }
    // PetscPrintf(PETSC_COMM_WORLD, "I_inv: I_n[0][0] %le I_inv[0][1] %le I_inv[0][2] %le I_inv[1][1] %le I_inv[1][2] %le I_inv[2][2] %le \n",FSinfo->I_inv[0][0],FSinfo->I_inv[0][1],FSinfo->I_inv[0][2],FSinfo->I_inv[1][1],FSinfo->I_inv[1][2],FSinfo->I_inv[2][2]);
    //  Mat_Maltiply(X,R,RT);
    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	X[i][j]=FSinfo->I_inv[i][0]*RT[0][j]+FSinfo->I_inv[i][1]*RT[1][j]+FSinfo->I_inv[i][2]*RT[2][j];
      }
    }
    //  Mat_Maltiply(X,FSinfo->I_inv,RT);

    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	I_n[i][j]=R[i][0]*X[0][j]+R[i][1]*X[1][j]+R[i][2]*X[2][j];
      }
    }
    //    Mat_Maltiply(I_n,R,X);
    //  PetscPrintf(PETSC_COMM_WORLD, "I_n: I_n[0][0] %le I_n[0][1] %le I_n[0][2] %le I_n[1][1] %le I_n[1][2] %le I_n[2][2] %le \n",I_n[0][0],I_n[0][1],I_n[0][2],I_n[1][1],I_n[1][2],I_n[2][2]);
    w_n[0]=I_n[0][0]*L_n[0]+I_n[0][1]*L_n[1]+I_n[0][2]*L_n[2];
    w_n[1]=I_n[1][0]*L_n[0]+I_n[1][1]*L_n[1]+I_n[1][2]*L_n[2];
    w_n[2]=I_n[2][0]*L_n[0]+I_n[2][1]*L_n[1]+I_n[2][2]*L_n[2];

    residue=0.0;

    for (i=0;i<4;i++){
      residue +=(q_n[i]-q_o[i])*(q_n[i]-q_o[i]);
    }

    // PetscPrintf(PETSC_COMM_WORLD, "w_n[0] %le w_n[1]  %le w_n[2] %le  \n",w_n[0],w_n[1],w_n[2]);
    residue=sqrt(residue/4);
    //  PetscPrintf(PETSC_COMM_WORLD, "residue is  %le  \n",residue);
    for (i=0;i<4;i++){
      q_o[i]=q_n[i];
    }
    for (i=0;i<3;i++){
      w_o[i]=w_n[i];
    }
  }

  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      FSinfo->R[i][j]=R[i][j];
    }
  }
 
  //  PetscPrintf(PETSC_COMM_WORLD, "Rotation Matrix: R[0][0] %le R[0][1] %le R[0][2] %le R[1][1] %le R[1][2] %le R[2][2] %le \n",FSinfo->R[0][0],FSinfo->R[0][1],FSinfo->R[0][2],FSinfo->R[1][1],FSinfo->R[1][2],FSinfo->R[2][2]);

  // PetscPrintf(PETSC_COMM_WORLD, "Inverse of Moment of Inertia: I_n[0][0]  %le I_n[0][1] %le I_n[0][2] %le I_n[1][1] %le I_n[1][2] %le I_n[2][2] %le \n",I_n[0][0],I_n[0][1],I_n[0][2],I_n[1][1],I_n[1][2],I_n[2][2]);


  for (i=0;i<4;i++){
    FSinfo->q[i]=q_n[i];
  }
  
 /*  FSinfo->S_ang_n[1]=w_r[0]; */
/*   FSinfo->S_ang_n[3]=w_r[1]; */
/*   FSinfo->S_ang_n[5]=w_r[2]; */

  FSinfo->S_ang_n[1]=w_n[0];
  FSinfo->S_ang_n[3]=w_n[1];
  FSinfo->S_ang_n[5]=w_n[2];

  // FSinfo->S_ang_n[0]=FSinfo->S_ang_r[0]+0.5*(FSinfo->S_ang_n[1]+FSinfo->S_ang_r[1])*dt;
  
  FSinfo->alpha[0]=(FSinfo->S_ang_n[1]-FSinfo->S_ang_r[1])/dt;
  FSinfo->alpha[1]=(FSinfo->S_ang_n[3]-FSinfo->S_ang_r[3])/dt;
  FSinfo->alpha[2]=(FSinfo->S_ang_n[5]-FSinfo->S_ang_r[5])/dt;

  for (i=0;i<3;i++){
    FSinfo->L_n[i]=L_n[i];
  }
  
  // output values
  PetscPrintf(PETSC_COMM_WORLD, "w_x %le w_y  %le w_z %le  \n",FSinfo->S_ang_n[1],FSinfo->S_ang_n[3],FSinfo->S_ang_n[5]);
  PetscPrintf(PETSC_COMM_WORLD, "a_x %le a_y  %le a_z %le  \n",FSinfo->alpha[0],FSinfo->alpha[1],FSinfo->alpha[2]);
  PetscPrintf(PETSC_COMM_WORLD, "q[0] %le q[1] %le q[2] %le q[3] %le \n",FSinfo->q[0],FSinfo->q[1],FSinfo->q[2],FSinfo->q[3]);
  return 0;
}

/* PetscErrorCode Mat_Maltiply(PetscReal C[3][3],PetscReal A[3][3],PetscReal B[3][3]) */
/* { */
/*   PetscInt i,j; */
/*   for (i=0;i<3;i++){ */
/*     for (j=0;j<3;j++){ */
/*       C[i][j]=A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j]; */
/*     } */
/*   } */
/*   return 0; */
/* } */
PetscErrorCode rotate_quarternion(IBMNodes *ibm, FSInfo *fsi, PetscReal dt)
{
  PetscInt n_v = ibm->n_v;
  PetscInt i,j;
  Cmpnts   p,a_c;
  PetscReal Rot[3][3];

 
  for (i=1;i<6;i+=2){
    fsi->S_new[i-1]=fsi->S_real[i-1]+0.5*(fsi->S_new[i]+fsi->S_real[i])*dt;
  }
  
  
  fsi->x_c=fsi->a_c[0]+fsi->S_new[0];
  fsi->y_c=fsi->a_c[1]+fsi->S_new[2];
  fsi->z_c=fsi->a_c[2]+fsi->S_new[4];
  
  
  a_c.x=fsi->a_c[0];
  a_c.y=fsi->a_c[1];
  a_c.z=fsi->a_c[2];

  PetscPrintf(PETSC_COMM_WORLD, " initial center of mass: a_c.x %le a_c.y %le a_c.z %le \n",a_c.x,a_c.y,a_c.z);
 
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      Rot[i][j]=fsi->R[i][j];
      //  PetscPrintf(PETSC_COMM_WORLD, "R[%d][%d]= %le \n",i,j,Rot[i][j]);
    }
  }
 
  for (i=0; i<n_v; i++) {
   
    p.x = ibm->x_bp0[i]-fsi->a_c[0];
    p.y = ibm->y_bp0[i]-fsi->a_c[1];
    p.z = ibm->z_bp0[i]-fsi->a_c[2];
  
    ibm->x_bp[i] =Rot[0][0]*p.x+Rot[0][1]*p.y+Rot[0][2]*p.z+fsi->x_c;
    ibm->y_bp[i] =Rot[1][0]*p.x+Rot[1][1]*p.y+Rot[1][2]*p.z+fsi->y_c;
    ibm->z_bp[i] =Rot[2][0]*p.x+Rot[2][1]*p.y+Rot[2][2]*p.z+fsi->z_c;
  }


  //  PetscPrintf(PETSC_COMM_WORLD, "fsi[%d].S_new[0] %le \n",NumberOfBodies+particle_added-1,fsi[NumberOfBodies+particle_added-1].S_new[0]);


  return(0);
}

PetscErrorCode Particle_Duplicate_Indentify (UserCtx *user,IBMNodes *ibm, FSInfo *fsi, PetscInt *X, PetscInt *Z)
{
  PetscReal     L_x,L_y,L_z; 
  //  PetscReal	*x_bp =ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscInt      ibi,i,n_v,z_mode=0,x_mode=0;

 
  for (ibi=0; ibi< NumberOfBodies; ibi++){
    
    xbp_min = 1.e23; xbp_max = -1.e23;
    ybp_min = 1.e23; ybp_max = -1.e23;
    zbp_min = 1.e23; zbp_max = -1.e23; 

    n_v = ibm[ibi].n_v;
   
    for(i=0; i<n_v; i++) {
    
      xbp_min = PetscMin(xbp_min, ibm[ibi].x_bp[i]);
      xbp_max = PetscMax(xbp_max, ibm[ibi].x_bp[i]);
    
      ybp_min = PetscMin(ybp_min, ibm[ibi].y_bp[i]);
      ybp_max = PetscMax(ybp_max, ibm[ibi].y_bp[i]);
    
      zbp_min = PetscMin(zbp_min, ibm[ibi].z_bp[i]);
      zbp_max = PetscMax(zbp_max, ibm[ibi].z_bp[i]);
    }
  
    ibm[ibi].x_min = xbp_min-0.1; ibm[ibi].x_max = xbp_max+0.1;
    ibm[ibi].y_min = ybp_min-0.1; ibm[ibi].y_max = ybp_max+0.1;
    ibm[ibi].z_min = zbp_min-0.1; ibm[ibi].z_max = zbp_max+0.1;

    PetscOptionsGetReal(PETSC_NULL, "-L_x", &L_x, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-L_y", &L_y, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-L_z", &L_z, PETSC_NULL);
    
    if (ibm[ibi].z_max < L_z && ibm[ibi].z_min > 0.0 ) z_mode=0;
    else if (ibm[ibi].z_min < 0.0 && ibm[ibi].z_max > 0.0) z_mode=1;
    else if (ibm[ibi].z_min < L_z && ibm[ibi].z_max > L_z) z_mode=2;
    else z_mode=3;

  
    if ( ibm[ibi].x_max < L_x && ibm[ibi].x_min > 0.0 ) x_mode=0;
    else if(ibm[ibi].x_min < 0.0 && ibm[ibi].x_max > 0.0) x_mode=1;
    else if(ibm[ibi].x_min < L_x && ibm[ibi].x_max > L_x) x_mode=2;
    else x_mode=3;
   
    if (z_mode==1 && fsi[ibi].pbc[2]==0) Z[ibi]=1;
    if (z_mode==2 && fsi[ibi].pbc[2]==0) Z[ibi]=2;

    if (x_mode==1 && fsi[ibi].pbc[0]==0) X[ibi]=1;
    if (x_mode==2 && fsi[ibi].pbc[0]==0) X[ibi]=2;

    if (z_mode==3 && (fsi[ibi].pbc[2]==1 || fsi[ibi].pbc[2]==2)) Z[ibi]=3;
    if (x_mode==3 && (fsi[ibi].pbc[0]==1 || fsi[ibi].pbc[0]==2)) X[ibi]=3;

    fsi[ibi].pbc[2]=z_mode;
    fsi[ibi].pbc[0]=x_mode;

    PetscPrintf(PETSC_COMM_WORLD,"pbc[0] %d pbcs[1] %d pbc[2] %d  \n",fsi[ibi].pbc[0],fsi[ibi].pbc[1],fsi[ibi].pbc[2]);
  }
 

  
  return(0);
}


PetscErrorCode Particle_Clone (IBMNodes *ibm1, FSInfo *fsi1,IBMNodes *ibm2, FSInfo *fsi2,PetscInt X, PetscInt Z, PetscReal dt,PetscInt ibi)

{
 
  Cmpnts     omega_c, a_c;
  PetscReal  rx,ry,rz,l_x,l_y,l_z;
  PetscInt i,j;

  PetscOptionsGetReal(PETSC_NULL, "-L_x", &l_x, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_y", &l_y, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-L_z", &l_z, PETSC_NULL);       

  if (Z==2){
    ibm_duplicate(ibm1,ibm2, 0.0, 0.0,-l_z);
      
    fsi1->a_c[0]= fsi2->a_c[0];
    fsi1->a_c[1]= fsi2->a_c[1];
    fsi1->a_c[2]= fsi2->a_c[2]-l_z;
  }
  else if (Z==1){
    ibm_duplicate(ibm1,ibm2, 0.0, 0.0,l_z);
         
    fsi1->a_c[0]= fsi2->a_c[0];
    fsi1->a_c[1]= fsi2->a_c[1];
    fsi1->a_c[2]= fsi2->a_c[2]+l_z;
  }
      
  else if (X==2){
	ibm_duplicate(ibm1,ibm2, 0.0, 0.0,-l_x);
	
	fsi1->a_c[0]= fsi2->a_c[0]-l_x;
	fsi1->a_c[1]= fsi2->a_c[1];
	fsi1->a_c[2]= fsi2->a_c[2];
  }
  else if (X==1){
    ibm_duplicate(ibm1,ibm2, 0.0, 0.0,l_x);
    
    fsi1->a_c[0]= fsi2->a_c[0]+l_x;
    fsi1->a_c[1]= fsi2->a_c[1];
    fsi1->a_c[2]= fsi2->a_c[2];
	
  }
  else   PetscPrintf(PETSC_COMM_WORLD, "Particle duplication error!!!!!!!!!!!\n");
      
  for (i=0;i<3;i++){
    fsi1->pbc[i]=3;
    fsi1->S_ang_n[2*i+1]= fsi2->S_ang_n[2*i+1];
    fsi1->alpha[i]=fsi2->alpha[i];
    fsi1->acc[i]=fsi2->acc[i];
    fsi1->S_new[2*i+1]=fsi2->S_new[2*i+1];
    fsi1->S_real[2*i]=fsi2->S_real[2*i];
    fsi1->S_real[2*i+1]=fsi2->S_real[2*i+1];
    fsi1->L_n[i]=fsi2->L_n[i];
  }

  fsi1->mu_s=fsi1->mu_s;
  fsi1->x_c=fsi1->a_c[0]+fsi1->S_new[0];
  fsi1->y_c=fsi1->a_c[1]+fsi1->S_new[2];
  fsi1->z_c=fsi1->a_c[2]+fsi1->S_new[4];
	 
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      fsi1->R[i][j]=fsi2->R[i][j];
    }
  }
  
  rotate_quarternion(ibm1,fsi1,dt);
  
  calc_ibm_normal(ibm1);
  
  omega_c.x=fsi1->S_ang_n[1];
  omega_c.y=fsi1->S_ang_n[3];
  omega_c.z=fsi1->S_ang_n[5];
  
  a_c.x=fsi1->x_c;
  a_c.y=fsi1->y_c;
  a_c.z=fsi1->z_c;
  
  for (i=0; i<ibm1->n_v; i++) {
    
    rx = ibm1->x_bp[i]-a_c.x;
    ry = ibm1->y_bp[i]-a_c.y;
    rz = ibm1->z_bp[i]-a_c.z;
    ibm1->u[i].x =   (rz*omega_c.y-omega_c.z*ry)+fsi1->S_new[1]  ;
    ibm1->u[i].y =   (rx*omega_c.z-omega_c.x*rz)+fsi1->S_new[3]  ;
    ibm1->u[i].z =   (ry*omega_c.x-omega_c.y*rx)+fsi1->S_new[5]  ;
    
  }

  if (fsi2->clone<1.e-6){
    fsi2->clone=0.001;
    fsi1->clone=1+(ibi+1.)/1000;
  }else if (fsi2->clone < 1.0){

    fsi2->clone=0.003;
    fsi1->clone=2+(ibi+1.)/1000;
  }else fsi1->clone=3+(ibi+1.)/1000;

 
  return(0);
}
