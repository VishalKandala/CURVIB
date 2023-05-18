#include "variables.h"

extern PetscInt block_number, NumberOfBodies,blank,tiout;


PetscErrorCode platelet_activation(UserCtx *user,IBMNodes *ibm,PetscInt ti,PetscBool rstart_flg,PetscInt immersed){

  DM             da = user->da, fda = user->fda;
 Cmpnts               ***ucont,***ucontl;

 DMDALocalInfo       info = user->info;
 PetscInt        xs = info.xs, xe = info.xs + info.xm;
 PetscInt         ys = info.ys, ye = info.ys + info.ym;
 PetscInt        zs = info.zs, ze = info.zs + info.zm;
 PetscInt        mx = info.mx, my = info.my, mz = info.mz;
 PetscMPIInt        size,rank;

 Cmpnts            ***icsi, ***jeta, ***kzet,***coor;

 PetscReal         dt=user->dt,RE=user->ren;
 PetscReal         eta=RE*0.0;
 PetscInt          nx=user->IM,ny=user->JM,nz=user->KM;
 PetscInt          i,j,k,itter;
 PetscInt        lxs, lxe, lys, lye, lzs, lze;
 PetscReal         ep=.05*dt;


 PetscReal         ***laj,***Fe,***Fee;
 PetscReal        ***MArr,***MArrnew,***shrrarr,***nvert;
 Mat A;
 KSP ksp;
 Vec RHS;
 struct Cmpntsgdb  ***cent;

 PetscInt ibi;
 PetscInt bi=user->_this;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

///////////////////////////////////////////////////////////
MPI_Comm_size(PETSC_COMM_WORLD,&size);
MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

DMDAVecGetArray(fda, user->lUcont, &ucont);
DMDAVecGetArray(da, user->lAj, &laj);
DMDAVecGetArray(da, user->shrr, &shrrarr);
DMDAVecGetArray(fda, user->Cent, &cent);
DMDAVecGetArray(da, user->lNvert, &nvert);
//////////////////////////////////////////////////////////Initialize

 if(ti==0){
VecSet(user->Qnew,0.0);
/* DMDAVecGetArray(da,user->Qnew,&MArr); */
/*  if(user->_this==0){ */
/*  for(k=zs;k<ze;k++){ */
/*    for(j=ys;j<ye;j++){ */
/*      for(i=xs;i<xe;i++){ */
/*        if(k<10) */
/*      MArr[k][j][i]=1.0; */
   
/*      } */
/*    } */
/*  } */
/*  } */
/*  DMDAVecRestoreArray(da,user->Qnew,&MArr); */
 }
 ///////////////////

 if(immersed&&user->_this==0){
   for (ibi=0;ibi<NumberOfBodies;ibi++) {
      ibm_interpolation_Q(&user[0], &ibm[ibi], ibi,0);
   }
 }

 ///////////////////
 DMGlobalToLocalBegin(da,user->Qnew,INSERT_VALUES,user->Ql);
 DMGlobalToLocalEnd(da,user->Qnew,INSERT_VALUES,user->Ql);
 //////////////////
VecDuplicate(user->Qnew,&RHS);
VecSet(RHS,0.0);
DMDAVecGetArray(da,user->Ql,&MArr);
DMDAVecGetArray(da,RHS,&MArrnew);

 ///////////////////////////////////////////////////////////////RHS
  for(k=zs;k<ze;k++){
   for(j=ys;j<ye;j++){
     for(i=xs;i<xe;i++){
       if(i==0){
	 MArrnew[k][j][i]=0.0;
       }else if(i==mx-1){
	 MArrnew[k][j][i]=0.0;
       }else if(j==0){
	 MArrnew[k][j][i]=0.0;
       }else if(j==my-1){
	 MArrnew[k][j][i]=0.0;
       }else if(k==0){
	 if(user->_this==1){
	   MArrnew[k][j][i]=MArr[k][j][i];
	 }
       }else if(k==mz-1){
	 MArrnew[k][j][i]=0.0;
       }else{
	 if(nvert[k][j][i]==0){
	   MArrnew[k][j][i]= MArr[k][j][i]-0.5*(dt*laj[k][j][i])*(0.5*(MArr[k][j][i+1]+MArr[k][j][i])*ucont[k][j][i].x-0.5*(MArr[k][j][i]+MArr[k][j][i-1])*ucont[k][j][i-1].x+0.5*(MArr[k][j+1][i]+MArr[k][j][i])*ucont[k][j][i].y-0.5*(MArr[k][j][i]+MArr[k][j-1][i])*ucont[k][j-1][i].y+0.5*(MArr[k+1][j][i]+MArr[k][j][i])*ucont[k][j][i].z-0.5*(MArr[k][j][i]+MArr[k-1][j][i])*ucont[k-1][j][i].z)+.01*(dt*laj[k][j][i])*(shrrarr[k][j][i]);
	 }else if(nvert[k][j][i]==1){
	   MArrnew[k][j][i]= MArr[k][j][i];
	 }else{
	   MArrnew[k][j][i]=0.0;
	 }
       }
     }
   }
 }
 DMDAVecRestoreArray(da,RHS,&MArrnew);
 DMDAVecRestoreArray(da,user->Ql,&MArr);

 ///////////////////////////////////////////////////////////////

/*  //Create Matrix */
/* DMCreateMatrix(da,MATAIJ,&A); */
/* MatZeroEntries(A); */

  //Mat
  PetscInt M;
  PetscInt N;
  VecGetLocalSize(user->Qnew,&M);
  VecGetSize(user->Qnew,&N);
  MatCreateAIJ(MPI_COMM_WORLD,M,M,N,N,30,PETSC_NULL,30,PETSC_NULL,&A);
  MatZeroEntries(A);


 //Load Mtrix

 PetscReal val=0.0;
 PetscInt row,col;
 PetscReal am1,a,ap1,bm1,b,bp1, cm1,c,cp1;

 for(k=lzs;k<lze;k++)
   for(j=lys;j<lye;j++)
     for(i=lxs;i<lxe;i++)
       {
	 if(nvert[k][j][i]==0){

     am1=1.0,a=-2.0,ap1=1.0,bm1=1.0,b=-2.0,bp1=1.0,cm1=1.0,c=-2.0,cp1=1.0;

     row = Gidx(i, j, k, user);
     col=row;
     val=(1.0+0.5*(0.5*ucont[k][j][i].x-0.5*ucont[k][j][i-1].x+0.5*ucont[k][j][i].y-0.5*ucont[k][j-1][i].y+0.5*ucont[k][j][i].z-0.5*ucont[k-1][j][i].z)*(dt*laj[k][j][i])+ep*eta-ep*(a+b+c));
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);

     col=Gidx(i+1, j, k, user);
     val=(0.5*(0.5*ucont[k][j][i].x)*(dt*laj[k][j][i])-eta-ep*ap1);
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);

     col=Gidx(i-1, j, k, user);
     val=(0.5*(-0.5*ucont[k][j][i-1].x)*(dt*laj[k][j][i])-eta-ep*am1);
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);

     col=Gidx(i, j+1, k, user);
     val=(0.5*(0.5*ucont[k][j][i].y)*(dt*laj[k][j][i])-eta-ep*bp1);
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);

     col=Gidx(i, j-1, k, user);
     val=(0.5*(-0.5*ucont[k][j-1][i].y)*(dt*laj[k][j][i])-eta-ep*bm1);
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);

     col=Gidx(i, j, k+1, user);
     val=(0.5*(0.5*ucont[k][j][i].z)*(dt*laj[k][j][i])-eta-ep*cp1);
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);

     col=Gidx(i, j, k-1, user);
     val=(0.5*(-0.5*ucont[k-1][j][i].z)*(dt*laj[k][j][i])-eta-ep*cm1);
     MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
	 }else{
	  row = Gidx(i, j, k, user);
	  col=row;
	  val=1.0;
	  MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
	 }

       }


    /////////////////////////////////////////////////////////////////boundary condition

 for(k=zs;k<ze;k++)
   for(j=ys;j<ye;j++)
     for(i=xs;i<xe;i++)
       {
     row = Gidx(i, j, k, user);

     if(i==0){
       col=Gidx(i, j, k, user);
       val=1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
     }else if(i==mx-1){
       col=Gidx(i, j, k, user);
       val=1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
     }else if (j==0){
       col=Gidx(i, j, k, user);
       val=1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
     }else if(j==my-1){
       col=Gidx(i, j, k, user);
       val=1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
     }else if (k==0){
       col=Gidx(i, j, k, user);
       val=1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
     }else if(k==mz-1){
       col=Gidx(i, j, k, user);
       val=1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
       //////////////
       col=Gidx(i, j, k-1, user);
       val=-1.0;
       MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);
     }
       }
 ///////////////////////////////////////////////////////////////

 MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
 MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
 KSPCreate(MPI_COMM_WORLD,&ksp);
//End Load Matrix


//Solve KSP
KSPAppendOptionsPrefix(ksp, "pl_");
KSPSetFromOptions(ksp);
KSPSetOperators(ksp,A,A);
KSPSetUp(ksp);
KSPSolve(ksp,RHS,user->Qnew);

 
//////////////////////////////////////
 if (ti == (ti/tiout) * tiout){
 char filen[50];
 PetscViewer viewer;
   sprintf(filen, "Qfield%5.5d_%1.1d.dat", ti, user->_this);
   PetscViewerBinaryOpen(PETSC_COMM_WORLD,filen,FILE_MODE_WRITE,&viewer);
   VecView(user->shrr,viewer);
   PetscViewerDestroy(&viewer);
 }
   ///////////////////////////////////
 DMDAVecRestoreArray(da,user->lAj,&laj);
 DMDAVecRestoreArray(fda,user->lUcont,&ucont);
 DMDAVecRestoreArray(fda,user->Cent,&cent);
 DMDAVecRestoreArray(da,user->shrr,&shrrarr);
 DMDAVecRestoreArray(da, user->lNvert, &nvert);
  //////////////////////////////////Destroy
  VecDestroy(&(RHS));
  MatDestroy(&(A));
  KSPDestroy(&(ksp));

   PetscPrintf(PETSC_COMM_WORLD,"End of the platelet file:%d\n",ti);

return(0);

}



PetscErrorCode shear(UserCtx *user,PetscInt ti){

  DM            da = user->da, fda = user->fda;
  DMDALocalInfo    info = user->info;
  PetscInt    xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt    zs = info.zs, ze = info.zs + info.zm;
  PetscInt    mx = info.mx, my = info.my, mz = info.mz;

  PetscInt    lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k, ii,jj,kk;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;//, ***l2;

  PetscInt   njac, nrot;
  PetscReal  a[3][3], v[3][3], d[3], II[4],x[3];

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;


  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);

  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->shrr, &q);
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
    a[0][0]=1./user->ren*s11;
    a[0][1]=1./user->ren*s12;
    a[0][2]=1./user->ren*s13;
    a[1][0]=1./user->ren*s21;
    a[1][1]=1./user->ren*s22;
    a[1][2]=1./user->ren*s23;
    a[2][0]=1./user->ren*s31;
    a[2][1]=1./user->ren*s32;
    a[2][2]=1./user->ren*s33;

    // invariants
    II[0]=1.;II[1]=0.;II[2]=0.;II[3]=0.;

    smax=sqrt(1.0/3.0)*sqrt(a[0][0]*a[0][0]+a[1][1]*a[1][1]+a[2][2]*a[2][2]-a[0][0]*a[1][1]-a[1][1]*a[2][2]-a[0][0]*a[2][2]+3.0*(a[0][1]*a[0][1]+a[1][2]*a[1][2]+a[0][2]*a[0][2]));

    q[k][j][i] = (smax);

    if(test<fabs(ucat[k][j][i].z)) test=fabs(ucat[k][j][i].z);
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD,"MAX ucat:%f\n",test);

   DMDAVecRestoreArray(fda, user->lUcat, &ucat);
   DMDAVecRestoreArray(fda, user->lCsi, &csi);
   DMDAVecRestoreArray(fda, user->lEta, &eta);
   DMDAVecRestoreArray(fda, user->lZet, &zet);

   DMDAVecRestoreArray(da, user->lAj, &aj);
   DMDAVecRestoreArray(da, user->lNvert, &nvert);
   DMDAVecRestoreArray(da, user->shrr, &q);
     // DMDAVecRestoreArray(da, user->Phi, &l2);
  return 0;
}




PetscErrorCode interpolate(UserCtx *user,PetscInt ti)
{

  DM             da = user->da, fda = user->fda;
  Cmpnts             ***ucont,***ucont1,***ucont12;
  
  DMDALocalInfo       info = user->info;
  PetscInt        xs = info.xs, xe = info.xs + info.xm;
  PetscInt         ys = info.ys, ye = info.ys + info.ym;
  PetscInt        zs = info.zs, ze = info.zs + info.zm;
  PetscInt        mx = info.mx, my = info.my, mz = info.mz;
  PetscMPIInt        size,rank;
  
  PetscInt          nx=user->IM,ny=user->JM,nz=user->KM;
  PetscInt          i,j,k,itter;
  PetscInt        lxs, lxe, lys, lye, lzs, lze;

  struct Cmpntsgdb  ***cent;
  PetscInt STEPS=10;
  PetscInt iter;

  PetscViewer	viewer;
  char filen[80];
  PetscInt N;
  Vec Ucont1,Ucont12;
  ///////////////////////////////////
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ///////////////////////////////////
  VecDuplicate(user->Ucont, &Ucont1);
  VecDuplicate(user->Ucont, &Ucont12);
  ///////////////////////////////////
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, Ucont1, &ucont1);
  DMDAVecGetArray(fda, Ucont12, &ucont12);
  ///////////////////////////////////reading next step
  sprintf(filen, "vfield%5.5d_%1.1d.dat", ti+STEPS-1, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecGetSize(Ucont1, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad((Ucont1),viewer);
  PetscViewerDestroy(&viewer);

  ///////////////////////////////////
  for(iter=0;iter<STEPS-1;iter++){

    DMDAVecGetArray(fda,Ucont12, &ucont12);
    
   for(k=zs;k<ze;k++){
     for(j=ys;j<ye;j++){
       for(i=xs;i<xe;i++){
	 ucont12[k][j][i].x=ucont[k][j][i].x+(iter+1)*((ucont1[k][j][i].x-ucont[k][j][i].x)/(STEPS-1.0));
	 ucont12[k][j][i].y=ucont[k][j][i].y+(iter+1)*((ucont1[k][j][i].y-ucont[k][j][i].y)/(STEPS-1.0));
	 ucont12[k][j][i].z=ucont[k][j][i].z+(iter+1)*((ucont1[k][j][i].z-ucont[k][j][i].z)/(STEPS-1.0));
       }
     }
   }
    DMDAVecRestoreArray(fda,Ucont12,&ucont12);
   ///////////////////////////////////writing interval STEPS
   
   sprintf(filen, "vfield%5.5d_%1.1d.dat", ti+iter, user->_this);
   PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
   VecView(Ucont12, viewer);
   PetscViewerDestroy(&viewer);
  }
   
   ///////////////////////////
   DMDAVecRestoreArray(fda,user->Ucont,&ucont);
   DMDAVecRestoreArray(fda,Ucont1,&ucont1);
   ///////////////////////////
   VecDestroy(&(Ucont1));VecDestroy(&(Ucont12));
   
   return 0;
}



PetscErrorCode ibm_interpolation_Q(UserCtx *user,
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

  PetscInt nbn;
  PetscInt nbnumber = user->ibmnumber;
  PetscInt i,j,k;
  PetscReal sb, sc;
  PetscInt	ip1, ip2, ip3, jp1, jp2, jp3, kp1, kp2, kp3;
  PetscReal sk1, sk2, sk3, cv1, cv2, cv3, phia, phic;
  Cmpnts	***ucat, ***lucat;
  PetscReal	***nvert, ***nvert_o, ***q, ***ql;
  PetscReal cs1, cs2, cs3;
  PetscInt	ni;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  IBMInfo *ibminfo;
  Vec Ql;

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
	
 
  DMDAVecGetArray(da, user->Qnew, &q);

  DMDAVecGetArray(da,user->Ql,&ql);
  DMDAVecGetArray(fda,user->lCent,&coor);
  current = user->ibmlist[ibi].head;
		
    while (current) {
      IBMInfo *ibminfo = &current->ibm_intp;
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
      Cmpnts Ua, Uc;  
    
      PetscReal Ua_n, Ua_nold;
     
      nfx=coor[k][j][i].x-ibminfo->pmin.x;
      nfy=coor[k][j][i].y-ibminfo->pmin.y;
      nfz=coor[k][j][i].z-ibminfo->pmin.z;
      double dr=sqrt(nfx*nfx+nfy*nfy+nfz*nfz);	
      nfx=nfx/dr;
      nfy=nfy/dr;
      nfz=nfz/dr;
 
     
      cv1 = ql[kp1][jp1][ip1];
      cv2 = ql[kp2][jp2][ip2];
      cv3 = ql[kp3][jp3][ip3];
	
      q[k][j][i] = (cv1 * sk1 + cv2 * sk2 + cv3 * sk3);
    }

    DMDAVecRestoreArray(fda,user->lCent,&coor);	

    DMDAVecRestoreArray(da, user->Qnew, &q);
    DMDAVecRestoreArray(da, user->Ql, &ql);
  ////////////////////

  }

  /////////////////////////////
  DMGlobalToLocalBegin(da,user->Qnew,INSERT_VALUES,user->Ql);
  DMGlobalToLocalEnd(da,user->Qnew,INSERT_VALUES,user->Ql);

  return(0);
}
