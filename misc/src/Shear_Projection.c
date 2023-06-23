//Description of the code :

What the code does is to calculate the velocity gradient at all nodes: fluid + IB

Then it projects the velocity gradient on the solid surface using linear interpolation of the velocity gradient tensor.

First, the calculation is done at every processor, thus the local portion of the grid computes the velocity gradient on local machine.
Only the local part of the surface mesh is computed.

Then all the processors send their velocity gradient to the ROOT processor.

The root processor then assemble all the information. It subsequently do averaging on each cell element and assemble it to 
have the full shear stress field on the surface mesh.

It is done in real time.

The ROOT processor is the one who keeps the whole shear stress field. If you want to send the shear stress field to every other processors 
please use MPI function BroadCast. I don't think you might need to do so. 

If you want to change the interpolation scheme other than the solid surface. Change the sb, sc in the IBM_Shear_Stress subroutine.
//-------------------------------


//
To use the code please cooperate the following subroutine
in places you want to calculate the shear stress.

for (bi = 0; bi < block_number; bi++)
{
	     IBM_Shear_Stress(&user[bi], user[bi].ibm->IBMAssemble);
             Shear_Projection(&user[bi], user[bi].ibm->IBMAssemble);

             PetscBarrier(PETSC_NULL);
           //Write IBM_Shear Stress
             Write_IBM_Shear_Stress(user[bi].ibm->IBMAssemble);
}


Please add to "control.dat" the following "-IBM_Shear 3"

//-------------------------------Add to the "variables.h"--------------------------------------------
//Tensor form for the shear
typedef struct {
PetscScalar   Txx, Txy, Txz, Tyx, Tyy, Tyz, Tzx, Tzy, Tzz;
} Tensor;

//----------------------------------*************************************_---------------------------
//------------------------------- Subroutines--------------------------------------------------------


PetscErrorCode  IBM_Shear_Stress(UserCtx *user, IBMNodes *ibm)
{
DA	        da = user->da, fda = user->fda;
DALocalInfo	info = user->info;
PetscInt	xs = info.xs, xe = info.xs + info.xm;
PetscInt  	ys = info.ys, ye = info.ys + info.ym;
PetscInt	zs = info.zs, ze = info.zs + info.zm;
PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
PetscInt        lxs, lxe, lys, lye, lzs, lze;
PetscReal       rei= 1./user->ren;
PetscInt	i, j, k, elmt;
Cmpnts          ***ucat,***cent;
PetscReal       ***nvert;
PetscReal       dwdz,dwdy,dwdx;
PetscReal       dvdz,dvdy,dvdx;
PetscReal       dudz,dudy,dudx;
PetscReal       A_x, A_y, A_z, Total_Ax, Total_Ay, Total_Az;
PetscReal       Cs_x, Cs_y , Cs_z; 
Tensor          Strain_Rate, Strain_Rate1, Strain_Rate2,Strain_Rate3,Strain_Surface, Strain_Fluid;
PetscReal       nfx,nfy,nfz;
Cmpnts           ***csi,***eta,***zet;
PetscReal       csi1,csi2,csi3;
PetscReal       eta1,eta2,eta3;
PetscReal       zet1,zet2,zet3;
PetscReal       ***iaj,***jaj,***kaj;
IBMInfo         *ibminfo;
IBMListNode     *current;
PetscInt        rank, processors;
PetscInt        ip1,ip2,ip3;
PetscInt        jp1,jp2,jp3;
PetscInt        kp1,kp2,kp3;
PetscReal       sb,sc;
PetscReal       sk1,sk2,sk3;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

MPI_Comm_rank(PETSC_COMM_WORLD, &rank); 
MPI_Comm_size(PETSC_COMM_WORLD, &processors);


DAVecGetArray(fda, user->lCent, &cent);
DAVecGetArray(fda, user->lUcat, &ucat);
DAVecGetArray(da, user->lNvert, &nvert);
DAVecGetArray(fda, user->lCsi, &csi);
DAVecGetArray(fda, user->lEta, &eta);
DAVecGetArray(fda, user->lZet, &zet);
DAVecGetArray(da, user->lIAj, &iaj);
DAVecGetArray(da, user->lJAj, &jaj);
DAVecGetArray(da, user->lKAj, &kaj);

  current            = user->ibmlist.head;

  while (current) 
  {
    ibminfo = &current->ibm_intp;
    current = current->next;	

    //current IBM Node
    i = ibminfo->ni; 
    j= ibminfo->nj; 
    k = ibminfo->nk;

    //Point to cell
    elmt = ibminfo->cell;

    // normal 
    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt]; 

    sb = ibminfo->d_s; 
    sc = sb + ibminfo->d_i;


    if (i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze)
	{


	  Velocity_Gradient(user, i,j,k, &Strain_Rate);

	  //Find 3 fluid points
	  ip1 = ibminfo->i1; jp1 = ibminfo->j1; kp1 = ibminfo->k1;
	  ip2 = ibminfo->i2; jp2 = ibminfo->j2; kp2 = ibminfo->k2;
	  ip3 = ibminfo->i3; jp3 = ibminfo->j3; kp3 = ibminfo->k3;
	  sk1  = ibminfo->cr1; sk2 = ibminfo->cr2; sk3 = ibminfo->cr3;
 
	  Velocity_Gradient(user, ip1,jp1,kp1, &Strain_Rate1);
	  Velocity_Gradient(user, ip2,jp2,kp2, &Strain_Rate2);
	  Velocity_Gradient(user, ip3,jp3,kp3, &Strain_Rate3);

	  Strain_Fluid.Txx =   Strain_Rate1.Txx * sk1 + Strain_Rate2.Txx * sk2 + Strain_Rate3.Txx * sk3;
	  Strain_Fluid.Txy =   Strain_Rate1.Txy * sk1 + Strain_Rate2.Txy * sk2 + Strain_Rate3.Txy * sk3;
	  Strain_Fluid.Txz =   Strain_Rate1.Txz * sk1 + Strain_Rate2.Txz * sk2 + Strain_Rate3.Txz * sk3;
	  Strain_Fluid.Tyy =   Strain_Rate1.Tyy * sk1 + Strain_Rate2.Tyy * sk2 + Strain_Rate3.Tyy * sk3;
	  Strain_Fluid.Tzz =   Strain_Rate1.Tzz * sk1 + Strain_Rate2.Tzz * sk2 + Strain_Rate3.Tzz * sk3;
	  Strain_Fluid.Tyz =   Strain_Rate1.Tyz * sk1 + Strain_Rate2.Tyz * sk2 + Strain_Rate3.Tyz * sk3;

          //Interpolate on the surface element Strain Rate
          Strain_Surface.Txx = Strain_Fluid.Txx + (Strain_Rate.Txx - Strain_Fluid.Txx) * sc / (sc - sb);
          Strain_Surface.Txy = Strain_Fluid.Txy + (Strain_Rate.Txy - Strain_Fluid.Txy) * sc / (sc - sb);
          Strain_Surface.Txz = Strain_Fluid.Txz + (Strain_Rate.Txz - Strain_Fluid.Txz) * sc / (sc - sb);
          Strain_Surface.Tyy = Strain_Fluid.Tyy + (Strain_Rate.Tyy - Strain_Fluid.Tyy) * sc / (sc - sb);
          Strain_Surface.Tzz = Strain_Fluid.Tzz + (Strain_Rate.Tzz - Strain_Fluid.Tzz) * sc / (sc - sb);
          Strain_Surface.Tyz = Strain_Fluid.Tyz + (Strain_Rate.Tyz - Strain_Fluid.Tyz) * sc / (sc - sb);
	  //-----------------------------------------------

          ibminfo->Shear.Txx = Strain_Surface.Txx;
          ibminfo->Shear.Txy = Strain_Surface.Txy;
          ibminfo->Shear.Txz = Strain_Surface.Txz;
          ibminfo->Shear.Tyy = Strain_Surface.Tyy;
          ibminfo->Shear.Tzz = Strain_Surface.Tzz;
          ibminfo->Shear.Tyz = Strain_Surface.Tyz;

	}//End of inside domain

  }//End of IBM Nodes Loop

/*   Restore Working arrays */

  DAVecRestoreArray(fda, user->lUcat, &ucat);
  DAVecRestoreArray(da, user->lNvert, &nvert);
  DAVecRestoreArray(fda, user->lCsi, &csi);
  DAVecRestoreArray(fda, user->lEta, &eta);
  DAVecRestoreArray(fda, user->lZet, &zet);
  DAVecRestoreArray(da, user->lIAj, &iaj);
  DAVecRestoreArray(da, user->lJAj, &jaj);
  DAVecRestoreArray(da, user->lKAj, &kaj);
  DAVecRestoreArray(fda, user->lCent, &cent);


  
  return 0;
}

PetscErrorCode  Velocity_Gradient(UserCtx *user, PetscInt i, PetscInt j, PetscInt k, Tensor *Vel_Gradient)
{
DA	        da = user->da, fda = user->fda;
DALocalInfo	info = user->info;
PetscInt	xs = info.xs, xe = info.xs + info.xm;
PetscInt  	ys = info.ys, ye = info.ys + info.ym;
PetscInt	zs = info.zs, ze = info.zs + info.zm;
PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
PetscInt        lxs, lxe, lys, lye, lzs, lze;
PetscReal       rei= 1./user->ren;
Cmpnts          ***ucat,***cent;
PetscReal       ***nvert;
PetscReal       dwdz=0,dwdy=0,dwdx=0;
PetscReal       dvdz=0,dvdy=0,dvdx=0;
PetscReal       dudz=0,dudy=0,dudx=0;
PetscReal       Txx,Tyy,Tzz;
PetscReal       Tyx,Tzx,Tzy;
Cmpnts           ***csi,***eta,***zet;
PetscReal       csi1,csi2,csi3;
PetscReal       eta1,eta2,eta3;
PetscReal       zet1,zet2,zet3;
PetscReal       ***iaj,***jaj,***kaj;
PetscReal       ***p;
PetscReal       dw,dv,du;
Vec             Ucat;
PetscInt        mode;
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


DAVecGetArray(fda, user->lCent, &cent);
DAVecGetArray(fda, user->lUcat, &ucat);
DAVecGetArray(da, user->lNvert, &nvert);
DAVecGetArray(fda, user->lICsi, &csi);
DAVecGetArray(fda, user->lJEta, &eta);
DAVecGetArray(fda, user->lKZet, &zet);
DAVecGetArray(da, user->lIAj, &iaj);
DAVecGetArray(da, user->lJAj, &jaj);
DAVecGetArray(da, user->lKAj, &kaj);
DAVecGetArray(da, user->P, &p);


	dwdz = 0.;
	dvdz = 0.;
	dudz = 0.;

	dwdy = 0.;
	dvdy = 0.;
	dudy = 0.;

	dwdx = 0.;
	dvdx = 0.;
	dudx = 0.;

        du = 0;
        dv = 0;  
        dw = 0;


		zet1 = 0.25*(zet[k][j][i].x+zet[k-1][j][i].x)*	(kaj[k][j][i]+kaj[k-1][j][i]);
		zet2 = 0.25*(zet[k][j][i].y+zet[k-1][j][i].y)*  (kaj[k][j][i]+kaj[k-1][j][i]);
		zet3 = 0.25*(zet[k][j][i].z+zet[k-1][j][i].z)*  (kaj[k][j][i]+kaj[k-1][j][i]);


		//Finding the differenting scheme
		mode = -1;


		if (k+2 < lze+2 && k < mz-3 && mode == -1)        
		  {
		    if ( nvert[k+1][j][i]<=1 && nvert[k+2][j][i]<=1)
		      mode = 0;
		  } 

		if (k-2>=lzs-2 && zs > 2 && mode == -1)
		  {
		    if ( nvert[k-1][j][i]<=1 && nvert[k-2][j][i]<=1)
		      mode = 1;
		  } 
           
		if (k-1>=lzs-2 && k+1<lze+2 && k > 1 && k < mz-2 && mode == -1)
		  {
		    
		    if (nvert[k+1][j][i] +  nvert[k-1][j][i]<=1) //Central
		      mode =2;
		  }


		if  (k+1<lze+2 && k <mz-2 && mode == -1)
		  {
		    if (nvert[k+1][j][i]<=1)
		      mode = 3;
		    
		   } 
		

		if (k-1>=lzs-2  && k > 1 && mode == -1)
		  {
		if  (nvert[k-1][j][i]<=1)
		  mode = 4;
		
		  }


	   //Now do the real differentiation
	   // 0 : 2nd forward
	   // 1 : 2nd backward
	   // 2 : Central
	   // 3 : 1st forward
	   // 4 : 1st backward

	   switch (mode) {

	   case(0): {
	       dw = (-3*ucat[k][j][i].z + 4 * ucat[k+1][j][i].z - ucat[k+2][j][i].z)/2;
	       dv = (-3*ucat[k][j][i].y + 4 * ucat[k+1][j][i].y - ucat[k+2][j][i].y)/2;
	       du = (-3*ucat[k][j][i].x + 4 * ucat[k+1][j][i].x - ucat[k+2][j][i].x)/2;	     
              break;
	   }

	   case(1): {
		 dw = (3*ucat[k][j][i].z - 4 * ucat[k-1][j][i].z + ucat[k-2][j][i].z)/2;
		 dv = (3*ucat[k][j][i].y - 4 * ucat[k-1][j][i].y + ucat[k-2][j][i].y)/2;
		 du = (3*ucat[k][j][i].x - 4 * ucat[k-1][j][i].x + ucat[k-2][j][i].x)/2;
              break;	     
	   }

	   case(2): {
		   dw = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2;
		   dv = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2;
		   du = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2;
              break;
	   }
	   case(3):{
		     dw = (ucat[k+1][j][i].z - ucat[k][j][i].z);
		     dv = (ucat[k+1][j][i].y - ucat[k][j][i].y);
		     du = (ucat[k+1][j][i].x - ucat[k][j][i].x);
              break;

	   }
	   case(4):{
	               dw = ucat[k][j][i].z - ucat[k-1][j][i].z;
		       dv = ucat[k][j][i].y - ucat[k-1][j][i].y;
		       du = ucat[k][j][i].x - ucat[k-1][j][i].x;
              break;
		
	     }
	   }




	   dwdz = dw*zet3;
	   dvdz = dv*zet3;
	   dudz = du*zet3;

	   dwdy = dw*zet2;
	   dvdy = dv*zet2;
	   dudy = du*zet2;

	   dwdx = dw*zet1;
	   dvdx = dv*zet1;
	   dudx = du*zet1;


	       eta1 = 0.25*(eta[k][j][i].x+eta[k][j-1][i].x)*    (jaj[k][j][i]+jaj[k][j-1][i]);
	       eta2 = 0.25*(eta[k][j][i].y+eta[k][j-1][i].y)*    (jaj[k][j][i]+jaj[k][j-1][i]);
	       eta3 = 0.25*(eta[k][j][i].z+eta[k][j-1][i].z)*    (jaj[k][j][i]+jaj[k][j-1][i]);


	     du = 0;
	     dv = 0;
	     dw = 0;

	     //Finding the differenting scheme
	     mode = -1;

	     if (j+2 < lye+2 && j < my -3 && mode == -1)
	       {

		 if ( nvert[k][j+1][i] <=1 &&  nvert[k][j+2][i] <=1)
		   mode = 0;
	       }
	  

	     if (j -2 >= lys-2 && j > 2 && mode == -1)
	       {
		 if (nvert[k][j-1][i] <=1 && nvert[k][j-2][i] <=1)
		   mode = 1;
	       }
	    

	     if (j -1 >= lys-2 && j+1 <lye+2 && j > 1 && j < mz-2 &&  mode == -1)
	       {
		 if (nvert[k][j+1][i] + nvert[k][j-1][i]<=1)  //Central
		   mode = 2;

	       }
	
	     if (j+1 <lye+2  && j < mz-2 && mode == -1)
	       { 
		 if (nvert[k][j+1][i]<=1) 
		   mode = 3;
	       }
	       
	     if (j -1 >=lys-2 && j > 1 && mode == -1)
	       {
		 if (nvert[k][j-1][i]<=1 ) 	  
		   mode = 4;

	       }
	  
	   //Now do the real differentiation
	   // 0 : 2nd forward
	   // 1 : 2nd backward
	   // 2 : Central
	   // 3 : 1st forward
	   // 4 : 1st backward
	     switch (mode) {

	     case(0): {
	    dw = (-3*ucat[k][j][i].z + 4 * ucat[k][j+1][i].z - ucat[k][j+2][i].z)/2;
	    dv = (-3*ucat[k][j][i].y + 4 * ucat[k][j+1][i].y - ucat[k][j+2][i].y)/2;
	    du = (-3*ucat[k][j][i].x + 4 * ucat[k][j+1][i].x - ucat[k][j+2][i].x)/2.;
              break;
	     }

	     case(1):{
       	      dw = (3*ucat[k][j][i].z - 4 * ucat[k][j-1][i].z + ucat[k][j-2][i].z)/2;
	      dv = (3*ucat[k][j][i].y - 4 * ucat[k][j-1][i].y + ucat[k][j-2][i].y)/2;
	      du = (3*ucat[k][j][i].x - 4 * ucat[k][j-1][i].x + ucat[k][j-2][i].x)/2.;	       
              break;
	     }
	     case(2):{
	
		dw = (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2;
		dv = (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2;
		du = (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.;
              break;
	     }
	     case(3):{
	
	       dw = (ucat[k][j+1][i].z - ucat[k][j][i].z);
	       dv = (ucat[k][j+1][i].y - ucat[k][j][i].y);
	       du = (ucat[k][j+1][i].x - ucat[k][j][i].x);
              break;

	     }
	       case(4):{
	
		 dw = (ucat[k][j][i].z - ucat[k][j-1][i].z);
		 dv = (ucat[k][j][i].y - ucat[k][j-1][i].y);
		 du = (ucat[k][j][i].x - ucat[k][j-1][i].x);
              break;
	       }
	     }
	     

	dwdz += dw*eta3;
	dvdz += dv*eta3;
	dudz += du*eta3;

	dwdy += dw*eta2;
	dvdy += dv*eta2;
	dudy += du*eta2;

	dwdx += dw*eta1;
	dvdx += dv*eta1;
	dudx += du*eta1;



	du = 0;
	dv = 0;
	dw = 0;

	    csi1 = 0.25*(csi[k][j][i].x+csi[k][j][i-1].x)* (iaj[k][j][i]+iaj[k][j][i-1]);
	    csi2 = 0.25*(csi[k][j][i].y+csi[k][j][i-1].y)* (iaj[k][j][i]+iaj[k][j][i-1]);
	    csi3 = 0.25*(csi[k][j][i].z+csi[k][j][i-1].z)* (iaj[k][j][i]+iaj[k][j][i-1]);

	    mode = -1;

	    if (i+2 < lxe+2  && i < mx -3 && mode == -1)
	      {

		if (nvert[k][j][i+1]<=1 && nvert[k][j][i+2]<=1)
		  mode =0;
	      }


	    if (i-2>=lxs-2 &&  i > 2 && mode == -1)
	      {
		if (nvert[k][j][i-1]<=1 && nvert[k][j][i-2]<=1) 
		  mode = 1;
	      }


	    if (i-1 >= lxs-2 && i+1 < lxe+2 &&  i > 1 && i < mx-2 && mode == -1)
	      {
		if (nvert[k][j][i+1]+ nvert[k][j][i-1]<=1) 
		  mode = 2;
	      }

	    if (i+1 < lxe+2 && i < mx-2 && mode == -1)
	      {
		if (nvert[k][j][i+1]<=1)
		  mode = 3;
	      }
	      
	    if (i-1>=lxs-2  && i > 1 && mode == -1)
	      {
		if (nvert[k][j][i-1]<=1) 
		  mode = 4;
	      }
	      
	    //Now do the real differentiation
	   // 0 : 2nd forward
	   // 1 : 2nd backward
	   // 2 : Central
	   // 3 : 1st forward
	   // 4 : 1st backward

	     switch (mode) {

	     case(0): {
	       dw = (-3*ucat[k][j][i].z + 4 * ucat[k][j][i+1].z -  ucat[k][j][i+2].z)/2;
	       dv = (-3*ucat[k][j][i].y + 4 * ucat[k][j][i+1].y -  ucat[k][j][i+2].y)/2;
	       du = (-3*ucat[k][j][i].x + 4 * ucat[k][j][i+1].x -  ucat[k][j][i+2].x)/2;
              break;
	     }

	     case(1):{

	       dw = (3*ucat[k][j][i].z - 4 * ucat[k][j][i-1].z +  ucat[k][j][i-2].z)/2;
	       dv = (3*ucat[k][j][i].y - 4 * ucat[k][j][i-1].y +  ucat[k][j][i-2].y)/2;
	       du = (3*ucat[k][j][i].x - 4 * ucat[k][j][i-1].x +  ucat[k][j][i-2].x)/2;
              break;
	     }
	     case(2):{
	
	       dw = (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2;
               dv = (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2;
   	       du = (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2;
              break;
	     }
	     case(3):{
		  dw = ucat[k][j][i+1].z - ucat[k][j][i].z;
		  dv = ucat[k][j][i+1].y - ucat[k][j][i].y;
		  du = ucat[k][j][i+1].x - ucat[k][j][i].x;
              break;
	     }
	     case(4):{
	
		  dw = ucat[k][j][i].z - ucat[k][j][i-1].z;
		  dv = ucat[k][j][i].y - ucat[k][j][i-1].y;
		  du = ucat[k][j][i].x - ucat[k][j][i-1].x;
              break;
	     }
	     }
	  	      

	dwdz += dw*csi3;
	dvdz += dv*csi3;
	dudz += du*csi3;

	dwdy += dw*csi2;
	dvdy += dv*csi2;
	dudy += du*csi2;

	dwdx += dw*csi1;
	dvdx += dv*csi1;
	dudx += du*csi1;

 
	Tzz = rei * (dwdz + dwdz);
	Tyy = rei * (dvdy + dvdy);
	Txx = rei * (dudx + dudx);
	Tzy = rei * (dwdy + dvdz);
	Tzx = rei * (dwdx + dudz);
	Tyx = rei * (dvdx + dudy);

        Vel_Gradient->Txx = Txx;
        Vel_Gradient->Txy = Tyx;
        Vel_Gradient->Txz = Tzx;
        Vel_Gradient->Tyy = Tyy;
        Vel_Gradient->Tzz = Tzz;
        Vel_Gradient->Tyz = Tzy;


/*   Restore Working arrays */

DAVecRestoreArray(fda,user->lUcat, &ucat);
DAVecRestoreArray(da, user->lNvert, &nvert);
DAVecRestoreArray(fda, user->lICsi, &csi);
DAVecRestoreArray(fda, user->lJEta, &eta);
DAVecRestoreArray(fda, user->lKZet, &zet);
DAVecRestoreArray(da, user->lIAj, &iaj);
DAVecRestoreArray(da, user->lJAj, &jaj);
DAVecRestoreArray(da, user->lKAj, &kaj);
DAVecRestoreArray(fda, user->lCent, &cent);
DAVecRestoreArray(da, user->P, &p);
  return 0;
}

PetscInt IB_Nodes_Gross(UserCtx *user)
{
IBMInfo         *ibminfo;
IBMListNode     *current;
PetscInt        Number_Of_IB_Nodes= 0;

  //Point to the first node
  current = user->ibmlist.head;
  Number_Of_IB_Nodes= 0;

  while (current) 
  {
    ibminfo = &current->ibm_intp;
    current = current->next;	

    Number_Of_IB_Nodes++;
  }//End of all IBM List nodes

  return Number_Of_IB_Nodes;

} 
PetscErrorCode Shear_Projection(UserCtx *user, IBMNodes *ibm)
{
IBMInfo         *ibminfo;
IBMListNode     *current;
PetscInt        *current_element;
DA	        da = user->da, fda = user->fda;
PetscReal       *Txx_a,*Tyy_a,*Tzz_a;
PetscReal       *Tzy_a,*Tzx_a,*Tyx_a;
PetscInt        *nIB;

PetscInt        rank, processors;
PetscInt        Number_Of_IB_Nodes= 0;
int             *rec,*dis;

PetscReal       Txx,Tyy,Tzz;
PetscReal       Tzy,Tzx,Tyx;

PetscReal       Cs_x, Cs_y, Cs_z;
PetscReal       A_x, A_y, A_z;

PetscReal       *Cs_x_a, *Cs_y_a, *Cs_z_a;
PetscReal       *Cs_x_buffer, *Cs_y_buffer, *Cs_z_buffer;

PetscReal      *A_x_a,*A_y_a,*A_z_a;
PetscReal      *A_x_buffer,*A_y_buffer, *A_z_buffer;

PetscInt        elmt,i,j;

PetscInt        *element_buffer;
PetscReal       *Txx_buffer,*Tyy_buffer,*Tzz_buffer;
PetscReal       *Tzy_buffer,*Tzx_buffer,*Tyx_buffer;
PetscInt         Total_IB_Nodes;
PetscInt         start=0,*stride;
MPI_Datatype     stype, ltype;
Cmpnts           ***cent;

 DAVecGetArray(fda, user->lCent, &cent);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); 
  MPI_Comm_size(PETSC_COMM_WORLD, &processors);


  //Now create the buffer on ROOT to receive the sending information
  //Sum up all the size of IB Nodes on ALL processors
  PetscMalloc(processors*sizeof(PetscInt),&rec);
  PetscMalloc(processors*sizeof(PetscInt),&dis);
  PetscMalloc(processors*sizeof(PetscInt),&stride);

  Number_Of_IB_Nodes= IB_Nodes_Gross(user);

  //Find total number of IB Nodes in all processors
  MPI_Allreduce(&Number_Of_IB_Nodes, &Total_IB_Nodes,     1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  //Allocate the stride which is possible for allocation of information
  MPI_Gather(&Number_Of_IB_Nodes, 1, MPI_INT,  stride, 1, MPI_INT,   0,MPI_COMM_WORLD);

  //Broadcast the information of IB Nodes for all machines
  MPI_Bcast(stride, processors, MPI_INT, 0,MPI_COMM_WORLD);

 //Allocate the rec and dis
  start = 0;

  for (i=0; i< processors; i++)
  {
    dis[i] = start;
    rec[i] = stride[i];
    start = start + stride[i];
  }

  // Here it is the part that send and receive the information from other processors
  // To send all the information to ROOT
  //Allocate the memory for shear stress of IB object
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscInt),  &current_element);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Txx_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Tyx_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Tzx_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Tyy_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Tzz_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Tzy_a);


   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Cs_x_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Cs_y_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &Cs_z_a);

   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &A_x_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &A_y_a);
   PetscMalloc(Number_Of_IB_Nodes*sizeof(PetscReal), &A_z_a);


   PetscMalloc(Total_IB_Nodes*sizeof(PetscInt),  &element_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Txx_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Tyx_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Tzx_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Tyy_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Tzz_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Tzy_buffer);

   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Cs_x_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Cs_y_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &Cs_z_buffer);

   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &A_x_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &A_y_buffer);
   PetscMalloc(Total_IB_Nodes*sizeof(PetscReal), &A_z_buffer);


 //Point to the first node
  current = user->ibmlist.head;
  i = 0;

  while (current) 
  {
     ibminfo = &current->ibm_intp;
     current = current->next;	

    //Point to cell
     elmt = ibminfo->cell;

     if (i < Number_Of_IB_Nodes)
       {
         current_element[i] = elmt;

	 //Stress components
	 Txx_a[i] = ibminfo->Shear.Txx;
	 Tyx_a[i] = ibminfo->Shear.Txy;
	 Tzx_a[i] = ibminfo->Shear.Txz;
	 Tyy_a[i] = ibminfo->Shear.Tyy;
	 Tzz_a[i] = ibminfo->Shear.Tzz;
	 Tzy_a[i] = ibminfo->Shear.Tyz;

	 //Traction components
	 Cs_x_a[i] = ibminfo->Cs_x;
	 Cs_y_a[i] = ibminfo->Cs_y;
	 Cs_z_a[i] = ibminfo->Cs_z;

	 //Area components
	 A_x_a[i] = ibminfo->Area_x;
	 A_y_a[i] = ibminfo->Area_y;
	 A_z_a[i] = ibminfo->Area_z;

       }


     i++;
  }//End of all IBM List nodes

  //Send number of element first
  MPI_Type_contiguous(Number_Of_IB_Nodes,  MPI_INT, &stype);
  MPI_Type_commit(&stype);

  ///Now send the information to ROOT
  MPI_Gatherv(current_element, 1, stype,   element_buffer,  rec, dis,   MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Type_contiguous(Number_Of_IB_Nodes,  MPI_DOUBLE, &ltype);
  MPI_Type_commit(&ltype);

  MPI_Gatherv(Txx_a,         1,  ltype,  Txx_buffer,             rec, dis ,  MPI_DOUBLE,  0, MPI_COMM_WORLD);

  MPI_Gatherv(Tyx_a,         1,  ltype,  Tyx_buffer,             rec, dis,   MPI_DOUBLE,  0, MPI_COMM_WORLD);

  MPI_Gatherv(Tzx_a,         1,  ltype,  Tzx_buffer,             rec, dis,   MPI_DOUBLE , 0, MPI_COMM_WORLD);

  MPI_Gatherv(Tyy_a,         1,  ltype,  Tyy_buffer,             rec, dis ,  MPI_DOUBLE , 0, MPI_COMM_WORLD);

  MPI_Gatherv(Tzz_a,         1,  ltype,  Tzz_buffer,             rec, dis,   MPI_DOUBLE , 0, MPI_COMM_WORLD);

  MPI_Gatherv(Tzy_a,         1,  ltype,  Tzy_buffer,             rec, dis,   MPI_DOUBLE , 0, MPI_COMM_WORLD);

  //Traction gathering
    MPI_Gatherv(Cs_x_a,        1, ltype, Cs_x_buffer,       rec, dis, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
    MPI_Gatherv(Cs_y_a,        1, ltype, Cs_y_buffer,       rec, dis, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
    MPI_Gatherv(Cs_z_a,        1, ltype, Cs_z_buffer,       rec, dis, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);

  //Area gathering
    MPI_Gatherv(A_x_a,        1, ltype, A_x_buffer,       rec, dis, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
    MPI_Gatherv(A_y_a,        1, ltype, A_y_buffer,       rec, dis, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
    MPI_Gatherv(A_z_a,        1, ltype, A_z_buffer,       rec, dis, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);

  // if it is ROOT then process the received information
  // Wait for all processors to proceed further
   
   
 if (rank == 0)
   {

     PetscMalloc(ibm->n_elmt*sizeof(PetscInt), &nIB);
     //Clear the stress information
     for (i= 0; i< ibm->n_elmt; i++)
       {
	     ibm->Txx[i] = 0;
	     ibm->Txy[i] = 0;
	     ibm->Txz[i] = 0;
	     ibm->Tyy[i] = 0;
	     ibm->Tzz[i] = 0;
	     ibm->Tyz[i] = 0;
	     nIB[i] = 0;

	     //Traction vector
	     ibm->Traction_x[i] = 0;
	     ibm->Traction_y[i] = 0;
	     ibm->Traction_z[i] = 0;

	     //Area
	     ibm->Area_x[i] = 0;
	     ibm->Area_y[i] = 0;
	     ibm->Area_z[i] = 0;

       }


     //Project onto the surface
     for (i= 0; i < Total_IB_Nodes; i++)
      {
     
	elmt = element_buffer[i];

	 if (elmt < ibm->n_elmt)
	  {

	     Txx = Txx_buffer[i];
	     Tyx = Tyx_buffer[i];
	     Tzx = Tzx_buffer[i];
	     Tyy = Tyy_buffer[i];
	     Tzz = Tzz_buffer[i];
	     Tzy = Tzy_buffer[i];

	     //Traction components
	     Cs_x = Cs_x_buffer[i];
	     Cs_y = Cs_y_buffer[i];
	     Cs_z = Cs_z_buffer[i];


	     //Area
	     A_x = A_x_buffer[i];
	     A_y = A_y_buffer[i];
	     A_z = A_z_buffer[i];


       	    //Transform information into IB cell's information

	    if ( fabs(Txx) + fabs(Tyx) + fabs(Tzx) + fabs(Tyy) + fabs(Tzz) + fabs(Tzy) > 0)
	     {
	     ibm->Txx[elmt] =ibm->Txx[elmt] + Txx;
	     ibm->Txy[elmt] =ibm->Txy[elmt] + Tyx;
	     ibm->Txz[elmt] =ibm->Txz[elmt] + Tzx;
	     ibm->Tyy[elmt] =ibm->Tyy[elmt] + Tyy;
	     ibm->Tzz[elmt] =ibm->Tzz[elmt] + Tzz;
	     ibm->Tyz[elmt] =ibm->Tyz[elmt] + Tzy;

	     nIB[elmt] = nIB[elmt] + 1;
	     }

	  }//End of element     	     
      }//End of for all IB Nodes

 
     PetscInt     shear_method;
     PetscTruth   flg;    

      PetscOptionsGetInt(PETSC_NULL, "-IBM_Shear", &shear_method, &flg);

     //Average the shear stress

     for (i = 0;i< ibm->n_elmt; i++)
       {

	 //Average over the same element
       if (shear_method == 3)
        {
	 if (nIB[i] > 0)
	   {
	     ibm->Txx[i] = ibm->Txx[i]/nIB[i];
     	     ibm->Txy[i] = ibm->Txy[i]/nIB[i];
             ibm->Txz[i] = ibm->Txz[i]/nIB[i];
	     ibm->Tyy[i] = ibm->Tyy[i]/nIB[i];
             ibm->Tzz[i] = ibm->Tzz[i]/nIB[i];
  	     ibm->Tyz[i] = ibm->Tyz[i]/nIB[i];

	   }
	}//End of shear method =3 -> Calculate from IB nodes


       }//End of all element

          for (i = 0;i< ibm->n_elmt; i++)
     {

        if (shear_method == 3)
        {
          PetscReal nx,ny,nz,tx,ty,tz,normal_stress;

          nx = ibm->nf_x[i];
          ny = ibm->nf_y[i];
          nz = ibm->nf_z[i];


	     //Traction vector compute
	  tx = ibm->Txx[i]*nx + ibm->Txy[i]*ny + ibm->Txz[i]*nz;
          ty = ibm->Txy[i]*nx + ibm->Tyy[i]*ny + ibm->Tyz[i]*nz;
          tz = ibm->Txz[i]*nx + ibm->Tyz[i]*ny + ibm->Tzz[i]*nz;

          normal_stress = tx*nx  +  ty*ny + tz*nz;

     
          ibm->Traction_x[i] = tx - normal_stress*nx;
          ibm->Traction_y[i] = ty - normal_stress*ny;
          ibm->Traction_z[i] = tz - normal_stress*nz;

        }
     }//End of shear_method == 3




       PetscFree(nIB);    

   }//End of rank == 0



    PetscFree(current_element);
    PetscFree(Txx_a);
    PetscFree(Tyx_a);
    PetscFree(Tzx_a);
    PetscFree(Tyy_a);
    PetscFree(Tzz_a);
    PetscFree(Tzy_a);

    PetscFree(element_buffer);
    PetscFree(Txx_buffer);
    PetscFree(Tyx_buffer);
    PetscFree(Tzx_buffer);
    PetscFree(Tyy_buffer);
    PetscFree(Tzz_buffer);
    PetscFree(Tzy_buffer);


    PetscFree(Cs_x_a);
    PetscFree(Cs_y_a);
    PetscFree(Cs_z_a);


    PetscFree(Cs_x_buffer);
    PetscFree(Cs_y_buffer);
    PetscFree(Cs_z_buffer);

    DAVecRestoreArray(fda, user->lCent, &cent);

  return 0;
}
//Write out the stress tensor T{ij}
PetscErrorCode Write_IBM_Shear_Stress(IBMNodes *ibm)
{

PetscInt        i;
FILE            *f;
PetscInt        rank;
char            filen[80];


  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (rank == 0) 
   {

     PetscInt     shear_method;
     PetscTruth   flg;    

     PetscOptionsGetInt(PETSC_NULL, "-IBM_Shear", &shear_method, &flg);


     if (shear_method == 3)
       {
          sprintf(filen,"IBM_Shear_%4.4d.dat",ti);

	  f = fopen(filen,"w");     

	  for (i= 0; i<ibm->n_elmt; i++)
	    { 
		 
	      PetscReal Tract_mag;

	      Tract_mag = sqrt(ibm->Traction_x[i]*ibm->Traction_x[i]+ 
			       ibm->Traction_y[i]*ibm->Traction_y[i]+ 
			       ibm->Traction_z[i]*ibm->Traction_z[i]);

	     fprintf(f,"%le  %le  %le %le\n",ibm->Traction_x[i],ibm->Traction_y[i],ibm->Traction_z[i],Tract_mag);


	    }

	  fclose(f);
       }//End of shear method 3

     if(shear_method == 4)
       {
	 //Print out the Twx component

         PetscPrintf(PETSC_COMM_WORLD,"Tzx is printing.......\n");

         sprintf(filen,"IBM_Shear_%4.4d.dat",ti);
         f = fopen(filen,"w");

         PetscReal Radius =1;

         for (i= 0; i<ibm->n_elmt; i++)
          {
	    fprintf(f,"%le %le %le %le  %le  %le\n",ibm->Txx[i], ibm->Txy[i], ibm->Txz[i], ibm->Tyy[i], ibm->Tzz[i],ibm->Tyz[i]);
          }

         fclose(f);
       }

  
   }//End of rank =0



return 0;
}