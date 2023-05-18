static char help[] = "Testing programming!";


#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

#ifdef TECIO
#include "TECIO.h"
#endif

int main(int argc, char **argv)
{

  PetscInt   njac, nrot;
  PetscReal  a[3][3], v[3][3], d[3];


  PetscInitialize(&argc, &argv, (char *)0, help);

  a[0][0]=1;
  a[0][1]=1;
  a[0][2]=2;
  a[1][0]=2;
  a[1][1]=2;
  a[1][2]=2;
  a[2][0]=2;
  a[2][1]=2;
  a[2][2]=1;

  njac = 3;
  
  d[0]=0.;
  d[1]=0.;
  d[2]=0.;

  PetscPrintf(PETSC_COMM_WORLD, "jacobi \n");
  
  jacobi(a, njac, &d, v, &nrot);
  PetscPrintf(PETSC_COMM_WORLD, "eig %d %le %le %le %le %le %le\n",nrot, d[0],d[1],d[2],v[0][0],v[0][2],v[2][0]);
  PetscPrintf(PETSC_COMM_WORLD, "A %d %le %le %le \n",nrot, a[0][0],a[0][1],a[0][2]);
  PetscPrintf(PETSC_COMM_WORLD, "A %d %le %le %le \n",nrot, a[1][0],a[1][1],a[1][2]);
  PetscPrintf(PETSC_COMM_WORLD, "A %d %le %le %le \n",nrot, a[2][0],a[2][1],a[2][2]);

  PetscPrintf(PETSC_COMM_WORLD, "v %d %le %le %le \n",nrot, v[0][0],v[0][1],v[0][2]);
  PetscPrintf(PETSC_COMM_WORLD, "v %d %le %le %le \n",nrot, v[1][0],v[1][1],v[1][2]);
  PetscPrintf(PETSC_COMM_WORLD, "v %d %le %le %le \n",nrot, v[2][0],v[2][1],v[2][2]);

  PetscPrintf(PETSC_COMM_WORLD, "sort \n");
  eigsrt(&d, v, njac, njac);
  
  PetscPrintf(PETSC_COMM_WORLD, "eig %le %le %le \n", d[0],d[1],d[2]);

  PetscFinalize();

}


//!--------------------------------------------------------------------
PetscErrorCode jacobi(PetscReal a[3][3],PetscInt n, 
		      PetscReal *d,PetscReal v[3][3],
		      PetscInt *nrot) 
//!--------------------------------------------------------------------  
{
 PetscPrintf(PETSC_COMM_WORLD, "jacobi \n");

  PetscInt  NMAX=500;  
  PetscInt  i, ip, iq, j;
  PetscReal c, g, h, s, sm,t,tau,theta,tresh,b[NMAX],z2[NMAX];

  PetscPrintf(PETSC_COMM_WORLD, "jacobi n %d \n", n);
 
  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) {
      v[ip][iq]=0.;
      v[ip][ip]=1.;
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "jacobi \n");
 
  for (ip=0; ip<n; ip++) {
    b[ip]=a[ip][ip];
    d[ip]=b[ip];
    z2[ip]=0.;
  }

  PetscPrintf(PETSC_COMM_WORLD, "jacobi \n");
 
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

    PetscPrintf(PETSC_COMM_WORLD, "jacobi %le %le \n", sm, tresh);
       
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

//!--------------------------------------------------------------------
PetscErrorCode eigsrt(PetscReal *d, PetscReal v[3][3], 
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
       
  PetscPrintf(PETSC_COMM_WORLD, "eig %le %le %le \n", d[0],d[1],d[2]);

  }

  return 0;
}



