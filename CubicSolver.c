static char help[] = "Testing programming!";


#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

#ifdef TECIO
#include "TECIO.h"
#endif

//!--------------------------------------------------------------------
PetscErrorCode SolveCubic(PetscReal I[4],PetscReal *x) 
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
  
  p = ( 3*I[2]/I[0] - pow((I[1]/I[0]), 2) ) / 3.;
  q = ( 2*pow((I[1]/I[0]),3) - 9*I[1]*I[2]/I[0]/I[0] + 27*I[3]/I[0] ) / 27.;

  D = pow(p/3,3) + pow(q/2,2);

  //  PetscPrintf(PETSC_COMM_WORLD, "solve %le %le %le\n", D,p,q);

  if (D<0.) {
    phi = acos(-q/2/sqrt(pow(fabs(p),3)/27));
    y1 =  2 * sqrt(fabs(p)/3) * cos(phi/3);
    y2 = -2 * sqrt(fabs(p)/3) * cos((phi+pi)/3);
    y3 = -2 * sqrt(fabs(p)/3) * cos((phi-pi)/3);
  } else if (fabs(D)<1e-6) {

    u = -pow(fabs(q)/2,1./3.);
    if (q<0) u=-u;

    //    PetscPrintf(PETSC_COMM_WORLD, "solve u %le\n", u);
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

  PetscPrintf(PETSC_COMM_WORLD, "solve \n");

  x[0] = y1 - I[1]/I[0]/3;//  PetscPrintf(PETSC_COMM_WORLD, "solve 1 \n");

  x[1] = y2 - I[1]/I[0]/3;//  PetscPrintf(PETSC_COMM_WORLD, "solve 2 \n");

  x[2] = y3 - I[1]/I[0]/3;

  PetscPrintf(PETSC_COMM_WORLD, "solve \n");
  
  return(0);
}

int main(int argc, char **argv)
{
  PetscReal  I[4],x[3];
  PetscErrorCode imag;
 
  PetscInitialize(&argc, &argv, (char *)0, help);
 
  PetscPrintf(PETSC_COMM_WORLD, "solve \n");

  I[0]=1;
  I[1]=2;
  I[2]=0;
  I[3]=-1;

  imag=SolveCubic(I,&x);

  PetscPrintf(PETSC_COMM_WORLD, "solve \n");
  
  PetscPrintf(PETSC_COMM_WORLD, "sol %le %le %le \n", x[0],x[1],x[2]);

  PetscFinalize();


}
