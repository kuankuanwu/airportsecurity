#include "mex.h"


double mrg32k3a( double *s10, double *s11, double *s12, double *s20, double *s21, 
                 double *s22 );

void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	/* more C code ... */
	double s10 = 1., s11 = 2., s12 = 3., s20 = 4., s21 = 5.;    /* Used in mrg32k3a function */
	double s22;
	s22 = mxGetScalar(prhs[0]);
	plhs[0] = mxCreateDoubleScalar(mrg32k3a(&s10,&s11,&s12,&s20,&s21,&s22));
	plhs[1] = mxCreateDoubleScalar(s22);
}

double mrg32k3a( double *s10, double *s11, double *s12, double *s20, double *s21, 
                 double *s22 )
/**************************************************************************
c     double precision function mrg32k3a( s22 )
c     Reference: Pierre L'Ecuyer, "Good parameters and implementations 
c     for combined multiple recursive random number generators", 
c     Operations Research 47 (1999), pp. 159--164.
c     The state of the generator is (s10,s11,s12,s20,s21,s22).
c     If s22 is negative, then the other five seeds are initialized.
c
c     input:  s22          a positive integer less than m2.
c     output: mrg32k3a     a pseudorandom number in (0,1).
c
c     Original Fortran code: Bruce Schmeiser and Raghu Pasupathy, June 2012.
c     Converted to C by Huifen Chen, August 20, 2012.
***************************************************************************/
{

      int    k;
      double p1, p2, rn;
      double norm, m1, m2, a12, a13n, a21, a23n;
             norm = 2.328306549295728E-10;
             m1   = 4294967087.;
             m2   = 4294944443.;
             a12  = 1403580.;
             a13n = 810728.;
             a21  = 527612.;
             a23n = 1370589.;

//    printf("s10, s11, s12, s20, s21, s22 = %lf %lf %lf %lf %lf %lf\n", 
//           *s10,*s11,*s12,*s20,*s21,*s22);
      

//    ...initialize (if s22 is not positive)
      if (*s22 <= 0.) 
      {
         *s22 = -*s22;
         *s10 = *s22 + 1;
         if (*s10 >= m1)  *s10 =  1;
         *s11 = *s22 + 2;
         if (*s11 >= m1)  *s11 =  2;
         *s12 = *s22 + 3;
         if (*s12 >= m1)  *s12 =  3;
         *s20 = *s22 + 4;
         if (*s20 >= m2)  *s20 =  4;
         *s21 = *s22 + 5;
         if (*s21 >= m2)  *s21 =  5;
      }

//    ...component 1
      p1 = a12 * (*s11) - a13n * (*s10);
      k  = (int) ( p1 / m1 );
      p1 = p1 - k * m1;
      if (p1 < 0.) p1 = p1 + m1;
      *s10 = *s11;
      *s11 = *s12;
      *s12 = p1;
//    ...component 2
      p2 = a21 * (*s22) - a23n * (*s20);
      k  = (int) ( p2 / m2 );
      p2 = p2 - k * m2;
      if (p2 < 0.) p2 = p2 + m2;
      *s20 = *s21;
      *s21 = *s22;
      *s22 = p2;
//    ...combination
      if (p1 <= p2) 
      {
         rn = (p1 - p2 + m1) * norm;
      } else { 
         rn = (p1 - p2) * norm;    
      }

      return(rn);
}