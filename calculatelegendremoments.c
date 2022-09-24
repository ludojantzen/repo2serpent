/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatelegendremoments.c                     */
/*                                                                           */
/* Created:       2017/04/06 (VVa)                                           */
/* Last modified: 2017/04/06 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Calculates the first seven Legendre moments for a scattering */
/*              cosine distribution.                                         */
/*                                                                           */
/* Comments: - The mu-distribution must be tabulated and use linear-linear   */
/*             interpolation.                                                */
/*           - This might be able to partly utilize polynomiallegendre.c     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "CalculateLegendreMoments:"

/*****************************************************************************/

void CalculateLegendreMoments(long np, long l0, double LMS[7])
{
  /* pointer l0 to the first cosine value (d1) i.e. after to the "np" value*/
  double unoInt;
  double muInt;
  double mu2Int;
  double mu3Int;
  double mu4Int;
  double mu5Int;
  double mu6Int;
  double mu7Int;
  double d1, d2;
  double p1, p2;
  /*
  double c1, c2;
  */
  double d12, d13, d14, d15, d16, d17, d18, d19;
  double d22, d23, d24, d25, d26, d27, d28;
  long k;
  double m;


  unoInt = 0.0;
  muInt = 0.0;
  mu2Int = 0.0;
  mu3Int = 0.0;
  mu4Int = 0.0;
  mu5Int = 0.0;
  mu6Int = 0.0;
  mu7Int = 0.0;

  for (k = 0; k < (np - 1); k++)
    {
      /* Cosine values */
      d1 = RDB[l0 + k];
      d2 = RDB[l0 + k + 1];

      /* PDF values */

      p1 = RDB[l0 + np + k];
      p2 = RDB[l0 + np + k + 1];

      /* CDF values */
      /*
      c1 = RDB[l0 + 2*np + k];
      c2 = RDB[l0 + 2*np + k + 1];
      */

      /* Slope of pdf over the cosine interval */

      m = (p2 - p1)/(d2 - d1);

      /* Powers of the limiting cosines */

      d12 = d1*d1;
      d13 = d12*d1;
      d14 = d13*d1;
      d15 = d14*d1;
      d16 = d15*d1;
      d17 = d16*d1;
      d18 = d17*d1;
      d19 = d18*d1;

      d22 = d2*d2;
      d23 = d22*d2;
      d24 = d23*d2;
      d25 = d24*d2;
      d26 = d25*d2;
      d27 = d26*d2;
      d28 = d27*d2;

      /* Integrals of (x^N times linear pdf) over the intervals */

      /*integral_d1^d2 x (p1+(x-d1) m) dx =
       * 1/6  (d1^3 m-3 d1^2 p1-3 d1 d2^2 m+d2^2 (2 d2 m+3 p1)) */
      muInt += (1.0/6.0)*  (d13*m - 3.0*d12*p1 - 3.0*d1*d22*m + d22*(2.0*d2*m + 3.0*p1));

      /*integral_d1^d2 x^2 (p1+(x-d1) m) dx =
       * 1/12 (d1^4 m-4 d1^3 p1-4 d1 d2^3 m+d2^3 (3 d2 m+4 p1)) */
      mu2Int += (1.0/12.0)*(d14*m - 4.0*d13*p1 - 4.0*d1*d23*m + d23*(3.0*d2*m + 4.0*p1));

      /*integral_d1^d2 x^3 (p1+(x-d1) m) dx =
       * 1/20 (d1^5 m-5 d1^4 p1-5 d1 d2^4 m+d2^4 (4 d2 m+5 p1)) */
      mu3Int += (1.0/20.0)*(d15*m - 5.0*d14*p1 - 5.0*d1*d24*m + d24*(4.0*d2*m + 5.0*p1));

      /*integral_d1^d2 x^4 (p1+(x-d1) m) dx =
       * 1/30 (d1^6 m-6 d1^5 p1-6 d1 d2^5 m+d2^5 (5 d2 m+6 p1)) */
      mu4Int += (1.0/30.0)*(d16*m - 6.0*d15*p1 - 6.0*d1*d25*m + d25*(5.0*d2*m + 6.0*p1));

      /*integral_d1^d2 x^5 (p1+(x-d1) m) dx =
       * 1/42 (d1^7 m-7 d1^6 p1-7 d1 d2^6 m+d2^6 (6 d2 m+7 p1)) */
      mu5Int += (1.0/42.0)*(d17*m - 7.0*d16*p1 - 7.0*d1*d26*m + d26*(6.0*d2*m + 7.0*p1));

      /*integral_d1^d2 x^6 (p1+(x-d1) m) dx =
       * 1/56 (d1^8 m-8 d1^7 p1-8 d1 d2^7 m+d2^7 (7 d2 m+8 p1)) */
      mu6Int += (1.0/56.0)*(d18*m - 8.0*d17*p1 - 8.0*d1*d27*m + d27*(7.0*d2*m + 8.0*p1));

      /*integral_d1^d2 x^7 (p1+(x-d1) m) dx =
       * 1/72 (d1^9 m-9 d1^8 p1-9 d1 d2^8 m+d2^8 (8 d2 m+9 p1)) */
      mu7Int += (1.0/72.0)*(d19*m - 9.0*d18*p1 - 9.0*d1*d28*m + d28*(8.0*d2*m + 9.0*p1));

      /* integral_d1^d2 1 (p1+(x-d1) m) dx =
       * 1/2 (d1-d2) (d1 m-d2 m-2 p1) */
      unoInt += (1.0/2.0)*(d1-d2)*(d1*m - d2*m - 2*p1);

    }

  /* Check the value of total integral of the linear-interpolated PDF */
  /* I.e. the value of CDF at 1.0 */

  CheckValue(FUNCTION_NAME, "unoInt", "", unoInt, 0.95, 1.05);

  /* Values of Legendre polynomials */

  LMS[0] = muInt;
  LMS[1] =  0.5*(3.0*mu2Int - 1.0);
  LMS[2] =  0.5*(5.0*mu3Int -   3.0*muInt);
  LMS[3] =     (35.0*mu4Int -  30.0*mu2Int + 3.0)/8.0;
  LMS[4] =     (63.0*mu5Int -  70.0*mu3Int +  15*muInt)/8.0;
  LMS[5] =    (231.0*mu6Int - 315.0*mu4Int + 105*mu2Int - 5.0)/16.0;
  LMS[6] =    (429.0*mu7Int - 693.0*mu5Int + 315*mu3Int - 35.0*muInt)/16.0;
}

/*****************************************************************************/
