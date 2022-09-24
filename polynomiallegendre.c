/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : polynomiallegendre.c                           */
/*                                                                           */
/* Created:       2017/02/20 (BWe)                                           */
/* Last modified: 2018/02/23 (BWe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Evaluate the Legendre polynomial                             */
/*                                                                           */
/* Comments: Vectorized version for calculating Legendre polynomials. The    */
/*           array 'vector' must be pre-allocated in the calling function    */
/*           with a size of n+1.                                             */
/*           Previously, values for n <= 10 were included in a switch-case   */
/*           block. However, testing revealed that this provided no          */
/*           measurable performance gains or losses whatsoever.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PolynomialLegendre:"

#define MAX_DIRECT_CALCULATION  12

#define Save(k, value)  (PutPrivateData(vectorPtr + k, value, threadID))
#define Load(k)         (GetPrivateData(vectorPtr + k, threadID))

/*****************************************************************************/

void PolynomialLegendre(long n, double x, long vectorPtr, long threadID, long type)
{
  long k;
  const double x2 = x * x;

  /* Use direct formula to efficiently evaluate the polynomials for n <= 12   */
  /*                                                                          */
  /* The performance benefit diminishes for higher n. It is expected that     */
  /* the cost of the direct calculation nears that of the recurrence relation */
  /* in the neighborhood of n == 15, although this theory is untested due to  */
  /* only implementing the direct calculations up to n == 12.                 */
  /*                                                                          */
  /* If you want to calculate the higher-order Legendre Coefficients and      */
  /* code them in then be my guest.                                           */

  switch (type)
  {
    case FET_CALCULATION_STANDARD:
    {
      switch (n)
      {
        default:
        case MAX_DIRECT_CALCULATION:  /* 12 */
          Save(12, ((((((676039 * x2 - 1939938) * x2 + 2078505) * x2 - 1021020) * x2 + 225225) * x2 - 18018) * x2 + 231) / 1024);

        case 11:
          Save(11, (((((88179 * x2 - 230945) * x2 + 218790) * x2 - 90090) * x2 + 15015) * x2 - 693) * x / 256);

        case 10:
          Save(10, (((((46189 * x2 - 109395) * x2 + 90090) * x2 - 30030) * x2 + 3465) * x2 - 63) / 256);

        case 9:
          Save(9, ((((12155 * x2 - 25740) * x2 + 18018) * x2 - 4620) * x2 + 315) * x / 128);

        case 8:
          Save(8, ((((6435 * x2 - 12012) * x2 + 6930) * x2 - 1260) * x2 + 35) / 128);

        case 7:
          Save(7, (((429 * x2 - 693) * x2 + 315) * x2 - 35) * x / 16);

        case 6:
          Save(6, (((231 * x2 - 315) * x2 + 105) * x2 - 5) / 16);

        case 5:
          Save(5, ((63 * x2 - 70) * x2 + 15) * x / 8);

        case 4:
          Save(4, ((35 * x2 - 30) * x2 + 3) / 8);

        case 3:
          Save(3, (5 * x2 - 3) * x / 2);

        case 2:
          Save(2, (3 * x2 - 1) / 2);

        case 1:
          Save(1, x);

        case 0:
          Save(0, 1);
      }

      /* Evaluate any remaining polynomials.                                  */
      /* The original recurrence relation is:                                 */
      /*       (2 * k - 1) * x * L_(k-1) - (k - 1) * L_(k-2)                  */
      /* L_k = ---------------------------------------------                  */
      /*                        k                                             */

      for (k = MAX_DIRECT_CALCULATION + 1; k <= n; ++k)
        Save(k, ((k + k - 1) * x * Load(k - 1) - (k - 1) * Load(k - 2)) / (double)k);

      break;
    }

    case FET_CALCULATION_ORTHONORMAL:
    {
      switch (n)
      {
        default:
        case MAX_DIRECT_CALCULATION:  /* 12 */
          Save(12, ((((((676039 * x2 - 1939938) * x2 + 2078505) * x2 - 1021020) * x2 + 225225) * x2 - 18018) * x2 + 231) / 1024
                   * 12.5);

        case 11:
          Save(11, (((((88179 * x2 - 230945) * x2 + 218790) * x2 - 90090) * x2 + 15015) * x2 - 693) * x / 256
                   * 11.5);

        case 10:
          Save(10, (((((46189 * x2 - 109395) * x2 + 90090) * x2 - 30030) * x2 + 3465) * x2 - 63) / 256
                   * 10.5);

        case 9:
          Save(9, ((((12155 * x2 - 25740) * x2 + 18018) * x2 - 4620) * x2 + 315) * x / 128
                  * 9.5);

        case 8:
          Save(8, ((((6435 * x2 - 12012) * x2 + 6930) * x2 - 1260) * x2 + 35) / 128
                  * 8.5);

        case 7:
          Save(7, (((429 * x2 - 693) * x2 + 315) * x2 - 35) * x / 16
                  * 7.5);

        case 6:
          Save(6, (((231 * x2 - 315) * x2 + 105) * x2 - 5) / 16
                  * 6.5);

        case 5:
          Save(5, ((63 * x2 - 70) * x2 + 15) * x / 8
                  * 5.5);

        case 4:
          Save(4, ((35 * x2 - 30) * x2 + 3) / 8
                  * 4.5);

        case 3:
          Save(3, (5 * x2 - 3) * x / 2
                  * 3.5);

        case 2:
          Save(2, (3 * x2 - 1) / 2
                  * 2.5);

        case 1:
          Save(1, x
                  * 1.5);

        case 0:
          Save(0, 1
                  * 0.5);
      }

      /* Evaluate any remaining polynomials.                                  */
      /*                                                                      */
      /* For this FET we are using the orthonormalized version of the         */
      /* polynomials, so each polynomial L_k is multiplied by:                */
      /*    (2 * k + 1)                                                       */
      /*    -----------      essentially:    k + 0.5                          */
      /*         2                                                            */
      /* Reversing this in the standard polynomials and implementing for the  */
      /* current polynomial results in the orthonormalized recurrence:        */
      /*       (2 * k + 1)   /                 (k - 1)             \          */
      /* L_k = ----------- * | x * L_(k-1) - ----------- * L_(k-2) |          */
      /*            k        \               (2 * k - 3)           /          */
      /*                                                                      */
      /* The options are 1) to use this form, or 2) to not apply the          */
      /* orthonormalization at first, and then loop through all the values in */
      /* a second loop and then apply the orthonormalization.                 */

      for (k = MAX_DIRECT_CALCULATION + 1; k <= n; ++k)
        Save(k, ((k + k + 1) / (double)k)
                * (x * Load(k - 1) - ((k - 1) / (double)(k + k - 3)) * Load(k - 2)));

      break;
    }
  }
}
