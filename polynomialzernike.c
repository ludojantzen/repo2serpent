/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : polynomialzernike.c                            */
/*                                                                           */
/* Created:       2017/02/20 (BWe)                                           */
/* Last modified: 2018/02/23 (BWe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Evaluate the Zernike polynomial                              */
/*                                                                           */
/* Comments: - Returns the Zernike polynomial for (n, m) evaluated at rho    */
/*             and phi. This form is normalized so that all evaluations are  */
/*             within the range [-1, 1].                                     */
/*                                                                           */
/*           - Merged several input files as local functions (JLe)           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PolynomialZernike:"

#define MAX_DIRECT_CALCULATION  10

#define Save(j, value)  (PutPrivateData(vectorPtr + j, value, threadID))
#define Load(j)         (GetPrivateData(vectorPtr + j, threadID))

/* Local prototypes */

long ConvertIndicesZernikeDoubleToSingle(long, long);

/*****************************************************************************/

void PolynomialZernike(long order, double rho, double phi, long vectorPtr, long threadID, long type)
{
  long j, n, m, q, a;
  double H1, H2, H3;
  const double rho2 = rho * rho;
  const double rho4 = rho2 * rho2;

#ifdef DEBUG
  /* Ensure order >= 0 */

  if (order < 0)
    Die(FUNCTION_NAME, "Negative value for Zernike order n %ld\n", order);
#endif /* DEBUG */

  switch (type)
  {
    case FET_CALCULATION_STANDARD:
    {
      if (rho <= ZERO)
      {
        /* An edge case on the recurrence relations will blow up if rho is very */
        /* small or zero. Evaluate using known formulas for rho == 0            */

        /* Set all array values to 0.0                                          */

        for (n = 0; n <= ConvertIndicesZernikeDoubleToSingle(order, order); ++n)
          Save(n, 0.0);

        for (n = 0; n <= order; n += 2)
        {
          j = ConvertIndicesZernikeDoubleToSingle(n, 0);

          if ((n / 2) % 2 != 0)
            Save(j, -1);
          else
            Save(j, 1);
        }

        return;
      }

      /* Calculate all the m >= 0 radial polynomials first */

      switch (order)
      {
        default:
        case MAX_DIRECT_CALCULATION:  /* 10 */
          /* n = 10, m = 10 */
          Save(65, rho4 * rho4 * rho2);
          /* n = 10, m = 8 */
          Save(64, (10 * rho2 - 9) * rho4 * rho4);
          /* n = 10, m = 6 */
          Save(63, ((45 * rho2 - 72) * rho2 + 28) * rho4 * rho2);
          /* n = 10, m = 4 */
          Save(62, (((120 * rho2 - 252) * rho2 + 168) * rho2 - 35) * rho4);
          /* n = 10, m = 2 */
          Save(61, ((((210 * rho2 - 504) * rho2 + 420) * rho2 - 140) * rho2 + 15) * rho2);
          /* n = 10, m = 0 */
          Save(60, (((((252 * rho2 - 630) * rho2 + 560) * rho2 - 210) * rho2 + 30) * rho2 - 1));

        case 9:
          /* n = 9, m = 9 */
          Save(54, rho4 * rho4 * rho);
          /* n = 9, m = 7 */
          Save(53, (9 * rho2 - 8) * rho4 * rho2 * rho);
          /* n = 9, m = 5 */
          Save(52, ((36 * rho2 - 56) * rho2 + 21) * rho4 * rho);
          /* n = 9, m = 3 */
          Save(51, (((84 * rho2 - 168) * rho2 + 105) * rho2 - 20) * rho2 * rho);
          /* n = 9, m = 1 */
          Save(50, ((((126 * rho2 - 280) * rho2 + 210) * rho2 - 60) * rho2 + 5) * rho);

        case 8:
          /* n = 8, m = 8 */
          Save(44, rho4 * rho4);
          /* n = 8, m = 6 */
          Save(43, (8 * rho2 - 7) * rho4 * rho2);
          /* n = 8, m = 4 */
          Save(42, ((28 * rho2 - 42) * rho2 + 15) * rho4);
          /* n = 8, m = 2 */
          Save(41, (((56 * rho2 - 105) * rho2 + 60) * rho2 - 10) * rho2);
          /* n = 8, m = 0 */
          Save(40, ((((70 * rho2 - 140) * rho2 + 90) * rho2 - 20) * rho2 + 1));

        case 7:
          /* n = 7, m = 7 */
          Save(35, rho4 * rho2 * rho);
          /* n = 7, m = 5 */
          Save(34, (7 * rho2 - 6) * rho4 * rho);
          /* n = 7, m = 3 */
          Save(33, ((21 * rho2 - 30) * rho2 + 10) * rho2 * rho);
          /* n = 7, m = 1 */
          Save(32, (((35 * rho2 - 60) * rho2 + 30) * rho2 - 4) * rho);

        case 6:
          /* n = 6, m = 6 */
          Save(27, rho4 * rho2);
          /* n = 6, m = 4 */
          Save(26, (6 * rho2 - 5) * rho4);
          /* n = 6, m = 2 */
          Save(25, ((15 * rho2 - 20) * rho2 + 6) * rho2);
          /* n = 6, m = 0 */
          Save(24, (((20 * rho2 - 30) * rho2 + 12) * rho2 - 1));

        case 5:
          /* n = 5, m = 5 */
          Save(20, rho4 * rho);
          /* n = 5, m = 3 */
          Save(19, (5 * rho2 - 4) * rho2 * rho);
          /* n = 5, m = 1 */
          Save(18, ((10 * rho2 - 12) * rho2 + 3) * rho);

        case 4:
          /* n = 4, m = 4 */
          Save(14, rho4);
          /* n = 4, m = 2 */
          Save(13, (4 * rho2 - 3) * rho2);
          /* n = 4, m = 0 */
          Save(12, ((6 * rho2 - 6) * rho2 + 1));

        case 3:
          /* n = 3, m = 3 */
          Save(9, rho2 * rho);
          /* n = 3, m = 1 */
          Save(8, (3 * rho2 - 2) * rho);

        case 2:
          /* n = 2, m = 2 */
          Save(5, rho2);
          /* n = 3, m = 0 */
          Save(4, (2 * rho2 - 1));

        case 1:
          /* n = 1, m = 1 */
          Save(2, rho);

        case 0:
          /* n = 0, m = 0 */
          Save(0, 1);
      }

      for (n = MAX_DIRECT_CALCULATION + 1; n <= order; ++n)
      {
        /* Calculate for m == n */

        j = ConvertIndicesZernikeDoubleToSingle(n, n);
        Save(j, pow(rho, n)
                * (n + n + 2) / M_PI);

        /* Calculate for m = n - 2 */

        j--; /* j = (n, n - 2) */
        Save(j, n * Load(j + 1) - (n + 1) * Load(j - (n + n)));
        /*               └───┴(n, n)             └────────┴(n - 2, n - 2) */

        /* Calculate any remaining polynomials */

        for (q = n; q >= 4; q -= 2)
        {
          H3 = (-4 * (q - 2) * (q - 3)) / (double)((n + q - 2) * (n - q + 4));
          H2 = (H3 * (n + q) * (n - q + 2)) / (4 * (q - 1)) + (q - 2);
          H1 = q * (q - 1) / 2 - q * H2 + (H3 * (n + q + 2) * (n - q)) / 8;
          j--; /* j = (n, q - 4) */
          Save(j,  H1 * Load(j + 2) + (H2 + H3 / rho2) * Load(j + 1));
          /*                 └───┴(n, q)                      └───┴(n, q - 2) */
        }
      }

      break;
    }

    case FET_CALCULATION_ORTHONORMAL:
    {
      if (rho <= ZERO)
      {
        /* An edge case on the recurrence relations will blow up if rho is very */
        /* small or zero. Evaluate using known formulas for rho == 0            */

        /* Set all array values to 0.0                                          */

        for (n = 0; n <= ConvertIndicesZernikeDoubleToSingle(order, order); ++n)
          Save(n, 0.0);

        for (n = 0; n <= order; n += 2)
        {
          j = ConvertIndicesZernikeDoubleToSingle(n, 0);

          if ((n / 2) % 2 != 0)
            Save(j, -1
                    * (n + 1) / M_PI);
          else
            Save(j, 1
                    * (n + 1) / M_PI);
        }

        return;
      }

      /* Calculate all the m >= 0 radial polynomials first */

      switch (order)
      {
        default:
        case MAX_DIRECT_CALCULATION:  /* 10 */
          /* n = 10, m = 10 */
          Save(65, rho4 * rho4 * rho2
                   * 22 / M_PI);
          /* n = 10, m = 8 */
          Save(64, (10 * rho2 - 9) * rho4 * rho4
                   * 22 / M_PI);
          /* n = 10, m = 6 */
          Save(63, ((45 * rho2 - 72) * rho2 + 28) * rho4 * rho2
                   * 22 / M_PI);
          /* n = 10, m = 4 */
          Save(62, (((120 * rho2 - 252) * rho2 + 168) * rho2 - 35) * rho4
                   * 22 / M_PI);
          /* n = 10, m = 2 */
          Save(61, ((((210 * rho2 - 504) * rho2 + 420) * rho2 - 140) * rho2 + 15) * rho2
                   * 22 / M_PI);
          /* n = 10, m = 0 */
          Save(60, (((((252 * rho2 - 630) * rho2 + 560) * rho2 - 210) * rho2 + 30) * rho2 - 1)
                   * 11 / M_PI);

        case 9:
          /* n = 9, m = 9 */
          Save(54, rho4 * rho4 * rho
                   * 20 / M_PI);
          /* n = 9, m = 7 */
          Save(53, (9 * rho2 - 8) * rho4 * rho2 * rho
                   * 20 / M_PI);
          /* n = 9, m = 5 */
          Save(52, ((36 * rho2 - 56) * rho2 + 21) * rho4 * rho
                   * 20 / M_PI);
          /* n = 9, m = 3 */
          Save(51, (((84 * rho2 - 168) * rho2 + 105) * rho2 - 20) * rho2 * rho
                   * 20 / M_PI);
          /* n = 9, m = 1 */
          Save(50, ((((126 * rho2 - 280) * rho2 + 210) * rho2 - 60) * rho2 + 5) * rho
                   * 20 / M_PI);

        case 8:
          /* n = 8, m = 8 */
          Save(44, rho4 * rho4
                   * 18 / M_PI);
          /* n = 8, m = 6 */
          Save(43, (8 * rho2 - 7) * rho4 * rho2
                   * 18 / M_PI);
          /* n = 8, m = 4 */
          Save(42, ((28 * rho2 - 42) * rho2 + 15) * rho4
                   * 18 / M_PI);
          /* n = 8, m = 2 */
          Save(41, (((56 * rho2 - 105) * rho2 + 60) * rho2 - 10) * rho2
                   * 18 / M_PI);
          /* n = 8, m = 0 */
          Save(40, ((((70 * rho2 - 140) * rho2 + 90) * rho2 - 20) * rho2 + 1)
                   * 9 / M_PI);

        case 7:
          /* n = 7, m = 7 */
          Save(35, rho4 * rho2 * rho
                   * 16 / M_PI);
          /* n = 7, m = 5 */
          Save(34, (7 * rho2 - 6) * rho4 * rho
                   * 16 / M_PI);
          /* n = 7, m = 3 */
          Save(33, ((21 * rho2 - 30) * rho2 + 10) * rho2 * rho
                   * 16 / M_PI);
          /* n = 7, m = 1 */
          Save(32, (((35 * rho2 - 60) * rho2 + 30) * rho2 - 4) * rho
                   * 16 / M_PI);

        case 6:
          /* n = 6, m = 6 */
          Save(27, rho4 * rho2
                   * 14 / M_PI);
          /* n = 6, m = 4 */
          Save(26, (6 * rho2 - 5) * rho4
                   * 14 / M_PI);
          /* n = 6, m = 2 */
          Save(25, ((15 * rho2 - 20) * rho2 + 6) * rho2
                   * 14 / M_PI);
          /* n = 6, m = 0 */
          Save(24, (((20 * rho2 - 30) * rho2 + 12) * rho2 - 1)
                   * 7 / M_PI);

        case 5:
          /* n = 5, m = 5 */
          Save(20, rho4 * rho
                   * 12 / M_PI);
          /* n = 5, m = 3 */
          Save(19, (5 * rho2 - 4) * rho2 * rho
                   * 12 / M_PI);
          /* n = 5, m = 1 */
          Save(18, ((10 * rho2 - 12) * rho2 + 3) * rho
                   * 12 / M_PI);

        case 4:
          /* n = 4, m = 4 */
          Save(14, rho4
                   * 10 / M_PI);
          /* n = 4, m = 2 */
          Save(13, (4 * rho2 - 3) * rho2
                   * 10 / M_PI);
          /* n = 4, m = 0 */
          Save(12, ((6 * rho2 - 6) * rho2 + 1)
                   * 5 / M_PI);

        case 3:
          /* n = 3, m = 3 */
          Save(9, rho2 * rho
                  * 8 / M_PI);
          /* n = 3, m = 1 */
          Save(8, (3 * rho2 - 2) * rho
                  * 8 / M_PI);

        case 2:
          /* n = 2, m = 2 */
          Save(5, rho2
                  * 6 / M_PI);
          /* n = 3, m = 0 */
          Save(4, (2 * rho2 - 1)
                  * 3 / M_PI);

        case 1:
          /* n = 1, m = 1 */
          Save(2, rho
                  * 4 / M_PI);

        case 0:
          /* n = 0, m = 0 */
          Save(0, 1
                  * 1 / M_PI);
      }

      for (n = MAX_DIRECT_CALCULATION + 1; n <= order; ++n)
      {
        /* Calculate for m == n */

        j = ConvertIndicesZernikeDoubleToSingle(n, n);
        Save(j, pow(rho, n)
                * (n + n + 2) / M_PI);

        /* Calculate for m = n - 2 */

        j--; /* j = (n, n - 2) */
        Save(j, n * Load(j + 1) - (n + 1) * Load(j - (n + n)));
        /*               └───┴(n, n)             └────────┴(n - 2, n - 2) */

        /* Calculate any remaining polynomials */

        for (q = n; q >= 4; q -= 2)
        {
          H3 = (-4 * (q - 2) * (q - 3)) / (double)((n + q - 2) * (n - q + 4));
          H2 = (H3 * (n + q) * (n - q + 2)) / (4 * (q - 1)) + (q - 2);
          H1 = q * (q - 1) / 2 - q * H2 + (H3 * (n + q + 2) * (n - q)) / 8;
          j--; /* j = (n, q - 4) */
          if (q == 4) /* m == 0 */
            Save(j, (H1 * Load(j + 2) + (H2 + H3 / rho2) * Load(j + 1))
                    * 0.5);
          else
            Save(j,  H1 * Load(j + 2) + (H2 + H3 / rho2) * Load(j + 1));
          /*                   └───┴(n, q)                      └───┴(n, q - 2) */
        }
      }

      break;
    }
  }

  /* Fill out negative ranks and apply the azimuthal components.  */
  /* Nothing needs to be done for the ranks m = 0.                */

  j = 0;
  for (n = 1; n <= order; ++n)
  {
    /* Work from the outside edges of the triangle inward */

    j += n; /* j = (n, -m) */
    for (m = 0, q = a = n; m < q; ++m, --q, a -= 2)
    {
      Save(j + m, Load(j + q) * sin(a * phi));
      /*   └───┴(n, -m)└───┴(n, m) */
      Save(j + q, Load(j + q) * cos(a * phi));
      /*   └───┴(n, m) └───┴(n, m) */
    }
  }

#ifdef TEST_TO_ENSURE_THAT_EVERYTHING_IS_COMPUTED_PROPERLY
  for (n = 0, j = 0; n <= order; ++n)
  {
    for (m = 0; m <= order * 2 + 1; ++m)
      if (m > order - n && m < order + n && m % 2 != n % 2)
      {
        printf(" %+6.6f", Load(j++));
      }
      else
        printf("          ");

    printf("\n");
  }
  Die(FUNCTION_NAME, "Zernike test");
#endif
}
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/* Description: Converts the double Zernike polynomial indices (n, m) to the */
/*              single index (j)                                             */
/*                                                                           */
/* Comments: Algorithm from: http://dx.doi.org/10.1080/09500340.2011.554896  */
/*                                                                           */
/*****************************************************************************/

long ConvertIndicesZernikeDoubleToSingle(long n, long m)
{
#ifdef DEBUG
  /* Ensure n >= 0 */

  if (n < 0)
    Die(FUNCTION_NAME, "Negative value for Zernike index n %ld\n", n);

  /* Ensure |m| <= n */

  if (labs(m) > n)
    Die(FUNCTION_NAME, "Zernike index m (%ld) further from zero than index n (%ld)\n", m, n);
#endif /* DEBUG */

  /* Convert (n, m) indices to the index j */

  return (n * (n + 2) + m) / 2;
}
