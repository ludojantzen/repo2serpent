/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : trapzreal.c                                    */
/*                                                                           */
/* Created:       2016/10/21 (TKa)                                           */
/* Last modified: 2017/03/02 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Trapezoidal integration functions.                           */
/*                                                                           */
/* Comments: - The same interpolation modes (types) are used as in           */
/*             ENDFInterp().                                                 */
/*           - Histogram integration added 3.2.2017 / 2.1.29 (JLe)           */
/*           - Lisäsin tähän option ottaa tuon koko kumulatiivisen funktion  */
/*             talteen (tarvitaan jatkuvien lähdejakaumien prosessoinnissa). */
/*             Pointterin paikalla voi välittää NULL:in.                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"


/*****************************************************************************/

void TrapzRealCum(const double *x, const double *y, double *F, long n,
                  long type) {
  /* Calculates cumulative integral F for function y(x). The first element
   * in F is set to zero.
   * */
  static char * const FUNCTION_NAME = "TrapzRealCum:";
  long i;
  double lxr, c;

  /* Check size */
  if (n < 2)
    Die(FUNCTION_NAME, "Number of elements = %ld less than two", n);

  F[0] = 0.0;

  /* Check type */
  switch (type) {
    case 1:
      {
        /* Histogram (JLe / 2.3.2017 / 2.1.29) */

        for (i = 1; i < n; i++)
          F[i] = F[i-1] + (x[i] - x[i-1])*y[i-1];

        break;
      }
    case 0:
    case 2:
      {
        /* Lin-lin */

        for (i = 1; i < n; i++) {
          F[i] = F[i-1] + 0.5*(x[i] - x[i-1])*(y[i] + y[i-1]);
        }

        break;
      }
    case 3:
      {
        /* Lin-log */

        for (i = 1; i < n; i++) {

          /* Check values */
          if (x[i-1] == x[i])
            Die(FUNCTION_NAME, "x[%ld] is equal to x[%ld] = %E (type %ld)",
                i-1, i, x[i], type);
          if (x[i-1] == 0.0)
            Die(FUNCTION_NAME, "x[%ld] is zero (type %ld)", i-1, type);

          F[i] = F[i-1] + y[i]*x[i] - y[i-1]*x[i-1]
                 + (y[i] - y[i-1])*(x[i-1] - x[i])/log(x[i]/x[i-1]);
        }

        break;
      }
    case 4:
      {
        /* Log-lin */

        for (i = 1; i < n; i++) {

          /* Check values */
          if (y[i-1] == 0.0)
            Die(FUNCTION_NAME, "y[%ld] is zero (type %ld)", i-1, type);

          if (y[i-1] == y[i])
            F[i] = F[i-1] + y[i-1]*(x[i] - x[i-1]);
          else
            F[i] = F[i-1] + (y[i] - y[i-1])*(x[i] - x[i-1])/log(y[i]/y[i-1]);
        }

        break;
      }
    case 5:
      {

        /* Log-log */

        /* Check first values */
        if (x[0] == 0.0)
          Die(FUNCTION_NAME, "x is zero (type %ld)", type);
        if (y[0] == 0.0)
          Die(FUNCTION_NAME, "y is zero (type %ld)", type);

        for (i = 1; i < n; i++) {

          /* Check values */
          if (x[i-1] == x[i])
            Die(FUNCTION_NAME, "x[%ld] is equal to x[%ld] = %E (type %ld)",
                i-1, i, x[i], type);
          if (x[i] == 0.0)
            Die(FUNCTION_NAME, "x[%ld] is zero (type %ld)", i, type);
          if (y[i] == 0.0)
            Die(FUNCTION_NAME, "y[%ld] is zero (type %ld)", i, type);

          lxr = log(x[i]/x[i-1]);
          c = log(y[i]/y[i-1])/lxr;

          if (c == -1)
            F[i] = F[i-1] + y[i-1]*x[i-1]*lxr;
          else
            F[i] = F[i-1] + (y[i]*x[i] - y[i-1]*x[i-1])/(c + 1.0);
        }

        break;
      }
    default:
      Die(FUNCTION_NAME, "Invalid interpolation mode %ld", type);
    }

}

/***************************************************************************/


/***************************************************************************/

double TrapzReal(const double *x, const double *y, long n, double *F0, 
                 long type) {
  /* Trapezoidal integration of function y(x) between x[0] and x[n-1].
   * See the mode definitions in TrapzRealCum().
   * */
  static char * const FUNCTION_NAME = "TrapzReal:";
  double Fend;
  double *F;

  /* Check size */
  if (n < 2)
    Die(FUNCTION_NAME, "Number of elements = %ld less than two", n);

  /* Allocate memory for the cumulative integral */

  if (F0 == NULL)
    F = (double *)Mem(MEM_ALLOC, n, sizeof(double));
  else
    F = F0;

  /* Calculate cumulative integral */

  TrapzRealCum(x, y, F, n, type);

  /* Store the last value */

  Fend = F[n-1];

  /* Free memory or */

  if (F0 == NULL)
    Mem(MEM_FREE, F);

  /* Return integral */

  return Fend;
}

/***************************************************************************/



