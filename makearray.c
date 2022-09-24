/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : makearray.c                                    */
/*                                                                           */
/* Created:       2010/12/15 (JLe)                                           */
/* Last modified: 2015/05/03 (TKa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates a linear- or log-interpolated array of doubles       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "MakeArray:"

/*****************************************************************************/

double *MakeArray(double min, double max, long N, long mode)
{
  static double *val;
  long n;

  /* Allocate memory for data */

  val = (double *)Mem(MEM_ALLOC, N, sizeof(double));

  /* Check mode */

  if (mode == 1)
    {
      /* Linear interpolation */
      
      for (n = 0; n < N; n++)
        val[n] = ((double)n)/((double)N - 1.0)*(max - min) + min;
    }
  else if (mode == 2)
    {
      /* Log interpolation */

      for (n = 0; n < N; n++)
        val[n] = min*pow(max/min, (double)n/((double)N - 1.0));
    }
  else
    Die(FUNCTION_NAME, "Invalid mode %ld", mode);

  /* Return pointer to data */

  return val;
}

/*****************************************************************************/
