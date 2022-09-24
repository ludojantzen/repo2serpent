/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printcoevals.c                                 */
/*                                                                           */
/* Created:       2016/02/22 (JLe)                                           */
/* Last modified: 2016/02/22 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Prints results in .coe output                                */
/*                                                                           */
/* Comments: - Errors are not printed at the moment                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintCoeVals:"

/*****************************************************************************/

void PrintCoeVals(FILE *fp, long loc0)
{
  long dim, ptr, nv, n, N[4], i0, i1, i2, i3;
  double val, err;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  
  /* Get dimension */

  if ((dim = (long)RDB[loc0 + SCORE_DIM]) > 4)
    Die(FUNCTION_NAME, "Increase array size");
  
  CheckValue(FUNCTION_NAME, "dim", "", dim, 1, 4);

  /* Pointer to list of maximum values */
  
  ptr = (long)RDB[loc0 + SCORE_PTR_NMAX];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Reset dimension vector */

  for (n = 0; n < 4; n++)
    N[n] = 1;

  /* Loop over dimensions, get sizes and number of values */

  nv = 1;
  for (n = 0; n < dim; n++)
    {
      N[n] = (long)RDB[ptr++];
      nv = nv*N[n];
    }

  /* Print variable name */

  fprintf(fp, "%s %ld", GetText(loc0 + SCORE_PTR_NAME), nv);

  /* Avoid compiler warning */

  val = -1.0;
  err = -1.0;

  /* Loop over dimensions */

  for (i0 = 0; i0 < N[0]; i0++)
  for (i1 = 0; i1 < N[1]; i1++)
  for (i2 = 0; i2 < N[2]; i2++)
  for (i3 = 0; i3 < N[3]; i3++)
    {
      if (dim == 1)
        {
          val = Mean(loc0, i0);
          err = RelErr(loc0, i0);
        }
      else if (dim == 2)
        {
          val = Mean(loc0, i0, i1);
          err = RelErr(loc0, i0, i1);
        }
      else if (dim == 3)
        {
          val = Mean(loc0, i0, i1, i2);
          err = RelErr(loc0, i0, i1, i2);
        }
      else if (dim == 4)
        {
          val = Mean(loc0, i0, i1, i2, i3);
          err = RelErr(loc0, i0, i1, i2, i3);
        }

      /* Print */

      if ((long)RDB[DATA_COEF_CALC_INCLUDE_ERRORS] == (double)NO)
        fprintf(fp, " %12.5E", val);
      else if ((long)RDB[DATA_COEF_CALC_INCLUDE_ERRORS] == (double)YES)
        {
          /* Check magnitude of errors */

          if (err < 0.00010)
            fprintf(fp, " %12.5E %7.1E", val, err);
          else
            fprintf(fp, " %12.5E %7.5f", val, err);
        }
      else
        Die(FUNCTION_NAME, "Invalid mode");
    }

  /* Newline */

  fprintf(fp, "\n");
}

/*****************************************************************************/
