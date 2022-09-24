/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printmvar.c                                    */
/*                                                                           */
/* Created:       2014/06/27 (JLe)                                           */
/* Last modified: 2014/08/15 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Prints results in Matlab format variables                    */
/*                                                                           */
/* Comments: - Used for printing Serpent 2 style results                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintMVar:"

/*****************************************************************************/

void PrintMVar(FILE *fp, long param)
{
  long loc0, dim, ptr, nv, n, N[4], i0, i1, i2, i3;
  char tmpstr[MAX_STR];
  double val, err;

  /* Get pointer */

  loc0 = (long)RDB[param];
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

  sprintf(tmpstr, "%-26s(idx, [1: %3ld])", GetText(loc0 + SCORE_PTR_NAME),2*nv);
  fprintf(fp, "%-40s = [", tmpstr);

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

      if (err < 0.00010)
        fprintf(fp, " %12.5E %7.1E", val, err);
      else
        fprintf(fp, " %12.5E %7.5f", val, err);
    }

  fprintf(fp, " ];\n");
}

/*****************************************************************************/
