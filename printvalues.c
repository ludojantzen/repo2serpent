/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printvalues.c                                  */
/*                                                                           */
/* Created:       2011/03/13 (JLe)                                           */
/* Last modified: 2012/04/13 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Formatted printing of mean and statistical error             */
/*                                                                           */
/* Comments: - Used for printing Serpent 1 style output                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintValues:"

/*****************************************************************************/

void PrintValues(FILE *fp, char *name, long param, long nv1, long nv2, long iu,
                 long n0, long m0)
{
  long n, m, ptr;
  char tmpstr[MAX_STR];
  double err;

  if (nv1 < 1)
    return;

  /* Get pointer to data */

  ptr = (long)RDB[param];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Print title */

  if (name != NULL)
    {
      if (nv2 < 1)
        sprintf(tmpstr, "%-26s(idx, [1: %3ld])", name, 2*nv1);
      else
        sprintf(tmpstr, "%-26s(idx, [1: %3ld])", name, 2*nv1*nv2);
    }
  else
    {
      if (nv2 < 1)
        sprintf(tmpstr, "%-26s(idx, [1: %3ld])", 
                GetText(ptr + SCORE_PTR_NAME), 2*nv1);
      else
        sprintf(tmpstr, "%-26s(idx, [1: %3ld])", 
                GetText(ptr + SCORE_PTR_NAME), 2*nv1*nv2);
    }

  fprintf(fp, "%-40s = [", tmpstr);

  /* Print values */

  if (iu < 0)
    {
      if (nv2 < 1)
        {
          /* Loop over values */
          
          for (n = 0; n < nv1; n++)
            {          
              /* Get relative statistical error */

              err = RelErr(ptr, n0 + n); 

              /* Print */

              if (err < 0.00010)
                fprintf(fp, " %12.5E %7.1E", Mean(ptr, n0 + n), err);
              else
                fprintf(fp, " %12.5E %7.5f", Mean(ptr, n0 + n), err);
            }
        }
      else
        {
          /* Loop over values */

          for (n = 0; n < nv1; n++)
            for (m = 0; m < nv2; m++)
              {
                /* Get relative statistical error */

                err = RelErr(ptr, n0 + n, m0 + m);

              /* Print */
                
                if (err < 0.00010)
                  fprintf(fp, " %12.5E %7.1E", Mean(ptr, n0 + n, m0 + m), err);
                else
                  fprintf(fp, " %12.5E %7.5f", Mean(ptr, n0 + n, m0 + m), err);
              }
        }
    }
  else
    Die(FUNCTION_NAME, "iu = %ld", iu);

  fprintf(fp, " ];\n");
}

/*****************************************************************************/
