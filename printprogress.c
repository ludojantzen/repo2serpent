/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printprogress.c                                */
/*                                                                           */
/* Created:       2012/12/30 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Prints completed fraction to output in routines handling     */
/*              materials.                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintProgress:"

/*****************************************************************************/

void PrintProgress(long mat, long type)
{
  double frac;

  /* Check OpenMP thread number and MPI task id */

  if ((OMP_THREAD_NUM > 0) || (mpiid > 0))
    return;

  /* Special */

  if (type == 0)
    {
      /* Print */

      fprintf(outp, "   0%% complete\n");

      /* Reset previous */

      WDB[DATA_PRINT_PREV_COMPLETE] = 0.0;

      /* Exit subroutine */

      return;
    }
  else if (type == 100)
    {
      /* Print */
      
      fprintf(outp, " 100%% complete\n\n");

      /* Exit subroutine */

      return;
    }

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Avoid compiler warning */

  frac = -1.0;

  /* Calculate fraction */

  if (type == 1)
    frac = 100.0*RDB[mat + MATERIAL_PROC_IDX]/RDB[DATA_N_MATERIALS];
  else if (type == 2)
    frac = 100.0*RDB[mat + MATERIAL_BURN_IDX]/RDB[DATA_N_BURN_MATERIALS];
  else
    Die(FUNCTION_NAME, "Invalid print type");

  /* Compare fraction to previously completed */

  if (frac - RDB[DATA_PRINT_PREV_COMPLETE] > 5.0)
    {
      /* Print */

      if ((frac > 4.5) && (frac < 95.5))
        fprintf(outp, " %3.0f%% complete\n",  frac);

      /* Remember last value */

      WDB[DATA_PRINT_PREV_COMPLETE] = frac;
    }
}

/*****************************************************************************/
