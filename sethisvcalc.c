/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sethisvcalc.c                                  */
/*                                                                           */
/* Created:       2020/03/05 (JLe)                                           */
/* Last modified: 2020/03/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Sets historiy variations                                     */
/*                                                                           */
/* Comments: - Called right after ReadInput() to invoke changes before       */
/*             anything else is done.                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetHisvCalc:"

/*****************************************************************************/

long SetHisvCalc(long idx, long cyc)
{
  long loc0, loc1, loc2, bra, dep, n, j, ptr, ns;

  /* Pointer to history variations */

  if ((loc0 = (long)RDB[DATA_PTR_HISV0]) < VALID_PTR)
    return idx;

  /* Set restart file writing on */

  WDB[DATA_WRITE_RESTART_FILE] = (double)YES;

  /* Check if first pass */

  if (idx == -1)
    {
      /* Reset restart file reading */

      WDB[DATA_READ_RESTART_FILE] = (double)NO;
      WDB[DATA_RESTART_READ_CONTINUE] = (double)NO;

      /* Pointer to variations */

      loc1 = (long)RDB[loc0 + HISV_PTR_VAR];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Set break point */

      WDB[DATA_HISV_BREAK_POINT] = RDB[loc1 + HISV_VAR_BU];

      /* Set continue flag */

      WDB[DATA_MORE_COEF_CALC] = (double)YES;

      /* Exit subroutine */

      return idx + 1;
    }

  /* Print */

  fprintf(outp, "Re-configuring model for history variation...\n");

  /* Set restart file reading */

  WDB[DATA_READ_RESTART_FILE] = (double)YES;
  WDB[DATA_RESTART_READ_CONTINUE] = (double)YES;

  /* Reset count */

  j = 0;

  /* Loop over break points points */

  loc1 = (long)RDB[loc0 + HISV_PTR_VAR];
  while (loc1 > VALID_PTR)
    {
      /* Number of variations */

      ns = (long)RDB[loc1 + HISV_VAR_N_VAR];
      CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 100000);

      /* Pointer to data */

      loc2 = (long)RDB[loc1 + HISV_VAR_PTR_VAR];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Loop over variations */

      for (n = 0; n < ns; n++)
        {
          /* Find match */

          bra = (long)RDB[DATA_PTR_BRA0];
          while (bra > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(bra + DEP_BRA_PTR_NAME, loc2 + n))
                break;

              /* Next */

              bra = NextItem(bra);
            }

          /* Check */

          if (bra < VALID_PTR)
            Die(FUNCTION_NAME, "Branch is not defined");

          /* Invoke branch */

          InvokeBranch(bra);
        }

      /* Check count */

      if (j++ == idx)
        break;

      /* Next */

      loc1 = NextItem(loc1);
    }

  /* Check */

  if (loc1 < VALID_PTR)
    Die(FUNCTION_NAME, "WTF?");

  /* Avoid compiler warning */

  n = -1;
  ns = -1;

  /* Reset count */

  j = 0;

  /* Re-configure depletion history */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Get number of steps */

      ns = (long)RDB[dep + DEP_HIS_N_STEPS];
      CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 1000000);

      /* Pointer to steps */

      ptr = (long)RDB[dep + DEP_HIS_PTR_STEPS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over steps and find match */

      for (n = 0; n < ns; n++)
        {
          /* Check */

          if (j++ == cyc)
            break;

          /* Adjust step counter */

          WDB[dep + DEP_HIS_N_STEPS] = RDB[dep + DEP_HIS_N_STEPS] - 1.0;
        }

      /* Check */

      if (n < ns)
        {
          /* Put pointer to next step */

          WDB[dep + DEP_HIS_PTR_STEPS] = (double)(ptr + n);

          /* Put pointer to depletion history */

          WDB[DATA_BURN_PTR_DEP] = (double)dep;

          /* Break loop */

          break;
        }

      /* Next */

      dep = NextItem(dep);
    }

  /* Put burnup step */

  WDB[DATA_BURN_STEP] = (double)cyc;

  /* Check */

  if ((n == ns) && (NextItem(loc1) > VALID_PTR))
    Error(loc0, "Error in break points");

  /* Check if more break points */

  if ((loc1 = NextItem(loc1)) > VALID_PTR)
    {
      /* Set break point */

      WDB[DATA_HISV_BREAK_POINT] = RDB[loc1 + HISV_VAR_BU];

      /* Set continue flag */

      WDB[DATA_MORE_COEF_CALC] = (double)YES;
    }
  else
    {
      /* Reset continue flag */

      WDB[DATA_MORE_COEF_CALC] = (double)NO;
    }

  /* OK */

  fprintf(outp, "OK.\n\n");

  /* Return updated index */

  return idx + 1;
}

/*****************************************************************************/
