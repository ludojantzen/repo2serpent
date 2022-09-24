/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : branchfrac.c                                   */
/*                                                                           */
/* Created:       2016/08/12 (JLe)                                           */
/* Last modified: 2016/08/12 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Returns isomeric branching fraction to ground state          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BranchFrac:"

/*****************************************************************************/

double BranchFrac(long rea, long RFS, double E, long id)
{
  double f, F0, F1, F;
  long ptr, erg, i, ne;

  /* Check reaction pointer and state number */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  CheckValue(FUNCTION_NAME, "RFS", "", RFS, 0, 1);

  /* Get pointer to branch fraction data */

  if ((ptr = (long)RDB[rea + REACTION_PTR_ISO_BRA]) < VALID_PTR)
    {
      /* Check state */

      if (RFS == 0)      
        return 1.0;
      else
        return 0.0;
    }

  /* Check if fixed value is given */

  if ((F = RDB[ptr + BRA_LIST_FIX_FRAC]) >= 0.0)
    {
      /* Check state */

      if (RFS == 0)      
        return F;
      else
        return 1.0 - F;
    }

  /* Pointer to energy array */

  erg = (long)RDB[ptr + BRA_LIST_PTR_E];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Get interpolation factor */

  if ((f = GridFactor(erg, E, id)) < 0.0)
    {
      /* Check state */

      if (RFS == 0)      
        return 1.0;
      else
        return 0.0;
    }

  /* Get number of energy points */

  ne = (long)RDB[ptr + BRA_LIST_NP];
  CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000000);

  /* Separate integer and decimal parts of interpolation factor */

  i = (long)f;
  f = f - (double)i;

  /* Check */

  CheckValue(FUNCTION_NAME, "i", "", i, 0, ne - 1);

  /* Pointer to factor */

  ptr = (long)RDB[ptr + BRA_LIST_PTR_F];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get tabulated values */
          
  F0 = RDB[ptr + i];
  CheckValue(FUNCTION_NAME, "F", "", F0, 0.0, 1.0);

  F1 = RDB[ptr + i + 1];
  CheckValue(FUNCTION_NAME, "F", "", F1, 0.0, 1.0);
              
  /* Interpolate */
          
  if (i == ne - 1)
    F = (1.0 - f)*F0;
  else
    F = f*(F1 - F0) + F0;

  /* Check state */

  if (RFS != 0)
    F = 1.0 - F;

  /* Check value */

  CheckValue(FUNCTION_NAME, "F", "", F, 0.0, 1.0);

  /* Return value */

  return F;
}

/*****************************************************************************/
