/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : interpolatefisse.c                             */
/*                                                                           */
/* Created:       2017/11/24 (RTu)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates temporary fission energy deposition array on         */
/*              unionized energy grid                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InterpolateFissE:"

/*****************************************************************************/

void InterpolateFissE(double *fisse, long rea, long id)
{
  long ptr, ne, n;
  const double *E;

  /* Check mode and pointers */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO)
    Die(FUNCTION_NAME, "Macroscopic cross sections not reconstructed");
  else if (fisse == NULL)
    Die(FUNCTION_NAME, "Null pointer for fisse");
  
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to energy grid */

  ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of energy points */

  ne = (long)RDB[ptr + ENERGY_GRID_NE];

  /* Get pointer to energy array */

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  E = &RDB[ptr];

  /***************************************************************************/

  /***** Construct fission energy deposition data ****************************/

  /* Loop over energy array */

  for (n = 0; n < ne; n++)
    fisse[n] = FissE(rea, E[n], id);
}

/*****************************************************************************/
