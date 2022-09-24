/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : adjustenergygrid.c                             */
/*                                                                           */
/* Created:       2010/12/09 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Reconfigures energu grid structure                           */
/*                                                                           */
/* Comments: - Tää muuttaa gridin samalla yhdelle tasolle                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AdjustEnergyGrid:"

/*****************************************************************************/

void AdjustEnergyGrid(long erg, long ne, const double *E)
{
  long ptr, n;
  double Emin, Emax;

  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Get pointer to data */

  ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Compare size to allocated and clear data */

  if (ne > (long)RDB[erg + ENERGY_GRID_ALLOC_NE])
    Die(FUNCTION_NAME, "New grid size exceeds maximum");
  else
    memset(&WDB[ptr], 0.0, 
	   (long)RDB[erg + ENERGY_GRID_ALLOC_NE]*sizeof(double));

  /* Check ascending order */

  for (n = 1; n < ne; n++)
    if (E[n] <= E[n - 1])
      Die(FUNCTION_NAME, "Energy array is not in ascending order");

  /* Put values */

  for (n = 0; n < ne; n++)
    WDB[ptr + n] = E[n];

  /* Adjust grid size */

  WDB[erg + ENERGY_GRID_NE] = (double)ne;

  /* Get minimum and maximum energy */

  Emin = E[0];
  Emax = E[ne - 1];

  /* Store values */

  WDB[erg + ENERGY_GRID_EMIN] = Emin;
  WDB[erg + ENERGY_GRID_EMAX] = Emax;

  /* Store log values */

  WDB[erg + ENERGY_GRID_LOG_EMIN] = log(Emin);
  WDB[erg + ENERGY_GRID_LOG_EMAX] = log(Emax);

  /* Reset pointers */

  WDB[erg + ENERGY_GRID_PTR_LOW] = NULLPTR;
  WDB[erg + ENERGY_GRID_PTR_HIGH] = NULLPTR;
  WDB[erg + ENERGY_GRID_PTR_BINS] = NULLPTR;

  /* Reset mid-point energy and number of bins */

  WDB[erg + ENERGY_GRID_EMID] = -1.0;
  WDB[erg + ENERGY_GRID_NB] = -1.0;
}

/*****************************************************************************/
