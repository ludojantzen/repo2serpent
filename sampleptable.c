/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sampleptable.c                                 */
/*                                                                           */
/* Created:       2015/11/07 (JLe)                                           */
/* Last modified: 2018/11/05 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Samples unresolved resonance probability table               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SamplePTable:"

/*****************************************************************************/

double SamplePTable(long urs, double E, long mt, long id)
{
  long n, ptr, erg, ne, M, k0, k1, type;
  double f, r, f0, f1, rnd;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(urs)", DATA_ARRAY, urs);

  /* Check if random number is already sampled */
  
  if ((rnd = TestValuePair(urs + URES_PTR_RND, E, id)) == -INFTY)
    {
      /* Sample new value */
      
      rnd = RandF(id);
      
      /* Store value */
      
      StoreValuePair(urs + URES_PTR_RND, E, rnd, id);
    }
  
  /* Store value for checking */
  
  StoreValuePair(urs + URES_PTR_RND_CHK, E, rnd, id);
  
  /* Get pointer to energy grid */
      
  erg = (long)RDB[urs + URES_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
  
  /* Get number of grid points */
  
  ne = (long)RDB[erg + ENERGY_GRID_NE];
  
  /* Get interval */
  
  n = GridSearch(erg, E);
  CheckValue(FUNCTION_NAME, "n", "", n, 0.0, ne - 2);
      
  /* Get number of probabilities */
  
  M = (long)RDB[urs + URES_NP];
  
  /* Pointer to probability data */
  
  ptr = (long)RDB[urs + URES_PTR_PROB];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Find intervals */
  
  for (k0 = 0; k0 < M; k0++)
    if (rnd < RDB[ptr + n*M + k0])
      break;
  
  for (k1 = 0; k1 < M; k1++)
    if (rnd < RDB[ptr + (n + 1)*M + k1])
      break;
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "k0", "", rnd, 0.0, RDB[ptr + n*M + k0]);
  CheckValue(FUNCTION_NAME, "k1", "", rnd, 0.0, RDB[ptr + (n + 1)*M + k1]);
  
  /* Get pointer to energy grid data */
  
  ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA]; 
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get interpolation type */
  
  type = (long)RDB[urs + URES_INT];
  
  /* Avoid compiler warning */
  
  r = -1.0;
  
  /* Calculate interpolation factor */
  
  if (type == 2)
    r = (E - RDB[ptr + n])/(RDB[ptr + n + 1] - RDB[ptr + n]);
  else if (type == 5)
    r = log(E/RDB[ptr + n])/log(RDB[ptr + n + 1]/RDB[ptr + n]);
  else
    Die(FUNCTION_NAME, "Invalid interpolation type");
  
  /* Check factor */
  
  CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);
  
  /* Pointer to table */
  
  ptr = (long)RDB[urs + URES_PTR_FACT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get values */
  
  f0 = RDB[ptr + n*M + k0];
  f1 = RDB[ptr + (n + 1)*M + k1];
  
  /* Avoid compiler warning */

  f = -1.0;

  /* Check interpolation scheme and interpolate */
  
  if (type == 2)
    f = r*(f1 - f0) + f0;
  else if (type == 5)
    {
      /* KERMA can have negative URES factors in EDEP_MODE_NEUTRON_PHOTON */
      
      if (((long)RDB[DATA_EDEP_MODE] == EDEP_MODE_NEUTRON_PHOTON) && (mt == 301))
        Die(FUNCTION_NAME, "Invalid interpolation type for KERMA URES factors");
      
      /* Check non-zero values */
      
      if ((f0 <= 0.0) || (f1 <= 0.0))
        return 1.0;
      else
        f = exp(r*log(f1/f0) + log(f0));
    }

  /* Pointer to maximum factors */
  
  ptr = (long)RDB[urs + URES_PTR_MAXF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Check */
  
  if (((long)RDB[DATA_EDEP_MODE] != EDEP_MODE_NEUTRON_PHOTON) || (mt != 301))
    {
      CheckValue(FUNCTION_NAME, "f0", "", f0, -10.0, RDB[ptr + n]);
      CheckValue(FUNCTION_NAME, "f1", "", f1, -10.0, RDB[ptr + n]);
      CheckValue(FUNCTION_NAME, "f", "", f, -10.0, RDB[ptr + n]);
    }
  
  /* Return value */

  return f;
}

/*****************************************************************************/
