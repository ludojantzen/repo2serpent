/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gridsearch.c                                   */
/*                                                                           */
/* Created:       2010/12/09 (JLe)                                           */
/* Last modified: 2018/01/18 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Finds energy interval from grid structure                    */
/*                                                                           */
/* Comments: - Yläraja otetaan väkisin mukaan vaikka sen pitäisi             */
/*             periaatteessa olla ulkopuolella.                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GridSearch:"

/*****************************************************************************/

long GridSearch(long erg, double E)
{
  double Emin, Emax, Emid, logE;
  long ptr, ne, idx, i0, nb, erg0, id;

  /* Check pointer and value */

  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Get minimum and maximum energy */

  Emin = RDB[erg + ENERGY_GRID_EMIN];
  CheckValue(FUNCTION_NAME, "Emin", "", Emin, 0.0, INFTY);

  Emax = RDB[erg + ENERGY_GRID_EMAX];
  CheckValue(FUNCTION_NAME, "Emax", "", Emax, 0.0, INFTY);

  /* Compare to limits */

  if ((E < Emin) || (E > Emax))
    return -1;
    
  /***************************************************************************/

  /***** Double-indexing *****************************************************/

  if ((ptr = (long)RDB[erg + ENERGY_GRID_PTR_DIX_IDX]) > VALID_PTR)
    {
      /* Pointer to unionized energy grid */

      erg0 = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(erg0)", DATA_ARRAY, erg0);

      /* Get OpenMP thread number */

      id = OMP_THREAD_NUM;

      /* Number of bins */

      ne = (long)RDB[erg + ENERGY_GRID_NE];
      
      /* Get index from grid factor (this allows access to previously */
      /* searched value) */

      if ((idx = (long)GridFactor(erg0, E, id)) < 0)
        Die(FUNCTION_NAME, "Unionized grid search failed");

      /* Get index from table */

      idx = (long)RDB[ptr + idx];
      CheckValue(FUNCTION_NAME, "idx", "", idx, 0, ne - 1);

      /* Check energy */

      ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckValue(FUNCTION_NAME, "E", "", E, RDB[ptr + idx], RDB[ptr + idx + 1]);

      /* Return index */

      return idx;
    }

  /***************************************************************************/

  /***** Sub-intervals *******************************************************/
  
  /* Get number of bins */
  
  while ((nb = (long)RDB[erg + ENERGY_GRID_NB]) > 0)
    {

#ifdef DEBUG

      if (RDB[erg + ENERGY_GRID_EMID] > 0.0)
        Die(FUNCTION_NAME, "Mid energy with bins");

#endif

      /* Check grid type */

      if ((long)RDB[erg + ENERGY_GRID_TYPE] == GRID_TYPE_LOG)
        {
          /* Take log of E */

          logE = log(E);
  
          /* Get log of minimum and maximum energy */
          
          Emin = RDB[erg + ENERGY_GRID_LOG_EMIN];
          Emax = RDB[erg + ENERGY_GRID_LOG_EMAX];
          
          /* Check energy */

          CheckValue(FUNCTION_NAME, "logE", "", logE, Emin, Emax);
          
          /* Calculate bin index */
          
          idx = (long)((double)nb*(logE - Emin)/(Emax - Emin));
        }
      else
        {
          /* Get minimum and maximum energy */
          
          Emin = RDB[erg + ENERGY_GRID_EMIN];
          Emax = RDB[erg + ENERGY_GRID_EMAX];
          
          /* Check energy */
          
          CheckValue(FUNCTION_NAME, "E", "", E, Emin, Emax);
          
          /* Calculate bin index */
          
          idx = (long)((double)nb*(E - Emin)/(Emax - Emin));
        }

      /* Check value */

      CheckValue(FUNCTION_NAME, "idx1", "", idx, -1, nb);

      /* Log-function may cause numerical problems */
      
      if (idx < 0)
        idx = 0;
      else if (idx > nb - 1)
        idx = nb - 1;

      /* Pointer to bins */
      
      ptr = (long)RDB[erg + ENERGY_GRID_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);      
      
      /* pointer to grid */
      
      erg = (long)RDB[ptr + idx];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);      
    }
  
  /***************************************************************************/

  /***** Binary-tree search **************************************************/
      
  while ((Emid = RDB[erg + ENERGY_GRID_EMID]) > 0.0)
    {
      if (E < Emid)
        erg = (long)RDB[erg + ENERGY_GRID_PTR_LOW];
      else
        erg = (long)RDB[erg + ENERGY_GRID_PTR_HIGH];
      
      /* Check pointer */
      
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
    }
  
  /***************************************************************************/

  /***** Find index in interval **********************************************/

  /* Check limits (JLe / 19.12.2017 / 2.1.30: Toi alaraja menee jostain   */
  /* syystä fotoneilla rikki xsplotterissa jonkun pyöristysvirheen takia) */

  CheckValue(FUNCTION_NAME, "E", "", E, 0.999999*RDB[erg + ENERGY_GRID_EMIN],
             RDB[erg + ENERGY_GRID_EMAX]);
  
  /* Get number of points */
  
  ne = (long)RDB[erg + ENERGY_GRID_NE];
  CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);

  /* Get pointer to data */

  ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get index to beginning */

  i0 = (long)RDB[erg + ENERGY_GRID_I0];
  CheckValue(FUNCTION_NAME, "i0", "", i0, 0, MAX_EGRID_NE - 2);

  /* Add value to pointer */

  ptr = ptr + i0;
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Find interval */
  
  if ((idx = SearchArray(&RDB[ptr], E, ne)) < 0)
    return -1;
        
  /* Check */

#ifdef DEBUG

  if (idx > ne - 2)
    Die(FUNCTION_NAME, "Search error 1");
  else if ((E < RDB[ptr + idx]) || (E > RDB[ptr + idx + 1]))
    Die(FUNCTION_NAME, "Search error 2 %ld %ld (%E %E %E %E)", idx, ne, 
        RDB[ptr + idx], E, RDB[ptr + idx + 1], RDB[ptr + idx + 2]);
  else if (E == RDB[ptr + idx + 1])
    Warn(FUNCTION_NAME, "Coincident energy point? (E = %E)", E);
  
#endif

  /* Return index to grid */
  
  return i0 + idx;

  /***************************************************************************/
}

/*****************************************************************************/
