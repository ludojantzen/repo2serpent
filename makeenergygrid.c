/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : makeenergygrid.c                               */
/*                                                                           */
/* Created:       2010/12/09 (JLe)                                           */
/* Last modified: 2014/05/21 (JLe)                                           */
/* Version:       2.1.21                                                     */
/*                                                                           */
/* Description: Allocates memory and stores energy array in a search grid    */
/*              based on a binary tree                                       */
/*                                                                           */
/* Comments: - Toi rekursiivisuus ei oikein toimi pienillä grideillä         */
/*             (stack overflow) --> disabloitu toistaiseksi.                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MakeEnergyGrid:"

/*****************************************************************************/

long MakeEnergyGrid(long ne, long i0, long lvl, long ptr, 
                    const double *E, long mode)
{
  double Emin, Emax, E0, mean, med;
  long erg, idx, loc0, loc1, n, nb, imin, imax, type;

  /* Check if first level */

  if (lvl == 0)
    {
      /* Check index */

      if (i0 != 0)
        Die(FUNCTION_NAME, "i0 = %ld", i0);

      /* Check pointer */

      if (ptr > 0)
        Die(FUNCTION_NAME, "ptr = %ld", ptr);

      /* Check ascending order */

      for (n = 0; n < ne - 1; n++)
        if (E[n] > E[n + 1])
          Die(FUNCTION_NAME, "Values are not in ascending order");

      /* Allocate memory for data */

      ptr = ReallocMem(DATA_ARRAY, ne);

      /* Copy data */

      memcpy(&WDB[ptr], E, ne*sizeof(double));
    }
  else if (lvl > 10000)
    Die(FUNCTION_NAME, "Recursion error");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Allocate memory for node */

  erg = ReallocMem(DATA_ARRAY, ENERGY_GRID_BLOCK_SIZE);

  /* Put values */

  WDB[erg + ENERGY_GRID_NE] = (double)ne;
  WDB[erg + ENERGY_GRID_I0] = (double)i0;
  WDB[erg + ENERGY_GRID_PTR_DATA] = (double)ptr;

  /* Put allocated size */

  if (lvl == 0)
    WDB[erg + ENERGY_GRID_ALLOC_NE] = (double)ne;

  /* Check interpoltation mode */

  if ((mode != EG_INTERP_MODE_LIN) && (mode != EG_INTERP_MODE_LOG))
    Die(FUNCTION_NAME, "Invalid interpolation mode %ld", mode);

  /* Put mode */

  WDB[erg + ENERGY_GRID_INTERP_MODE] = (double)mode;

  /* Get minimum and maximum energy */

  Emin = E[i0];
  Emax = E[i0 + ne - 1];

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

  /* Allocate memory for previous values */

  AllocValuePair(erg + ENERGY_GRID_PTR_PREV_VAL);

  /* Calculate mean of minimum and maximum */

  mean = (Emin + Emax)/2.0;

  /* Median */

  med = E[i0 + (long)((double)ne/2.0)];

  /* Set grid type (viimeinen ehto lisätty 21.5.2014 / 2.1.21, estämään */
  /* log(0) myöhemmin. */

  if ((ne < 5) || (fabs(med/mean - 1.0) < 0.20) || (Emin == 0.0))
    type = GRID_TYPE_LIN;
  else
    type = GRID_TYPE_LOG;

  /* Put type */

  WDB[erg + ENERGY_GRID_TYPE] = (double)type;

  /* TÄÄ !!!!!!!!!!!!!!!!!!!!!! */
  /*
  return erg;
  */
  /******************************/

  /* Check if more levels are needed */
  
  if ((ne < 100) || (Emax/Emin < 10.0))
    return erg;

  /* Number of bins */

  nb = (long)((double)ne/10.0);

  /* Check type */

  if (nb < 100)
    {
      /***** Binary tree structure *******************************************/

      /* Find mid-point */

      if (type == GRID_TYPE_LIN)
        {      
          /* Linear grid */

          if ((idx = SearchArray(&E[i0], (Emin + Emax)/2.0, ne)) < 0)
            Die(FUNCTION_NAME, "idx < 0");
        }
      else
        {      
          /* Log grid */
          
          if ((idx = SearchArray(&E[i0], exp((log(Emin) + log(Emax))/2.0), 
                                 ne)) < 0)
            Die(FUNCTION_NAME, "idx < 0");
        }

      /* Check index */

      if ((idx == 0) || (idx == ne - 1))
        return erg;

      /* Put value */
      
      WDB[erg + ENERGY_GRID_EMID] = E[i0 + idx];
      
      /* Call recursively for upper and lower half */
      
      loc0 = MakeEnergyGrid(idx + 1, i0, lvl + 1, ptr, E, mode);
      WDB[erg + ENERGY_GRID_PTR_LOW] = (double)loc0;
      
      loc0 = MakeEnergyGrid(ne - idx, i0 + idx, lvl + 1, ptr, E, mode);
      WDB[erg + ENERGY_GRID_PTR_HIGH] = (double)loc0;
      
      /***********************************************************************/
    }
  else
    {
      /***** Bins ************************************************************/
      
      /* Allocate memory for bins */

      loc1 = ReallocMem(DATA_ARRAY, nb);
      
      /* Put pointer */
      
      WDB[erg + ENERGY_GRID_PTR_BINS] = (double)loc1;
      
      /* Put number of bins */
      
      WDB[erg + ENERGY_GRID_NB] = (double)nb;
      
      /* Loop over bins */
      
      imin = 0;
      
      for (n = 0; n < nb; n++)
        {
          /* Calculate energy for lower bin */

          if (type == GRID_TYPE_LIN)
            {
              /* Linear grid */
              
              E0 = ((double)n/((double)nb))*(Emax - Emin) + Emin;
            }
          else
            {
              /* Log grid */
              
              E0 = exp(((double)n/((double)nb))*(log(Emax) - log(Emin)) 
                       + log(Emin));
            }

          /* Adjust for numerical errors */
          
          if (E0 < Emin)
            E0 = Emin;
          
          /* Search index */

          if (E0 == Emin)
            imin = 0;
          else if ((imin = SearchArray(&E[i0], E0, ne)) < 0)
            Die(FUNCTION_NAME, "imin < 0: %E %E %E, %ld", E[i0], E0,
                E[i0 + ne - 1], ne);
          
          /* Calculate energy for upper bin */
          
          if (type == GRID_TYPE_LIN)
            {
              /* Linear grid */
              
              E0 = ((double)(n + 1)/((double)nb))*(Emax - Emin) + Emin;
            }
          else
            {
              /* Log grid */
              
              E0 = exp(((double)(n + 1)/((double)nb))*(log(Emax) - log(Emin)) 
                       + log(Emin));
            }
          
          /* Adjust for numerical errors */
          
          if (E0 > Emax)
            E0 = Emax;
          
          /* Search index */
          
          if (E0 == Emax)
            imax = ne;
          else if ((imax = SearchArray(&E[i0], E0, ne)) < 0)
            Die(FUNCTION_NAME, "imax < 0: %E %E %E, %ld", E[i0], E0,
                E[i0 + ne - 1], ne);
          
          /* Adjust upper bin */
          
          if (imax > ne - 2)
            imax = ne - 1;
          else
            
            imax = imax + 1;

          /* Check indexes */
          
          if (imin >= imax)
            Die(FUNCTION_NAME, "imin = %ld, imax = %ld", imin, imax);
          
          if (!((E0 > -INFTY) && (E0 < INFTY)))
            printf("%ld %ld %ld %ld %ld %ld %E %E\n", type, n, nb, i0, imin, imax, Emin, Emax);

          /* Check energy */
          
          CheckValue(FUNCTION_NAME, "E0", "", E0, E[i0 + imin], E[i0 + imax]);
          
          /* Call recursively */
          
          loc0 = MakeEnergyGrid(imax - imin + 1, i0 + imin, lvl + 1, ptr, E,
                                mode);
          WDB[loc1++] = (double)loc0;              
        }

      /***********************************************************************/
    }

  /* Return pointer to grid */
  
  return erg;        
}

/*****************************************************************************/
