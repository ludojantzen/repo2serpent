/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculateuresmajorants.c                       */
/*                                                                           */
/* Created:       2011/10/05 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates nuclide-wise majorants in ures region             */
/*                                                                           */
/* Comments: - Used for calculating other majorants                          */
/*           - Nuclide-wise energy grid changed to unionized grid to correct */
/*             a bug in optimization modes 1 and 3 (15.2.2016 / 2.1.25)      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateUresMajorants:"

/*****************************************************************************/

void CalculateUresMajorants(long nuc)
{
  long erg, ptr, pte, rea, i, urs, mp, m, n, nu0, nup, i0, np, ne, e0;
  double *tot, *xs;
  const double *f;

  /* Check nuclide pointer */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check ures flag */

  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED))
    return;
  
  /* Pointer to energy grid */

  ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of grid points */

  ne = (long)RDB[ptr + ENERGY_GRID_NE];

  /* Pointer to data */
  
  e0 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(e0)", DATA_ARRAY, e0);

  /* Allocate memory for temporary arrays (workarray ei toimi tässä */
  /* eikä energiavektorissa koska muistia varataan myöhemmin) */

  tot = (double *)Mem(MEM_ALLOC, ne, sizeof(double));
  xs = (double *)Mem(MEM_ALLOC, ne, sizeof(double));

  /* Get pointer to total xs */

  rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to cross section data array */
          
  ptr = (long)RDB[rea + REACTION_PTR_XS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Pointer to corresponding energy array */

  pte = (long)RDB[rea + REACTION_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);
  
  pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

  /* First point in xs data and number of points */
          
  i0 = (long)RDB[rea + REACTION_XS_I0];
  np = (long)RDB[rea + REACTION_XS_NE];          

  /* Interpolate data into grid */          

  InterpolateData(&RDB[e0], tot, ne, &RDB[pte + i0], &RDB[ptr], np, 0, 
                  NULL, NULL, NO);

  /* Reset indexes */

  nu0 = -1;
  nup = -1;

  /* Loop over partials */

  for (i = 0; i < 3; i++)
    {
      /* Get reaction pointer */

      if (i == 0)
        rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
      else if (i == 1)
        rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
      else if (i == 2)
        rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];

      /* Check pointer */

      if (rea < VALID_PTR)
        continue;

      /* Clear data */

      memset(xs, 0.0, ne*sizeof(double));

      /* Pointer to cross section data array */
          
      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to corresponding energy array */

      pte = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);
      
      pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

      /* First point in xs data and number of points */
          
      i0 = (long)RDB[rea + REACTION_XS_I0];
      np = (long)RDB[rea + REACTION_XS_NE];          
      
      /* Interpolate data into grid */          
      
      InterpolateData(&RDB[e0], xs, ne, &RDB[pte + i0], &RDB[ptr], np, 0,
                      NULL, NULL, NO);
      
#ifdef DEBUG

      /* Check data */

      for (n = 0; n < ne; n++)
        if (xs[n] < 0.0)
          Die(FUNCTION_NAME, "Negative xs (%s mt = %ld)",
              GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT]);

#endif

      /* Get pointer to probability table data */

      urs = (long)RDB[rea + REACTION_PTR_URES];
      CheckPointer(FUNCTION_NAME, "(urs)", DATA_ARRAY, urs);

      /* Get pointer to ures energy grid */
  
      erg = (long)RDB[urs + URES_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg1)", DATA_ARRAY, erg);

      /* Number of points */

      mp = (long)RDB[erg + ENERGY_GRID_NE];

      /* Pointer to data */

      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(erg2)", DATA_ARRAY, erg);

      /* Pointer to maximum factors */
    
      ptr = (long)RDB[urs + URES_PTR_MAXF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Direct pointer (tää on OK, koska pointteria ei käytetä enää */
      /* alempana olevan muistinvarauksen jälkeen) */

      f = &RDB[ptr];

      /* Check for non-positive values */

      for (m = 0; m < mp; m++)
        if (f[m] < 0.0)
          Die(FUNCTION_NAME, "Negative maximum factor (%s, mt = %ld)",
              GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT]);
      
      /* Reset index */

      m = -1;
      nu0 = -1;
      nup = 0;

      /* Loop ures region */

      for (n = 0; n < ne; n++)
        if ((RDB[e0 + n] >= RDB[rea + REACTION_URES_EMIN]) &&
            (RDB[e0 + n] <= RDB[rea + REACTION_URES_EMAX]))
          {
            /* Put first point */

            if (nu0 == -1)
              nu0 = n;
            
            /* Add to counter */

            nup++;

            /* Find correct interval */

            while ((RDB[erg + m + 1] <= RDB[e0 + n]) && (m < mp - 2))
              m++;

#ifdef DEBUG

            ptr = (long)RDB[urs + URES_PTR_EGRID];
            if ((m != GridSearch(ptr, RDB[e0 + n])) &&
                (RDB[e0 + n] < RDB[rea + REACTION_URES_EMAX]))
              Die(FUNCTION_NAME, "Error in ures grid index search %ld %ld",
                  m, GridSearch(ptr, RDB[e0 + n]));

#endif
            /* Adjust or replace value */

            if ((long)RDB[urs + URES_IFF] == 0)
              {
                /* Check non-zero values */
                
                if (f[m] > 0.0)
                  {
                    tot[n] = tot[n] - xs[n] + f[m];
                    xs[n] = f[m];
                  }
              }
            else
              {
                /* Adjust */ 

                tot[n] = tot[n] + (f[m] - 1)*xs[n];
                xs[n] = xs[n]*f[m];
              }

            /* Check non-positive */

            if (tot[n] <= 0.0)
              Die(FUNCTION_NAME, "Non-positive values");
          }
    
      /* Allocate memory for data */

      ptr = ReallocMem(DATA_ARRAY, nup);

      /* Put data */

      memcpy(&WDB[ptr], &xs[nu0], nup*sizeof(double));

      /* Put pointers and values */

      WDB[rea + REACTION_PTR_URES_MAX] = (double)ptr;
      WDB[rea + REACTION_URES_MAX_N0] = (double)nu0;
      WDB[rea + REACTION_URES_MAX_NP] = (double)nup;    
    }

  /* Check index */

  if (nu0 < 0)
    Die(FUNCTION_NAME, "No ures channels");

#ifdef DEBUG

  /* Check one more time */

  for (n = 0; n < nup; n++)
    if (tot[nu0 + n] <= 0.0)
      Die(FUNCTION_NAME, "Non-positive values in total (%s)", 
          GetText(nuc + NUCLIDE_PTR_NAME));

#endif

  /* Allocate memory for data */

  ptr = ReallocMem(DATA_ARRAY, nup);
  
  /* Put data */
  
  memcpy(&WDB[ptr], &tot[nu0], nup*sizeof(double));

  /* Pointer to total */

  rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  
  /* Put pointers and values */
  
  WDB[rea + REACTION_PTR_URES_MAX] = (double)ptr;
  WDB[rea + REACTION_URES_MAX_N0] = (double)nu0;
  WDB[rea + REACTION_URES_MAX_NP] = (double)nup;    

  /* Free allocated memory */

  Mem(MEM_FREE, tot);
  Mem(MEM_FREE, xs);
}

/*****************************************************************************/
