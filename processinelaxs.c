/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processinelaxs.c                               */
/*                                                                           */
/* Created:       2016/08/09 (JLe)                                           */
/* Last modified: 2016/08/09 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Process total inelastic cross sections (mainly for burnup    */
/*              calculation and production of isomeric states)               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessInelaXS:"

/*****************************************************************************/

void ProcessInelaXS(long nuc)
{
  long rea, rea0, loc0, ptr, mt, ne0, ne1, i0, i1, n;
  double Emin;

  /***************************************************************************/

  /***** Calculate sum if not pre-calculated *********************************/

  /* Find inelastic reaction */

  rea0 = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea0 > VALID_PTR)
    {
      /* Check mt */

      if (((long)RDB[rea0 + REACTION_MT] == 4) && 
          ((long)RDB[rea0 + REACTION_TYPE] == REACTION_TYPE_SPECIAL))
        break;

      /* Next */

      rea0 = NextItem(rea0);
    }

  /* Check (ei voi käyttää VALID_PTR) */

  if (rea0 > VALID_PTR)
    if (RDB[rea0 + REACTION_PTR_EGRID] < 0.0)
      {
        /* Reset pointer and minimum energy */

        rea = -1;
        Emin = INFTY;

        /* Find first inelastic reaction */

        ptr = (long)RDB[nuc + NUCLIDE_PTR_REA];
        while (ptr > VALID_PTR)
          {
            /* Get mt */

            mt = (long)RDB[ptr + REACTION_MT];
            
            /* Check mt */
            
            if ((mt > 50) && (mt < 92))
              {
                /* Compare energies */

                if (RDB[ptr + REACTION_EMIN] < Emin)
                  {
                    /* Update pointer and energy */

                    rea = ptr;
                    Emin = RDB[ptr + REACTION_EMIN];
                  }
              }
            
            /* Next */

            ptr = NextItem(ptr);
          }

        /* Check reaction */

        if (rea < VALID_PTR)
          Die(FUNCTION_NAME, "No inelastic levels found");

        /* Copy pointers and data */
        
        WDB[rea0 + REACTION_PTR_EGRID] = RDB[rea + REACTION_PTR_EGRID];
        WDB[rea0 + REACTION_XS_I0] = RDB[rea + REACTION_XS_I0];
        WDB[rea0 + REACTION_XS_NE] = RDB[rea + REACTION_XS_NE];

        /* Allocate memory for previous value */

        AllocValuePair(rea0 + REACTION_PTR_PREV_XS);

        /* Get index to first */
        
        i0 = (long)RDB[rea0 + REACTION_XS_I0];

        /* Get number of points */

        ptr = (long)RDB[rea + REACTION_PTR_EGRID];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
        ne0 = (long)RDB[ptr + ENERGY_GRID_NE] - i0;

         /* Allocate memory for cross sections */

        loc0 = ReallocMem(DATA_ARRAY, ne0);
        WDB[rea0 + REACTION_PTR_XS] = (double)loc0;

        /* Sum over reactions */

        rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
        while (rea > VALID_PTR)
          {
            /* Get mt */

            mt = (long)RDB[rea + REACTION_MT];
            
            /* Check mt */
            
            if ((mt > 50) && (mt < 92))
              {
                /* Index to first, number of points and pointer to data */

                i1 = (long)RDB[rea + REACTION_XS_I0];
                CheckValue(FUNCTION_NAME, "i1", "", i1, i0, i0 + ne0);

                ne1 = (long)RDB[rea + REACTION_XS_NE];
                CheckValue(FUNCTION_NAME, "ne1", "", ne1, 2, ne0);

                ptr = (long)RDB[rea + REACTION_PTR_XS];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                for (n = 0; n < ne1; n++)
                  WDB[loc0 + i1 - i0 + n] = RDB[loc0 + i1 - i0 + n] 
                    + RDB[ptr + n];
              }
            
            /* Next */

            rea = NextItem(rea);
          }                
      }

  /***************************************************************************/
  
  /***** Handle branch to isomeric state *************************************/

  /* Find inelastic reaction */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check mt */

      if (((long)RDB[rea + REACTION_MT] == 4) && 
          ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_TRA_BRANCH))
        break;

      /* Next */

      rea = NextItem(rea);
    }

  /* Check pointer */

  if (rea > VALID_PTR)
    {
      /* Pointer to parent */
      
      rea0 = (long)RDB[rea + REACTION_PTR_BRANCH_PARENT];
      CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

      /* Copy pointers and other data */

      WDB[rea + REACTION_PTR_EGRID] = RDB[rea0 + REACTION_PTR_EGRID];
      WDB[rea + REACTION_PTR_XS] = RDB[rea0 + REACTION_PTR_XS];
      WDB[rea + REACTION_XS_I0] = RDB[rea0 + REACTION_XS_I0];
      WDB[rea + REACTION_XS_NE] = RDB[rea0 + REACTION_XS_NE];
      WDB[rea + REACTION_PTR_PREV_XS] = RDB[rea0 + REACTION_PTR_PREV_XS];
    }

  /***************************************************************************/
}

/*****************************************************************************/
