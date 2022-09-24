/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processlattices.c                              */
/*                                                                           */
/* Created:       2010/10/10 (JLe)                                           */
/* Last modified: 2015/09/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Creates lattice surfaces, etc.                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessLattices:"

/*****************************************************************************/

void ProcessLattices()
{
  long lat, ptr, reg1, reg2, n;
  double r1, r2;

  /* Loop over lattices */

  lat = (long)RDB[DATA_PTR_L0];
  while (lat > VALID_PTR)
    {
      /* Check type */
      
      switch ((long)RDB[lat + LAT_TYPE])
        {
        case LAT_TYPE_CLU:
          {
            /* Get pointer to rings */
              
            ptr = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                   
            /* Loop over rings */ 

            for (n = 0; n < (long)RDB[lat + LAT_N_RINGS]; n++)
              {
                /* First ring */

                reg1 = ListPtr(ptr, n);
                r1 = RDB[reg1 + RING_RAD];

                /* Second ring */

                if (n < (long)RDB[lat + LAT_N_RINGS] - 1)
                  {
                    reg2 = ListPtr(ptr, n + 1);
                    r2 = RDB[reg2 + RING_RAD];
                  }
                else
                  r2 = INFTY;

                /* Check values */

                if (r1 < 0.0)
                  Error(lat, "Invalid ring radius %E\n", r1);

                if (r2 < 0.0)
                  Error(lat, "Invalid ring radius %E\n", r2);

                if (r1 > r2)
                  Error(lat, "Ring radii must be in ascending order");

                /* Put limit */
                
                WDB[reg1 + RING_RLIM] = (r1 + r2)/2.0;
              }

            /* Break case */

            break;
          }
        case LAT_TYPE_XPRISM:
        case LAT_TYPE_YPRISM:
          {
            /* Check that x- and y-pitch are the same */

            if (RDB[lat + LAT_PITCHX] != RDB[lat + LAT_PITCHY])
              Error(lat, "x- and y-pitches must be the same");

            /* Break case */

            break;
          }
        }

      /* Allocate memory for collision counter */

      AllocValuePair(lat + LAT_COL_CELL_IDX);

      /* Next lattice */
      
      lat = NextItem(lat);
    }
}

/*****************************************************************************/
