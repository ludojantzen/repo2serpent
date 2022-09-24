/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : icmidx.c                                       */
/*                                                                           */
/* Created:       2013/09/15 (JLe)                                           */
/* Last modified: 2013/10/25 (JLe)                                           */
/* Version:       2.1.16                                                     */
/*                                                                           */
/* Description: Calculates surface index for ICM                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ICMIdx:"

/* Surfaces */

#define WEST  0
#define SOUTH 1
#define EAST  2
#define NORTH 3

/*****************************************************************************/

long ICMIdx(long surf, double x, double y, double z, double *u0, double *v0, 
            double *w0, double *u1, double *v1, double *w1, double *u2, 
            double *v2, double *w2)

{
  long ptr, i, j, k, type, N, idx, is;
  double p, f;
  const double *param;
  
  /* Get pointer to parameters */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get parameters */

  param = &RDB[ptr];

  /* Avoid compiler warning */
  
  idx = -1;
  is = -1;

  /* Get type */

  type = (long)RDB[surf + SURFACE_TYPE];

  switch (type)
    {
    case SURF_SQC:
      {
        /* Square prism (toi indeksointi kiertää W-reunan */
        /* N-nurkasta vastapäivään) */

        x = x - param[0];
        y = y - param[1];

        /* Get pitch */

        p = 2.0*param[2] - 5.0*EXTRAP_L;

        /* Get lattice indexes for surface */

        GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_SQC);
 
        /* Get major index */
        
        if ((i == 0) && (j == 0))
          {
            /* Error */

            Die(FUNCTION_NAME, "Error in indexes");
          }
        else if ((i < 0) && (j == 0))
          {
            /* West face */
            
            idx = WEST;
          }
        else if ((j < 0) && (i == 0))
          {
            /* South face */

            idx = SOUTH;
          }
        else if ((i > 0) && (j == 0))
          {
            /* East face */

            idx = EAST;
          }
        else if ((j > 0) && (i == 0))
          {
            /* North face */

            idx = NORTH;
          }
        else if ((i < 0) && (j < 0))
          {
            /* South-West corner */

            if (fabs(x) > fabs(y))
              idx = WEST;
            else
              idx = SOUTH;
          }
        else if ((j < 0) && (i > 0))
          {
            /* South-East corner */

            if (fabs(x) > fabs(y))
              idx = EAST;
            else
              idx = SOUTH;
          }
        else if ((i > 0) && (j > 0))
          {
            /* North-East corner */

            if (fabs(x) > fabs(y))
              idx = EAST;
            else
              idx = NORTH;
          }
        else if ((j > 0) && (i < 0))
          {
            /* North-West corner */

            if (fabs(x) > fabs(y))
              idx = WEST;
            else
              idx = NORTH;
          }
        else
          {
            /* This should not happen anymore */
            
            Die(FUNCTION_NAME, "Not possible");
          }

        /* Modify pitch to get minor indexes right (kerroin > 2 */
        /* numeriikan takia) */

        p = 2.01*param[2];

        /* Avoid compiler warning */

        f = -1.0;

        /* Get minor indexes and surface normals */
        
        if (idx == WEST)
          {
            /* West face */

            f = (0.5*p - y)/p;

            /* Surface vectors */

            *u0 = -1.0;
            *v0 =  0.0;
            *w0 =  0.0;
            
            *u1 =  0.0;
            *v1 = -1.0;
            *w1 =  0.0;

            *u2 =  0.0;
            *v2 =  0.0;
            *w2 =  1.0;
          }
        else if (idx == SOUTH)
          {
            /* South face */
            
            f = (0.5*p + x)/p;

            /* Surface vectors */
            
            *u0 =  0.0;
            *v0 =  1.0;
            *w0 =  0.0;

            *u1 = -1.0;
            *v1 =  0.0;
            *w1 =  0.0;

            *u2 =  0.0;
            *v2 =  0.0;
            *w2 =  1.0;
          }
        else if (idx == EAST)
          {
            /* East face */
            
            f = (0.5*p + y)/p;

            /* Surface vectors */

            *u0 =  1.0;
            *v0 =  0.0;
            *w0 =  0.0;

            *u1 =  0.0;
            *v1 =  1.0;
            *w1 =  0.0;

            *u2 =  0.0;
            *v2 =  0.0;
            *w2 =  1.0;
          }
        else if (idx == NORTH)
          {
            /* North face */
            
            f = (0.5*p - x)/p;

            /* Surface vectors */
            
            *u0 =  0.0;
            *v0 = -1.0;
            *w0 =  0.0;

            *u1 =  1.0;
            *v1 =  0.0;
            *w1 =  0.0;

            *u2 =  0.0;
            *v2 =  0.0;
            *w2 =  1.0;
          }
        else
          {
            /* This should not happen anymore */
            
            Die(FUNCTION_NAME, "Not possible either");
          }

        /* Check factor */

        CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

        /* Get number of sub-segments */

        N = (long)RDB[DATA_ICM_NSUB];

        /* Get pointer to sub-segments */
        
        ptr = (long)RDB[DATA_ICM_PTR_SUB];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get index */
        
        is = SearchArray(&RDB[ptr], f, N + 1);
        CheckValue(FUNCTION_NAME, "is", "", is, 0, N - 1);
        
        /* Check main index one more time */

        CheckValue(FUNCTION_NAME, "idx", "", idx, 0, 3);

        /* Convert */

        idx = idx*N + is;

        /* Break loop */

        break;
      }
    default:
      {
        /* Assume that all other surfaces have a single segment */
        /* (tämä toimii jos joku haluaa laskea yksinkertaisia   */
        /* albedoja) */

        idx = 0;

        /* Break loop */

        break;
      }
    }

  if (idx < 0)
    Die(FUNCTION_NAME, "WTF?");

  /* Return index */

  return idx;
}

/*****************************************************************************/
