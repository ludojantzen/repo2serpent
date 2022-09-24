/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : surfacevol.c                                   */
/*                                                                           */
/* Created:       2011/07/03 (JLe)                                           */
/* Last modified: 2016/10/25 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates volume inside surface                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SurfaceVol:"

/*****************************************************************************/

double SurfaceVol(long surf)
{
  long ptr, type, params;
  double vol, x, y, z, r, h, r1, r2, a1, a2, l, d;

  /* Check surface pointer */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Pointer to parameter list */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  
  /* Check number of parameters */

  if((type != SURF_INF) && (ptr < 0))
    Die(FUNCTION_NAME, "Surface %s has no parameters",
        GetText(surf + SURFACE_PTR_NAME));

  /* Reset volume */

  vol = -1.0;

  /* Get number of parameters */

  params = (long)RDB[surf + SURFACE_N_PARAMS];

  /***************************************************************************/

  switch(type)
    {
      /***********************************************************************/

    case SURF_CYL:
    case SURF_CYLX:
    case SURF_CYLY:
    case SURF_CYLZ:
      {
        /* Get radius */

        r = RDB[ptr + 2];
        
        /* Get height */

        if (params == 5)
          h = RDB[ptr + 4] - RDB[ptr + 3];
        else
          h = 1.0;

        /* Calculate volume */

        vol = PI*r*r*h;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_PAD:
      {
        /* Get parameters */

        r1 = RDB[ptr + 2];
        r2 = RDB[ptr + 3];
        a1 = RDB[ptr + 4] - (long)(RDB[ptr + 4]/360)*360.0;
        a2 = RDB[ptr + 5] - (long)(RDB[ptr + 5]/360)*360.0;

        /* Calculate volume */

        vol = PI*(r2*r2 - r1*r1)*fabs(a1 - a2)/360.0;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_SPH:
      {
        /* Get radius */

        r = RDB[ptr + 3];

        /* Calculate volume */

        vol = PI*r*r*r*4.0/3.0;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_CUBE:
      {
        /* Get radius */

        r = RDB[ptr + 3];

        /* Calculate volume */

        vol = r*r*r*8.0;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_CUBOID:
      {
        /* Get dimensions */
        
        x = RDB[ptr + 1] - RDB[ptr];
        y = RDB[ptr + 3] - RDB[ptr + 2];
        z = RDB[ptr + 5] - RDB[ptr + 4];

        /* Calculate volume */

        vol = x*y*z;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_RECT:
      {
        /* Get dimensions */
        
        x = RDB[ptr + 1] - RDB[ptr];
        y = RDB[ptr + 3] - RDB[ptr + 2];

        /* Calculate volume */

        vol = x*y;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_PX:
    case SURF_PY:
    case SURF_PZ:
    case SURF_PLANE:
    case SURF_MPLANE:
    case SURF_QUADRATIC:
    case SURF_CKX:
    case SURF_CKY:
    case SURF_CKZ:
    case SURF_X:
    case SURF_Y:
    case SURF_Z:
    case SURF_TORX:
    case SURF_TORY:
    case SURF_TORZ:
      {
        /* Surfaces contain no volume (or not implemented) */

        vol = 0.0;

        /* Break case */

        break;
      }      

      /**********************************************************************/

    case SURF_CONE:
      {
        /* Get parameters */

        r = RDB[ptr + 3];
        h = RDB[ptr + 4];

        /* Calculate volume */

        vol = PI*r*r*h/3.0;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_INF:
      {
        /* Infinite volume */

        vol = INFTY;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_SQC:
      {
        /* Get radius */

        r = RDB[ptr + 2];

        /* Calculate volume */

        vol = r*r*4.0;

        /* Check definition of rounded corners */

        if (params == 4)
          {
            Die(FUNCTION_NAME, "Rounded corners");
          }

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_HEXYC:
    case SURF_HEXXC:
      {
        /* Get radius */

        r = RDB[ptr + 2];

        /* Calculate volume */

        vol = 2.0*SQRT3*r*r;

        /* Check definition of rounded corners */

        if (params == 4)
          Die(FUNCTION_NAME, "Rounded corners");
        
        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_HEXYPRISM:
    case SURF_HEXXPRISM:
      {
        /* Get radius */

        r = RDB[ptr + 2];

        /* Get height */

        if (params == 5)
          h = RDB[ptr + 4] - RDB[ptr + 3];
        else
          h = 1.0;

        /* TBD */

        /* Calculate volume */

        vol = 2.0*SQRT3*r*r*h;

        /* Break case */

        break;
      }
      
      /**********************************************************************/
      
    case SURF_CROSS:
      {
        /* Get parameters */

        l = RDB[ptr + 2];
        d = RDB[ptr + 3];

        /* Check rounded corners */

        if (params == 5)
          Die (FUNCTION_NAME, "Rounded corners");

        /* Calculate volume */

        vol = 4.0*l*l - 4.0*(l - d)*(l - d);

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_SVC:
      {
        /* TBD */

        Die (FUNCTION_NAME, "Cannot calculate volume for surface type %ld",
             type);

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_DODE:
      {
        /* Get radius */

        r = RDB[ptr + 2];

        /* Calculate volume */

        vol = 12.0*(2.0 - SQRT3)*r*r;

        /* Break case */

        break;
      }
        
      /**********************************************************************/
          
    case SURF_OCTA:
      {
        /* Get radius */

        r = RDB[ptr + 2];

        /* Calculate volume */

        vol = 8.0*(SQRT2 - 1.0)*r*r;

        /* Break case */

        break;
      }
      
      /**********************************************************************/
          
    case SURF_ASTRA:
      {
        /* TBD */

        Die (FUNCTION_NAME, "Cannot calculate volume for surface type %ld",
             type);

        /* Break case */

        break;
      }
    
      /***********************************************************************/        

    default:
      {
        Die (FUNCTION_NAME, "Cannot calculate volume for surface type %ld",
             type);

        break;
      }
    }

  /***************************************************************************/
  
  /* Return volume */

  return vol;
}

/*****************************************************************************/
