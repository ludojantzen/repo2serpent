/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : usersurf.c                                     */
/*                                                                           */
/* Created:       2012/03/30 (JLe)                                           */
/* Last modified: 2012/12/24 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: User-defined surface routines                                */
/*                                                                           */
/* Comments: - This subroutine is intended to simplify the development of    */
/*             new surface types. The routine consists of three parts,       */
/*             described below:                                              */
/*                                                                           */
/*             1) Calculating the maximum dimensions of the surface          */
/*             2) Testing if point (x,y,z) is inside or outside the surface  */
/*             3) Calculating the shortest distance to the surface from      */
/*                point (x,y,z) in direction (u,v,w)                         */
/*                                                                           */
/*             Surface parameters are passed into vector "params" of size    */
/*             np from input card:                                           */
/*                                                                           */
/*             surf <name> usr <p_1> <p_2> ... <p_np>                        */
/*                                                                           */
/*             i.e. surface type "usr", with an arbitrary number of          */
/*             parameters.                                                   */
/*                                                                           */
/*             The example coding handles a cuboid defined by planes at      */
/*             x0, x1, y0, y1, z0, z1.                                       */
/*                                                                           */
/*           - This is the only routine that needs to be touched when adding */
/*             a new surface type.                                           */
/*                                                                           */
/*           - Notice that installing new updates will overwrite any user    */
/*             specifide coding, so be sure to make backups before updating. */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UserSurf:"

/*****************************************************************************/

void UserSurf(long part, long np, const double *params, double *xmin, 
              double *xmax, double *ymin, double *ymax, double *zmin, 
              double *zmax, long *in, double *dmin, double x, double y, 
              double z, double u, double v, double w)
{
  double d;

  /***************************************************************************/

  /***** Common part (do not touch) ******************************************/

  /* Check parameters pointer */

  if ((np > 0) && (params == NULL))
    Die(FUNCTION_NAME, "Pointer error");

  /***************************************************************************/

  /***** User-defined part ***************************************************/

  /* Variable np gives the number of surface parameters read from the   */
  /* input. Array params gives the parameters. Notice that in C the     */
  /* indexing is zero-based, so the first given parameter is params[0], */
  /* and the last parameter is params[np - 1]. */

  if (part == 1)
    {
      /***********************************************************************/

      /***** Calculate maximum surface dimensions ****************************/

      /* This part determines the maximum dimensions of the surface. The   */
      /* procedure is not mandatory, but if the dimensions are not set,    */
      /* the surface cannot be used as the outer boundary of the geometry. */
      /* If it is, the sampling of initial source  will most likely fail   */
      /* in criticality source calculation, or other problems may occur in */
      /* source and plotter routines. */

      /* For the example cuboidal surface the maximum dimensions are */
      /* simply: */

      *xmin = params[0];
      *xmax = params[1];
      *ymin = params[2];
      *ymax = params[3];
      *zmin = params[4];
      *zmax = params[5];

      /***********************************************************************/
    }
  else if (part == 2)
    {
      /***********************************************************************/

      /***** Test if point is inside surface *********************************/

      /* This part tests if a point (x,y,z) is inside the surface or not */
      /* and sets variable *in = YES or *in = NO accordingly. This part  */
      /* must always be defined. */

      /* Tests for the example cuboidal surface */

      if (x < params[0])
        *in = NO;
      else if (x > params[1])
        *in = NO;
      else if (y < params[2])
        *in = NO;
      else if (y > params[3])
        *in = NO;
      else if (z < params[4])
        *in = NO;
      else if (z > params[5])
        *in = NO;
      else
        *in = YES;

      /***********************************************************************/
    }
  else if (part == 3)
    {
      /***********************************************************************/

      /***** Calculate minimum distance to surface ***************************/

      /* This part calculates the shortest distance to the surface from   */
      /* point (x,y,z) in direction (u,v,w) and stores the value in       */
      /* variable *dmin. The plotter routine may produce some results     */
      /* even if this part is not defined, but the tracking routine will  */
      /* fail. Notice that the distance doesn't have to be the absolute   */
      /* minimum, it just must never be longer than the absolute minimum. */
      /* In the cuboidal surface example case the routine may return the  */
      /* distance to planes that extend beyond the actual boundaries when */
      /* the particle is outside, but this is not a problem because the   */
      /* actual boundary is never missed. Basically the subroutine        */
      /* sees the cuboid as:                                              */
      /*                                                                  */
      /*   |      |                                                       */
      /* --+------+--            +------+                                 */
      /*   |      |    and not:  |      |                                 */
      /* --+------+--            +------+                                 */
      /*   |      |                                                       */
      /*                                                                  */
      /* The routine should set variable *dmin to INFTY if the surface is */
      /* not within line-of-sight. */

      /* Reset the minimum distance */

      *dmin = INFTY;

      /* Planes on x-axis */

      if (u != 0.0)
        {
          if (((d = -(x - params[0])/u) > 0.0) && (d < *dmin))
            *dmin = d;
          if (((d = -(x - params[1])/u) > 0.0) && (d < *dmin))
            *dmin = d;
        }

      /* Planes on y-axis */
        
      if (v != 0.0)
        {
          if (((d = -(y - params[2])/v) > 0.0) && (d < *dmin))
            *dmin = d;
          if (((d = -(y - params[3])/v) > 0.0) && (d < *dmin))
            *dmin = d;          
        }

      /* Planes on z-axis */

      if (w != 0.0)
        {
          if (((d = -(z - params[4])/w) > 0.0) && (d < *dmin))
            *dmin = d;
          if (((d = -(z - params[5])/w) > 0.0) && (d < *dmin))
            *dmin = d;
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid part");

  /***************************************************************************/
}

/*****************************************************************************/
