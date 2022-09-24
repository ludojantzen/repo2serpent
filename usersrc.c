/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : usersrc.c                                      */
/*                                                                           */
/* Created:       2012/04/03 (JLe)                                           */
/* Last modified: 2016/07/27 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: User-defined source routine                                  */
/*                                                                           */
/* Comments: - This subroutine is intended to be used for complicated user-  */
/*             defined source routines. The routine sets coordinates,        */
/*             direction, energy ,weight and time for sampled source         */
/*             particles. The routine is called for source definitions using */
/*             the "si" entry:                                               */
/*                                                                           */
/*             src <name> si <np> <v_1> <v_2> ... <v_np>                     */
/*                                                                           */
/*             where <np> is the number of parameters and <v_i> are values   */
/*             passed to the subroutine. The routine can be used in          */
/*             combination with other source parameters.                     */
/*                                                                           */
/*           - Notice that installing new updates will overwrite any user    */
/*             specifide coding, so be sure to make backups before updating. */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UserSrc:"

/*****************************************************************************/

void UserSrc(long src, double *x, double *y, double *z, double *u,  double *v, 
             double *w, double *E, double *wgt, double *t, long id)
{

  /* The entire source code is commented out in order to suppress */
  /* unnecessary compiler warnings about unused variables, etc.   */

#ifdef REMOVE_THIS_IFDEF_STATEMENT_WHEN_USED

  long ptr, np;
  double x0, y0, z0, xmin, xmax, ymin, ymax, zmin, zmax, u0, v0, w0, E0, rnd;
  const double *params;

  /***************************************************************************/

  /***** Common part (don't touch) *******************************************/

  /* Check source pointer */

  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Point source coordinates (from "sp x0 y0 z0"): */

  x0 = RDB[src + SRC_X0];
  y0 = RDB[src + SRC_Y0];
  z0 = RDB[src + SRC_Z0];

  /* Limits (from "sx xmin xmax", "sy ymin ymax", "sz zmin zmax"): */

  xmin = RDB[src + SRC_XMIN];
  xmax = RDB[src + SRC_XMAX];
  ymin = RDB[src + SRC_YMIN];
  ymax = RDB[src + SRC_YMAX];
  zmin = RDB[src + SRC_ZMIN];
  zmax = RDB[src + SRC_ZMAX];

  /* Monodirectional direction cosines (from "sd u0 v0 w0"): */

  u0 = RDB[src + SRC_U0];
  v0 = RDB[src + SRC_V0];
  w0 = RDB[src + SRC_W0];

  /* Energy (from "se E0"): */

  E0 = RDB[src + SRC_E];

  /* Get pointer to user-data */
  
  ptr = (long)RDB[src + SRC_PTR_USR];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of input values */

  np = (long)RDB[ptr + SRC_USR_NP];

  /* Get pointer to values */

  ptr = (long)RDB[ptr + SRC_USR_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  params = &RDB[ptr];

  /***************************************************************************/

  /***** User-defined part ***************************************************/
  
  /* Particle coordinates, direction, energy, weight and time are set  */
  /* in the following for monodirectional, monoenergetic point source. */
  /* The number of parameters entered after "si" is stored in variable */
  /* np, and the values in array params. Notice that in C the indexing */
  /* is zero-based, so the first given parameter is params[0], and the */
  /* last parameter is params[np - 1]. Other variables defined in the  */
  /* source card are:                                                  */
  /*                                                                   */
  /* x0, y0, z0  -- Point source coordinates (values following "sp")   */
  /* xmin, xmax  -- Boundaries on the x-axis (values following "sx")   */
  /* ymin, ymax  -- Boundaries on the y-axis (values following "sy")   */
  /* zmin, zmax  -- Boundaries on the z-axis (values following "sz")   */
  /* u0, v0, w0  -- Direction cosines (values following "sd")          */
  /* E0          -- Source energy (value following "se")               */
  /*                                                                   */
  /* The values of x, y, z, u, v, w, E wgt and t, passed as function   */
  /* arguments are initially set based on the other entries in the     */
  /* source card. This subroutine is called after all other routines.  */
  /*                                                                   */
  /* Notice that altough these values can be used in the subroutine,   */
  /* they are not automatically assigned to the corresponding          */
  /* variables of the source particle.                                 */

  /* Override coordinates (point-source at origin) */

  *x = 0.0;
  *y = 0.0;
  *z = 0.0;

  /* Sample uniformly distributed random number on the unit interval */
  /* (just to give an example, the value is not used for anything).  */

  rnd = RandF(id);

  /* Override direction cosines. The vector length must be normalized   */
  /* to 1.0. Subroutine IsotropicDirection() samples a random direction */
  /* that is automatically normalized. */

  IsotropicDirection(u, v, w, id);

  /* Override energy (in MeV) */

  if (1 == 2)
    {
      /* Monoenergetic source */
      
      *E = 1.0;
    }
  else
    {
      /* Maxwellian spectrum */

      if (np > 0)
        *E = MaxwellEnergy(params[0], id);
      else
        Die(FUNCTION_NAME, "kT not provided");
    }

  /* Override weight to unity */

  *wgt = 1.0;

  /* Override time */

  *t = 0.0;     

#endif

  /***************************************************************************/
}

/*****************************************************************************/
