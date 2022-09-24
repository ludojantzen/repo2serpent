/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : userdf.c                                       */
/*                                                                           */
/* Created:       2012/05/09 (JLe)                                           */
/* Last modified: 2018/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: User-defined density factor for multi-physics interface      */
/*                                                                           */
/* Comments: - This subroutine is intended to be used for user-defined       */
/*             density and temperature distributions with the multi-physics  */
/*             interface (type 3).                                           */
/*                                                                           */
/*             The function arguments are:                                   */
/*                                                                           */
/*             mat     : Pointer to the material where the collision occurs  */
/*             f       : The density factor, i.e. the ratio of local density */
/*                       to the value given in the main input file (variable */
/*                       adens). The factor must be between 0.0 and 1.0.     */
/*             T       : Local temperature. The value must be between Tmin   */
/*                       and Tmax, or zero, if the value is not used         */
/*             x, y, z : Coordinates where the collision occurs relative to  */
/*                       the root universe (usually universe "0")            */
/*             np      : Number of parameters in the interface file (type 3) */
/*             params  : Vector of parameters                                */
/*                                                                           */
/*             Additional variables:                                         */
/*                                                                           */
/*             adens   : Atomic density of the interface material, as given  */
/*                       in the main input file                              */
/*             Tmin    : Minimum allowed temperature                         */
/*             Tmin    : Maximum allowed temperature                         */
/*                                                                           */
/*           - Notice that indexing in C is zero-based, so the first and     */
/*             last input parameters are params[0] and params[np - 1],       */
/*             respectively.                                                 */
/*                                                                           */
/*           - Notice that installing new updates will overwrite any user    */
/*             specifide coding, so be sure to make backups before updating. */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UserIFC:"

/*****************************************************************************/

void UserIFC(long mat, double *f, double *T, double x, double y, double z, 
             double t, long np, const double *params)
{
  long n;
  double adens, Tmin, Tmax;

  /***************************************************************************/

  /***** Get additional variables (don't touch) ******************************/

  /* Atomic density given in material card */

  adens = RDB[mat + MATERIAL_ADENS];

  /* Minimum and maximum temperature */

  Tmin = RDB[mat + MATERIAL_TMS_TMIN];
  Tmax = RDB[mat + MATERIAL_TMS_TMAX];

  /* Reset state variables */

  *f = 1.0;
  *T = 0.0;

  /***************************************************************************/
  
  /***** Transient behaviour for M&C 2013 paper ******************************/

#ifdef mmmmmmmmmmm
  double t0, t1, f0, f1;

  /* Get parameters */

  t0 = params[0];
  t1 = params[1];
  f0 = params[2];
  f1 = params[3];

  /* Compare to limits */

  if (t < t0)
    return f0;
  else if (t > t1)
    return f1;
  else 
    f = (t - t0)/(t1 - t0)*(f1 - f0) + f0;
  
  /* Exit */

  return f;

#endif
  
  /***************************************************************************/

  /* Polynomial dependence axial density distribution as an example */

  *f = params[0];
  for (n = 1; n < np; n++)
    *f = *f + params[n]*pow(z, (double)n);

  /* Distribution for the example case in NSE paper: */

  if (z < params[0])
    *f = 1.0;
  else
    *f = params[2] - (1.0 - params[3])*
      ((z - params[0])/(params[1] - params[0]))*
      ((z - params[0])/(params[1] - params[0]));
}

/*****************************************************************************/
