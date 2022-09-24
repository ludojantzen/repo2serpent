/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testunisym.c                                   */
/*                                                                           */
/* Created:       2019/10/08 (JLe)                                           */
/* Last modified: 2019/10/08 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Tests for universe symmetry based on real reflection or      */
/*              translation.                                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestUniSym:"

/*****************************************************************************/

long TestUniSym(double x, double y, double z)
{
  long uni, sym;
  double x0, y0, z0, u, v, w;


  /* Get pointer to root universe */

  uni =  (long)RDB[DATA_PTR_ROOT_UNIVERSE];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get pointer to symmetry */

  if ((sym = (long)RDB[uni + UNIVERSE_PTR_SYM]) < VALID_PTR)
    return NO;

  /* Check if true symmetry */

  if ((long)RDB[sym + SYMMETRY_COORD_TRANS] == NO)
    return NO;

  /* Remember original coordinates */

  x0 = x;
  y0 = y;
  z0 = z;

  /* Put direction cosines */

  u = 0;
  v = 0;
  w = 0;

  /* Call transformation */

  UniSym(sym, &x, &y, &z, &u, &v, &w);

  /* Check */

  if ((x == x0) && (y == y0) && (z == z0))
    return NO;
  else
    return YES;
}

/*****************************************************************************/
