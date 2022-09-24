/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gettemp.c                                      */
/*                                                                           */
/* Created:       2012/01/11 (JLe)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns material temperature at given location               */
/*                                                                           */
/* Comments: - Multifysiikkarajapinnasta haettava lämpötila siirrettiin      */
/*             yhdessä tiheyskertoimen kanssa ifcpoint.c -aliohjelmaan       */
/*             18.3.2013 / 2.1.13.                                           */
/*                                                                           */
/*           - Noi muuttujat pitäis jotenkin pystyä kuitenkin välittämään    */
/*             tonne IFCPoint():lle tässä ja GetTemp():ssä.                  */
/*                                                                           */
/*           - Return value is zero if no temperature is given.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetTemp:"

/*****************************************************************************/

double GetTemp(long mat, long id)
{
  double T, f, x, y, z;
  long ptr, uni;

  /***************************************************************************/

  /***** Temperature from interface ******************************************/

  /* Reset values */

  f = 1.0;
  T = -1.0;

  /* Get density factor from multi-physics interface */

  IFCPoint(mat, &f, &T, -1.0, id);

  /* Check value */

  if (T > 0.0)
    return T;

  /***************************************************************************/

  /***** Temperature from material *******************************************/

  /* Get Temperature */

  if ((mat > VALID_PTR) && (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE))
    {
      /* Use minimum temperature */

      T = RDB[mat + MATERIAL_TMS_TMIN];

      /* Check with maximum */

      if (T != RDB[mat + MATERIAL_TMS_TMAX])
        {
          /* Give error, but get coordinates first */

          /* Get collision universe */

          ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
          uni = (long)GetPrivateData(ptr, id);
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Get coordinates */

          ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          x = GetPrivateData(ptr, id);

          ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          y = GetPrivateData(ptr, id);

          ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          z = GetPrivateData(ptr, id);

          Die(FUNCTION_NAME,
              "On-the-fly temperature treatment set for material %s, but "
              "could not find temperature at coordinates (%E %E %E).",
              GetText(mat + MATERIAL_PTR_NAME), x,y,z);
        }
    }
  else
    T = 0.0;

  /* Return value */

  return T;

  /***************************************************************************/
}

/*****************************************************************************/
