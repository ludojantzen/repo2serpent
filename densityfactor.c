/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : densityfactor.c                                */
/*                                                                           */
/* Created:       2012/02/14 (JLe)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns material density factor at a given location          */
/*                                                                           */
/* Comments: - Multifysiikkarajapinnasta haettava kerroin siirrettiin        */
/*             yhdessä lämpötilan kanssa ifcpoint.c -aliohjelmaan            */
/*             18.3.2013 / 2.1.13.                                           */
/*                                                                           */
/*           - Noi muuttujat pitäis jotenkin pystyä kuitenkin välittämään    */
/*             tonne IFCPoint():lle tässä ja GetTemp():ssä.                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DensityFactor:"

/*****************************************************************************/

double DensityFactor(long mat, double x, double y, double z, double t, long id)
{
  double f, T;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return 1.0;

  /* Reset values */

  f = 1.0;
  T = 0.0;

  /* Get density factor from multi-physics interface */

  IFCPoint(mat, &f, &T, -1.0, id);

  /* Multiply by global factor */

  f = f*RDB[DATA_GLOBAL_DF];

  /* Check value (non-set values are now set to 1.0) */

  if ((f < 0.0) || (f > 1.0))
    {
      /* Check plotter mode */

      if ((long)RDB[DATA_PLOTTER_MODE] == (double)YES)
        return -1.0;
      else
        return 1.0;
    }

  /* Return value */

  return f;
}

/*****************************************************************************/
