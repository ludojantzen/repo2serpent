/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movest.c                                       */
/*                                                                           */
/* Created:       2012/10/05 (JLe)                                           */
/* Last modified: 2020/04/09 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Moves particle forward using surface-tracking                */
/*                                                                           */
/* Comments: - WhereAmI() pitää nyt kutsua juuri ennen NearestBoundary():a,  */
/*             koska jälkimmäinen funktio tarvitsee tallennetut suunta-      */
/*             kosinit, jotka on saattaneet muuttua törmäyksessä. Tässä on   */
/*             potentiaalinen ongelmakohta, jos tuota dataa aletaan käyttää  */
/*             myös muuhun tarkoitukseen, ja se unohdetaan päivittää samalla */
/*             tavalla. --> voisi miettiä jonkun paremman ratkaisun. Voi     */
/*             myös olla että noita kutsuja WhereAmI():hin tehdään           */
/*             tarpeettoman monta.                                           */
/*                                                                           */
/*           - Toi uusi moodi ei siirrä pistettä pinnan yli vaan lähelle     */
/*             sitä, minkä jälkeen seuraava piste arvotaan DT:llä. Tuon      */
/*             pitäisi toimia paremmin STL-geometriatyypin kanssa, ja        */
/*             nopeuttaa erityisesti fuusiogeometrioissa missä on isoja      */
/*             vakuumialueita joihin DT jää helposti jumiin.                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveST:"

#define STL_SAFE_DISTANCE 0.5

/*****************************************************************************/

long MoveST(long part, double totxs, double minxs, long *cell, double *xs0,
            double *x, double *y, double *z, double *l0, double u, double v,
            double w, long id)
{
  long ptr, type, cell0, stl;
  double d, xs, l, x0, y0, z0;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "", u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", w, -1.0, 1.0);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Add to ST fraction counter */

  ptr = (long)RDB[RES_ST_TRACK_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, 1.0, ptr, id, 2 - type);

  /***** Sensitivity calculations ********************/
  /* This will force more collisions                 */
  /* 50 % of collisions are rejected in collision.c  */

  /* TODO: Add this to the extra majorants instead (somehow) */

  if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
    totxs *= 2;

  /*****************************/

  /* Get cross section */

  if (totxs < minxs)
    xs = minxs;
  else
    xs = totxs;

  /* Check cross section and sample path length */

  CheckValue(FUNCTION_NAME, "xs1", "", xs, ZERO, INFTY);
  l = -log(RandF(id))/xs;

  /* This call must be made to get the stored direction cosines right */
  /* (may have changed in a collision or when symmetries are invoked) */

  cell0 = WhereAmI(*x, *y, *z, u, v, w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell0);

  /* Get distance to nearest boundary */

  d = NearestBoundary(id);

  x0 = *x;
  y0 = *y;
  z0 = *z;

  /* Check STL mode */

  ptr = (long)RDB[DATA_ST_USE_STL_MODE];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  if ((stl = (long)GetPrivateData(ptr, id)) == YES)
    {
      /* Check domain decomposition */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        Die(FUNCTION_NAME, "DD may not work with STL geometries");

      /* STL mode, reset flag */

      PutPrivateData(ptr, NO, id);

      /* Adjust surface distance */

      if ((d = d - STL_SAFE_DISTANCE) < 0.0)
        d = 0.0;
    }

  /* Compare distances */

  if (l < d)
    {
      /***********************************************************************/

      /***** Move to collision site ******************************************/

      /* Remember initial coordinates */

      x0 = *x;
      y0 = *y;
      z0 = *z;

      /* Move particle to collision site */

      *x = *x + u*l;
      *y = *y + v*l;
      *z = *z + w*l;

      /* Find location */

      *cell = WhereAmI(*x, *y, *z, u, v, w, id);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

      /* Tätä tarvitaan hoitamaan geometriavirheet UMSH-geometrioissa */
      /* (Mikähän tääkin on? 22.7.2015) */

      if ((long)RDB[DATA_PTR_UMSH0] > VALID_PTR)
        if ((*cell < VALID_PTR) ||
            ((long)RDB[*cell + CELL_PTR_MAT] < VALID_PTR))
          {
            /* Set cross section and path length */

            *xs0 = -1.0;
            *l0 = l;

            /* Check distance */

            CheckValue(FUNCTION_NAME, "l0", "1", *l0, 0.0, INFTY);

            /* Return virtual collision */

            return TRACK_END_VIRT;
          }

      /* Check with previous (tää antoi aiheettomia varoituksia universumi- */
      /* symmetrioiden kanssa ennen kuin symmetriarajapinta lisättiin       */
      /* nearestboundary.c:hen). JLe: Error tolerance 9.10.2016 / 2.1.28.   */

      if (cell0 > VALID_PTR)
        if (*cell != cell0)
          {
            /* Enforce delta-tracking for next collision */

            ptr = (long)RDB[DATA_DT_ENFORCE_NEXT_TRACK];
            CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
            PutPrivateData(ptr, YES, id);

            /* Add count */

#ifdef OPEN_MP
#pragma omp atomic
#endif
            WDB[DATA_UNDEF_POS_COUNT] += 1.0;

            /* Check and print warning or error */

            if (((long)RDB[DATA_IGNORE_UNDEFINED_CELLS] == NO) ||
                (((long)RDB[DATA_MAX_UNDEF_POS] > -1) &&
                 ((long)RDB[DATA_UNDEF_POS_COUNT] >
                  (long)RDB[DATA_MAX_UNDEF_POS])))
              Die(FUNCTION_NAME,
                   "Surface crossing missed between [%1.2E, %1.2E, %1.2E] and [%1.2E %1.2E %1.2E] (d = %E, l = %E)",
                   x0, y0, z0, *x, *y, *z, d, l);
            else
              Warn(FUNCTION_NAME,
                   "Surface crossing missed between [%1.2E, %1.2E, %1.2E] and [%1.2E %1.2E %1.2E] (d = %E, l = %E)",
                   x0, y0, z0, *x, *y, *z, d, l);

            /* Reset total xs to enforce virtual collision */

            totxs = 0.0;
          }

      /* Check cross section (could be from previous time step) */

      CheckValue(FUNCTION_NAME, "xs2", "", xs, ZERO, INFTY);

      /* Set cross section and path length */

      *xs0 = xs;
      *l0 = l;

      /* Check distance */

      CheckValue(FUNCTION_NAME, "l0", "2", *l0, 0.0, INFTY);

      /* Sample between real and virtual collision */

      if (RandF(id) < totxs/minxs)
        {
          /* Check material pointer */

	  if ((long)RDB[*cell + CELL_PTR_MAT] < VALID_PTR)
	    {
	      /* NOTE: This has happened in some STL geometries, possibly */
	      /* due to overlapping solids? (JLe / 8.11.2019 / 2.1.32).   */

	      if (stl == NO)
		TrackingError(TRACK_ERR_NO_MATERIAL, -1, -1, -1, id);
	      else
		{
		  Note(0,
			 "Collision in void cell %s at (%1.5E, %1.5E, %1.5E)",
		       GetText(*cell + CELL_PTR_NAME), *x, *y, *z);

		  /* Set cross section and path length */

		  *xs0 = -1.0;
		  *l0 = l;

		  /* Return virtual collision */

		  return TRACK_END_VIRT;
		}
	    }

          /* Return collision */

          return TRACK_END_COLL;
        }
      else
        {
          /* Return virtual */

          return TRACK_END_VIRT;
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Move to surface crossing ****************************************/

      /* Check STL mode */

      if (stl == YES)
        {
          /* Check if moved */

          if (d == 0.0)
            {
              /* Copy cell pointer */

              *cell = cell0;
            }
          else
            {
              /* Move particle near surface */

              *x = *x + u*d;
              *y = *y + v*d;
              *z = *z + w*d;

              /* Find location (tarvitaankohan tätä muuhun kuin */
              /* tohon testiin?) */

              *cell = WhereAmI(*x, *y, *z, u, v, w, id);
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

              /* Check */

              if (cell0 > VALID_PTR)
                if (*cell != cell0)
                  Warn(FUNCTION_NAME,
                       "Possible overlap at (%1.2E, %1.2E, %1.2E)",
                       *x, *y, *z);
            }

          /* Set cross section and path length */

          *xs0 = -1.0;
          *l0 = d;

          /* Check distance */

          CheckValue(FUNCTION_NAME, "l0", "4", *l0, 0.0, INFTY);

          /* Enforce delta-tracking for next collision */

          ptr = (long)RDB[DATA_DT_ENFORCE_NEXT_TRACK];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          PutPrivateData(ptr, YES, id);

          /* Return virtual collision */

          return TRACK_END_VIRT;
        }
      else
        {
          /* Conventional surface-tracking, move particle over surface */

          *x = *x + u*(d + EXTRAP_L);
          *y = *y + v*(d + EXTRAP_L);
          *z = *z + w*(d + EXTRAP_L);

          /* Set cell search list option */
          /*
          ptr = (long)RDB[DATA_CELL_SEARCH_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          PutPrivateData(ptr, (double)CELL_SEARCH_LIST_SURF, id);
          */
          /* Find location */

          *cell = WhereAmI(*x, *y, *z, u, v, w, id);
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

          /* Set cell search list option */
          /*
          ptr = (long)RDB[DATA_CELL_SEARCH_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          PutPrivateData(ptr, (double)CELL_SEARCH_LIST_CELL, id);
          */
          /* Set cross section and path length */

          *xs0 = -1.0;
          *l0 = d + EXTRAP_L;

          /* Check distance */

          CheckValue(FUNCTION_NAME, "l0", "3", *l0, 0.0, INFTY);

          /* Return surface crossing */

          return TRACK_END_SURF;
        }

      /***********************************************************************/
    }

  /* Avoid compiler warning */

  return -1;
}

/*****************************************************************************/
