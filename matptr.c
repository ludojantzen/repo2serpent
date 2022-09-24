/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : matptr.c                                       */
/*                                                                           */
/* Created:       2012/05/10 (JLe)                                           */
/* Last modified: 2019/04/04 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns material pointer for depletion zones                 */
/*                                                                           */
/* Comments: - Reunat lasketaan alueeseen mukaan koska muuten tää saattaa    */
/*             feilata pintalähteen kanssa.                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MatPtr:"

/*****************************************************************************/

long MatPtr(long mat, long id)
{
  long ptr, lvl0, lvl, div, lst, lst0, idx, nx, ny, nz, nrad, nseg, reg, uni;
  long i, j, k, fpe, n, sz;
  double x, y, z, t, r, phi;

  /* Check null pointer */

  if (mat < VALID_PTR)
    return mat;

  /* Check divisor */

  if ((div = (long)RDB[mat + MATERIAL_PTR_DIV]) < VALID_PTR)
    return mat;

  /* Check type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
    Die(FUNCTION_NAME, "Error");

  /***************************************************************************/

  /***** Division into zones *************************************************/

  /* Pointer to material list */

  lst0 = (long)RDB[div + DIV_PTR_MAT_LIST];
  CheckPointer(FUNCTION_NAME, "(lst0)", DATA_ARRAY, lst0);

  /* Check division */

  if ((long)RDB[div + DIV_SEP] == YES)
    {
      /* Check if divided at top level */
      /*
      if ((long)RDB[div + DIV_SEP_LVL] == 0)
      */

      /* Tää ei nyt jostain syystä pelaa */

      if (1 == 2)
        {
          /* Get top level zone index */

          ptr = (long)RDB[DATA_PTR_ZONE_IDX];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          idx = (long)GetPrivateData(ptr, id);

          /* Find region */

          if ((lst = SeekList(lst0, DIV_MAT_LIST_ZONE_IDX, (double)idx,
                              SORT_MODE_ASCEND)) < VALID_PTR)
            Die(FUNCTION_NAME, "Unable to find divide zone");
        }
      else
        {
          /* Index to highest level where material is located */

          i = (long)RDB[DATA_GEOM_LEVELS] - (long)RDB[div + DIV_LVL_MAX] - 1;

          /* Get level pointer */

          lvl0 = (long)RDB[DATA_PTR_LVL0];
          CheckPointer(FUNCTION_NAME, "(lvl0)", DATA_ARRAY, lvl0);
          lvl0 = LastItem(lvl0);

          /* Loop to correct level */

          for (j = 0; j < (long)RDB[div + DIV_SEP_LVL] + i; j++)
            if ((lvl0 = PrevItem(lvl0)) < VALID_PTR)
              Die(FUNCTION_NAME, "Error in level");

          /* Reset pointer and index */

          lst = -1;
          idx = -1;

          /* Loop over levels from here to down */

          while (lvl0 > VALID_PTR)
            {
              /* Pointer to private data */

              lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
              CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

              /* Get index */

              if ((idx = (long)GetPrivateData(lvl + LVL_PRIV_ZONE_IDX, id))
                  > -1)
                {
                  /* New (potentially) faster search option added */
                  /* 5.4.2019 (2.1.31 / JLe). */

                  if (1 == 2)
                    {
                      /* Find region */

                      lst = SeekList(lst0, DIV_MAT_LIST_ZONE_IDX, (double)idx,
                                     SORT_MODE_ASCEND);
                    }
                  else
                    {
                      /* Pointer to search list */

                      ptr = (long)RDB[div + DIV_PTR_SEARCH_LIST];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      /* Get list size */

                      sz = (long)RDB[div + DIV_SEARCH_LIST_SZ];
                      CheckValue(FUNCTION_NAME, "sz", "", sz, 1, INFTY);

                      /* Search array and get pointer */

                      if ((n = SeekArray(&RDB[ptr], idx, sz)) > -1)
                        lst = ListPtr(lst0, n);
                      else
                        lst = -1;
                    }

                  /* Check */

                  if (lst > VALID_PTR)
                    break;
                }

              /* Pointer to previous level */

              lvl0 = PrevItem(lvl0);
            }

          /* Check pointer */

          if (lst < VALID_PTR)
           {
             /* Check top level. NOTE: ei ymmärrä ihan tarkkaan mistä */
             /* nää errorit tulee, mutta tällä voi joskus korjaantua. */
             /* (JLe 22.12.2016 / 2.1.26) */

             if (((long)RDB[div + DIV_SEP_LVL] + 1) !=
                 ((long)RDB[DATA_GEOM_LEVELS] -
                  (long)RDB[div + DIV_LVL_MAX] - 1))
               Error(div,
                     "Error in level parameter, try replacing \"sep %ld\" with \"sep %ld\"",
                     (long)RDB[div + DIV_SEP_LVL] + 1,
                     (long)RDB[div + DIV_SEP_LVL] + 2);
             else
               Die(FUNCTION_NAME, "Unable to find divide zone");
           }
        }

      /* Check index */

      if (idx != (long)RDB[lst + DIV_MAT_LIST_ZONE_IDX])
        Die(FUNCTION_NAME, "WTF?");
    }
  else
    lst = lst0;

  /***************************************************************************/

  /***** Division into sub-regions *******************************************/

  /* Get universe pointer */

  uni = (long)RDB[lst + DIV_MAT_LIST_PTR_UNIV];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get pointer to sub-regions */

  reg = (long)RDB[lst + DIV_MAT_LIST_PTR_REG];
  CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

  /* Get bin sizes */

  nx = (long)RDB[div + DIV_NX];
  ny = (long)RDB[div + DIV_NY];
  nz = (long)RDB[div + DIV_NZ];
  nrad = (long)RDB[div + DIV_NRAD];
  nseg = (long)RDB[div + DIV_NSEG];

  /* Check for single-valued */

  if (nx*ny*nz*nrad*nseg == 1)
    idx = 0;
  else
    {
      /* Get local coordinates */

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
      x = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
      y = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
      z = GetPrivateData(ptr, id);

      /* Get time */

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_T];
      t = GetPrivateData(ptr, id);

      /* Coordinate change to cold conditions */

      if ((fpe = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
        CoordExpans(fpe, &x, &y, &z, t, 1);

      /* Reset indexes */

      i = 0;
      j = 0;
      k = 0;

      /* Check type */

      if ((nx > 1) || (ny > 1))
        {
          /*******************************************************************/

          /***** Cartesian division ******************************************/

          /* x-coordinate */

          if (nx > 1)
            {
              /* Calculate normalized coordinate */

              x = (x - RDB[div + DIV_XMIN])/
                (RDB[div + DIV_XMAX] - RDB[div + DIV_XMIN]);

              /* Check mode */

              if (((long)RDB[div + DIV_LIMS_CHECK] == YES) &&
                  ((x < 0.0) || (x > 1.0)))
                Error(div, "X-dimension doesn't cover region");
              else if (x == 1.0)
                i = nx - 1;
              else if (x < 0.0)
                i = -1;
              else
                i = (long)(x*nx);
            }

          /* y-coordinate */

          if (ny > 1)
            {
              /* Calculate normalized coordinate */

              y = (y - RDB[div + DIV_YMIN])/
                (RDB[div + DIV_YMAX] - RDB[div + DIV_YMIN]);

              /* Check mode */

              if (((long)RDB[div + DIV_LIMS_CHECK] == YES) &&
                  ((y < 0.0) || (y > 1.0)))
                Error(div, "Y-dimension doesn't cover region");
              else if (y == 1.0)
                j = ny - 1;
              else if (y < 0.0)
                j = -1;
              else
                j = (long)(y*ny);
            }

          /* z-coordinate */

          if (nz > 1)
            {
              /* Calculate normalized coordinate */

              z = (z - RDB[div + DIV_ZMIN])/
                (RDB[div + DIV_ZMAX] - RDB[div + DIV_ZMIN]);

              /* Check mode */

              if (((long)RDB[div + DIV_LIMS_CHECK] == YES) &&
                  ((z < 0.0) || (z > 1.0)))
                Error(div, "Z-dimension doesn't cover region");
              else if (z == 1.0)
                k = nz - 1;
              else if (z < 0.0)
                k = -1;
              else
                k = (long)(z*nz);
            }

          /* Calculate index */

          if ((i < 0) || (i > nx - 1) || (j < 0) || (j > ny - 1) ||
              (k < 0) || (k > nz - 1))
            idx = nx*ny*nz;
          else
            idx = i + j*nx + k*nx*ny;

          CheckValue(FUNCTION_NAME, "idx", "", idx, 0, nx*ny*nz);
        }
      else
        {
          /*******************************************************************/

          /***** Axial coordinates *******************************************/

          /* z-coordinate */

          if (nz > 1)
            {
              /* Calculate normalized coordinate */

              z = (z - RDB[div + DIV_ZMIN])/
                (RDB[div + DIV_ZMAX] - RDB[div + DIV_ZMIN]);

              /* Check mode */

              if (((long)RDB[div + DIV_LIMS_CHECK] == YES) &&
                  ((z < 0.0) || (z > 1.0)))
                Error(div, "Axial dimension doesn't cover region");
              else if (z == 1.0)
                i = nz - 1;
              else if (z < 0.0)
                i = -1;
              else
                i = (long)(z*nz);
            }

          /* Radial coordinate */

          if (nrad > 1)
            {
              /* Calculate square radius */

              r = x*x + y*y;

              /* Calculate portion of volume inside this radius */

              r = (r - RDB[div + DIV_RMIN]*RDB[div + DIV_RMIN])/
                (RDB[div + DIV_RMAX]*RDB[div + DIV_RMAX] -
                 RDB[div + DIV_RMIN]*RDB[div + DIV_RMIN]);

              /* Check mode */

              if (((long)RDB[div + DIV_LIMS_CHECK] == YES) &&
                  ((r < 0.0) || (r > 1.0)))
                Error(div, "Radial dimension doesn't cover region");
              else if (r == 1.0)
                j = nrad - 1;
              else if (r < 0.0)
                j = -1;
              else
                j = (long)(r*nrad);
            }

          /* Angular segment */

          if (nseg > 1)
            {
              /* Calculate angle */

              phi = PolarAngle(x, y);

              /* Add tilt */

              phi = phi + RDB[div + DIV_SEG0];

              /* Normalize */

              phi = phi/(2.0*PI);
              phi = phi - (double)((long)phi);

              /* Check phi */

              CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 1.0);

              /* Calculate index */

              k = (long)(phi*nseg);
            }

          /* Calculate index */

          if ((i < 0) || (i > nz - 1) || (j < 0) || (j > nrad - 1))
            idx = nz*nrad*nseg;
          else
            idx = i + j*nz + k*nz*nrad;

          CheckValue(FUNCTION_NAME, "idx", "", idx, 0, nz*nrad*nseg - 1);

          /*******************************************************************/
        }
    }

  /* Get pointer */

  mat = (long)RDB[reg + idx];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /***************************************************************************/

  /* Return pointer */

  return mat;
}

/*****************************************************************************/
