/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cellvolumes.c                                  */
/*                                                                           */
/* Created:       2012/08/22 (JLe)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Calculates cell volumes                                      */
/*                                                                           */
/* Comments: - Works for a very limited cell types only                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CellVolumes:"

/*****************************************************************************/

void CellVolumes()
{
  long cell, loc0, loc1, surf, ptr, n, type, prnt;
  long tetlist, ntet, tet, i;
  double zmin, zmax, r1, r2, vol;

  /***************************************************************************/

  /***** Main geometry cell types ********************************************/

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Check if volume was already calculated and that cell consists */
      /* of intersections only */

      if ((RDB[cell + CELL_VOLUME] > 0.0) ||
          ((loc0 = (long)RDB[cell + CELL_PTR_SURF_INSC]) < VALID_PTR))
        {
          /* Next cell */

          cell = NextItem(cell);

          /* Cycle loop */

          continue;
        }

      /* Reset min and max dimensions */

      zmin = 0.0;
      zmax = 0.0;
      r1 = 0.0;
      r2 = 0.0;

      /* Reset type */

      type = -1;

      /* Loop over intersection list */

      n = 0;
      while ((loc1 = ListPtr(loc0, n++)) > 0)
        {
          /* Pointer to surface */

          surf = (long)RDB[loc1 + CELL_INSC_PTR_SURF];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Check surface type */

          if ((long)RDB[surf + SURFACE_TYPE] == SURF_PZ)
            {
              /* Z-plane. Get pointer to parameter list. */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Compare to min and max */

              if ((long)RDB[loc1 + CELL_INSC_SIDE] < 0)
                zmax = RDB[ptr];
              else if ((long)RDB[loc1 + CELL_INSC_SIDE] > 0)
                zmin = RDB[ptr];
              else
                Die(FUNCTION_NAME, "Surface error");

              /* Set type */

              if (type == SURF_SPH)
                type = -1;
              else
                type = SURF_CYL;
            }
          else if ((long)RDB[surf + SURFACE_TYPE] == SURF_CYL)
            {
              /* Cylinder. Get pointer to parameter list. */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Compare to min and max */

              if ((long)RDB[loc1 + CELL_INSC_SIDE] < 0)
                r2 = RDB[ptr + 2];
              else if ((long)RDB[loc1 + CELL_INSC_SIDE] > 0)
                r1 = RDB[ptr + 2];
              else
                Die(FUNCTION_NAME, "Surface error");

              /* Check if cylinder is truncated (12.4.2012 / 1.1.18) */

              if ((long)RDB[surf + SURFACE_N_PARAMS] > 3)
                {
                  /* Set type to error */

                  type = -2;

                  /* Break loop */

                  break;
                }

              /* Set type */

              if (type == SURF_SPH)
                type = -1;
              else
                type = SURF_CYL;
            }
          else if ((long)RDB[surf + SURFACE_TYPE] == SURF_SPH)
            {
              /* Sphere. Get pointer to parameter list. */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Compare to min and max */

              if ((long)RDB[loc1 + CELL_INSC_SIDE] < 0)
                r2 = RDB[ptr + 3];
              else if ((long)RDB[loc1 + CELL_INSC_SIDE] > 0)
                r1 = RDB[ptr + 3];
              else
                Die(FUNCTION_NAME, "Surface error");

              /* Set type */

              if (type == SURF_CYL)
                type = -1;
              else
                type = SURF_SPH;
            }
        }

      /* Reset volume */

      vol = 0.0;

      /* Check type */

      if (type == SURF_CYL)
        {
          /* Calculate cross-sectional area */

          if (r2 > 0.0)
            vol = PI*(r2*r2 - r1*r1);

          /* Multiply by height */

          if (zmax > zmin)
            vol = vol*(zmax - zmin);
        }
      else if (type != -2)
        {
          /* Calculate volume */

          if (r2 > 0.0)
            vol = 4.0/3.0*PI*(r2*r2*r2 - r1*r1*r1);
        }

      /* Set value */

      if (vol > 0.0)
        WDB[cell + CELL_VOLUME] = vol;

      /* Next cell */

      cell = NextItem(cell);
    }

  /***************************************************************************/

  /***** Tetrahedral cells in unstructured mesh geometries *******************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_UMSH0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to interface structure */

      loc1 = (long)RDB[loc0 + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      tetlist = (long)RDB[loc1 + IFC_PTR_TET_LIST];
      CheckPointer(FUNCTION_NAME, "(tetlist)", DATA_ARRAY, tetlist);

      ntet = (long)RDB[loc1 + IFC_NC];

      for (i = 0; i < ntet; i++)
        {
          /* Get pointer to tet */

          tet = (long)RDB[tetlist + i];
          CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

          /* Get pointer to parent */

          prnt = (long)RDB[tet + TET_PTR_PARENT];
          CheckPointer(FUNCTION_NAME, "(prnt)", DATA_ARRAY, prnt);

          /* Pointer to geometry cell */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell/tet)", DATA_ARRAY, cell);

          /* Add volume to cell */

          WDB[cell + CELL_VOLUME] += TetraVol(tet);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
