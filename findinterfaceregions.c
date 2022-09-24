/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findinterfaceregions.c                         */
/*                                                                           */
/* Created:       2012/02/15 (JLe)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Recursive algorithm that determines output regions for       */
/*              multi-physics interface                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindInterfaceRegions:"

/*****************************************************************************/

void FindInterfaceRegions(long loc0, long  uni, long lvl, long recu, double x0,
                          double y0, double z0, long idx0)
{
  long loc1, loc2, nst, ptr, reg, n, cell, lat, lst, mat0, mat, nx, ny;
  long nr, i0, j0, i, j, type, surf, idx, idx1, ifcmat;
  double x, y, z, pitch, wdth, phi, rad;

  /* Check flag */

  if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == NO)
    return;

  /* Get pointer to output material */

  ifcmat = (long)RDB[loc0 + IFC_PTR_OUTPUT_MAT];
  CheckPointer(FUNCTION_NAME, "(ifcmat)", DATA_ARRAY, ifcmat);

  /* Update level pointer or get pointer to first level and universe */

  if (lvl < VALID_PTR)
    {
      lvl = (long)RDB[DATA_PTR_LVL0];
      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
    }
  else
    lvl = NextItem(lvl);

  /* Check level and universe pointers */

  CheckPointer(FUNCTION_NAME, "(lvl)", DATA_ARRAY, lvl);
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check infinite loop */

  if (recu++ > 1000)
    Die(FUNCTION_NAME, "Infinite geometry loop involving universe %s",
        GetText(uni + UNIVERSE_PTR_NAME));

  /* Coordinate transformation to local origin */

  if ((ptr = (long)RDB[uni + UNIVERSE_PTR_TRANS]) > VALID_PTR)
    {
      x0 = x0 + RDB[ptr + TRANS_X0];
      y0 = y0 + RDB[ptr + TRANS_Y0];
      z0 = z0 + RDB[ptr + TRANS_Z0];
    }

  /* Check symmetries */

  if ((long)RDB[uni + UNIVERSE_PTR_SYM] > VALID_PTR)
    Error(0, "Interface output doesn't work with symmetries (to be fixed)");

  /* Check universe type */

  switch((long)RDB[uni + UNIVERSE_TYPE])
    {
    case UNIVERSE_TYPE_NEST:
      {
        /***** Nest universe *************************************************/

        /* Pointer to nest */

        nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
        CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

        /* Get pointer to regions */

        reg = (long)RDB[nst + NEST_PTR_REGIONS];
        CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

        /* Loop over regions */

        while (reg > VALID_PTR)
          {
            /* Get region index and update global index */

            idx = (long)RDB[reg + NEST_REG_IDX];
            idx1 = idx0 + ((long)RDB[lvl + LVL_ZONE_IDX_MULT])*idx;

            /* Check fill pointer */

            if ((uni = (long)RDB[reg + NEST_REG_PTR_FILL]) > VALID_PTR)
              {
                /* Filled region, call recursively */

                FindInterfaceRegions(loc0, uni, lvl, recu, x0, y0, z0, idx1);
              }

            /* Next region */

            reg = NextItem(reg);
          }

        /* Get pointer to outermost region */

        reg = (long)RDB[nst + NEST_PTR_REGIONS];
        reg = LastItem(reg);

        /* Get pointer to cell */

        cell = (long)RDB[reg + NEST_REG_PTR_CELL];
        CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

        /* Check that the region is not the only one and get pointer */
        /* to material */

        if ((PrevItem(reg) > VALID_PTR) &&
            ((mat0 = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR))
          {
            /* Reset pointer */

            mat = -1;

            /* Compare name */

            if (mat0 == ifcmat)
              mat = mat0;

            /* Check if material was divided for burnup calculation */

            if ((ptr = (long)RDB[mat0 + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
              if (ptr == ifcmat)
                mat = mat0;

            /* Check match */

            if (mat > VALID_PTR)
              {
                /* Get pointer to surface */

                surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT];
                CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

                /* Check type */

                if ((long)RDB[surf + SURFACE_TYPE] != SURF_CYL)
                  Error(loc0, "Interface allowed only with cylindrical nests");

                /* Pointer to parameter list */

                ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                /* Get radius */

                rad = RDB[ptr + 2];

                /* Create new structure */

                loc1 = NewItem(loc0 + IFC_PTR_OUT, IFC_OUT_LIST_BLOCK_SIZE);

                /* Put data */

                WDB[loc1 + IFC_OUT_X0] = x0;
                WDB[loc1 + IFC_OUT_Y0] = y0;
                WDB[loc1 + IFC_OUT_R] = rad;
                WDB[loc1 + IFC_OUT_PTR_IFC] = (double)loc0;

                /* Add scoring regions */

                reg = (long)RDB[nst + NEST_PTR_REGIONS];
                while (reg > VALID_PTR)
                  {
                    /* Exclude last */

                    if (NextItem(reg) < VALID_PTR)
                      break;

                    /* Get pointer to cell */

                    cell = (long)RDB[reg + NEST_REG_PTR_CELL];
                    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

                    /* Pointer to material */

                    mat = (long)RDB[cell + CELL_PTR_MAT];

                    /* Check material pointer and fissile flag */

                    if (mat > VALID_PTR)
                      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT)
                        {
                          /* Get region index global index */

                          idx = (long)RDB[reg + NEST_REG_IDX];
                          idx1 = idx0
                            + ((long)RDB[lvl + LVL_ZONE_IDX_MULT])*idx;

                          /* Create new structure */

                          loc2 = NewItem(loc0 + IFC_PTR_SCORE,
                                         IFC_SCORE_LIST_BLOCK_SIZE);

                          /* Put region index */

                          WDB[loc2 + IFC_SCORE_REG_IDX] = (double)idx1;

                          /* Reset stat index */

                          WDB[loc2 + IFC_SCORE_STAT_IDX] = -1.0;

                          /* Put pointers */

                          WDB[loc2 + IFC_SCORE_PTR_OUT] = (double)loc1;
                          WDB[loc1 + IFC_OUT_PTR_SCORE] = (double)loc2;
                        }

                    /* Next region */

                    reg = NextItem(reg);
                  }
              }
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_CELL:
      {
        /***** Cell universe *************************************************/

        /* Pointer to cell list */

        lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
        CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

        /* Loop over cell list */

        while (lst > VALID_PTR)
          {
            /* Get region index and update global index */

            idx = (long)RDB[lst + CELL_LIST_REG_IDX];
            idx1 = idx0 + ((long)RDB[lvl + LVL_ZONE_IDX_MULT])*idx;

            /* Pointer to cell */

            cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Check fill pointer */

            if ((uni = (long)RDB[cell + CELL_PTR_FILL]) > 0)
              {
                /* Filled region, call recursively */

                FindInterfaceRegions(loc0, uni, lvl, recu, x0, y0, z0, idx1);
              }

            /* Next */

            lst = NextItem(lst);
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_LATTICE:
      {
        /***** Lattice universe **********************************************/

        /* Pointer to lattice */

        lat = (long)RDB[uni + UNIVERSE_PTR_LAT];
        CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

        /* Check type */

        if ((long)RDB[lat + LAT_TYPE] == LAT_TYPE_CLU)
          {
            /***** Circular array ********************************************/

            /* Get pointer to rings */

            reg = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

            /* Reset index */

            idx = 0;

            /* Loop over rings */

            while (reg > VALID_PTR)
              {
                /* Pointer to items */

                ptr = (long)RDB[reg + RING_PTR_FILL];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                /* Number of sectors */

                nr = (long)RDB[reg + RING_N_SEC];

                /* Sector width */

                wdth = 2.0*PI/((double)nr);

                /* Loop over sectors */

                for (n = 0; n < (long)RDB[reg + RING_N_SEC]; n++)
                  {
                    /* Get global index and update local */

                    idx1 = idx0 + ((long)RDB[lvl + LVL_ZONE_IDX_MULT])*idx;
                    idx++;

                    /* Sector center angle */

                    phi = (double)n*wdth + RDB[reg + RING_TILT];

                    /* Adjust */

                    while (phi < 0.0)
                      phi = phi + 2.0*PI;

                    while (phi >= 2.0*PI)
                      phi = phi - 2.0*PI;

                    /* Transfer co-ordinates */

                    x = x0 + RDB[reg + RING_RAD]*cos(phi);
                    y = y0 + RDB[reg + RING_RAD]*sin(phi);
                    z = z0;

                    /* Pointer to universe */

                    uni = (long)RDB[ptr + n];
                    CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

                    /* Call recursively */

                    FindInterfaceRegions(loc0, uni, lvl, recu, x, y, z, idx1);
                  }

                /* Next ring */

                reg = NextItem(reg);
              }

            /*****************************************************************/
          }
        else
          {
            /***** Simple types **********************************************/

            /* Get parameters */

            nx = (long)RDB[lat + LAT_NX];
            ny = (long)RDB[lat + LAT_NY];
            pitch = RDB[lat + LAT_PITCH];
            type = (long)RDB[lat + LAT_TYPE];

            /* Transfer coordinates */

            x0 = x0 + RDB[lat + LAT_ORIG_X0];
            y0 = y0 + RDB[lat + LAT_ORIG_Y0];

            /* If even number of cells, shift origin by pitch/2 */

            x0 = x0 + (1 - (double)(nx % 2))*0.5*pitch;
            y0 = y0 + (1 - (double)(ny % 2))*0.5*pitch;

            /* Get pointer to lattice cells */

            ptr = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

            /* Reset index */

            idx = 0;

            /* Loop over items */

            for (j0 = 0; j0 < ny; j0++)
              for (i0 = 0; i0 < nx; i0++)
                {
                  /* Get global index and update local */

                  idx1 = idx0 + ((long)RDB[lvl + LVL_ZONE_IDX_MULT])*idx;
                  idx++;

                  /* Transfer to centered indexing system */

                  i = i0 - (long)((double)nx/2.0);
                  j = j0 - (long)((double)ny/2.0);

                  /* Avoid compiler warning */

                  x = 0.0;
                  y = 0.0;
                  z = z0;

                  /* Calculate local coordinates */

                  if (type == LAT_TYPE_S)
                    {
                      x = x0 + (double)i*pitch;
                      y = y0 + (double)j*pitch;
                    }
                  else if (type == LAT_TYPE_HX)
                    {
                      x = x0 + ((double)i + COS60*(double)j)*pitch;
                      y = y0 + (double)j*SIN60*pitch;
                    }
                  else if (type == LAT_TYPE_HY)
                    {
                      x = x0 + (double)j*SIN60*pitch;
                      y = y0 + ((double)i + COS60*(double)j)*pitch;
                    }
                  else
                    Die (FUNCTION_NAME, "Unsupported lattice type %ld", type);

                  /* Index to lattice element */

                  n = i0 + nx*j0;

                  /* Call recursively */

                  if ((uni = (long)RDB[ptr + n]) > VALID_PTR)
                    FindInterfaceRegions(loc0, uni, lvl, recu, x, y, z, idx1);
                }

            /*****************************************************************/
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    }
}

/*****************************************************************************/
