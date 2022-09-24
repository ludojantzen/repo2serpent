/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nearestboundary.c                              */
/*                                                                           */
/* Created:       2010/10/10 (JLe)                                           */
/* Last modified: 2020/01/15 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Finds distance to the nearest boundary surface               */
/*                                                                           */
/* Comments: - 2014/03/13 added time and angular support for fuelperf. ifc   */
/*                                                                           */
/*           - Täällä on aika pitkiä pätkiä interface- ja umsh-juttua,       */
/*             jotka voisi ehkä joskus tulevaisuudessa laittaa omiksi        */
/*             aliohjelmikseen samaan tyyliin kuin NearestSTLSurf().         */
/*           - Tarviikohan tolle polttoainerajapinnalle laskea noita etäi-   */
/*             syyksiä lämpötilajakauman eri aksiaalikerroksiin ja radiaali  */
/*             noodeihin, lähinnä varmaan jos joku TMS-majorantti lasketaan  */
/*             erikseen eri alueilla   (?)                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestBoundary:"

/*****************************************************************************/

double NearestBoundary(long id)
{
  long lvl0, lvl, reg, cell, pbd, pbl, surf, type, n, np, ptr, loc0, ltype;
  long nbhr, uni, ncol, ang, i, j, loc1, tet, pt0, pt1, pt2, uni0, surf0;
  double min, d, x, y, z, u, v, w, y2, z2, params[MAX_SURFACE_PARAMS];
  double t, phi, phi2, min0;
  const long *face[4];
  static const long face1[3] = {0,2,1};
  static const long face2[3] = {1,2,3};
  static const long face3[3] = {0,3,2};
  static const long face4[3] = {0,1,3};

  /* Reset variables */

  min = INFTY;
  min0 = -1.0;
  uni0 = -1;
  surf0 = -1;

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Disabled for now */

      break;

      /* Reset surface pointer */

      ptr = (long)RDB[uni + UNIVERSE_PTR_NEAREST_SURF];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      PutPrivateData(ptr, NULLPTR, id);

      uni = NextItem(uni);
    }


  /* Pointer to first level */

  lvl0 = (long)RDB[DATA_PTR_LVL0];
  CheckPointer(FUNCTION_NAME, "(lvl0)", DATA_ARRAY, lvl0);

  /* Loop over levels */

  while (lvl0 > VALID_PTR)
    {
      /* Pointer to private data */

      lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
      CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

      /* Get coordinates */

      x = GetPrivateData(lvl + LVL_PRIV_X, id);
      y = GetPrivateData(lvl + LVL_PRIV_Y, id);
      z = GetPrivateData(lvl + LVL_PRIV_Z, id);

      /* Get direction cosines */

      u = GetPrivateData(lvl + LVL_PRIV_U, id);
      v = GetPrivateData(lvl + LVL_PRIV_V, id);
      w = GetPrivateData(lvl + LVL_PRIV_W, id);

      /* Check coordinates and direction cosines */

      CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "u", "", u, -1.0, 1.0);
      CheckValue(FUNCTION_NAME, "v", "", v, -1.0, 1.0);
      CheckValue(FUNCTION_NAME, "w", "", w, -1.0, 1.0);

      /* Pointer to universe */

      uni = (long)GetPrivateData(lvl + LVL_PRIV_PTR_UNIV, id);
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Get time */

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_T];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      t = GetPrivateData(ptr, id);

      /* Check for symmetry */

      if ((ptr = (long)RDB[uni + UNIVERSE_PTR_SYM]) > VALID_PTR)
        if ((d = SymmetryBoundary(ptr, x, y, z, u, v, w)) < min)
          min = d;

      /* Check for divisor (NOTE: tänne ei pitäisi mennä transport-vaiheessa */
      /* ollenkaan, sillä dt pakotetaan päälle asettamalla MATERIAL_DT_MODE  */
      /* DT_MAT_FORCE:ksi dividezone.c:ssä. VolumesMC() rutiinissa tänne     */
      /* mennään, sillä kutsu tulee suoraan. (15.1.2020 / 2.1.32 / JLE)      */

      if ((ptr = (long)GetPrivateData(lvl + LVL_PRIV_PTR_DIV, id)) > VALID_PTR)
        if ((d = DivSubDistance(ptr, uni, x, y, z, u, v, w, t)) < min)
          min = d;

      /* Get level type */

      ltype = (long)GetPrivateData(lvl + LVL_PRIV_TYPE, id);

      /* Check type */

      switch (ltype)
        {
        case UNIVERSE_TYPE_NEST:
          {
            /*****************************************************************/

            /***** Nest ******************************************************/

            /* Check for fuel performance interface */

            if ((loc0 = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
              {
                /* First check nest boundaries */
                /* Check nest region surfaces */

                /* Pointer to nest region */

                reg = (long)GetPrivateData(lvl + LVL_PRIV_PTR_NEST_REG, id);
                CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

                /* First surface */

                if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_IN]) > VALID_PTR)
                  {
                    /* Get type */

                    type = (long)RDB[surf + SURFACE_TYPE];

                    /* Die on non-cylindrical */

                    if(type != SURF_CYL)
                      Die(FUNCTION_NAME,"Non-cylindrical surface");

                    /* Get number of parameters */

                    np = (long)RDB[surf + SURFACE_N_PARAMS];

                    /* Pointer to parameter list */

                    ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    /* Copy surface parameters */

                    memcpy(&params,&RDB[ptr],np*sizeof(double));

                    y2 = 0.0;
                    z2 = z;

                    /* Cold to hot expansion of cylinder radius */
                    /* (params[2]) */

                    CoordExpans(loc0,&params[2],&y2,&z2,t,2);

                    /* Get distance */

                    d = SurfaceDistance(-1, params, type, np, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "1", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      min = d;
                  }

                /* Second surface */

                if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT]) >
                    VALID_PTR)
                  {
                    /* Get type */

                    type = (long)RDB[surf + SURFACE_TYPE];

                    /* Die on non-cylindrical */

                    if(type != SURF_CYL)
                      Die(FUNCTION_NAME,"Non-cylindrical surface");

                    /* Get number of parameters */

                    np = (long)RDB[surf + SURFACE_N_PARAMS];

                    /* Pointer to parameter list */

                    ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    /* Copy surface parameters */

                    memcpy(&params,&RDB[ptr],np*sizeof(double));

                    y2 = 0.0;
                    z2 = z;

                    /* Cold to hot expansion of cylinder radius */
                    /* (params[2]) */

                    CoordExpans(loc0,&params[2],&y2,&z2,t,2);

                    /* Get distance */

                    d = SurfaceDistance(-1, params, type, np, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "2", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      min = d;
                  }

                /* Find axial zone */

                loc1 = (long)RDB[loc0 + IFC_FUEP_PTR_AX];

                while (loc1 > VALID_PTR)
                  {
                    /* Compare coordinates */

                    if ((z >= RDB[loc1 + IFC_FUEP_AX_ZMIN]) &&
                        (z < RDB[loc1 + IFC_FUEP_AX_ZMAX]))
                      {

                        /* Break loop */

                        break;
                      }

                    /* Next */

                    loc1 = NextItem(loc1);
                  }

                /* Break case if not found*/

                if (loc1 < VALID_PTR)
                  break;

                /* Only check axial boundaries in 3D calculation */

                if ((long)RDB[DATA_GEOM_DIM] == 3)
                  {
                    /* Distance to lower boundary */

                    params[0] = RDB[loc1 + IFC_FUEP_AX_ZMIN];
                    d = SurfaceDistance(-1, params, SURF_PZ, 1, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "3", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      min = d;

                    /* Distance to upper boundary */

                    params[0] = RDB[loc1 + IFC_FUEP_AX_ZMAX];
                    d = SurfaceDistance(-1, params, SURF_PZ, 1, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "4", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      min = d;
                  }

                /* Get the correct angular segment */

                ang = (long)RDB[loc1 + IFC_FUEP_AX_PTR_ANG];

                /* Get polar angle */

                phi = PolarAngle(x,y);

                while (ang > VALID_PTR)
                  {
                    /* Rotate if needed */

                    if(phi > 2.0*PI+RDB[ang + IFC_FUEP_ANG_AMIN])
                      phi2 = phi - 2.0*PI;
                    else
                      phi2 = phi;

                    /* Compare coordinates */

                    if ((phi2 >= RDB[ang + IFC_FUEP_ANG_AMIN]) &&
                        (phi2 < RDB[ang + IFC_FUEP_ANG_AMAX]))
                      break;

                    /* Next */

                    ang = NextItem(ang);
                  }

                /* Break case if not found */

                if (ang < VALID_PTR)
                  break;

                /* Calculate distance to angular zone boundaries            */
                /* Tää ei huomioi että rajoittava taso on vain puoliääretön */
                /* (Tulee ylimääräisiä pysäytyksiä, jos kulmasegmenttejä    */
                /*  on pariton määrä)                                       */

                params[0] = RDB[ang + IFC_FUEP_ANG_CMIN];
                params[1] = 1.0;
                params[2] = 0.0;
                params[3] = 0.0;

                d = SurfaceDistance(-1, params, SURF_PLANE, 4, x, y, z,
                                    u, v, w, id);
                CheckValue(FUNCTION_NAME, "d", "5", d, 0.0, INFTY);

                /* Compare to minimum */

                if (d < min)
                  min = d;

                params[0] = RDB[ang + IFC_FUEP_ANG_CMAX];

                d = SurfaceDistance(-1, params, SURF_PLANE, 4, x, y, z,
                                    u, v, w, id);
                CheckValue(FUNCTION_NAME, "d", "6", d, 0.0, INFTY);

                /* Compare to minimum */

                if (d < min)
                  min = d;

              }
            else
              {
                /* Pointer to nest region */

                reg = (long)GetPrivateData(lvl + LVL_PRIV_PTR_NEST_REG, id);
                CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

                /* First surface */

                if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_IN]) > VALID_PTR)
                  {
                    /* Get type */

                    type = (long)RDB[surf + SURFACE_TYPE];

                    /* Get number of parameters */

                    np = (long)RDB[surf + SURFACE_N_PARAMS];

                    /* Pointer to parameter list */

                    ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    /* Get distance */

                    d = SurfaceDistance(-1, &RDB[ptr], type, np, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "10", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      min = d;
                  }

                /* Second surface */

                if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT]) > VALID_PTR)
                  {
                    /* Get type */

                    type = (long)RDB[surf + SURFACE_TYPE];

                    /* Get number of parameters */

                    np = (long)RDB[surf + SURFACE_N_PARAMS];

                    /* Pointer to parameter list */

                    ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    /* Get distance */

                    d = SurfaceDistance(-1, &RDB[ptr], type, np, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "11", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      min = d;
                  }
              }

            /* Break case */

            break;

            /*****************************************************************/
          }

        case UNIVERSE_TYPE_LATTICE:
          {
            /*****************************************************************/

            /***** Lattice ***************************************************/

            /* Check pointer to lattice */

            if ((long)GetPrivateData(lvl + LVL_PRIV_PTR_LAT, id) < VALID_PTR)
              Die(FUNCTION_NAME, "No lattice pointer");

            /* Get surface type */

            type = (long)GetPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, id);

            /* Get number of surface parameters */

            np = (long)GetPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, id);

            /* Get parameters (noiden kerrointen pitää olla peräkkäin) */

            for (n = 0; n < np; n++)
              params[n] = GetPrivateData(lvl + LVL_PRIV_LAT_SURF_C0 + n, id);

            /* Get distance */

            d = SurfaceDistance(-1, params, type, np, x, y, z, u, v, w, id);
            CheckValue(FUNCTION_NAME, "d", "12", d, 0.0, INFTY);

            /* Compare to minimum */

            if (d < min)
              min = d;

            /* Check if type is vertical stack */

            if (type == SURF_PZ)
              {
                /* Put parameter for second surface */

                params[0] = GetPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, id);

                /* Get distance (ei tarkisteta, reunimmaiset pinnat voi */
                /* olla +/- INFTY) */

                d = SurfaceDistance(-1, params, type, np, x, y, z, u, v, w, id);

                /* Compare to minimum */

                if (d < min)
                  min = d;
              }

            /* Break case */

            break;

            /*****************************************************************/
          }

        case UNIVERSE_TYPE_CELL:
          {
            /*****************************************************************/

            /***** Cell ******************************************************/

            /* Pointer to cell */

            cell = (long)GetPrivateData(lvl + LVL_PRIV_PTR_CELL, id);
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Pointer to surface list */

            loc0 = (long)RDB[cell + CELL_PTR_SURF_LIST];
            CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

            /* Loop over list */

            while ((surf = (long)RDB[loc0++]) > VALID_PTR)
              {
                /* Get type */

                type = (long)RDB[surf + SURFACE_TYPE];

                /* Check infinite */

                if (type != SURF_INF)
                  {
                    /* Get number of parameters */

                    np = (long)RDB[surf + SURFACE_N_PARAMS];

                    /* Pointer to parameter list */

                    ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    /* Get distance */

                    d = SurfaceDistance(surf, &RDB[ptr], type, np, x, y, z,
                                        u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "14", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      {
                        /* Put minimum distance */

                        min = d;

                        /* Put universe values */

                        min0 = min;
                        uni0 = uni;
                        surf0 = surf;
                      }
                  }
              }

            /* Break case */

            break;

            /*****************************************************************/
          }

        case UNIVERSE_TYPE_UMSH:
          {
            /*****************************************************************/

            /***** Unstructured mesh based geometry **************************/

            /* Pointer to cell */

            if ((cell = (long)GetPrivateData(lvl + LVL_PRIV_PTR_CELL, id))
                > VALID_PTR)
              {
                /* Get current tetra-cell (set in whereami.c) */

                ptr = (long)RDB[cell + CELL_PTR_PREV_TET];

                if ((tet = (long)GetPrivateData(ptr, id)) < VALID_PTR)
                  Die(FUNCTION_NAME, "UMSH tet was not stored");

                /* Reset pointer to neighbour */

                nbhr = -1;

                /* Store point list for easy access (faceX are constant lists) */

                face[0] = face1;
                face[1] = face2;
                face[2] = face3;
                face[3] = face4;

                /* Loop over tet faces */

                for (i = 0; i < 4; i++)
                  {
                    /* Get pointers to points */

                    pt0 = (long)RDB[tet + TET_POINTS + face[i][0]];
                    pt1 = (long)RDB[tet + TET_POINTS + face[i][1]];
                    pt2 = (long)RDB[tet + TET_POINTS + face[i][2]];

                    /* Copy points to params */

                    for (j = 0; j < 3; j++)
                      {
                        /* Store coordinates to params */

                        params[0 + j] = RDB[pt0 + j];
                        params[3 + j] = RDB[pt1 + j];
                        params[6 + j] = RDB[pt2 + j];
                      }

                    /* Get distance */

                    d = SurfaceDistance(-1, params, SURF_PLANE, 9,
                                        x, y, z, u, v, w, id);
                    CheckValue(FUNCTION_NAME, "d", "15", d, 0.0, INFTY);

                    /* Compare to minimum */

                    if (d < min)
                      {
                        min = d;

                        /* Get neighbour */

                        nbhr = (long)RDB[tet + TET_NEIGHBOURS + i];
                      }
                  }

                /* Get collision count */

                ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
                ncol = (long)GetPrivateData(ptr, id);

                /* Store value */

                StoreValuePair(uni + UNIVERSE_PTR_NEXT_CELL, (double)ncol,
                               (double)nbhr, id);
              }
            else
              {
                /* Point is in background universe, get distance to search */
                /* mesh boundaries. */

                /* Pointer to geometry */

                loc0 = (long)GetPrivateData(lvl + LVL_PRIV_PTR_UMSH, id);
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Pointer to interaface structure */

                loc0 = (long)RDB[loc0 + UMSH_PTR_IFC];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Get minimum distance */

                d = NearestUMSHSurf(loc0, x, y, z, u, v, w, id);

                /* Compare to minimum */

                if (d < min)
                  min = d;
              }

            /* Break case */

            break;

            /*****************************************************************/
          }

        case UNIVERSE_TYPE_STL:
          {
            /*****************************************************************/

            /***** STL based geometry ****************************************/

            /* Pointer to geometry */

            loc0 = (long)GetPrivateData(lvl + LVL_PRIV_PTR_STL, id);
            CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

            /* Get minimum distance */

            d = NearestSTLSurf(loc0, x, y, z, u, v, w, id);

            /* Compare to minimum */

            if (d < min)
              min = d;

            /* Break case */

            break;

            /*****************************************************************/
          }

        case UNIVERSE_TYPE_PBED:
          {
            /*****************************************************************/

            /***** Explicit stochastic geometry ******************************/

            /* Check direct pointer to pebble */

            if ((pbl = (long)GetPrivateData(lvl + LVL_PRIV_PTR_PEBBLE, id)) >
                VALID_PTR)
              {
                /* Put surface parameters */

                params[0] = RDB[pbl + PEBBLE_X0];
                params[1] = RDB[pbl + PEBBLE_Y0];
                params[2] = RDB[pbl + PEBBLE_Z0];
                params[3] = RDB[pbl + PEBBLE_RAD];

                if ((x - params[0])*(x - params[0]) +
                    (y - params[1])*(y - params[1]) +
                    (z - params[2])*(z - params[2]) > params[3]*params[3])
                  Die(FUNCTION_NAME, "not inside");

                /* Get distance */

                d = SurfaceDistance(-1, params, SURF_SPH, 4, x, y, z,
                                    u, v, w, id);
                CheckValue(FUNCTION_NAME, "d", "16", d, 0.0, 2.0*params[3]);
              }
            else
              {
                /* Pointer to PB geometry */

                pbd = (long)GetPrivateData(lvl + LVL_PRIV_PTR_PBED, id);
                CheckPointer(FUNCTION_NAME, "(pbd)", DATA_ARRAY, pbd);

                /* Get minimum distance */

                d = NearestPBSurf(pbd, x, y, z, u, v, w, id);
              }

            /* Compare to minimum */

            if (d < min)
              min = d;

            /* Break case */

            break;

            /*****************************************************************/
          }

        default:
          {
            /* Invalid type */

            Die(FUNCTION_NAME, "Invalid universe type");
          }
        }

      /* Break loop if level is last */

      if ((long)GetPrivateData(lvl + LVL_PRIV_LAST, id) == YES)
              break;

      /* Next level */

      lvl0 = NextItem(lvl0);
    }

  /* Check if shortest distance was in a cell universe (disabled) */

  if (1 == 2)
  if (min0 == min)
    {
      /* Check universe and surface pointer */

      CheckPointer(FUNCTION_NAME, "(uni0)", DATA_ARRAY, uni0);
      CheckPointer(FUNCTION_NAME, "(surf0)", DATA_ARRAY, surf0);

      /* Put surface pointer */

      ptr = (long)RDB[uni0 + UNIVERSE_PTR_NEAREST_SURF];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      PutPrivateData(ptr, surf0, id);
    }

  /* Return shortest distance */

  return min;
}

/*****************************************************************************/
