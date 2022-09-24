/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : wwimportance.c                                 */
/*                                                                           */
/* Created:       2015/10/08 (JLe)                                           */
/* Last modified: 2020/04/13 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Return the importance for weight windows                     */
/*                                                                           */
/* Comments: - Tässä luupataan kaikkien yli ja palautetaan kumulatiivinen    */
/*             arvo.                                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WWImportance:"

/*****************************************************************************/

double WWImportance(long type, double x0, double y0, double z0, double u,
                    double v, double w, double E, long itype)
{
  long wwd, msh, rmx, erg, ptr, loc0, loc1, ne, n, j, nmax, ng;
  double f, p, c, x, y, z;

  /* Pointer to weight window structure */

  wwd = (long)RDB[DATA_PTR_WWD0];
  CheckPointer(FUNCTION_NAME, "(wwd)", DATA_ARRAY, wwd);

  /* Loop over weight window structures */

  while (wwd > VALID_PTR)
    {
      /* Check particle type */

      if (((long)RDB[wwd + WWD_PARTICLE_TYPE] > 0) &&
          ((long)RDB[wwd + WWD_PARTICLE_TYPE] != type))
        {
          /* Pointer to next */

          wwd = NextItem(wwd);

          /* Cycle loop */

          continue;
        }

      /* Avoid compiler warning */

      f = -1.0;

      /* Check type */

      if ((long)RDB[wwd + WWD_TYPE] == WWD_MESH_TYPE_MCNP)
        {
          /*******************************************************************/

          /***** MCNP Weight window file *************************************/

          /* Get pointer to mesh */

          msh = (long)RDB[wwd + WWD_PTR_MESH];

          /* Check */

          if (msh < VALID_PTR)
            {
              /* Pointer to next */

              wwd = NextItem(wwd);

              /* Cycle loop */

              continue;
            }

          /* Pointer to data */

          if ((ptr = MeshPtr(msh, x0, y0, z0)) > VALID_PTR)
            {
              /* Pointer to structure */

              loc0 = (long)RDB[ptr];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Number of energy groups */

              ne = (long)RDB[wwd + WWD_NE];
              CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

              /* Pointer to array */

              erg = (long)RDB[wwd + WWD_PTR_ERG];
              CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

              /* Get energy group (maksimi ei oo mukana --> ongelmia) */

              if ((n = SearchArray(&RDB[erg], E, ne + 1)) > -1)
                {
                  /* Check value */

                  CheckValue(FUNCTION_NAME, "n", "", n, 0, ne - 1);

                  /* Pointer to importance vector */

                  ptr = (long)RDB[loc0 + WWD_MESH_PTR_IMP];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Get importance */

                  f = RDB[ptr + n];
                }

              /***************************************************************/
            }
        }
      else if (((long)RDB[wwd + WWD_TYPE] == WWD_MESH_TYPE_ITER) ||
               ((long)RDB[wwd + WWD_TYPE] == WWD_MESH_TYPE_SSS))
        {
          /*******************************************************************/

          /***** Serpent type mesh *******************************************/

          /* Ignore collisions */

          if ((itype == WWMESH_COLL) && ((long)RDB[DATA_PLOTTER_MODE] == NO))
            return 0.0;

          /* Get pointer to mesh (not set for first iteration) */

          if ((msh = (long)RDB[wwd + WWD_PTR_MESH]) < VALID_PTR)
            f = 0.0;
          else
            {
              /* Pointer to data */

              if ((ptr = MeshPtr(msh, x0, y0, z0)) < VALID_PTR)
                {
                  /* Pointer to next */

                  wwd = NextItem(wwd);

                  /* Cycle loop */

                  continue;
                }

              /* Pointer to structure */

              loc0 = (long)RDB[ptr];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Pointer to reponse matrix solution */

              rmx = (long)RDB[wwd + WWD_PTR_RMX];
              CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

              /* Get number of energy groups */

              ng = (long)RDB[rmx + RMX_NG];
              CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

              /* Get pointer to few-group structure */

              if ((ptr = (long)RDB[rmx + RMX_PTR_EGRID]) > VALID_PTR)
                {
                  /* Get group index */

                  if ((n = GridSearch(ptr, E)) < 0)
                    Error(wwd,
                          "Particle energy %1.2E beyond energy boundaries", E);
                }
              else
                n = 0;

              /* Check energy group */

              CheckValue(FUNCTION_NAME, "n", "", n, 0, ng - 1);

              /* Check simulation mode */

              if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
                  (1 != 2))
                {
                  /* Pointer to source importances */

                  ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Ignore importances for source cells */

                  if (RDB[ptr + n] != 0.0)
                    return 0.0;
                }

              /* Check mode */

              if ((long)RDB[wwd + WWD_BOUNDS_TYPE] == WW_MESH_BOUNDS_SEP)
                {
                  /************************************************************/

                  /***** Segment-wise importances *****************************/

                  /* Check call type */

                  if (itype == WWMESH_SRC)
                    ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                  else if (itype == WWMESH_COLL)
                    ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
                  else if (itype == WWMESH_BOUND)
                    {
                      /* Move to previous point on track */

                      x = x0 - 2.0*EXTRAP_L*u;
                      y = y0 - 2.0*EXTRAP_L*v;
                      z = z0 - 2.0*EXTRAP_L*w;

                      /* Pointer to data */

                      if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
                        Die(FUNCTION_NAME, "WTF1?");

                      /* Pointer to structure */

                      loc1 = (long)RDB[ptr];
                      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                      /* Compare */

                      if (loc1 == loc0)
                        Die(FUNCTION_NAME, "WTF2?");

                      /* Get maximum number of neighbours */

                      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                      /* Find neighbour cell */

                      ptr = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
                      for (j = 0; j < nmax; j++)
                        {
                          /* Check pointer */

                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          /* Compare */

                          if ((long)RDB[ptr + RMX_CELL_BOUND_PTR_CELL] == loc1)
                            break;

                          /* Next */

                          ptr = NextItem(ptr);
                        }

                      /* Pointer to current importances */

                      if (j < nmax)
                        ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR]
                          + j*ng;
                      else
                        ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];

                      /* Check zero */

                      if (RDB[ptr + n] == 0.0)
                        ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
                    }
                  else
                    Die(FUNCTION_NAME, "Invalid mode %ld", itype);

                  /************************************************************/
                }
              else if ((long)RDB[wwd + WWD_BOUNDS_TYPE] == WW_MESH_BOUNDS_AVG)
                {
                  /************************************************************/

                  /***** Cell-averaged importances ****************************/

                  /* Check call type (collisions ingnored for now) */

                  if (itype == WWMESH_SRC)
                    ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                  else
                    ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];

                  /************************************************************/
                }
              else
                Die(FUNCTION_NAME, "Error in bounds type");

              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

              /* Get importance (n = energy group) */

              f = fabs(RDB[ptr + n]);
            }

          /* Check */

          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Unsupported mesh type");

      /* Get coefficient and exponential */

      if (((c = RDB[wwd + WWD_MULT]) > 0.0) &&
          ((p = RDB[wwd + WWD_POW]) > 0.0))
        f = c*pow(f, p);

      /* Scale by normalization factor */

      f = f*RDB[wwd + WWD_NORM_FACT];
      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

      /* Return */

      if (f > 0.0)
        return f;

      /* Next */

      wwd = NextItem(wwd);
    }

  /* Return zero */

  return 0.0;
}

/*****************************************************************************/
