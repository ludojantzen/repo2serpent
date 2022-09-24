/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoregc.c                                      */
/*                                                                           */
/* Created:       2011/05/04 (JLe)                                           */
/* Last modified: 2019/04/24 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores reaction rates needed for group constant generation   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreGC:"

/*****************************************************************************/

void ScoreGC(double flx, double tot, double capt, double fiss, double fissE,
             double nsf, long mat, long part, double E, double spd, double wgt,
             double u, double v, double w, long id)
{
  long gcu, ptr, loc0, iso, ng, ntot, ncol, uni, dng;
  double lambda, f, xs;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Score micro-depletion data */

  ScoreMicroDep(flx, mat, E, wgt, id);

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Check for multiple levels */

  if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
    {
      /* Single level, get pointer */

      if ((gcu = (long)TestValuePair(DATA_GCU_PTR_UNI, (double)ncol, id))
          < VALID_PTR)
        return;
    }
  else
    {
      /* Multiple levels, get pointer to list */

      gcu = (long)RDB[DATA_PTR_GCU0];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
    }

  /****************************************************************************/

  /***** Score data **********************************************************/

  /* Loop over universes */

  while (gcu > VALID_PTR)
    {
      /* Check multi-level mode */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == YES)
        {
          /* Pointer to universe */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Check collision */

          if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id) < 0.0)
            {
              /* Next universe */

              gcu = NextItem(gcu);

              /* Cycle loop */

              continue;
            }
        }

      /***********************************************************************/

      /***** MORA data *******************************************************/

      if ((loc0 = (long)RDB[gcu + GCU_PTR_MORA]) > VALID_PTR)
        {
          /* Get pointer to energy grid */

          ptr = (long)RDB[loc0 + MORA_PTR_EG];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Total number of groups */

          ntot = (long)RDB[loc0 + MORA_N_EG];

          /* Find group */

          if ((ng = GridSearch(ptr, E)) > -1)
            {
              /* Flux */

              ptr = (long)RDB[loc0 + MORA_PTR_FLX];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(flx, wgt, ptr, id, ng);
              AddBuf1D(flx, wgt, ptr, id, ntot);

              /* Score total reaction rate */

              ptr = (long)RDB[loc0 + MORA_PTR_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(tot, wgt, ptr, id, ng);

              /* Score fission rate */

              ptr = (long)RDB[loc0 + MORA_PTR_FISS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(fiss, wgt, ptr, id, ng);

              /* Score capture rate */

              ptr = (long)RDB[loc0 + MORA_PTR_CAPT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(capt, wgt, ptr, id, ng);

              /* Score fission energy production rate */

              ptr = (long)RDB[loc0 + MORA_PTR_KAPPA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(fissE, wgt, ptr, id, ng);
            }
        }

      /************************************************************************/

      /***** Micro-group data ************************************************/

      /* Get pointer to microgroup energy grid */

      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Number of groups */

      ntot = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

      /* Get group index */

      if ((ng = GridSearch(ptr, E)) > -1)
        {
          /* Convert index */

          ng = ntot - ng - 1;
          CheckValue(FUNCTION_NAME, "ng2", "", ng, 0, ntot - 1);

          /* Put values */

          ptr = (long)RDB[gcu + GCU_MICRO_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng, wgt*flx, id);

          ptr = (long)RDB[gcu + GCU_MICRO_INV_V];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng, wgt*flx/spd, id);

          /* Check material pointer */

          if (mat > VALID_PTR)
            {
              ptr = (long)RDB[gcu + GCU_MICRO_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*tot, id);

              ptr = (long)RDB[gcu + GCU_MICRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*(capt + fiss), id);

              ptr = (long)RDB[gcu + GCU_MICRO_FISS];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*fiss, id);

              ptr = (long)RDB[gcu + GCU_MICRO_FISSE];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*fissE, id);

              ptr = (long)RDB[gcu + GCU_MICRO_NSF];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*nsf, id);

              /* Flux in fissile materials */

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_FISS_FLX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  AddPrivateRes(ptr + ng, wgt*flx, id);
                }
            }

          /* Eddington factors (Andrew Hall 2016/3/31) */

          if ((long)RDB[DATA_OPTI_EDDINGTON_CALC] == YES)
            {
              ptr = RDB[gcu + GCU_MICRO_FLX_XX];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx*u*u, id);

              ptr = RDB[gcu + GCU_MICRO_FLX_XY];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx*u*v, id);

              ptr = RDB[gcu + GCU_MICRO_FLX_XZ];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx*u*w, id);

              ptr = RDB[gcu + GCU_MICRO_FLX_YY];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx*v*v, id);

              ptr = RDB[gcu + GCU_MICRO_FLX_YZ];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx*v*w, id);

              ptr = RDB[gcu + GCU_MICRO_FLX_ZZ];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx*w*w, id);
            }
        }

      /***********************************************************************/

      /***** Transport-correction ********************************************/

      /* Check material pointer */

      if (mat > VALID_PTR)
        {
          /* Transport correction (NOTE: Täällä tohon vaikutusalaan   */
          /* lisätään total tai f*total, riippuen siitä, käytetäänkö  */
          /* korjausta vai ei. P1-komponentti vähennetään vastaavasti */
          /* scorescattering.c:ssä. */

          if ((ptr = (long)RDB[mat + MATERIAL_PTR_TRANSP_CORR]) < VALID_PTR)
            {
              /* No correction for material, use total */

              xs = tot;
            }
          else if (E < RDB[ptr + TRANSP_CORR_EMIN])
            {
              /* Below limit, use total */

              xs = tot;
            }
          else
            {
              /* Get transport correction */

              f = TransportCorrection(mat, E, id);

              /* Check if nuclide list is given */

              if ((loc0 = (long)RDB[ptr + TRANSP_CORR_PTR_ISO]) < VALID_PTR)
                {
                  /* No list, use corrected total */

                  xs = f*tot;
                }
              else
                {
                  /* Reset total */

                  xs = tot;

                  /* Loop over list */

                  while ((iso = (long)RDB[loc0++]) > VALID_PTR)
                    {
                      /* Pointer to nuclide data */

                      ptr = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

                      /* Pointer to total xs */

                      ptr = (long)RDB[ptr + NUCLIDE_PTR_TOTXS];
                      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);

                      /* Add corrected total, subtract uncorrected total */

                      xs = xs + (f - 1.0)*RDB[iso + COMPOSITION_ADENS]*
                        flx*MicroXS(ptr, E, id);
                    }
                }
            }

          /* Score */

          if ((long)RDB[DATA_PTR_TRC0] > VALID_PTR)
            {
              ptr = (long)RDB[gcu + GCU_MICRO_TRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*xs, id);

              ptr = (long)RDB[gcu + GCU_MICRO_TRC_FLX];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ng, wgt*flx, id);
            }
        }

      /***********************************************************************/

      /***** Beta-eff using Meulekamp's method *******************************/

      /* Check material pointer */

      if (mat > VALID_PTR)
        {
          /* Get delayed neutron group */

          if ((dng = (long)RDB[part + PARTICLE_DN_GROUP]) > 0)
            {
              /* Score fission rate */

              ptr = (long)RDB[gcu + GCU_MEULEKAMP_BETA_EFF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(fiss, wgt, ptr, id, 0);
              AddBuf1D(fiss, wgt, ptr, id, dng);

              /* Get decay constant */

              lambda = RDB[part + PARTICLE_DN_LAMBDA];

              /* Score lambda */

              ptr = (long)RDB[gcu + GCU_MEULEKAMP_LAMBDA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(fiss*lambda, wgt, ptr, id, 0);
              AddBuf1D(fiss*lambda, wgt, ptr, id, dng);
            }

          /* Score total fission rate */

          ptr = (long)RDB[gcu + GCU_MEULEKAMP_TOT_FISS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(fiss, wgt, ptr, id, 0);
        }

      /***********************************************************************/

      /* Next universe */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
        break;
      else
        gcu = NextItem(gcu);
    }

  /***************************************************************************/
}

/*****************************************************************************/
