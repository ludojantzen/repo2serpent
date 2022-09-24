/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorefission.c                                 */
/*                                                                           */
/* Created:       2011/06/11 (JLe)                                           */
/* Last modified: 2018/10/01 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores fission parameters                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreFission:"

/*****************************************************************************/

void ScoreFission(long mat, long rea, double tnu, double dnu, double lambda,
                  long dng, double E0, double E1, double wgt0, double wgt1,
                  long idx0, long idx1, long id)
{
  long ptr, gcu, ntot, ng, nuc, ncol, loc0, fmx, uni;
  long norm;
  double fE;

  /* Check material, reaction and universe pointers */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  ncol = (long)GetPrivateData(ptr, id);

  /* Check neutron type */

  /****************************************************************************/

  /***** Incident neutron *****************************************************/

  if (E0 > 0.0)
    {
      /************************************************************************/

      /***** Fission matrix ***************************************************/

      /* Check matrix indexes */

      if ((idx0 > -1) && (idx1 > -1))
        {
          /* Get pointer to matrix */

          fmx = (long)RDB[DATA_PTR_FMTX];
          CheckPointer(FUNCTION_NAME, "(fmx)", DATA_ARRAY, fmx);

          /* Get pointer to data */

          ptr = (long)RDB[fmx + FMTX_PTR_MTX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Score total */

          AddBuf(tnu, wgt0, ptr, id, -1, 0, idx0, idx1);

          /* Score prompt or delayed */

          if (dng == 0)
            AddBuf(tnu, wgt0, ptr, id, -1, 1, idx0, idx1);
          else
            AddBuf(dnu, wgt0, ptr, id, -1, 2, idx0, idx1);
        }

      /************************************************************************/

      /***** Common parameters for incident neutron ***************************/

      /* Get nuclide pointer */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check ZAI */

      if (((long)RDB[nuc + NUCLIDE_ZAI] == 922330) ||
          ((long)RDB[nuc + NUCLIDE_ZAI] == 922350) ||
          ((long)RDB[nuc + NUCLIDE_ZAI] == 942390) ||
          ((long)RDB[nuc + NUCLIDE_ZAI] == 942410))
        {
          /* Score analog fissile loss rate */

          ptr = (long)RDB[RES_ANA_CONV_RATIO];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt0, ptr, id, 2);
        }

      /* Score particle balance */

      ptr = (long)RDB[RES_N_BALA_LOSS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_LOSS_FISS, 0);
      AddBuf(wgt0, 1.0, ptr, id, -1, BALA_N_LOSS_FISS, 1);

      /* Total and nuclide-wise rates */

      ptr = (long)RDB[RES_ANA_FISS_FRAC];
      CheckPointer(FUNCTION_NAME, "(ptr 4)", DATA_ARRAY, ptr);
      AddBuf(1.0, wgt0, ptr, id, -1, 0, 0);

      if ((long)RDB[nuc + NUCLIDE_ZAI] == 902320)
        AddBuf(1.0, wgt0, ptr, id, -1, 1, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 922330)
        AddBuf(1.0, wgt0, ptr, id, -1, 2, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 922350)
        AddBuf(1.0, wgt0, ptr, id, -1, 3, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 922380)
        AddBuf(1.0, wgt0, ptr, id, -1, 4, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 942390)
        AddBuf(1.0, wgt0, ptr, id, -1, 5, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 942400)
        AddBuf(1.0, wgt0, ptr, id, -1, 6, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 942410)
        AddBuf(1.0, wgt0, ptr, id, -1, 7, 0);

      /* Score analog reaction rate if not implicit mode */

      if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == NO)
        {
          /* Total reaction rate */

          ptr = (long)RDB[RES_TOT_NEUTRON_RR];
          CheckPointer(FUNCTION_NAME, "(ptr 5)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt0, ptr, id, 0);

          /* Total fission rate */

          ptr = (long)RDB[RES_TOT_FISSRATE];
          CheckPointer(FUNCTION_NAME, "(ptr 6)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt0, ptr, id, 0);

          /* Check material burn flag */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            AddBuf1D(1.0, wgt0, ptr, id, 1);
          else
            AddBuf1D(1.0, wgt0, ptr, id, 2);

          /* Total fission neutron production rate */

          ptr = (long)RDB[RES_TOT_NSF];
          CheckPointer(FUNCTION_NAME, "(ptr 6)", DATA_ARRAY, ptr);
          AddBuf1D(tnu, wgt0, ptr, id, 0);

          /* Check material burn flag */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            AddBuf1D(tnu, wgt0, ptr, id, 1);
          else
            AddBuf1D(tnu, wgt0, ptr, id, 2);

          /* Get fission energy deposition (JLe: muutettiin 9.3.2016 / 2.1.25) */

          /*
            fE = RDB[rea + REACTION_Q]*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ;
          */

          fE = FissE(rea, E0, id);
          CheckValue(FUNCTION_NAME, "fE", "", fE/MEV, 0.0, 300.0);

          /* Total fission energy production rate */

          ptr = (long)RDB[RES_TOT_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr 6)", DATA_ARRAY, ptr);
          AddBuf1D(fE, wgt0, ptr, id, 0);

          /* Check material burn flag */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            AddBuf1D(fE, wgt0, ptr, id, 1);
          else
            AddBuf1D(fE, wgt0, ptr, id, 2);

          /* Pointer to normalization */

          if ((norm = (long)RDB[mat + MATERIAL_PTR_NORM]) > VALID_PTR)
            {
              /* Fission rate */

              ptr = (long)RDB[norm + NORM_PTR_FISSRATE];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(1.0, wgt0, ptr, id, 0);

              /* Fission neutron production rate */

              ptr = (long)RDB[norm + NORM_PTR_NSF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(tnu, wgt0, ptr, id, 0);

              /* Fission energy */

              ptr = (long)RDB[norm + NORM_PTR_FISSE];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(fE, wgt0, ptr, id, 0);
            }
        }

      /* Fission rate by delayed neutrons */

      if (dng > 0)
        {
          /* Fission rate for analog delayed neutron k-eff */

          ptr = (long)RDB[RES_ANA_KEFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(tnu, wgt0, ptr, id, 2);
        }
      else
        {
          /* Fission rate for analog prompt neutron k-eff */

          ptr = (long)RDB[RES_ANA_KEFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(tnu, wgt0, ptr, id, 1);
        }

      /* Average lethargy of neutron causing fission */

      ptr = (long)RDB[RES_ANA_ALF];
      CheckPointer(FUNCTION_NAME, "(ptr) ANA_ALF", DATA_ARRAY, ptr);

      AddBuf1D(log(RDB[DATA_NEUTRON_EMAX]/E0), wgt0, ptr, id, 0);

      /* Average energy of neutron causing fission */

      ptr = (long)RDB[RES_ANA_AFGE];
      CheckPointer(FUNCTION_NAME, "(ptr) ANA_AFGE", DATA_ARRAY, ptr);

      AddBuf1D(E0, wgt0, ptr, id, 0);

      /* Absorption and production for six-factor formula */

      if (E0 < 0.625E-6)
        {
          ptr = (long)RDB[RES_SIX_FF_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Total and fuel absorption */

          AddBuf1D(1.0, wgt0, ptr, id, 0);
          AddBuf1D(1.0, wgt0, ptr, id, 1);

          /* Neutron production */

          ptr = (long)RDB[RES_SIX_FF_ETA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(tnu, wgt0, ptr, id, 0);
        }
      else
        {
          /* Fast absorption */

          ptr = (long)RDB[RES_SIX_FF_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt0, ptr, id, 2);
        }

      /* Fast fission factor */

      ptr = (long)RDB[RES_SIX_FF_EPSILON];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(tnu, wgt0, ptr, id, 0);

      /************************************************************************/

      /***** Group constant generation ****************************************/

      /* Check active cycle and corrector step */

      if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
          (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
           ((long)RDB[DATA_B1_BURNUP_CORR] == NO)))
        return;

      /* Check that group constants */

      if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
        return;

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

              if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id)
                  < 0.0)
                {
                  /* Next universe */

                  gcu = NextItem(gcu);

                  /* Cycle loop */

                  continue;
                }
            }

          /*******************************************************************/

          /***** Mora data ***************************************************/

          /* Check pointer to MORA data */

          if ((loc0 = (long)RDB[gcu + GCU_PTR_MORA]) > VALID_PTR)
            {
              /* Get pointer to energy grid */

              ptr = (long)RDB[loc0 + MORA_PTR_EG];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Find group */

              if ((ng = GridSearch(ptr, E0)) > -1)
                {
                  /* Score nubar */

                  if (dng == 0)
                    ptr = (long)RDB[loc0 + MORA_PTR_PNU];
                  else
                    ptr = (long)RDB[loc0 + MORA_PTR_DNU];

                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddBuf1D(tnu, wgt0, ptr, id, ng);
                }
            }

          /*******************************************************************/

          /* Next universe */

          if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
            break;
          else
            gcu = NextItem(gcu);
        }

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Emitted neutron *****************************************************/

  else if (E1 > 0.0)
    {
      /************************************************************************/

      /***** Chi in infinite spectrum *****************************************/

      /* Check active cycle and corrector step */

      if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
          (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
           ((long)RDB[DATA_B1_BURNUP_CORR] == NO)))
        return;

      /* Analog beta-zero and lambda */

      if (dng > 0)
        {
          ptr = (long)RDB[RES_FWD_ANA_BETA_ZERO];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(wgt1, wgt0, ptr, id, 0);
          AddBuf1D(wgt1, wgt0, ptr, id, dng);

          ptr = (long)RDB[RES_FWD_ANA_LAMBDA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(lambda, wgt0, ptr, id, 0);
          AddBuf1D(lambda, wgt0, ptr, id, dng);
        }

      /* Check that group constants are calculated */

      if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
        return;

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

              if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id)
                  < 0.0)
                {
                  /* Next universe */

                  gcu = NextItem(gcu);

                  /* Cycle loop */

                  continue;
                }
            }

          /*******************************************************************/

          /***** MORA data ***************************************************/

          /* Check pointer to MORA data */

          if ((loc0 = (long)RDB[gcu + GCU_PTR_MORA]) > VALID_PTR)
            {
              /* Get pointer to energy grid */

              ptr = (long)RDB[loc0 + MORA_PTR_EG];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Number of groups */

              ntot = (long)RDB[loc0 + MORA_N_EG];

              /* Find group */

              if ((ng = GridSearch(ptr, E1)) > -1)
                {
                  /* Score chi */

                  if (dng == 0)
                    ptr = (long)RDB[loc0 + MORA_PTR_CHIP];
                  else
                    ptr = (long)RDB[loc0 + MORA_PTR_CHID];

                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddBuf1D(wgt1, wgt0, ptr, id, ng);
                  AddBuf1D(wgt1, wgt0, ptr, id, ntot);
                }
            }

          /********************************************************************/

          /***** Micro-group chi **********************************************/

          /* Get pointer to microgroup energy grid */

          ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Number of groups */

          ntot = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

          /* Get few-group index of emitted neutron */

          if ((ng = GridSearch(ptr, E1)) > -1)
            {
              /* Check index */

              CheckValue(FUNCTION_NAME, "ng", "", ng, 0, ntot - 1);

              /* Total chi */

              ptr = (long)RDB[gcu + GCU_MICRO_CHIT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + ntot - ng - 1, wgt0*wgt1, id);

              /* Prompt and delayed */

              if (dng == 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_CHIP];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  AddPrivateRes(ptr + ntot - ng - 1, wgt0*wgt1, id);
                }
              else
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_CHID];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  AddPrivateRes(ptr + ntot - ng - 1, wgt0*wgt1, id);
                }
            }

          /*******************************************************************/

          /* Next universe */

          if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
            break;
          else
            gcu = NextItem(gcu);
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Both energies are zero");
}

/******************************************************************************/
