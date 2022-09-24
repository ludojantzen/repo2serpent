/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearmicrogroupxs.c                            */
/*                                                                           */
/* Created:       2013/03/04 (JLe)                                           */
/* Last modified: 2019/11/06 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Resets micro-group cross sections needed for B1 and P1       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearMicroGroupXS:"

/*****************************************************************************/

void ClearMicroGroupXS()
{
  long ptr, nmg, gcu, adf, ppw, loc0, alb, n, m, i;

  /* Check if group constants are generated */

  if((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check if access to RES2 array is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");

  /* Reduce private results (to reset thread-wise data) */

  ReducePrivateRes();

  /* Avoid compiler warning */

  i = 0;

  /* Reset normalization factor */

  WDB[DATA_MICRO_CALC_NORM] = 0.0;

  /* Reset batch counter */

  WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;

  if ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP)
    {
      /* Get pointer to micro-group structure */

      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr0)", DATA_ARRAY, ptr);

      /* Number of micro- and macro-groups */

      nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Get micro-group reaction rates */

          ptr = (long)RDB[gcu + GCU_MICRO_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_FISS_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr2)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr3)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr5)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_FISS];
          CheckPointer(FUNCTION_NAME, "(ptr4)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_NSF];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_INV_V];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_CHIT];
          CheckPointer(FUNCTION_NAME, "(ptr7)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_CHIP];
          CheckPointer(FUNCTION_NAME, "(ptr8)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_CHID];
          CheckPointer(FUNCTION_NAME, "(ptr9)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*sizeof(double));

          /* Scattering matrixes */

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT0];
          CheckPointer(FUNCTION_NAME, "(ptr10)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP0];
          CheckPointer(FUNCTION_NAME, "(ptr11)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT1];
          CheckPointer(FUNCTION_NAME, "(ptr12)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP1];
          CheckPointer(FUNCTION_NAME, "(ptr13)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT2];
          CheckPointer(FUNCTION_NAME, "(ptr14)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP2];
          CheckPointer(FUNCTION_NAME, "(ptr15)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT3];
          CheckPointer(FUNCTION_NAME, "(ptr16)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP3];
          CheckPointer(FUNCTION_NAME, "(ptr17)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT4];
          CheckPointer(FUNCTION_NAME, "(ptr18)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP4];
          CheckPointer(FUNCTION_NAME, "(ptr19)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT5];
          CheckPointer(FUNCTION_NAME, "(ptr20)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP5];
          CheckPointer(FUNCTION_NAME, "(ptr21)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT6];
          CheckPointer(FUNCTION_NAME, "(ptr22)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP6];
          CheckPointer(FUNCTION_NAME, "(ptr23)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT7];
          CheckPointer(FUNCTION_NAME, "(ptr24)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP7];
          CheckPointer(FUNCTION_NAME, "(ptr25)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmg*nmg*sizeof(double));

          /* Transport corrected xs */

          if ((long)RDB[DATA_PTR_TRC0] > VALID_PTR)
            {
              ptr = (long)RDB[gcu + GCU_MICRO_TRC];
              CheckPointer(FUNCTION_NAME, "(ptr25)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_TRC_FLX];
              CheckPointer(FUNCTION_NAME, "(ptr25)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));
            }

          /* Check poison calculation */

          if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
            {
              ptr = (long)RDB[gcu + GCU_MICRO_I135_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr26)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_XE135_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr27)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_XE135M_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr27)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM147_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM148_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM148M_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM149_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_SM149_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr29)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_I135_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr30)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_XE135_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr31)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_XE135M_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr31)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM147_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM148_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM148M_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_PM149_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_SM149_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr33)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_XE135_MACRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr35)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_XE135M_MACRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr35)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_SM149_MACRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr37)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));
            }

          /* Eddington factors (Andrew Hall 2016/3/31) */

          if ((long)RDB[DATA_OPTI_EDDINGTON_CALC] == YES)
            {
              ptr = (long)RDB[gcu + GCU_MICRO_FLX_XX];
              CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_FLX_XY];
              CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_FLX_XZ];
              CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_FLX_YY];
              CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_FLX_YZ];
              CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_FLX_ZZ];
              CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));
            }

          /* Discontinuity factors  */

          if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) > VALID_PTR)
            {
              /* Get number of vertices and corners */

              n = (long)RDB[adf + ADF_NSURF];
              m = (long)RDB[adf + ADF_NCORN];

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              if (m > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr39)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, nmg*m*sizeof(double));
                }

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_CELL_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr40)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr41)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr42)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr43)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr44)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              if (m > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr45)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, nmg*m*sizeof(double));

                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_OUT_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr46)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, nmg*m*sizeof(double));
                }

              /* Sign moments of discontinuity factors */

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SGN_SURF_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr47)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SGN_SURF_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr48)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SGN_SURF_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr49)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmg*n*sizeof(double));
            }

          /* Pin-power calculation */

          if ((ppw = (long)RDB[gcu + GCU_PTR_PPW]) > VALID_PTR)
            {
              /* Get number of macro energy groups */

              m = (long)RDB[DATA_ERG_FG_NG];

              /* Get number of pins */

              if ((n = (long)RDB[ppw + PPW_NP]) > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_PPW_POW];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, (m + 1)*n*sizeof(double));

                  ptr = (long)RDB[gcu + GCU_MICRO_PPW_XYZ];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, 3*(m + 1)*n*sizeof(double));
                }
            }

          /* Albedo currents */

          if ((alb = (long)RDB[gcu + GCU_PTR_ALB]) > VALID_PTR)
            {
              /* Get number of macro energy groups */

              m = (long)RDB[DATA_ERG_FG_NG];

              /* Get number of surfaces */

              if ((n = (long)RDB[alb + ALB_NSURF]) > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ALB_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, m*n*sizeof(double));

                  ptr = (long)RDB[gcu + GCU_MICRO_ALB_OUT_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, m*m*n*n*sizeof(double));
                }
            }

          /* Micro-depletion */

          if ((loc0 = (long)RDB[gcu + GCU_PTR_MDEP]) > VALID_PTR)
            {
              /* Get number of macro energy groups */

              m = (long)RDB[DATA_ERG_FG_NG];

              /* Get number of reactions */

              if ((n = (long)RDB[loc0 + MDEP_N_REA]) > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_MICRO_DEP_XS];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);
                  memset(&RES2[ptr], 0.0, m*n*sizeof(double));
                }
            }

          /* Cumulative migration method diffusion coefficients */

          ptr = (long)RDB[gcu + GCU_BUF_CMM_CUMU_R2];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, 5*nmg*sizeof(double));

          /* Next universe */

          gcu = NextItem(gcu);
        }
    }

#ifdef DEBUG

  else
    {
      /* Get pointer to micro-group structure */

      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr0)", DATA_ARRAY, ptr);

      /* Number of micro- and macro-groups */

      nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Get micro-group reaction rates */

          ptr = (long)RDB[gcu + GCU_MICRO_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero1");

          ptr = (long)RDB[gcu + GCU_MICRO_FISS_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr2)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero2");

          ptr = (long)RDB[gcu + GCU_MICRO_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr3)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero3");

          ptr = (long)RDB[gcu + GCU_MICRO_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr5)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero4");

          ptr = (long)RDB[gcu + GCU_MICRO_FISS];
          CheckPointer(FUNCTION_NAME, "(ptr4)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero5");

          ptr = (long)RDB[gcu + GCU_MICRO_NSF];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero6");

          ptr = (long)RDB[gcu + GCU_MICRO_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero7");

          ptr = (long)RDB[gcu + GCU_MICRO_INV_V];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero8");

          ptr = (long)RDB[gcu + GCU_MICRO_CHIT];
          CheckPointer(FUNCTION_NAME, "(ptr7)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero9");

          ptr = (long)RDB[gcu + GCU_MICRO_CHIP];
          CheckPointer(FUNCTION_NAME, "(ptr8)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero10");

          ptr = (long)RDB[gcu + GCU_MICRO_CHID];
          CheckPointer(FUNCTION_NAME, "(ptr9)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero11");

          /* Scattering matrixes */

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT0];
          CheckPointer(FUNCTION_NAME, "(ptr10)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero12");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP0];
          CheckPointer(FUNCTION_NAME, "(ptr11)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero13");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT1];
          CheckPointer(FUNCTION_NAME, "(ptr12)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero14");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP1];
          CheckPointer(FUNCTION_NAME, "(ptr13)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero15");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT2];
          CheckPointer(FUNCTION_NAME, "(ptr14)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero16");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP2];
          CheckPointer(FUNCTION_NAME, "(ptr15)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero17");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT3];
          CheckPointer(FUNCTION_NAME, "(ptr16)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero18");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP3];
          CheckPointer(FUNCTION_NAME, "(ptr17)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero19");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT4];
          CheckPointer(FUNCTION_NAME, "(ptr18)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero20");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP4];
          CheckPointer(FUNCTION_NAME, "(ptr19)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero21");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT5];
          CheckPointer(FUNCTION_NAME, "(ptr20)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero22");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP5];
          CheckPointer(FUNCTION_NAME, "(ptr21)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero23");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT6];
          CheckPointer(FUNCTION_NAME, "(ptr22)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero24");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP6];
          CheckPointer(FUNCTION_NAME, "(ptr23)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero25");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATT7];
          CheckPointer(FUNCTION_NAME, "(ptr24)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero26");

          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP7];
          CheckPointer(FUNCTION_NAME, "(ptr25)", RES2_ARRAY, ptr);

          for (i = 0; i < nmg*nmg; i++)
            if (RES2[ptr++] != 0.0)
              Die(FUNCTION_NAME, "Nonzero27");

          /* Check poison calculation */

          if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
            {
              ptr = (long)RDB[gcu + GCU_MICRO_I135_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr26)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero28");

              ptr = (long)RDB[gcu + GCU_MICRO_XE135_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr27)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero29");

              ptr = (long)RDB[gcu + GCU_MICRO_XE135M_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr27)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero29");

              ptr = (long)RDB[gcu + GCU_MICRO_PM147_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero30");

              ptr = (long)RDB[gcu + GCU_MICRO_PM148_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero30");

              ptr = (long)RDB[gcu + GCU_MICRO_PM148M_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero30");

              ptr = (long)RDB[gcu + GCU_MICRO_PM149_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr28)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero30");

              ptr = (long)RDB[gcu + GCU_MICRO_SM149_YIELD];
              CheckPointer(FUNCTION_NAME, "(ptr29)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero31");

              ptr = (long)RDB[gcu + GCU_MICRO_I135_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr30)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero32");

              ptr = (long)RDB[gcu + GCU_MICRO_XE135_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr31)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero33");

              ptr = (long)RDB[gcu + GCU_MICRO_XE135M_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr31)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero33");

              ptr = (long)RDB[gcu + GCU_MICRO_PM147_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero34");

              ptr = (long)RDB[gcu + GCU_MICRO_PM148_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero34");

              ptr = (long)RDB[gcu + GCU_MICRO_PM148M_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero34");

              ptr = (long)RDB[gcu + GCU_MICRO_PM149_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr32)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero34");

              ptr = (long)RDB[gcu + GCU_MICRO_SM149_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr33)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero35");

              ptr = (long)RDB[gcu + GCU_MICRO_XE135_MACRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr35)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero37");

              ptr = (long)RDB[gcu + GCU_MICRO_XE135M_MACRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr35)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero37");

              ptr = (long)RDB[gcu + GCU_MICRO_SM149_MACRO_ABS];
              CheckPointer(FUNCTION_NAME, "(ptr37)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero39");
            }

          /* Discontinuity factors  */

          if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) > VALID_PTR)
            {
              /* Get number of vertices and corners */

              n = (long)RDB[adf + ADF_NSURF];
              m = (long)RDB[adf + ADF_NCORN];

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg*n; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero40");

              if (m > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr39)", RES2_ARRAY, ptr);

                  for (i = 0; i < nmg*m; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero41");
                }

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_CELL_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr40)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero42");

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr41)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg*n; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero43");

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr42)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg*n; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero44");

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr43)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg*n; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero45");

              ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr44)", RES2_ARRAY, ptr);

              for (i = 0; i < nmg*n; i++)
                if (RES2[ptr++] != 0.0)
                  Die(FUNCTION_NAME, "Nonzero46");

              if (m > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr45)", RES2_ARRAY, ptr);

                  for (i = 0; i < nmg*m; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero47");

                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_OUT_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr46)", RES2_ARRAY, ptr);

                  for (i = 0; i < nmg*m; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero48");
                }
            }

          /* Pin-power calculation */

          if ((ppw = (long)RDB[gcu + GCU_PTR_PPW]) > VALID_PTR)
            {
              /* Get number of macro energy groups */

              m = (long)RDB[DATA_ERG_FG_NG];

              /* Get number of pins */

              if ((n = (long)RDB[ppw + PPW_NP]) > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_PPW_POW];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);

                  for (i = 0; i < (m + 1)*n; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero49");

                  ptr = (long)RDB[gcu + GCU_MICRO_PPW_XYZ];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);

                  for (i = 0; i < 3*(m + 1)*n; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero50");
                }
            }

          /* Albedo currents */

          if ((alb = (long)RDB[gcu + GCU_PTR_ALB]) > VALID_PTR)
            {
              /* Get number of macro energy groups */

              m = (long)RDB[DATA_ERG_FG_NG];

              /* Get number of surfaces */

              if ((n = (long)RDB[alb + ALB_NSURF]) > 0)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ALB_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);

                  for (i = 0; i < m*n; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero51");

                  ptr = (long)RDB[gcu + GCU_MICRO_ALB_OUT_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr38)", RES2_ARRAY, ptr);

                  for (i = 0; i < m*m*n*n; i++)
                    if (RES2[ptr++] != 0.0)
                      Die(FUNCTION_NAME, "Nonzero52");
                }
            }

          /* Next universe */

          gcu = NextItem(gcu);
        }
    }

#endif
}

/*****************************************************************************/
