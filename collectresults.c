/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectresults.c                               */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2019/08/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Collects results after cycle or batch                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectResults:"

/*****************************************************************************/

void CollectResults()
{
  long ptr, mat, i, n, m, ng, n0, n1, n2, nz, gcu, pbd, nuc, rea, loc0, loc1;
  long uni, nst, nr, k, pts, nc, np, nseg, ma0, ms0, g0, g1;
  long nmus, nmua, ma1, ms1, ng0, ng1, loc2;
  double tot, nuxn, capt, scatt, fiss, leak, val, flx, tots, gent, keff;
  double norm, sum, fE, div, fmass, wgt0, wgt1, beta, invv, nsf;
  double kcp, cut;

  /***************************************************************************/

  /***** Get parallel data ***************************************************/

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Collect MPI parallel data */

  CollectBuf();

  /* Avoid compiler warning */

  sum = 0.0;
  tots = 0.0;

  /***************************************************************************/

  /***** Collision and reaction sampling *************************************/

  /* Loop over particle types */

  for (n = 0; n < 2; n++)
    {
      /* Source sampling efficiency */

      ptr = (long)RDB[RES_SRC_SAMPLING_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      if ((div = BufVal(ptr, n + 2)) > 0.0)
        {
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);
        }

      /* Average number of tracking loops */

      ptr = (long)RDB[RES_AVG_TRACK_LOOPS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);

      /* Fraction of failed loops */

      if ((div = BufN(ptr, n)) > 0.0)
        {
          val = BufVal(ptr, n + 2);
          AddStat(val/div, ptr, n + 2);
        }

      /* Mean ww source splitting */

      ptr = (long)RDB[RES_SRC_WW_SPLIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);

      /* Mean ww source sampling efficiency */

      ptr = (long)RDB[RES_SRC_WW_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      if ((div = BufVal(ptr, n, 1)) > 0.0)
        {
          val = BufVal(ptr, n, 0);
          AddStat(1.0 - val/div, ptr, n, 0);
        }

      /* Mean source weight */

      ptr = (long)RDB[RES_SRC_MEAN_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);

      /* Mean ww weight balance */

      ptr = (long)RDB[RES_WW_BALA_ROULETTE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);

      ptr = (long)RDB[RES_WW_BALA_SPLIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);

      /* Get total number of sampled tracks (ST + DT) */

      ptr = (long)RDB[RES_ST_TRACK_FRAC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, n);

      ptr = (long)RDB[RES_DT_TRACK_FRAC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div + BufVal(ptr, n);

      /* Calculate fractions */

      if (div > 0.0)
        {
          ptr = (long)RDB[RES_ST_TRACK_FRAC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);

          ptr = (long)RDB[RES_DT_TRACK_FRAC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);
        }

      /* Efficiency of delta-tracking */

      ptr = (long)RDB[RES_DT_TRACK_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, n) + BufVal(ptr, n + 2);

      if (div > 0.0)
        {
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);
        }

      /* Total collision efficiency */

      ptr = (long)RDB[RES_TOT_COL_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, n) + BufVal(ptr, n + 2);

      if (div > 0.0)
        {
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);
        }

      /* Reaction sampling efficiency */

      ptr = (long)RDB[RES_REA_SAMPLING_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, n) + BufVal(ptr, n + 2);

      if (div > 0.0)
        {
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);
        }

      /* Failure rate */

      ptr = (long)RDB[RES_REA_SAMPLING_FAIL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, n);

      if (div > 0.0)
        AddStat(val/div, ptr, n);

      /* Efficiency of interface collision rejection */

      ptr = (long)RDB[RES_IFC_COL_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);

      if (val > 0.0)
        AddStat(val, ptr, n);

      /* Average number of tracks, collisions and surface crossings per */
      /* history (NOTE: noi on per starting particle, eli varsinkin     */
      /* fotoneilla ja neutronien external source -laskussa pitää ottaa */
      /* huomioon se, että sekundääriset lasketaan summaan mukaan).     */

      if ((div = RDB[DATA_CYCLE_BATCH_SIZE]) > 0.0)
        {
          ptr = (long)RDB[RES_AVG_TRACKS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);

          ptr = (long)RDB[RES_AVG_SURF_CROSS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);

          ptr = (long)RDB[RES_AVG_REAL_COL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);

          ptr = (long)RDB[RES_AVG_VIRT_COL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/div, ptr, n);
        }

      /* Minimum macroscopic XS */

      ptr = (long)RDB[RES_MIN_MACROXS];
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);
    }

  /* TMS sampling efficiency */

  ptr = (long)RDB[RES_TMS_SAMPLING_EFF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  if ((div = BufVal(ptr, 1)) > 0)
    {
      val = BufVal(ptr, 0);
      AddStat(val/div, ptr, 0);
    }

  /* TMS fail statistics */

  ptr = (long)RDB[RES_TMS_FAIL_STAT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  if ((div = BufVal(ptr, 0)) > 0)
    {
      val = BufVal(ptr, 1);
      AddStat(val/div, ptr, 1);

      val = BufVal(ptr, 2);
      AddStat(val/div, ptr, 2);

      val = BufVal(ptr, 3);
      AddStat(val/div, ptr, 3);
    }

  /* STL ray test failures */

  if ((long)RDB[DATA_PTR_STL0] > VALID_PTR)
    {
      /* Get pointer */

      ptr = (long)RDB[RES_STL_RAY_TEST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get total */

      if ((div = BufVal(ptr, 0)) > 0.0)
        {
          val = BufVal(ptr, 1);
          AddStat(val/div, ptr, 0);

          val = BufVal(ptr, 2);
          AddStat(val/div, ptr, 1);

          val = BufVal(ptr, 3);
          AddStat(val/div, ptr, 2);

          val = BufVal(ptr, 4);
          AddStat(val/div, ptr, 3);

          val = BufVal(ptr, 5);
          AddStat(val/div, ptr, 4);
        }
    }

  /* Response matrix calculation */

  if ((long)RDB[DATA_RMTX_CALC] == YES)
    {
      /* Fail rate */

      ptr = (long)RDB[RES_RMX_CURR_SEARCH_FAIL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 1)) > 0.0)
        {
          val = BufVal(ptr, 0);
          AddStat(val/div, ptr, 0);
        }
    }

  /* Photon sampling */

  if ((long)RDB[DATA_PHOTON_PRODUCTION] != NO)
    {
      /* Fail rate */

      ptr = (long)RDB[RES_PHOTON_SAMPLING_FAIL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 1)) > 0.0)
        {
          val = BufVal(ptr, 0);
          AddStat(val/div, ptr, 0);
        }
    }

  /* Particle balance for domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Particles sent to limbo */

      ptr = (long)RDB[RES_DD_BALA_IN];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      for (n = 0; n < mpitasks; n++)
        for (m = 0; m < 2; m++)
          {
            val = BufVal(ptr, n, m);
            AddStat(val, ptr, n, m);
          }

      /* Particles retrived from limbo */

      ptr = (long)RDB[RES_DD_BALA_OUT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      for (n = 0; n < mpitasks; n++)
        for (m = 0; m < 2; m++)
          {
            val = BufVal(ptr, n, m);
            AddStat(val, ptr, n, m);
          }
    }

  /***************************************************************************/

  /***** Particle balance for photons ****************************************/

  /* Check transport mode */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Particle source */

      ptr = (long)RDB[RES_G_BALA_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over particles and weight */

      for (m = 0; m < 3; m++)
        {
          /* Loop over terms */

          tot = 0.0;
          for (n = 0; n < 6; n++)
            {
              val = BufVal(ptr, n, m);
              AddStat(val, ptr, n, m);
              tot = tot + val;
            }

          /* Total */

          AddStat(tot, ptr, BALA_G_SRC_TOT, m);
        }

      /* Particle loss */

      ptr = (long)RDB[RES_G_BALA_LOSS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over particles and weight */

      for (m = 0; m < 2; m++)
        {
          /* Loop over terms */

          tot = 0.0;
          for (n = 0; n < 4; n++)
            {
              val = BufVal(ptr, n, m);
              AddStat(val, ptr, n, m);
              tot = tot + val;
            }

          /* Total */

          AddStat(tot, ptr, BALA_G_LOSS_TOT, m);
        }
    }

  /***************************************************************************/

  /***** Normalized total reaction rates for photons *************************/

  /* Check transport mode */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Get normalization factor for photon source */

      norm = NormCoef(PARTICLE_TYPE_GAMMA);
      CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);

      /* Store value */

      ptr = (long)RDB[RES_NORM_COEF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(norm, ptr, 1);

      /* Reset total loss rate */

      tot = 0.0;

      /* Total leak rate */

      ptr = (long)RDB[RES_TOT_PHOTON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      tot = tot + val;
      AddStat(norm*val, ptr, 0);

      /* Total energy cut-off rate */

      ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      tot = tot + val;
      AddStat(norm*val, ptr, 0);

      /* Photoelectric capture rate */

      ptr = (long)RDB[RES_PHOTOELE_CAPT_RATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      tot = tot + val;
      AddStat(norm*val, ptr, 0);

      /* Pair production capture rate */

      ptr = (long)RDB[RES_PAIRPROD_CAPT_RATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      tot = tot + val;
      AddStat(norm*val, ptr, 0);

      /* Total loss rate */

      ptr = (long)RDB[RES_TOT_PHOTON_LOSSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(norm*tot, ptr, 0);

      /* Total reaction rate */

      ptr = (long)RDB[RES_TOT_PHOTON_RR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      AddStat(norm*val, ptr, 0);

      /* Total source rate */

      ptr = (long)RDB[RES_TOT_PHOTON_SRCRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      AddStat(norm*val, ptr, 0);

      /* Total flux */

      ptr = (long)RDB[RES_TOT_PHOTON_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      AddStat(norm*val, ptr, 0);

      /* Total heating rate */

      ptr = (long)RDB[RES_TOT_PHOTON_HEATRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);
      AddStat(norm*val, ptr, 0);

      /* Analog reaction rate estimates */

      if ((long)RDB[DATA_ANA_RR_PCALC] != ARR_MODE_NONE)
        {
          /* Loop over nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
                {
                  /* Loop over reactions */

                  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
                  while (rea > VALID_PTR)
                    {
                      /* Get pointer */

                      if ((ptr = (long)RDB[rea + REACTION_PTR_ANA_RATE]) >
                          VALID_PTR)
                        {
                          /* Add to statistics */

                          val = BufVal(ptr, 0);
                          AddStat(val*norm, ptr, 0);
                        }

                      /* Next reaction */

                      rea = NextItem(rea);
                    }
                }

              /* Next nuclide */

              nuc = NextItem(nuc);
            }
        }

      /* Lifetime */

      ptr = (long)RDB[RES_ANA_PHOTON_LIFETIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, 0);
      AddStat(val, ptr, 0);

      /* Sampled decay source */

      if ((long)RDB[DATA_USE_DECAY_SRC] == YES)
        {
          /* Loop over materials */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Sampled source */

              if ((ptr = (long)RDB[mat + MATERIAL_SAMPLED_DECAY_SRC])
                  > VALID_PTR)
                {
                  val = BufVal(ptr, 1);
                  AddStat(norm*val, ptr, 1);
                }

              /* Next material */

              mat = NextItem(mat);
            }
        }
    }

  /* Check neutron transport mode */

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
    return;

  /* Get normalization factor for neutron source */

  norm = NormCoef(PARTICLE_TYPE_NEUTRON);
  CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);

  /* Store value */

  ptr = (long)RDB[RES_NORM_COEF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddStat(norm, ptr, 0);

  /* Sampled decay source */

  if ((long)RDB[DATA_USE_DECAY_SRC] == YES)
    {
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Sampled source */

          if ((ptr = (long)RDB[mat + MATERIAL_SAMPLED_DECAY_SRC]) > VALID_PTR)
            {
              val = BufVal(ptr, 0);
              AddStat(norm*val, ptr, 0);
            }

          /* Next material */

          mat = NextItem(mat);
        }
    }

  /***************************************************************************/

  /***** Particle balance for neutrons ***************************************/

  /* Particle source */

  ptr = (long)RDB[RES_N_BALA_SRC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over particles and weight */

  for (m = 0; m < 2; m++)
    {
      /* Loop over terms */

      tot = 0.0;
      for (n = 0; n < 4; n++)
        {
          val = BufVal(ptr, n, m);
          AddStat(val, ptr, n, m);
          tot = tot + val;
        }

      /* Total */

      AddStat(tot, ptr, BALA_N_SRC_TOT, m);
    }

  /* Particle loss */

  ptr = (long)RDB[RES_N_BALA_LOSS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over particles and weight */

  for (m = 0; m < 2; m++)
    {
      /* Loop over terms */

      tot = 0.0;
      for (n = 0; n < 5; n++)
        {
          val = BufVal(ptr, n, m);
          AddStat(val, ptr, n, m);
          tot = tot + val;
        }

      /* Total */

      AddStat(tot, ptr, BALA_N_LOSS_TOT, m);
    }

  /***************************************************************************/

  /***** Common integral parameters for neutrons *****************************/

  /* Neutron production term */

  ptr = (long)RDB[RES_TOT_NSF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  nsf = BufVal(ptr, 0);

  /* Fission term */

  ptr = (long)RDB[RES_TOT_FISSRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  fiss = BufVal(ptr, 0);

  /* Total source weights */

  ptr = (long)RDB[RES_INI_SRC_WGT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  wgt0 = BufVal(ptr, 0);

  ptr = (long)RDB[RES_NEW_SRC_WGT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  wgt1 = BufVal(ptr, 0);

  /* Check initial weight */

  CheckValue(FUNCTION_NAME, "wgt0", "", wgt0, 0.0, INFTY);

  /* Check mode */

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN))
    {
      /* Collision estimate of k-eff */

      ptr = (long)RDB[RES_COL_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(1.0 - wgt0/(wgt0 + nsf), ptr, 0);

      /* Analog estimate of k-eff */

      ptr = (long)RDB[RES_ANA_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(1.0 - wgt0/(wgt0 + wgt1), ptr, 0);

      /* Source multiplication */

      ptr = (long)RDB[RES_SRC_MULT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat((wgt1 + wgt0)/wgt0, ptr, 0);

      /* Mean number of generations */

      ptr = (long)RDB[RES_MEAN_NGEN];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, 0);
      AddStat(val, ptr, 0);

      /* Loop over generations */

      for (i = 0; i < MAX_EXT_K_GEN; i++)
        if (wgt0 > 0.0)
          {
            /* Get new source weight */

            ptr = (long)RDB[RES_NEW_SRC_WGT];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            wgt1 = BufVal(ptr, i + 1);

            /* Score k-eff */

            ptr = (long)RDB[RES_EXT_K];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            AddStat(wgt1/wgt0, ptr, i);

            /* Update weight */

            wgt0 = wgt1;
          }

      /* Mean prompt chain length */

      ptr = (long)RDB[RES_PROMPT_CHAIN_LENGTH];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, 0);
      AddStat(val, ptr, 0);

      /* Prompt generation fractions */

      if ((ptr = (long)RDB[RES_PROMPT_GEN_CUMU]) > VALID_PTR)
        if ((div = BufVal(ptr, 0)) > 0.0)
          {
            /* Cumulative fractions */

            for (n = 0; n < (long)RDB[DATA_MAX_PROMPT_CHAIN_LENGTH]; n++)
              {
                val = BufVal(ptr, n);
                AddStat(val/div, ptr, n);
              }

            /* Times */

            ptr = (long)RDB[RES_PROMPT_GEN_TIMES];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

            for (n = 0; n < (long)RDB[DATA_MAX_PROMPT_CHAIN_LENGTH]; n++)
              {
                val = BufMean(ptr, n);
                AddStat(val, ptr, n);
              }
          }

      /* Prompt generation populations */

      if ((ptr = (long)RDB[RES_PROMPT_GEN_POP]) > VALID_PTR)
        if ((div = BufVal(ptr, 0)) > 0.0)
          for (n = 0; n < (long)RDB[DATA_MAX_PROMPT_CHAIN_LENGTH]; n++)
            {
              val = BufVal(ptr, n);
              AddStat(val/div, ptr, n);
            }
    }
  else
    {
      /* Analog estimate of k-eff */

      ptr = (long)RDB[RES_ANA_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Total */

      val = BufMean(ptr, 0);
      AddStat(val, ptr, 0);

      /* Prompt */

      val = BufVal(ptr, 1)/RDB[DATA_CRIT_POP]/((double)mpitasks);
      AddStat(val, ptr, 1);

      /* Delayed */

      val = BufVal(ptr, 2)/RDB[DATA_CRIT_POP]/((double)mpitasks);
      AddStat(val, ptr, 2);

      /* Check Wielandt shift */

      if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
        {
          /* Backtransformed collision estimate of k-eff */

          kcp = (RDB[DATA_WIELANDT_KEFF]*nsf/wgt0)/
            (RDB[DATA_WIELANDT_KEFF] + nsf/wgt0);

          /* Store backtransformed keff */

          ptr = (long)RDB[RES_COL_KEFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(kcp, ptr, 0);

          /* K and P */

          ptr = (long)RDB[RES_WIELANDT_K];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);
          AddStat(val, ptr, 0);

          ptr = (long)RDB[RES_WIELANDT_P];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);
          AddStat(val, ptr, 0);
        }
      else
        {
          /* Collision estimate of k-eff */

          ptr = (long)RDB[RES_COL_KEFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(nsf/wgt0, ptr, 0);
        }
    }

  /* Total capture rate */

  ptr = (long)RDB[RES_TOT_CAPTRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  capt = BufVal(ptr, 0);

  /* Scattering production rate */

  ptr = (long)RDB[RES_TOT_INLPRODRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  nuxn = BufVal(ptr, 0);

  /* Leak term */

  ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  leak = BufVal(ptr, 0);

  ptr = (long)RDB[RES_ALB_NEUTRON_LEAKRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  leak = leak + BufVal(ptr, 0);

  /* Implicit estimate of k-eff */

  ptr = (long)RDB[RES_IMP_KEFF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  keff = nsf/(capt + fiss - nuxn + leak + ZERO);
  AddStat(keff, ptr, 0);

  /* Implicit estimate of k-inf */

  ptr = (long)RDB[RES_IMP_KINF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddStat(nsf/(capt + fiss - nuxn + ZERO), ptr, 0);

  /* Albedo */

  ptr = (long)RDB[RES_GEOM_ALBEDO];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddStat(RDB[DATA_GEOM_ALBEDO1], ptr, 0);
  AddStat(RDB[DATA_GEOM_ALBEDO2], ptr, 1);
  AddStat(RDB[DATA_GEOM_ALBEDO3], ptr, 2);

  /* Total flux */

  ptr = (long)RDB[RES_TOT_NEUTRON_FLUX];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  div = BufVal(ptr, 0);

  /* 1/v */

  if (div > 0.0)
    {
      ptr = (long)RDB[RES_TOT_RECIPVEL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      invv = BufVal(ptr, 0);
      AddStat(invv/div, ptr, 0);
    }

  /* Total photon production rate */

  ptr = (long)RDB[RES_TOT_PHOTON_PRODRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  val = BufVal(ptr, 0);
  AddStat(norm*val, ptr, 0);

  val = BufVal(ptr, 1);
  AddStat(norm*val, ptr, 1);

  /***************************************************************************/

  /***** Average energy causing fission **************************************/

  /* Average lethargy of neutron causing fission (analog estimator) */

  ptr = (long)RDB[RES_ANA_ALF];
  CheckPointer(FUNCTION_NAME, "(ptr) ANA_ALF", DATA_ARRAY, ptr);

  val = BufMean(ptr, 0);
  AddStat(val, ptr, 0);

  /* Energy corresponding to average lethargy of neutron causing fission */

  ptr = (long)RDB[RES_ANA_EALF];
  CheckPointer(FUNCTION_NAME, "(ptr) ANA_EALF", DATA_ARRAY, ptr);

  if (val > 0.0)
    AddStat(RDB[DATA_NEUTRON_EMAX]*exp(-val), ptr, 0);

  /* Average energy of neutron causing fission (analog estimator) */

  ptr = (long)RDB[RES_ANA_AFGE];
  CheckPointer(FUNCTION_NAME, "(ptr) ANA_AFGE", DATA_ARRAY, ptr);
  val = BufMean(ptr, 0);

  AddStat(val, ptr, 0);

  /* Check fission rate */

  if (fiss > 0.0)
    {
      /* Average lethargy of neutron causing fission (implicit estimator) */

      ptr = (long)RDB[RES_IMP_ALF];
      CheckPointer(FUNCTION_NAME, "(ptr) IMP_ALF", DATA_ARRAY, ptr);

      val = BufVal(ptr, 0);
      AddStat(val/fiss, ptr, 0);

      /* Energy corresponding to average lethargy of neutron causing fission */

      ptr = (long)RDB[RES_IMP_EALF];
      CheckPointer(FUNCTION_NAME, "(ptr) IMP_EALF", DATA_ARRAY, ptr);

      AddStat(RDB[DATA_NEUTRON_EMAX]*exp(-val/fiss), ptr, 0);

      /* Average energy of neutron causing fission (implicit estimator) */

      ptr = (long)RDB[RES_IMP_AFGE];
      CheckPointer(FUNCTION_NAME, "(ptr) IMP_AFGE", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);

      AddStat(val/fiss, ptr, 0);
    }

  /***************************************************************************/

  /***** Forward-weighted time constants *************************************/

  /* Implicit lifetime estimators */

  if (nsf > 0.0)
    {
      /* Analog beta zero */

      ptr = (long)RDB[RES_FWD_ANA_BETA_ZERO];
      CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

      /* Loop over precursor groups */

      for (n = 0; n < (long)RDB[DATA_PRECURSOR_GROUPS] + 1; n++)
        {
          val = BufVal(ptr, n);
          AddStat(val/nsf, ptr, n);
        }

      /* Analog decay constant */

      ptr = (long)RDB[RES_FWD_ANA_LAMBDA];
      CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

      /* Loop over precursor groups */

      for (n = 0; n < (long)RDB[DATA_PRECURSOR_GROUPS] + 1; n++)
        {
          val = BufMean(ptr, n);
          AddStat(val, ptr, n);
        }
    }

  /***************************************************************************/

  /***** Adjoint-weighted time constants *************************************/

  /* Get total fission rate */

  ptr = (long)RDB[RES_TOT_FISSRATE];
  CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
  fiss = BufVal(ptr, 0);

  /* Get analog k-eff */

  keff = RDB[DATA_CYCLE_KEFF];

  /* Check */

  if (fiss > 0.0)
    {
      /* Meulekamp beta-eff and lambda */

      for (n = 0; n < (long)RDB[DATA_PRECURSOR_GROUPS] + 1; n++)
        {
          /* Beta-eff */

          ptr = (long)RDB[RES_ADJ_MEULEKAMP_BETA_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);
          val = BufVal(ptr, n);
          AddStat(val/fiss, ptr, n);

          /* Lambda */

          if (val > 0.0)
            {
              ptr = (long)RDB[RES_ADJ_MEULEKAMP_LAMBDA];
              CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n)/val;
              AddStat(val, ptr, n);
            }
        }

      /* IFP estimators */

      if ((np = (long)RDB[DATA_IFP_CHAIN_LENGTH]) > 0)
        {
          /* Loop over chain */

          for (m = 0; m < np; m++)
            {
              /* Avoid compiler warning */

              gent = 0.0;
              beta = 0.0;

              /* Loop over types (t/p/d) */

              for (n = 0; n < 3; n++)
                {
                  /* Lifetime */

                  ptr = (long)RDB[RES_ADJ_IFP_LIFETIME];
                  CheckPointer(FUNCTION_NAME, "(ptr8)", DATA_ARRAY, ptr);
                  val = BufMean(ptr, n, m);
                  AddStat(val, ptr, n, m);

                  if (m == 0)
                    {
                      ptr = (long)RDB[RES_ADJ_NAUCHI_LIFETIME];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddStat(val, ptr, n);
                    }

                  /* Generation time */

                  ptr = (long)RDB[RES_ADJ_IFP_GEN_TIME];
                  CheckPointer(FUNCTION_NAME, "(ptr7)", DATA_ARRAY, ptr);
                  AddStat(val/keff, ptr, n, m);

                  if (m == 0)
                    {
                      ptr = (long)RDB[RES_ADJ_NAUCHI_GEN_TIME];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddStat(val/keff, ptr, n);
                    }

                  /* Remember prompt value */

                  if (n == 1)
                    gent = val/keff;
                }

              /* Loop over precursor groups */

              for (n = 0; n < (long)RDB[DATA_PRECURSOR_GROUPS] + 1; n++)
                {
                  /* Analog beta-eff */

                  ptr = (long)RDB[RES_ADJ_IFP_ANA_BETA_EFF];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  val = BufVal(ptr, n, m);
                  div = BufVal(ptr, 9, m);

                  if (div > 0.0)
                    {
                      /* Score */

                      AddStat(val/div, ptr, n, m);

                      /* Score NFP value */

                      if (m == 0)
                        {
                          ptr = (long)RDB[RES_ADJ_NAUCHI_BETA_EFF];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                          AddStat(val/div, ptr, n);
                        }
                    }

                  /* Analog Lambda */

                  if (val > 0.0)
                    {
                      ptr = (long)RDB[RES_ADJ_IFP_ANA_LAMBDA];
                      CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);
                      val = BufVal(ptr, n, m)/val;
                      AddStat(val, ptr, n, m);

                      /* Score NFP value */

                      if (m == 0)
                        {
                          ptr = (long)RDB[RES_ADJ_NAUCHI_LAMBDA];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                          AddStat(val, ptr, n);
                        }
                    }

                  /* Implicit beta-eff */

                  ptr = (long)RDB[RES_ADJ_IFP_IMP_BETA_EFF];
                  CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);
                  val = BufVal(ptr, n, m);
                  AddStat(val/fiss, ptr, n, m);

                  /* Remember total value */

                  if (n == 0)
                    beta = val/fiss;

                  /* Implicit Lambda */

                  if (val > 0.0)
                    {
                      ptr = (long)RDB[RES_ADJ_IFP_IMP_LAMBDA];
                      CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);
                      val = BufVal(ptr, n, m)/val;
                      AddStat(val, ptr, n, m);
                    }
                }

              /* Rossi alpha */

              if (gent > 0.0)
                {
                  ptr = (long)RDB[RES_ADJ_IFP_ROSSI_ALPHA];
                  CheckPointer(FUNCTION_NAME, "(ptr10)", DATA_ARRAY, ptr);
                  AddStat(-beta/gent, ptr, m);
                }
            }
        }

      /* PERT estimators */

      if ((np = (long)RDB[DATA_IFP_CHAIN_LENGTH]) > 0)
        {
          /* Avoid compiler warning */

          gent = 0.0;
          beta = 0.0;

          /* Loop over types (t/p/d) */

          for (n = 0; n < 3; n++)
            {
              /* Lifetime */

              ptr = (long)RDB[RES_ADJ_PERT_LIFETIME];
              CheckPointer(FUNCTION_NAME, "(ptr11)", DATA_ARRAY, ptr);
              val = BufMean(ptr, n);
              AddStat(val, ptr, n);

              /* Generation time */

              ptr = (long)RDB[RES_ADJ_PERT_GEN_TIME];
              CheckPointer(FUNCTION_NAME, "(ptr12)", DATA_ARRAY, ptr);
              AddStat(val/keff, ptr, n);

              /* Remember prompt value (not calculated, use total) */

              if (n == 0)
                gent = val/keff;
            }

          /* Beta-eff */

          ptr = (long)RDB[RES_ADJ_PERT_BETA_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr13)", DATA_ARRAY, ptr);

          /* Loop over precursor groups */

          for (n = 0; n < (long)RDB[DATA_PRECURSOR_GROUPS] + 1; n++)
            {
              val = BufMean(ptr, n);
              AddStat(val, ptr, n);

              /* Remember total value */

              if (n == 0)
                beta = val;
            }

          /* Rossi alpha */

          if (gent > 0.0)
            {
              ptr = (long)RDB[RES_ADJ_PERT_ROSSI_ALPHA];
              CheckPointer(FUNCTION_NAME, "(ptr14)", DATA_ARRAY, ptr);
              AddStat(-beta/gent, ptr, 0);
            }
        }
    }

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Get total fission rate */

      ptr = (long)RDB[gcu + GCU_MEULEKAMP_TOT_FISS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fiss = BufVal(ptr, 0);

      /* Check */

      if (fiss > 0.0)
        {
          /* Meulekamp beta-eff and lambda */

          for (n = 0; n < (long)RDB[DATA_PRECURSOR_GROUPS] + 1; n++)
            {
              /* Beta-eff */

              ptr = (long)RDB[gcu + GCU_MEULEKAMP_BETA_EFF];
              CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n);
              AddStat(val/fiss, ptr, n);

              /* Lambda */

              if (val > 0.0)
                {
                  ptr = (long)RDB[gcu + GCU_MEULEKAMP_LAMBDA];
                  CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);
                  val = BufVal(ptr, n)/val;
                  AddStat(val, ptr, n);
                }
            }
        }

      /* Next universe */

      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Misc. time constants ************************************************/

  /* Loop over total, prompt and delayed values */

  for (n = 0; n < 3; n++)
    {
      /* Slowing-down time */

      ptr = (long)RDB[RES_ANA_SLOW_TIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);

      /* Thermal neutron lifetime */

      ptr = (long)RDB[RES_ANA_THERM_TIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufMean(ptr, n);
      AddStat(val, ptr, n);
    }

  /* Get source weight */

  ptr = (long)RDB[RES_INI_SRC_WGT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  div = BufVal(ptr, 0);

  /* Total delayed neutron fraction */

  if (nsf > 0.0)
    {
      ptr = (long)RDB[RES_FWD_ANA_BETA_ZERO];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      beta = BufVal(ptr, 0)/nsf;
    }
  else
    beta = 0.0;

  /* Total neutron thermalization fraction */

  ptr = (long)RDB[RES_ANA_THERM_TIME];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  val = BufWgt(ptr, 0);

  ptr = (long)RDB[RES_ANA_THERM_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddStat(val/div, ptr, 0);

  /* Prompt neutron thermalization fraction */

  ptr = (long)RDB[RES_ANA_THERM_TIME];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  val = BufWgt(ptr, 1);

  ptr = (long)RDB[RES_ANA_THERM_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddStat(val/((1.0 - beta)*div), ptr, 1);

  /* Delayed neutron thermalization fraction */

  if ((beta*div) > 0.0)
    {
      ptr = (long)RDB[RES_ANA_THERM_TIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufWgt(ptr, 2);

      ptr = (long)RDB[RES_ANA_THERM_FRAC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(val/(beta*div), ptr, 2);
    }

  /* Analog estimate of delayed neutron emission time */

  ptr = (long)RDB[RES_ANA_DELAYED_EMTIME];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  val = BufMean(ptr, 0);
  AddStat(val, ptr, 0);

  /* Average number of collision */

  ptr = (long)RDB[RES_ANA_MEAN_NCOL];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  val = BufMean(ptr, 0);
  AddStat(val, ptr, 0);

  val = BufMean(ptr, 1);
  AddStat(val, ptr, 1);

  /***************************************************************************/

  /***** Fission matrix ******************************************************/

  /* Get pointer to matrix data */

  if ((loc0 = (long)RDB[DATA_PTR_FMTX]) > VALID_PTR)
    {
      /* Get size */

      ng = (long)RDB[loc0 + FMTX_SIZE];

      /* Get pointer to source term */

      pts = (long)RDB[loc0 + FMTX_PTR_SRC];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      /* Get pointer to matrix */

      ptr = (long)RDB[loc0 + FMTX_PTR_MTX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over groups and source regions */

      for (i = 0; i < 3; i++)
        {
#ifdef OPEN_MP
#pragma omp parallel private(n, m, div, val)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            for (n = 0; n < ng; n++)
              {
                /* Source term */

                div = BufVal(pts, i, n);

                /* Loop over target regions score */

                for (m = 0; m < ng; m++)
                  {
                    /* Get value */

                    val = BufVal(ptr, i, n, m);

                    /* Check divisor and score */

                    if (div > 0.0)
                      AddStat(val/div, ptr, i, n, m);
                  }
              }
          }
        }
    }

  /***************************************************************************/

  /***** Multi-physics interface *********************************************/

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];
  while (loc0 > VALID_PTR)
    {
      /* Return if inactive cycles and no coupled calculation */

      if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) &&
          (RDB[DATA_RUN_CC] == NO))
        return;

      /* Check flag */

      if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == NO)
        {
          /* Next interface */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Check type */

      if (((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FUEP) ||
          ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FPIP))
        {
          /* Loop over fuel performance code interfaces */

          loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
          while (loc1 > VALID_PTR)
            {
              /* Pointer to universe */

              uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
              CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

              /* Get number of nests from nest count or FINIX definition */

              div = RDB[nst + NEST_COUNT];

              if ((ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FINIX]) > VALID_PTR)
                if (RDB[ptr + FINIX_N_RODS] > 0.0)
                  div = RDB[ptr + FINIX_N_RODS];

              /* Get pointer to statistics */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

              loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
              CheckPointer(FUNCTION_NAME, "loc2", DATA_ARRAY, loc2);

#ifdef OPEN_MP
#pragma omp parallel private(n, m, k, val)
#endif
              {
                /* Loop over regions */

#ifdef OPEN_MP
#pragma omp for
#endif
                for (n = 0; n < (long)RDB[loc2 + FUEP_NZ]; n++)
                  for (m = 0; m < (long)RDB[loc2 + FUEP_NA]; m++)
                    for (k = 0; k < (long)RDB[loc2 + FUEP_NR]; k++)
                      {
                        val = BufVal(ptr, n, m, k)/div;
                        AddStat(norm*val, ptr, n, m, k);
                      }
              }

              /* Get pointer to fast flux */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FLUX];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

              loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];
              CheckPointer(FUNCTION_NAME, "loc2", DATA_ARRAY, loc2);

#ifdef OPEN_MP
#pragma omp parallel private(n, m, k, val)
#endif
              {

                /* Loop over regions */

#ifdef OPEN_MP
#pragma omp for
#endif
                for (n = 0; n < (long)RDB[loc2 + FUEP_NZ]; n++)
                  for (m = 0; m < (long)RDB[loc2 + FUEP_NA]; m++)
                    for (k = 0; k < (long)RDB[loc2 + FUEP_NR]; k++)
                      {
                        val = BufVal(ptr, n, m, k)/div;
                        AddStat(norm*val, ptr, n, m, k);
                      }
              }

              /* Next */

              loc1 = NextItem(loc1);
            }
        }
      else if (((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FET_DENSITY)
               || ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FET_TEMP))
        {
          /* FET Interface */

          loc1 = (long)RDB[loc0 + IFC_FET_OUTPUT_PARAMS_PTR];
          CheckPointer(FUNCTION_NAME, "loc1", DATA_ARRAY, loc1);

          ptr = (long)RDB[loc0 + IFC_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          CollectFET(&RDB[loc1], ptr, (long)RDB[loc0 + IFC_STAT_NREG]);
        }

      else
        {
          /* Get pointer to common statistics */

          ptr = (long)RDB[loc0 + IFC_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          /* Get number of axial bins */

          if ((nz = (long)RDB[loc0 + IFC_NZ]) > 0)
            {
              /* Get number of radial bins */

              nr = (long)RDB[loc0 + IFC_NR];
              CheckValue(FUNCTION_NAME, "nr", "", nr, 1, 1000000);
            }
          else
            nr = 0;

#ifdef OPEN_MP
#pragma omp parallel private(m, n, k, val)
#endif
          {

            /* Loop over regions */

#ifdef OPEN_MP
#pragma omp for
#endif
            for (m = 0; m < (long)RDB[loc0 + IFC_STAT_NREG]; m++)
              {
                /* Check Number of axial bins */

                if (nz > 0)
                  {
                    /* Three-dimensional array, loop over bins */

                    for (n = 0; n < nz; n++)
                      for (k = 0; k < nr; k++)
                        {
                          val = BufVal(ptr, m, n, k);
                          AddStat(norm*val, ptr, m, n, k);
                        }
                  }
                else
                  {
                    /* Single bin */

                    val = BufVal(ptr, m);
                    AddStat(norm*val, ptr, m);
                  }
              }
          }
        }

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /* Check active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /***************************************************************************/

  /***** Normalized total reaction rates for neutrons ************************/

  /* Total cut-off rate (no binned values) */

  ptr = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  cut = BufVal(ptr, 0);
  AddStat(norm*cut, ptr, 0);

  /* Total leak rate (no binned values) */

  ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  leak = BufVal(ptr, 0);
  AddStat(norm*leak, ptr, 0);

  /* Albedo leak rate */

  ptr = (long)RDB[RES_ALB_NEUTRON_LEAKRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  val = BufVal(ptr, 0);
  AddStat(norm*val, ptr, 0);

  /* Add to leak term */

  leak = leak + val;

  /* Total reaction rate */

  ptr = (long)RDB[RES_TOT_NEUTRON_RR];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  val = BufVal(ptr, 0);
  AddStat(norm*val, ptr, 0);

  /* Loop over total, burnable and non-burnable rates */

  for (i = 0; i < 3; i++)
    {
      /* Total fission rate */

      ptr = (long)RDB[RES_TOT_FISSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fiss = BufVal(ptr, i);
      AddStat(norm*fiss, ptr, i);

      /* Fission nubar */

      if (fiss > 0.0)
        {
          ptr = (long)RDB[RES_TOT_NSF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          nsf = BufVal(ptr, i);

          ptr = (long)RDB[RES_TOT_NUBAR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(nsf/fiss, ptr, i);
        }
      else
        nsf = 0.0;

      /* Mean fission energy */

      if (fiss > 0.0)
        {
          ptr = (long)RDB[RES_TOT_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          fE = BufVal(ptr, i);
          AddStat(fE/fiss/MEV, ptr, i);
        }
      else
        fE = 0.0;

      /* Total power */

      ptr = (long)RDB[RES_TOT_NEUTRON_POWER];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(norm*fE, ptr, i);

      /* Total power density */

      if (i == 0)
        fmass = RDB[DATA_INI_FMASS];
      else if (i == 1)
        fmass = RDB[DATA_INI_BURN_FMASS];
      else
        fmass = RDB[DATA_INI_FMASS] - RDB[DATA_INI_BURN_FMASS];

      if (fmass > 0.0)
        {
          ptr = (long)RDB[RES_TOT_POWDENS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(1E-6*norm*fE/fmass, ptr, i);
        }

      /* Total generation rate */

      ptr = (long)RDB[RES_TOT_GENRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(norm*nsf, ptr, i);

      /* Total capture rate */

      ptr = (long)RDB[RES_TOT_CAPTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      capt = BufVal(ptr, i);
      AddStat(norm*capt, ptr, i);

      /* Total absorption rate */

      ptr = (long)RDB[RES_TOT_ABSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(norm*(capt + fiss), ptr, i);

      /* Total loss rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_LOSSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(norm*(leak + capt + fiss), ptr, i);

      /* Total flux */

      ptr = (long)RDB[RES_TOT_NEUTRON_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, i);
      AddStat(norm*val, ptr, i);

      /* Total source rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, i);
      AddStat(norm*val, ptr, i);
    }

  /***************************************************************************/

  /***** Analog reaction rate estimators for neutrons ************************/

  /* Conversion ratio */

  ptr = (long)RDB[RES_ANA_CONV_RATIO];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  val = BufVal(ptr, 1);
  div = BufVal(ptr, 2);

  if (div > 0.0)
    AddStat(val/div, ptr, 0);

  /* Fission fractions */

  ptr = (long)RDB[RES_ANA_FISS_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  div = BufVal(ptr, 0, 0);

  if (div > 0.0)
    for (n = 1; n < 8; n++)
      {
        val = BufVal(ptr, n, 0);
        AddStat(norm*val, ptr, n, 0);
        AddStat(val/div, ptr, n, 1);
      }

  /* Capture fractions */

  ptr = (long)RDB[RES_ANA_CAPT_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  div = BufVal(ptr, 0, 0);

  if (div > 0.0)
    for (n = 1; n < 10; n++)
      {
        val = BufVal(ptr, n, 0);
        AddStat(norm*val, ptr, n, 0);
        AddStat(val/div, ptr, n, 1);
      }

  /***************************************************************************/

  /***** Six-factor formula **************************************************/

  if (1 != 2)
    {
      /* Reset k-eff */

      keff = 1.0;

      /* Fast fission factor */

      ptr = (long)RDB[RES_SIX_FF_ETA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 0)) > 0.0)
        {
          ptr = (long)RDB[RES_SIX_FF_EPSILON];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0)/div;

          /* Score */

          AddStat(val, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val;
        }

      /* Thermal fission factor */

      ptr = (long)RDB[RES_SIX_FF_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 1)) > 0.0)
        {
          ptr = (long)RDB[RES_SIX_FF_ETA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0)/div;

          /* Score */

          AddStat(val, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val;
        }
      else
        keff = 0.0;

      /* Thermal utilization factor */

      ptr = (long)RDB[RES_SIX_FF_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 0)) > 0.0)
        {
          val = BufVal(ptr, 1)/div;

          /* Score */

          ptr = (long)RDB[RES_SIX_FF_F];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val;
        }

      /* Resonance escape probability */

      ptr = (long)RDB[RES_INI_SRC_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, 0);

      ptr = (long)RDB[RES_SIX_FF_LF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div - BufVal(ptr, 0);

      ptr = (long)RDB[RES_SIX_FF_LT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div - BufVal(ptr, 0);

      if (div > 0.0)
        {
          /* Fraction of neutrons absorbed in thermal region */

          ptr = (long)RDB[RES_SIX_FF_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0)/div;

          /* Score */

          ptr = (long)RDB[RES_SIX_FF_P];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val;
        }

      /* Add to six-factor k-inf */

      ptr = (long)RDB[RES_SIX_FF_KINF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(keff, ptr, 0);

      /* Fast non-leakage probability */

      ptr = (long)RDB[RES_INI_SRC_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, 0);

      ptr = (long)RDB[RES_SIX_FF_LF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);

      AddStat(1.0 - val/div, ptr, 0);

      /* Multiply to k-eff */

      keff = keff*(1.0 - val/div);

      /* Thermal non-leakage probability */

      if ((div = div - val) > 0.0)
        {
          ptr = (long)RDB[RES_SIX_FF_LT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);

          AddStat(1.0 - val/div, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*(1.0 - val/div);
        }

      /* Add to six-factor k-eff */

      ptr = (long)RDB[RES_SIX_FF_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(keff, ptr, 0);
    }
  else
    {
      /* Reset k-eff */

      keff = 1.0;

      /* Total number of new fission neutrons */

      ptr = (long)RDB[RES_SIX_FF_ETA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 0)) > 0.0)
        {
          /* Number of new fission neutrons produced in thermal fission */

          ptr = (long)RDB[RES_SIX_FF_EPSILON];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);

          /* Score fast fission factor */

          AddStat(val/div, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val/div;
        }

      /* Number of emitted fission neutrons */

      ptr = (long)RDB[RES_INI_SRC_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, 0);

      /* Number of escaped fast neutrons */

      ptr = (long)RDB[RES_SIX_FF_LF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div - BufVal(ptr, 0);

      if (div > 0.0)
        {
          /* Number of escaped thermal neutrons */

          ptr = (long)RDB[RES_SIX_FF_LT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);

          /* Number of absorbed thermal neutrons */

          ptr = (long)RDB[RES_SIX_FF_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = val + BufVal(ptr, 0);

          /* Score resoance escape probability */

          ptr = (long)RDB[RES_SIX_FF_P];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val/div, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val/div;
        }

      /* Number of escaped thermal neutrons */

      ptr = (long)RDB[RES_SIX_FF_LT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, 0);

      /* Number of absorbed thermal neutrons */

      ptr = (long)RDB[RES_SIX_FF_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div + BufVal(ptr, 0);

      if (div > 0.0)
        {
          /* Number of thermal neutrons absorbed in fuel */

          val = BufVal(ptr, 1);

          /* Score thermal utilization */

          ptr = (long)RDB[RES_SIX_FF_F];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val/div, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val/div;
        }

      /* Number of thermal neutrons absorbed in fuel */

      ptr = (long)RDB[RES_SIX_FF_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if ((div = BufVal(ptr, 1)) > 0.0)
        {
          /* Number of new fission neutrons from thermal fission */

          ptr = (long)RDB[RES_SIX_FF_ETA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);

          /* Score thermal fission factor */

          AddStat(val/div, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*val/div;
        }
      else
        keff = 0.0;

      /* Add to six-factor k-inf */

      ptr = (long)RDB[RES_SIX_FF_KINF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(keff, ptr, 0);

      /* Number of emitted fission neutrons */

      ptr = (long)RDB[RES_INI_SRC_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, 0);

      /* Number of escaped fast neutrons */

      ptr = (long)RDB[RES_SIX_FF_LF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      val = BufVal(ptr, 0);

      /* Score fast non-leakage probability */

      AddStat(1.0 - val/div, ptr, 0);

      /* Multiply to k-eff */

      keff = keff*(1.0 - val/div);

      /* Number of emitted fission neutrons */

      ptr = (long)RDB[RES_INI_SRC_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = BufVal(ptr, 0);

      /* Number of escaped fast neutrons */

      ptr = (long)RDB[RES_SIX_FF_LF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div - BufVal(ptr, 0);

      /* Number of absorbed fast neutrons */

      ptr = (long)RDB[RES_SIX_FF_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      div = div - BufVal(ptr, 2);

      if (1 == 2)
      if (div  > 0.0)
        {
          /* Number of escaped thermal neutrons */

          ptr = (long)RDB[RES_SIX_FF_LT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          val = BufVal(ptr, 0);

          /* Score thermal non-leakage probability */

          AddStat(1.0 - val/div, ptr, 0);

          /* Multiply to k-eff */

          keff = keff*(1.0 - val/div);
        }

      /* Add to six-factor k-eff */

      ptr = (long)RDB[RES_SIX_FF_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(keff, ptr, 0);
    }

  /***************************************************************************/

  /***** Interface current method ********************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_ICM0];
  while (loc0 > VALID_PTR)
    {
      /* Number of energy groups, segment and mu bins */

      ng0 = (long)RDB[DATA_ICM_NG0];
      ng1 = (long)RDB[DATA_ICM_NG1];
      nseg = (long)RDB[DATA_ICM_NSEG];
      nmua = (long)RDB[DATA_ICM_NMU1];
      nmus = (long)(RDB[DATA_ICM_NMU0]*RDB[DATA_ICM_NMU2]);
      np = (long)RDB[loc0 + ICM_NP];

      /* Loop over data */

      for (n0 = 0; n0 < nseg; n0++)
        for (ma0 = 0; ma0 < nmua; ma0++)
          for (ms0 = 0; ms0 < nmus; ms0++)
            for (g0 = 0; g0 < ng0; g0++)
              {
                /* Get inward current */

                ptr = (long)RDB[loc0 + ICM_RES_CURR0];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                div = BufVal(ptr, n0, ma0, ms0, g0);

                /* Check */

                if (div > 0.0)
                  {
                    /* Coupling coefficients */

                    ptr = (long)RDB[loc0 + ICM_RES_CC1];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    for (n1 = 0; n1 < nseg; n1++)
                      for (ma1 = 0; ma1 < nmua; ma1++)
                        for (ms1 = 0; ms1 < nmus; ms1++)
                          for (g1 = 0; g1 < ng0; g1++)
                            {
                              val = BufVal(ptr, n0, ma0, ms0, g0,
                                           n1, ma1, ms1, g1);
                              AddStat(val/div, ptr, n0, ma0, ms0, g0,
                                      n1, ma1, ms1, g1);
                            }

                    ptr = (long)RDB[loc0 + ICM_RES_CC2];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    for (n1 = 0; n1 < nseg; n1++)
                      for (ma1 = 0; ma1 < nmua; ma1++)
                        for (ms1 = 0; ms1 < nmus; ms1++)
                          for (g1 = 0; g1 < ng0; g1++)
                            {
                              val = BufVal(ptr, n0, ma0, ms0, g0,
                                           n1, ma1, ms1, g1);
                              AddStat(val/div, ptr, n0, ma0, ms0, g0,
                                      n1, ma1, ms1, g1);
                            }

                    /* Assembly flux reconstruction factors */

                    if ((ptr = (long)RDB[loc0 + ICM_RES_AFLX1]) > VALID_PTR)
                      for (g1 = 0; g1 < ng1; g1++)
                        {
                          val = BufVal(ptr, n0, ma0, ms0, g0, g1);
                          AddStat(norm*val/div, ptr, n0, ma0, ms0, g0, g1);
                        }

                    if ((ptr = (long)RDB[loc0 + ICM_RES_AFLX2]) > VALID_PTR)
                      for (g1 = 0; g1 < ng1; g1++)
                        {
                          val = BufVal(ptr, n0, ma0, ms0, g0, g1);
                          AddStat(norm*val/div, ptr, n0, ma0, ms0, g0, g1);
                        }

                    /* Assembly fiss and src rate reconstruction factors */

                    ptr = (long)RDB[loc0 + ICM_RES_AFISS1];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    ptr = (long)RDB[loc0 + ICM_RES_ASRC1];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    ptr = (long)RDB[loc0 + ICM_RES_AFISS2];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    ptr = (long)RDB[loc0 + ICM_RES_ASRC2];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    /* Assembly absorption rate reconstruction factors */

                    ptr = (long)RDB[loc0 + ICM_RES_AABS1];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    ptr = (long)RDB[loc0 + ICM_RES_AABS2];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    /* Assembly power reconstruction factors */

                    ptr = (long)RDB[loc0 + ICM_RES_APOW1];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    ptr = (long)RDB[loc0 + ICM_RES_APOW2];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    /* Leak rate reconstruction factors */

                    ptr = (long)RDB[loc0 + ICM_RES_LEAK1];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    ptr = (long)RDB[loc0 + ICM_RES_LEAK2];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    val = BufVal(ptr, n0, ma0, ms0, g0);
                    AddStat(norm*val/div, ptr, n0, ma0, ms0, g0);

                    /* Check number of pins */

                    if (np > 0)
                      {
                        /* Pin flux reconstruction factors */

                        if ((ptr = (long)RDB[loc0 + ICM_RES_PFLX1]) > VALID_PTR)
                          for (g1 = 0; g1 < ng1; g1++)
                            for (i = 0; i < np; i++)
                              {
                                val = BufVal(ptr, n0, ma0, ms0, g0, g1, i);
                                AddStat(norm*val/div, ptr, n0, ma0, ms0, g0,
                                        g1, i);
                              }

                        if ((ptr = (long)RDB[loc0 + ICM_RES_PFLX2]) > VALID_PTR)
                          for (g1 = 0; g1 < ng1; g1++)
                            for (i = 0; i < np; i++)
                              {
                                val = BufVal(ptr, n0, ma0, ms0, g0, g1, i);
                                AddStat(norm*val/div, ptr, n0, ma0, ms0, g0,
                                        g1, i);
                              }

                        /* Pin power reconstruction factors */

                        ptr = (long)RDB[loc0 + ICM_RES_PPOW1];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        for (i = 0; i < np; i++)
                          {
                            val = BufVal(ptr, n0, ma0, ms0, g0, i);
                            AddStat(norm*val/div, ptr, n0, ma0, ms0, g0, i);
                          }

                        ptr = (long)RDB[loc0 + ICM_RES_PPOW2];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        for (i = 0; i < np; i++)
                          {
                            val = BufVal(ptr, n0, ma0, ms0, g0, i);
                            AddStat(norm*val/div, ptr, n0, ma0, ms0, g0, i);
                          }
                      }
                  }
              }

      /* Next set */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** MORA cross sections *************************************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Pointer to MORA data */

      if ((loc0 = (long)RDB[gcu + GCU_PTR_MORA]) > VALID_PTR)
        {
          /* Number of energy groups and cosine bins */

          ng = (long)RDB[loc0 + MORA_N_EG];
          nc = (long)RDB[loc0 + MORA_N_COS];

          /* Get total flux */

          ptr = (long)RDB[loc0 + MORA_PTR_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          flx = BufVal(ptr, ng);

          /* Get total prompt and delayed chi */

          ptr = (long)RDB[loc0 + MORA_PTR_CHIP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          wgt0 = BufVal(ptr, ng);

          ptr = (long)RDB[loc0 + MORA_PTR_CHID];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          wgt1 = BufVal(ptr, ng);

          /* Loop over energy groups */

          for (n = 0; n < ng; n++)
            {
              ptr = (long)RDB[loc0 + MORA_PTR_FLX];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              div = BufVal(ptr, n);

              if (flx > 0.0)
                AddStat(div/flx, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              tot = BufVal(ptr, n );

              if (div > 0.0)
                AddStat(tot/div, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_CAPT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              capt = BufVal(ptr, n );

              if (div > 0.0)
                AddStat(capt/div, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_FISS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              fiss = BufVal(ptr, n );

              if (div > 0.0)
                AddStat(fiss/div, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_KAPPA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n );

              if (fiss > 0.0)
                AddStat(val/fiss, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_PNU];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n );

              if (fiss > 0.0)
                AddStat(val/fiss, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_DNU];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n );

              if (fiss > 0.0)
                AddStat(val/fiss, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_CHIP];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n );

              if (wgt0 > 0.0)
                AddStat(val/wgt0, ptr, n);

              ptr = (long)RDB[loc0 + MORA_PTR_CHID];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              val = BufVal(ptr, n );

              if (wgt1 > 0.0)
                AddStat(val/wgt1, ptr, n);

              scatt = tot - fiss - capt;

              for (m = 0; m < ng; m++)
                for (i = 0; i < nc; i++)
                  if (div > 0.0)
                    {
                      ptr = (long)RDB[loc0 + MORA_PTR_SCATTP];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      val = BufVal(ptr, m, n, i);
                      AddStat(val/div, ptr, m, n, i);

                      ptr = (long)RDB[loc0 + MORA_PTR_SCATTW];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      val = BufVal(ptr, m, n, i);
                      AddStat(val/div, ptr, m, n, i);
                    }
            }
        }

      /* Next universe */

      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Microscopic reaction rates  for neutrons (analog estimators) ********/

  /* Check mode */

  if ((long)RDB[DATA_ANA_RR_NCALC] != ARR_MODE_NONE)
    {
      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)
            {
              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Get pointer */

                  if ((ptr = (long)RDB[rea + REACTION_PTR_ANA_RATE]) >
                      VALID_PTR)
                    {
                      /* Add to statistics */

                      val = BufVal(ptr, 0);
                      AddStat(val*norm, ptr, 0);
                    }

                  /* Next reaction */

                  rea = NextItem(rea);
                }
            }

          /* Next nuclide */

          nuc = NextItem(nuc);
        }
    }

  /***************************************************************************/

  /***** Transmutation cross sections for burnup calculation *****************/

  /* Check burnup mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
    {
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if domain decomposition is in use */

          if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Check material burn-flag */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            {
              /* Check for leakage spectrum correction (tää tarvitaan, */
              /* koska reaktionopeuksia ei skoorata ennen kuin toi     */
              /* rutiini on ajettu kerran. */

              if (((long)RDB[DATA_B1_BURNUP_CORR] == NO) ||
                  ((long)RDB[DATA_OPTI_GC_CALC] == NO) ||
                  ((long)RDB[DATA_MICRO_CALC_BATCH_NUM] > 0))
                {
                  /* Total flux for normalization */

                  ptr = (long)RDB[mat + MATERIAL_PTR_BURN_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  val = BufVal(ptr, 0);
                  AddStat(val*norm, ptr, 0);
                }
            }

          /* Next material */

          mat = NextItem(mat);
        }
    }

  /***************************************************************************/

  /***** Core power distributions ********************************************/

  /* Check if distribution is given */

  if (RDB[DATA_CORE_PDE_DEPTH] > 0.0)
    {
      /* Core level */

      if ((ptr = (long)RDB[DATA_CORE_PDE_PTR_RES0]) > VALID_PTR)
        {

#ifdef OPEN_MP
#pragma omp parallel private(n0, val)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
              {
                val = BufVal(ptr, n0);
                AddStat(val*norm, ptr, n0);
              }
          }
        }

      /* Pin level */

      if ((ptr = (long)RDB[DATA_CORE_PDE_PTR_RES1]) > VALID_PTR)
        {

#ifdef OPEN_MP
#pragma omp parallel private(n0, n1, val)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif

            for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
              for (n1 = 0; n1 < (long)RDB[DATA_CORE_PDE_N1]; n1++)
                {
                  val = BufVal(ptr, n0, n1);
                  AddStat(val*norm, ptr, n0, n1);
                }
          }
        }

      /* Region level */

      if ((ptr = (long)RDB[DATA_CORE_PDE_PTR_RES2]) > VALID_PTR)
        {
          if ((long)RDB[DATA_CORE_PDE_N1] > 0)
            {
#ifdef OPEN_MP
#pragma omp parallel private(n0, n1, n2, val)
#endif
              {

#ifdef OPEN_MP
#pragma omp for
#endif
                for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
                  for (n1 = 0; n1 < (long)RDB[DATA_CORE_PDE_N1]; n1++)
                    for (n2 = 0; n2 < (long)RDB[DATA_CORE_PDE_N2]; n2++)
                      {
                        val = BufVal(ptr, n0, n1, n2);
                        AddStat(val*norm, ptr, n0, n1, n2);
                      }
              }
            }
          else
            {

#ifdef OPEN_MP
#pragma omp parallel private(n0, n2, val)
#endif
              {

#ifdef OPEN_MP
#pragma omp for
#endif
                for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
                  for (n2 = 0; n2 < (long)RDB[DATA_CORE_PDE_N2]; n2++)
                    {
                      val = BufVal(ptr, n0, 0, n2);
                      AddStat(val*norm, ptr, n0, 0, n2);
                    }
              }
            }
        }
    }

  /***************************************************************************/

  /***** Pebble-bed power distributions **************************************/

  /* Loop over pbed geometries */

  pbd = (long)RDB[DATA_PTR_PB0];
  while (pbd > VALID_PTR)
    {
      /* Get pointer to distribution */

      if ((ptr = (long)RDB[pbd + PBED_PTR_POW]) > VALID_PTR)
        {
#ifdef OPEN_MP
#pragma omp parallel private(i, val)
#endif
          {

            /* Loop over pebbles */

#ifdef OPEN_MP
#pragma omp for
#endif
            for (i = 0; i < (long)RDB[pbd + PBED_N_PEBBLES]; i++)
              {
                /* Get power and add to statistics */

                val = BufVal(ptr, i);
                AddStat(val*norm, ptr, i);
              }
          }
        }

      /* Next */

      pbd = NextItem(pbd);
    }

  /***************************************************************************/

  /* Iterable nuclide reaction rates */

  /* Absorption rate */

  ptr = (long)RDB[RES_TOT_ITER_NUC_ABSRATE];

  val = BufVal(ptr,0);

  AddStat(val*norm, ptr, 0);

  /***************************************************************************/
}

/*****************************************************************************/
