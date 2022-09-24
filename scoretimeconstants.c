/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoretimeconstants.c                           */
/*                                                                           */
/* Created:       2012/09/05 (JLe)                                           */
/* Last modified: 2017/06/07 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Scores various analog time constants                         */
/*                                                                           */
/* Comments: - Scored at the end of particle lifetime after absorption,      */
/*             escape or cut-off.                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreTimeConstants:"

/*****************************************************************************/

void ScoreTimeConstants(double t, double wgt, long part, long trk, long id)
{
  long ptr, loc0, n, ng;
  double t0, tt, td, ti, wb, wl, a, c, d, lambda;

  /* Get emission time */

  t0 = RDB[part + PARTICLE_T0];

  /* Check times */

  CheckValue(FUNCTION_NAME, "t", "", t, t0, INFTY);

  /***************************************************************************/

  /***** Analog time constants ***********************************************/

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] == PARTICLE_TYPE_GAMMA)
    {
      /* Score photon lifetime */

      ptr = (long)RDB[RES_ANA_PHOTON_LIFETIME];
      AddBuf1D(t - t0, wgt, ptr, id, 0);
    }
  else
    {
      /* Score delayed neutron emission time */

      if ((td = RDB[part + PARTICLE_TD]) > 0.0)
        {
          ptr = (long)RDB[RES_ANA_DELAYED_EMTIME];
          AddBuf1D(td, wgt, ptr, id, 0);
        }

      /* Score slowing-down time and thermal neutron lifetime */

      if ((tt = RDB[part + PARTICLE_TT]) > 0.0)
        {
          /* Check */

          CheckValue(FUNCTION_NAME, "tt", "", tt, t0, t);

          /* Score total slowing-down time */

          ptr = (long)RDB[RES_ANA_SLOW_TIME];
          AddBuf1D(tt - t0, wgt, ptr, id, 0);

          /* Score prompt or delayed */

          if (td == 0.0)
            AddBuf1D(tt - t0, wgt, ptr, id, 1);
          else
            AddBuf1D(tt - t0, wgt, ptr, id, 2);

          /* Score total thermal neutron lifetime */

          ptr = (long)RDB[RES_ANA_THERM_TIME];
          AddBuf1D(t - tt, wgt, ptr, id, 0);

          /* Score prompt or delayed */

          if (td == 0.0)
            AddBuf1D(t - tt, wgt, ptr, id, 1);
          else
            AddBuf1D(t - tt, wgt, ptr, id, 2);
        }
    }

  /***************************************************************************/

  /***** IFP time constants **************************************************/

#ifdef OLD_IFP

  /* Analog IFP estimator of beta-eff */

  if ((loc0 = (long)RDB[part + PARTICLE_PTR_FISS_PROG]) > VALID_PTR)
    {
      /* Loop over progenies */

      for (n = 0; n < (long)RDB[DATA_IFP_CHAIN_LENGTH]; n++)
        {
          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Get pointer to stats */

          ptr = (long)RDB[RES_ADJ_IFP_ANA_BETA_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Score total */

          AddBuf(1.0, wgt, ptr, id, -1, 9, n);

          /* Get delayed neutron group */

          if ((ng = (long)RDB[loc0 + FISS_PROG_DN_GROUP]) > 0)
            {
              /* Score beta-eff */

              AddBuf(1.0, wgt, ptr, id, -1, 0, n);
              AddBuf(1.0, wgt, ptr, id, -1, ng, n);

              /* Get decay constant */

              lambda = RDB[loc0 + FISS_PROG_LAMBDA];

              /* Score lambda */

              ptr = (long)RDB[RES_ADJ_IFP_ANA_LAMBDA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf(lambda, wgt, ptr, id, -1, 0, n);
              AddBuf(lambda, wgt, ptr, id, -1, ng, n);
            }

          /* Pointer to next */

          loc0 = NextItem(loc0);
        }
    }

  /* Check track type */

  if ((trk == TRACK_END_FISS) &&
      ((loc0 = (long)RDB[part + PARTICLE_PTR_FISS_PROG]) > VALID_PTR))
    {
      /* Get pointer to tally */

      ptr = (long)RDB[RES_ADJ_IFP_LIFETIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over progenies */

      for (n = 0; n < (long)RDB[DATA_IFP_CHAIN_LENGTH]; n++)
        {
          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Get time */

          if ((ti = RDB[loc0 + FISS_PROG_LIFETIME]) > 0.0)
            {
              /* Score total lifetime */

              AddBuf(ti, wgt, ptr, id, -1, 0, n);

              /* Score prompt or delayed */

              if ((long)RDB[loc0 + FISS_PROG_DN_GROUP] == 0)
                AddBuf(ti, wgt, ptr, id, -1, 1, n);
              else
                AddBuf(ti, wgt, ptr, id, -1, 2, n);
            }

          /* Pointer to next */

          loc0 = NextItem(loc0);
        }

      /* Beta-eff and lifetime by perturbation */

      a = RDB[DATA_PERT_VAR_A];
      c = RDB[DATA_PERT_VAR_C];
      d = 0.0;

      if ((long)RDB[part + PARTICLE_DN_GROUP] > 0)
        wb = wgt*(1.0 + a);
      else
        wb = wgt;

      wl = wgt*(1.0 + c*t);

      /* Loop over progenies */

      loc0 = (long)RDB[part + PARTICLE_PTR_FISS_PROG];
      while (loc0 > VALID_PTR)
        {
          /* Perturbed calculation of delayed neutron fraction */

          if ((long)RDB[loc0 + FISS_PROG_DN_GROUP] > 0)
            wb = wb*(1.0 + a);

          /* Perturbed calculation of lifetime */

          wl = wl*(1.0 + c*RDB[loc0 + FISS_PROG_LIFETIME]);

          /* Check count */

          if ((d = d + 1.0) == RDB[DATA_PERT_N_BATCH])
            break;

          /* Next progeny */

          loc0 = NextItem(loc0);
        }

      /* Calculate exponential */

      d = 1.0/(d + 1.0);

      /* Beta and lambda */

      wb = (pow(wb/wgt, d) - 1.0)/a;
      wl = (pow(wl/wgt, d) - 1.0)/c;

      /* Score */

      ptr = (long)RDB[RES_ADJ_PERT_BETA_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(wb, wgt, ptr, id, 0);

      ptr = (long)RDB[RES_ADJ_PERT_LIFETIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(wl, wgt, ptr, id, 0);
    }

#else

  /* Analog IFP estimator of beta-eff */

  n = 0;

  /* Loop over events */

  loc0 = (long)RDB[part + PARTICLE_PTR_EVENTS];
  while (loc0 > VALID_PTR)
    {
      /* Check if event is fission */

      if ((long)RDB[loc0 + EVENT_TYPE] == EVENT_TYPE_FISS)
        {
          /* Get pointer to stats */

          ptr = (long)RDB[RES_ADJ_IFP_ANA_BETA_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Score total */

          AddBuf(1.0, wgt, ptr, id, -1, 9, n);

          /* Get delayed neutron group */

          if ((ng = (long)RDB[loc0 + EVENT_DN_GROUP]) > 0)
            {
              /* Score beta-eff */

              AddBuf(1.0, wgt, ptr, id, -1, 0, n);
              AddBuf(1.0, wgt, ptr, id, -1, ng, n);

              /* Get decay constant */

              lambda = RDB[loc0 + EVENT_LAMBDA];

              /* Score lambda */

              ptr = (long)RDB[RES_ADJ_IFP_ANA_LAMBDA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf(lambda, wgt, ptr, id, -1, 0, n);
              AddBuf(lambda, wgt, ptr, id, -1, ng, n);
            }

          /* Update index */

          n++;
        }

      /* Check index */

      if (n == (long)RDB[DATA_IFP_CHAIN_LENGTH])
        break;

      /* Pointer to next event */

      loc0 = NextItem(loc0);
    }

  /* Check track type */

  if (trk == TRACK_END_FISS)
    {
      /* Get pointer to tally */

      ptr = (long)RDB[RES_ADJ_IFP_LIFETIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over events */

      n = 0;

      /* Get pointer to last event (should be fission) */

      loc0 = (long)RDB[part + PARTICLE_PTR_EVENTS];

      while (loc0 > VALID_PTR)
        {
          /* Check if event is fission */

          if ((long)RDB[loc0 + EVENT_TYPE] == EVENT_TYPE_FISS)
            {
              /* Get time */

              if ((ti = RDB[loc0 + EVENT_LIFETIME]) > 0.0)
                {
                  /* Score total lifetime */

                  AddBuf(ti, wgt, ptr, id, -1, 0, n);

                  /* Score prompt or delayed */

                  if ((long)RDB[loc0 + EVENT_DN_GROUP] == 0)
                    AddBuf(ti, wgt, ptr, id, -1, 1, n);
                  else
                    AddBuf(ti, wgt, ptr, id, -1, 2, n);
                }

              /* Update index */

              n++;
            }

          /* Check index */

          if (n == (long)RDB[DATA_IFP_CHAIN_LENGTH])
            break;

          /* Pointer to next event */

          loc0 = NextItem(loc0);
        }

      /* Beta-eff and lifetime by perturbation */

      a = RDB[DATA_PERT_VAR_A];
      c = RDB[DATA_PERT_VAR_C];
      d = 0.0;

      if ((long)RDB[part + PARTICLE_DN_GROUP] > 0)
        wb = wgt*(1.0 + a);
      else
        wb = wgt;

      wl = wgt*(1.0 + c*t);

      /* Loop over events */

      n = 0;

      loc0 = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (loc0 > VALID_PTR)
        {
          /* Check if event is fission */

          if ((long)RDB[loc0 + EVENT_TYPE] == EVENT_TYPE_FISS)
            {
              /* Perturbed calculation of delayed neutron fraction */

              if ((long)RDB[loc0 + EVENT_DN_GROUP] > 0)
                wb = wb*(1.0 + a);

              /* Perturbed calculation of lifetime */

              wl = wl*(1.0 + c*RDB[loc0 + EVENT_LIFETIME]);

              /* Check count */

              if ((d = d + 1.0) == RDB[DATA_PERT_N_BATCH])
                break;

              /* Update index */

              n++;
            }

          /* Check index */

          if (n == (long)RDB[DATA_IFP_CHAIN_LENGTH])
            break;

          /* Pointer to next event */

          loc0 = NextItem(loc0);
        }

      /* Calculate exponential */

      d = 1.0/(d + 1.0);

      /* Beta and lambda */

      wb = (pow(wb/wgt, d) - 1.0)/a;
      wl = (pow(wl/wgt, d) - 1.0)/c;

      /* Score */

      ptr = (long)RDB[RES_ADJ_PERT_BETA_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(wb, wgt, ptr, id, 0);

      ptr = (long)RDB[RES_ADJ_PERT_LIFETIME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(wl, wgt, ptr, id, 0);
    }

#endif

  /***************************************************************************/
}

/*****************************************************************************/
