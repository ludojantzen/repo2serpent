/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eventstosensitivity.c                          */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2018/09/17 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores    particle events to sensitivities after particle    */
/*              track has been terminated                                    */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EventsToSensitivity:"

/*****************************************************************************/

void EventsToSensitivity(long part, double wgt, long detstat, double detval, long id)
{
  long loc0, loc1, loc2, label, gen, n, keffptr, buffptr, beffptr, beffgptr;
  long leffptr, lambdaptr, lambdagptr, voidptr, dngroups, mingen, maxg;
  long imat, izai, irea, iene, imu, mat0, zai0;
  long nmat, nzai, nrea, nene, nmu;
  long zaiarr, matarr, reaarr, itotnubar, itotchi;
  double val, hitmiss, lifes, nacolls, lambda;

#ifndef OLD_IFP

  double t;
  long ng;

#endif

  dngroups = 0;
  lifes = 0.0;
  nacolls = 0.0;

  /* No need to score sensitivities during inactive cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Get pointer to sensitivity block or return */

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /*************************************************************/
  /* Get prompt lifetime, delayed neutron group and number of  */
  /* collisions in void fraction material from ancestor        */
  /*************************************************************/

  if (detstat < VALID_PTR)
    GetIFPAncestorData(part, &lifes, &dngroups, &lambda, &nacolls);

  /**********************************************************************/

  maxg = (long)RDB[DATA_SENS_LAST_GEN];

  /* Get pointers to sensitivity results */

  /* Multiplication factor */

  keffptr = (long)RDB[RES_ADJ_PERT_KEFF_SENS];
  CheckPointer(FUNCTION_NAME, "(keffptr)", DATA_ARRAY, keffptr);

  /* Total buffer (tracks total weight in all scored generations) */

  buffptr = (long)RDB[RES_SENS_BUFF];
  CheckPointer(FUNCTION_NAME, "(buffptr)", DATA_ARRAY, buffptr);

  /* Delayed neutron fraction */

  beffptr = (long)RDB[RES_ADJ_PERT_BEFF_SENS];
  CheckPointer(FUNCTION_NAME, "(beffptr)", DATA_ARRAY, beffptr);

  /* Delayed neutron decay constant */

  lambdaptr = (long)RDB[RES_ADJ_PERT_LAMBDA_SENS];
  CheckPointer(FUNCTION_NAME, "(lambdaptr)", DATA_ARRAY, lambdaptr);

  /* Get pointers to group-wise delayed neutron statistics */

  if (dngroups > 0)
    {
      switch (dngroups)
        {
        case 1:
          beffgptr = RES_ADJ_PERT_BEFF_G1_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G1_SENS;
          break;
        case 2:
          beffgptr = RES_ADJ_PERT_BEFF_G2_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G2_SENS;
          break;
        case 3:
          beffgptr = RES_ADJ_PERT_BEFF_G3_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G3_SENS;
          break;
        case 4:
          beffgptr = RES_ADJ_PERT_BEFF_G4_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G4_SENS;
          break;
        case 5:
          beffgptr = RES_ADJ_PERT_BEFF_G5_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G5_SENS;
          break;
        case 6:
          beffgptr = RES_ADJ_PERT_BEFF_G6_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G6_SENS;
          break;
        case 7:
          beffgptr = RES_ADJ_PERT_BEFF_G7_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G7_SENS;
          break;
        case 8:
          beffgptr = RES_ADJ_PERT_BEFF_G8_SENS;
          lambdagptr = RES_ADJ_PERT_LAMBDA_G8_SENS;
          break;
        default:
          beffgptr = NULLPTR;
          lambdagptr = NULLPTR;
          Die(FUNCTION_NAME, "WTF!");
        }

      /* Get pointer to group-wise statistics */

      beffgptr = (long)RDB[beffgptr];
      lambdagptr = (long)RDB[lambdagptr];
    }
  else
    {
      beffgptr = NULLPTR;
      lambdagptr = NULLPTR;
    }

  /* Prompt neutron generation time */

  leffptr = (long)RDB[RES_ADJ_PERT_LEFF_SENS];
  CheckPointer(FUNCTION_NAME, "(leffptr)", DATA_ARRAY, leffptr);

  /* Void reactivity coefficient */

  voidptr = (long)RDB[RES_ADJ_PERT_VOID_SENS];
  CheckPointer(FUNCTION_NAME, "(voidptr)", DATA_ARRAY, voidptr);

  /***********************************************************************/
  /*************** Process events to sensitivities ***********************/
  /***********************************************************************/

  if ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_EVENT)
    {
      /* Get pointer to the most recent event */

      loc1 = (long)RDB[part + PARTICLE_PTR_EVENTS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check if last event experienced by neutron was a fission */
      /* to figure out generation numbering */

      if ((long)RDB[loc1 + EVENT_TYPE] == EVENT_TYPE_FISS)
        gen = -1;
      else
        gen = 0;

      /* Loop over event history */

      while (loc1 > VALID_PTR)
        {
          /* Check if reaction is a fission */

          if ((long)RDB[loc1 + EVENT_TYPE] == EVENT_TYPE_FISS)
            {
              /* Increment generation due to fission */

              gen++;

              /* Check generation */

              if (gen == maxg)
                break;
            }

          /************************/
          /* Score stuff **********/
          /************************/

          /* Get label */

          label = (long)RDB[loc1 + EVENT_LABEL];

          if (label == 0)
            {
              /* Event not included in sensitivity calculation */
              /* Get next event and cycle loop */

              loc1 = NextItem(loc1);
              continue;
            }

          /* Make the label positive so that it can be used as an index */

          if (label < 0)
            {
              hitmiss = -1.0;
              label = -label;
            }
          else
            hitmiss = 1.0;

          /* Get event value and mulitply with hit/miss */

          val = RDB[loc1 + EVENT_VAL]*hitmiss;

          /************************************/
          /* Get minimum generation to score  */
          /************************************/

          /* Get minimum ancestor estimator to score when summing: gen N doesn't */
          /* contribute to sums with a lower number of latent generations... */

          if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_SCORE_FLAG_HIS)
            mingen = gen;
          else
            mingen = (long)RDB[DATA_SENS_LAST_GEN] - 1;

          /* Detector sensitivities don't score k-eff, leff etc. */

          if (detstat > VALID_PTR)
            for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
              ScoreSensLabel(detstat, n, label, val*detval, wgt, id);
          else
            {
              /* K-eff sensitivities */

              ScoreSensLabel(keffptr, gen, label, val, wgt, id);

              /* Score total weight in all generations */

              for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                ScoreSensLabel(buffptr, n, label, val, wgt, id);

              /* Check whether ancestor was a delayed neutron and score */
              /* beta-eff and ell-eff sensitivities */

              if (dngroups > 0)
                {
                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_BEFF)
                    {
                      /* Score total beff */

                      for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                        ScoreSensLabel(beffptr, n, label, val, wgt, id);

                      /* Score group-wise beff */

                      for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                        ScoreSensLabel(beffgptr, n, label, val, wgt, id);

                    }

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LAMBDA)
                    {
                      /* Score total lambda */

                      for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                        ScoreSensLabel(lambdaptr, n, label, lambda*val, wgt, id);

                      /* Score group-wise lambda */

                      for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                        ScoreSensLabel(lambdagptr, n, label, lambda*val, wgt, id);

                    }
                }
              else if (lifes > 0.0)
                {
                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LEFF)
                    for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                      ScoreSensLabel(leffptr, n, label, lifes*val, wgt, id);
                }

              /* Number of collisions in Sodium (Na) used for void fraction sensitivity? */

              if (nacolls != 0.0)
                if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_VOID)
                  for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                    ScoreSensLabel(voidptr, n, label, nacolls*val, wgt, id);
            }

          /* Pointer to next event */

          loc1 = NextItem(loc1);
        }

    }
  else if ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_DIRECT)
    {
      /* Get pointer to the 0-gen block */

      loc1 = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get sizes of sub bins */

      nmat = (long)RDB[loc0 + SENS_N_MAT];
      nzai = (long)RDB[loc0 + SENS_N_ZAI];
      nrea = (long)RDB[loc0 + SENS_N_PERT] + 1;
      nene = (long)RDB[loc0 + SENS_N_ENE] + 1;
      nmu  = (long)RDB[loc0 + SENS_N_MU];

      imu = 0;

      /* Get pointer to index array */

      matarr = (long)RDB[loc0 + SENS_PTR_MAT_INDICES];
      CheckPointer(FUNCTION_NAME, "(matarr)", DATA_ARRAY, matarr);

      /* Get index to first material that should be dealt with */

      mat0 = (long)RDB[matarr + SENS_NON_MAT_IDX];

      /* Get pointer to index array */

      zaiarr = (long)RDB[loc0 + SENS_PTR_ZAI_INDICES];
      CheckPointer(FUNCTION_NAME, "(zaiarr)", DATA_ARRAY, zaiarr);

      /* Get index to first zai that should be dealt with */

      zai0 = (long)RDB[zaiarr + SENS_NON_ZAI_IDX];

      /* Get pointer to index array */

      reaarr = (long)RDB[loc0 + SENS_PTR_PERT_INDICES];
      CheckPointer(FUNCTION_NAME, "(reaarr)", DATA_ARRAY, reaarr);

      /* Get indices for total nubar and total chi as they can be */
      /* skipped here */

      if ((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
        itotnubar = (long)RDB[reaarr + NUBAR_TOT_IDX];
      else
        itotnubar = -1;

      if ((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
        itotchi = (long)RDB[reaarr + CHI_TOT_IDX];
      else
        itotchi = -1;

      /* Get index to first zai that should be dealt with */

      zai0 = (long)RDB[zaiarr + SENS_NON_ZAI_IDX];
      /* Reset generation */

      gen = 0;

      /* Loop over all generations */

      while (loc1 > VALID_PTR)
        {
          /************************************/
          /* Get minimum generation to score  */
          /************************************/

          /* Get minimum ancestor estimator to score when summing: gen N doesn't */
          /* contribute to sums with a lower number of latent generations... */

          if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_SCORE_FLAG_HIS)
            mingen = gen;
          else
            mingen = (long)RDB[DATA_SENS_LAST_GEN] - 1;

          /* Get pointer to data */

          loc2 = (long)RDB[loc1 + SENS_EBLOCK_PTR_DATA];

          /*********************************************************************/
          /* Loop over all labels, this loop takes most of the simulation time */
          /* so reducing the number of loop cycles saves a lot of time         */
          /*********************************************************************/

          /* First loop is over materials  */

          for (imat = mat0; imat < nmat; imat++)
            {
              /* Second loop over ZAIs  */

              for (izai = zai0; izai < nzai; izai++)
                {
                  /* Third loop over reactions */

                  for (irea = 1; irea < nrea; irea++)
                    {
                      /* Skip total nubar and total chi */

                      if ((irea == itotnubar) || (irea == itotchi))
                        continue;

                      /* Final loop over energy bins */

                      for (iene = 0; iene < nene; iene++)
                        {
                          /* Calculate bin idx */

                          label = 1 + imat*nzai*nrea*nene*nmu +
                            izai*nrea*nene*nmu +
                            irea*nene*nmu +
                            iene*nmu + imu;

                          /* Get label value */

                          val = RDB[loc2 + label];

                          /* Skip zero-valued labels */

                          if (val == 0)
                            continue;

                          /* Detector sensitivities */

                          if (detstat > VALID_PTR)
                            {

                              for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                ScoreSensLabel(detstat, n, label, val*detval, wgt, id);

                              /* Do not score others when called from coldet.c */

                              continue;
                            }

                          /* K-eff sensitivities */

                          ScoreSensLabel(keffptr, gen, label, val, wgt, id);

                          /* Score total weight in all generations */

                          for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                            ScoreSensLabel(buffptr, n, label, val, wgt, id);

                          /* Check whether ancestor was a delayed neutron and score */
                          /* beta-eff and ell-eff sensitivities */

                          if (dngroups > 0)
                            {
                              if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_BEFF)
                                {
                                  /* Score sensitivity for total beff */

                                  for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                    ScoreSensLabel(beffptr, n, label, val, wgt, id);

                                  /* Score sensitivity for specific group */

                                  for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                    ScoreSensLabel(beffgptr, n, label, val, wgt, id);

                                }

                              if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LAMBDA)
                                {
                                  /* Score total lambda */

                                  for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                    ScoreSensLabel(lambdaptr, n, label, lambda*val, wgt, id);

                                  /* Score group-wise lambda */

                                  for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                    ScoreSensLabel(lambdagptr, n, label, lambda*val, wgt, id);

                                }
                            }
                          else if (lifes > 0.0)
                            {
                              if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LEFF)
                                for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                  ScoreSensLabel(leffptr, n, label, val*lifes, wgt, id);
                            }

                          /* Number of collisions in Sodium (Na) used for void fraction sensitivity? */

                          if (nacolls != 0.0)
                            if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_VOID)
                              for (n = mingen; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
                                ScoreSensLabel(voidptr, n, label, nacolls*val, wgt, id);
                        }
                    }
                }
            }

          /* Pointer to next event block */

          loc1 = NextItem(loc1);

          /* Increment generation */

          gen++;

          /* Check generation */

          if (gen == maxg)
            break;
        }
    }

  /*****************************************/
  /* score total weights in zero locations */
  /*****************************************/

  if (detstat < VALID_PTR)
    {
      /* Total weight for keff buffer */

      for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
        AddBuf(0.0, wgt, keffptr, id, -1, n, 0);

      /* Total weight in all generations */

      for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
        AddBuf(0.0, wgt, buffptr, id, -1, n, 0);

      /* beta-eff and ell-eff */

      if (dngroups > 0)
        {
          /* Score total weight for total beff buffer */
          /* = IFP estimate for beta-eff */

          for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
            AddBuf(0.0, wgt, beffptr, id, -1, n, 0);

          /* Score total weight for group-wise beff buffer */
          /* = IFP estimate for group-wise beta-eff */

          for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
            AddBuf(0.0, wgt, beffgptr, id, -1, n, 0);

          /* Score total weight for total lambda buffer */
          /* = IFP estimate for lambda */

          for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
            AddBuf(0.0, wgt*lambda, lambdaptr, id, -1, n, 0);

          /* Score total weight for group-wise lambda buffer */
          /* = IFP estimate for group-wise lambda */

          for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
            AddBuf(0.0, wgt*lambda, lambdagptr, id, -1, n, 0);
        }
      else if (lifes > 0.0)
        for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
          AddBuf(0.0, wgt*lifes, leffptr, id, -1, n, 0);

      /* Void coefficient */
      /* Note: nacolls will score negative value */
      /* = IFP estimate for void reactivity coefficient */

      if (nacolls != 0.0)
        for (n = 0; n < (long)RDB[DATA_SENS_LAST_GEN]; n++)
          AddBuf(nacolls, wgt, voidptr, id, -1, n, 0);
    }

}
