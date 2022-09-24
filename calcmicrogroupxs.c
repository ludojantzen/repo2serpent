/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calcmicrogroupxs.c                             */
/*                                                                           */
/* Created:       2013/03/04 (JLe)                                           */
/* Last modified: 2019/08/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates micro-group cross sections needed for B1 and P1   */
/*                                                                           */
/* Comments: - Calculation of microscopic poison cross sections was changed  */
/*             5.3.2016 / 2.1.26 (JLe)                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalcMicroGroupXS:"

/*****************************************************************************/

void CalcMicroGroupXS()
{
  long ptr, nmg, gcu, n, m;
  double *flx, *tot, *fiss, *abs, *nsf, *fissE, *invv, *chit, *chip, *chid;
  double *s0, *sp0, *s1, *sp1, *s2, *sp2, *s3, *sp3, *s4, *sp4, *s5, *sp5;
  double *s6, *sp6, *s7, *sp7, *fflx;
  double *Xep, *Xemp, *Ip, *Smp, *Pm7p, *Pm8p, *Pm8mp, *Pm9p, *Ia, *Xea, *Xema;
  double *trc, *tflx;
  double *Pm7a, *Pm8a, *Pm8ma, *Pm9a, *Sma, *Xeam, *Xemam, *Smam, sum;

  /* Check active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check if group constants are generated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check if in corrector calculation */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    {
      /* This subroutine performs some checks in debug mode */

      ClearMicroGroupXS();

      /* Exit */

      return;
    }

  /* Add to normalization coefficient */

  WDB[DATA_MICRO_CALC_NORM] = RDB[DATA_MICRO_CALC_NORM]
    + NormCoef(PARTICLE_TYPE_NEUTRON);

  /* Check if batch is full */

  if (RDB[DATA_MICRO_CALC_BATCH_COUNT] < RDB[DATA_MICRO_CALC_BATCH_SIZE])
    return;
  else if (RDB[DATA_MICRO_CALC_BATCH_COUNT] > RDB[DATA_MICRO_CALC_BATCH_SIZE])
    Die(FUNCTION_NAME, "Batch overflow");

  /* Update batch number */

  WDB[DATA_MICRO_CALC_BATCH_NUM] = RDB[DATA_MICRO_CALC_BATCH_NUM] + 1.0;

  /* Get pointer to micro-group structure */

  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of micro- and macro-groups */

  nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

  /* Reduce results buffer */

  WDB[DATA_RES2_REDUCED] = (double)NO;
  ReducePrivateRes();

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /***********************************************************************/

      /***** Get micro-group reaction rates **********************************/

      /* Flux */

      ptr = (long)RDB[gcu + GCU_MICRO_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      flx = &RES2[ptr];

      /* Flux in fissile zones */

      ptr = (long)RDB[gcu + GCU_MICRO_FISS_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fflx = &RES2[ptr];

      /* Reaction rates */

      ptr = (long)RDB[gcu + GCU_MICRO_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      tot = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_FISS];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fiss = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      abs = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      nsf = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_FISSE];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fissE = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_INV_V];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      invv = &RES2[ptr];

      /* Fission spectra */

      ptr = (long)RDB[gcu + GCU_MICRO_CHIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      chit = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_CHIP];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      chip = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_CHID];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      chid = &RES2[ptr];

      /* Scattering matrixes */

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT0];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s0 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP0];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp0 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT1];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s1 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP1];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp1 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s2 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp2 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT3];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s3 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP3];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp3 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT4];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s4 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP4];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp4 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT5];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s5 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP5];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp5 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT6];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s6 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP6];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp6 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT7];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s7 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP7];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp7 = &RES2[ptr];

      /* Transport correction */

      if ((long)RDB[DATA_PTR_TRC0] > VALID_PTR)
        {
          ptr = (long)RDB[gcu + GCU_MICRO_TRC];
          CheckPointer(FUNCTION_NAME, "(ptr25)", RES2_ARRAY, ptr);
          trc = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_TRC_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr25)", RES2_ARRAY, ptr);
          tflx = &RES2[ptr];
        }
      else
        {
          trc = NULL;
          tflx = NULL;
        }

      /* Check poison calculation */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
        {
          ptr = (long)RDB[gcu + GCU_MICRO_I135_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Ip = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_XE135_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Xep = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_XE135M_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Xemp = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM147_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm7p = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM148_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm8p = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM148M_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm8mp = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM149_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm9p = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_SM149_YIELD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Smp = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_I135_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Ia = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_XE135_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Xea = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_XE135M_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Xema = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM147_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm7a = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM148_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm8a = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM148M_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm8ma = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_PM149_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Pm9a = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_SM149_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Sma = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_XE135_MACRO_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Xeam = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_XE135M_MACRO_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Xemam = &RES2[ptr];

          ptr = (long)RDB[gcu + GCU_MICRO_SM149_MACRO_ABS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          Smam = &RES2[ptr];
        }
      else
        {
          /* Avoid compiler warning */

          Ip = NULL;
          Xep = NULL;
          Xemp = NULL;
          Pm7p = NULL;
          Pm8p = NULL;
          Pm8mp = NULL;
          Pm9p = NULL;
          Smp = NULL;
          Ia = NULL;
          Xea = NULL;
          Xema = NULL;
          Pm7a = NULL;
          Pm8a = NULL;
          Pm8ma = NULL;
          Pm9a = NULL;
          Sma = NULL;
          Xeam = NULL;
          Xemam = NULL;
          Smam = NULL;
        }

      /***********************************************************************/

      /***** Calculate cross sections ****************************************/

      /* Cross sections */

      for (n = 0; n < nmg; n++)
        {
          /* Check zero flux */

          if (flx[n] > 0.0)
            {
              /* Calculate cross sections */

              tot[n] = tot[n]/flx[n];
              fiss[n] = fiss[n]/flx[n];
              abs[n] = abs[n]/flx[n];
              nsf[n] = nsf[n]/flx[n];
              fissE[n] = fissE[n]/flx[n];
              invv[n] = invv[n]/flx[n];

              if (trc != NULL)
                trc[n] = trc[n]/flx[n];
            }
          else
            {
              tot[n] = 0.0;
              fiss[n] = 0.0;
              abs[n] = 0.0;
              nsf[n] = 0.0;
              fissE[n] = 0.0;
              invv[n] = 0.0;

              if (trc != NULL)
                trc[n] = 0.0;
            }
        }

      /* Fission spectra */

      sum = 0.0;
      for (n = 0; n < nmg; n++)
        sum = sum + chit[n];

      for (n = 0; n < nmg; n++)
        if (sum > 0.0)
          chit[n] = chit[n]/sum;

      sum = 0.0;
      for (n = 0; n < nmg; n++)
        sum = sum + chip[n];

      for (n = 0; n < nmg; n++)
        if (sum > 0.0)
          chip[n] = chip[n]/sum;

      sum = 0.0;
      for (n = 0; n < nmg; n++)
        sum = sum + chid[n];

      for (n = 0; n < nmg; n++)
        if (sum > 0.0)
          chid[n] = chid[n]/sum;

      /* Scattering matrixes */

      for (n = 0; n < nmg; n++)
        for (m = 0; m < nmg; m++)
          {
            if (flx[m] > 0.0)
              {
                s0[m*nmg + n] = s0[m*nmg + n]/flx[m];
                sp0[m*nmg + n] = sp0[m*nmg + n]/flx[m];
                s1[m*nmg + n] = s1[m*nmg + n]/flx[m];
                sp1[m*nmg + n] = sp1[m*nmg + n]/flx[m];
                s2[m*nmg + n] = s2[m*nmg + n]/flx[m];
                sp2[m*nmg + n] = sp2[m*nmg + n]/flx[m];
                s3[m*nmg + n] = s3[m*nmg + n]/flx[m];
                sp3[m*nmg + n] = sp3[m*nmg + n]/flx[m];
                s4[m*nmg + n] = s4[m*nmg + n]/flx[m];
                sp4[m*nmg + n] = sp4[m*nmg + n]/flx[m];
                s5[m*nmg + n] = s5[m*nmg + n]/flx[m];
                sp5[m*nmg + n] = sp5[m*nmg + n]/flx[m];
                s6[m*nmg + n] = s6[m*nmg + n]/flx[m];
                sp6[m*nmg + n] = sp6[m*nmg + n]/flx[m];
                s7[m*nmg + n] = s7[m*nmg + n]/flx[m];
                sp7[m*nmg + n] = sp7[m*nmg + n]/flx[m];
              }
            else
              {
                s0[m*nmg + n] = 0.0;
                sp0[m*nmg + n] = 0.0;
                s1[m*nmg + n] = 0.0;
                sp1[m*nmg + n] = 0.0;
                s2[m*nmg + n] = 0.0;
                sp2[m*nmg + n] = 0.0;
                s3[m*nmg + n] = 0.0;
                sp3[m*nmg + n] = 0.0;
                s4[m*nmg + n] = 0.0;
                sp4[m*nmg + n] = 0.0;
                s5[m*nmg + n] = 0.0;
                sp5[m*nmg + n] = 0.0;
                s6[m*nmg + n] = 0.0;
                sp6[m*nmg + n] = 0.0;
                s7[m*nmg + n] = 0.0;
                sp7[m*nmg + n] = 0.0;
              }
          }

      /* Poison cross sections */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
        for (n = 0; n < nmg; n++)
          {
            /* Check zero fission rate */

            if (flx[n]*fiss[n] > 0.0)
              {
                Ip[n] = Ip[n]/(flx[n]*fiss[n]);
                Xep[n] = Xep[n]/(flx[n]*fiss[n]);
                Xemp[n] = Xemp[n]/(flx[n]*fiss[n]);
                Pm7p[n] = Pm7p[n]/(flx[n]*fiss[n]);
                Pm8p[n] = Pm8p[n]/(flx[n]*fiss[n]);
                Pm8mp[n] = Pm8mp[n]/(flx[n]*fiss[n]);
                Pm9p[n] = Pm9p[n]/(flx[n]*fiss[n]);
                Smp[n] = Smp[n]/(flx[n]*fiss[n]);
              }
            else
              {
                Ip[n] = 0.0;
                Xep[n] = 0.0;
                Xemp[n] = 0.0;
                Pm7p[n] = 0.0;
                Pm8p[n] = 0.0;
                Pm8mp[n] = 0.0;
                Pm9p[n] = 0.0;
                Smp[n] = 0.0;
              }

            /* Check zero fissile flux */
            /*
            if (fflx[n] > 0.0)
              {
                Ia[n] = Ia[n]/fflx[n];
                Xea[n] = Xea[n]/fflx[n];
                Pma[n] = Pma[n]/fflx[n];
                Sma[n] = Sma[n]/fflx[n];
              }
            else
              {
                Ia[n] = 0.0;
                Xea[n] = 0.0;
                Pma[n] = 0.0;
                Sma[n] = 0.0;
              }
            */
            /* Check zero flux */

            if (flx[n] > 0.0)
              {
                Ia[n] = Ia[n]/flx[n];
                Xea[n] = Xea[n]/flx[n];
                Xema[n] = Xema[n]/flx[n];
                Pm7a[n] = Pm7a[n]/flx[n];
                Pm8a[n] = Pm8a[n]/flx[n];
                Pm8ma[n] = Pm8ma[n]/flx[n];
                Pm9a[n] = Pm9a[n]/flx[n];
                Sma[n] = Sma[n]/flx[n];

                Xeam[n] = Xeam[n]/flx[n];
                Xemam[n] = Xemam[n]/flx[n];
                Smam[n] = Smam[n]/flx[n];
              }
            else
              {
                Ia[n] = 0.0;
                Xea[n] = 0.0;
                Xema[n] = 0.0;
                Pm7a[n] = 0.0;
                Pm8a[n] = 0.0;
                Pm8ma[n] = 0.0;
                Pm9a[n] = 0.0;
                Sma[n] = 0.0;

                Xeam[n] = 0.0;
                Xemam[n] = 0.0;
                Smam[n] = 0.0;
              }

            /* Ratio of fissile to total flux */

            if (flx[n] > 0.0)
              fflx[n] = fflx[n]/flx[n];
            else
              fflx[n] = 0.0;
          }

      /***********************************************************************/

      /* Next universe */

      gcu = NextItem(gcu);
    }

  /* Do micro-group calculation */

  MicroCalc();

  /* Clear data */

  ClearMicroGroupXS();
}

/*****************************************************************************/
