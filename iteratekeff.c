/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratekeff.c                                  */
/*                                                                           */
/* Created:       2013/10/02 (JLe)                                           */
/* Last modified: 2020/03/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: K-eff iteration                                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterateKeff:"

/*****************************************************************************/

void IterateKeff()
{
  long ptr, mode, skip, ncyc, idx, fix, inuc, iso;
  double nsf, fiss, capt, nuxn, leak, L0, L1, val, f, f0, newf, adens;
  double keff, val0, C0, C1, relax;

  /* Check mode */

  if ((mode = (long)RDB[DATA_ITER_MODE]) == ITER_MODE_NONE)
    return;

  /* Get fix mode */

  fix = (long)RDB[DATA_ITER_FIX];

  /* Number of cycles and actual number of skip cycles (setoptimization.c) */

  ncyc = (long)RDB[DATA_ITER_NCYC];
  idx = (long)RDB[DATA_CYCLE_IDX];

  if (fix == YES)
    skip = (long)((RDB[DATA_CRIT_SKIP] - RDB[DATA_ITER_NCYC])/2.0);
  else
    skip = (long)(RDB[DATA_CRIT_SKIP] - RDB[DATA_ITER_NCYC]);

  /* Check cycles */

  if ((idx < skip) || ((fix == YES) && (idx > skip + ncyc)))
    return;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Collect MPI parallel data */

  CollectBuf();

  /* Check mode */

  if (mode == ITER_MODE_ALBEDO)
    {
      /***********************************************************************/

      /***** Albedo iteration ************************************************/

      /* Get k-eff */

      keff = RDB[DATA_ITER_KEFF];
      CheckValue(FUNCTION_NAME, "keff", "", keff, 0.1, 2.5);

      /* Fission nubar */

      ptr = (long)RDB[RES_TOT_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      nsf = BufVal(ptr, 0);

      /* Fission term */

      ptr = (long)RDB[RES_TOT_FISSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fiss = BufVal(ptr, 0);

      /* Total capture rate */

      ptr = (long)RDB[RES_TOT_CAPTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      capt = BufVal(ptr, 0);

      /* Scattering production rate */

      ptr = (long)RDB[RES_TOT_INLPRODRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      nuxn = BufVal(ptr, 0);

      /* Physical leakage rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      leak = BufVal(ptr, 0);

      /* Get previous albedo leakage rate */

      ptr = (long)RDB[RES_ALB_NEUTRON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      L0 = BufVal(ptr, 0);

      /* Calculate estimate for new albedo leakage rate */

      L1 = nsf/keff - capt - fiss - leak + nuxn;

      /* Avoid compiler warning */

      val = -1.0;

      /* Get previous value */

      if ((val0 = RDB[DATA_ITER_VAL]) < 0.0)
        {
          /* Not set, use initial guess */

          ptr = (long)RDB[RES_ANA_KEFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          if (Mean(ptr, 0) > 1.0)
            val = 0.999;
          else
            val = 1.001;
        }
      else if (L0 != 0.0)
        {
          /* Calculate new */

          val = (val0 - 1.0)*L1/L0 + 1.0;

          /* Add to statistics */

          ptr = (long)RDB[RES_ITER_VAL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          /* Fix value for last iteration */

          if ((fix == YES) && (idx > skip + ncyc))
            val = Mean(ptr, 0);
        }
      else
        Die(FUNCTION_NAME, "L0 == 0");

      /* Check (mustakaan reuna ei välttämättä anna tarpeeksi vuotoa) */

      if (val < 0.0)
        Error(0, "Albedo iteration failed");

      /* Put value */

      WDB[DATA_ITER_VAL] = val;

      /* Put albedos */

      if ((f = RDB[DATA_ITER_ALB_F1]) > 0.0)
        WDB[DATA_GEOM_ALBEDO1] = val*f;

      if ((f = RDB[DATA_ITER_ALB_F2]) > 0.0)
        WDB[DATA_GEOM_ALBEDO2] = val*f;

      if ((f = RDB[DATA_ITER_ALB_F3]) > 0.0)
        WDB[DATA_GEOM_ALBEDO3] = val*f;

      /***********************************************************************/
    }
  else if (mode == ITER_MODE_NUCLIDE)
    {
      /**********************************************************************/
      /* Figure out factor that the atomic densities have to be scaled with */
      /**********************************************************************/

      /* Get target k-eff */

      keff = RDB[DATA_ITER_KEFF];
      CheckValue(FUNCTION_NAME, "keff", "", keff, 0.1, 2.5);

      /* Fission nubar */

      ptr = (long)RDB[RES_TOT_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      nsf = BufVal(ptr, 0);

      /* Fission term */

      ptr = (long)RDB[RES_TOT_FISSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fiss = BufVal(ptr, 0);

      /* Total capture rate */

      ptr = (long)RDB[RES_TOT_CAPTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      capt = BufVal(ptr, 0);

      /* Scattering production rate */

      ptr = (long)RDB[RES_TOT_INLPRODRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      nuxn = BufVal(ptr, 0);

      /* Physical leakage rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      leak = BufVal(ptr, 0);

      /* Get current nuclide absorption rate */

      ptr = (long)RDB[RES_TOT_ITER_NUC_ABSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      C0 = BufVal(ptr, 0);

      /* Calculate estimate for required absorption rate */

      C1 = nsf/keff - capt - fiss - leak + nuxn + C0;

      /* Additional relaxation coefficient */

      relax = 0.3;

      /* Calculate scaling factor */

      if (C1 < 0)
        {
          /* The iterable nuclides should produce neutrons instead */
          /* of absorbing to reach the required keff */

          Warn(FUNCTION_NAME, "The k-eff iteration with nuclide concentrations"
               " may fail due to insufficient absorption by the adjustable "
               "nuclides.");

          /* Non-relaxed factor used for scoring the estimate */

          f0 = 0.0;

          /* Relaxed factor used for updating concentrations */

          f = 1.0 + relax*(-1.0);

        }
      else if (C0 != 0.0)
        {
          /* Calculate estimate for scaling factor */
          /* The atomic density is to be scaled with 1.0 + f*/

          f = C1/C0 - 1.0;

          /* Non-relaxed factor used for scoring the estimate */

          f0 = 1.0 + f;

          /* Relaxed factor used for updating concentrations */

          f = 1.0 + relax*f;

        }
      else
        {
          /* No captures to iterable nuclides in the previous cycle */
          /* Increase the concentrations tenfold (but relax) */

          /* Non-relaxed factor used for scoring the estimate */

          f0 = 9.0;

          /* Relaxed factor used for updating concentrations */

          f = 1.0 + relax*(9.0);

        }

      /* Use the unrelaxed value for scoring */

      val = RDB[DATA_ITER_VAL]*f0;

      /* Add to statistics */

      ptr = (long)RDB[RES_ITER_VAL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(val, ptr, 0);

      /* Calculate the relaxed value */

      newf = RDB[DATA_ITER_VAL]*f;

      /* Fix value for last iteration */

      if ((fix == YES) && (idx > skip + ncyc))
        {
          /* Get absolute scaling factor */

          newf = Mean(ptr, 0);

          /* Calculate relative scaling factor to be used now to */
          /* end up with the absolute one */

          f = newf/RDB[DATA_ITER_VAL];
        }

      /* Store value (absolute) */

      WDB[DATA_ITER_VAL] = newf;

      /******************************/
      /* Scale the atomic densities */
      /******************************/

      inuc = (long)RDB[DATA_PTR_ITER_NUC0];
      CheckPointer(FUNCTION_NAME, "(inuc)", DATA_ARRAY, inuc);

      while (inuc > VALID_PTR)
        {

          /* Loop over compositions to scale atomic densities */

          ptr = (long)RDB[inuc + ITER_NUCLIDE_PTR_COMPOSITION_LIST];

          while ((iso = (long)RDB[ptr]) > VALID_PTR)
            {

              /* Scale atomic density with factor f */

              adens = RDB[iso + COMPOSITION_ADENS];
              WDB[iso + COMPOSITION_ADENS] = adens*f;

              /* Get next composition */

              ptr++;
            }

          /* Get next nuclide */

          inuc = NextItem(inuc);
        }

      /**********************************************/
      /* Update maximum concentrations for majorant */
      /**********************************************/

      ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
      while (ptr > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + MAJORANT_EXTRA_TYPE] ==
              MAJORANT_EXTRA_NUCLIDE_ITER)
            WDB[ptr + MAJORANT_EXTRA_FRAC] = RDB[ptr + MAJORANT_EXTRA_FRAC]*f;

          /* Pointer to next */

          ptr = NextItem(ptr);
        }

    }
  else if (mode == ITER_MODE_USER)
    {
      /**********************************************************************/

      /***** User-defined iteration *****************************************/

      /* Call user-defined subroutine */

      UserIter(2);

      /**********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid iteration mode");
}

/*****************************************************************************/
