/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nfkerma.c                                      */
/*                                                                           */
/* Created:       2018/03/22 (RTu)                                           */
/* Last modified: 2018/11/09 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates non-fission KERMA                                 */
/*                                                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NFKERMA:"

/*****************************************************************************/

void NFKERMA(long nuc, double *fiss_kerma)
{
  long i, j, ptr, rea, I0, NE, ace, NES, n, m, mt, N, M, IFF, ptr0, c0, c1, c2;
  double nfiss_xs, t, f, k, nfiss_kerma;
  double *tot_xs, *fiss_xs, *XSS, *NXS, *JXS, *E, *E0;
  double *inf_fiss_xs, *inf_fiss_kerma, *inf_tot_kerma, *tot_kerma;
  long fiss_ptrs[] = {-1, -1, -1, -1, -1};

  /* Get pointers to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
  CheckPointer(FUNCTION_NAME, "ace", ACE_ARRAY, ace);

  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];
  NXS = &ACE[(long)ACE[ace + ACE_PTR_NXS]];
  JXS = &ACE[(long)ACE[ace + ACE_PTR_JXS]];

  /* Number of energy points */

  NES = (long)NXS[2];

  /* Allocate memory for total fission cross section */

  fiss_xs = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

  /* Find fission reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

  while (rea > VALID_PTR)
    {
      /* Get mt */

      mt = RDB[rea + REACTION_MT];

      /* Store pointer to fission XS data (MTs 18, 19, 20, 21 and 38) */

      if ((mt > 17) && (mt < 22))
        {
          fiss_ptrs[mt - 18] = rea;

          if (mt == 18)
            break;
        }
      else if (mt == 38)
        {
          fiss_ptrs[4] = rea;
        }

      /* Next reaction */

      rea = NextItem(rea);
    }

  /* Calculate total fission cross section */

  for (j = 0; j < 5; j++)
    {
      rea = fiss_ptrs[j];
      if (rea > VALID_PTR)
        {
          /* Pointer to xs data */

          ptr = (long)RDB[rea + REACTION_PTR_XS];

          /* Index to first point and number of points */

          I0 = (long)RDB[rea + REACTION_XS_I0];
          NE = (long)RDB[rea + REACTION_XS_NE];

          /* Check value */

          CheckValue(FUNCTION_NAME, "NE", "", I0 + NE, 2, NES);

          /* Loop over points and add to sum */

          for (i = 0; i < NE; i++)
            {
              /* Add value to total */

              fiss_xs[I0 + i] = fiss_xs[I0 + i] + XSS[ptr + i];
            }

          /* Total fission cross section MT18 found, break loop */

          if (j == 0)
            break;
        }
    }

  /* Total cross section array */

  tot_xs = &XSS[(long)JXS[0] - 1 + NES];

  /* KERMA reaction pointer */

  rea = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS];
  CheckPointer(FUNCTION_NAME, "rea", ACE_ARRAY, rea);

  /* KERMA xs array */

  ptr0 = (long)RDB[rea + REACTION_PTR_XS];
  CheckPointer(FUNCTION_NAME, "ptr0", ACE_ARRAY, ptr0);

  tot_kerma = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

   for (i = 0; i < NES; i++)
     tot_kerma[i] = XSS[ptr0 + i]*tot_xs[i];

  /* Modify KERMA ures data for non-fission KERMA */

  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED)
    {
      /* Pointer to UNR block */

      if ((ptr = (long)JXS[22] - 1) < 0)
        Die(FUNCTION_NAME, "Ures flag set but ACE block not found");

      /* Number of incident energies */

      N = (long)XSS[ptr];

      /* Number of probabilities */

      M = (long)XSS[ptr + 1];

      /* Factors flag */

      IFF = (long)XSS[ptr + 5];

      /* Get pointer to incident energies */

      E = &XSS[ptr + 6];

      /* Update pointer to the start of probability tables */

      ptr = ptr + 6 + N;

      /* Get pointer to ACE energy grid */

      E0 = &XSS[(long)JXS[0] - 1];

      /* Interpolate fission cross section to URES energies */

      inf_fiss_xs = (double *)Mem(MEM_ALLOC, N, sizeof(double));

      n = InterpolateData(E, inf_fiss_xs, N, E0, fiss_xs,
                          NES, 0, NULL, NULL, YES);

      if (n > 0)
        Warn(FUNCTION_NAME,
             "%ld negative fission cross sections (%s)",
             n, GetText(nuc + NUCLIDE_PTR_NAME));

      /* Interpolate fission KERMA to URES energies */

      inf_fiss_kerma = (double *)Mem(MEM_ALLOC, N, sizeof(double));

      n = InterpolateData(E, inf_fiss_kerma, N, E0, fiss_kerma,
                          NES, 0, NULL, NULL, YES);

      if (n > 0)
        Warn(FUNCTION_NAME,
             "%ld negative fission kerma cross sections (%s)",
             n, GetText(nuc + NUCLIDE_PTR_NAME));

      /* Interpolate total KERMA to URES energies */

      inf_tot_kerma = (double *)Mem(MEM_ALLOC, N, sizeof(double));

      n = InterpolateData(E, inf_tot_kerma, N, E0, tot_kerma,
                          NES, 0, NULL, NULL, YES);

      if (n > 0)
        Warn(FUNCTION_NAME,
             "%ld negative kerma cross sections (%s)",
             n, GetText(nuc + NUCLIDE_PTR_NAME));

      /* Init counters */

      c0 = 0;
      c1 = 0;
      c2 = 0;

      /* Loop over incident energies and bins */

      for (n = 0; n < N; n++)
        for (m = 0; m < M; m++)
          {
            /* Get total cross section URES value */

            t = XSS[ptr + n*6*M + M + m];

            /* Get fission cross section URES value */

            f = XSS[ptr + n*6*M + 3*M + m];

            /* Get original KERMA URES value */

            k = XSS[ptr + n*6*M + 5*M + m];

            /* Convert to non-fission KERMA URES value */

            if (fabs(t) < ZERO)
              {
                c0++;
                k = 0.0;
              }
            else
              {
                if (IFF == 1)
                  {
                    nfiss_kerma = inf_tot_kerma[n]-inf_fiss_kerma[n];
                    if (fabs(nfiss_kerma) < ZERO)
                      {
                        c1++;
                        k = 0.0;
                      }
                    else
                      {
                        k = (k*inf_tot_kerma[n]-f/t*inf_fiss_kerma[n])/nfiss_kerma;
                      }
                  }
                else
                  {
                    if (fabs(inf_fiss_xs[n]) < ZERO)
                      {
                        c2++;
                      }
                    else
                      k = k - f/t*inf_fiss_kerma[n]/inf_fiss_xs[n];
                  }
              }

            XSS[ptr + n*6*M + 5*M + m] = k;
          }

      if (c0 > 0)
        Warn(FUNCTION_NAME, "%ld zero total URES cross sections/factors in data (%s)",
             c0, GetText(nuc + NUCLIDE_PTR_NAME));

      if (c1 > 0)
        Warn(FUNCTION_NAME, "%ld zero non-fission KERMAs in data (%s)",
             c1, GetText(nuc + NUCLIDE_PTR_NAME));

      if (c2 > 0)
        Warn(FUNCTION_NAME, "%ld zero fission cross sections in data (%s)",
             c2, GetText(nuc + NUCLIDE_PTR_NAME));

      /* Free temporary arrays */

      Mem(MEM_FREE, inf_fiss_xs);
      Mem(MEM_FREE, inf_fiss_kerma);
      Mem(MEM_FREE, inf_tot_kerma);
    }

  /* Convert to non-fission kerma */

  n = 0;
  m = 0;

  for (i = 0; i < NES; i++)
    {
      nfiss_xs = tot_xs[i] - fiss_xs[i];
      if (fabs(nfiss_xs) < ZERO)
        {
          m++;
          XSS[ptr0 + i] = 0.0;
        }
      else
        {
          XSS[ptr0 + i] = (tot_kerma[i] - fiss_kerma[i])/nfiss_xs;
          if (XSS[ptr0 + i] < 0.0)
            n++;
        }
    }

  if (n > 0)
    Warn(FUNCTION_NAME,
         "%ld negative non-fission KERMAs in data (%s)",
         n, GetText(nuc + NUCLIDE_PTR_NAME));
  if (m > 0)
    Warn(FUNCTION_NAME,
         "%ld zero non-fission cross sections in data (%s)",
         m, GetText(nuc + NUCLIDE_PTR_NAME));

  /* Free temporary array */

  Mem(MEM_FREE, fiss_xs);
  Mem(MEM_FREE, tot_kerma);

}

/*****************************************************************************/
