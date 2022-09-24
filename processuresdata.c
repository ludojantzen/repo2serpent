/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processuresdata.c                              */
/*                                                                           */
/* Created:       2011/01/08 (JLe)                                           */
/* Last modified: 2018/10/01 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes unresolved resonance probability table data        */
/*                                                                           */
/* Comments: - Serpent 1 hakee ensimmäisen energiapisteen etukäteen          */
/*             generoidusta listasta, mikä voi nopeuttaa sämpläystä. Samaa   */
/*             voisi harkita tähän.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessUresData:"

/*****************************************************************************/

void ProcessUresData(long nuc)
{
  long L, N, M, INT, ILF, IOA, IFF, loc0, neg, n, m, i, j, rea, rea0, ptr;
  long urs, pte, ptp, ace, JXS[32];
  double *E0, f, fmax, *XSS, sum;

  /* Check ures flag */

  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED))
    return;

  /* Check TMS */

  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
    {
      /* Reset ures flag */

      ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_USED);

      /* Exit */

      return;
    }

  /* Pointer to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];

  /* Read data to JXS array */

  ptr = (long)ACE[ace + ACE_PTR_JXS];

  for (n = 0; n < 32; n++)
    JXS[n] = (long)ACE[ptr++];

  /* Get pointer to XSS array */

  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

  /* Get pointer to UNR block */

  if ((L = JXS[22] - 1) < 0)
    Die(FUNCTION_NAME, "Ures flag set but ACE block not found");

  /***************************************************************************/

  /***** Read pointers and flags *********************************************/

  /* Read number of incident energies */

  N = (long)XSS[L];

  /* Read number of probabilities */

  M = (long)XSS[L + 1];

  /* Read interpolation parameter (2 = lin-lin, 5 = log-log) */

  INT = (long)XSS[L + 2];

  /* Check interpolation scheme, (2 = lin-lin, 5 = log-log) */

  if ((INT != 2) && (INT != 5))
    Die(FUNCTION_NAME,
        "Unrecognized URES interpolation scheme: INT = %ld (nuclide %s)",
        INT, GetText(nuc + NUCLIDE_PTR_NAME));

  /* XSS[L + 3] is ILF (Inelastic competition flag), XSS[L + 4] is IOA */
  /* (other absorption flag). Neither is used here. */

  ILF = (long)XSS[L + 3];
  IOA = (long)XSS[L + 4];

#ifdef DEBUG

  if (ILF > -1)
    Warn(FUNCTION_NAME, "ILF = %ld\n", ILF);

  if (IOA > -1)
    Warn(FUNCTION_NAME, "IOA = %ld\n", IOA);

#endif

  /* Other factors flag (0 = cross sections, 1 = factors of smooth) */

  IFF = (long)XSS[L + 5];

  /* Get pointer to incident energies */

  E0 = &XSS[L + 6];

  if (1 != 2)
    {
      /* Check points and discard data if negative values are found. This  */
      /* is the case for some natural elements and the values seem to make */
      /* sense without the minus sign. (27.9.2009) */

      for (n = 0; n < N; n++)
        if (E0[n] < 0.0)
          {
            /* Print warning */

            Note(0, "Unable to process ures data for %s, tables discarded",
                 GetText(nuc + NUCLIDE_PTR_NAME));

            /* Reset ures options */

            ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_AVAIL);
            ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_USED);

            /* Update counters */

            WDB[DATA_URES_AVAIL] = RDB[DATA_URES_AVAIL] - 1.0;
            WDB[DATA_URES_USED] = RDB[DATA_URES_USED] - 1.0;

            /* Exit subroutine */

            return;
          }
    }
  else
    {
      /* Or take absolute values (28.9.2009). PURR sets negative energies */
      /* for materials consisting of multiple isotopes with overlapping   */
      /* ures regions. */

      for (n = 0; n < N; n++)
        E0[n] = fabs(E0[n]);
    }

  /* Set energy limits for nuclide */

  WDB[nuc + NUCLIDE_URES_EMIN] = E0[0];
  WDB[nuc + NUCLIDE_URES_EMAX] = E0[N - 1];

  /* Check limits */

  CheckValue(FUNCTION_NAME, "Emin", "", E0[0], 1E-7, 20.0);
  CheckValue(FUNCTION_NAME, "Emax", "", E0[N - 1], E0[0], 20.0);

  /* Compare to global boundaries */

  if (RDB[nuc + NUCLIDE_URES_EMIN] < RDB[DATA_URES_EMIN])
    WDB[DATA_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];

  if (RDB[nuc + NUCLIDE_URES_EMAX] > RDB[DATA_URES_EMAX])
    WDB[DATA_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMAX];

  /* Allocate memory for random number (used for all reaction modes of */
  /* nuclides with same ZAI) */

  if ((long)RDB[nuc + NUCLIDE_PTR_URES_RND] < VALID_PTR)
    {
      /* Allocate memory */

      AllocValuePair(nuc + NUCLIDE_PTR_URES_RND);

      /* Loop over nuclides */

      ptr = (long)RDB[DATA_PTR_NUC0];
      while (ptr > VALID_PTR)
        {
          /* Check ZAI and copy pointer */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == (long)RDB[ptr + NUCLIDE_ZAI])
            WDB[ptr + NUCLIDE_PTR_URES_RND] = RDB[nuc + NUCLIDE_PTR_URES_RND];

          /* Next */

          ptr = NextItem(ptr);
        }
    }

  /* Update pointer */

  L = L + 6 + N;

  /***************************************************************************/

  /***** Store reaction data *************************************************/

  /* Make energy grid */

  pte = MakeEnergyGrid(N, 0, 0, -1, E0, EG_INTERP_MODE_LIN);

  /* Allocate memory for probabilities */

  ptp = ReallocMem(DATA_ARRAY, N*M);

  /* Store probabilities */

  i = 0;
  for (n = 0; n < N; n++)
    {
      for (m = 0; m < M; m++)
        {
          WDB[ptp + i++] = XSS[L + n*6*M + m];

          /* Check ascending order */

          if (m > 0)
            if (XSS[L + n*6*M + m] < XSS[L + n*6*M + m - 1])
              Die(FUNCTION_NAME, "Probabilities not in ascending order");
        }

      /* Check last point */

      if (XSS[L + n*6*M + M - 1] != 1.0)
        Die(FUNCTION_NAME, "Probability distribution not normalized");
    }

  /* Loop over reaction modes */

  for (j = 0; j < 5; j++)
    {
      /* Get pointer to reaction data */

      rea = -1;

      if(j == 0)
        rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      else if(j == 1)
        rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
      else if (j == 2)
        rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
      else if (j == 3)
        rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
      else if (j == 4)
        rea = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS];

      /* Check pointer */

      if (rea < VALID_PTR)
        continue;

      /* Set minimum and maximum energy */

      WDB[rea + REACTION_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];
      WDB[rea + REACTION_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMAX];

      /* Allocate memory for previous unadjusted value */

      AllocValuePair(rea + REACTION_PTR_PREV_URES_XS0);

      /* Allocate memory for ures data block */

      urs = NewItem(rea + REACTION_PTR_URES, URES_BLOCK_SIZE);

      /* Allocate memory for previous factor */

      AllocValuePair(urs + URES_PTR_PREV_FACT);

      /* Set pointer to grid */

      WDB[urs + URES_PTR_EGRID] = (double)pte;

      /* Set number of probabilities, pointer to grid and interpolation */

      WDB[urs + URES_NP] = (double)M;
      WDB[urs + URES_PTR_PROB] = (double)ptp;
      WDB[urs + URES_INT] = (double)INT;

      /* Link pointer to random number */

      WDB[urs + URES_PTR_RND] = RDB[nuc + NUCLIDE_PTR_URES_RND];

      /* Allocate memory for rnd check */

      AllocValuePair(urs + URES_PTR_RND_CHK);

      /* Store interpolation flag */

      WDB[urs + URES_IFF] = (double)IFF;

      /* Allocate memory for factors */

      ptr = ReallocMem(DATA_ARRAY, N*M);

      /* Set pointer */

      WDB[urs + URES_PTR_FACT] = (double)ptr;

      /* Store points */

      i = 0;
      for (n = 0; n < N; n++)
        for (m = 0; m < M; m++)
          {
            /* Get value */

            f = XSS[L + n*6*M + m + (j + 1)*M];

            /* Convert MEVs to Joules in KERMA data if necessary */

            if ((IFF == 0) && (j == 4))
              f = f*MEV;

            /* Check */

            if ((IFF == 0) && ((f < -10.0) || (f > MAX_XS)))
              {
                /* Some 0K nuclides have NaN in ures data */

                if (isnan(f))
                  {
                    /* Print warning */

                    Warn(FUNCTION_NAME,
                         "Nuclide %s has NaN's in ures ACE data",
                         GetText(nuc + NUCLIDE_PTR_NAME));

                    /* Reset ures flag */

                    ResetOption(nuc + NUCLIDE_TYPE_FLAGS,
                                NUCLIDE_FLAG_URES_USED);

                    /* Exit */

                    return;
                  }
                else
                  Die(FUNCTION_NAME, "Invalid table value %1.5E (%s IFF = 0)",
                      f, GetText(nuc + NUCLIDE_PTR_NAME));
              }
            else if (IFF == 1)
              {
                if ((f < -2E-1) || (f > 10000.0))
                  {
                    /* Some 0K nuclides have NaN in ures data */

                    if (isnan(f))
                      {
                        /* Print warning */

                        Warn(FUNCTION_NAME,
                             "Nuclide %s has NaN's in ures ACE data",
                             GetText(nuc + NUCLIDE_PTR_NAME));

                        /* Reset ures flag */

                        ResetOption(nuc + NUCLIDE_TYPE_FLAGS,
                                    NUCLIDE_FLAG_URES_USED);

                        /* Exit */

                        return;
                      }
                    /* Large factors are possible for non-fission KERMA in EDEP_MODE_NEUTRON_PHOTON */
                    else if (j == 4 && (long)RDB[DATA_EDEP_MODE] == EDEP_MODE_NEUTRON_PHOTON)
                      Warn(FUNCTION_NAME,
                           "Abnormal table value %.6E (%s IFF = 1)",
                           f,GetText(nuc + NUCLIDE_PTR_NAME));
                    else
                      Die(FUNCTION_NAME,
                          "Invalid table value %1.5E (%s IFF = 1)", f,
                          GetText(nuc + NUCLIDE_PTR_NAME));
                  }
                else if (f < 0.0)
                  Warn(FUNCTION_NAME,
                       "Negative table value %1.5E (%s IFF = 1)", f,
                       GetText(nuc + NUCLIDE_PTR_NAME));
              }

            /* Add point */

            WDB[ptr + i++] = f;
          }
    }

  /****************************************************************************/

  /***** Check and adjust totals (IFF = 0 only) *******************************/

  if (IFF == 0)
    {
      /* Reset number of negative points */

      neg = 0.0;

      /* Loop over points */

      for (i = 0; i < N*M; i++)
        {
          /* Reset sum */

          sum = 0.0;

          /* Loop over other reaction modes */

          for (j = 0; j < 3; j++)
            {
              /* Get pointer to reaction data */

              rea = -1;

              if(j == 0)
                rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
              else if (j == 1)
                rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
              else if (j == 2)
                rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];

              /* Check pointer */

              if (rea < VALID_PTR)
                continue;

              /* Pointer to ures data */

              urs = (long)RDB[rea + REACTION_PTR_URES];
              CheckPointer(FUNCTION_NAME, "(urs)", DATA_ARRAY, urs);

              /* Pointer to factors */

              ptr  = (long)RDB[urs + URES_PTR_FACT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Check negative */

              if (RDB[ptr + i] < 0.0)
                {
                  /* Adjust to zero */

                  WDB[ptr + i] = 0.0;

                  /* Add count */

                  neg++;
                }

              /* Add to sum */

              sum = sum + RDB[ptr + i];
            }

          /* Pointer to total cross section */

          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Pointer to ures data */

          urs = (long)RDB[rea + REACTION_PTR_URES];
          CheckPointer(FUNCTION_NAME, "(urs)", DATA_ARRAY, urs);

          /* Pointer to factors */

          ptr = (long)RDB[urs + URES_PTR_FACT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Put total */

          WDB[ptr + i] = sum;
        }

      /* Print warning */

      if (neg > 0)
        Warn(FUNCTION_NAME, "%ld negative values in ures data", neg);
    }

  /****************************************************************************/

  /***** Calculate maximum fractions ******************************************/

  /* Tää siirrettiin tuolta ylempää 26.8.2016 / 2.1.27 (JLe) */

  /* Loop over reaction modes */

  for (j = 0; j < 5; j++)
    {
      /* Get pointer to reaction data */

      rea = -1;

      if(j == 0)
        rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      else if(j == 1)
        rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
      else if (j == 2)
        rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
      else if (j == 3)
        rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
      else if (j == 4)
        rea = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS];

      /* Check pointer */

      if (rea < VALID_PTR)
        continue;

      /* Pointer to data */

      urs = (long)RDB[rea + REACTION_PTR_URES];
      CheckPointer(FUNCTION_NAME, "(urs)", DATA_ARRAY, urs);

      /* Pointer to table */

      ptr = (long)RDB[urs + URES_PTR_FACT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Allocate memory for maximum factors */

      loc0 = ReallocMem(DATA_ARRAY, N);

      /* Set pointer */

      WDB[urs + URES_PTR_MAXF] = (double)loc0;

      /* Loop over values */

      for (n = 0; n < N; n++)
        {
          fmax = RDB[loc0 + n];
          for (m = 0; m < M; m++)
            {
              /* Compare to value in table m at energy n */

              if ((f =  RDB[ptr + n*M + m]) > fmax)
                fmax = f;

              /* Compare to value in table m at energy n + 1 */

              if (n < N - 1)
                if ((f = RDB[ptr + (n + 1)*M + m]) > fmax)
                  fmax = f;

              /* Compare to value in table m at energy n - 1 */
              /* 14.3.2018: This was in the older implementation for some reason. */
              /*
              if (n > 0)
                if ((f = RDB[ptr + (n - 1)*M + m]) > fmax)
                  fmax = f;
              */
            }
          WDB[loc0 + n] = fmax;
        }
    }

  /***************************************************************************/

  /***** Sum of absorptions **************************************************/

  /* Pointer to absorption cross section */

  rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to (n,gamma) */

  rea0 = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

  /* Allocate memory for previous unadjusted value */
  /*
  AllocValuePair(rea + REACTION_PTR_PREV_URES_XS0);
  */
  /* Copy data */

  WDB[rea + REACTION_PTR_URES] = RDB[rea0 + REACTION_PTR_URES];
  WDB[rea + REACTION_URES_EMIN] = RDB[rea0 + REACTION_URES_EMIN];
  WDB[rea + REACTION_URES_EMAX] = RDB[rea0 + REACTION_URES_EMAX];
  WDB[rea + REACTION_URES_MAX_N0] = RDB[rea0 + REACTION_URES_MAX_N0];
  WDB[rea + REACTION_URES_MAX_NP] = RDB[rea0 + REACTION_URES_MAX_NP];


  WDB[rea + REACTION_PTR_PREV_URES_XS0] = RDB[rea0 + REACTION_PTR_PREV_URES_XS0];

  /* Check for non-(n,gamma) capture reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check implicit capture, ty and mt */

      if (((long)RDB[DATA_OPT_IMPL_CAPT] == YES) &&
          ((long)RDB[rea + REACTION_TY] == 0) &&
          ((long)RDB[rea + REACTION_MT] != 102) &&
          ((long)RDB[rea + REACTION_MT] > 99) &&
          ((long)RDB[rea + REACTION_MT] < 200))
        Warn(FUNCTION_NAME, "Non-(n,gamma) reaction mode (%s %ld)",
             GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT]);

      /* Next reaction */

      rea = NextItem(rea);
    }

  /***************************************************************************/
}

/*****************************************************************************/
