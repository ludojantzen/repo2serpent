/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processnubardata.c                             */
/*                                                                           */
/* Created:       2010/02/05 (JLe)                                           */
/* Last modified: 2020/01/29 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes ENDF nubar data                                    */
/*                                                                           */
/* Comments: - Adopted from energydistribution.c in Serpent 1.1.12           */
/*                                                                           */
/*           - NOTE: Aseta fissile flagi täällä                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessNubarData:"

/*****************************************************************************/

void ProcessNubarData(long rea)
{
  long nuc, ace, ptr, mt, loc0, loc1, fisse_ptr, n, i, ZAI;
  long L, LNU, NR, NP, NG, NXS[16], JXS[32];
  double lambda, P, *XSS;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  if (((mt < 18) || (mt > 21)) && (mt != 38))
    {
      /* Reset pointers */

      WDB[rea + REACTION_PTR_TNUBAR] = NULLPTR;
      WDB[rea + REACTION_PTR_DNUBAR] = NULLPTR;
      WDB[rea + REACTION_PTR_PREC_LIST] = NULLPTR;

      /* Exit subroutine */

      return;
    }

  /* Pointer to nuclide data */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Link pointers if not total or first-chance fission */
  /* NOTE: Tää olettaa että noi on järjestyksessä */

  if ((mt != 18) && (mt != 19))
    {
      /* Get pointer to primary fission channel */

      ptr = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Link pointers */

      WDB[rea + REACTION_PTR_TNUBAR] = RDB[ptr + REACTION_PTR_TNUBAR];
      WDB[rea + REACTION_PTR_DNUBAR] = RDB[ptr + REACTION_PTR_DNUBAR];
      WDB[rea + REACTION_PTR_PREC_LIST] = RDB[ptr + REACTION_PTR_PREC_LIST];

      /* Exit subroutine */

      return;
    }

  /* Pointer to fission energy deposition data */

  fisse_ptr = (long)RDB[nuc + NUCLIDE_PTR_FISSE_DATA];

  /* Pointer to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];

  /* Read data to NXS array */

  ptr = (long)ACE[ace + ACE_PTR_NXS];

  for (n = 0; n < 16; n++)
    NXS[n] = (long)ACE[ptr++];

  /* Read data to JXS array */

  ptr = (long)ACE[ace + ACE_PTR_JXS];

  for (n = 0; n < 32; n++)
    JXS[n] = (long)ACE[ptr++];

  /* Get pointer to XSS array */

  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

  /***************************************************************************/

  /***** Total nubar data ****************************************************/

  /* Check existence and type of nubar data */

  if ((L = JXS[1] - 1) < 0)
    Die(FUNCTION_NAME, "Nuclide %s, mt %ld has no prompt fission nubar data",
        GetText(nuc + NUCLIDE_PTR_NAME), mt);

  /* If XSS[L] > 0, total nubar data is assumed. Otherwise get */
  /* pointer to total nubar data. */

  if ((long)XSS[L] == 0)
    Die(FUNCTION_NAME, " XSS[L] = 0");
  else if ((long)XSS[L] < 0)
    L = (long)JXS[1] + (long)fabs(XSS[JXS[1] - 1]);

  /* Allocate memory for data block */

  loc0 = NewItem(rea + REACTION_PTR_TNUBAR, NUBAR_BLOCK_SIZE);

  /* Get data type */

  LNU = (long)XSS[L];

  /* Put type */

  WDB[loc0 + NUBAR_DATA_TYPE] = (double)LNU;

  /* Allocate memory for previous value */

  AllocValuePair(loc0 + NUBAR_PTR_PREV_VAL);

  /* Some libraries have erroneous data (L = L - 1). Try this. */

  if ((LNU != 1) && (LNU != 2))
    {
      /* Tää on joku ylimääräinen? */

      Die(FUNCTION_NAME, " LNU = %ld", LNU);

      L = L - 1;
      LNU = (long)XSS[L];
    }

  /* Check data type */

  if (LNU == 1)
    {
      /***********************************************************************/

      /***** Polynomial data *************************************************/

      /* Get number of coefficients */

      NP = (long)XSS[L + 1];

      /* Allocate memory for data */

      loc1 = ReallocMem(DATA_ARRAY, NP + 1);

      /* Put pointer */

      WDB[loc0 + NUBAR_PTR_POLY_DATA] = (double)loc1;

      /* Put number of points */

      WDB[loc1++] = (double)NP;

      /* Read coefficients */

      for (n = 0; n < NP; n++)
        WDB[loc1++] = XSS[L + 2 + n];

      /* Check if fission energy deposition data is defined */
      /* Note: This is now at E = 0, but maybe it should be at the threshold */
      /*       for fertile nuclides */

      if (fisse_ptr > VALID_PTR)
        WDB[fisse_ptr + FISSE_DATA_NUBAR0] = XSS[L + 2];

      /***********************************************************************/
    }
  else if (LNU == 2)
    {
      /***********************************************************************/

      /***** Tabular data ****************************************************/

      /* Check number of interpolation regions */

      if ((NR = (long)XSS[L + 1]) != 0)
        {
          /* Calculate number of non linear-linear regions */

          NP = 0;

          for (n = 0; n < NR; n++)
            {
              if (((long)XSS[L + 2 + NR + n] != 1) &&
                  ((long)XSS[L + 2 + NR + n] != 2))
                Die(FUNCTION_NAME, "Error in nubar data (%s mt %ld)",
                    GetText(nuc + NUCLIDE_PTR_NAME), mt);
              else if ((long)XSS[L + 2 + NR + n] == 1)
                NP++;
            }

          if (NP > 0)
            Warn(FUNCTION_NAME, "Nubar data converted to lin-lin (%s mt %ld)",
                 GetText(nuc + NUCLIDE_PTR_NAME), mt);
        }

      if (NR > 0)
        Die(FUNCTION_NAME, "NR = %ld %s", NR,
            GetText(nuc + NUCLIDE_PTR_NAME));

      /* Get number of energy points */

      NP = (long)XSS[L + 2*NR + 2];

      /* Create energy grid */

      loc1 = MakeEnergyGrid(NP, 0, 0, -1, &XSS[L + 2*NR + 3],
                            EG_INTERP_MODE_LIN);

      WDB[loc0 + NUBAR_PTR_EGRID] = (double)loc1;

      /* Allocate memory for points */

      loc1 = ReallocMem(DATA_ARRAY, NP);

      /* Put pointer */

      WDB[loc0 + NUBAR_PTR_PTS] = (double)loc1;

      /* Copy data */

      memcpy(&WDB[loc1], &XSS[L + 2*NR + 3 + NP], NP*sizeof(double));

      /* Check values (may be zero) */

      for (n = 0; n < NP; n++)
        {
          /* NOTE: Tähän tarvii noi kaarisulkeet */

          CheckValue(FUNCTION_NAME, "tnu", "", RDB[loc1 + n], 0.0, 10.0);
        }

      /* Check if fission energy deposition data is defined */

      if (fisse_ptr > VALID_PTR)
        WDB[fisse_ptr + FISSE_DATA_NUBAR0] = RDB[loc1];
    }

  /* Exit if delayed neutron data doesn't exist */

  if ((long)JXS[23] - 1 < 0)
    return;

  /***************************************************************************/

  /***** Precursor data ******************************************************/

  /* Get ZAI */

  ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

  /* Check that global group number is set */

  if ((long)RDB[DATA_PRECURSOR_GROUPS] == 0)
    Die(FUNCTION_NAME, "Number of precursor groups not set");

  /* Get number of precursor groups */

  if ((NG = NXS[7]) > MAX_PRECURSOR_GROUPS)
    Die(FUNCTION_NAME, "Too many precursor groups (%ld, max = %d)",
        NG, MAX_PRECURSOR_GROUPS);
  else if (NG != (long)RDB[DATA_PRECURSOR_GROUPS])
    {
      /* Multiple precursor group structures in calculation */

#ifdef DEBUG

      /* Print warning for all in debug mode */

      Note(0, "Delayed neutron data discarded for %s (mismatch in groups)",
           GetText(nuc + NUCLIDE_PTR_NAME));
#else

      /* Check significant actinides */

      if ((ZAI == 922330) || (ZAI == 922340) || (ZAI == 922350) ||
          (ZAI == 922360) || (ZAI == 922380) || (ZAI == 932370) ||
          (ZAI == 942380) || (ZAI == 942390) || (ZAI == 942400) ||
          (ZAI == 942410) || (ZAI == 942420) || (ZAI == 952410) ||
          (ZAI == 952420) || (ZAI == 952421) || (ZAI == 952430) ||
          (ZAI == 902320) || (ZAI == 922320))
        Note(0, "Delayed neutron data discarded for %s (mismatch in groups)",
             GetText(nuc + NUCLIDE_PTR_NAME));

#endif

      /* Skip data */

      return;
    }

  /* Get pointer to precursor data */

  L = (long)JXS[24] - 1;

  /* Loop over precursor groups */

  for (n = 0; n < NG; n++)
    {
      /* Create item */

      loc0 = NewItem(rea + REACTION_PTR_PREC_LIST, PREC_BLOCK_SIZE);

      /* Put index */

      WDB[loc0 + PREC_IDX] = (double)(n + 1);

      /* Get decay constant for group and convert to seconds */

      if ((lambda = XSS[L]/SHAKE) == 0.0)
        Warn(FUNCTION_NAME, "Zero lambda");

      CheckValue(FUNCTION_NAME, "lambda", "", lambda, 0.0, 11.0);

      /* Set value decay constant */

      WDB[loc0 + PREC_LAMBDA] = lambda;

      /* Process energy distribution */

      ProcessEDistributions(rea, loc0);

      /* Get number of interpolation regions. The interpolation  */
      /* parameters are not read, which means that linear-linear */
      /* interpolation is always used. */

      NR = (long)XSS[L + 1];

      /* Get number of points */

      NP = (long)XSS[L + 2 + 2*NR];
      CheckValue(FUNCTION_NAME, "NP", " (prec)", NP, 0, 500);

      /* Create energy grid */

      loc1 = MakeEnergyGrid(NP, 0, 0, -1, &XSS[L + 3 + 2*NR],
                            EG_INTERP_MODE_LIN);

      WDB[loc0 + PREC_PTR_EGRID] = (double)loc1;

      /* Allocate memory for probabilities */

      loc1 = ReallocMem(DATA_ARRAY, NP);

      /* Set pointer */

      WDB[loc0 + PREC_PTR_PTS] = (double)loc1;

      /* Loop over probabilities */

      for (i = 0; i < NP; i++)
        {
          /* Get value */

          P = XSS[L + 3 + 2*NR + NP + i];
          CheckValue(FUNCTION_NAME, "P", " (prec)", P, 1E-4, 0.9);

          /* Set value */

          WDB[loc1++] = P;
        }

      /* Update pointer */

      L = L + 3 + 2*NR + 2*NP;
    }

  /* Close list */

  loc0 = (long)RDB[rea + REACTION_PTR_PREC_LIST];
  CloseList(loc0);

  /****************************************************************************/

  /***** Delayed nubar data **************************************************/

  /* Get pointer and check if data exists */

  if ((L = (long)JXS[23] - 1) == 0)
    Die(FUNCTION_NAME, "JXS[23] - 1 = 0");
  else if (L < 0)
    return;

  /* Check data type (must always be tabular) */

  if ((long)XSS[L] != 2)
    Die(FUNCTION_NAME, "Invalid delayed nubar data type (%s mt %ld)",
        GetText(nuc + NUCLIDE_PTR_NAME), mt);

  /* Check number of interpolation regions */

  if ((long)XSS[L + 1] != 0)
    Die(FUNCTION_NAME, "Interpolation in delayed nubar data (%s mt %ld)",
        GetText(nuc + NUCLIDE_PTR_NAME), mt);

  /* Allocate memory for data block */

  loc0 = NewItem(rea + REACTION_PTR_DNUBAR, NUBAR_BLOCK_SIZE);

  /* Put type */

  WDB[loc0 + NUBAR_DATA_TYPE] = 2.0;

  /* Allocate memory for previous value */

  AllocValuePair(loc0 + NUBAR_PTR_PREV_VAL);

  /* Get number of energy points */

  NP = (long)XSS[L + 2];

  /* Create energy grid */

  loc1 = MakeEnergyGrid(NP, 0, 0, -1, &XSS[L + 3], EG_INTERP_MODE_LIN);
  WDB[loc0 + NUBAR_PTR_EGRID] = (double)loc1;

  /* Allocate memory for points */

  loc1 = ReallocMem(DATA_ARRAY, NP);

  /* Put pointer */

  WDB[loc0 + NUBAR_PTR_PTS] = (double)loc1;

  /* Copy data */

  memcpy(&WDB[loc1], &XSS[L + 3 + NP], NP*sizeof(double));

  /* Check values */

  for (n = 0; n < NP; n++)
    {
      /* NOTE: Tähän tarvii noi kaarisulkeet */

      CheckValue(FUNCTION_NAME, "dnu", "", RDB[loc1 + n], 0.0, 0.11);
    }

  /****************************************************************************/
}

/*****************************************************************************/
