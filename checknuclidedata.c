/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checknuclidedata.c                             */
/*                                                                           */
/* Created:       2010/09/13 (JLe)                                           */
/* Last modified: 2019/11/13 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Performs various checs on nuclide data and prints some       */
/*              output.                                                      */
/*                                                                           */
/* Comments: - Lisää tarkistukset hajoamis branchingien ja fission yieldien  */
/*             summille.                                                     */
/*           - Toi stable nuclide decay modes -tarkistus ei oo ihan pätevä   */
/*             (Zn-70 on ENDFB/VII:ssa stabiili mutta sillä on beta-         */
/*             hajoamismoodi, näitä on paljon muitakin)                      */
/*           - Tää rutiini löytää vielä aika paljon bugeja.                  */
/*           - Osa Warn-käskyistä korvattu warningeilla että ajo menee läpi. */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckNuclideData:"

/*****************************************************************************/

void CheckNuclideData()
{
  long n, m, nuc, ptr, table, ace, i, j, rea, mt, yld, loc0, tgt, tot, nsab;
  long ZAI;
  double t;
  char outfile[MAX_STR], sym[MAX_STR], name[MAX_STR];
  FILE *fp;

  fprintf(outp, "Checking data and printing output...\n");

  /***************************************************************************/

  /***** Check basic data ****************************************************/

  /* Loop over nuclides */

  n = 0;
  while ((nuc = ListPtr((long)RDB[DATA_PTR_NUC0], n++)) > VALID_PTR)
    {
      /* Check used-flag */

      if (!((long)RDB[nuc + NUCLIDE_OPTIONS] & OPT_USED))
        Warn(FUNCTION_NAME, "Used flag not set for %s",
            GetText(nuc + NUCLIDE_PTR_NAME));

      /* Check type */

      if ((long)RDB[nuc + NUCLIDE_TYPE] < 1)
        Warn(FUNCTION_NAME, "Nuclide %s has no type",
            GetText(nuc + NUCLIDE_PTR_NAME));

      /* Check ZAI, ZA, Z, A and I */

      CheckValue(FUNCTION_NAME, "ZAI", "", RDB[nuc + NUCLIDE_ZAI], 10000,
                 1200000);
      CheckValue(FUNCTION_NAME, "ZA", "", RDB[nuc + NUCLIDE_ZA], 1000,
                 120000);
      CheckValue(FUNCTION_NAME, "Z", "", RDB[nuc + NUCLIDE_Z], 1, 119);
      CheckValue(FUNCTION_NAME, "A", "", RDB[nuc + NUCLIDE_A], 0, 290);
      CheckValue(FUNCTION_NAME, "I", "", RDB[nuc + NUCLIDE_I], 0, 3);

      /* Check that combinations match */

      if (RDB[nuc + NUCLIDE_ZAI] - 10000.0*RDB[nuc + NUCLIDE_Z] -
          10.0*RDB[nuc + NUCLIDE_A] - RDB[nuc + NUCLIDE_I] != 0.0)
        Warn(FUNCTION_NAME, "Mismatch in ZAI (%s)",
            GetText(nuc + NUCLIDE_PTR_NAME));

      if (RDB[nuc + NUCLIDE_ZA] - 1000.0*RDB[nuc + NUCLIDE_Z] -
          RDB[nuc + NUCLIDE_A] != 0)
        Warn(FUNCTION_NAME, "Mismatch in ZA (%s)",
            GetText(nuc + NUCLIDE_PTR_NAME));

      /* Check AW and AWR */

      CheckValue(FUNCTION_NAME,"AW", "", RDB[nuc + NUCLIDE_AW], 0.9, 290.0);
      CheckValue(FUNCTION_NAME,"AWR", "", RDB[nuc + NUCLIDE_AWR], 0.9,
                 290.0);

      /* Check ratio */

      if (fabs(RDB[nuc + NUCLIDE_AW]/RDB[nuc + NUCLIDE_AWR] - M_NEUTRON)
          > 1E-2)
        Warn(FUNCTION_NAME, "%1.1f %% mismatch between AW and AWR (%s)",
             GetText(nuc + NUCLIDE_PTR_NAME),
             fabs(RDB[nuc + NUCLIDE_AW]/RDB[nuc + NUCLIDE_AWR] - M_NEUTRON));

      /* Lambda and number of decay modes */

      if (RDB[nuc + NUCLIDE_LAMBDA] == 0.0)
        {
          /* Stable nuclide */

          if ((long)RDB[nuc + NUCLIDE_N_DECAY_REA] != 0)
            Warn(FUNCTION_NAME, "Stable nuclide %s has %d decay modes",
                 GetText(nuc + NUCLIDE_PTR_NAME),
                 (long)RDB[nuc + NUCLIDE_N_DECAY_REA]);
        }
      else
        {
          /* Radioactive nuclide */

          if ((long)RDB[nuc + NUCLIDE_N_DECAY_REA] == 0)
            Warn(FUNCTION_NAME, "Radioactive nuclide %s has no decay modes",
                GetText(nuc + NUCLIDE_PTR_NAME));
        }

      /* Nuclide type and transport modes */

      if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY) &&
          (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TRANSMU_DATA)))
        {
          /* Decay nuclide without transmutation data */

          if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA] != 0)
            Die(FUNCTION_NAME, "Decay nuclide %s has %d reaction channels",
                GetText(nuc + NUCLIDE_PTR_NAME),
                (long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA]);
          if ((long)RDB[nuc + NUCLIDE_N_SPECIAL_REA] != 0)
            Die(FUNCTION_NAME,"Decay nuclide %s has %d special cross sections",
                GetText(nuc + NUCLIDE_PTR_NAME),
                (long)RDB[nuc + NUCLIDE_N_SPECIAL_REA]);
        }
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
        {
          /* Dosimetry type nuclide */

          if ((long)(RDB[nuc + NUCLIDE_N_TRANSPORT_REA] +
                     RDB[nuc + NUCLIDE_N_SPECIAL_REA]) == 0)
            Warn(FUNCTION_NAME, "Dosimetry nuclide %s has no reaction data",
                GetText(nuc + NUCLIDE_PTR_NAME));
        }
      else
        {
          /* Transport or SAB type nuclide */

          if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA] == 0)
            Warn(FUNCTION_NAME, "Nuclide %s has no reaction channels",
                GetText(nuc + NUCLIDE_PTR_NAME));
        }

      /* Check transport data */

      if (((long)RDB[DATA_BURN_DECAY_CALC] == NO) &&
          (RDB[DATA_N_TRANSPORT_NUCLIDES] == 0.0) &&
          ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES))
        Error(0, "No neutron data in neutron transport problem");

      if ((RDB[DATA_N_PHOTON_NUCLIDES] == 0.0) &&
          ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES))
        Error(0, "No photon data in photon transport problem");

      /* Check flags */

      if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
           NUCLIDE_FLAG_TRANSPORT_DATA)
          && ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA] == 0))
        Warn(FUNCTION_NAME, "Transport flag without reaction data (%s)",
            GetText(nuc + NUCLIDE_PTR_NAME));

      if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
           NUCLIDE_FLAG_DOSIMETRY_DATA)
          && ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA] +
              (long)RDB[nuc + NUCLIDE_N_SPECIAL_REA] == 0))
        Warn(FUNCTION_NAME, "Dosimetry flag without reaction data (%s)",
            GetText(nuc + NUCLIDE_PTR_NAME));

      if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DECAY_DATA)
          && ((long)RDB[nuc + NUCLIDE_N_DECAY_REA] == 0)
          && (RDB[nuc + NUCLIDE_LAMBDA] > 0.0))
        Warn(FUNCTION_NAME, "Decay flag without reaction data (%s)",
            GetText(nuc + NUCLIDE_PTR_NAME));

      /* Check that nuclide originated from somewhere */

      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
        {
          /* Involved in depletion calculation */

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DP) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_AP) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BP) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SRC) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FP))
            Die(FUNCTION_NAME, "Nuclide %s has no origin %ld",
                GetText(nuc + NUCLIDE_PTR_NAME), nuc);
        }
      else
        {
          /* No depletion */

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SRC))
            Die(FUNCTION_NAME, "Nuclide %s has no origin",
                GetText(nuc + NUCLIDE_PTR_NAME));
        }
    }

  /***************************************************************************/

  /***** Check for duplicates ************************************************/

  /* Sort list by ZAI */

  nuc = (long)RDB[DATA_PTR_NUC0];
  SortList(nuc, NUCLIDE_ZAI, SORT_MODE_ASCEND);

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over remaining with same ZAI */

      ptr = NextItem(nuc);
      while (ptr > VALID_PTR)
        {
          /* Compare ZAI */

          if (RDB[nuc + NUCLIDE_ZAI] != RDB[ptr + NUCLIDE_ZAI])
            break;

          /* Compare library ID and temperature and check type */

          if ((RDB[nuc + NUCLIDE_TEMP] == RDB[ptr + NUCLIDE_TEMP]) &&
              CompareStr(nuc + NUCLIDE_PTR_LIB_ID, ptr + NUCLIDE_PTR_LIB_ID) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA) &&
              !((long)RDB[ptr + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BA) &&
              !((long)RDB[ptr + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BA) &&
              (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) ==
               ((long)RDB[ptr + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)))
            Die(FUNCTION_NAME, "Duplicate found ZAI %ld id = %s T = %1.1fT\n",
                (long)RDB[nuc + NUCLIDE_ZAI], GetText(nuc + NUCLIDE_PTR_LIB_ID),
                RDB[nuc + NUCLIDE_TEMP]);

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next */

      nuc = NextItem(nuc);
    }

  /* Sort list by type */

  nuc = (long)RDB[DATA_PTR_NUC0];
  SortList(nuc, NUCLIDE_TYPE, SORT_MODE_ASCEND);

  /***************************************************************************/

  /***** Check daughter temperatures and id's ********************************/

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Get pointer to target */

          if (((tgt = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR) &&
              (tgt != (long)RDB[DATA_PTR_NUCLIDE_LOST]))
            {
              /* Check library ID */

              if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
                    & NUCLIDE_FLAG_SAB_DATA) &&
                  !((long)RDB[tgt + NUCLIDE_TYPE_FLAGS]
                    & NUCLIDE_FLAG_SAB_DATA))
                if (!CompareStr(nuc + NUCLIDE_PTR_LIB_ID,
                                tgt + NUCLIDE_PTR_LIB_ID))
                  Warn(FUNCTION_NAME, "Mismatch in library ID");

              /* Check TMS flags */

              if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) &&
                  ((long)RDB[tgt + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
                {
                  /* Check temperature */

                  if (RDB[nuc + NUCLIDE_TEMP] < RDB[tgt + NUCLIDE_TEMP])
                    Die(FUNCTION_NAME,
                        "Mismatch in temperature %s --> %s %E %E (1)",
                        GetText(nuc + NUCLIDE_PTR_NAME),
                        GetText(tgt + NUCLIDE_PTR_NAME),
                        RDB[nuc + NUCLIDE_TEMP], RDB[tgt + NUCLIDE_TEMP]);
                }
              else if
                (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) &&
                 !((long)RDB[tgt + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
                {
                  /* Check temperature */

                  if (RDB[nuc + NUCLIDE_TEMP] != RDB[tgt + NUCLIDE_TEMP])
                    Die(FUNCTION_NAME,
                        "Mismatch in temperature %s --> %s %E %E (2)",
                        GetText(nuc + NUCLIDE_PTR_NAME),
                        GetText(tgt + NUCLIDE_PTR_NAME),
                        RDB[nuc + NUCLIDE_TEMP], RDB[tgt + NUCLIDE_TEMP]);
                }
              else
                Die(FUNCTION_NAME, "Error in TMS flags");
            }
          else if ((yld = (long)RDB[rea +  REACTION_PTR_FISSY]) > VALID_PTR)
            {
              /* Get pointer to distribution */

              if ((yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR]) < 0)
                Warn(FUNCTION_NAME, "Pointer error");

              /* Loop over distribution */

              m = 0;
              while ((loc0 = ListPtr(yld, m++)) > VALID_PTR)
                {
                  /* Get pointer to target */

                  if (((tgt = (long)RDB[loc0 + FY_PTR_TGT]) > VALID_PTR) &&
                      (tgt != (long)RDB[DATA_PTR_NUCLIDE_LOST]))
                    {
                      /* Check library ID */

                      if (!CompareStr(nuc + NUCLIDE_PTR_LIB_ID,
                                      tgt + NUCLIDE_PTR_LIB_ID))
                        Warn(FUNCTION_NAME, "Mismatch in library ID");

              /* Check TMS flags */

              if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) &&
                  ((long)RDB[tgt + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
                {
                  /* Check temperature */

                  if (RDB[nuc + NUCLIDE_TEMP] < RDB[tgt + NUCLIDE_TEMP])
                    Die(FUNCTION_NAME,
                        "Mismatch in temperature %s --> %s %E %E",
                        GetText(nuc + NUCLIDE_PTR_NAME),
                        GetText(tgt + NUCLIDE_PTR_NAME),
                        RDB[nuc + NUCLIDE_TEMP], RDB[tgt + NUCLIDE_TEMP]);
                }
              else if
                (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) &&
                 !((long)RDB[tgt + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
                {
                  /* Check temperature */

                  if (RDB[nuc + NUCLIDE_TEMP] != RDB[tgt + NUCLIDE_TEMP])
                    Die(FUNCTION_NAME,
                        "Mismatch in temperature %s --> %s %E %E",
                        GetText(nuc + NUCLIDE_PTR_NAME),
                        GetText(tgt + NUCLIDE_PTR_NAME),
                        RDB[nuc + NUCLIDE_TEMP], RDB[tgt + NUCLIDE_TEMP]);
                }
              else
                Die(FUNCTION_NAME, "Error in TMS flags");
                    }
                }
            }

          /* Next */

          rea = NextItem(rea);
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***************************************************************************/

  /* Reset table counter */

  table = 0;

  /* Get output file */

  sprintf(outfile, "%s.out", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for writing */

  if ((fp = fopen(outfile, "a")) == NULL)
    Warn(FUNCTION_NAME, "Unable to open file for writing");

  /***************************************************************************/

  /***** Print summary data **************************************************/

  table++;

  if (1 != 2)
    {
      fprintf(fp, "\n --- Table %2ld: Summary of nuclide data: \n\n", table);
      fprintf(fp, " Data for %ld nuclides included in calculation:\n\n",
              (long)RDB[DATA_N_TOT_NUCLIDES]);

      if (RDB[DATA_N_TRANSPORT_NUCLIDES] > 0.0)
        fprintf(fp, " - %ld Transport nuclides\n",
                (long)RDB[DATA_N_TRANSPORT_NUCLIDES]);
      else
        fprintf(fp, " - No transport nuclides\n");

      if (RDB[DATA_N_DOSIMETRY_NUCLIDES] > 0.0)
        fprintf(fp, " - %ld Dosimetry nuclides\n",
                (long)RDB[DATA_N_DOSIMETRY_NUCLIDES]);
      else
        fprintf(fp, " - No dosimetry nuclides\n");

      if (RDB[DATA_N_DECAY_NUCLIDES] > 0.0)
        fprintf(fp, " - %ld Decay nuclides\n",
                (long)RDB[DATA_N_DECAY_NUCLIDES]);
      else
        fprintf(fp, " - No decay nuclides\n");

      if (RDB[DATA_N_TRANSPORT_REA] > 0.0)
        fprintf(fp, " - %ld reaction channels\n",
                (long)RDB[DATA_N_TRANSPORT_REA]);
      else
        fprintf(fp, " - No reaction channels\n");

      if (RDB[DATA_N_SPECIAL_REA] > 0.0)
        fprintf(fp, " - %ld special reactions\n",
                (long)RDB[DATA_N_SPECIAL_REA]);
      else
        fprintf(fp, " - No special reactions\n");

      if (RDB[DATA_N_TRANSMUTATION_REA] > 0.0)
        fprintf(fp, " - %ld transmutation reactions\n",
                (long)RDB[DATA_N_TRANSMUTATION_REA]);
      else
        fprintf(fp, " - No transmutation reactions\n");

      if (RDB[DATA_N_TRANSPORT_BRANCH] > 0.0)
        fprintf(fp, " - %ld transmutation branch reactions\n",
                (long)RDB[DATA_N_TRANSPORT_BRANCH]);
      else
        fprintf(fp, " - No transmutation branch reactions\n");

      if (RDB[DATA_N_DECAY_REA] > 0.0)
        fprintf(fp, " - %ld decay reactions\n",
                (long)RDB[DATA_N_DECAY_REA]);
      else
        fprintf(fp, " - No decay reactions\n");

      if (RDB[DATA_N_DECAY_BRANCH] > 0.0)
        fprintf(fp, " - %ld decay branch reactions\n",
                (long)RDB[DATA_N_DECAY_BRANCH]);
      else
        fprintf(fp, " - No decay branch reactions\n");

      if (RDB[DATA_N_DEAD_PATH] > 0.0)
        fprintf(fp, " - %ld lost transmutation paths\n\n",
                (long)RDB[DATA_N_DEAD_PATH]);
      else
        fprintf(fp, " - All reaction linked to target\n\n");

      fprintf(fp, " Columns:\n\n");
      fprintf(fp, "  1. Index\n");
      fprintf(fp, "  2. Level (first occurence in transmutation paths)\n");
      fprintf(fp, "  3. Name (nuclide identifier data)\n");
      fprintf(fp, "  4. Library ID\n");
      fprintf(fp, "  5. Primary data type: TRA = transport\n");
      fprintf(fp, "                        DOS = dosimetry\n");
      fprintf(fp, "                        SAB = thermal scattering\n");
      fprintf(fp, "                        DEC = only decay data\n");
      fprintf(fp, "                        STR = only structural data (mass)\n");
      fprintf(fp, "                        PHO = photon interaction data\n");
      fprintf(fp, "  6. Z (atomic number)\n");
      fprintf(fp, "  7. A (isotope number)\n");
      fprintf(fp, "  8. I (isomeric state)\n");
      fprintf(fp, "  9. Atomic weight\n");
      fprintf(fp, " 10. Temperature (in Kelvin)\n");
      fprintf(fp, " 11. Decay constant (1/s)\n");
      fprintf(fp, " 12. Flags: I - nuclide is present in initial composition\n");
      fprintf(fp, "            B - nuclide is involved in burnup calculation\n");
      fprintf(fp, "            A - nuclide is produced as activation product\n");
      fprintf(fp, "            D - nuclide is produced as decay product\n");
      fprintf(fp, "            F - nuclide is produced as fission product\n");
      fprintf(fp, "            B - nuclide is produced in decay or transmutation branching\n");
      fprintf(fp, "            T - nuclide has neutron transport cross section data\n");
      fprintf(fp, "            S - nuclide has S(a,b) data\n");
      fprintf(fp, "            I - nuclide has isomeric branching reactions\n");
      fprintf(fp, "            F - nuclide is fissile\n");
      fprintf(fp, "            P - nuclide is a delayed neutron precursor\n");
      fprintf(fp, "            D - nuclide has delayed neutron data\n");
      fprintf(fp, "            U - nuclide has uses ures ptable data\n");
      fprintf(fp, "            T - nuclide cross sections are adjusted by TMS\n");
      fprintf(fp, "            D - nuclide has dosimetry cross section data\n");
      fprintf(fp, "            R - nuclide has radioactive decay data\n");
      fprintf(fp, "            N - nuclide has neutron-induced fission yield data\n");
      fprintf(fp, "            S - nuclide has spontaneous fission yield data\n");
      fprintf(fp, "            P - nuclide has photon transport cross section data\n");
      fprintf(fp, "            T - nuclide has supplementary transmutation cross section data\n\n");

      /* Loop over nuclides */

      n = 0;
      while ((nuc = ListPtr((long)RDB[DATA_PTR_NUC0], n++)) > VALID_PTR)
        {
                /* Counter */

          fprintf(fp, " %4ld  ", n);

          /* Level */

          fprintf(fp, "%2ld  ", (long)RDB[nuc + NUCLIDE_PATH_LEVEL]);

          /* Name and library id */

          fprintf(fp, "%-12s  %3s  ",
                  GetText(nuc + NUCLIDE_PTR_NAME),
                  GetText(nuc + NUCLIDE_PTR_LIB_ID));


          /* Primary type */

          if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
            fprintf(fp, "TRA  ");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
            fprintf(fp, "DOS  ");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB)
            fprintf(fp, "SAB  ");
          else if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY) &&
                   ((long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE] > VALID_PTR))
            fprintf(fp, "DEC  ");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
            fprintf(fp, "STR  ");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
            fprintf(fp, "PHO  ");
          else
            Warn(FUNCTION_NAME, "Invalid type %ld for %s",
                (long)RDB[nuc + NUCLIDE_TYPE],
                GetText(nuc + NUCLIDE_PTR_NAME));

          /* Z, A and I */

          fprintf(fp, "%3ld %3ld %ld  ", (long)RDB[nuc + NUCLIDE_Z],
                  (long)RDB[nuc + NUCLIDE_A],
                  (long)RDB[nuc + NUCLIDE_I]);

          /* AW */

          fprintf(fp, "%9.5f  ", RDB[nuc + NUCLIDE_AW]);

          /* Temperature */

          fprintf(fp, "%6.1f  ", RDB[nuc + NUCLIDE_TEMP]);

          /* Decay constant */

          if ((ace = (long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE]) < VALID_PTR)
            fprintf(fp, "        N/A  ");
          else if (ACE[ace + ACE_LAMBDA] == 0.0)
            fprintf(fp, "     stable  ");
          else
            fprintf(fp, "%11.5E  ", RDB[nuc + NUCLIDE_LAMBDA]);

          /* Flags */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL)
            fprintf(fp, "I");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
            fprintf(fp, "B");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_AP)
            fprintf(fp, "A");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DP)
            fprintf(fp, "D");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FP)
            fprintf(fp, "F");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BP)
            fprintf(fp, "B");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_TRANSPORT_DATA)
            fprintf(fp, "T");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_SAB_DATA)
            fprintf(fp, "S");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_BRA_DATA)
            fprintf(fp, "I");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_FISSILE)
            fprintf(fp, "F");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_DELNU_PREC)
            fprintf(fp, "P");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_ACE_PREC_GROUPS] > 0)
            fprintf(fp, "D");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_URES_USED)
            fprintf(fp, "U");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_TMS)
            fprintf(fp, "T");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_DOSIMETRY_DATA)
            fprintf(fp, "D");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_DECAY_DATA)
            fprintf(fp, "R");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_NFY_DATA)
            fprintf(fp, "N");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SFY_DATA)
            fprintf(fp, "S");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_PHOTON_DATA)
            fprintf(fp, "P");
          else
            fprintf(fp, "-");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TRANSMU_DATA)
            fprintf(fp, "T");
          else
            fprintf(fp, "-");

          /* Name and symbol */

          sprintf(sym, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 1));
          sprintf(name, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 2));

          fprintf(fp, "  %s (%s)", name, sym);

          /* Newline */

          fprintf(fp, "\n");
        }
    }

  /***************************************************************************/

  /***** Print reaction and decay data ***************************************/

  table++;

  if (1 != 2)
    {
      fprintf(fp, "\n --- Table %2ld: Reaction and decay data: \n",  table);

      /* Loop over nuclides */

      n = 0;
      while ((nuc = ListPtr((long)RDB[DATA_PTR_NUC0], n++)) > VALID_PTR)
        {
          /* Name and symbol */

          sprintf(sym, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 1));
          sprintf(name, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 2));

          /* Print header data */

          fprintf(fp, "\n Nuclide %4ld / %ld : %s -- %s (%s)\n\n", n,
                  (long)RDB[DATA_N_TOT_NUCLIDES],
                  GetText(nuc + NUCLIDE_PTR_NAME), name, sym);

          fprintf(fp, " Pointers                     : %ld %ld\n",
                  nuc, (long)RDB[nuc + NUCLIDE_PTR_ACE]);

          fprintf(fp, " Primary type                 : ");

          if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
            fprintf(fp, "Transport\n");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
            fprintf(fp, "Dosimetry\n");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB)
            fprintf(fp, "Thermal scattering\n");
          else if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY) &&
                   ((long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE] > VALID_PTR))
            fprintf(fp, "Decay\n");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
            fprintf(fp, "Decay (no data)\n");
          else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
            fprintf(fp, "Photon interaction\n");
          else
            Warn(FUNCTION_NAME, "Invalid type %ld for %s",
                (long)RDB[nuc + NUCLIDE_TYPE],
                GetText(nuc + NUCLIDE_PTR_NAME));

          fprintf(fp, " Nuclide ZAI                  : %ld\n",
                  (long)RDB[nuc + NUCLIDE_ZAI]);

          fprintf(fp, " Library ID                   : %s\n",
                  GetText(nuc + NUCLIDE_PTR_LIB_ID));

          if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_ACE]) > 0)
            WDB[DATA_DUMMY] = ACE[ptr + ACE_PTR_FILE];
          else
            WDB[DATA_DUMMY] = -1.0;

          if ((long)RDB[DATA_DUMMY] > VALID_PTR)
            fprintf(fp, " Cross section data file      : %s\n",
                    GetText(DATA_DUMMY));
          else
            fprintf(fp, " Cross section data file      : N/A\n");

          fprintf(fp, " Atomic weight (AW)           : %1.5f\n",
                  RDB[nuc + NUCLIDE_AW]);
          fprintf(fp, " Atomic weight ratio (AWR)    : %1.5f\n",
                  RDB[nuc + NUCLIDE_AWR]);
          fprintf(fp, " Temperature                  : %1.1f Kelvin\n",
                  RDB[nuc + NUCLIDE_TEMP]);
          fprintf(fp, " XS temperature               : %1.1f Kelvin\n",
                  RDB[nuc + NUCLIDE_XS_TEMP]);

          fprintf(fp, " Minimum and maximum energy   : %1.5E %1.5E\n",
                  RDB[nuc + NUCLIDE_EMIN], RDB[nuc + NUCLIDE_EMAX]);

          fprintf(fp, " Half-life                    :");

          if ((ace = (long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE]) < VALID_PTR)
            fprintf(fp, " N/A\n");
          else if (ACE[ace + ACE_LAMBDA] == 0.0)
            fprintf(fp, " stable\n");
          else if (RDB[nuc + NUCLIDE_LAMBDA] == 0.0)
            fprintf(fp, " stable beyond cutoff (Half-life > %s)\n",
                    TimeIntervalStr(RDB[DATA_DEP_HALF_LIFE_CUTOFF]));
          else
            {
              /* Convert to half-life */

              t = log(2.0)/RDB[nuc + NUCLIDE_LAMBDA];

              /* Print */

              if (t < 60.0)
                fprintf(fp, " %s \n", TimeIntervalStr(t));
              else
                fprintf(fp, " %1.2E seconds (%s)\n", t,
                        TimeIntervalStr(t));
            }

          if (RDB[nuc + NUCLIDE_SPEC_ING_TOX] > 0.0)
            fprintf(fp, " Specific ingestion toxicity  : %1.2E Sv/Bq\n",
                    RDB[nuc + NUCLIDE_SPEC_ING_TOX]);
          else
            fprintf(fp, " Specific ingestion toxicity  : N/A\n");

          if (RDB[nuc + NUCLIDE_SPEC_INH_TOX] > 0.0)
            fprintf(fp, " Specific inhalation toxicity : %1.2E Sv/Bq\n",
                    RDB[nuc + NUCLIDE_SPEC_INH_TOX]);
          else
            fprintf(fp, " Specific inhalation toxicity : N/A\n");

          /* Energy deposited per fission */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
            fprintf(fp, " Energy deposited per fission : %1.2f MeV\n",
                    RDB[nuc + NUCLIDE_FISSE]/MEV);

          /* Nubar for spontaneous fission */

          if (RDB[nuc + NUCLIDE_SF_NUBAR] > 0.0)
            fprintf(fp, " Spontaneous fission nubar    : %1.5f\n",
                    RDB[nuc + NUCLIDE_SF_NUBAR]);

          /* Newline */

          fprintf(fp, "\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
            {
              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL)
                fprintf(fp, " - Nuclide is present in initial composition\n");
              else
                fprintf(fp, " - Nuclide is not present in initial composition\n");

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_AP)
                fprintf(fp, " - Nuclide is produced in transmutation\n");
              else
                fprintf(fp, " - Nuclide is not produced in transmutation\n");

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DP)
                fprintf(fp, " - Nuclide is produced in decay\n");
              else
                fprintf(fp, " - Nuclide is not produced in decay\n");

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FP)
                fprintf(fp, " - Nuclide is produced in fission\n");
              else
                fprintf(fp, " - Nuclide is not produced in fission\n");
            }
          else
            fprintf(fp, " - Nuclide is not involved in burnup calculation\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_TRANSPORT_DATA)
            fprintf(fp, " - Nuclide has neutron transport data\n");
          else
            fprintf(fp, " - Nuclide has no neutron transport data\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_PHOTON_DATA)
            fprintf(fp, " - Nuclide has photon transport data\n");
          else
            fprintf(fp, " - Nuclide has no photon transport data\n");
          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA)
            {
              /* Loop over S(a,b) reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  if ((long)RDB[rea + REACTION_MT] == 1002)
                    fprintf(fp, " - Nuclide has elastic S(a,b) channel between %1.2E and %1.2E MeV\n", RDB[rea + REACTION_EMIN], RDB[rea + REACTION_EMAX]);

                  else if ((long)RDB[rea + REACTION_MT] == 1004)
                    fprintf(fp, " - Nuclide has inelastic S(a,b) channel between %1.2E and %1.2E MeV\n", RDB[rea + REACTION_EMIN], RDB[rea + REACTION_EMAX]);

                  /* Next */

                  rea = NextItem(rea);
                }
            }
          else
            fprintf(fp, " - Nuclide has no S(a,b) data\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BRA_DATA)
            fprintf(fp, " - Nuclide has isomeric branching data\n" );
          else
            fprintf(fp, " - Nuclide has no isomeric branching data\n" );

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
            fprintf(fp, " - Nuclide is fissile\n" );
          else
            fprintf(fp, " - Nuclide is not fissile\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DELNU_PREC)
            fprintf(fp, " - Nuclide is a delayed neutron precursor\n" );
          else
            fprintf(fp, " - Nuclide is not a delayed neutron precursor\n");

          if ((long)RDB[nuc + NUCLIDE_ACE_PREC_GROUPS] > 0)
            fprintf(fp,
                    " - Nuclide has delayed neutron data in %ld precursor groups\n", (long)RDB[nuc + NUCLIDE_ACE_PREC_GROUPS]);
          else
            fprintf(fp, " - Nuclide has no delayed neutron data\n");

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                NUCLIDE_FLAG_URES_AVAIL))
            fprintf(fp, " - Nuclide has no ures ptable data\n");
          else if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                NUCLIDE_FLAG_URES_USED))
            fprintf(fp, " - Nuclide has ures ptable data but it is not used\n");
          else if ((RDB[nuc + NUCLIDE_URES_EMIN] < 0.0) && (RDB[nuc + NUCLIDE_URES_EMAX] < 0.0))
            fprintf(fp, " - Nuclide has unsupported ures ptable data that is not used\n");
          else
            fprintf(fp, " - Nuclide uses ures ptable data between %1.5E and %1.5E MeV\n", RDB[nuc + NUCLIDE_URES_EMIN], RDB[nuc + NUCLIDE_URES_EMAX]);

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
            fprintf(fp, " - Nuclide temperature is adjusted with TMS between %1.1f and %1.1f Kelvin\n", RDB[nuc + NUCLIDE_TMS_MIN_TEMP],
                    RDB[nuc + NUCLIDE_TMS_MAX_TEMP]);
          else
            fprintf(fp, " - Nuclide temperature is not adjusted with TMS\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_DOSIMETRY_DATA)
            fprintf(fp, " - Nuclide has dosimetry data\n");
          else
            fprintf(fp, " - Nuclide has no dosimetry data\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_DECAY_DATA)
            fprintf(fp, " - Nuclide has decay data\n");
          else
            fprintf(fp, " - Nuclide has no decay data\n");

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS]
              & NUCLIDE_FLAG_TRANSMU_DATA)
            fprintf(fp, " - Nuclide has supplementary transmutation cross section data\n");
          else
            fprintf(fp, " - Nuclide has no supplementary transmutation cross section data\n");

           if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA]) < VALID_PTR)
            fprintf(fp, " - Nuclide has no NFY data\n");
          else
            {
              if ((long)RDB[ptr + FISSION_YIELD_PARENT_ZAI] !=
                   (long)RDB[nuc + NUCLIDE_ZAI])
                fprintf(fp, " - Nuclide is using NFY data from %s\n",
                        ZAItoIso((long)RDB[ptr + FISSION_YIELD_PARENT_ZAI],
                                 1));
              else
                fprintf(fp, " - Nuclide has NFY data \n");
            }

          if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_SFY_DATA]) < VALID_PTR)
            fprintf(fp, " - Nuclide has no SFY data\n\n");
          else
            fprintf(fp, " - Nuclide has SFY data\n\n");

          /* Modes etc. */

          if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA] > 0)
            fprintf(fp,  " - %ld reaction channels\n",
                    (long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA]);
          else
            fprintf(fp,  " - No reaction channels\n");

          if ((long)RDB[nuc + NUCLIDE_N_SPECIAL_REA] > 0)
            fprintf(fp,  " - %ld special reactions\n",
                    (long)RDB[nuc + NUCLIDE_N_SPECIAL_REA]);
          else
            fprintf(fp,  " - No special reactions\n");

          if ((long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA] > 0)
            fprintf(fp,  " - %ld transmutation reactions\n",
                    (long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA]);
          else
            fprintf(fp,  " - No transmutation reactions\n");

          if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] > 0)
            fprintf(fp,  " - %ld transmutation branch reactions\n",
                    (long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH]);
          else
            fprintf(fp,  " - No transmutation branch reactions\n");

          if ((long)RDB[nuc + NUCLIDE_N_DECAY_REA] > 0)
            fprintf(fp,  " - %ld decay reactions\n",
                    (long)RDB[nuc + NUCLIDE_N_DECAY_REA]);
          else
            fprintf(fp,  " - No decay reactions\n");

          if ((long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH] > 0)
            fprintf(fp,  " - %ld decay branch reactions\n",
                    (long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH]);
          else
            fprintf(fp,  " - No decay branch reactions\n");

          if ((long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_PATH] > 0)
            fprintf(fp,  " - %ld transmutation paths\n",
                    (long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_PATH]);
          else
            fprintf(fp,  " - No transmutation paths\n");

          if ((long)RDB[nuc + NUCLIDE_N_DEAD_PATH] > 0)
            fprintf(fp,  " - %ld lost transmutation paths\n",
                    (long)RDB[nuc + NUCLIDE_N_DEAD_PATH]);
          else if ((long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA] +
                   (long)RDB[nuc + NUCLIDE_N_DECAY_REA] > 0)
            fprintf(fp,  " - All reaction paths linked to target\n");

          /* Print transport reactions */

          if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA] > 0)
            {
              if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA]  == 1)
                fprintf(fp, "\n 1 reaction channel:\n\n");
              else
                fprintf(fp, "\n %ld reactions channel:\n\n",
                        (long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA]);

              /* Reset count */

              i = 0;
              nsab = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

              while (rea > VALID_PTR)
                {
                  /* Get mt */

                  mt = (long)RDB[rea + REACTION_MT];

                  /* add to number of S(a,b) reactions */

                  if ((mt == 1002) || (mt == 1004) ||
                      (mt == 2002) || (mt == 2004))
                    nsab++;

                  /* Check type */

                  if ((long)RDB[rea + REACTION_TYPE] ==
                      REACTION_TYPE_PARTIAL)
                    {
                      /* Counter */

                      fprintf(fp, " %3ld  ", 1 + i++);

                      /* MT */

                      fprintf(fp, "  MT = %6ld  ", mt);

                      /* Q-value */

                      fprintf(fp, "Q = %12.5E MeV  ",
                              RDB[rea + REACTION_Q]);

                      /* Minimum energy */

                      fprintf(fp, "Emin = %11.5E MeV  ",
                              RDB[rea + REACTION_EMIN]);

                      /* Fraction */

                      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)
                        fprintf(fp, "frac = %11.5E  ",
                                RDB[rea + REACTION_BR]);

                      /* Target */

                      if ((ptr = (long)RDB[rea + REACTION_PTR_FISSY])
                          > VALID_PTR)
                        fprintf(fp, "NFY (E = %11.5E MeV)     ",
                                RDB[ptr + FISSION_YIELD_E]);
                      else if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
                        {
                          /* Get target ZAI */

                          ZAI = ReactionTargetZAI(rea);

                          if ((ptr = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
                            {
                              if (ptr != RDB[DATA_PTR_NUCLIDE_LOST])
                                fprintf(fp, "Product nuclide = %-10s    ", GetText(ptr + NUCLIDE_PTR_NAME));
                              else if (ZAI > 0)
                                fprintf(fp, "Product nuclide = lost: %-7s ", ZAItoIso(ZAI, 1));
                              else
                                Die(FUNCTION_NAME, "Horrible error 1");
                            }
                          else if (ZAI > 0)
                            fprintf(fp, "Product nuclide = cut: %-8s ", ZAItoIso(ZAI, 1));
                          else if (ZAI == 0)
                            fprintf(fp, "Product nuclide = self          ");
                          else if ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] < 0)
                            fprintf(fp, "Product nuclide = FP            ");
                          else
                            Die(FUNCTION_NAME, "Horrible error 2 (%ld)", ZAI);

                        }

                      /* Description */

                      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
                        fprintf(fp, ": %s", ReactionMT(mt, YES));
                      else
                        fprintf(fp, ": %s", ReactionMT(mt, NO));

                      /* Newline */

                      fprintf(fp, "\n");
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Check S(a,b) flag and reaction count */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA)
                {
                  if (nsab == 0)
                    Warn(FUNCTION_NAME, "Nuclide %s has S(a,b) flag but no data",
                        GetText(nuc + NUCLIDE_PTR_NAME));
                }
              else
                {
                  if (nsab > 0)
                    Warn(FUNCTION_NAME, "Nuclide %s has S(a,b) data but no flat",
                        GetText(nuc + NUCLIDE_PTR_NAME));
                }

              /* Check count */

              if (i != (long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA])
                Warn(FUNCTION_NAME, "Mismatch in count (%s) 1",
                    GetText(nuc + NUCLIDE_PTR_NAME));
            }

          /* Print additional transport branches */

          if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] > 0)
            {
              if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] == 1)
                fprintf(fp, "\n 1 additional transport branch:\n\n");
              else
                fprintf(fp, "\n %ld additional transport branches:\n\n",
                        (long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH]);

              /* Reset count */

              i = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

              while (rea > VALID_PTR)
                {
                  /* Get mt */

                  mt = (long)RDB[rea + REACTION_BRANCH_MT];

                  /* Check type */

                  if ((long)RDB[rea + REACTION_TYPE] ==
                      REACTION_TYPE_TRA_BRANCH)
                    {
                      /* Counter */

                      fprintf(fp, " %3ld  ", 1 + i++);

                      /* MT */

                      fprintf(fp, "  MT = %6ld  ", mt);

                      /* Q-value */

                      fprintf(fp, "Q = %12.5E MeV  ",
                              RDB[rea + REACTION_Q]);

                      /* Minimum energy */

                      fprintf(fp, "Emin = %11.5E MeV  ",
                              RDB[rea + REACTION_EMIN]);

                      /* Check if fission distribution */

                      /* Fraction */

                      fprintf(fp, "frac = %11.5E  ",
                              RDB[rea + REACTION_BR]);
                      if ((ptr = (long)RDB[rea + REACTION_PTR_FISSY])
                          > VALID_PTR)
                        fprintf(fp, "NFY (E = %11.5E MeV)     ",
                                RDB[ptr + FISSION_YIELD_E]);
                      else if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
                        {
                          /* Get target ZAI */

                          ZAI = ReactionTargetZAI(rea);

                          if ((ptr = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
                            {
                              if (ptr != RDB[DATA_PTR_NUCLIDE_LOST])
                                fprintf(fp, "Product nuclide = %-10s    ", GetText(ptr + NUCLIDE_PTR_NAME));
                              else if (ZAI > 0)
                                fprintf(fp, "Product nuclide = lost: %-7s ", ZAItoIso(ZAI, 1));
                              else
                                Die(FUNCTION_NAME, "Horrible error 3");
                            }
                          else if (ZAI > 0)
                            fprintf(fp, "Product nuclide = cut: %-8s ", ZAItoIso(ZAI, 1));
                          else if (ZAI == 0)
                            fprintf(fp, "Product nuclide = self          ");
                          else
                            Die(FUNCTION_NAME, "Horrible error 4");
                        }

                      /* Description */

                      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
                        fprintf(fp, ": %s", ReactionMT(mt, YES));
                      else
                        fprintf(fp, ": %s", ReactionMT(mt, NO));

                      /* Newline */

                      fprintf(fp, "\n");
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Check count */

              if (i != (long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH])
                Warn(FUNCTION_NAME, "Mismatch in count (%s) 2",
                    GetText(nuc + NUCLIDE_PTR_NAME));
            }

          /* Print special reactions */

          if ((long)RDB[nuc + NUCLIDE_N_SPECIAL_REA] > 0)
            {
              if ((long)RDB[nuc + NUCLIDE_N_SPECIAL_REA] == 1)
                fprintf(fp, "\n 1 special reaction:\n\n");
              else
                fprintf(fp, "\n %ld special reactions:\n\n",
                        (long)RDB[nuc + NUCLIDE_N_SPECIAL_REA]);

              /* Reset count */

              i = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

              while (rea > VALID_PTR)
                {
                  /* Get mt */

                  mt = (long)RDB[rea + REACTION_MT];

                  /* Check type */

                  if ((long)RDB[rea + REACTION_TYPE] ==
                      REACTION_TYPE_SPECIAL)
                    {
                      /* Counter */

                      fprintf(fp, " %3ld  ", 1 + i++);

                      /* Check threshold cut-off */

                      /* MT */

                      fprintf(fp, "  MT = %6ld  ", mt);

                      /* Q-value */

                      fprintf(fp, "Q = %12.5E MeV  ",
                              RDB[rea + REACTION_Q]);

                      /* Minimum energy */

                      fprintf(fp, "Emin = %11.5E MeV  ",
                              RDB[rea + REACTION_EMIN]);

                      /* Description */

                      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
                        fprintf(fp, ": %s", ReactionMT(mt, YES));
                      else
                        fprintf(fp, ": %s", ReactionMT(mt, NO));

                      /* Newline */

                      fprintf(fp, "\n");
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Check count */

              if (i != (long)RDB[nuc + NUCLIDE_N_SPECIAL_REA])
                Warn(FUNCTION_NAME, "Mismatch in count (%s) 3",
                    GetText(nuc + NUCLIDE_PTR_NAME));
            }


          /* Print decay reactions */

          if ((long)RDB[nuc + NUCLIDE_N_DECAY_REA] > 0)
            {
              if ((long)RDB[nuc + NUCLIDE_N_DECAY_REA] == 1)
                fprintf(fp, "\n 1 decay reaction:\n\n");
              else
                fprintf(fp, "\n %ld decay reactions:\n\n",
                        (long)RDB[nuc + NUCLIDE_N_DECAY_REA]);

              /* Reset count */

              i = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

              while (rea > VALID_PTR)
                {
                  /* Get mt */

                  mt = (long)RDB[rea + REACTION_MT];

                  /* Check type */

                  if ((long)RDB[rea + REACTION_TYPE] ==
                      REACTION_TYPE_DECAY)
                    {
                      /* Counter */

                      fprintf(fp, " %3ld  ", 1 + i++);

                      /* MT */

                      fprintf(fp, "RTYP = %ld", mt - 10000);

                      if ((long)RDB[rea + REACTION_RTYP2] > 10000)
                        fprintf(fp, "%ld",
                                (long)RDB[rea + REACTION_RTYP2] - 10000);
                      else
                        fprintf(fp, " ");
                      if ((long)RDB[rea + REACTION_RTYP3] > 10000)
                        fprintf(fp, "%ld",
                                (long)RDB[rea + REACTION_RTYP3] - 10000);
                      else
                        fprintf(fp, " ");
                      if ((long)RDB[rea + REACTION_RTYP4] > 10000)
                        fprintf(fp, "%ld",
                                (long)RDB[rea + REACTION_RTYP4] - 10000);
                      else
                        fprintf(fp, " ");
                      if ((long)RDB[rea + REACTION_RTYP5] > 10000)
                        fprintf(fp, "%ld",
                                (long)RDB[rea + REACTION_RTYP5] - 10000);
                      else
                        fprintf(fp, " ");

                      /* Q-value */

                      fprintf(fp, " Q = %12.5E MeV  ",
                              RDB[rea + REACTION_Q]);

                      /* Minimum energy */
                      /*
                      fprintf(fp, "Emin = %11.5E MeV  ",
                              RDB[rea + REACTION_EMIN]);
                      */

                      /* Fraction */

                      fprintf(fp, "frac = %11.5E  ",
                              RDB[rea + REACTION_BR]);

                      /* Target */

                      if ((ptr = (long)RDB[rea + REACTION_PTR_FISSY])
                          > VALID_PTR)
                        fprintf(fp, "SFY (E = %11.5E MeV)       ",
                                RDB[ptr + FISSION_YIELD_E]);
                      else if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
                        {
                          /* Get target ZAI */

                          ZAI = ReactionTargetZAI(rea);

                          if ((ptr = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
                            {
                              if (ptr != RDB[DATA_PTR_NUCLIDE_LOST])
                                fprintf(fp, "Product nuclide = %-10s    ", GetText(ptr + NUCLIDE_PTR_NAME));
                              else if (ZAI > 0)
                                fprintf(fp, "Product nuclide = lost: %-7s ", ZAItoIso(ZAI, 1));
                              else if (ZAI == -2)
                                fprintf(fp, "SF without yield                ");
                              else
                                Die(FUNCTION_NAME, "Horrible error 5");
                            }
                          else if (ZAI > 0)
                            fprintf(fp, "Product nuclide = cut: %-8s ", ZAItoIso(ZAI, 1));
                          else if (ZAI == 0)
                            fprintf(fp, "Product nuclide = self          ");
                          else if (ZAI != FISSION_YIELD_TYPE_SFY)
                            Die(FUNCTION_NAME, "Horrible error 6 %s %ld %ld", GetText(nuc + NUCLIDE_PTR_NAME), mt, ZAI);
                          else if (ZAI == -2)
                            fprintf(fp, "SF without yield                ");
                          else
                            Die(FUNCTION_NAME, "Horrible error 6b %s %ld %ld", GetText(nuc + NUCLIDE_PTR_NAME), mt, ZAI);
                        }

                      /* Description */

                      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
                        fprintf(fp, ": %s", ReactionMT(mt, YES));
                      else
                        fprintf(fp, ": %s", ReactionMT(mt, NO));

                      /* Check additional modes */

                      if ((long)RDB[rea + REACTION_RTYP2] != 10000)
                        fprintf(fp, " + %s",
                                ReactionMT((long)RDB[rea + REACTION_RTYP2], NO));
                      if ((long)RDB[rea + REACTION_RTYP3] != 10000)
                        fprintf(fp, " + %s",
                                ReactionMT((long)RDB[rea + REACTION_RTYP3], NO));
                      if ((long)RDB[rea + REACTION_RTYP4] != 10000)
                        fprintf(fp, " + %s",
                                ReactionMT((long)RDB[rea + REACTION_RTYP4], NO));
                      if ((long)RDB[rea + REACTION_RTYP5] != 10000)
                        fprintf(fp, " + %s",
                                ReactionMT((long)RDB[rea + REACTION_RTYP5], NO));

                      /* Newline */

                      fprintf(fp, "\n");
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Check count */

              if (i != (long)RDB[nuc + NUCLIDE_N_DECAY_REA])
                Warn(FUNCTION_NAME, "Mismatch in count (%s) 4",
                    GetText(nuc + NUCLIDE_PTR_NAME));
            }

          /* Print additional decay branches */

          if ((long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH] > 0)
            {
              if ((long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH] == 1)
                fprintf(fp, "\n 1 additional decay branch:\n\n");
              else
                fprintf(fp, "\n %ld additional decay branches:\n\n",
                        (long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH]);

              /* Reset count */

              i = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

              while (rea > VALID_PTR)
                {
                  /* Get mt */

                  mt = (long)RDB[rea + REACTION_BRANCH_MT];

                  /* Check type */

                  if ((long)RDB[rea + REACTION_TYPE] ==
                      REACTION_TYPE_DEC_BRANCH)
                    {
                      /* Counter */

                      fprintf(fp, " %3ld  ", 1 + i++);

                      /* MT */

                      fprintf(fp, "RTYP = %-6ld  ", mt - 10000);

                      /* Q-value */

                      fprintf(fp, "Q = %12.5E MeV  ",
                              RDB[rea + REACTION_Q]);

                      /* Minimum energy */
                      /*
                      fprintf(fp, "Emin = %11.5E MeV  ",
                              RDB[rea + REACTION_EMIN]);
                      */


                      /* Fraction */

                      fprintf(fp, "frac = %11.5E  ",
                              RDB[rea + REACTION_BR]);

                      /* Check if fission distribution */

                      if ((ptr = (long)RDB[rea + REACTION_PTR_FISSY])
                          > VALID_PTR)
                        fprintf(fp, "SFY (E = %11.5E MeV)       ",
                                RDB[ptr + FISSION_YIELD_E]);
                      else if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
                        {
                          /* Get target ZAI */

                          ZAI = ReactionTargetZAI(rea);

                          if ((ptr = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
                            {
                              if (ptr != RDB[DATA_PTR_NUCLIDE_LOST])
                                fprintf(fp, "Product nuclide = %-10s    ", GetText(ptr + NUCLIDE_PTR_NAME));
                              else if (ZAI > 0)
                                fprintf(fp, "Product nuclide = lost: %-7s ", ZAItoIso(ZAI, 1));
                              else
                                Die(FUNCTION_NAME, "Horrible error 7");
                            }
                          else if (ZAI > 0)
                            fprintf(fp, "Product nuclide = cut: %-8s ", ZAItoIso(ZAI, 1));
                          else if (ZAI == 0)
                            fprintf(fp, "Product nuclide = self          ");
                          else
                            Die(FUNCTION_NAME, "Horrible error 8");
                        }

                      /* Description */

                      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
                        fprintf(fp, ": %s", ReactionMT(mt, YES));
                      else
                        fprintf(fp, ": %s", ReactionMT(mt, NO));

                      /* Newline */

                      fprintf(fp, "\n");
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Check count */

              if (i != (long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH])
                Warn(FUNCTION_NAME, "Mismatch in count (%s) 5",
                    GetText(nuc + NUCLIDE_PTR_NAME));
            }
        }
    }

  /***************************************************************************/

  /***** Print fission yield data ********************************************/

  table++;

  if ((1 != 2) && ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES))
    {
      fprintf(fp, "\n --- Table %2ld: Fission yield data: \n", table);

      /* Loop over nuclides */

      n = 0;
      while ((nuc = ListPtr((long)RDB[DATA_PTR_NUC0], n++)) > VALID_PTR)
        {
          /* Check that either NFY or SFY data is given */

          if (((long)RDB[nuc + NUCLIDE_PTR_NFY_DATA] > VALID_PTR) ||
              ((long)RDB[nuc + NUCLIDE_PTR_SFY_DATA] > VALID_PTR))
            {
              /* Name and symbol */

              sprintf(sym, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 1));
              sprintf(name, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 2));

              fprintf(fp, "\n Nuclide %s -- %s (%s) ",
                      GetText(nuc + NUCLIDE_PTR_NAME), name, sym);

              /***************************************************************/

              /***** Neutron-induced fission *********************************/

              /* Pointer to distribution */

              if ((yld = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA]) > VALID_PTR)
                {
                  /* Check that yield is own */

                  if ((long)RDB[yld + FISSION_YIELD_PARENT_ZAI] ==
                      (long)RDB[nuc + NUCLIDE_ZAI])
                    {
                      fprintf(fp, "\n");

                      /* Reset count */

                      j = 0;

                      /* Loop over distributions */

                      while (yld > VALID_PTR)
                        {
                          /* Get number of fission products */

                          tot = (long)RDB[yld + FISSION_YIELD_NFP];

                          fprintf(fp,
                                  "\n Neutron-induced fission:\n\n");

                          fprintf(fp, " Interpolation region     : %ld / %ld\n", 1 + j++,
                                  (long)RDB[nuc + NUCLIDE_NFY_NE]);
                          fprintf(fp, " Interpolation energy     : %1.5E MeV\n",
                              RDB[yld + FISSION_YIELD_E]);
                          fprintf(fp, " Number of products       : %ld / %ld\n", tot,
                                  (long)RDB[yld + FISSION_YIELD_ORIG_NFP]);
                          fprintf(fp, " Cumulative yield cut-off : %1.5E\n\n",
                                  RDB[DATA_DEP_FP_YIELD_CUTOFF]);

                          fprintf(fp, " Columns:\n\n");
                          fprintf(fp, " 1. Index\n");
                          fprintf(fp, " 2. Product ZAI\n");
                          fprintf(fp, " 3. Product name\n");
                          fprintf(fp, " 4. Independent yield\n");
                          fprintf(fp, " 5. Cumulative yield\n");
                          fprintf(fp, " 6. Target in data\n\n");

                          /*
                          printf("i%lde%ld = zeros(200,1);\n", (long)RDB[nuc + NUCLIDE_ZAI], j);
                          printf("i%lde%ldE = %E;\n", (long)RDB[nuc + NUCLIDE_ZAI], j, RDB[yld + FISSION_YIELD_E]);
                          */
                          /* Pointer to distribution */

                          ptr = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];

                          /* Reset count */

                          i = 0;

                          /* Loop over distribution */

                          while ((loc0 = ListPtr(ptr, i++)) > VALID_PTR)
                            {
                              /* Get pointer to target */

                              if ((tgt = (long)RDB[loc0 + FY_PTR_TGT]) > VALID_PTR)
                                fprintf(fp, " %4ld  %6ld  %-7s  %11.5E  %11.5E  %s\n",
                                        i, (long)RDB[loc0 + FY_TGT_ZAI],
                                        ZAItoIso((long)RDB[loc0 + FY_TGT_ZAI],1),
                                        RDB[loc0 + FY_INDEPENDENT_FRAC],
                                        RDB[loc0 + FY_CUMULATIVE_FRAC],
                                        GetText(tgt + NUCLIDE_PTR_NAME));
                              else
                                fprintf(fp, " %4ld  %6ld  %-7s  %11.5E  %11.5E  (path not processed)\n",
                                        i, (long)RDB[loc0 + FY_TGT_ZAI],
                                        ZAItoIso((long)RDB[loc0 + FY_TGT_ZAI],1),
                                        RDB[loc0 + FY_INDEPENDENT_FRAC],
                                        RDB[loc0 + FY_CUMULATIVE_FRAC]);

                              /*
                              m = (long)(RDB[loc0 + FY_TGT_ZAI]/10.0);
                              m = m - 1000*((long)(((double)m)/1000.0));
                              printf("i%lde%ld(%ld) = i%lde%ld(%ld) + %E;\n", (long)RDB[nuc + NUCLIDE_ZAI], j, m, (long)RDB[nuc + NUCLIDE_ZAI], j, m, RDB[loc0 + FY_INDEPENDENT_FRAC]);
                              */
                              /* Tää printtaa prekursorit */
                              /*
                              if ((j == 1) && (RDB[loc0 + FY_INDEPENDENT_FRAC] > 1E-4))
                                if ((long)RDB[tgt + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DELNU_PREC)
                                  printf("%E %E %E %% %10s %% delnu\n",
                                         RDB[loc0 + FY_INDEPENDENT_FRAC],
                                         RDB[loc0 + FY_CUMULATIVE_FRAC],
                                         LOG2/RDB[tgt + NUCLIDE_LAMBDA],
                                         ZAItoIso((long)RDB[tgt + NUCLIDE_ZAI],1)
                                         );
                              */
                            }

                          /* Check count */

                          if (i - 1 != tot)
                            Warn(FUNCTION_NAME, "Mismatch in count 6");

                          /* Next */

                          yld = NextItem(yld);
                        }

                      /* Check count */

                      if (j != (long)RDB[nuc + NUCLIDE_NFY_NE])
                        Warn(FUNCTION_NAME, "Mismatch in count 7");
                    }
                  else
                    fprintf(fp, "is using NFY data from %s\n",
                            ZAItoIso((long)RDB[yld + FISSION_YIELD_PARENT_ZAI], 1));
                }

              /***************************************************************/

              /***** Spontaneous fission *************************************/

              /* Pointer to distribution */

              if ((yld = (long)RDB[nuc + NUCLIDE_PTR_SFY_DATA]) > VALID_PTR)
                {
                  /* Check pointer */

                  if (yld > VALID_PTR)
                    {
                      /* Get number of fission products */

                      tot = (long)RDB[yld + FISSION_YIELD_NFP];

                      fprintf(fp,
                              "\n Spontaneous fission:\n\n");

                      fprintf(fp, " Decay constant           : %1.5E\n",
                              RDB[nuc + NUCLIDE_LAMBDA]);
                      fprintf(fp, " SF decay fraction        : %1.5E\n",
                              RDB[nuc + NUCLIDE_SF_BR]);
                      fprintf(fp, " Number of products       : %ld / %ld\n", tot,
                              (long)RDB[yld + FISSION_YIELD_ORIG_NFP]);
                      fprintf(fp, " Cumulative yield cut-off : %1.5E\n\n",
                              RDB[DATA_DEP_FP_YIELD_CUTOFF]);

                      fprintf(fp, " Columns:\n\n");
                      fprintf(fp, " 1. Index\n");
                      fprintf(fp, " 2. Product ZAI\n");
                      fprintf(fp, " 3. Product name\n");
                      fprintf(fp, " 4. Independent yield\n");
                      fprintf(fp, " 5. Cumulative yield\n");
                      fprintf(fp, " 6. Target in data\n\n");

                      /* Pointer to distribution */

                      ptr = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];

                      /* Reset count */

                      i = 0;

                      /* Loop over distribution */

                      while ((loc0 = ListPtr(ptr, i++)) > VALID_PTR)
                        {
                          /* Get pointer to target */

                          if ((tgt = (long)RDB[loc0 + FY_PTR_TGT]) > VALID_PTR)
                            fprintf(fp, " %4ld  %6ld  %-7s  %11.5E  %11.5E  %s\n",
                                    i, (long)RDB[loc0 + FY_TGT_ZAI],
                                    ZAItoIso((long)RDB[loc0 + FY_TGT_ZAI],1),
                                    RDB[loc0 + FY_INDEPENDENT_FRAC],
                                    RDB[loc0 + FY_CUMULATIVE_FRAC],
                                    GetText(tgt + NUCLIDE_PTR_NAME));
                          else
                            fprintf(fp, " %4ld  %6ld  %-7s  %11.5E  %11.5E  (path not processed)\n",
                                    i, (long)RDB[loc0 + FY_TGT_ZAI],
                                    ZAItoIso((long)RDB[loc0 + FY_TGT_ZAI],1),
                                    RDB[loc0 + FY_INDEPENDENT_FRAC],
                                    RDB[loc0 + FY_CUMULATIVE_FRAC]);
                        }

                      /* Check count */

                      if (i - 1 != tot)
                        Warn(FUNCTION_NAME, "Mismatch in count 8");
                    }
                }

              /***************************************************************/
            }
        }
    }

  /***************************************************************************/

  /***** Print dead chains and lost data *************************************/

  table++;

  if ((1 != 2) && ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES))
    {
      fprintf(fp, "\n --- Table %2ld: Lost transmutation paths: \n\n", table);

      if (RDB[DATA_N_DEAD_PATH] > 0.0)
        fprintf(fp, " %ld lost transmutation paths:\n\n",
                (long)RDB[DATA_N_DEAD_PATH]);
      else
        fprintf(fp, " No lost transmutation paths.\n\n");

      /* Reset count */

      i = 0;

      /* Loop over nuclides */

      n = 0;
      while ((nuc = ListPtr((long)RDB[DATA_PTR_NUC0], n++)) > VALID_PTR)
        {
          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Get mt */

              mt = (long)RDB[rea + REACTION_MT];

              /* Get pointer to target */

              ptr = (long)RDB[rea + REACTION_PTR_TGT];
              if (ptr == (long)RDB[DATA_PTR_NUCLIDE_LOST])
                {
                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_DECAY)
                    {
                      j = mt - 10000;

                      if ((long)RDB[rea + REACTION_RTYP2] > 10000)
                        j = 10*j + (long)RDB[rea + REACTION_RTYP2] - 10000;
                      if ((long)RDB[rea + REACTION_RTYP3] > 10000)
                        j = 10*j + (long)RDB[rea + REACTION_RTYP3] - 10000;
                      if ((long)RDB[rea + REACTION_RTYP4] > 10000)
                        j = 10*j + (long)RDB[rea + REACTION_RTYP4] - 10000;
                      if ((long)RDB[rea + REACTION_RTYP5] > 10000)
                        j = 10*j + (long)RDB[rea + REACTION_RTYP5] - 10000;

                      if ((long)RDB[rea + REACTION_TGT_ZAI] == -2)
                        fprintf(fp, " %5ld  %10s  RTYP = %3ld           frac = %11.5E  missing target = %7s  :  SF without yield\n", 1 + i++,
                                GetText(nuc + NUCLIDE_PTR_NAME), j, RDB[rea + REACTION_BR],
                                ZAItoIso((long)RDB[rea + REACTION_TGT_ZAI], 1));
                      else
                        fprintf(fp, " %5ld  %10s  RTYP = %3ld           frac = %11.5E  missing target = %7s  :  %s\n", 1 + i++,
                                GetText(nuc + NUCLIDE_PTR_NAME), j, RDB[rea + REACTION_BR],
                                ZAItoIso((long)RDB[rea + REACTION_TGT_ZAI], 1), ReactionMT(mt, NO));
                    }
                  else
                    fprintf(fp, " %5ld  %10s    MT = %3ld           frac = %11.5E  missing target = %7s  :  %s\n", 1 + i++,
                            GetText(nuc + NUCLIDE_PTR_NAME), mt, RDB[rea + REACTION_BR],
                            ZAItoIso((long)RDB[rea + REACTION_TGT_ZAI], 1), ReactionMT(mt, NO));
                }

              /* Next */

              rea = NextItem(rea);
            }

          /* Loop over NFY yields */

          yld = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA];

          while (yld > VALID_PTR)
            {
              /* Check that yield is own */

              if ((long)RDB[yld + FISSION_YIELD_PARENT_ZAI] ==
                  (long)RDB[nuc + NUCLIDE_ZAI])
                {
                  /* Pointer to distribution */

                  ptr = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];

                  /* Loop over yield */

                  j = 0;
                  while ((loc0 = ListPtr(ptr, j++)) > VALID_PTR)
                    {
                      /* Get pointer to target */

                      tgt = (long)RDB[loc0 + FY_PTR_TGT];

                      if (tgt == (long)RDB[DATA_PTR_NUCLIDE_LOST])
                        fprintf(fp, " %5ld  %10s  Fission              frac = %11.5E  missing target = %7s\n", 1 + i++,
                                GetText(nuc + NUCLIDE_PTR_NAME), RDB[loc0 + FY_INDEPENDENT_FRAC],
                                ZAItoIso((long)RDB[loc0 + FY_TGT_ZAI], 1));
                    }
                }

              /* Next */

              yld = NextItem(yld);
            }

          /* Loop over NFY yields */

          yld = (long)RDB[nuc + NUCLIDE_PTR_SFY_DATA];

          while (yld > VALID_PTR)
            {
              /* Check that yield is own */

              if ((long)RDB[yld + FISSION_YIELD_PARENT_ZAI] ==
                  (long)RDB[nuc + NUCLIDE_ZAI])
                {
                  /* Pointer to distribution */

                  ptr = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];

                  /* Loop over yield */

                  j = 0;
                  while ((loc0 = ListPtr(ptr, j++)) > VALID_PTR)
                    {
                      /* Get pointer to target */

                      tgt = (long)RDB[loc0 + FY_PTR_TGT];

                      if (tgt == (long)RDB[DATA_PTR_NUCLIDE_LOST])
                        fprintf(fp, " %5ld  %10s  Spontaneous fission  frac = %11.5E  missing target = %7s\n", 1 + i++,
                                GetText(nuc + NUCLIDE_PTR_NAME), RDB[loc0 + FY_INDEPENDENT_FRAC],
                                ZAItoIso((long)RDB[loc0 + FY_TGT_ZAI], 1));
                    }
                }

              /* Next */

              yld = NextItem(yld);
            }
        }

      /* Check count */

      if (i != (long)RDB[DATA_N_DEAD_PATH])
        Warn(FUNCTION_NAME, "Mismatch in count 9 (%ld %ld)",
             i, (long)RDB[DATA_N_DEAD_PATH]);
    }

  /***************************************************************************/

  /* Close file */

  fclose(fp);

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
