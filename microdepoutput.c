/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : microdepoutput.c                               */
/*                                                                           */
/* Created:       2017/01/01 (JLe)                                           */
/* Last modified: 2019/11/15 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints micro-depletion output                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MicroDepOutout:"

/*****************************************************************************/

void MicroDepOutput()
{
  long loc0, loc1, loc2, ptr, mat, ntot, nr, i, j, gcu, uni, ZAI, MT, nuc, yld;
  long rea, n, I, FY[100], tgt;
  double E;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check micro-depletion data */

  if ((long)RDB[DATA_PTR_MDEP0] < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

 /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check branch index */

  if ((long)RDB[DATA_COEF_CALC_IDX] < 0)
    sprintf(tmpstr, "%s_mdx%ld.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_BURN_STEP]);
  else
    sprintf(tmpstr, "%s_mdx%ldb%ld.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_COEF_CALC_BU_IDX],
            (long)RDB[DATA_COEF_CALC_IDX]);

  /* Open file for writing */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /***************************************************************************/

  /***** Print fission yields ************************************************/

  /* Reset count */

  nr = 0;

  /* Loop over micro depletion data and collect fissionable actinides */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over material data */

      loc1 = (long)RDB[loc0 + MDEP_PTR_MAT];
      while (loc1 > VALID_PTR)
        {
          /* Loop over reactions */

          loc2 = (long)RDB[loc1 + MDEP_MAT_PTR_REA];
          while (loc2 > VALID_PTR)
            {
              /* Pointer to reaction */

              rea = (long)RDB[loc2 + MDEP_REA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Get mt */

              MT = (long)RDB[rea + REACTION_MT];

              /* Pointer to nuclide */

              nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Add fission reactions to list */

              if (((MT > 17) && (MT < 22)) || (MT == 38))
                if ((yld = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA]) > VALID_PTR)
                  {
                    /* Check if already added */

                    for (n = 0; n < nr; n++)
                      if ((long)RDB[FY[n] + NUCLIDE_ZAI] ==
                          (long)RDB[nuc + NUCLIDE_ZAI])
                        break;

                    /* check */

                    if (n == nr)
                      FY[nr++] = nuc;
                  }

              /* Check */

              if (nr == 100)
                Die(FUNCTION_NAME, "Array overflow");

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      loc0 = NextItem(loc0);
    }

  /* Loop over list */

  for (n = 0; n < nr; n++)
    {
      /* Get pointer to nuclide */

      nuc = FY[n];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Pointer to distribution */

      if ((yld = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA]) > VALID_PTR)
        {
          /* Check that yield is own */

          if ((long)RDB[yld + FISSION_YIELD_PARENT_ZAI] !=
              (long)RDB[nuc + NUCLIDE_ZAI])
            fprintf(fp, "%% NOTE: %ld is using yield from %ld\n\n",
                    (long)RDB[nuc + NUCLIDE_ZAI],
                    (long)RDB[yld + FISSION_YIELD_PARENT_ZAI]);

          /* Loop over distributions */

          while (yld > VALID_PTR)
            {
              /* Get energy and index */

              E = RDB[yld + FISSION_YIELD_E];
              j = (long)RDB[yld + FISSION_YIELD_IDX];

              /* Print energy */

              fprintf(fp, "NFY_%ld_%ldE = %1.5E ;\n",
                      (long)RDB[nuc + NUCLIDE_ZAI], j, E);

              /* Print yields */

              fprintf(fp, "NFY_%ld_%ld = [\n",
                      (long)RDB[nuc + NUCLIDE_ZAI], j);

              /* Loop over distribution */

              ptr = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];
              while (ptr > VALID_PTR)
                {
                  /* Get ZAI */

                  if ((long)RDB[ptr + FY_PTR_TGT] ==
                      (long)RDB[DATA_PTR_NUCLIDE_LOST])
                    ZAI = -1;
                  else
                    ZAI = (long)RDB[ptr + FY_TGT_ZAI];

                  /* Print */

                  fprintf(fp, "%6ld %11.5E  %11.5E %% %s\n", ZAI,
                          RDB[ptr + FY_INDEPENDENT_FRAC],
                          RDB[ptr + FY_CUMULATIVE_FRAC],
                          ZAItoIso((long)RDB[ptr + FY_TGT_ZAI],1));

                  /* Next */

                  ptr = NextItem(ptr);
                }

              /* Close bracket */

              fprintf(fp, "];\n\n");

              /* Next */

              yld = NextItem(yld);
            }
        }
    }

  /***************************************************************************/

  /***** Print decay data ****************************************************/

  fprintf(fp, "%% Decay data, columns:\n");
  fprintf(fp, "%% 1. ZAI\n");
  fprintf(fp, "%% 2. decay constant (1/s)\n");
  fprintf(fp, "%% 3. specific decay energy (J)\n");
  fprintf(fp, "%% 4. reaction type\n");
  fprintf(fp, "%% 5. branch fraction\n");
  fprintf(fp, "%% 6. product ZAI\n\n");

  fprintf(fp, "dec = [\n");

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check that nuclide is included in burnup calculation */

      if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP))
        {
          /* Next nuclide */

          nuc = NextItem(nuc);

          /* Cycle loop */

          continue;
        }

      /* Loop over decay reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_DECAY)
            {
              /* Next reaction */

              rea = NextItem(rea);

              /* Cycle loop */

              continue;
            }

          /* Print ZAI */

          ZAI = (long)RDB[nuc + NUCLIDE_ZAI];
          fprintf(fp, "%6ld  ", ZAI);

          /* Print decay constant */

          fprintf(fp, "%11.5E  ", RDB[nuc + NUCLIDE_LAMBDA]);

          /* Print specific decay energy */

          fprintf(fp, "%11.5E  ", RDB[nuc + NUCLIDE_DECAY_E]*MEV);

          /* Print mt */

          MT = (long)RDB[rea + REACTION_MT];
          n = MT - 10000;

          /* Add partials */

          if ((long)RDB[rea + REACTION_RTYP2] > 10000)
            n = 10*n + (long)RDB[rea + REACTION_RTYP2] - 10000;
          if ((long)RDB[rea + REACTION_RTYP3] > 10000)
            n = 10*n + (long)RDB[rea + REACTION_RTYP3] - 10000;
          if ((long)RDB[rea + REACTION_RTYP4] > 10000)
            n = 10*n + (long)RDB[rea + REACTION_RTYP4] - 10000;
          if ((long)RDB[rea + REACTION_RTYP5] > 10000)
            n = 10*n + (long)RDB[rea + REACTION_RTYP5] - 10000;

          fprintf(fp, "%5ld  ", n);

          /* Print fraction */

          fprintf(fp, "%11.5E  ", RDB[rea + REACTION_BR]);

          /* Get target */

          tgt = ReactionTargetZAI(rea);

          /* Print target */

          if (tgt > 0)
            fprintf(fp, "%6ld", tgt);
          else
            fprintf(fp, "%6d", 0);

          /* Print comment */

          sprintf(tmpstr, "%s", ZAItoIso(tgt, 1));

          if (tgt > 0)
            {
              /* Check partial */

              if (n != MT - 10000)
                fprintf(fp, "  %% %s to %s (multiple reactions)\n",
                        ZAItoIso(ZAI, 1), tmpstr);
              else
                fprintf(fp, "  %% %s %s to %s\n", ZAItoIso(ZAI, 1),
                        ReactionMT(MT, NO), tmpstr);
            }
          else if (tgt == -2)
            fprintf(fp, "  %% %s spontaneous fission\n", ZAItoIso(ZAI, 1));
          else
            Die(FUNCTION_NAME, "tgt = %ld (%ld RTYP %ld)", tgt, ZAI, n);

          /* Next reaction */

          rea = NextItem(rea);
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  fprintf(fp, "];\n\n");

  /***************************************************************************/

  /***** Print cross sections ************************************************/

  /* Loop over micro depletion data */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to gcu data */

      gcu = (long)RDB[loc0 + MDEP_PTR_GCU];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

      /* Pointer to universe */

      uni = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Print universe */

      fprintf(fp, "%% Universe: %s\n", GetText(uni + UNIVERSE_PTR_NAME));

      /* Print materials */

      fprintf(fp, "%% materials:");

      /* Reset previous */

      ptr = -1;

      /* Loop over data */

      loc1 = (long)RDB[loc0 + MDEP_PTR_MAT];
      while (loc1 > VALID_PTR)
        {
          /* Material pointer */

          mat = (long)RDB[loc1 + MDEP_MAT_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
            mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];

          /* Print name */

          if (mat != ptr)
            fprintf(fp, " %s", GetText(mat + MATERIAL_PTR_NAME));

          /* remember previous */

          ptr = mat;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Newline */

      fprintf(fp, "\n");

      /* Number of energy groups */

      ntot = (long)RDB[DATA_ERG_FG_NG];
      CheckValue(FUNCTION_NAME, "ntot", "", ntot, 1, 100000);

      /* Number of reactions */

      nr = (long)RDB[loc0 + MDEP_N_REA];
      CheckValue(FUNCTION_NAME, "nr", "", nr, 1, 100000);

      /* Pointer to search keys */

      loc1 = (long)RDB[loc0 + MDEP_PTR_KEY];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Pointer to flux */

      ptr = (long)RDB[gcu + GCU_INF_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Print flux ratio */

      fprintf(fp, "\nFLUX_%s = [", GetText(uni + UNIVERSE_PTR_NAME));

      /* Loop over groups and print results */

      for (j = 0; j < ntot; j++)
        fprintf(fp, " %12.5E %1.5f", Mean(ptr, j)/RDB[loc0 + MDEP_VOLUME],
                RelErr(ptr, j));

      /* Close bracket */

      fprintf(fp, " ];\n");

      /* Pointer to results */

      ptr = (long)RDB[gcu + GCU_RES_FG_MICRO_DEP_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to average densities */

      loc2 = RDB[gcu + GCU_RES_FG_MICRO_DEP_ADENS];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Print cross sections */

      fprintf(fp, "\nXS_%s = [\n", GetText(uni + UNIVERSE_PTR_NAME));

      /* Loop over reactions */

      for (i = 0; i < nr; i++)
        {
          /* Get ZAI, MT and I */

          ZAI = (long)(RDB[loc1 + i]/10000.0);
          MT = (long)((RDB[loc1 + i] - (double)ZAI*10000.0)/10.0);
          I = (long)RDB[loc1 + i] - 10000*ZAI - 10*MT;

          /* Adjust I if fission */

          if (((MT == 18) || (MT == 19)) && (I > 0))
            I = I - 2;

          /* Print ZAI, MT and branch */

          fprintf(fp, "%6ld  %3ld  %ld  ", ZAI, MT, I);

          /* Print average density */

          fprintf(fp, "%11.5E %1.5f  ", Mean(loc2, i, j), RelErr(loc2, i, j));

          /* Loop over groups and print results */

          for (j = 0; j < ntot; j++)
            fprintf(fp, "%11.5E %1.5f  ", Mean(ptr, i, j), RelErr(ptr, i, j));

          fprintf(fp, "\n");
        }

      /* Close bracket */

      fprintf(fp, "];\n\n");

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /* Close file and exit */

  fclose(fp);
}

/*****************************************************************************/
