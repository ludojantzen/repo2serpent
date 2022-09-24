/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : linkreations.c                                 */
/*                                                                           */
/* Created:       2011/03/04 (JLe)                                           */
/* Last modified: 2020/05/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Links reactions to sources and detectors, etc.               */
/*                                                                           */
/* Comments: - Toi fissiotuotejakauman energiariippuvuus on hoidettu nyt     */
/*             tosi hölmösti. Siinä odotetaan että energia osuu tietylle     */
/*             välille, mikä ei välttämättä pidä paikkansa.                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LinkReactions:"

/*****************************************************************************/

void LinkReactions()
{
  long det, src, ptr, mat, rea, nuc, iso, mt, loc0, loc1, loc2, n, key, i, zai;
  long tst;
  long RFS, idx, ify;

  /***************************************************************************/

  /***** Detectors ***********************************************************/

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Loop over reaction bins */

      ptr = (long)RDB[det + DET_PTR_RBINS];
      while (ptr > VALID_PTR)
        {
          /* Get material pointer */

          mat = (long)RDB[ptr + DET_RBIN_PTR_MAT];

          /* Get mt */

          mt = (long)RDB[ptr + DET_RBIN_MT];

          /* Check for TLE type */

          if ((loc0 = (long)RDB[det + DET_PTR_SBINS]) > VALID_PTR)
            {
              /* Check for surface flux detector */

              if ((long)RDB[loc0 + DET_SBIN_SURF_NORM] == -2)
                {
                  /* Check response */

                  if ((mat < VALID_PTR) && (mt != 0) &&
                      (mt != MT_USER_DEFINED))
                    Error(det,
                          "Void response material not allowed with surface flux detector");
                  else if ((mt ==  MT_SOURCE_RATE) ||
                           (mt ==  MT_MACRO_RECOILE) ||
                           (mt ==  MT_MACRO_HEATPHOTANA) ||
                           (mt ==  MT_LEAK_RATE))
                    Error(det, "Response type %ld not allowed with surface flux detector", mt);
                }
              else if ((mt != 0.0) && (mt != MT_USER_DEFINED) &&
                       ((long)RDB[loc0 + DET_SBIN_TYPE] !=
                        SUPERDET_TYPE_TLEFLUX))
                Error(det,
                      "Response not allowed with surface current detector");
            }

          /* Reset reaction pointer */

          rea = -1;

          /* Check mt */

          if (mt == MT_MACRO_TOTXS)
            {
              if (mat > VALID_PTR)
                {
                  if ((long)RDB[det + DET_PARTICLE] == PARTICLE_TYPE_NEUTRON)
                    rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
                  else
                    rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
                }
            }
          else if (mt == MT_MACRO_ABSXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_ABSXS];
            }
          else if (mt == MT_MACRO_ELAXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_ELAXS];
            }
          else if (mt == MT_MACRO_INLPRODXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_INLPXS];
            }
          else if (mt == MT_MACRO_FISSXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_FISSXS];
            }
          else if (mt == MT_MACRO_HEATXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_HEATTXS];
            }
          else if (mt == MT_MACRO_PHOTXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_PHOTPXS];
            }
          else if (mt == MT_MACRO_PROTPXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_PROTPXS];
            }
          else if (mt == MT_MACRO_DEUTPXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_DEUTPXS];
            }
          else if (mt == MT_MACRO_TRITPXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_TRITPXS];
            }
          else if (mt == MT_MACRO_HE3PXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_HE3PXS];
            }
          else if (mt == MT_MACRO_HE4PXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_HE4PXS];
            }
          else if (mt == MT_MACRO_FISSE)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_FISSE];
            }
          else if (mt == MT_MACRO_NSF)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_NSF];
            }
          else if (mt == MT_MACRO_RECOILE)
            {
              rea = -1;
            }
          else if (mt == MT_MACRO_HEATPHOTANA)
            {
              rea = -1;
            }
          else if (mt == MT_SOURCE_RATE)
            {
              rea = -1;
            }
          else if (mt == MT_LEAK_RATE)
            {
              rea = -1;
            }
          else if (mt == MT_MACRO_MAJORANT)
            {
              rea = -1;
            }
          else if (mt == MT_NEUTRON_DENSITY)
            {
              rea = -1;
            }
          else if ((mt <= MT_PRIMARY_LIVE_SOURCE) &&
                   (mt >= MT_SECONDARY_DN_SOURCE_G8))
            {
              rea = -1;
            }
          else if (mt == MT_USER_DEFINED)
            {
              rea = -1;
            }
          else if (mt == MT_PHOTON_PULSE_HEIGHT)
            {
              rea = -1;
            }
          else if ((mt == MT_PHOTON_DOSE) || ((mt > -249) && (mt < -200)))
            {
              rea = -1;
            }
          else if ((mt >= MT_ELECTRON_AUGER) && (mt <= MT_ELECTRON_PE))
            {
              rea = -1;
            }
          else if (mt == MT_MACRO_HEATTOT)
            {
              rea = -1;
            }
          else if (mt == MT_MACRO_TOTPHOTXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
            }
          else if (mt == MT_MACRO_TMP_MAJORANTXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS];
            }
          else if (mt == MT_MACRO_HEATPHOTXS)
            {
              if (mat > VALID_PTR)
                rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS];
            }
          else if (mt < 0)
            Error(det, "MT %ld not allowed in response function", mt);
          else if (mt > 0)
            {
              /* Check material pointer */

              if (mat < VALID_PTR)
                Error(det, "Response function with mt %ld must be associated with a material", mt);

              /* Pointer to composition */

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

              /* Check number of nuclides */

              if (NextItem(iso) > VALID_PTR)
                Error(det,  "Material %s used in response must consist of single nuclide", GetText(mat + MATERIAL_PTR_NAME));

              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check total and total absorption */

              if (mt == 1)
                rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              else if (mt == 101)
                rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
              else
                {
                  /* Find matching mt */

                  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
                  while (rea > VALID_PTR)
                    {
                      /* Compare */

                      if ((long)RDB[rea + REACTION_MT] == mt)
                        break;

                      /* Pointer to next */

                      rea = NextItem(rea);
                    }

                  /* Check pointer */

                  if (rea < VALID_PTR)
                    {
                      fprintf(stdout,
                              "\nNuclide %s has the following reactions:\n\n",
                              GetText(nuc + NUCLIDE_PTR_NAME));

                      if ((long)RDB[nuc + NUCLIDE_PTR_TOTXS] > VALID_PTR)
                        fprintf(stdout, "MT %-6ld : %s\n", (long)1,
                                ReactionMT(1, YES));
                      if ((long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS] > VALID_PTR)
                        fprintf(stdout, "MT %-6ld : %s\n", (long)101,
                                ReactionMT(101, YES));

                      /* Loop over reactions */

                      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
                      while (rea > VALID_PTR)
                        {
                          /* Check type */

                          if (((long)RDB[rea + REACTION_TYPE] ==
                               REACTION_TYPE_PARTIAL) ||
                              ((long)RDB[rea + REACTION_TYPE] ==
                               REACTION_TYPE_SPECIAL))
                            {
                              if ((long)RDB[rea + REACTION_RFS] == 0)
                                fprintf(stdout, "MT %-6ld : %s\n",
                                        (long)RDB[rea + REACTION_MT],
                                        ReactionMT((long)RDB[rea + REACTION_MT], YES));
                              else
                                fprintf(stdout, "MT %-6ld : %s to isomeric state\n",                                         (long)(10*RDB[rea + REACTION_MT] + 1),
                                        ReactionMT((long)RDB[rea + REACTION_MT], YES));
                            }

                          /* Next reaction */

                          rea = NextItem(rea);
                        }

                      Error(det,
                            "Reaction mt %ld not found for response function",
                            (long)RDB[ptr + DET_RBIN_MT]);
                    }
                }

              /* Check branching to isomeric state */

              if (((long)RDB[ptr + DET_RBIN_RFS] > 0) &&
                  ((long)RDB[rea + REACTION_PTR_ISO_BRA] < VALID_PTR))
                Error(det, "Reaction mt %ld has no isomeric branching data",
                      (long)RDB[ptr + DET_RBIN_MT]);
            }

          /* Check pointer */

          if ((mt != 0) && (mt != MT_PHOTON_DOSE) && (rea < VALID_PTR) &&
              (mat > VALID_PTR))
            {
              if ((mt == MT_MACRO_RECOILE) || (mt == MT_MACRO_HEATPHOTANA))
                Error(det, "Material entry must be void with mt %ld", mt);
              else if ((long)RDB[mat + MATERIAL_DIV_TYPE] ==
                       MAT_DIV_TYPE_PARENT)
                Error(det,
                      "Reaction mt %ld not allowed with divided/burnable material %s",
                     mt, GetText(mat + MATERIAL_PTR_NAME));
              else if ((mt != MT_MACRO_FISSXS) &&
                       (mt != MT_MACRO_FISSE) &&
                       (mt != MT_MACRO_PROTPXS) &&
                       (mt != MT_MACRO_DEUTPXS) &&
                       (mt != MT_MACRO_TRITPXS) &&
                       (mt != MT_MACRO_HE3PXS) &&
                       (mt != MT_MACRO_HE4PXS))
                Error(det, "Reaction mt %ld not found for response function",
                      mt);
              else
                Note(det, "Reaction mt %ld not found for response function",
                      mt);
            }

          /* Set pointer */

          WDB[ptr + DET_RBIN_PTR_REA] = (double)rea;

          /* Next bin */

          ptr = NextItem(ptr);
        }

      /* Next detector */

      det = NextItem(det);
    }

  /***************************************************************************/

  /***** Sources *************************************************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Pointer to nuclide */

      if ((nuc = (long)RDB[src + SRC_PTR_XSDATA]) > VALID_PTR)
        {
          /* Get mt */

          mt = (long)RDB[src + SRC_PTR_REA];

          /* Find matching mt */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          if ((rea = SeekList(rea, REACTION_MT, (double)mt, NO)) < VALID_PTR)
            {
              fprintf(stdout,
                      "\nNuclide %s has the following reactions:\n\n",
                      GetText(nuc + NUCLIDE_PTR_NAME));

              if ((long)RDB[nuc + NUCLIDE_PTR_TOTXS] > VALID_PTR)
                fprintf(stdout, "MT %-6ld : %s\n", (long)1, ReactionMT(1, YES));
              if ((long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS] > VALID_PTR)
                fprintf(stdout, "MT %-6ld : %s\n", (long)101, ReactionMT(101, YES));

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Check type */

                  if (((long)RDB[rea + REACTION_TYPE] ==
                       REACTION_TYPE_PARTIAL) ||
                      ((long)RDB[rea + REACTION_TYPE] ==
                       REACTION_TYPE_SPECIAL))
                    fprintf(stdout, "MT %-6ld : %s\n",
                            (long)RDB[rea + REACTION_MT],
                            ReactionMT((long)RDB[rea + REACTION_MT], YES));

                  /* Next reaction */

                  rea = NextItem(rea);
                }

              Error(src, "Reaction mt %ld not found for source distribution",
                    mt);
            }

          /* Put pointer */

          WDB[src + SRC_PTR_REA] = (double)rea;

          /* Check that reaction has distribution data */

          if ((long)RDB[rea + REACTION_PTR_ERG] < VALID_PTR)
            Error(src, "Reaction mt %ld has no distribution data", mt);
        }

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Micro-depletion stuff **********************************************/

  /* Loop over definitions */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to reaction search keys */

      ptr = (long)RDB[loc0 + MDEP_PTR_KEY];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Link compositions */

      loc1 = (long)RDB[loc0 + MDEP_PTR_MAT];
      while (loc1 > VALID_PTR)
        {
          mat = (long)RDB[loc1 + MDEP_MAT_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Add to volume */

          if (RDB[mat + MATERIAL_VOLUME] == 0.0)
            Error(loc0, "Material %s volume is zero",
                  GetText(mat + MATERIAL_PTR_NAME));
          else
            WDB[loc0 + MDEP_VF] = RDB[loc0 + MDEP_VF]
              + RDB[mat + MATERIAL_VOLUME];

          /* Loop over composition and link reactions */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Check type */

                  if (((long)RDB[rea + REACTION_TYPE] !=
                       REACTION_TYPE_PARTIAL) &&
                      ((long)RDB[rea + REACTION_TYPE] !=
                       REACTION_TYPE_TRA_BRANCH) &&
                      ((long)RDB[rea + REACTION_MT] > 999))
                    {
                      /* Next reaction */

                      rea = NextItem(rea);

                      /* Cycle loop */

                      continue;
                    }

                  /* Calculate search key */

                  key = ((long)RDB[nuc + NUCLIDE_ZAI])*1000 +
                    (long)RDB[rea + REACTION_MT];

                  /* Get isomeric state flag */

                  RFS = (long)RDB[rea + REACTION_RFS];

                  /* Get fission yield index */

                  ify = -1;
                  if ((loc2 = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
                    ify = (long)RDB[loc2 + FISSION_YIELD_IDX];

                  /* Search from list */

                  n = 0;
                  while ((long)RDB[ptr + n] > VALID_PTR)
                    {
                      /* Separate branch identifier */

                      i = (long)(RDB[ptr + n]/10.0);
                      i = (long)RDB[ptr + n] - 10*i;

                      /* Separate mt */

                      mt = (long)((RDB[ptr + n] - (double)i)/10000.0);
                      mt = (long)(RDB[ptr + n] - 10000*mt)/10;

                      /* Separate ZAI */

                      zai = (long)((RDB[ptr + n] - 10*mt - i)/10000.0);

                      /* Create test id */

                      if ((mt == 452) || (mt == 455) || (mt == 456) ||
                          (mt == 458))
                        {
                          /* Nubar or kappa, use fission */

                          tst = 1000*zai + 18;
                        }
                      else
                        tst = 1000*zai + mt;

                      /* Compare */
                      /*
                      if ((long)(RDB[ptr + n]/10.0) == key)
                      */

                      if (tst == key)
                        {
                          /* Check branch identifiers */

                          if ((i == 0) ||
                              ((i == 1) && (RFS == 0)) ||
                              ((i == 2) && (RFS == 1)) ||
                              /*
                              ((i == 3) && (fabs(IE/2.53E-08 - 1.0) < 0.01)) ||
                              ((i == 4) && (fabs(IE/5.00E-01 - 1.0) < 0.01)) ||
                              ((i == 5) && (fabs(IE/1.40E+01 - 1.0) < 0.01)))
                              ((i == 3) && (IE <= 1E-6)) ||
                              ((i == 4) && ((IE > 1E-6) && (IE <= 1.0))) ||
                              ((i == 5) && ((IE > 1.0))))
                              */
                              ((i > 2) && (i - 2 == ify)))
                            {
                              /* Allocate memory */

                              loc2 = NewItem(loc1 + MDEP_MAT_PTR_REA,
                                             MDEP_REA_BLOCK_SIZE);

                              /* Put minimum energy */

                              if (RDB[rea + REACTION_FISSY_IE0] <
                                  RDB[rea + REACTION_EMIN])
                                WDB[loc2 + MDEP_REA_EMIN] =
                                  RDB[rea + REACTION_EMIN];
                              else
                                WDB[loc2 + MDEP_REA_EMIN] =
                                  RDB[rea + REACTION_FISSY_IE0];

                              /* Put mt */

                              WDB[loc2 + MDEP_REA_MT] = (double)mt;

                              /* Exclude branches from density calculation */

                              if ((long)RDB[rea + REACTION_TYPE] ==
                                  REACTION_TYPE_TRA_BRANCH)
                                WDB[loc2 + MDEP_REA_PARTIAL_IDX] = 1.0;

                              /* Put other data */

                              WDB[loc2 + MDEP_REA_IDX] = (double)n;
                              WDB[loc2 + MDEP_REA_PTR_ISO] = (double)iso;
                              WDB[loc2 + MDEP_REA_PTR_REA] = (double)rea;
                            }
                        }

                      /* Next */

                      n++;
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Next */

              iso = NextItem(iso);
            }

          /* Loop over composition and link total fission */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Reset partial index */

              idx = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Check type */

                  if (((long)RDB[rea + REACTION_TYPE] !=
                       REACTION_TYPE_PARTIAL) &&
                      ((long)RDB[rea + REACTION_TYPE] !=
                       REACTION_TYPE_TRA_BRANCH) &&
                      ((long)RDB[rea + REACTION_MT] > 999))
                    {
                      /* Next reaction */

                      rea = NextItem(rea);

                      /* Cycle loop */

                      continue;
                    }

                  /* Skip everything but partial fission */

                  if (((long)RDB[rea + REACTION_MT] != 19) &&
                      ((long)RDB[rea + REACTION_MT] != 20) &&
                      ((long)RDB[rea + REACTION_MT] != 21) &&
                      ((long)RDB[rea + REACTION_MT] != 38))
                    {
                      /* Next reaction */

                      rea = NextItem(rea);

                      /* Cycle loop */

                      continue;
                    }

                  /* Calculate search key */

                  key = ((long)RDB[nuc + NUCLIDE_ZAI])*1000 + 18;

                  /* Get isomeric state flag */

                  RFS = (long)RDB[rea + REACTION_RFS];

                  /* Get fission yield index */

                  ify = -1;
                  if ((loc2 = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
                    ify = (long)RDB[loc2 + FISSION_YIELD_IDX];

                  /* Search from list */

                  n = 0;
                  while ((long)RDB[ptr + n] > VALID_PTR)
                    {
                      /* Separate branch identifier */

                      i = (long)(RDB[ptr + n]/10.0);
                      i = (long)RDB[ptr + n] - 10*i;

                      /* Separate mt */

                      mt = (long)((RDB[ptr + n] - (double)i)/10000.0);
                      mt = (long)(RDB[ptr + n] - 10000*mt)/10;

                      /* Separate ZAI */

                      zai = (long)((RDB[ptr + n] - 10*mt - i)/10000.0);

                      /* Create test id */

                      if ((mt == 452) || (mt == 455) || (mt == 456) ||
                          (mt == 458))
                        {
                          /* Nubar or kappa, use fission */

                          tst = 1000*zai + 18;
                        }
                      else
                        tst = 1000*zai + mt;

                      /* Compare */
                      /*
                      if ((long)(RDB[ptr + n]/10.0) == key)
                      */

                      if (tst == key)
                        {
                          /* Check branch identifiers */

                          if ((i == 0) ||
                              ((i == 1) && (RFS == 0)) ||
                              ((i == 2) && (RFS == 1)) ||
                              /*
                              ((i == 3) && (fabs(IE/2.53E-08 - 1.0) < 0.01)) ||
                              ((i == 4) && (fabs(IE/5.00E-01 - 1.0) < 0.01)) ||
                              ((i == 5) && (fabs(IE/1.40E+01 - 1.0) < 0.01)))
                              ((i == 3) && (IE <= 1E-6)) ||
                              ((i == 4) && ((IE > 1E-6) && (IE <= 1.0))) ||
                              ((i == 5) && ((IE > 1.0))))
                              */
                              ((i > 2) && (i - 2 == ify)))
                            {
                              /* Allocate memory */

                              loc2 = NewItem(loc1 + MDEP_MAT_PTR_REA,
                                             MDEP_REA_BLOCK_SIZE);

                              /* Put minimum energy */

                              if (RDB[rea + REACTION_FISSY_IE0] <
                                  RDB[rea + REACTION_EMIN])
                                WDB[loc2 + MDEP_REA_EMIN] =
                                  RDB[rea + REACTION_EMIN];
                              else
                                WDB[loc2 + MDEP_REA_EMIN] =
                                  RDB[rea + REACTION_FISSY_IE0];

                              /* Put mt */

                              WDB[loc2 + MDEP_REA_MT] = (double)mt;

                              /* Put index */

                              WDB[loc2 + MDEP_REA_PARTIAL_IDX] =
                                (double)(idx++);

                              /* Put other data */

                              WDB[loc2 + MDEP_REA_IDX] = (double)n;
                              WDB[loc2 + MDEP_REA_PTR_ISO] = (double)iso;
                              WDB[loc2 + MDEP_REA_PTR_REA] = (double)rea;
                            }
                        }

                      /* Next */

                      n++;
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Next */

              iso = NextItem(iso);
            }

          /* Loop over composition and link total capture */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Reset partial index */

              idx = 0;

              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Check type */

                  if (((long)RDB[rea + REACTION_TYPE] !=
                       REACTION_TYPE_PARTIAL) &&
                      ((long)RDB[rea + REACTION_TYPE] !=
                       REACTION_TYPE_TRA_BRANCH) &&
                      ((long)RDB[rea + REACTION_MT] > 999))
                    {
                      /* Next reaction */

                      rea = NextItem(rea);

                      /* Cycle loop */

                      continue;
                    }

                  /* Skip everything but partial captures */

                  if (((long)RDB[rea + REACTION_MT] < 102) ||
                      ((long)RDB[rea + REACTION_MT] > 118))
                    {
                      /* Next reaction */

                      rea = NextItem(rea);

                      /* Cycle loop */

                      continue;
                    }

                  /* Calculate search key */

                  key = ((long)RDB[nuc + NUCLIDE_ZAI])*1000 + 101;

                  /* Search from list */

                  n = 0;
                  while ((long)RDB[ptr + n] > VALID_PTR)
                    {
                      /* Compare */

                      if ((long)(RDB[ptr + n]/10.0) == key)
                        {
                          /* Allocate memory */

                          loc2 = NewItem(loc1 + MDEP_MAT_PTR_REA,
                                         MDEP_REA_BLOCK_SIZE);

                          /* Put minimum energy */

                          if (RDB[rea + REACTION_FISSY_IE0] <
                              RDB[rea + REACTION_EMIN])
                            WDB[loc2 + MDEP_REA_EMIN] =
                              RDB[rea + REACTION_EMIN];
                          else
                            WDB[loc2 + MDEP_REA_EMIN] =
                              RDB[rea + REACTION_FISSY_IE0];

                          /* Put mt */

                          WDB[loc2 + MDEP_REA_MT] = 101.0;

                          /* Put index */

                          WDB[loc2 + MDEP_REA_PARTIAL_IDX] = (double)(idx++);

                          /* Put other data */

                          WDB[loc2 + MDEP_REA_IDX] = (double)n;
                          WDB[loc2 + MDEP_REA_PTR_ISO] = (double)iso;
                          WDB[loc2 + MDEP_REA_PTR_REA] = (double)rea;
                        }

                      /* Next */

                      n++;
                    }

                  /* Next */

                  rea = NextItem(rea);
                }

              /* Next */

              iso = NextItem(iso);
            }

          /* Sort list (is not created if there are no reactions) */

          if ((loc2 = (long)RDB[loc1 + MDEP_MAT_PTR_REA]) > VALID_PTR)
            SortList(loc2, MDEP_REA_EMIN, SORT_MODE_ASCEND);

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Number of reactions */

      n = (long)RDB[loc0 + MDEP_N_REA];
      CheckValue(FUNCTION_NAME, "n", "", n, 0, 1000000);

      /* Allocate memory for average densities */

      ptr = ReallocMem(DATA_ARRAY, n);
      WDB[loc0 + MDEP_PTR_AVG_ADENS] = (double)ptr;

      ptr = ReallocMem(DATA_ARRAY, n);
      WDB[loc0 + MDEP_PTR_BTCH_AVG_ADENS] = (double)ptr;

      /* Next definition */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
