/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readrestartfile.c                              */
/*                                                                           */
/* Created:       2014/01/23 (JLe)                                           */
/* Last modified: 2020/06/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads material compositions from a restart file before       */
/*              running the simulation.                                      */
/*                                                                           */
/* Comments: RESTART_CHECK    -- Check compatibility with materials          */
/*           RESTART_OVERRIDE -- Override material compositions (used in     */
/*                               burnup mode)                                */
/*           RESTART_REPLACE  -- Replace initial compositions (used in       */
/*                               transport mode)                             */
/*                                                                           */
/*           - Replace mode works only with photon data (neutron transport   */
/*             calculation requires full nuclide names to be stored and a    */
/*             revised file format)                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadRestartFile:"

/*****************************************************************************/

void ReadRestartFile(long mode)
{
  long sz, n, nnuc, zai, mat, iso, nuc, *ptr, ok, idx, i;
  double pt, bu, days, days0, adens, mdens, mbu, d0, closest, sum;
  char tmpstr[MAX_STR], fname[MAX_STR], ZAI[MAX_STR];
  FILE *fp;

  /* Check if file is read */

  if ((long)RDB[DATA_READ_RESTART_FILE] == NO)
    return;

  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    Die(FUNCTION_NAME, "Domain decomposition in use");

  /* File name */

  if ((long)RDB[DATA_RESTART_READ_PTR_FNAME] > VALID_PTR)
    sprintf(fname, "%s", GetText(DATA_RESTART_READ_PTR_FNAME));
  else
    sprintf(fname, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    {
      /* Check continue mode */

      if ((long)RDB[DATA_RESTART_READ_CONTINUE] == YES)
        {
          /* Start from beginning */

          fprintf(outp, "No restart file, starting from beginning...\n\n");

          /* Exit */

          return;
        }
      else
        Error(0, "Restart file \"%s\" does not exist", fname);
    }

  /* Print */

  if (mode == RESTART_CHECK)
    fprintf(outp, "Checking restart file \"%s\"...\n", fname);
  else
    fprintf(outp,
            "Reading material compositions from restart file \"%s\":\n\n",
            fname);

 /***************************************************************************/

  /***** Find point **********************************************************/

  /* Check if index is given */

  if ((idx = (long)RDB[DATA_RESTART_READ_IDX]) > 0)
    {
      /* Reset point */

      pt = -1;
    }
  else if ((long)RDB[DATA_RESTART_READ_CONTINUE] == YES)
    {
      /* Continue from previous, set values */

      idx = -1;
      pt = -1;
    }
  else
    {
      /* Get point */

      if ((pt = RDB[DATA_RESTART_READ_POINT]) == 0.0)
        Error(0, "Restart to zero burnup");
    }

  /* Reset closest or point if index is given */

  closest = INFTY;
  d0 = -1.0;

  /* NOTE: Tää oli asetettu nollaksi 2.1.29 ja aikaisemmat. Ei kuitenkaan */
  /* toimi niin vaan löytää väärän pisteen (2.7.2017 / 2.1.30 / JLe) */

  i = 1;

  /* Read loop */

  while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
    {
      /* Read name */

      if ((sz = fread(tmpstr, sizeof(char), n, fp)) == 0)
        Error(0, "Error in restart file");

      /* Put EOF */

      tmpstr[n] = '\0';

      /* Read nominal burnup and time */

      if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      /* Read number of nuclides, atomic and mass density and burnup */

      if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      /* Check values (mdens is not used) */

      CheckValue(FUNCTION_NAME, "nnuc", "", nnuc, 1, 10000);
      CheckValue(FUNCTION_NAME, "adens", "", adens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mdens", "", mdens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mbu", "", mbu, 0.0, 1000.0);

      /* Tää lisättiin 30.11.2017 / 2.1.30. Haku toimi aikaisemmin   */
      /* sattumalta kun filessä oleva ensimmäinen piste oli nollassa */
      /* ja d0 alustettiin nollaksi. Tämän pitäisi tehdä sama vaikka */
      /* eka piste on jotain muuta. (JLe) */

      if (d0 < 0.0)
        d0 = days;

      /* Check index was given */

      if ((idx > 0) && (idx == i))
        {
          /* Set point */

          pt = -days;
          closest = days;

          /* Break loop */

          break;
        }

      /* Check closest */

      if (pt > 0.0)
        {
          if (fabs(bu - pt) < fabs(closest - pt))
            closest = bu;
        }
      else if (pt < 0.0)
        {
          if (fabs(days + pt) < fabs(closest + pt))
            closest = days;
        }

      /* Skip nuclide data */

      n = (2*nnuc)*sizeof(double);
      fseek(fp, n, SEEK_CUR);

      /* Update counter */

      if (d0 != days)
        {
          i++;
          d0 = days;
        }
    }

  /* Check if continued from last */

  if ((long)RDB[DATA_RESTART_READ_CONTINUE] == YES)
    {
      /* Set points */

      pt = -days;
      closest = days;
      idx = i;
    }

  /* Rewind file */

  rewind(fp);

  /* Check point */

  if (((pt > 0.0) && (fabs(closest/pt - 1.0) > 0.001)) ||
      ((pt < 0.0) && (fabs(-closest/pt - 1.0) > 0.001)) ||
      ((idx > 0) && (idx != i)))
    {
      /* Point not found, print points */

      fprintf(stdout,
              "\nRestart file \"%s\" contains the following burnup points:\n\n",
              fname);

      /* Reset days */

      d0 = -1;
      i = 0;

      /* Read loop */

      while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
        {
          /* Read name */

          if ((sz = fread(tmpstr, sizeof(char), n, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");

          /* Put EOF */

          tmpstr[n] = '\0';

          /* Read data */

          if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");
          if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");
          if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");
          if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");
          if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");
          if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
            Die(FUNCTION_NAME, "Error in restart file");

          /* Print */

          if (days != d0)
            {
              /* Print burnup and time */

              fprintf(stdout, "%2ld : BU = %7.3f MWd/kgu, time = %1.5E days",
                      i++, bu, days);

              /* Identify closest point */

              if ((idx == 0) && (((pt > 0) && (bu == closest)) ||
                                 ((pt < 0) && (days == closest))))
                fprintf(stdout, " (closest point)\n");
              else
                fprintf(stdout, "\n");
            }

          /* Set days */

          d0 = days;

          /* Skip material data */

          n = (2*nnuc)*sizeof(double);
          fseek(fp, n, SEEK_CUR);
        }

      if (idx > 0)
        Error(0, "Burnup point %ld not found in restart file", idx - 1);
      else if (pt > 0.0)
        Error(0, "Burnup point %1.3f MWd/kgU not found in restart file", pt);
      else
        Error(0, "Burnup point %1.3f days not found in restart file", -pt);
    }
  else if ((pt > 0.0) && (fabs(closest/pt - 1.0) > 1E-6))
    Note(0, "No exact match in restart file, using %1.2f MWd/kgU", closest);
  else if ((pt < 0.0) && (fabs(-closest/pt - 1.0) > 1E-6))
    Note(0, "No exact match in restart file, using %s",
         TimeIntervalStr(closest*60*60*24));

  /* Put point to closest value */

  if (pt > 0)
    pt = closest;
  else
    pt = -closest;

  /***************************************************************************/

  /***** Read data ***********************************************************/

  /* Reset variables */

  ok = NO;
  days = 0;
  idx = 0;

  /* Read loop */

  while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
    {
      /* Remember days */

      days0 = days;

      /* Read name */

      if ((sz = fread(tmpstr, sizeof(char), n, fp)) == 0)
        Error(0, "Error in restart file");

      /* Put EOF */

      tmpstr[n] = '\0';

      /* Read nominal burnup and time */

      if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      /* Update index if new day */

      if (days != days0)
        idx++;

      /* Read number of nuclides, atomic and mass density and burnup */

      if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      /* Check values (mdens is not used) */

      CheckValue(FUNCTION_NAME, "nnuc", "", nnuc, 1, 10000);
      CheckValue(FUNCTION_NAME, "adens", "", adens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mdens", "", mdens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mbu", "", mbu, 0.0, 1000.0);

      /* Check point */

      if (((pt > 0.0) && (fabs(bu/pt - 1.0) > 1E-12)) ||
          ((pt < 0.0) && (fabs(-days/pt - 1.0) > 1E-12)) ||
          ((pt == 0.0) && (days != 0.0)))
        {
          /* Skip material */

          n = (2*nnuc)*sizeof(double);
          fseek(fp, n, SEEK_CUR);
        }
      else
        {
          /* Set flag */

          ok = YES;

          /* Check check mode */

          if (mode == RESTART_CHECK)
            {
              /* Skip material */

              n = (2*nnuc)*sizeof(double);
              fseek(fp, n, SEEK_CUR);

              /* Cycle loop */

              continue;
            }

          /* Check continue mode */

          if (((long)RDB[DATA_RESTART_READ_CONTINUE] == YES) &&
              ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR))
            {
              /* Put index */

              WDB[DATA_BURN_STEP] = idx;
            }

          /* Set burnup and irradiation time */

          WDB[DATA_BURN_CUM_BURNUP] = bu;
          WDB[DATA_BURN_CUM_BURNTIME0] = days0*24.0*60.0*60.0;
          WDB[DATA_BURN_CUM_BURNTIME] = days*24.0*60.0*60.0;

          /* Find material */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Compare */

              if (!strcmp(tmpstr, GetText(mat + MATERIAL_PTR_NAME)))
                break;

              /* Next material */

              mat = NextItem(mat);
            }

          /* Check mode */

          if (mode == RESTART_REPLACE)
            {
              /***************************************************************/

              /***** Replace entire material data ****************************/

              /* Check material pointer and create new if not found */

              if (mat < VALID_PTR)
                {
                  /* Check neutron transport mode */

                  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
                    Error(0, "No materials corresponding to restart file");
                  else
                    mat = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);
                }
              else if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
                       ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES))
                Error(mat,
                      "Material %s must be burnable to work with restart file",
                      GetText(mat + MATERIAL_PTR_NAME));

              /* Put parameter name and file name */

              WDB[mat + PARAM_PTR_NAME] = (double)PutText("mat");
              WDB[mat + PARAM_PTR_FNAME] = (double)PutText(fname);

              /* Put material name and atomic density */

              WDB[mat + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);
              WDB[mat + MATERIAL_ADENS] = adens;

              /* Reset pointer to composition */

              WDB[mat + MATERIAL_PTR_COMP] = -1;

              /* Loop over composition */

              for (n = 0; n < nnuc; n++)
                {
                  /* Read ZAI and atomic density */

                  if ((sz = fread(&zai, sizeof(long), 1, fp)) == 0)
                    Error(0, "Error in restart file");

                  if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
                    Error(0, "Error in restart file");

                  /* Check if Xenon is set to zero */

                  if ((long)RDB[DATA_RESTART_READ_ZERO_XE] == YES)
                    if ((zai == 541350) || (zai == 541351))
                      adens = 0.0;

                  /* Check if Samarium is set to zero */

                  if ((long)RDB[DATA_RESTART_READ_ZERO_SM] == YES)
                    if (zai == 621490)
                      adens = 0.0;

                  /* Skip lost at beginning */

                  if (n > 0)
                    {
                      /* Create new item */

                      iso = NewItem(mat + MATERIAL_PTR_COMP,
                                    COMPOSITION_BLOCK_SIZE);

                      /* Put name */

                      sprintf(ZAI, "%ld", zai);
                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        (double)PutText(ZAI);

                      /* Put atomic denisty */

                      WDB[iso + COMPOSITION_ADENS] = adens;
                    }
                }

              /***************************************************************/
            }
          else
            {
              /***************************************************************/

              /***** Replace material composition ****************************/

              /* Check pointer (JLe 10.6.2020 / 2.1.32: jostain syystä     */
              /* writedepfile.c kirjoittaa kompositiot div-kortilla        */
              /* jaetuille materiaaleille, vaikka burn-flagi ei olisikaan  */
              /* päällä. Data on parentille silloin nollaa, mikä aiheuttaa */
              /* errorin tuolla alempana). */

              if ((mat > VALID_PTR) &&
                  ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
                {
                  /* Put density (JLe 16.6.2017 / 2.1.30: tässä oli ennen    */
                  /* bugi joka aiheutti sen, että branch-laskussa käytettiin */
                  /* aina materiaalikortissa annettua tiheyttä. Skaalauksen  */
                  /* tarkoitus on ottaa branch-kortilla määritellyt tiheyden */
                  /* muutokset huomioon.) */

                  /* JLe 11.1.2019: ja toinen bugi, joka asettaa väärän    */
                  /* tiheyden sellaisissa restart-laskuissa, jotka alkavat */
                  /* muualta kuin nollapalamasta. Kertoimen asettaminen    */
                  /* siirrettiin invokebranch.c:hen */

                  if (RDB[mat + MATERIAL_RESTART_ADENS_F] > 0.0)
                    WDB[mat + MATERIAL_ADENS] =
                      RDB[mat + MATERIAL_RESTART_ADENS_F]*adens;
                  else
                    WDB[mat + MATERIAL_ADENS] = adens;

                  /* Check */

                  CheckValue(FUNCTION_NAME, "adens", "",
                             RDB[mat + MATERIAL_ADENS], 0.0, INFTY);

                  /* Put burnup */

                  WDB[mat + MATERIAL_BURNUP] = mbu;

                  /* Get pointer to composition list */

                  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
                      < VALID_PTR)
                    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

                  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

                  /* Allocate memory for temporary nuclide pointers */

                  n = ListSize(iso);
                  ptr = (long *)Mem(MEM_ALLOC, n, sizeof(double));

                  /* Store negative ZAI's in composition vector */

                  n = 0;

                  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
                      < VALID_PTR)
                    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

                  while (iso > VALID_PTR)
                    {
                      /* Pointer to nuclide */

                      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                      /* Store pointer */

                      ptr[n++] = nuc;

                      /* Replace pointer with ZAI */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        RDB[nuc + NUCLIDE_ZAI];

                      /* Reset density */

                      WDB[iso + COMPOSITION_ADENS] = 0.0;

                      /* Pointer to next */

                      iso = NextItem(iso);
                    }
                }
              else
                {
                  /* Not found, Skip material */

                  n = (2*nnuc)*sizeof(double);
                  fseek(fp, n, SEEK_CUR);

                  /* Reset nuclide count to avoid loop */

                  nnuc = -1;
                  iso = -1;
                  ptr = NULL;

                  /* Cycle loop */

                  continue;
                }

              /* Read data */

              for (n = 0; n < nnuc; n++)
                {
                  /* Read ZAI and atomic density */

                  if ((sz = fread(&zai, sizeof(long), 1, fp)) == 0)
                    Error(0, "Error in restart file");

                  if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
                    Error(0, "Error in restart file");

                  /* Check if Xenon is set to zero */

                  if ((long)RDB[DATA_RESTART_READ_ZERO_XE] == YES)
                    if ((zai == 541350) || (zai == 541351))
                      adens = 0.0;

                  /* Check if Samarium is set to zero */

                  if ((long)RDB[DATA_RESTART_READ_ZERO_SM] == YES)
                    if (zai == 621490)
                      adens = 0.0;

                  /* Pointer to composition vector */

                  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
                      < VALID_PTR)
                    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

                  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

                  /* Find nuclide */

                  if ((iso = SeekList(iso, COMPOSITION_PTR_NUCLIDE,
                                      (double)zai, SORT_MODE_ASCEND))
                      > VALID_PTR)
                    {
                      /* Put density */

                      WDB[iso + COMPOSITION_ADENS] = adens;
                    }
                  else if (zai > 0)
                    {
                      /* Tuo järjestetty haku feilaa jos pointterit ei    */
                      /* jostain syystä olekaan järjestyksessä. Yritetään */
                      /* uusiksi */

                      if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
                          < VALID_PTR)
                        iso = (long)RDB[mat + MATERIAL_PTR_COMP];

                      /* Search unsorted */

                      if ((iso = SeekList(iso, COMPOSITION_PTR_NUCLIDE,
                                          (double)zai, NO)) > VALID_PTR)
                        WDB[iso + COMPOSITION_ADENS] = adens;
                      else
                        Note(0, "Nuclide %ld not found in composition", zai);
                    }
                }

              /* Renormalize densities */

              sum = 0.0;

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Add to total */

                  sum = sum + RDB[iso + COMPOSITION_ADENS];

                  /* Next */

                  iso = NextItem(iso);
                }

              /* Check */

              if (sum == 0.0)
                Die(FUNCTION_NAME, "Zero density");

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Add to total */

                  WDB[iso + COMPOSITION_ADENS] =
                    RDB[iso + COMPOSITION_ADENS]*RDB[mat + MATERIAL_ADENS]/sum;

                  /* Next */

                  iso = NextItem(iso);
                }

              /* Restore pointers */

              if (ptr != NULL)
                {
                  /* Loop over composition */

                  n = 0;

                  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
                      < VALID_PTR)
                    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

                  while (iso > VALID_PTR)
                    {
                      /* Put pointer */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)ptr[n++];

                      /* Pointer to next */

                      iso = NextItem(iso);
                    }

                  /* Free memory */

                  Mem(MEM_FREE, ptr);
                }

              /* Calculate mass density */

              mdens = 0.0;
              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Pointer to nuclide */

                  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Add to total */

                  mdens = mdens + RDB[iso + COMPOSITION_ADENS]*
                    RDB[nuc + NUCLIDE_AW]/N_AVOGADRO;

                  /* Next */

                  iso = NextItem(iso);
                }

              /* Put value */

              WDB[mat + MATERIAL_MDENS] = mdens;

              /***************************************************************/
            }

          /*******************************************************************/
        }
    }

  /* Check ok flag */

  if (ok == NO)
    Die(FUNCTION_NAME, "WTF?");

  /***************************************************************************/

  /* Check continue mode */

  if (((long)RDB[DATA_RESTART_READ_CONTINUE] == YES) &&
      ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR))
    {
      /* Put total number of steps */

      WDB[DATA_BURN_TOT_STEPS] = RDB[DATA_BURN_TOT_STEPS]
        + RDB[DATA_BURN_STEP];
    }

  /* Store starting point */

  WDB[DATA_RESTART_START_POINT] = RDB[DATA_BURN_CUM_BURNTIME];

  /* Calculate burnups, etc */

  SumDivCompositions();

  /* Close file */

  fclose(fp);

  /* Print */

  if (mode != RESTART_CHECK)
    {
      fprintf(outp, "- %1.2f MWd/kgU burnup\n", RDB[DATA_BURN_CUM_BURNUP]);

      if (days < 1E4)
        fprintf(outp, "- %1.2f days irradiation time\n\n",
                RDB[DATA_BURN_CUM_BURNTIME]/(24.0*60.0*60.0));
      else
        fprintf(outp, "- %s irradiation time\n\n",
                TimeIntervalStr(RDB[DATA_BURN_CUM_BURNTIME]));
    }
  else
    fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
