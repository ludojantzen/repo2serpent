/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writedepfile.c                                 */
/*                                                                           */
/* Created:       2011/05/24 (JLe)                                           */
/* Last modified: 2020/03/07 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Writes binary file containing data from burnup calculation   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteDepFile:"

/*****************************************************************************/

void WriteDepFile()
{
  long nuc, mat, nnuc, nmat, iso, n, m, count, rad;
  char tmpstr[MAX_STR];
  double adens, val, bu, days;
  FILE *fp;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /***************************************************************************/

  /***** Write binary work file for restarts *********************************/

  /* Get burnup and days */

  bu = RDB[DATA_BURN_CUM_BURNUP];
  days = RDB[DATA_BURN_CUM_BURNTIME]/24.0/60.0/60.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Write data */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) ||
          ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT))
        StoreComposition(mat, bu, days);

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Write binary data for depletion output ******************************/

  /* File name */

  sprintf(tmpstr, "%s.dep", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for writing */

  if (((long)RDB[DATA_BURN_STEP] == 0) &&
      (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    fp = fopen(tmpstr, "w");
  else
    fp = fopen(tmpstr, "a");

  /* Check pointer */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Get number of nuclides */

  nnuc = (long)RDB[DATA_TOT_NUCLIDES];

  /* Calculate number of materials */

  nmat = 0;

  mat = (long)RDB[DATA_PTR_M0];
  while(mat > VALID_PTR)
    {
      /* Check output flag */

      if ((long)RDB[mat + MATERIAL_BURN_PRINT_OUTPUT] == YES)
        nmat++;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Check first step */

  if (((long)RDB[DATA_BURN_STEP] == 0) &&
      (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    {
      /***********************************************************************/

      /***** Nuclide data ****************************************************/

      /* Write number of materials */

      fwrite(&nmat, sizeof(long), 1, fp);

      /* Write number of nuclides */

      fwrite(&nnuc, sizeof(long), 1, fp);

      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Write ZAI */

          n = (long)RDB[nuc + NUCLIDE_ZAI];
          fwrite(&n, sizeof(long), 1, fp);

          /* Write Z */

          n = (long)RDB[nuc + NUCLIDE_Z];
          fwrite(&n, sizeof(long), 1, fp);

          /* Write atomic weight */

          val = RDB[nuc + NUCLIDE_AW];
          fwrite(&val, sizeof(double), 1, fp);

          /* Decay constant */

          val = RDB[nuc + NUCLIDE_LAMBDA];
          fwrite(&val, sizeof(double), 1, fp);

          /* Decay heat */

          val = RDB[nuc + NUCLIDE_DECAY_E]*RDB[nuc + NUCLIDE_LAMBDA];
          fwrite(&val, sizeof(double), 1, fp);

          /* Spontaneous fission rate */

          val = RDB[nuc + NUCLIDE_SF_BR]*RDB[nuc + NUCLIDE_LAMBDA];
          fwrite(&val, sizeof(double), 1, fp);

          /* Reset photon emission rate */

          val = 0.0;

          /* Loop over radiations */

          rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
          while (rad > VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[rad + NUCLIDE_RAD_TYPE] == PARTICLE_TYPE_GAMMA)
                {
                  /* Photon emission rate */

                  val = RDB[rad + NUCLIDE_RAD_SPEC_I]*
                    RDB[nuc + NUCLIDE_LAMBDA];

                  /* break loop */

                  break;
                }

              /* Next */

              rad = NextItem(rad);
            }

          /* Write photon emission rate */

          fwrite(&val, sizeof(double), 1, fp);

          /* Specific ingestion toxicity */

          if ((val = RDB[nuc + NUCLIDE_SPEC_ING_TOX]) < 0.0)
            val = 0.0;

          fwrite(&val, sizeof(double), 1, fp);

          /* Specific inhalation toxicity */

          if ((val = RDB[nuc + NUCLIDE_SPEC_INH_TOX]) < 0.0)
            val = 0.0;

          fwrite(&val, sizeof(double), 1, fp);

          /* Next nuclide */

          nuc = NextItem(nuc);
        }

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Reset some nuclide-wise quantities **********************************/

  /* NOTE: These are used for top lists */

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR )
    {
      WDB[nuc + NUCLIDE_TOP_MASS] = 0.0;
      WDB[nuc + NUCLIDE_TOP_ACTIVITY] = 0.0;
      WDB[nuc + NUCLIDE_TOP_SF] = 0.0;
      WDB[nuc + NUCLIDE_TOP_GSRC] = 0.0;
      WDB[nuc + NUCLIDE_TOP_DECAY_HEAT] = 0.0;
      WDB[nuc + NUCLIDE_TOP_ING_TOX] = 0.0;
      WDB[nuc + NUCLIDE_TOP_INH_TOX] = 0.0;

      /* Next */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Material-wise data **************************************************/

  /* Reset index */

  m = 0;

  mat = (long)RDB[DATA_PTR_M0];
  while(mat > VALID_PTR)
    {
      /* Check output flag */

      if ((long)RDB[mat + MATERIAL_BURN_PRINT_OUTPUT] == NO)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check domain decomposition */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
          ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT))
        Die(FUNCTION_NAME, "Trying to write divide compositions with DD");

      /* Calculate number of nuclides */

      count = 0;

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Get atomic density */

          adens = RDB[iso + COMPOSITION_ADENS];
          CheckValue(FUNCTION_NAME, "adens", "", adens, 0.0, INFTY);

          /* Check that value is non-zero */

          if (adens > ZERO)
            count++;

          /* Next nuclide in composition */

          iso = NextItem(iso);
        }

      /* Write burnup (nää vie tilaa mutta yksinkertaistaa lukemista) */

      val = RDB[DATA_BURN_CUM_BURNUP];
      fwrite(&val, sizeof(double), 1, fp);

      val = RDB[DATA_BURN_CUM_REAL_BURNUP];
      fwrite(&val, sizeof(double), 1, fp);

      /* Write burn time */

      val = RDB[DATA_BURN_CUM_BURNTIME]/24.0/60.0/60.0;
      fwrite(&val, sizeof(double), 1, fp);

      /* Write number of nuclides */

      fwrite(&count, sizeof(long), 1, fp);

      /* Write material volume */

      val = RDB[mat + MATERIAL_VOLUME];
      fwrite(&val, sizeof(double), 1, fp);

      /* Write material burnup */

      val = RDB[mat + MATERIAL_BURNUP];
      fwrite(&val, sizeof(double), 1, fp);

      /* Write length of material name and name */

      sprintf(tmpstr, "%s", GetText(mat + MATERIAL_PTR_NAME));
      n = strlen(tmpstr);

      fwrite(&n, sizeof(long), 1, fp);
      fwrite(tmpstr, sizeof(char), n, fp);

      /* Write index */

      fwrite(&m, sizeof(long), 1, fp);

      /* Write burnup step */

      n = (long)RDB[DATA_BURN_STEP];

      /* In CI calculations burn step is a bit different */
      /*
      if(RDB[DATA_BURN_SIE] == (double)YES)
        {
          n = (long)RDB[DATA_BURN_STEP_TOT];
        }
      */
      fwrite(&n, sizeof(long), 1, fp);

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

          /* Nuclide index */

          n = (long)RDB[nuc + NUCLIDE_IDX];

          /* Atomic density */

          adens = RDB[iso + COMPOSITION_ADENS];
          CheckValue(FUNCTION_NAME, "adens", "", adens, 0.0, INFTY);

          /* Check that value is non-zero */

          if (adens > ZERO)
            {
              /* Write index */

              fwrite(&n, sizeof(long), 1, fp);

              /* Write composition */

              fwrite(&adens, sizeof(double), 1, fp);

              /* Convert atomic density to number of nuclides */

              adens = adens*RDB[mat + MATERIAL_VOLUME];

              /* Add to totals */

              WDB[nuc + NUCLIDE_TOP_MASS] = RDB[nuc + NUCLIDE_TOP_MASS]
                + adens*RDB[nuc + NUCLIDE_AW]/N_AVOGADRO;

              WDB[nuc + NUCLIDE_TOP_ACTIVITY] = RDB[nuc + NUCLIDE_TOP_ACTIVITY]
                + adens*RDB[nuc + NUCLIDE_LAMBDA]/BARN;

              WDB[nuc + NUCLIDE_TOP_SF] = RDB[nuc + NUCLIDE_TOP_SF]
                + adens*RDB[nuc + NUCLIDE_SF_BR]*RDB[nuc + NUCLIDE_LAMBDA]/BARN;

              WDB[nuc + NUCLIDE_TOP_DECAY_HEAT] =
                WDB[nuc + NUCLIDE_TOP_DECAY_HEAT]
                + adens*RDB[nuc + NUCLIDE_DECAY_E]*
                RDB[nuc + NUCLIDE_LAMBDA]*MEV/BARN;

              WDB[nuc + NUCLIDE_TOP_ING_TOX] = RDB[nuc + NUCLIDE_TOP_ING_TOX]
                + adens*RDB[nuc + NUCLIDE_SPEC_ING_TOX]*
                RDB[nuc + NUCLIDE_LAMBDA]/BARN;

              WDB[nuc + NUCLIDE_TOP_INH_TOX] =
                WDB[nuc + NUCLIDE_TOP_INH_TOX]
                + adens*RDB[nuc + NUCLIDE_SPEC_INH_TOX]*
                RDB[nuc + NUCLIDE_LAMBDA]/BARN;

              /* Loop over radiations */

              rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
              while (rad > VALID_PTR)
                {
                  /* Check type */

                  if ((long)RDB[rad + NUCLIDE_RAD_TYPE] == PARTICLE_TYPE_GAMMA)
                    {
                      /* Photon emission rate */

                      WDB[nuc + NUCLIDE_TOP_GSRC] = RDB[nuc + NUCLIDE_TOP_GSRC]
                        + adens*RDB[rad + NUCLIDE_RAD_SPEC_I]*
                        RDB[nuc + NUCLIDE_LAMBDA]/BARN;

                      /* break loop */

                      break;
                    }

                  /* Next */

                  rad = NextItem(rad);
                }
            }

          /* Next nuclide in composition */

          iso = NextItem(iso);
        }

      /* Update index */

      m++;

      /* Write burnup step */

      n = (long)RDB[DATA_BURN_STEP];

      /* In CI calculations burn step is a bit different */
      /*
      if(RDB[DATA_BURN_SIE] == (double)YES)
        {
          n = (long)RDB[DATA_BURN_STEP_TOT];
        }
      */
      fwrite(&n, sizeof(long), 1, fp);

      /* Next material */

      mat = NextItem(mat);
    }

  /* Check index */

  if (m == 0)
    Die(FUNCTION_NAME, "No material data written");

  /* Close file */

  fclose(fp);

  /***************************************************************************/

  /***** Print material-wise burnups in a separate file **********************/

  if ((1 == 2) && ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP))
    {
      /* Open file for writing */

      sprintf(tmpstr, "%s_bu%ld.m", GetText(DATA_PTR_INPUT_FNAME),
              (long)RDB[DATA_BURN_STEP]);

      if ((fp = fopen(tmpstr, "w"))  == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Print data */

      fprintf(fp, "dat = [\n");

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if sub-zone */

          if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
            {
              /* Calculate FIMA */

              if (RDB[mat + MATERIAL_FIMA_ADENS0] > 0.0)
                val = 1.0 - RDB[mat + MATERIAL_FIMA_ADENS]/
                  RDB[mat + MATERIAL_FIMA_ADENS0];
              else
                val = 0.0;

              /* Write mass, burnup and FIMA */

              fprintf(fp, "%1.5E %1.5E %1.5E %1.5E %% %s\n",
                      RDB[mat + MATERIAL_MASS],
                      RDB[mat + MATERIAL_INI_FMASS],
                      RDB[mat + MATERIAL_BURNUP], val,
                      GetText(mat + MATERIAL_PTR_NAME));
            }

          /* Next */

          mat = NextItem(mat);
        }

      fprintf(fp, "];\n");

      /* Close file */

      fclose(fp);
    }

  /***************************************************************************/
}

/*****************************************************************************/
