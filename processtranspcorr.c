/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processtranspcorr.c                            */
/*                                                                           */
/* Created:       2016/06/11 (JLe)                                           */
/* Last modified: 2019/10/21 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes transport corrections for in-scattering            */
/*              approximation.                                               */
/*                                                                           */
/* Comments: - Tätä yritettiin toteuttaa siten että jaetuille materiaaleille */
/*             luodaan omat structuret, jne. Ongelmaksi tuli kuitenkin se,   */
/*             että niillä ei vielä tässä vaiheessa ole noita composition    */
/*             listejä. Eli linkkaus pitää tehdä jossain muualla. Nyt toi    */
/*             antaa errorin jos materiaali on jaettu tai burnable.          */
/*             (JLE / 28.12.2016 / 2.1.28)                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessTranspCorr:"

/*****************************************************************************/

void ProcessTranspCorr()
{
  long loc0, loc1, mat, mat0, gcu, ptr, ng, nt, n, ZAI, iso, nuc;
  FILE *fp;

  /* Check if corrections are defined */

  if ((long)RDB[DATA_PTR_TRC0] < VALID_PTR)
    return;

  /* Check if group constant generation is on */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /***************************************************************************/

  /***** Process data ********************************************************/

  /* Loop over corrections */

  loc0 = (long)RDB[DATA_PTR_TRC0];
  while (loc0 > VALID_PTR)
    {
      /***********************************************************************/

      /***** Link to material data *******************************************/

      /* Reset counter */

      n = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if not parent (tää liittyy tohon jakoon) */

          if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
            mat0 = mat;

          /* Compare */

          if (CompareStr(mat0 + MATERIAL_PTR_NAME, loc0 + TRANSP_CORR_PTR_MAT))
            {
              /* Link pointer */

              WDB[mat + MATERIAL_PTR_TRANSP_CORR] = (double)loc0;

              /* Nää tarkistukset estää errorin myöhemmin (composition */
              /* listaa ei oo luotu) */

              if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
                Error(loc0, "Method does not work with divided materials");

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                Error(loc0, "Method does not work with burnable materials");

              /* Increase counter */

              n++;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check count */

      if (n == 0)
        Error(loc0,
              "Material %s assigned with transport correction not defined",
              GetText(loc0 + TRANSP_CORR_PTR_MAT));

      /***********************************************************************/

      /****** Read data ******************************************************/

      /* Test format */

      WDB[DATA_DUMMY] = RDB[loc0 + TRANSP_CORR_FNAME];
      TestDOSFile(GetText(DATA_DUMMY));

      /* Open file for reading */

      fp = OpenDataFile(loc0 + TRANSP_CORR_FNAME, "transport correction data");

      /* Read number of energy groups */

      if (fscanf(fp, "%ld", &ng) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      /* Put value */

      WDB[loc0 + TRANSP_CORR_NG] = (double)ng;

      /* Allocate memory for boundaries */

      ptr = ReallocMem(DATA_ARRAY, ng + 1);
      WDB[loc0 + TRANSP_CORR_PTR_ENE] = (double)ptr;

      /* Read data */

      for (n = 0; n < ng + 1; n++)
        if (fscanf(fp, "%lf", &WDB[ptr + n]) == EOF)
          Die(FUNCTION_NAME, "fscanf error");

      /* Check negative points */

      for (n = 0; n < ng + 1; n++)
        if (RDB[ptr + n] < 0.0)
          Error(loc0, "Negative energy in file %s",
                GetText(loc0 + TRANSP_CORR_FNAME));

      /* Check order */

      for (n = 1; n < ng + 1; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Error(loc0, "Energies in file %s not in ascending order",
                GetText(loc0 + TRANSP_CORR_FNAME));

      /* Read number of temperatures */

      if (fscanf(fp, "%ld", &nt) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      /* Put value */

      WDB[loc0 + TRANSP_CORR_NT] = (double)nt;

      /* Check if multiple temperatures */

      if (nt > 1)
        {
          /* Allocate memory for temperatures */

          ptr = ReallocMem(DATA_ARRAY, nt);
          WDB[loc0 + TRANSP_CORR_PTR_TEMP] = (double)ptr;

          /* Read data */

          for (n = 0; n < nt; n++)
            if (fscanf(fp, "%lf", &WDB[ptr + n]) == EOF)
              Die(FUNCTION_NAME, "fscanf error");

          /* Check order */

          for (n = 1; n < nt; n++)
            if (RDB[ptr + n] <= RDB[ptr + n - 1])
              Error(loc0, "Temperatures in file %s not in ascending order",
                    GetText(loc0 + TRANSP_CORR_FNAME));
        }
      else
        nt = 1;

      /* Allocate memory for data */

      ptr = ReallocMem(DATA_ARRAY, ng*nt);
      WDB[loc0 + TRANSP_CORR_PTR_DATA] = (double)ptr;

      /* Read data */

      for (n = 0; n < ng*nt; n++)
        if (fscanf(fp, "%lf", &WDB[ptr + n]) == EOF)
          Die(FUNCTION_NAME, "fscanf error");

      /* Check negative points */

      for (n = 0; n < ng*nt; n++)
        if (RDB[ptr + n] < 0.0)
          Error(loc0, "Negative point in file %s",
                GetText(loc0 + TRANSP_CORR_FNAME));

      /* Check values */

      for (n = 0; n < ng*nt; n++)
        if ((RDB[ptr + n] < 0.0) || (RDB[ptr + n] > 100.0))
          Error(loc0, "Invalid value %E in transport correction curve",
                RDB[ptr + n]);

      /* Close file */

      fclose(fp);

      /***********************************************************************/

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Link to nuclide data ************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to transport correction */

      if ((loc0 = (long)RDB[mat + MATERIAL_PTR_TRANSP_CORR]) < VALID_PTR)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Count size of nuclide list */

      n = 0;
      if ((ptr = (long)RDB[loc0 + TRANSP_CORR_PTR_ISO]) > VALID_PTR)
        while ((ZAI = (long)RDB[ptr++]) > 0)
          n++;

      /* Check */

      if (n == 0)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Create duplicate */

      loc0 = DuplicateItem(loc0);
      WDB[mat + MATERIAL_PTR_TRANSP_CORR] = (double)loc0;

      /* Allocate memory for new list */

      loc1 = ReallocMem(DATA_ARRAY, n + 1);

      /* Pointer to old list */

      ptr = (long)RDB[loc0 + TRANSP_CORR_PTR_ISO];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Put new pointer */

      WDB[loc0 + TRANSP_CORR_PTR_ISO] = (double)loc1;

      /* Loop over ZAI */

      while ((ZAI = (long)RDB[ptr]) > 0)
        {
          /* Loop over composition and match ZAI */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Compare ZAI */

              if (ZAI == (long)RDB[nuc + NUCLIDE_ZAI])
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check if found and store pointer */

          if (iso < VALID_PTR)
            Error(loc0, "nuclide %ld not found in material %s", ZAI,
                  GetText(mat + MATERIAL_PTR_NAME));
          else
            WDB[loc1++] = (double)iso;

          /* Next */

          ptr++;
        }

      /* Add null terminator */

      WDB[loc1] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Allocate memory for additional cross sections ***********************/

  /* Get pointer to micro-group structure */

  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr0)", DATA_ARRAY, ptr);

  /* Number of groups */

  ng = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

  /* Get pointer to macro-group structure */

  ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of groups */

  n = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

  /* Loop over gcu structures */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Allocate memory for micro-group data */

      ptr = AllocPrivateData(ng, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_TRC] = (double)ptr;

      ptr = AllocPrivateData(ng, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_TRC_FLX] = (double)ptr;

      /* Allocate memory for results */

      ptr = NewStat("TRC_TRANSPXS", 1, n);
      WDB[gcu + GCU_TRC_TRANSPXS] = (double)ptr;

      ptr = NewStat("TRC_DIFFCOEF", 1, n);
      WDB[gcu + GCU_TRC_DIFFCOEF] = (double)ptr;

      /* Remember last value for coefoutput.c */

      WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

      /* Next */

      gcu = NextItem(gcu);
    }

  /***************************************************************************/
}

/*****************************************************************************/
