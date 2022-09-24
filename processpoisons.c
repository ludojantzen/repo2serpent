/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processpoisons.c                               */
/*                                                                           */
/* Created:       2012/12/03 (JLe)                                           */
/* Last modified: 2020/04/16 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Adds data for equilibriun Xe-135 calculation and poison      */
/*              production                                                   */
/*                                                                           */
/* Comments: - Tässä on nyt ihan järjettömän monimutkainen viritelmä         */
/*             tekemässä suht. yksinkertaista asiaa.                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPoisons:"

/*****************************************************************************/

void ProcessPoisons()
{
  long mat, mat0, iso, nuc, ptr, add, rea, yld, ZA, I, TMS;
  double f;

  /* Check poison calculation flag */

  if (((long)RDB[DATA_OPTI_POISON_CALC] == NO) &&
      ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] != YES) &&
      ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] != YES) &&
      ((long)RDB[DATA_PTR_XENON_MAT_LIST] < VALID_PTR) &&
      ((long)RDB[DATA_PTR_SAMARIUM_MAT_LIST] < VALID_PTR))
    return;

  /* Check fission yield library */

  if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
    if ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] < VALID_PTR)
      Error(0, "Fission yield library is required for poison cross sections");

  if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] > -1)
    if ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] < VALID_PTR)
      Error(0, "Fission yield library is required for equilibirum Xe-135 calculation");

  if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] > -1)
    if ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] < VALID_PTR)
      Error(0, "Fission yield library is required for equilibirum Sm-149 calculation");

  /* Check decay data library */

  if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] > -1)
    if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] < VALID_PTR)
      Error(0, "Decay library is required for equilibirum Xe-135 calculation");

  if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] > -1)
    if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] < VALID_PTR)
      Error(0, "Decay library is required for equilibirum Sm-149 calculation");

  fprintf(outp, "Setting fission product poison data...\n\n");

  /***************************************************************************/

  /***** Flag materials for equilibrium Xe-135 calculation *******************/

  /* Get flag */

  if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] < 0)
    add = NO;
  else if ((long)RDB[DATA_PTR_XENON_MAT_LIST] < VALID_PTR)
    add = YES;
  else
    add = !((long)RDB[DATA_XENON_EQUILIBRIUM_MODE]);

  /* Loop over materials and set or reset flags */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if divided */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        {
          /* For divided materials the parent is fissile and the divided is */
          /* physical at this point (JLe 30.5.2015 / 2.1.24) */

          if (((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
              ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
            WDB[mat0 + MATERIAL_XENON_EQUIL_CALC] = (double)add;
        }
      else if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
               ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        WDB[mat + MATERIAL_XENON_EQUIL_CALC] = (double)add;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Loop over listed materials and set flags */

  if ((ptr = (long)RDB[DATA_PTR_XENON_MAT_LIST]) > VALID_PTR)
    while ((long)RDB[ptr] > VALID_PTR)
      {
        /* Find material */

        mat = (long)RDB[DATA_PTR_M0];
        while (mat > VALID_PTR)
          {
            /* Check if material was produced by division */

            if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
              mat0 = mat;

            /* Check */

            if (CompareStr(mat0 + MATERIAL_PTR_NAME, ptr))
              {
                /* Check pointer */

                if (mat < VALID_PTR)
                  Error(0, "Material %s not defined", GetText(ptr));
                else if (!((long)RDB[mat0 + MATERIAL_OPTIONS] &
                           OPT_FISSILE_MAT))
                  Error(0, "Material %s is not fissile", GetText(ptr));
                else
                  WDB[mat0 + MATERIAL_XENON_EQUIL_CALC] =
                    RDB[DATA_XENON_EQUILIBRIUM_MODE];
              }

            /* Next material */

            mat = NextItem(mat);
          }

        /* Next */

        ptr++;
      }

  /***************************************************************************/

  /***** Flag materials for equilibrium Sm-149 calculation *******************/

  /* Get flag */

  if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] < 0)
    add = NO;
  else if ((long)RDB[DATA_PTR_SAMARIUM_MAT_LIST] < VALID_PTR)
    add = YES;
  else
    add = !((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE]);

  /* Loop over materials and set or reset flags */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if divided */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        {
          /* For divided materials the parent is fissile and the divided is */
          /* physical at this point (JLe 30.5.2015 / 2.1.24) */

          if (((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
              ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
            WDB[mat0 + MATERIAL_SAMARIUM_EQUIL_CALC] = (double)add;
        }
      else if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
               ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        WDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] = (double)add;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Loop over listed materials and set flags */

  if ((ptr = (long)RDB[DATA_PTR_SAMARIUM_MAT_LIST]) > VALID_PTR)
    while ((long)RDB[ptr] > VALID_PTR)
      {
        /* Find material */

        mat = (long)RDB[DATA_PTR_M0];
        while (mat > VALID_PTR)
          {
            /* Check if material was produced by division */

            if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
              mat0 = mat;

            /* Check */

            if (CompareStr(mat0 + MATERIAL_PTR_NAME, ptr))
              {
                /* Check pointer */

                if (mat < VALID_PTR)
                  Error(0, "Material %s not defined", GetText(ptr));
                else if (!((long)RDB[mat0 + MATERIAL_OPTIONS] &
                           OPT_FISSILE_MAT))
                  Error(0, "Material %s is not fissile", GetText(ptr));
                else
                  WDB[mat0 + MATERIAL_SAMARIUM_EQUIL_CALC] =
                    RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE];
              }

            /* Next material */

            mat = NextItem(mat);
          }

        /* Next */

        ptr++;
      }

  /***************************************************************************/

  /***** Add I-135 data ******************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if (((long)RDB[mat0 + MATERIAL_XENON_EQUIL_CALC] == NO) &&
          ((long)RDB[DATA_OPTI_POISON_CALC] == NO))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing I-135 data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 531350)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_I135_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 531350,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find I-135 data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_I135_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Add Xe-135 data *****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if (((long)RDB[mat0 + MATERIAL_XENON_EQUIL_CALC] == NO) &&
          ((long)RDB[DATA_OPTI_POISON_CALC] == NO))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing Xe-135 data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 541350)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_XE135_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 541350,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find Xe-135 data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_XE135_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Find nuclide in majorant extra list */

      ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
      while (ptr > VALID_PTR)
        {
          /* Check nuclide pointer */

          if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] == nuc)
            break;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Check if found */

      if (ptr < VALID_PTR)
        {
          /* New entry */

          ptr = NewItem(DATA_MAJORANT_PTR_EXTRA_XS, MAJORANT_EXTRA_BLOCK_SIZE);

          /* Put pointer and type */

          WDB[ptr + MAJORANT_EXTRA_PTR_NUC] = (double)nuc;
          WDB[ptr + MAJORANT_EXTRA_TYPE] =
            (double)MAJORANT_EXTRA_FP_POISON_ITER;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Add Xe-135m data ****************************************************/

  /* Check mode */

  if ((long)RDB[DATA_OPTI_POISON_CALC_XE135M] == YES)
    {
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* For divided materials the parent is fissile and the divided is */
          /* physical at this point (JLe 30.5.2015 / 2.1.24) */

          if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
            mat0 = mat;

          /* Check if physical and fissile */

          if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
              !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Check if added */

          if (((long)RDB[mat0 + MATERIAL_XENON_EQUIL_CALC] == NO) &&
              ((long)RDB[DATA_OPTI_POISON_CALC] == NO))
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Find existing Xe-135m data */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check ZAI */

              if ((long)RDB[nuc + NUCLIDE_ZAI] == 541351)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso > VALID_PTR)
            {
              /* Set pointer */

              WDB[mat0 + MATERIAL_PTR_XE135M_ISO] = (double)iso;
            }
          else
            {
              /* Avoid compiler warning */

              nuc = -1.0;

              /* Find first fissile isotope */

              iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Pointer to nuclide data */

                  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Check flag */

                  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                      NUCLIDE_FLAG_FISSILE)
                    break;

                  /* Next */

                  iso = NextItem(iso);
                }

              /* Check pointer */

              if (iso < VALID_PTR)
                Die(FUNCTION_NAME, "No fissile isotopes");

              /* Get TMS flag */

              TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

              /* Add nuclide */

              if ((nuc = AddNuclide(NULL, 541351,
                                    GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                    RDB[nuc + NUCLIDE_TEMP],
                                    NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
                Error(0,
                      "Unable to find Xe-135m data in cross section library");

              /* Add nuclide to composition */

              iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

              /* Put nuclide pointer */

              WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

              /* Set initial type flag */

              SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

              /* Set level */

              WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

              /* Set pointer */

              WDB[mat0 + MATERIAL_PTR_XE135M_ISO] = (double)iso;
            }

          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Set TMS temperatures */

          if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
            WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

          if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
            WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

          /* Find nuclide in majorant extra list */

          ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
          while (ptr > VALID_PTR)
            {
              /* Check nuclide pointer */

              if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] == nuc)
                break;

              /* Next */

              ptr = NextItem(ptr);
            }

          /* Check if found */

          if (ptr < VALID_PTR)
            {
              /* New entry */

              ptr = NewItem(DATA_MAJORANT_PTR_EXTRA_XS,
                            MAJORANT_EXTRA_BLOCK_SIZE);

              /* Put pointer and type */

              WDB[ptr + MAJORANT_EXTRA_PTR_NUC] = (double)nuc;
              WDB[ptr + MAJORANT_EXTRA_TYPE] =
                (double)MAJORANT_EXTRA_FP_POISON_ITER;
            }

          /* Next material */

          mat = NextItem(mat);
        }
    }

  /***************************************************************************/

  /***** Add Pm-147 data *****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == NO)
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing Pm-147 data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 611470)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM147_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 611470,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find Pm-147 data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM147_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Add Pm-148 data *****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == NO)
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing Pm-148 data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 611480)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM148_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 611480,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find Pm-148 data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM148_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Add Pm-148m data ****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == NO)
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing Pm-148m data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 611481)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM148M_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 611481,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find Pm-148m data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM148M_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Add Pm-149 data *****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if (((long)RDB[mat0 + MATERIAL_SAMARIUM_EQUIL_CALC] == NO) &&
          ((long)RDB[DATA_OPTI_POISON_CALC] == NO))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing Pm-149 data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 611490)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM149_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 611490,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find Pm-149 data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_PM149_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Add Sm-149 data *****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* For divided materials the parent is fissile and the divided is */
      /* physical at this point (JLe 30.5.2015 / 2.1.24) */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check if physical and fissile */

      if (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if added */

      if (((long)RDB[mat0 + MATERIAL_SAMARIUM_EQUIL_CALC] == NO) &&
          ((long)RDB[DATA_OPTI_POISON_CALC] == NO))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Find existing Sm-149 data */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check ZAI */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == 621490)
            break;

          /* Next */

          iso = NextItem(iso);
        }

      /* Check pointer */

      if (iso > VALID_PTR)
        {
          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_SM149_ISO] = (double)iso;
        }
      else
        {
          /* Avoid compiler warning */

          nuc = -1.0;

          /* Find first fissile isotope */

          iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check flag */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
                break;

              /* Next */

              iso = NextItem(iso);
            }

          /* Check pointer */

          if (iso < VALID_PTR)
            Die(FUNCTION_NAME, "No fissile isotopes");

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Add nuclide */

          if ((nuc = AddNuclide(NULL, 621490,
                                GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                RDB[nuc + NUCLIDE_TEMP],
                                NUCLIDE_TYPE_TRANSPORT, TMS)) < VALID_PTR)
            Error(0, "Unable to find Sm-149 data in cross section library");

          /* Add nuclide to composition */

          iso = NewItem(mat0 + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

          /* Set initial type flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

          /* Set level */

          WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;

          /* Set pointer */

          WDB[mat0 + MATERIAL_PTR_SM149_ISO] = (double)iso;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set TMS temperatures */

      if (RDB[mat0 + MATERIAL_TMS_TMAX] > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[mat0 + MATERIAL_TMS_TMAX];

      if (RDB[mat0 + MATERIAL_TMS_TMIN] < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[mat0 + MATERIAL_TMS_TMIN];

      /* Find nuclide in majorant extra list */

      ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
      while (ptr > VALID_PTR)
        {
          /* Check nuclide pointer */

          if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] == nuc)
            break;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Check if found */

      if (ptr < VALID_PTR)
        {
          /* New entry */

          ptr = NewItem(DATA_MAJORANT_PTR_EXTRA_XS, MAJORANT_EXTRA_BLOCK_SIZE);

          /* Put pointer and type */

          WDB[ptr + MAJORANT_EXTRA_PTR_NUC] = (double)nuc;
          WDB[ptr + MAJORANT_EXTRA_TYPE] =
            (double)MAJORANT_EXTRA_FP_POISON_ITER;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Find fission yields yields ******************************************/

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check fissile flag */

      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
        {
          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Pointer to fission yield */

              if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
                {
                  /* Reset yields */

                  WDB[rea + REACTION_XE135_YIELD] = 0.0;
                  WDB[rea + REACTION_XE135M_YIELD] = 0.0;
                  WDB[rea + REACTION_I135_YIELD] = 0.0;
                  WDB[rea + REACTION_SM149_YIELD] = 0.0;
                  WDB[rea + REACTION_PM147_YIELD] = 0.0;
                  WDB[rea + REACTION_PM148_YIELD] = 0.0;
                  WDB[rea + REACTION_PM148M_YIELD] = 0.0;
                  WDB[rea + REACTION_PM149_YIELD] = 0.0;

                  /* Pointer to distribution */

                  ptr = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Loop over distribution */

                  while (ptr > VALID_PTR)
                    {
                      /* Get ZA */

                      ZA = (long)(RDB[ptr + FY_TGT_ZAI]/10.0);
                      I = (long)RDB[ptr + FY_TGT_ZAI] - 10*ZA;

                      /* Independent yields for Xe-135 and Sm-149 */

                      f = RDB[ptr + FY_INDEPENDENT_FRAC];

                      if (ZA == 54135)
                        {
                          /* Check mode */

                          if ((long)RDB[DATA_OPTI_POISON_CALC_XE135M] == NO)
                            WDB[rea + REACTION_XE135_YIELD] =
                              RDB[rea + REACTION_XE135_YIELD] + f;
                          else if (I == 0)
                            WDB[rea + REACTION_XE135_YIELD] =
                              RDB[rea + REACTION_XE135_YIELD] + f;
                          else
                            WDB[rea + REACTION_XE135M_YIELD] =
                              RDB[rea + REACTION_XE135M_YIELD] + f;
                        }

                      if (ZA == 62149)
                        WDB[rea + REACTION_SM149_YIELD] =
                          RDB[rea + REACTION_SM149_YIELD] + f;

                      /* Independent yields for Pm-147, Pm-148 and Pm-148m */

                      if (ZA == 61147)
                        WDB[rea + REACTION_PM147_YIELD] =
                          RDB[rea + REACTION_PM147_YIELD] + f;

                      if ((ZA == 61148) && (I == 0))
                        WDB[rea + REACTION_PM148_YIELD] =
                          RDB[rea + REACTION_PM148_YIELD] + f;

                      if ((ZA == 61148) && (I > 0))
                        WDB[rea + REACTION_PM148M_YIELD] =
                          RDB[rea + REACTION_PM148M_YIELD] + f;

                      /* Cumulative yields for I-135 and Pm-149 */

                      f = RDB[ptr + FY_CUMULATIVE_FRAC];

                      if (ZA == 53135)
                        WDB[rea + REACTION_I135_YIELD] =
                          RDB[rea + REACTION_I135_YIELD] + f;

                      if (ZA == 61149)
                        WDB[rea + REACTION_PM149_YIELD] =
                          RDB[rea + REACTION_PM149_YIELD] + f;

                      /* Next */

                      ptr = NextItem(ptr);
                    }
                }

              /* Next reaction */

              rea = NextItem(rea);
            }
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Allocate memory for statistics **************************************/

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material was produced by division */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        mat0 = mat;

      /* Check fissile flag */

      if (((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
          ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          ptr = NewStat("MAT_I135_PROD_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_I135_PROD_RATE] = (double)ptr;

          ptr = NewStat("MAT_XE135_PROD_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_XE135_PROD_RATE] = (double)ptr;

          ptr = NewStat("MAT_PM149_PROD_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_PM149_PROD_RATE] = (double)ptr;

          ptr = NewStat("MAT_SM149_PROD_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_SM149_PROD_RATE] = (double)ptr;

          ptr = NewStat("MAT_I135_ABS_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_I135_ABS_RATE] = (double)ptr;

          ptr = NewStat("MAT_XE135_ABS_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_XE135_ABS_RATE] = (double)ptr;

          ptr = NewStat("MAT_PM149_ABS_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_PM149_ABS_RATE] = (double)ptr;

          ptr = NewStat("MAT_SM149_ABS_RATE", 1, 1);
          WDB[mat + MATERIAL_PTR_SM149_ABS_RATE] = (double)ptr;
        }

      /* Check equilibrium Xe-135 calculation */

      if ((long)RDB[mat0 + MATERIAL_XENON_EQUIL_CALC] == YES)
        {
          ptr = NewStat("MAT_I135_CONC", 1, 1);
          WDB[mat + MATERIAL_PTR_I135_CONC] = (double)ptr;

          ptr = NewStat("MAT_XE135_CONC", 1, 1);
          WDB[mat + MATERIAL_PTR_XE135_CONC] = (double)ptr;
        }

      /* Check equilibrium Sm-149 calculation */

      if ((long)RDB[mat0 + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
        {
          ptr = NewStat("MAT_PM149_CONC", 1, 1);
          WDB[mat + MATERIAL_PTR_PM149_CONC] = (double)ptr;

          ptr = NewStat("MAT_SM149_CONC", 1, 1);
          WDB[mat + MATERIAL_PTR_SM149_CONC] = (double)ptr;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /* Print newline */

  fprintf(outp, "OK.\n");
}

/*****************************************************************************/
