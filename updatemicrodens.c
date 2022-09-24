/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updatemicrodens.c                              */
/*                                                                           */
/* Created:       2019/10/23 (JLe)                                           */
/* Last modified: 2019/11/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates average nuclide densities for micro-depeltion.    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateMicroDens:"

/*****************************************************************************/

void UpdateMicroDens()
{
  long loc0, loc1, loc2, mat, iso, ptr, nrea, n, m, nbtch;
  double V;

  /* Get batch number */

  if ((long)RDB[DATA_CYCLE_IDX] == (long)RDB[DATA_CRIT_SKIP] - 1)
    nbtch = 0;
  else
    nbtch = (long)RDB[DATA_MICRO_CALC_BATCH_COUNT];

  /***************************************************************************/

  /***** Densities for micro-depletion xs ************************************/

  /* Loop over definitions */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Number of reactions */

      nrea = (long)RDB[loc0 + MDEP_N_REA];
      CheckValue(FUNCTION_NAME, "nrea", "", nrea, 0, 1000000);

      /* Pointer to batch average densities */

      ptr = (double)RDB[loc0 + MDEP_PTR_BTCH_AVG_ADENS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Multiply */

      for (n = 0; n < nrea; n++)
        WDB[ptr + n] = RDB[ptr + n]*((double)nbtch);

      /* Pointer to average densities */

      ptr = (double)RDB[loc0 + MDEP_PTR_AVG_ADENS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Reset data */

      memset(&WDB[ptr], 0.0, nrea*sizeof(double));

      /* Loop over mdep materials */

      loc1 = (long)RDB[loc0 + MDEP_PTR_MAT];
      while (loc1 > VALID_PTR)
        {
          /* Loop over materials */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Compare */

              if (((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] ==
                   (long)RDB[loc1 + MDEP_MAT_PTR_MAT]) ||
                  (((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
                   && ((long)RDB[loc1 + MDEP_MAT_PTR_MAT] == mat)))
                {
                  /* Get volume */

                  if ((V = RDB[mat + MATERIAL_VOLUME]) < ZERO)
                    Error(loc0, "Material %s volume is needed for calculation",
                          GetText(mat + MATERIAL_PTR_NAME));

                  /* Calculate ratio */

                  V = V/RDB[loc0 + MDEP_VOLUME];
                  CheckValue(FUNCTION_NAME, "V", "", V, ZERO, INFTY);

                  /* Check */

                  if (V > 1.0)
                    Error(loc0, "Material %s volume exceeds total",
                          GetText(mat + MATERIAL_PTR_NAME));

                  /* Loop over reactions */

                  loc2 = (long)RDB[loc1 + MDEP_MAT_PTR_REA];
                  while (loc2 > VALID_PTR)
                    {
                      /* Get index */

                      n = (long)RDB[loc2 + MDEP_REA_IDX];
                      CheckValue(FUNCTION_NAME, "n", "", n, 0,
                                 (long)RDB[loc0 + MDEP_N_REA] - 1);

                      /* Pointer to isotope */

                      iso = (long)RDB[loc2 + MDEP_REA_PTR_ISO];
                      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

                      /* Check partial index and add to average density */

                      if ((long)RDB[loc2 + MDEP_REA_PARTIAL_IDX] == 0)
                        {
                          /* Cycle-averaged */

                          ptr = (long)RDB[loc0 + MDEP_PTR_AVG_ADENS];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          WDB[ptr + n] = RDB[ptr + n] +
                            V*RDB[iso + COMPOSITION_ADENS];

                          /* Batch-averaged */

                          ptr = (long)RDB[loc0 + MDEP_PTR_BTCH_AVG_ADENS];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          WDB[ptr + n] = RDB[ptr + n] +
                            V*RDB[iso + COMPOSITION_ADENS];
                        }

                      /* Next */

                      loc2 = NextItem(loc2);
                    }
                }

              /* Next */

              mat = NextItem(mat);
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Batch-averaged */

      ptr = (long)RDB[loc0 + MDEP_PTR_BTCH_AVG_ADENS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (n = 0; n < nrea; n++)
        WDB[ptr + n] = RDB[ptr + n]/((double)nbtch + 1.0);

      /* Loop over mdep materials */

      loc1 = (long)RDB[loc0 + MDEP_PTR_MAT];
      while (loc1 > VALID_PTR)
        {
          /* Put partial densities */

          loc2 = (long)RDB[loc1 + MDEP_MAT_PTR_REA];
          while (loc2 > VALID_PTR)
            {
              /* Check if partial */

              if ((long)RDB[loc2 + MDEP_REA_PARTIAL_IDX] == 0)
                {
                  /* Pointer to next */

                  loc2 = NextItem(loc2);

                  /* Cycle loop */

                  continue;
                }

              /* Find matching nuclide that is not partial */

              ptr = (long)RDB[loc1 + MDEP_MAT_PTR_REA];
              while (ptr > VALID_PTR)
                {
                  /* Check */

                  if ((long)RDB[ptr + MDEP_REA_PARTIAL_IDX] == 0)
                    if ((long)RDB[ptr + MDEP_REA_PTR_ISO] ==
                        (long)RDB[loc2 + MDEP_REA_PTR_ISO])
                      break;

                  /* Next */

                  ptr = NextItem(ptr);
                }

              /* Check */

              if (ptr < VALID_PTR)
                Die(FUNCTION_NAME, "No match found");

              /* Get indexes */

              n = (long)RDB[loc2 + MDEP_REA_IDX];
              m = (long)RDB[ptr + MDEP_REA_IDX];

              /* Put density */

              ptr = (long)RDB[loc0 + MDEP_PTR_AVG_ADENS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              WDB[ptr + n] = RDB[ptr + m];

              /* Put batch-averaged density */

              ptr = (long)RDB[loc0 + MDEP_PTR_BTCH_AVG_ADENS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              WDB[ptr + n] = RDB[ptr + m];

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next definition */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Densities for poison xs *********************************************/

  /* Check poison calculation flag */

  if ((long)RDB[DATA_OPTI_POISON_CALC] == NO)
    return;

  /* Reset volume */

  WDB[DATA_POISON_XS_VF] = 0.0;

  /* Reset densities */

  WDB[DATA_POISON_XS_I135_ADENS] = 0.0;
  WDB[DATA_POISON_XS_XE135_ADENS] = 0.0;
  WDB[DATA_POISON_XS_XE135M_ADENS] = 0.0;
  WDB[DATA_POISON_XS_PM147_ADENS] = 0.0;
  WDB[DATA_POISON_XS_PM148_ADENS] = 0.0;
  WDB[DATA_POISON_XS_PM148M_ADENS] = 0.0;
  WDB[DATA_POISON_XS_PM149_ADENS] = 0.0;
  WDB[DATA_POISON_XS_SM149_ADENS] = 0.0;

  /* Reset batch-wise densities */

  WDB[DATA_POISON_XS_I135_BTCH_ADENS] =
    RDB[DATA_POISON_XS_I135_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_XE135_BTCH_ADENS] =
    RDB[DATA_POISON_XS_XE135_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_XE135M_BTCH_ADENS] =
    RDB[DATA_POISON_XS_XE135M_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_PM147_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM147_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_PM148_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM148_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_PM148M_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM148M_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_PM149_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM149_BTCH_ADENS]*((double)nbtch);
  WDB[DATA_POISON_XS_SM149_BTCH_ADENS] =
    RDB[DATA_POISON_XS_SM149_BTCH_ADENS]*((double)nbtch);

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check material type */

      if (((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
          && ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT)
          && ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Get volume */

          if ((V = RDB[mat + MATERIAL_VOLUME]) < ZERO)
            Error(0, "Material %s volume is needed for poison calculation",
                  GetText(mat + MATERIAL_PTR_NAME));

          /* Add to volume */

          WDB[DATA_POISON_XS_VF] = RDB[DATA_POISON_XS_VF] + V;

          /* Calculate ratio */

          V = V/RDB[DATA_POISON_XS_VOL];
          CheckValue(FUNCTION_NAME, "V", "", V, ZERO, INFTY);

          /* Check */

          if (V > 1.0)
            Error(0, "Material %s volume exceeds total",
                  GetText(mat + MATERIAL_PTR_NAME));

          /* Add to densities */

          if ((iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_I135_ADENS] =
                RDB[DATA_POISON_XS_I135_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_I135_BTCH_ADENS] =
                RDB[DATA_POISON_XS_I135_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_XE135_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_XE135_ADENS] =
                RDB[DATA_POISON_XS_XE135_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_XE135_BTCH_ADENS] =
                RDB[DATA_POISON_XS_XE135_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_XE135M_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_XE135M_ADENS] =
                RDB[DATA_POISON_XS_XE135M_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_XE135M_BTCH_ADENS] =
                RDB[DATA_POISON_XS_XE135M_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_PM147_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_PM147_ADENS] =
                RDB[DATA_POISON_XS_PM147_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_PM147_BTCH_ADENS] =
                RDB[DATA_POISON_XS_PM147_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_PM148_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_PM148_ADENS] =
                RDB[DATA_POISON_XS_PM148_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_PM148_BTCH_ADENS] =
                RDB[DATA_POISON_XS_PM148_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_PM148M_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_PM148M_ADENS] =
                RDB[DATA_POISON_XS_PM148M_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_PM148M_BTCH_ADENS] =
                RDB[DATA_POISON_XS_PM148M_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_PM149_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_PM149_ADENS] =
                RDB[DATA_POISON_XS_PM149_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_PM149_BTCH_ADENS] =
                RDB[DATA_POISON_XS_PM149_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }

          if ((iso = (long)RDB[mat + MATERIAL_PTR_SM149_ISO]) > VALID_PTR)
            {
              WDB[DATA_POISON_XS_SM149_ADENS] =
                RDB[DATA_POISON_XS_SM149_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
              WDB[DATA_POISON_XS_SM149_BTCH_ADENS] =
                RDB[DATA_POISON_XS_SM149_BTCH_ADENS]
                + V*RDB[iso + COMPOSITION_ADENS];
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Set batch-wise densities */

  WDB[DATA_POISON_XS_I135_BTCH_ADENS] =
    RDB[DATA_POISON_XS_I135_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_XE135_BTCH_ADENS] =
    RDB[DATA_POISON_XS_XE135_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_XE135M_BTCH_ADENS] =
    RDB[DATA_POISON_XS_XE135M_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_PM147_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM147_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_PM148_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM148_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_PM148M_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM148M_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_PM149_BTCH_ADENS] =
    RDB[DATA_POISON_XS_PM149_BTCH_ADENS]/((double)nbtch + 1.0);
  WDB[DATA_POISON_XS_SM149_BTCH_ADENS] =
    RDB[DATA_POISON_XS_SM149_BTCH_ADENS]/((double)nbtch + 1.0);

  /***************************************************************************/
}

/*****************************************************************************/
