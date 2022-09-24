/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : materialvolumes.c                              */
/*                                                                           */
/* Created:       2011/07/03 (JLe)                                           */
/* Last modified: 2019/06/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates material volumes from cell volumes                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MaterialVolumes:"

/*****************************************************************************/

void MaterialVolumes()
{
  long loc0, loc1, mat, mat0, cell, det, n, uni, ptr, idx, prnt, prev;
  long tet, ntet, tetlist, i;
  double vol;

  /***************************************************************************/

  /***** Init ****************************************************************/

  /* Reset volumes */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset volume, mass and count */

      WDB[mat + MATERIAL_VOLUME] = 0.0;
      WDB[mat + MATERIAL_MASS] = 0.0;
      WDB[mat + MATERIAL_VOL_COUNT] = 0.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /**** Materials given in mvol-card *****************************************/

  /* Reset previous */

  prev = -1;

  /* Loop over cards */

  ptr = (long)RDB[DATA_PTR_MVOL0];
  while (ptr > VALID_PTR)
    {
      /* Start from previous */

      if ((mat = prev) > VALID_PTR)
        {
          /* Loop over materials */

          while (mat > VALID_PTR)
            {
              /* Update previous */

              prev = mat;

              /* Check pointer to parent and get divisor index */

              if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
                idx = (long)RDB[mat + MATERIAL_DIV_ZONE_IDX];
              else
                {
                  mat0 = mat;
                  idx = 0;
                }

              /* Compare */

              if ((idx == (long)RDB[ptr + MVOL_REG_IDX]) &&
                  CompareStr(ptr + MVOL_PTR_MAT, mat0 + MATERIAL_PTR_NAME))
                {
                  /* Put volume */

                  WDB[mat + MATERIAL_VOLUME_GIVEN] = RDB[ptr + MVOL_VOL];

                  /* Break loop */

                  break;
                }

              /* Next material */

              mat = NextItem(mat);
            }
        }

      /* Check if found */

      if (mat < VALID_PTR)
        {
          /* Loop over all materials */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Update previous */

              prev = mat;

              /* Check pointer to parent and get divisor index */

              if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
                idx = (long)RDB[mat + MATERIAL_DIV_ZONE_IDX];
              else
                {
                  mat0 = mat;
                  idx = 0;
                }

              /* Compare */

              if ((idx == (long)RDB[ptr + MVOL_REG_IDX]) &&
                  CompareStr(ptr + MVOL_PTR_MAT, mat0 + MATERIAL_PTR_NAME))
                {
                  /* Put volume */

                  WDB[mat + MATERIAL_VOLUME_GIVEN] = RDB[ptr + MVOL_VOL];

                  /* Break loop */

                  break;
                }

              /* Next material */

              mat = NextItem(mat);
            }
        }

      /* Check pointer */

      if (mat < VALID_PTR)
        Note(0, "No match for %s, reg. %ld in mvol card",
             GetText(ptr + MVOL_PTR_MAT), (long)RDB[ptr + MVOL_REG_IDX]);

      /* Next */

      ptr = NextItem(ptr);
    }

  /***************************************************************************/

  /***** Materials in regular cells ******************************************/

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Get volume */

      vol = RDB[cell + CELL_VOLUME];

      /* Get universe pointer */

      uni = (long)RDB[cell + CELL_PTR_UNI];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Check universe type */

      if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_SUPER)
        {
          /* Get multiplicity */

          if (((n = (long)RDB[cell + CELL_VOL_COUNT]) == 0) &&
              ((long)RDB[DATA_PTR_REP0] < VALID_PTR))
            {
              /* Used-flagia ei oo asetettu, tms? */

              if ((long)RDB[DATA_PTR_STL0] > VALID_PTR)
                Warn(FUNCTION_NAME, "Tälle pitää tehdä jotain");
              else
                Die(FUNCTION_NAME,"Cell %s not included in volume calculation",
                    GetText(cell + CELL_PTR_NAME));
            }

          /* Add to material volume */

          if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
            {
              WDB[mat + MATERIAL_VOLUME] = RDB[mat + MATERIAL_VOLUME] +
                ((double)n)*vol;
              WDB[mat + MATERIAL_VOL_COUNT] =
                RDB[mat + MATERIAL_VOL_COUNT] + (double)n;
            }
        }

      /* Next item */

      cell = NextItem(cell);
    }

  /***************************************************************************/

  /***** Materials in ustructured mesh based geometries **********************/

  /* Loop over unstructured mesh based geometries */

  loc0 = (long)RDB[DATA_PTR_UMSH0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to interface structure */

      loc1 = (long)RDB[loc0 + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      tetlist = (long)RDB[loc1 + IFC_PTR_TET_LIST];
      CheckPointer(FUNCTION_NAME, "(tetlist)", DATA_ARRAY, tetlist);

      ntet = (long)RDB[loc1 + IFC_NC];

      for (i = 0; i < ntet; i++)
        {
          /* Get pointer to tet */

          tet = (long)RDB[tetlist + i];
          CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

          /* Get pointer to parent */

          prnt = (long)RDB[tet + TET_PTR_PARENT];
          CheckPointer(FUNCTION_NAME, "(prnt)", DATA_ARRAY, prnt);

          /* Pointer to geometry cell */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Get volume */

          vol = RDB[cell + CELL_VOLUME];

          /* Get multiplicity */

          if (((n = (long)RDB[cell + CELL_VOL_COUNT]) == 0) &&
              ((long)RDB[DATA_PTR_REP0] < VALID_PTR))
            Die(FUNCTION_NAME, "Cell %s not included in volume calculation",
                GetText(cell + CELL_PTR_NAME));

          /* Add to material volume */

          if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
            {
              WDB[mat + MATERIAL_VOLUME] =
                RDB[mat + MATERIAL_VOLUME] + vol/(double)n;
              WDB[mat + MATERIAL_VOL_COUNT] =
                RDB[mat + MATERIAL_VOL_COUNT] + 1.0/(double)n;
            }

        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Burnable materials *************************************************/

  /* Depletion zones and parent materials: In Serpent 1 type division  */
  /* the volumes of divided regions are calculated, and summed to give */
  /* the total volume. In the new approach the calculation gives the   */
  /* total volume, which is divided between the regions. */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_S1)
        {
          /* Get pointer to parent */

          ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Put volume */

          WDB[ptr + MATERIAL_VOLUME] = RDB[ptr + MATERIAL_VOLUME] +
            RDB[mat + MATERIAL_VOLUME];
        }
      else if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
        {
          /* Get pointer to parent */

          ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Put volume */

          WDB[mat + MATERIAL_VOLUME] =
            RDB[ptr + MATERIAL_VOLUME]/RDB[ptr + MATERIAL_DIV_N_TOT_ZONES];

          if (RDB[mat + MATERIAL_VOLUME_GIVEN] < 0.0)
            WDB[mat + MATERIAL_VOLUME_GIVEN] =
              RDB[ptr + MATERIAL_VOLUME_GIVEN]/
              RDB[ptr + MATERIAL_DIV_N_TOT_ZONES];
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Materials in continuous reprocessing ********************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material is associated with reprocessing */

      if ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] > VALID_PTR)
        if ((long)RDB[mat + MATERIAL_VOL_COUNT] == 0)
          WDB[mat + MATERIAL_VOL_COUNT] = 1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Materials in activation detectors ***********************************/

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
     /* Loop over activation bins */

      ptr = (long)RDB[det + DET_PTR_ABINS];
      while (ptr > VALID_PTR)
        {
          /* Get material pointer */

          mat = (long)RDB[ptr + DET_MBIN_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Check volume */

          if (RDB[mat + MATERIAL_VOLUME_GIVEN] < ZERO)
            Error(det, "Volume of material %s must be given",
                  GetText(mat + MATERIAL_PTR_NAME));

          /* Set count */

          WDB[mat + MATERIAL_VOL_COUNT] = 1.0;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next detector */

      det = NextItem(det);
    }

  /***************************************************************************/

  /***** Remaining stuff *****************************************************/

  /* Override user-given values */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if volume is given */

      if ((vol = RDB[mat + MATERIAL_VOLUME_GIVEN]) >= 0.0)
        {
          /* Check count */

          if ((long)RDB[mat + MATERIAL_VOL_COUNT] > 0)
            WDB[mat + MATERIAL_VOLUME] = vol;

          /* Check parent count */

          if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
            if ((long)RDB[mat0 + MATERIAL_VOL_COUNT] > 0)
              WDB[mat + MATERIAL_VOLUME] = vol;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Reset volumes of divided materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to parent and reset volume */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        WDB[mat0 + MATERIAL_VOLUME] = 0.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Calculate volumes of divided materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to parent and add to volume */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        WDB[mat0 + MATERIAL_VOLUME] = RDB[mat0 + MATERIAL_VOLUME] +
          RDB[mat + MATERIAL_VOLUME];

      /* Next material */

      mat = NextItem(mat);
    }

  /* Set initial mass flag for CalculateMasses() */

  WDB[DATA_BURN_CALC_INI_MASS] = (double)YES;

  /***************************************************************************/
}

/*****************************************************************************/
