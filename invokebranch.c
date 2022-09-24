/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : invokebranch.c                                 */
/*                                                                           */
/* Created:       2014/04/08 (JLe)                                           */
/* Last modified: 2020/03/06 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Invokes depletion branch in automated burnup sequence for    */
/*              group constant generation                                    */
/*                                                                           */
/* Comments: - Called right after ReadInput() to invoke changes before       */
/*             anything else is done.                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InvokeBranch:"

/*****************************************************************************/

void InvokeBranch(long loc0)
{
  long loc1, loc2, loc3, mat, cell, lat, nst, ptr;
  double sum, f, rep;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /***************************************************************************/

  /***** Replace materials ***************************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_REPLACE_MAT];
  while (loc1 > VALID_PTR)
    {
      /* Find first material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1,
                         mat + MATERIAL_PTR_NAME))
            break;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check pointer */

      if (mat < VALID_PTR)
        Die(FUNCTION_NAME, "Material %s is not defined",
              GetText(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1));

      /* Remove original */

      RemoveItem(mat);

      /* Find second material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT2,
                         mat + MATERIAL_PTR_NAME))
            break;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check pointer */

      if (mat < VALID_PTR)
        Die(FUNCTION_NAME, "Material %s is not defined",
            GetText(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT2));

      /* Override name */

      WDB[mat + MATERIAL_PTR_NAME] = RDB[loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1];

      /* Reset replaced flag */

      ResetOption(mat + MATERIAL_OPTIONS, OPT_REPLACED_MAT);

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/

  /***** Replace universes ***************************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_REPLACE_UNI];
  while (loc1 > VALID_PTR)
    {
      /* Reset count */

      rep = 0;

      /* Loop over cells */

      cell = (long)RDB[DATA_PTR_C0];
      while (cell > VALID_PTR)
        {
          /* Compare with first and second universe */

          if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1,
                         cell + CELL_PTR_UNI))
            {
              /* Copy pointer */

              ptr = cell;

              /* Pointer to next */

              cell = NextItem(cell);

              /* Remove cell from list */

              RemoveItem(ptr);

              /* Cycle loop */

              continue;
            }
          else if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2,
                              cell + CELL_PTR_UNI))
            {
              /* Replace pointer */

              WDB[cell + CELL_PTR_UNI] =
                RDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1];

              /* Add counter */

              rep++;
            }

          /* Next cell */

          cell = NextItem(cell);
        }

      /* Loop over lattices */

      lat = (long)RDB[DATA_PTR_L0];
      while (lat > VALID_PTR)
        {
          /* Compare with first and second universe */

          if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1,
                         lat + LAT_PTR_NAME))
            {
              /* Copy pointer */

              ptr = lat;

              /* Pointer to next */

              lat = NextItem(lat);

              /* remove lattice from list */

              RemoveItem(ptr);

              /* Cycle loop */

              continue;
            }
          else if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2,
                              lat + LAT_PTR_NAME))
            {
              /* Replace pointer */

              WDB[lat + LAT_PTR_NAME] =
                RDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1];

              /* Add counter */

              rep++;
            }

          /* next lattice */

          lat = NextItem(lat);
        }

      /* Loop over nests */

      nst = (long)RDB[DATA_PTR_NST0];
      while (nst > VALID_PTR)
        {
          /* Compare with first and second universe */

          if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1,
                         nst + NEST_PTR_NAME))
            {
              /* Copy pointer */

              ptr = nst;

              /* Pointer to next */

              nst = NextItem(nst);

              /* Remove nest from list */

              RemoveItem(ptr);

              /* Cycle loop */

              continue;
            }
          else if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2,
                              nst + NEST_PTR_NAME))
            {
              /* Replace pointer */

              WDB[nst + NEST_PTR_NAME] =
                RDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1];

              /* Add counter */

              rep++;
            }

          /* next nest */

          nst = NextItem(nst);
        }

      /* Check count */

      if (rep == 0)
        Error(loc0, "Universe %s not defined",
              GetText(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2));

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/

  /***** Changes in material states ******************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_STP];
  while (loc1 > VALID_PTR)
    {
      /* Find material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(loc1 + DEP_BRA_STP_PTR_MAT, mat + MATERIAL_PTR_NAME))
            break;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check pointer */

      if (mat < VALID_PTR)
        Die(FUNCTION_NAME, "Material %s is not defined",
              GetText(loc1 + DEP_BRA_STP_PTR_MAT));

      /* Set density adjustment factor (11.1.2019 / JLe) */

      if ((RDB[loc1 + DEP_BRA_STP_DENSITY] == -INFTY) &&
          (RDB[mat + MATERIAL_ADENS] != -INFTY))
        Error(loc0, "Sum over composition cannot be used with adjustment");
      else if ((RDB[loc1 + DEP_BRA_STP_DENSITY] != -INFTY) &&
               (RDB[mat + MATERIAL_ADENS] == -INFTY))
        Error(loc0, "Sum over composition cannot be used with adjustment");
      else
        WDB[mat + MATERIAL_RESTART_ADENS_F] = RDB[loc1 + DEP_BRA_STP_DENSITY]/
          RDB[mat + MATERIAL_ADENS];

      /* Check */

      if (RDB[mat + MATERIAL_RESTART_ADENS_F] < 0.0)
        Error(loc0, "Density must be given in consistent units");

      /* Override density and temperature */

      WDB[mat + MATERIAL_ADENS] = RDB[loc1 + DEP_BRA_STP_DENSITY];

      if (RDB[loc1 + DEP_BRA_STP_TEMP] > 0.0)
        {
          /* Reset Doppler temperature (11.1.2018 / 2.1.30 / JLe) */

          WDB[mat + MATERIAL_DOPPLER_TEMP] = -1.0;

          /* Set coefficient temperature */

          WDB[mat + MATERIAL_COEF_TEMP] = RDB[loc1 + DEP_BRA_STP_TEMP];
        }

      /* Check if sum is given */

      if (RDB[mat + MATERIAL_ADENS] == -INFTY)
        {
          /* Loop over composition and calculate sum */

          sum = 0.0;
          ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (ptr > VALID_PTR)
            {
              /* Add to sum */

              sum = sum + RDB[ptr + COMPOSITION_ADENS];

              /* Next */

              ptr = NextItem(ptr);
            }

          /* Put value */

          WDB[mat + MATERIAL_ADENS] = sum;
        }

      /* Check if material has S(a,b) data */

      loc3 = (long)RDB[mat + MATERIAL_PTR_SAB];
      while (loc3 > VALID_PTR)
        {
          /* Find match */

          loc2 = (long)RDB[loc1 + DEP_BRA_STP_PTR_SAB];
          while (loc2 > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(loc2 + DEP_BRA_STP_SAB_PTR_THERM,
                             loc3 + THERM_PTR_ALIAS))
                break;

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Check */

          if (loc2 < VALID_PTR)
            Die(FUNCTION_NAME, "Missing S(a,b) entries for %s",
                  GetText(loc3 + THERM_PTR_ALIAS));

          /* Next */

          loc3 = NextItem(loc3);
        }

      /* Override S(a,b) data */

      loc2 = (long)RDB[loc1 + DEP_BRA_STP_PTR_SAB];
      while (loc2 > VALID_PTR)
        {
          /* Check if material has associated S(a,b) data */

          if ((loc3 = (long)RDB[mat + MATERIAL_PTR_SAB]) < VALID_PTR)
            Die(FUNCTION_NAME,
                "Material %s is not associated with S(a,b) data",
                GetText(mat + MATERIAL_PTR_NAME));

          /* Check that names match */

          while (loc3 > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(loc2 + DEP_BRA_STP_SAB_PTR_THERM,
                             loc3 + THERM_PTR_ALIAS))
                break;

              /* Next */

              loc3 = NextItem(loc3);
            }

          /* Check pointer */

          if (loc3 < VALID_PTR)
            Die(FUNCTION_NAME,
                "Material %s is not associated with S(a,b) data %s",
                GetText(mat + MATERIAL_PTR_NAME),
                GetText(loc2 +  DEP_BRA_STP_SAB_PTR_THERM));

          /* Find match */

          loc3 = (long)RDB[DATA_PTR_T0];
          while (loc3 > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(loc2 + DEP_BRA_STP_SAB_PTR_THERM,
                             loc3 + THERM_PTR_ALIAS))
                break;

              /* Next */

              loc3 = NextItem(loc3);
            }

          /* Check pointer */

          if (loc3 < VALID_PTR)
            Die(FUNCTION_NAME, "S(a,b) library %s is not defined",
                  GetText(loc2 + DEP_BRA_STP_SAB_PTR_THERM));

          /* Override temperature and set library names */

          if (RDB[loc1 + DEP_BRA_STP_TEMP] > 0.0)
            {
              /* Put temperature */

              WDB[loc3 + THERM_T] = RDB[loc1 + DEP_BRA_STP_TEMP];

              /* Reset pointer */

              WDB[loc3 + THERM_PTR_SAB] = NULLPTR;

              /* Put data */

              ptr = NewItem(loc3 + THERM_PTR_SAB, SAB_BLOCK_SIZE);
              WDB[ptr + SAB_PTR_NAME] = RDB[loc2 + DEP_BRA_STP_SAB_PTR_LIB1];

              ptr = NewItem(loc3 + THERM_PTR_SAB, SAB_BLOCK_SIZE);
              WDB[ptr + SAB_PTR_NAME] = RDB[loc2 + DEP_BRA_STP_SAB_PTR_LIB2];

              /* Set interpolation mode */

              WDB[loc3 + THERM_INTERP_MODE] = (double)THERM_INTERP_MAKXSF;
            }

          /* Next */

          loc2 = NextItem(loc2);
        }

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/

  /***** Universe transformation *********************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_TRANS];
  while (loc1 > VALID_PTR)
    {
      /* Find transformation */

      loc2 = (long)RDB[DATA_PTR_TR0];
      while (loc2 > VALID_PTR)
        {
          /* Compare names */

          if (CompareStr(loc1 + DEP_BRA_TRANS_PTR_TRANS, loc2 + TRANS_PTR_NAME))
            break;

          /* Next */

          loc2 = NextItem(loc2);
        }

      /* Check pointer */

      if (loc2 < VALID_PTR)
        Die(FUNCTION_NAME, "Transformation %s is not defined",
            GetText(loc1 + DEP_BRA_TRANS_PTR_TRANS));

      /* Create copy (NOTE: Historiavariaatioissa luodaan kopiot koska   */
      /* transformaatioiden pitää kumuloitua. Tämä on ehkä se miten myös */
      /* muiden branchien pitäisi toimia, mutta if-lause lisättiin ettei */
      /* sotketa aikaisempaa käyttäytyämistä, 6.3.2020 / 2.1.32 / JLe).  */

      if ((long)RDB[DATA_PTR_HISV0] > VALID_PTR)
        loc2 = DuplicateItem(loc2);

      /* Put name */

      WDB[loc2 + TRANS_PTR_NAME] = RDB[loc1 + DEP_BRA_TRANS_PTR_UNI];

      /* Check type and copy pointer */

      if ((long)RDB[loc2 + TRANS_TYPE] == TRANS_TYPE_UNI)
        WDB[loc2 + TRANS_PTR_UNI] = RDB[loc2 + TRANS_PTR_NAME];
      else if ((long)RDB[loc2 + TRANS_TYPE] == TRANS_TYPE_FILL)
        WDB[loc2 + TRANS_PTR_CELL] = RDB[loc2 + TRANS_PTR_NAME];
      else if ((long)RDB[loc2 + TRANS_TYPE] == TRANS_TYPE_SURF)
        WDB[loc2 + TRANS_PTR_SURF] = RDB[loc2 + TRANS_PTR_NAME];
      else
        Die(FUNCTION_NAME, "Invalid transformation type %ld",
            (long)RDB[loc2 + TRANS_TYPE]);

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/

  /***** Other stuff *********************************************************/

  /* Replace universes for group constant generation */

  if ((loc1 = (long)RDB[loc0 + DEP_BRA_PTR_GCU]) > VALID_PTR)
    {
      /* Get pointer to existing, and check that list is single-valued */

      if ((ptr = (long)RDB[DATA_PTR_GCU0]) > VALID_PTR)
        if (ListSize(ptr) > 1)
          Error(loc0, "Replaced GCU-list must be single-valued");

      /* Replace pointer */

      WDB[DATA_PTR_GCU0] = (double)loc1;

      /* Loop over ADF's */

      ptr = (long)RDB[DATA_PTR_ADF0];
      while (ptr > VALID_PTR)
        {
          /* Replace universe name */

          WDB[ptr + ADF_PTR_GCU] = RDB[loc1 + GCU_PTR_UNIV];

          /* Pointer to next */

          ptr = NextItem(ptr);
        }

      /* Loop over pin power distributions */

      ptr = (long)RDB[DATA_PTR_PPW0];
      while (ptr > VALID_PTR)
        {
          /* Replace universe name */

          WDB[ptr + PPW_PTR_GCU] = RDB[loc1 + GCU_PTR_UNIV];

          /* Pointer to next */

          ptr = NextItem(ptr);
        }

      /* Loop over ALB's */

      ptr = (long)RDB[DATA_PTR_ALB0];
      while (ptr > VALID_PTR)
        {
          /* Replace universe name */

          WDB[ptr + ALB_PTR_GCU] = RDB[loc1 + GCU_PTR_UNIV];

          /* Pointer to next */

          ptr = NextItem(ptr);
        }
    }

  /* Put zero poison flags */

  WDB[DATA_RESTART_READ_ZERO_XE] = RDB[loc0 + DEP_BRA_ZERO_XE];
  WDB[DATA_RESTART_READ_ZERO_SM] = RDB[loc0 + DEP_BRA_ZERO_SM];

  /* Adjust normalization */

  if ((f = RDB[loc0 + DEP_BRA_NORM]) >= 0.0)
    {
      /* Loop over normalizations */

      loc1 = (long)RDB[DATA_PTR_NORM];
      while (loc1 > VALID_PTR)
        {
          /* Adjust */

          if (RDB[loc1 + NORM_POWER] > 0.0)
            WDB[loc1 + NORM_POWER] = RDB[loc1 + NORM_POWER]*f;
          else if (RDB[loc1 + NORM_POWDENS] > 0.0)
            WDB[loc1 + NORM_POWDENS] = RDB[loc1 + NORM_POWDENS]*f;
          else if (RDB[loc1 + NORM_GENRATE] > 0.0)
            WDB[loc1 + NORM_GENRATE] = RDB[loc1 + NORM_GENRATE]*f;
          else if (RDB[loc1 + NORM_FISSRATE] > 0.0)
            WDB[loc1 + NORM_FISSRATE] = RDB[loc1 + NORM_FISSRATE]*f;
          else if (RDB[loc1 + NORM_ABSRATE] > 0.0)
            WDB[loc1 + NORM_ABSRATE] = RDB[loc1 + NORM_ABSRATE]*f;
          else if (RDB[loc1 + NORM_LOSSRATE] > 0.0)
            WDB[loc1 + NORM_LOSSRATE] = RDB[loc1 + NORM_LOSSRATE]*f;
          else if (RDB[loc1 + NORM_FLUX] > 0.0)
            WDB[loc1 + NORM_FLUX] = RDB[loc1 + NORM_FLUX]*f;
          else if (RDB[loc1 + NORM_SRCRATE] > 0.0)
            WDB[loc1 + NORM_SRCRATE] = RDB[loc1 + NORM_SRCRATE]*f;
          else if (RDB[loc1 + NORM_SFRATE] > 0.0)
            WDB[loc1 + NORM_SFRATE] = RDB[loc1 + NORM_SFRATE]*f;

          /* Next */

          loc1 = NextItem(loc1);
        }
    }

  /***************************************************************************/
}

/*****************************************************************************/
