/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkcoefcalc.c                                */
/*                                                                           */
/* Created:       2014/08/18 (JLe)                                           */
/* Last modified: 2020/06/23 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Checks that coefficient matrix and branches are correctly    */
/*              defined before running the calculation.                      */
/*                                                                           */
/* Comments: - The purpose is to stop the calculation early, rather than     */
/*             terminating with an error later on.                           */
/*                                                                           */
/*           - TODO: Palamapisteiden ja universumien testaamiseen pitäis     */
/*                   keksiä jotain.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckCoefCalc:"

/*****************************************************************************/

void CheckCoefCalc()
{
  long loc0, loc1, loc2, ns, n, ptr, bra, mat, tra1, tra2, sab1, sab2;
  char fname[MAX_STR];
  double pt;
  FILE *fp;

  /* Check for coeffiecient, history variation and case matrix calculations */

  if (((long)RDB[DATA_PTR_COEF0] < VALID_PTR) &&
      ((long)RDB[DATA_PTR_HISV0] < VALID_PTR) &&
      ((long)RDB[DATA_CASEMTX_RUN] == NO))
    return;

  /***************************************************************************/

  /***** Process case matrix calculations ************************************/

  /* Check if case matrix calculation is on */

  if ((long)RDB[DATA_CASEMTX_RUN] == YES)
    {
      /* Check definitions */

      if ((loc0 = (long)RDB[DATA_PTR_CASEMTX0]) < VALID_PTR)
        Error(0, "No case matrix definitions");

      /* Check name passed from command line */

      if ((long)RDB[DATA_CASEMTX_RUN_PTR_NAME] < VALID_PTR)
        Die(FUNCTION_NAME, "Pointer error");

      /* Find match */

      while (loc0 > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(loc0 + CASEMTX_PTR_NAME, DATA_CASEMTX_RUN_PTR_NAME))
            break;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Check pointer */

      if (loc0 < VALID_PTR)
        Error(0, "No match for case matrix \"%s\"",
              GetText(DATA_CASEMTX_RUN_PTR_NAME));

      /* Get pointer to coefficient matrix structure */

      loc1 = (long)RDB[loc0 + CASEMTX_PTR_COEF];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Put pointers */

      WDB[DATA_PTR_CASEMTX0] = (double)loc0;
      WDB[DATA_PTR_COEF0] = (double)loc1;

      /* Put number of calculations */

      WDB[DATA_TOT_COEF_CALC] = RDB[loc0 + CASEMTX_TOT_COEF_CALC];
      WDB[DATA_COEF_CALC_TOT_RUNS] = RDB[loc0 + CASEMTX_COEF_CALC_TOT_RUNS];

      /* Check with input */

      if (RDB[DATA_CASEMTX_RUN_COE_IDX] > RDB[DATA_TOT_COEF_CALC])
        Error(loc0, "Invalid restart index %ld (max %ld)",
              (long)RDB[DATA_CASEMTX_RUN_COE_IDX],
              (long)RDB[DATA_TOT_COEF_CALC]);

      /* Number of histories */

      if ((ns = (long)RDB[loc0 + CASEMTX_N_HIS]) < 1)
        Die(FUNCTION_NAME, "Number of histories shouldn't be zero");

      /* Check with given index */

      if ((n = (long)RDB[DATA_CASEMTX_RUN_HIS_IDX]) > ns)
        Error(loc0, "Invalid history index %ld (max %ld)", n, ns);
      else if (n > 0)
        {
          /* Pointer to histories */

          loc1 = (long)RDB[loc0 + CASEMTX_PTR_HIS] + n - 1;
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Find match */

          bra = (long)RDB[DATA_PTR_BRA0];
          while (bra > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(bra + DEP_BRA_PTR_NAME, loc1))
                break;

              /* Next */

              bra = NextItem(bra);
            }

          /* Check */

          if (bra < VALID_PTR)
            Error(loc0, "Branch %s is not defined", GetText(loc1));

          /* Set used flag */

          SetOption(bra + DEP_BRA_OPTIONS, OPT_USED);

          /* Set pointer */

          WDB[loc0 + CASEMTX_PTR_HIS] = (double)bra;

          /* Put name of restart files */

          sprintf(fname, "%s_%s_h%ld.wrk", GetText(DATA_PTR_INPUT_FNAME),
                  GetText(loc0 + CASEMTX_PTR_NAME), n);
          WDB[DATA_RESTART_WRITE_PTR_FNAME] = (double)PutText(fname);
          WDB[DATA_RESTART_READ_PTR_FNAME] = (double)PutText(fname);

          /* Check that restart file exists if single run */

          if ((long)RDB[DATA_CASEMTX_RUN_COE_IDX] > 0)
            {
              if ((fp = fopen(fname, "r")) == NULL)
                Error(loc0, "History calculations must be run first");
              else
                fclose(fp);
            }
        }
      else
        Error(loc0, "Invalid history index %ld", n);
    }

  /***************************************************************************/

  /***** Check coefficient matrixes ******************************************/

  /* Loop over coefficient matrixes */

  loc0 = (long)RDB[DATA_PTR_COEF0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over branch matrix */

      loc1 = (long)RDB[loc0 + COEF_PTR_MTX];
      while (loc1 > VALID_PTR)
        {
          /* Number of branches */

          if ((ns = (long)RDB[loc1 + COEF_MTX_N_BRA]) < 1)
            Die(FUNCTION_NAME, "Number of branches shouldn't be zero");

          /* Loop over branches */

          loc2 = (long)RDB[loc1 + COEF_MTX_PTR_BRA];
          for (n = 0; n < ns; n++)
            {
              /* Find match */

              bra = (long)RDB[DATA_PTR_BRA0];
              while (bra > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(bra + DEP_BRA_PTR_NAME, loc2))
                    break;

                  /* Next */

                  bra = NextItem(bra);
                }

              /* Check */

              if (bra < VALID_PTR)
                Error(loc0, "Branch %s is not defined", GetText(loc2));

              /* Set used flag */

              SetOption(bra + DEP_BRA_OPTIONS, OPT_USED);

              /* Next branch in matrix */

              loc2++;
            }

          /* Pointer to next */

          loc1 = NextItem(loc1);
        }

      /* Next coefficient matrix */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check history variations ********************************************/

  /* Reset point */

  pt = 0.0;

  /* Loop over history variations */

  loc0 = (long)RDB[DATA_PTR_HISV0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over points */

      loc1 = (long)RDB[loc0 + HISV_PTR_VAR];
      while (loc1 > VALID_PTR)
        {
          /* Check burnup (nolla testataan readinput.c:ssä) */

          if (pt == 0.0)
            pt = RDB[loc1 + HISV_VAR_BU];
          else if (fabs(RDB[loc1 + HISV_VAR_BU]) < fabs(pt))
            Error(loc0, "Variation points must be given in ascending order");
          else if (pt*RDB[loc1 + HISV_VAR_BU] < 0.0)
            Error(loc0, "Variation points must either in burnup or time");

          /* Number of variations */

          if ((ns = (long)RDB[loc1 + HISV_VAR_N_VAR]) < 1)
            Die(FUNCTION_NAME, "Number of variations shouldn't be zero");

          /* Pointer to data */

          loc2 = (long)RDB[loc1 + HISV_VAR_PTR_VAR];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Loop over variations */

          for (n = 0; n < ns; n++)
            {
              /* Find match */

              bra = (long)RDB[DATA_PTR_BRA0];
              while (bra > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(bra + DEP_BRA_PTR_NAME, loc2 + n))
                    break;

                  /* Next */

                  bra = NextItem(bra);
                }

              /* Check */

              if (bra < VALID_PTR)
                Error(loc0, "Branch %s is not defined", GetText(loc2 + n));

              /* Set used flag */

              SetOption(bra + DEP_BRA_OPTIONS, OPT_USED);
            }

          /* Pointer to next */

          loc1 = NextItem(loc1);
        }

      /* Next coefficient matrix */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check branches ******************************************************/

  /* Pointer to branches */

  if ((bra = (long)RDB[DATA_PTR_BRA0]) < VALID_PTR)
    Error(0, "No branch definitions");

  /* Remove unused */

  RemoveFlaggedItems(bra, DEP_BRA_OPTIONS, OPT_USED, NO);

  /* Loop over branches */

  bra = (long)RDB[DATA_PTR_BRA0];
  while (bra > VALID_PTR)
    {
      /***********************************************************************/

      /***** Check that materials are defined ********************************/

      ptr = (long)RDB[bra + DEP_BRA_PTR_REPLACE_MAT];
      while (ptr > VALID_PTR)
        {
          /* Find first material */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Compare name */

              if (CompareStr(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT1,
                             mat + MATERIAL_PTR_NAME))
                break;

              /* Next material */

              mat = NextItem(mat);
            }

          /* Check pointer */

          if (mat < VALID_PTR)
            Error(bra, "Material %s is not defined",
                  GetText(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT1));

          /* Get pointer to S(a,b) data */

          sab1 = (long)RDB[mat + MATERIAL_PTR_SAB];

          /* Find second material */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Compare name */

              if (CompareStr(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT2,
                             mat + MATERIAL_PTR_NAME))
                break;

              /* Next material */

              mat = NextItem(mat);
            }

          /* Check pointer */

          if (mat < VALID_PTR)
            Error(bra, "Material %s is not defined",
                  GetText(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT2));

          /* Get pointer to S(a,b) data */

          sab2 = (long)RDB[mat + MATERIAL_PTR_SAB];

          /* Set flag */

          SetOption(mat + MATERIAL_OPTIONS, OPT_REPLACED_MAT);
          SetOption(mat + MATERIAL_OPTIONS, OPT_USED);

          /* Check that S(a,b) libraries match */

          while (sab1 > VALID_PTR)
            {
              /* Loop over second */

              loc1 = sab2;
              while (loc1 > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(sab1 + THERM_PTR_ALIAS,
                                 loc1 + THERM_PTR_ALIAS))
                    break;

                  /* Next */

                  loc1 = NextItem(loc1);
                }

              /* Check */

              if (loc1 < VALID_PTR)
              Error(bra,
                    "Materials in stp variation must be associated with the same therm cards");

              /* Next */

              sab1 = NextItem(sab1);
            }

          /* Next */

          ptr = NextItem(ptr);
        }

      /***********************************************************************/

      /***** Check transformations *******************************************/

      ptr = (long)RDB[bra + DEP_BRA_PTR_TRANS];
      while (ptr > VALID_PTR)
        {
          /* Reset pointer */

          tra1 = -1;

          /* Find transformation */

          tra2 = (long)RDB[DATA_PTR_TR0];
          while (tra2 > VALID_PTR)
            {
              /* Compare names */

              if (CompareStr(ptr + DEP_BRA_TRANS_PTR_TRANS,
                             tra2 + TRANS_PTR_NAME))
                {
                  /* Set pointer if not defined */

                  if (tra1 < VALID_PTR)
                    tra1 = tra2;
                  else
                    Error(bra, "Multiple transformations named %s",
                          GetText(ptr + DEP_BRA_TRANS_PTR_TRANS));
                }

              /* Next */

              tra2 = NextItem(tra2);
            }

          /* Check pointer */

          if (tra1 < VALID_PTR)
            Error(bra, "Transformation %s is not defined",
                  GetText(ptr + DEP_BRA_TRANS_PTR_TRANS));

          /* Next */

          ptr = NextItem(ptr);
        }

      /***********************************************************************/

      /***** Checks for state-point variation ********************************/

      ptr = (long)RDB[bra + DEP_BRA_PTR_STP];
      while (ptr > VALID_PTR)
        {
          /* Find material */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Compare name */

              if (CompareStr(ptr + DEP_BRA_STP_PTR_MAT,
                             mat + MATERIAL_PTR_NAME))
                break;

              /* Next material */

              mat = NextItem(mat);
            }

          /* Check pointer */

          if (mat < VALID_PTR)
            Error(bra, "Material %s is not defined",
                  GetText(ptr + DEP_BRA_STP_PTR_MAT));

          /* Check if material has S(a,b) data */

          sab2 = (long)RDB[mat + MATERIAL_PTR_SAB];
          while (sab2 > VALID_PTR)
            {
              /* Find match */

              sab1 = (long)RDB[ptr + DEP_BRA_STP_PTR_SAB];
              while (sab1 > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(sab1 + DEP_BRA_STP_SAB_PTR_THERM,
                                 sab2 + THERM_PTR_ALIAS))
                    break;

                  /* Next */

                  sab1 = NextItem(sab1);
                }

              /* Check */

              if (sab1 < VALID_PTR)
                Error(bra, "Missing entry for %s",
                      GetText(sab2 + THERM_PTR_ALIAS));

              /* Next */

              sab2 = NextItem(sab2);
            }

          /* Loop over data */

          sab1 = (long)RDB[ptr + DEP_BRA_STP_PTR_SAB];
          while (sab1 > VALID_PTR)
            {

              /* Check that names match */

              sab2 = (long)RDB[mat + MATERIAL_PTR_SAB];
              while (sab2 > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(sab1 + DEP_BRA_STP_SAB_PTR_THERM,
                                 sab2 + THERM_PTR_ALIAS))
                    break;

                  /* Next */

                  sab2 = NextItem(sab2);
                }

              /* Check pointer */

              if (sab2 < VALID_PTR)
                Note(bra,
                     "Material %s has no thermal scattering data %s",
                     GetText(mat + MATERIAL_PTR_NAME),
                     GetText(sab1 + DEP_BRA_STP_SAB_PTR_THERM));

              /* Next */

              sab1 = NextItem(sab1);
            }

          /* Next */

          ptr = NextItem(ptr);
        }

      /***********************************************************************/

      /* Next branch */

      bra = NextItem(bra);
    }
}

/*****************************************************************************/
