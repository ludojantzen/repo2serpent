/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : poisoneq.c                                     */
/*                                                                           */
/* Created:       2012/12/04 (JLe)                                           */
/* Last modified: 2019/11/15 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates equilibrium concentrations for Xe-135 and Sm-149  */
/*                                                                           */
/* Comments: - Tähän vois nyt ottaa I-135 ja Pm-149 kaappauksen mukaan       */
/*                                                                           */
/*           - TODO: Xe-135m käsittely omana nuklidinaan                     */
/*                                                                           */
/*           - Noi konsentraatiot pakotetaan nolliksi jos noi                */
/*             DATA_RESTART_READ_ZERO_ZE ja ..._SM -optiot on asetettu       */
/*             (liittyy branch-laskuihin ilman myrkkyjä).                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PoisonEq:"

/*****************************************************************************/

void PoisonEq()
{
  long mat, ptr, iso;
  double norm, sum1, sum2, sum3, dir, prec, abs, val, div, vol, old;

  /* Check modes */

  if (((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] < 0) &&
      ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] < 0))
    return;

  /* Reset maximum concentrations */

  ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
  while (ptr > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[ptr + MAJORANT_EXTRA_TYPE] ==
          MAJORANT_EXTRA_FP_POISON_ITER)
        WDB[ptr + MAJORANT_EXTRA_FRAC] = 0.0;

      /* Pointer to next */

      ptr = NextItem(ptr);
    }

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Get normalization factor */

  norm = NormCoef(PARTICLE_TYPE_NEUTRON);
  CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);

  /***************************************************************************/

  /***** Equilibrium xenon ***************************************************/

  if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] > -1)
    {
      /* Reset sums and divider */

      sum1 = 0.0;
      sum2 = 0.0;
      sum3 = 0.0;
      div = 0.0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Get volume */

          vol = RDB[mat + MATERIAL_VOLUME];

          /* Add to total volume (NOTE: tossa katsotaan pointteri, koska */
          /* DD:ssä XENON_EQUIL_CALC -flagi ei ole kaikilla päällä.) */

          if (((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT) &&
              ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT) &&
              ((long)RDB[mat + MATERIAL_PTR_XE135_PROD_RATE] > VALID_PTR))
            div = div + vol;

          /* Check equilibrium flag and division */

          if (((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == NO) ||
              ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) ||
              (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)) ||
              (vol < ZERO))
            {
              /* Next material */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Pointer to composition */

          iso = (long)RDB[mat + MATERIAL_PTR_XE135_ISO];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

          /* Get production rates */

          ptr = (long)RDB[mat + MATERIAL_PTR_XE135_PROD_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          dir = norm*BufVal(ptr, 0)/vol;

          ptr = (long)RDB[mat + MATERIAL_PTR_I135_PROD_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          prec = norm*BufVal(ptr, 0)/vol;

          /* Absorption rate */

          ptr = (long)RDB[mat + MATERIAL_PTR_XE135_ABS_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          abs = BARN*norm*BufVal(ptr, 0)/vol;

          /* Add to total absorption */

          sum1 = sum1 + norm*BufVal(ptr, 0);

          /* Equilibrium I-135 concentration */

          if ((long)RDB[DATA_RESTART_READ_ZERO_XE] == YES)
            val = 0.0;
          else
            val = BARN*prec/RDB[DATA_I135_DC];

          ptr = (long)RDB[mat + MATERIAL_PTR_I135_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          sum2 = sum2 + val*vol;

          /* Equilibrium Xe-135 concentration */

          if ((long)RDB[DATA_RESTART_READ_ZERO_XE] == YES)
            val = 0.0;
          else
            val = BARN*(dir + prec)/(abs + RDB[DATA_XE135_DC]);

          ptr = (long)RDB[mat + MATERIAL_PTR_XE135_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          sum3 = sum3 + val*vol;

          /* Get old concentration */

          old = RDB[iso + COMPOSITION_ADENS];

          /* Set new concentration */

          WDB[iso + COMPOSITION_ADENS] = val;

          /* Adjust atomic density */

          WDB[mat + MATERIAL_ADENS] = RDB[mat + MATERIAL_ADENS] - old + val;

          /* Find nuclide in majorant extra list */

          ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
          while (ptr > VALID_PTR)
            {
              /* Compare pointers */

              if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] ==
                  ((long)RDB[iso + COMPOSITION_PTR_NUCLIDE]))
                {
                  /* Compare fraction */

                  if (val > RDB[ptr + MAJORANT_EXTRA_FRAC])
                    WDB[ptr + MAJORANT_EXTRA_FRAC] = val;

                  /* Break loop */

                  break;
                }

              /* Pointer to next */

              ptr = NextItem(ptr);
            }

          /* Check pointer */

          if (ptr < VALID_PTR)
            Die(FUNCTION_NAME, "Nuclide not found in majorant extra list");

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check divisor */

      CheckValue(FUNCTION_NAME, "div", "", div, ZERO, INFTY);

      /* Total I-135 and Xe-135 concentrations */

      ptr = (long)RDB[RES_I135_EQUIL_CONC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(sum2/div, ptr, 0);

      ptr = (long)RDB[RES_XE135_EQUIL_CONC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(sum3/div, ptr, 0);

      /* Total Xe-135 absorption rate */

      ptr = (long)RDB[RES_XE135_ABSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(sum1, ptr, 0);
    }

  /***************************************************************************/

  /***** Equilibrium samarium ************************************************/

  if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] > -1)
    {
      /* Reset sums and divider */

      sum1 = 0.0;
      sum2 = 0.0;
      sum3 = 0.0;
      div = 0.0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Get volume */

          vol = RDB[mat + MATERIAL_VOLUME];

          /* Add to total volume (NOTE: tossa katsotaan pointteri, koska */
          /* DD:ssä SAMARIUM_EQUIL_CALC -flagi ei ole kaikilla päällä.) */

          if (((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT) &&
              ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT) &&
              ((long)RDB[mat + MATERIAL_PTR_SM149_PROD_RATE] > VALID_PTR))
            div = div + vol;

          /* Check equilibrium flag and division */

          if (((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == NO) ||
              ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) ||
              (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)) ||
              (vol < ZERO))
            {
              /* Next material */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Pointer to composition */

          iso = (long)RDB[mat + MATERIAL_PTR_SM149_ISO];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

          /* Get production rates */

          ptr = (long)RDB[mat + MATERIAL_PTR_SM149_PROD_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          dir = norm*BufVal(ptr, 0)/vol;

          ptr = (long)RDB[mat + MATERIAL_PTR_PM149_PROD_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          prec = norm*BufVal(ptr, 0)/vol;

          /* Absorption rate */

          ptr = (long)RDB[mat + MATERIAL_PTR_SM149_ABS_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          abs = BARN*norm*BufVal(ptr, 0)/vol;

          /* Add to total absorption */

          sum1 = sum1 + norm*BufVal(ptr, 0);

          /* Equilibrium Pm-149 concentration */

          if ((long)RDB[DATA_RESTART_READ_ZERO_SM] == YES)
            val = 0.0;
          else
            val = BARN*prec/RDB[DATA_PM149_DC];

          ptr = (long)RDB[mat + MATERIAL_PTR_PM149_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          sum2 = sum2 + val*vol;

          /* Equilibrium Sm-149 concentration */

          if ((long)RDB[DATA_RESTART_READ_ZERO_SM] == YES)
            val = 0.0;
          else if (abs > 0.0)
            val = BARN*(dir + prec)/abs;
          else
            val = 0.0;

          ptr = (long)RDB[mat + MATERIAL_PTR_SM149_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddStat(val, ptr, 0);

          sum3 = sum3 + val*vol;

          /* Get old concentration */

          old = RDB[iso + COMPOSITION_ADENS];

          /* Set new concentration */

          WDB[iso + COMPOSITION_ADENS] = val;

          /* Adjust atomic density */

          WDB[mat + MATERIAL_ADENS] = RDB[mat + MATERIAL_ADENS] - old + val;

          /* Find nuclide in majorant extra list */

          ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
          while (ptr > VALID_PTR)
            {
              /* Compare pointers */

              if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] ==
                  ((long)RDB[iso + COMPOSITION_PTR_NUCLIDE]))
                {
                  /* Compare fraction */

                  if (val > RDB[ptr + MAJORANT_EXTRA_PTR_NUC])
                    WDB[ptr + MAJORANT_EXTRA_PTR_NUC] = val;

                  /* Break loop */

                  break;
                }

              /* Pointer to next */

              ptr = NextItem(ptr);
            }

          /* Find nuclide in majorant extra list */

          ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
          while (ptr > VALID_PTR)
            {
              /* Compare pointers */

              if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] ==
                  ((long)RDB[iso + COMPOSITION_PTR_NUCLIDE]))
                {
                  /* Compare fraction */

                  if (val > RDB[ptr + MAJORANT_EXTRA_FRAC])
                    WDB[ptr + MAJORANT_EXTRA_FRAC] = val;

                  /* Break loop */

                  break;
                }

              /* Pointer to next */

              ptr = NextItem(ptr);
            }

          /* Check pointer */

          if (ptr < VALID_PTR)
            Die(FUNCTION_NAME, "Nuclide not found in majorant extra list");

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check divisor */

      CheckValue(FUNCTION_NAME, "div", "", div, ZERO, INFTY);

      /* Total Pm-149 and Sm-149 concentrations */

      ptr = (long)RDB[RES_PM149_EQUIL_CONC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(sum2/div, ptr, 0);

      ptr = (long)RDB[RES_SM149_EQUIL_CONC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(sum3/div, ptr, 0);

      /* Total Sm-149 absorption rate */

      ptr = (long)RDB[RES_SM149_ABSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(sum1, ptr, 0);
    }

  /* Update micro depletion and poison densities */

  UpdateMicroDens();

  /***************************************************************************/
}

/*****************************************************************************/
