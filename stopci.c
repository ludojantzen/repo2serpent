/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stopci.c                                       */
/*                                                                           */
/* Created:       2011/05/20 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Check for stopping criterion in corrector iteration          */
/*                                                                           */
/* Comments: Called from burnupcycle                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopCI:"

/*****************************************************************************/

void StopCI()
{

  long mat, dep;
  double sc, vol, ide, nom, denom, axs;

  /* If maximum number of iterations is given */

  if((long)RDB[DATA_BURN_CI_I] >= (long)RDB[DATA_BURN_CI_MAXI]-1)
    WDB[DATA_BURN_CI_LAST] = YES;
  else
    WDB[DATA_BURN_CI_LAST] = NO;

  /* The convergence criterion works at least for SIE, but maybe not for */
  /* other corrector iteration methods                                   */

  if ((long)RDB[DATA_BURN_SIE] == NO)
    return;

  /* if the previous step was a predictor step */
  if (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    {
      /* Store values from first predictor */
      if((long)RDB[DATA_BURN_STEP] == 0)
        {

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {

              WDB[mat + MATERIAL_CI_BOS_ABSXS] = RDB[mat + MATERIAL_CI_AVE_ABSXS];
              WDB[mat + MATERIAL_CI_AVE_ABSXS] = 0.0;

              WDB[mat + MATERIAL_CI_BOS_ABSXS2] = RDB[mat + MATERIAL_CI_AVE_ABSXS2];
              WDB[mat + MATERIAL_CI_AVE_ABSXS2] = 0.0;

              mat = NextItem(mat);
            }


        }
    }

  /* Loop over materials to calculate stopping criterion */

  mat = (long)RDB[DATA_PTR_M0];

  /* Reset nominator and denominator */

  nom = 0.0;
  denom = 0.0;

  while (mat > VALID_PTR)
    {

      /* Check burn flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /* Get material volume */
          vol = RDB[mat + MATERIAL_VOLUME];

          if (!(vol > 0))
            Die(FUNCTION_NAME, "Invalid material volume");

          /* Get material ideal deviation */
          ide = RDB[mat + MATERIAL_CI_IDE];

          if (ide < 0)
            Die(FUNCTION_NAME, "Invalid material ideal deviation");

          /* Get material average ABSXS */

          axs = RDB[mat + MATERIAL_CI_AVE_ABSXS];

          if (axs < 0)
            Die(FUNCTION_NAME, "Invalid material ABSXS");

          /* Add to numerator and denominator */

          nom = nom + vol*ide;

          denom = denom + vol*axs;

        }

      mat = NextItem(mat);
    }

  /* Calculate stopping criterion */

  sc = nom/denom;

  /* Reset stopping flag */

  WDB[DATA_BURN_CI_LAST] = NO;

  /* Check if ideal deviation is below tolerance */

  if(sc < RDB[DATA_BURN_CI_TOLER])
    {
      /* This was last iteration */

      WDB[DATA_BURN_CI_LAST] = YES;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {

          /* Store BOS absorption cross sections */
          /* Used in updatecistop.c              */

          WDB[mat + MATERIAL_CI_BOS_ABSXS] = RDB[mat + MATERIAL_CI_AVE_ABSXS];
          WDB[mat + MATERIAL_CI_AVE_ABSXS] = 0.0;

          WDB[mat + MATERIAL_CI_BOS_ABSXS2] = RDB[mat + MATERIAL_CI_AVE_ABSXS2];
          WDB[mat + MATERIAL_CI_AVE_ABSXS2] = 0.0;

          /* Get pointer to transmutation list */
          dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];

          while (dep > VALID_PTR)
            {

              /* Store microscopic average cross section for next step */

              WDB[dep + DEP_TRA_AV0] = RDB[dep + DEP_TRA_AV1];

              dep = NextItem(dep);
            }

          /* Get pointer to fission list */
          dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];

          while (dep > VALID_PTR)
            {

              /* Store microscopic average cross section for next step */

              WDB[dep + DEP_TRA_AV0] = RDB[dep + DEP_TRA_AV1];

              dep = NextItem(dep);
            }


          mat = NextItem(mat);
        }


    }

  /* Check if next iteration would exeed maximum iterations */

  if ((long)RDB[DATA_BURN_CI_I] + 1 >= (long)RDB[DATA_BURN_CI_MAXI])
    {
      /* This was last iteration */

      WDB[DATA_BURN_CI_LAST] = YES;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {

          /* Store BOS absorption cross sections */
          /* Used in updatecistop.c              */

          WDB[mat + MATERIAL_CI_BOS_ABSXS] = RDB[mat + MATERIAL_CI_AVE_ABSXS];
          WDB[mat + MATERIAL_CI_AVE_ABSXS] = 0.0;

          WDB[mat + MATERIAL_CI_BOS_ABSXS2] = RDB[mat + MATERIAL_CI_AVE_ABSXS2];
          WDB[mat + MATERIAL_CI_AVE_ABSXS2] = 0.0;

          /* Get pointer to transmutation list */
          dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];

          while (dep > VALID_PTR)
            {

              /* Store microscopic average cross section for next step */

              WDB[dep + DEP_TRA_AV0] = RDB[dep + DEP_TRA_AV1];

              dep = NextItem(dep);
            }

          /* Get pointer to fission list */
          dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];

          while (dep > VALID_PTR)
            {

              /* Store microscopic average cross section for next step */

              WDB[dep + DEP_TRA_AV0] = RDB[dep + DEP_TRA_AV1];

              dep = NextItem(dep);
            }

          mat = NextItem(mat);
        }
    }

}

/*****************************************************************************/
