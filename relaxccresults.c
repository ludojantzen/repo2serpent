/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : relaxccresults.c                               */
/*                                                                           */
/* Created:       2019/02/11 (VVa)                                           */
/* Last modified: 2019/02/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Relaxes power and transmutation cross sections for coupled   */
/*              calculations.                                                */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RelaxCCResults:"

/*****************************************************************************/

void RelaxCCResults()
{
  long id, mat;

  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == NO)
    return;

  /* Collect results from MPI tasks in dynamic mode                       */
  /* (In crit mode, this subroutine is never called from transportcycle ) */

  if (RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
    CollectParallelData();

  /* Calculate relaxation factor alpha */

  CalculateRelAlpha();

  /* Calculate population size for next iteration */

  CalculateRelPopSize();

  /* Relax the tallied power distribution */

  RelaxInterfacePower();

  /* Intra-step relaxation of tallied transmutation cross sections */
  /* We don't want to overwrite the corrector transmuxs on further */
  /* predictors */

  if ((RDB[DATA_BURN_SIE] == (double)NO) ||
      (!((RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) &&
         (RDB[DATA_BURN_STEP] != 0.0))))
    {
      /* Reduce private results */

      ReducePrivateRes();

      /* Get OpenMP id */

      id = OMP_THREAD_NUM;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];

      while (mat > VALID_PTR)
        {
          /* Check burn flag */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            {

              /* Calculate momentary transmutation XS*/

              CalculateTransmuXS(mat, id);

              /* Relax transmutation XS inside burnup step */

              RelaxTransmuXS(mat, id);
            }

          /* Next material */

          mat = NextItem(mat);
        }
    }
}
