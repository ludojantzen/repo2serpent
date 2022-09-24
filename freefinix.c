/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : freefinix.c                                    */
/*                                                                           */
/* Created:       2014/11/05 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Frees all finix arrays                                       */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "FreeFinix:"

/*****************************************************************************/

void FreeFinix()
{
  long fib;
  Rod *rod;
  Boundary_conditions *bc;
  Scenario *scenario;
  Results *results, *boiresults;
  Options *options;

  /* Check that some finix pins are defined */

  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  /* MPI-id 0 is the task that has allocated these */
  /* let it free them */

  if (mpiid > 0)
    return;

  /* Loop over finix blocks to free arrays from them */

  while (fib > VALID_PTR)
    {
      /* Get Finix pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      boiresults = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS_BOI]);
      bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);
      scenario = (Scenario *)((long)RDB[fib + FINIX_PTR_SCENARIO]);

      /* Free arrays */

      finix_rod_destruct(rod);
      finix_bc_destruct(bc, options);
      finix_scenario_destruct(scenario);
      finix_results_destruct(boiresults, options);
      finix_results_destruct(results, options);
      free(options);

      /* Next finix pin*/

      fib = NextItem(fib);
    }

}

#endif

/*****************************************************************************/
