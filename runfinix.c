/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : runfinix.c                                     */
/*                                                                           */
/* Created:       2013/03/27 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Runs Finix for steady state or transient                     */
/*                                                                           */
/* Comments:  - Transientti pitää vielä kirjoittaa                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"
#include "./FINIX/initial.h"
#include "./FINIX/transient.h"
#include "./FINIX/aux_functions.h"

#define FUNCTION_NAME "RunFinix:"

/*****************************************************************************/

void RunFinix(long fib, long fpe)
{
  long nz, nr, i, j;
  double tming, tmaxg;
  Rod *rod;
  Boundary_conditions *bc;
  Results *results, *boiresults;
  Options *options;
  char **err=NULL;

  /* Check steady state */

  if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {

      /* Criticality source simulation */

      /* Get Finix pointers for steady state */

      /* Get steady state Finix pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);

      /* Solve steady state solution */

      fprintf(outp, "Solving FINIX steady state for rod %s\n",
              GetText(fib + FINIX_PTR_UNI_NAME));

      err = finix_solve_initial_steady_state(rod, bc, results, options);

      /* Handle errors */

      if (err != NULL)
        {
          finix_printf_err(err);
          Die(FUNCTION_NAME, "Error in executing FINIX");
        }

      finix_free_err(&err);

    }
  else
    {

      /* Time dependent simulation */

      /* Get transport time bin limits */

      tming = RDB[DATA_TIME_CUT_TMIN];

      /* Get transport time bin limits */

      tmaxg = RDB[DATA_TIME_CUT_TMAX];

      /* Get Finix pointers for the time step*/

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      boiresults = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS_BOI]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);

      /* Set temperature distribution at beginning of step to equal that */
      /* of the end of previous step */

      nz = options->axial_nodes;
      nr = options->pellet_radial_nodes + options->clad_radial_nodes;

      for (i = 0 ; i < nz ; i++)
        for (j = 0 ; j < nr ; j++)
          results->temperature[i][j] = boiresults->temperature[i][j];

      /* Set up boundary conditions */

      /* Set beginning time */

      bc->time = tming;

      /* Set timestep length */

      bc->dt = tmaxg - tming;

      /* Boundary conditions should be updated before runfinix.c if needed */

      /* Run Finix for the current step */

      fprintf(outp,"Solving FINIX transient from %E to %E (%E s)\n", tming, tmaxg, tmaxg - tming);

      /* Solve transient */

      finix_solve_transient(bc->dt, rod, bc, results, options);

      /* Handle errors */

      if (err != NULL)
        {
          finix_printf_err(err);
          Die(FUNCTION_NAME, "Error in transient FINIX solution");
        }

      finix_free_err(&err);


    }
}

#endif

/*****************************************************************************/
