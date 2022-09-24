/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processfinix.c                                 */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Sets up FINIX data for temperature feedback                  */
/*              Runs all initialization FINIX routines and                   */
/*              prepares all arrays                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"
#include "./FINIX/finix_initialization.h"
#include "./FINIX/aux_functions.h"
#include "./FINIX/initial.h"
#include "./FINIX/finix_output.h"

#define FUNCTION_NAME "ProcessFinix:"

#define FILE_TYPE_ROD      0
#define FILE_TYPE_OPTIONS  1
#define FILE_TYPE_SCENARIO 2

/*****************************************************************************/

void ProcessFinix()
{
  long msh0, i,j, fib, npins;
  double val;
  Rod *rod;
  Boundary_conditions *bc;
  Scenario *scenario;
  Results *results, *boiresults;
  Options *options;
  char **err=NULL;

  /* Check that some finix pins are defined */

  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  /* Close list to allow omp-loops over it later */

  CloseList(fib);

  fprintf(outp, "Processing Finix pins\n");

  /* reset number of pins */

  npins = 0;

  /* Loop over FINIX definitions */

  while (fib > VALID_PTR)
    {

      /* Link power binning here */

      msh0 = (long)RDB[DATA_PTR_FINBIN0];

      while (msh0 > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(msh0 + MESH_PTR_NAME, fib + FINIX_PTR_POWMSH))
            break;

          /* Next possible binning mesh */

          msh0 = NextItem(msh0);
        }

      /* Check that power binning was found */

      if (msh0 < VALID_PTR)
        Error(fib, "Power binning \"%s\" not defined", GetText(fib + FINIX_PTR_POWMSH));

      /* Put direct pointer to mesh */

      WDB[fib + FINIX_PTR_POWMSH] = (double)msh0;

      /* Construct arrays */

      rod = finix_rod_construct();
      bc = finix_bc_construct();
      scenario = finix_scenario_construct();
      results = finix_results_construct();
      boiresults = finix_results_construct();
      options = finix_options_construct();

      /* Store pointers */

      WDB[fib + FINIX_PTR_ROD] = (double)((long)rod);
      WDB[fib + FINIX_PTR_BC] = (double)((long)bc);
      WDB[fib + FINIX_PTR_SCENARIO] = (double)((long)scenario);
      WDB[fib + FINIX_PTR_RESULTS] = (double)((long)results);
      WDB[fib + FINIX_PTR_RESULTS_BOI] = (double)((long)boiresults);
      WDB[fib + FINIX_PTR_OPTIONS] = (double)((long)options);

      /* Create finix_rod.inp for this rod */

      WriteFinixInputFile(fib, FILE_TYPE_ROD);

      /* Create finix_options.inp for this rod */

      WriteFinixInputFile(fib, FILE_TYPE_OPTIONS);

      /* Create finix_scenario.inp for this rod */

      WriteFinixInputFile(fib, FILE_TYPE_SCENARIO);

      /* Initialize data structures */

      err = finix_initialize_data_structures(rod, bc, scenario, results, options);

      /* Initialize data structures */

      err = finix_initialize_data_structures(rod, bc, scenario, boiresults, options);

      /* Handle errors */

      if (err != NULL)
        {
          finix_printf_err(err);
          Die(FUNCTION_NAME, "Error in initialization of FINIX data structures");
        }

      finix_free_err(&err);

      /* Put zero linear power to all axial segments */

      for (i = 0 ; i < options->axial_nodes ; i++)
        bc->linear_power[i] = 0.0;

      /* Set radial nodes to equal-distance, rather than equal area */

      for (i = 0 ; i < options->axial_nodes ; i++)
        {

          /* Pellet */
          for (j = 0; j < options->pellet_radial_nodes ; j++)
            {
              /* Calculate radial node position based on equal distance */

              val = rod->pellet_inner_radius
                + (double)j/(double)(options->pellet_radial_nodes - 1)*
                (rod->pellet_outer_radius - rod->pellet_inner_radius);

              /* EOI values */

              results->radial_node_position[i][j] = val;
              results->radial_node_position_cold[i][j] = val;

              /* BOI values */

              boiresults->radial_node_position[i][j] = val;
              boiresults->radial_node_position_cold[i][j] = val;

              /* Set burnups */
              /*
              boiresults->burnup[i][j] = TFBU[j];
              boiresults->burnup[i][j] = TFBU[j];
              */
            }

          /* Clad */

          for (; j < options->pellet_radial_nodes
                 + options->clad_radial_nodes ; j++)
            {
              /* Calculate radial node position based on equal distance */

              val = rod->clad_inner_radius
                + (double)(j - options->pellet_radial_nodes)/
                (double)(options->clad_radial_nodes - 1)*
                (rod->clad_outer_radius - rod->clad_inner_radius);

              /* EOI values */

              results->radial_node_position[i][j] = val;
              results->radial_node_position_cold[i][j] = val;

              /* BOI values */

              boiresults->radial_node_position[i][j] = val;
              boiresults->radial_node_position_cold[i][j] = val;
            }
        }

      /* Solve initial steady state (HZP) */

      err = finix_solve_initial_steady_state(rod, bc, results, options);

      /* Solve initial steady state (HZP) for BOI results */

      err = finix_solve_initial_steady_state(rod, bc, boiresults, options);

      /* Handle errors */

      if (err != NULL)
        {
          finix_printf_err(err);
          Die(FUNCTION_NAME, "Error in initialization of FINIX data structures");
        }

      finix_free_err(&err);

      npins++;

      /* Process next pin */

      fib = NextItem(fib);

    }

  /* Store total number of pins (needed for interface file) */

  fib = RDB[DATA_PTR_FIN0];

  while(fib > VALID_PTR)
    {
      WDB[fib + FINIX_N_PIN] = (double)npins;
      fib = NextItem(fib);
    }

  /* Set Density Factor usage on */

  WDB[DATA_USE_DENSITY_FACTOR] = (double)YES;

  /* Read initial conditions */

  ReadFinixIFC();

  /* Write pins to an interface file */

  WriteFinixIFC();

  /* Create interface data structure */

  CreateFinixIFC();

}

#endif

/*****************************************************************************/
