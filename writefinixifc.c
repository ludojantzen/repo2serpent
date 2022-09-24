/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writefinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Creates an IFC-template for the ReadInterface()              */
/*              for the initialization of the FINIX interface                */
/*                                                                           */
/* Comments:   -Tän vois periaatteessa kirjoittaa suoraan muistiinkin        */
/*              mutta lienee järkevintä luoda kaikki rajapinnat vasta        */
/*              ReadInterface():ssa                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "WriteFinixIFC:"

/*****************************************************************************/

void WriteFinixIFC()
{
  long fib, msh, nz, i, j, k, nnt;
  double zmin, zmax, nt;
  char tmpstr[MAX_STR];
  Rod *rod;
  Options *options;
  Results *results;
  FILE *fp;

  /* Print interface file name */

  sprintf(tmpstr,"./Finix.ifc");

  /* Try to open the interface file for writing */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Could not open file \"%s\" for writing", tmpstr);

  fib = (long)RDB[DATA_PTR_FIN0];

  /* Write the interface file header */

  /* Output filename*/
  sprintf(tmpstr,"Finixifcout.m");

  /* Interface type, outputfile, number of pins */

  fprintf(fp,"6 %s %ld\n", tmpstr, (long)RDB[fib + FINIX_N_PIN]);

  /* Loop over pins */

  while (fib > VALID_PTR)
    {
      /* Get FINIX pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      msh = (long)RDB[fib + FINIX_PTR_POWMSH];

      /* Pin universe */

      fprintf(fp,"%s\n",GetText(fib + FINIX_PTR_UNI_NAME));

      /* Power tally */

      fprintf(fp,"%ld %e %e 1 0 360 %ld %f %f\n", (long)RDB[msh + MESH_N2],
              RDB[msh + MESH_MIN2], RDB[msh + MESH_MAX2],
              (long)RDB[msh + MESH_N0], RDB[msh + MESH_MIN0],
              RDB[msh + MESH_MAX0]);


      /* Fast flux tally (only one radial bin) */

      fprintf(fp,"%ld %f %f 1 0 360 1 %f %f 1 15\n", (long)RDB[msh + MESH_N2],
              RDB[msh + MESH_MIN2], RDB[msh + MESH_MAX2], RDB[msh + MESH_MIN0],
              100*(rod->clad_outer_radius));

      /* Number of FINIX timesteps in each transport time step */

      nnt = 1;

      /* Get number of timesteps */

      if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
        nt = (double)(RDB[DATA_DYN_NB]*nnt);
      else
        nt = (double)1;

      /* Get number of axial zones */

      nz = (long)(options->axial_nodes);

      /* Write number of axial zones */

      fprintf(fp, "%ld\n", nz);

      /* Loop over the axial segments */

      for (i = 0 ; i < nz ; i++)
        {
          /* Get zmin and zmax of this segment */

          zmin = RDB[msh + MESH_MIN2] +
            i*(RDB[msh + MESH_MAX2]-
               RDB[msh + MESH_MIN2])/(double)nz;
          zmax = RDB[msh + MESH_MIN2] +
            (i+1)*(RDB[msh + MESH_MAX2]-
                   RDB[msh + MESH_MIN2])/(double)nz;

          /* Put limits of axial segment (only 1 angular zone) */

          fprintf(fp,"%f %f 1\n", zmin, zmax);

          /* Put limits of angular zone and number of radial nodes */

          if (rod->pellet_inner_radius > 0.0)
            {
              /* Print zoning */

              fprintf(fp,"0 360 %ld\n",
                      (long)(options->pellet_radial_nodes
                             + options->clad_radial_nodes + 1));

              /* Print central hole temperature*/

              fprintf(fp,"%E %E %E\n", 0.0, 0.0, results->temperature[i][0]);

            }
          else
            fprintf(fp,"0 360 %ld\n",
                    (long)(options->pellet_radial_nodes
                           + options->clad_radial_nodes));

          /* Print rest of the nodes */

          for(j=0; j < options->pellet_radial_nodes
                + options->clad_radial_nodes; j++)
            fprintf(fp,"%E %E %E\n", results->radial_node_position_cold[i][j]*100,
                    results->radial_node_position[i][j]*100, results->temperature[i][j]);

        }

      fib = NextItem(fib);
    }

  fclose(fp);
}

#endif

/*****************************************************************************/
