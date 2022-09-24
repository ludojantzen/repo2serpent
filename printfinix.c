/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printfinix.c                                   */
/*                                                                           */
/* Created:       2013/10/12 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Prints some interesting stuff from transient Finix           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "PrintFinix:"

/*****************************************************************************/

void PrintFinix()
{
  long fib, ifc, fpe, fin, ptr;
  long tb, i, j, l, nr, nz, tmax, tme;
  Boundary_conditions *bc;
  Results *results;
  Options *options;
  char *name, outfile[MAX_STR];
  FILE* fp;

  /* Return if no finix rods are defined */

  if((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Get pointer to interface block */

  ifc = (long)RDB[fib + FINIX_PTR_IFC];
  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Print output file name */

  sprintf(outfile, "%s_fin.m", GetText(DATA_PTR_INPUT_FNAME));

  /* Get time bin index  */

  tb = (long)RDB[DATA_DYN_TB];

  /* Open output file */

  if (tb == 0)
    fp = fopen(outfile, "w");
  else
    fp = fopen(outfile, "a");

  /* Get number of time bins */

  tmax = (long)RDB[DATA_DYN_NB];

  /* On the first time step print time data */

  if (tb == 0)
    {
      /* Print time information */

      fprintf(fp, "t = [");

      /* Get pointer to time bin structure */
      /* Does not exist in steady state    */

      if ((tme = (long)RDB[DATA_DYN_PTR_TIME_BINS]) < VALID_PTR)
        {
          /* Put time zero */

          fprintf(fp, "0.0];\n\n");
        }
      else
        {
          /* Get pointer to bins */

          tme = (long)RDB[tme + TME_PTR_BINS];
          CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

          /* Print time interval boundary times */

          for (l = 0; l < tmax + 1; l++)
            fprintf(fp, "%E ", RDB[tme + l]);

          /* Close the array */

          fprintf(fp, "];\n\n");
        }
    }

  /* Loop over pins */

  fpe = (long)RDB[ifc + IFC_PTR_FUEP];

  while (fpe > VALID_PTR)
    {

      /* Get pointer to finix block */

      fin = (long)RDB[fpe + IFC_FUEP_PTR_FINIX];

      /* Get FINIX option and bc pointers */

      options = (Options *)((long)RDB[fin + FINIX_PTR_OPTIONS]);
      bc = (Boundary_conditions *)((long)RDB[fin + FINIX_PTR_BC]);

      /* Get pointer to pin names */

      ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get first uni name and use it for this pin */

      name = GetText(ptr);

      /* At the first time interval initialize arrays */
      /* Makes it faster to run the .m file in MATLAB */

      if (tb == 0)
        {
          /* Get number of axial nodes */

          nz = options->axial_nodes;

          /* Get number of radial nodes */

          nr = options->pellet_radial_nodes + options->clad_radial_nodes;

          /* Print initialization command */

          fprintf(fp, "ROD_%s_T = zeros(%ld, %ld, %ld);\n", name, tmax + 1, nz, nr);
          fprintf(fp, "ROD_%s_P = zeros(%ld, %ld, %ld);\n", name, tmax + 1, nz, nr);
          fprintf(fp, "ROD_%s_RC = zeros(%ld, %ld, %ld);\n", name, tmax + 1, nz, nr);
          fprintf(fp, "ROD_%s_RH = zeros(%ld, %ld, %ld);\n", name, tmax + 1, nz, nr);
          fprintf(fp, "ROD_%s_HFlux = zeros(%ld, %ld);\n\n", name, tmax + 1, nz);
        }

      /* For first time interval also write the BOI data */

      if (tb == 0)
        l = 0;
      else
        l = 1;

      for (; l <= 1; l++)
        {

          /*****************************************/
          /* Print temperature array starting part */
          /*****************************************/

          fprintf(fp, "ROD_%s_T(%ld,:,:) = [\n", name, tb + l + 1);

          /* Get Finix pointers for the time step*/

          if (l == 0)
            results = (Results *)((long)RDB[fin + FINIX_PTR_RESULTS_BOI]);
          else
            results = (Results *)((long)RDB[fin + FINIX_PTR_RESULTS]);

          /* Loop over axial segments */

          for (i = 0; i < options->axial_nodes; i++)
            {
              /* Loop over radial nodes to print temperature */

              for ( j = 0 ; j < options->pellet_radial_nodes +
                      options->clad_radial_nodes ; j++)
                fprintf(fp,"%E ",results->temperature[i][j]);

              fprintf(fp,"\n");
            }
          fprintf(fp,"];\n\n");

          /*******************************************/
          /* Print power density array starting part */
          /*******************************************/

          fprintf(fp, "ROD_%s_P(%ld,:,:) = [\n", name, tb + l + 1);

          /* Loop over axial segments */

          for (i = 0; i < options->axial_nodes; i++)
            {
              /* Loop over radial nodes to print power density */

              for ( j = 0 ; j < options->pellet_radial_nodes +
                      options->clad_radial_nodes ; j++)
                fprintf(fp,"%E ",bc->linear_power[i]/
                        (PI*results->radial_node_position[i][options->pellet_radial_nodes]*
                         results->radial_node_position[i][options->pellet_radial_nodes])*
                        bc->power_distr[i][j]/1E6);

              fprintf(fp,"\n");
            }

          fprintf(fp,"];\n\n");

          /****************************************/
          /* Print hot radius array starting part */
          /****************************************/

          fprintf(fp,"ROD_%s_RH(%ld,:,:) = [\n", name, tb + l + 1);

          /* Loop over axial segments */

          for (i = 0; i < options->axial_nodes; i++)
            {
              /* Loop over radial nodes  */

              for ( j = 0 ; j < options->pellet_radial_nodes +
                      options->clad_radial_nodes ; j++)
                fprintf(fp,"%E ",results->radial_node_position[i][j]);

              fprintf(fp,"\n");
            }

          fprintf(fp,"];\n\n");

          /*****************************************/
          /* Print cold radius array starting part */
          /*****************************************/

          fprintf(fp, "ROD_%s_RC(%ld,:,:) = [\n", name, tb + l + 1);

          /* Loop over axial segments */

          for (i = 0; i < options->axial_nodes; i++)
            {
              /* Loop over radial nodes */

              for ( j = 0 ; j < options->pellet_radial_nodes +
                      options->clad_radial_nodes ; j++)
                fprintf(fp,"%E ",results->radial_node_position_cold[i][j]);

              fprintf(fp,"\n");
            }

          fprintf(fp,"];\n\n");

          /***************************************/
          /* Print heat flux array starting part */
          /***************************************/

          /* Open array */

          fprintf(fp, "ROD_%s_HFlux(%ld,:) = [\n", name, tb + l + 1);

          /* Loop over axial segments */

          for (i = 0; i < options->axial_nodes; i++)
            {
              fprintf(fp,"%E ",bc->heat_flux[i]);
            }

          fprintf(fp,"];\n\n");

        }

      /* Next rod */

      fpe = NextItem(fpe);
    }

  fclose(fp);

}

#endif

/*****************************************************************************/
