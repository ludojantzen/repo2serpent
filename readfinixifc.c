/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readfinixifc.c                                 */
/*                                                                           */
/* Created:       2015/03/11 (VVa)                                           */
/* Last modified: 2018/11/02 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads initial guess from interface file for FINIX            */
/*                                                                           */
/* Comments:   -If the interface file does not correspond to the FINIX pin   */
/*              definitions, it won't be read                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "ReadFinixIFC:"

/*****************************************************************************/

void ReadFinixIFC()
{
  long fib, i, j, k, maxk;
  long dummy, pins;
  long nz, na, nr;
  double zmin, zmax, amin, amax, rmin, rmax;
  double T, rh, rc;
  double dbl0, dbl1;
  char tmpstr[MAX_STR];
  Rod *rod;
  Options *options;
  Results *results, *boiresults;
  FILE *fp;

  /* Return if no initialization file is given */

  if ((long)RDB[DATA_PTR_FININIT_FNAME] < VALID_PTR)
    return;

  /* Print interface file name */

  sprintf(tmpstr, "%s", GetText(DATA_PTR_FININIT_FNAME));

  /* Try to open the interface file for writing */
  /* If file cannot be read, just return        */
  /* HZP initial conditions will be used        */

  if ((fp = fopen(tmpstr, "r")) == NULL)
    Die(FUNCTION_NAME, "Could not open initial conditions file \"%s\" for reading", tmpstr);

  fprintf(outp, "Reading initial conditions for FINIX pins from \"%s\"\n", tmpstr);

  /* Read interface type */

  if (fscanf(fp, "%ld", &dummy) == EOF)
    Die(FUNCTION_NAME, "Could not read interface type");

  /* Check that interface type is 6 */

  if (dummy != IFC_TYPE_FPIP)
    Die(FUNCTION_NAME, "Interface type was not %ld (was %ld)", IFC_TYPE_FPIP, dummy);

  /* Get output file name (not checked) */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "Could not read output filename");

  /* Get number of pins in interface */

  if (fscanf(fp, "%ld", &dummy) == EOF)
    Die(FUNCTION_NAME, "Could not read number of pins");

  /* Count number of finix pins */

  fib = (long)RDB[DATA_PTR_FIN0];
  pins = 0;

  while (fib > VALID_PTR)
    {
      /* Increment pin counter */

      pins++;

      /* Get next finix block*/

      fib = NextItem(fib);
    }

  /* Check that the number of pins in interface file corresponds */
  /* to the number of FINIX-pins in memory                       */

  if (dummy != pins)
    Die(FUNCTION_NAME, "Number of pins does not match: IFC: %ld, FINIX: %ld", dummy, pins);

  /* Read fuel rod information */

  for (i = 0; i < pins; i++)
    {
      /* Read pin universe */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      /* Find finix block from memory */

      fib = (long)RDB[DATA_PTR_FIN0];
      CheckPointer(FUNCTION_NAME, "(fib)", DATA_ARRAY, fib);

      /* Find correct pin */

      while (fib > VALID_PTR)
        {
          /* Compare ifc pin name to FINIX pin name */

          if(!strcmp(tmpstr, GetText(fib + FINIX_PTR_UNI_NAME)))
            break;

          /* Next FINIX pin */

          fib = NextItem(fib);

        }

      /* Check if found */

      if (fib < VALID_PTR)
        Die(FUNCTION_NAME, "Interface file contains unknown pin %s", tmpstr);

      /* Get FINIX pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      boiresults = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS_BOI]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);

      /* Read power output limits (not cheked) */

      if (fscanf(fp, "%ld %lf %lf %ld %lf %lf %ld %lf %lf",
                &nz, &zmin, &zmax, &na, &amin, &amax,
                &nr, &rmin, &rmax) == EOF)
        Die(FUNCTION_NAME, "Could not read power binning");

      /* Read flux output limits (not checked) */

      if (fscanf(fp, "%ld %lf %lf %ld %lf %lf %ld %lf %lf",
                &nz, &zmin, &zmax, &na, &amin, &amax,
                &nr, &rmin, &rmax) == EOF)
        Die(FUNCTION_NAME, "Could not read flux binning");

      /* Read flux energy limits (not checked) */

      if (fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
        Die(FUNCTION_NAME, "Could not read flux energy limits");

      /* Read number of axial zones */

      if (fscanf(fp, "%ld", &nz) == EOF)
        Die(FUNCTION_NAME, "Could not read number of axial zones");

      /* Check number of axial zones */

      if (nz != options->axial_nodes)
        Die(FUNCTION_NAME,
            "Number of axial zones doesn't match IFC: %ld vs FINIX: %d",
            nz, options->axial_nodes);

      /* Loop over axial zones */

      for (j = 0; j < options->axial_nodes ; j++)
        {
          /* Read axial limits (not checked) */

          if (fscanf(fp, "%lf %lf", &dbl0, &dbl1) == EOF)
            Die(FUNCTION_NAME, "Could not read axial limits");

          /* Read number of angular zones */

          if (fscanf(fp, "%ld", &dummy) == EOF)
            Die(FUNCTION_NAME, "Could not read number of angular zones");

          /* Get angular limits of the zone */

          if (fscanf(fp, "%lf %lf", &dbl0, &dbl1) == EOF)
            Die(FUNCTION_NAME, "Could not read angular limits");

          /* Read number of radial points */

          if (fscanf(fp, "%ld", &dummy) == EOF)
            Die(FUNCTION_NAME, "Could not read radial points");

          /* Calculate number of supposed radial points */

          if (rod->pellet_inner_radius > 0.0)
            {
              /* Rod with central hole */

              nr = options->pellet_radial_nodes
                + options->clad_radial_nodes + 1;

            }
          else
            {

              /* Rod without central hole */

              nr = options->pellet_radial_nodes
                + options->clad_radial_nodes;

            }

          /* Check number of radial points in file */

          if (nr != dummy)
            Die(FUNCTION_NAME,
                "Number of radial points doesn't match IFC: %ld FINIX: %ld\n",
                dummy, nr);

          /* Read and discard central hole data */

          if (rod->pellet_inner_radius > 0.0)
            if (fscanf(fp, "%lf %lf %lf", &rc, &rh, &T) == EOF)
              Die(FUNCTION_NAME, "Could not first radial node");

          maxk = options->pellet_radial_nodes + options->clad_radial_nodes;

          /* Loop over radial points and get hot radius and temperature */

          for (k = 0; k < maxk ; k++)
            {

              /* Read cold radius, hot radius and temperature */

              if (fscanf(fp, "%lf %lf %lf", &rc, &rh, &T) == EOF)
                Die(FUNCTION_NAME, "Could not read node data");

              /* Check cold radius */
              /* TODO: Ehkä joku suhteellinen vertailu */

              if (fabs(results->radial_node_position_cold[j][k]*100 - rc) > 1E-3)
                Die(FUNCTION_NAME, "Something wrong with radial node position %f vs %f\n",
                    results->radial_node_position_cold[j][k]*100, rc);

              /* Store hotradius */

              boiresults->radial_node_position[j][k] = rh/100.0;
              results->radial_node_position[j][k] = rh/100.0;

              /* Store temperature */

              boiresults->temperature[j][k] = T;
              results->temperature[j][k] = T;

            }

        }
    }

  fclose(fp);

  fprintf(outp, "Ok.\n");

  return;
}

#endif

/*****************************************************************************/
