/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printprecdet.c                                 */
/*                                                                           */
/* Created:       2015/05/11 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Prints contents of precursor detectors                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintPrecDet:"

/*****************************************************************************/

void PrintPrecDet()
{
  long loc0, loc1, ptr, stp, i, j, k, nt, ng, nm;
  char outfile[MAX_STR];
  FILE *fp;

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* Check if output filename is saved */

  if ((long)RDB[loc0 + PRECDET_PTR_OUT_FNAME] < VALID_PTR)
    return;

  /* Higher mpi-tasks do not print detector output */

  if (mpiid > 0)
    return;

  /****************************************************/
  /* Print out the <filename_base>.main file          */
  /****************************************************/

  /* File name */

  sprintf(outfile, "%s.main", GetText(loc0 + PRECDET_PTR_OUT_FNAME));

  /* Open file for writing */

  fp = fopen(outfile, "w");

  /* Get pointer to the detector that tallies the total number */
  /* of live neutron */

  loc1 = (long)RDB[loc0 + PRECDET_PTR_LIVE_DET];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  /* Pointer to statistics */

  ptr = (long)RDB[loc1 + DET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Print live neutron amount and relerr */

  fprintf(fp, "%E %E\n", Mean(ptr, 0, 0), RelErr(ptr, 0, 0));

  /************** Print second line (binning) **********************/

  /* Get number of time bins */

  nt = (long)RDB[loc0 + PRECDET_NT];

  /* Get number of group bins */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Calculate number of mesh bins */

  nm = 1;
  nm = nm*(long)RDB[ptr + MESH_N0];
  nm = nm*(long)RDB[ptr + MESH_N1];
  nm = nm*(long)RDB[ptr + MESH_N2];

  /* Print binning to line */

  fprintf(fp, "%ld %ld %ld %ld %ld\n", nt, ng,
          (long)RDB[ptr + MESH_N0],
          (long)RDB[ptr + MESH_N1],
          (long)RDB[ptr + MESH_N2]);

  /**** Print third line (current time) ****/

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    fprintf(fp, "0.0\n");
  else
    fprintf(fp, "%E\n", RDB[DATA_TIME_CUT_TMAX]);

  /* Get pointer to decay constant array */

  ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Print lambdas of different groups */

  for (j = 0; j < ng; j++)
    fprintf(fp, "%E ", RDB[ptr + j]);

  fprintf(fp,"\n");

  fclose(fp);

  /****************************************************/
  /* Print out the <filename_base>.prec file          */
  /****************************************************/

  /* File name */

  sprintf(outfile, "%s.prec", GetText(loc0 + PRECDET_PTR_OUT_FNAME));

  /* Open file for writing */

  fp = fopen(outfile, "w");

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  /* Get pointer to decay constant array */

  ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over bins to print stats (calculate stable concentrations by dividing with lambda) */

  for (i = 0; i < nt; i++)
    for (j = 0; j < ng; j++)
      for (k = 0; k < nm; k++)
        {
          /* Print stable values */
          /* In criticality source mode we'll divide by lambda as the tallied values */
          /* are production rates */

          if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
            fprintf(fp, "%ld %ld %ld %E %E\n", i, j, k, Mean(stp, i, j, k)/RDB[ptr + j], RelErr(stp, i, j, k));
          else
            fprintf(fp, "%ld %ld %ld %E %E\n", i, j, k, Mean(stp, i, j, k), RelErr(stp, i, j, k));

        }

  /* Close output file */

  fclose(fp);

}

/*****************************************************************************/
