/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : finalizeccprecdet.c                            */
/*                                                                           */
/* Created:       2016/06/08 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads final mesh-based precursor values from store file and  */
/*              writes them to the buffer and adds the statistic             */
/*                                                                           */
/* Comments:   -Needed because mesh-based precursors are cleared after each  */
/*             batch in coupled calculation mode                             */
/*             -Needed if we want to print out the final mesh based values   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FinalizeCCPrecDet:"

/*****************************************************************************/

void FinalizeCCPrecDet()
{
  long ptr, nt, nm, ng, nbin, maxb, j, k, tb;
  long nb, stp, loc0;
  long tmplong, flag;
  char tmpstr[MAX_STR];
  double val;
  FILE *fp;

  /***************************************************************************/

  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* This is only required in coupled transient calculations */

  if (RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN)
    return;

  /* Get current time bin index */
  /* Plus one since reading the EOI values */

  tb = (long)RDB[DATA_DYN_TB] + 1;

  /* Get number of time bins */

  nt = (long)RDB[loc0 + PRECDET_NT];

  /* Check current time bin index */
  /* Should be equal since it is the last bin */

  if (tb != nt - 1)
    Die(FUNCTION_NAME, "Current time interval %ld, number of time bins %ld", tb, nt);

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Calculate number of mesh bins */

  nm = 1;
  nm = nm*(long)RDB[ptr + MESH_N0];
  nm = nm*(long)RDB[ptr + MESH_N1];
  nm = nm*(long)RDB[ptr + MESH_N2];

  /* Group bins */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Calculate number of bins to save */

  nbin = nm*ng;

  /* Create store file name */
  /* Now reading EOI values */

  flag = (long)RDB[DATA_EOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storepmeshA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storepmeshB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Get index of batch */

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Open file for reading */

  if ((fp = fopen(tmpstr, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", tmpstr);

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  /* Get maximum number of batches */

#ifdef MPI_MODE1
  maxb = (long)RDB[DATA_SRC_BATCHES];
#else
  maxb = (long)(RDB[DATA_SRC_BATCHES]/(double)mpitasks);
#endif

  /* Loop over all batches */

  for (nb = 0; nb < maxb; nb++)
    {

      /* Clear buffers in order to add next batch */

      ClearBuf();

      /* Read batch number */

      if (fread(&tmplong, sizeof(long), 1, fp) != 1)
        Die(FUNCTION_NAME, "Could not read batch number from file %s", tmpstr);

      /* Check batch number */

      if (nb != tmplong)
        Die(FUNCTION_NAME, "Wrong batch number read from mesh file %ld (expected %ld)",
            tmplong, nb);

      /* Read number of bins */

      if (fread(&tmplong, sizeof(long), 1, fp) != 1)
        Die(FUNCTION_NAME, "Could not read number of bins from file %s", tmpstr);

      /* Check number of bins */

      if (nbin != tmplong)
        Die(FUNCTION_NAME, "Wrong number of bins read from mesh file %ld (expected %ld)",
            tmplong, nbin);

      /* Loop over bins to read bin values */

      for (j = 0; j < ng; j++)
        for (k = 0; k < nm; k++)
          {
            /* Read bin value from file */

            if (fread(&val, sizeof(double), 1, fp) != 1)
              Die(FUNCTION_NAME, "Could not read bin value from %s", tmpstr);

            /* Put bin value to the BOI bin of current time interval */
            /*            printf("Adding to buffer %E\n", val);*/
            AddBuf(val, 1.0, stp, 0, -1, tb, j, k);

          }

      /* Reduce buffers for adding statistics */

      ReduceBuffer();

      /* Add statistics from this batch */

      /* Loop over bins to read bin values */

      for (j = 0; j < ng; j++)
        for (k = 0; k < nm; k++)
          {
            /* Get bin value */

            val = BufVal(stp, tb, j, k);

            /* Add statistic */
            /*            printf("Adding value %E\n", val);*/
            AddStat(val, stp, tb, j, k);

          }

      /* Handle next bin or end for loop*/

    }

  /* Close mesh file */

  fclose(fp);

  /***************************************************************************/
}

/*****************************************************************************/
