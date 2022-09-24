/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : precursorstostore.c                            */
/*                                                                           */
/* Created:       2016/06/04 (VVa)                                           */
/* Last modified: 2018/11/02 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Stores dynamic mode precursors in the end of time interval   */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrecursorsToStore:"

/*****************************************************************************/

void PrecursorsToStore()
{
  long loc0, ptr, src, next, id, nb, tot, flag, tb, nt, ng, nm, nbin, j, k;
  long stp;
  double val;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check that precursors are used */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /*****************************************/
  /* Save mesh-based precursor populations */
  /*****************************************/

  /* Reduce scoring buffers */

  ReduceBuffer();

  /* Get current time bin index */
  /* Plus one since saving the EOI values */

  tb = (long)RDB[DATA_DYN_TB] + 1;

  /* Get number of time bins */

  nt = (long)RDB[loc0 + PRECDET_NT];

  /* Check current time bin index */

  if (tb >= nt)
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

  /* Create store file name*/

  flag = (long)RDB[DATA_EOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storepmeshA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storepmeshB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Get index of batch */

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Open file for writing */

  if (nb == 0)
    fp = fopen(tmpstr, "w");
  else
    fp = fopen(tmpstr, "a");

  /* Write batch number */

  fwrite(&nb, sizeof(long), 1, fp);

  /* Write number of bins */

  fwrite(&nbin, sizeof(long), 1, fp);

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  /* Loop over bins to write bin values to file */

  for (j = 0; j < ng; j++)
    for (k = 0; k < nm; k++)
      {
        /* Get value of the buffer */

        val = BufVal(stp, tb, j, k);

        /* Write current bin value */

        fwrite(&val, sizeof(double), 1, fp);

      }

  /* Close output file */

  fclose(fp);

  /* Return if point-wise precursors don't have to be stored */

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] != PREC_MODE_POINT)
    return;

  /***********************************************************/
  /* Point-wise precursors                                   */

  /* Initialize id */

  id = 0;

  /* Get pointer to precursor source definition */

  src = (long)RDB[DATA_PART_PTR_PSOURCE];
  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Add source size to total number of particles */

  tot = ListSize(src) - 1;

  /* Create store file name*/

  flag = (long)RDB[DATA_EOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storeprecA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storeprecB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Get index of batch */

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Open file for writing */

  if (nb == 0)
    fp = fopen(tmpstr, "w");
  else
    fp = fopen(tmpstr, "a");

  /* Write batch number */

  fwrite(&nb, sizeof(long), 1, fp);

  /* Write number of particles that follows */

  fwrite(&tot, sizeof(long), 1, fp);

  /* Check that source still contains the dummy */

  if ((long)RDB[src + PARTICLE_TYPE] != PARTICLE_TYPE_DUMMY)
    Die(FUNCTION_NAME, "Dummy particle missing from precursor source");

  /* Skip dummy in source */

  next = NextItem(src);

  while (next > VALID_PTR)
    {

      /* Get current */

      ptr = next;

      /* Get next one */

      next = NextItem(ptr);

      /* Write particle to file */

      fwrite(&RDB[ptr + LIST_DATA_SIZE], sizeof(double), (PARTICLE_BLOCK_SIZE -
                                                          LIST_DATA_SIZE), fp);

      /* Decrease number of particles that still needs to be saved */

      tot--;

      /* Remove item from source */

      RemoveItem(ptr);

      /* Put original to stack */

      ToStack(ptr, id++);

      /* Check id */

      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
        id = 0;

    }

  /* Close the store file */

  fclose(fp);

  /* Check that we stored as many particles as we were thinking to */

  if (tot != 0)
    Die(FUNCTION_NAME, "Stored %ld particles too few", tot);

}
