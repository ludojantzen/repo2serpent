/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hisval.c                                       */
/*                                                                           */
/* Created:       2011/05/18 (JLe)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Returns value stored in history array of scores              */
/*                                                                           */
/* Comments: - Toi pointterin haku rutiini on sama kaikille stat funktioille */
/*             syy on se että en tiedä miten vaihtuvapituinen argumentt-     */
/*             välitetään funktiosta toiseen. Sitä voisi miettiä vielä.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "HisVal:"

/*****************************************************************************/

double HisVal(long ptr, long cycle, ...)
{
  long i, idx, nmax, n, loc0, bins, dim, stp, stp0, stp1, skip;
  double val;
  va_list argp;
  va_start (argp, cycle);

  /***************************************************************************/

  /***** Get pointer to statistics *******************************************/

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "", DATA_ARRAY, ptr);

  /* Get dimension */

  dim = (long)RDB[ptr + SCORE_DIM];
  CheckValue(FUNCTION_NAME, "dim", "", dim, 1, 10);

  /* Reset bin index and maximum */

  idx = 0;
  nmax = 1;

  /* Pointer to list of maximum values */

  loc0 = (long)RDB[ptr + SCORE_PTR_NMAX];

  /* Loop over dimensions */

  for (n = 0; n < dim; n++)
    {
      /* Get number of bins */

      bins = (long)RDB[loc0++];

      /* Get bin index */

      i = va_arg(argp, int);

      /* Check index */

#ifdef DEBUG

      if ((i < 0) || (i > bins - 1))
        Die(FUNCTION_NAME,
            "Invalid bin index %ld in %s at level %ld (max = %ld)",
            i, GetText(ptr + SCORE_PTR_NAME), n, bins - 1);

#endif

      /* Update index */

      idx = idx + i*nmax;

      /* Update maximum */

      nmax = nmax*bins;
    }

#ifdef DEBUG

  /* Check index and maximum number of bins */

  if (idx > nmax - 1)
    Die(FUNCTION_NAME, "idx > nmax in %s", GetText(ptr + SCORE_PTR_NAME));
  if (idx > (long)RDB[ptr + SCORE_STAT_SIZE] - 1)
    Die(FUNCTION_NAME, "idx > max in %s", GetText(ptr + SCORE_PTR_NAME));
  if (nmax > (long)RDB[ptr + SCORE_STAT_SIZE])
    Die(FUNCTION_NAME, "nmax > max in %s", GetText(ptr + SCORE_PTR_NAME));

#endif

  /* Get pointer to history data */

  if ((stp = (long)RDB[ptr + SCORE_PTR_HIS]) < 1)
    return 0.0;

  /* Get direct pointers */

  stp0 = stp + (idx + nmax*(cycle - 1))*STAT_BLOCK_SIZE;
  stp1 = stp + (idx + nmax*cycle)*STAT_BLOCK_SIZE;

  /* Get number of inactive cycles if using fission source passing */

  if((RDB[DATA_USE_FSP] == (double)NO) ||
     ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
      && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)))
    {
      /* Set number of inactive batches */

      skip = (long)RDB[DATA_CRIT_SKIP];
    }
  else
    {
      /* Set number of inactive batches */

      skip = (long)RDB[DATA_FSP_CRIT_SKIP];
    }

  /****************************************************************************/

  /***** Access data **********************************************************/

  /* Get value (TODO: tässä pitää ottaa huomioon MPI taskit) */
  /* Use cycle-wise value for first inactive and first active cycle */

  if ((cycle == (long)RDB[DATA_CRIT_SKIP] - skip) ||
      (cycle == (long)RDB[DATA_CRIT_SKIP]))
    val = RES1[stp1 + STAT_X];
  else
    val = RES1[stp1 + STAT_X] - RES1[stp0 + STAT_X];

  /* Return value */

  return val/((double)mpitasks);

  /****************************************************************************/
}

/*****************************************************************************/
