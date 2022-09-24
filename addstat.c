/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addstat.c                                      */
/*                                                                           */
/* Created:       2010/11/10 (JLe)                                           */
/* Last modified: 2011/11/30 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Adds score to statistics                                     */
/*                                                                           */
/* Comments: - Toi pointterin haku rutiini on sama kaikille stat funktioille */
/*             syy on se että en tiedä miten vaihtuvapituinen argumentt-     */
/*             välitetään funktiosta toiseen. Sitä voisi miettiä vielä.      */
/*                                                                           */
/*           - Serpent 1:n addstat.c:ssä tulos tallennetaan myös bufferiin   */
/*             tässä ei sitä tehdä.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddStat:"

/*****************************************************************************/

void AddStat(double val, long ptr, ...)
{
  long i, idx, nmax, n, loc0, bins, dim, stp, his;
  va_list argp;
  va_start (argp, ptr);

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

  /* Get pointer to statistics */

  stp = (long)RDB[ptr + SCORE_PTR_DATA] + idx*STAT_BLOCK_SIZE;

  /****************************************************************************/

  /***** Access data **********************************************************/

#ifdef DEBUG

  /* Check value */

  if (!((val > -INFTY) && (val < INFTY)))
    Die(FUNCTION_NAME, "Invalid value %E in %s", val, 
	GetText(ptr + SCORE_PTR_NAME));

#endif

  /* Add value */

  RES1[stp + STAT_N] = RES1[stp + STAT_N] + 1.0;
  RES1[stp + STAT_X] = RES1[stp + STAT_X] + val;
  RES1[stp + STAT_X2] = RES1[stp + STAT_X2] + val*val;
  
  /* Add history data */

  if ((his = (long)RDB[ptr + SCORE_PTR_HIS]) > 0)
    {
      /* Get direct pointer */

      his = his + (idx + nmax*((long)RDB[DATA_CYCLE_IDX]))*STAT_BLOCK_SIZE;

      /* Add data */

      RES1[his + STAT_N] = RES1[stp + STAT_N];
      RES1[his + STAT_X] = RES1[stp + STAT_X];
      RES1[his + STAT_X2] = RES1[stp + STAT_X2];
    }

  /****************************************************************************/
}

/*****************************************************************************/
