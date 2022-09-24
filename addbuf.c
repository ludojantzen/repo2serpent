/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addbuf.c                                       */
/*                                                                           */
/* Created:       2010/11/10 (JLe)                                           */
/* Last modified: 2014/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Adds results to score buffer                                 */
/*                                                                           */
/* Comments: - Toi pointterin haku rutiini on sama kaikille stat funktioille */
/*             syy on se että en tiedä miten vaihtuvapituinen argumentt-     */
/*             välitetään funktiosta toiseen. Sitä voisi miettiä vielä.      */
/*                                                                           */
/*           - Indeksin voi antaa suoraan 1-dimensioisssa taulukoissa        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddBuf:"

/*****************************************************************************/

void AddBuf(double val, double wgt, long ptr, long id, long idx, ...)
{
  long i, nmax, n, loc0, bins, dim, sz;
  va_list argp;
  va_start (argp, idx);

  /***************************************************************************/

  /***** Get pointer to statistics *******************************************/

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "", DATA_ARRAY, ptr);

  /* Check id */

  CheckValue(FUNCTION_NAME, "id", "", id, 0, MAX_OMP_THREADS);

  /* Check direct index */

  if (idx < 0)
    {
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

      /* Check if access is allowed */

      if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
        Die(FUNCTION_NAME, "BUF array not ready for access");

      /* Check index and maximum number of bins */
      
      if (idx > nmax - 1)
        Die(FUNCTION_NAME, "idx > nmax in %s", GetText(ptr + SCORE_PTR_NAME));
      if (idx > (long)RDB[ptr + SCORE_STAT_SIZE] - 1)
        Die(FUNCTION_NAME, "idx > max in %s", GetText(ptr + SCORE_PTR_NAME));
      if (nmax > (long)RDB[ptr + SCORE_STAT_SIZE])
        Die(FUNCTION_NAME, "nmax > max in %s", GetText(ptr + SCORE_PTR_NAME));
      
#endif      
    }

  /****************************************************************************/

  /***** Access data **********************************************************/

#ifdef DEBUG

  /* Check value */

  if (!((val > -INFTY) && (val < INFTY)))
    Die(FUNCTION_NAME, "Invalid value %E in %s", val, 
        GetText(ptr + SCORE_PTR_NAME));

  /* Check weight */

  if (!((wgt > 0.0) && (wgt < INFTY)))
    Die(FUNCTION_NAME, "Invalid weight %E in %s", wgt, 
        GetText(ptr + SCORE_PTR_NAME));

#endif

  /* Get pointer to buffer */

  loc0 = (long)RDB[ptr + SCORE_PTR_BUF];
  CheckPointer(FUNCTION_NAME, "(loc0)", BUF_ARRAY, loc0);

  /* Get pointer to data */

  loc0 = loc0 + idx*BUF_BLOCK_SIZE;
  CheckPointer(FUNCTION_NAME, "(loc0)", BUF_ARRAY, loc0);

  /* Check that buffer is not reduced */
      
  if ((long)RDB[DATA_BUF_REDUCED] == YES)
    Die(FUNCTION_NAME, "Trying to add to reduced buffer");

  /* Check if shared or private */

  if ((long)RDB[DATA_OPTI_SHARED_BUF] == YES)
    {
      /* Shared buffer, put data */

#ifdef OPEN_MP
#pragma omp atomic
#endif
      BUF[loc0 + BUF_VAL] += wgt*val;

#ifdef OPEN_MP
#pragma omp atomic
#endif
      BUF[loc0 + BUF_WGT] += wgt;

#ifdef OPEN_MP
#pragma omp atomic
#endif
      BUF[loc0 + BUF_N] += 1.0;
    }
  else
    {
      /* Get buffer segment size */

      sz = (long)RDB[DATA_REAL_BUF_SIZE];
      loc0 = loc0 + id*sz;

      /* Add data */
      
      BUF[loc0 + BUF_VAL] +=  wgt*val;
      BUF[loc0 + BUF_WGT] +=  wgt;
      BUF[loc0 + BUF_N] += 1.0;
    }

  /****************************************************************************/
}

/*****************************************************************************/
