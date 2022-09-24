/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addbuf1d.c                                     */
/*                                                                           */
/* Created:       2014/04/04 (JLe)                                           */
/* Last modified: 2014/04/05 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Simplified version of addbuf.c for 1d results                */
/*                                                                           */
/* Comments: - Tarkoitus optimoida tai ainakin helpottaa profilointia.       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddBuf1D:"

/*****************************************************************************/

void AddBuf1D(double val, double wgt, long ptr, long id, long idx)
{
  long loc0, sz;

  /* Check pointer and id */

  CheckPointer(FUNCTION_NAME, "", DATA_ARRAY, ptr);
  CheckValue(FUNCTION_NAME, "id", "", id, 0, MAX_OMP_THREADS);

#ifdef DEBUG

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "BUF array not ready for access");

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
}

/*****************************************************************************/
