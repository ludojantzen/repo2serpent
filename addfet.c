/******************************************************************************/
/*                                                                            */
/* serpent 2 (beta-version) : addfet.c                                        */
/*                                                                            */
/* Created:       2018/01/31 (BWe)                                            */
/* Last modified: 2018/02/23 (BWe)                                            */
/* Version:       2.1.31                                                      */
/*                                                                            */
/* Description: Adds results to FET score buffer                              */
/*                                                                            */
/* Comments: - This routine is distilled from AddBuf(...) from assuming:      */
/*                  - idx      == -1                                          */
/*                  - dim      == 2                                           */
/*                  - va_arg0  == idx                                          */
/*                  - va_arg1  == rbin                                        */
/*                  - Private buffering                                       */
/*                                                                            */
/******************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddFET:"

/*****************************************************************************/

void AddFET(const double *const paramsPtr, long statPtr, long particleID,
            double expansionVal, double scoreVal, double weight, long idx,
            long threadID, long rbin)
{
  long buffBin, particleIdx, size, uncertaintyPtr;
  const double weightedVal = expansionVal * scoreVal * weight;

  buffBin = (long)RDB[statPtr + SCORE_PTR_NMAX];
  idx = idx + rbin*(long)RDB[buffBin];

  /* Get pointer to buffer */

  buffBin = (long)RDB[statPtr + SCORE_PTR_BUF];
  CheckPointer(FUNCTION_NAME, "(buffBin)", BUF_ARRAY, buffBin);

  /* Get pointer to data */

  buffBin = buffBin + idx*BUF_BLOCK_SIZE;
  CheckPointer(FUNCTION_NAME, "(buffBin)", BUF_ARRAY, buffBin);

  /* Check that buffer is not reduced */

  if ((long)RDB[DATA_BUF_REDUCED] == YES)
    Die(FUNCTION_NAME, "Trying to add to reduced buffer");

  size = (long)RDB[DATA_REAL_BUF_SIZE];
  buffBin = buffBin + threadID*size;

  /* Add data */

  BUF[buffBin + BUF_FET_VAL] += weightedVal;

  uncertaintyPtr = (long)paramsPtr[FET_PARAM_UNC_ARRAY_PTR];
  particleIdx = (long)paramsPtr[FET_PARAM_NCOEF_TOTAL];
  if (particleID == GetPrivateData(uncertaintyPtr + particleIdx, threadID))
  {
    /* Particles are the same, so contribute the new addition */

    BUF[buffBin + BUF_FET_SUM_SQ_SUM] += weightedVal*(2*BUF[buffBin + BUF_FET_SUM] + weightedVal);
    BUF[buffBin + BUF_FET_SUM] += weightedVal;
  }
  else
  {
    /* Particles are different, so start the sum over */

    BUF[buffBin + BUF_FET_SUM_SQ_SUM] += weightedVal*weightedVal;
    BUF[buffBin + BUF_FET_SUM] = weightedVal;

    /* Is this the last coefficient? */

    if (idx == particleIdx - 1)
    {
      /* Increment the particle count */

      buffBin += BUF_BLOCK_SIZE; /* The next bin is just one block size away */
      BUF[buffBin + BUF_FET_VAL] += 1;

      /* Store the new particle ID */

      PutPrivateData(uncertaintyPtr + particleIdx, particleID, threadID);
    }
  }
}
