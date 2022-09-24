/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findsenszaiindex.c                             */
/*                                                                           */
/* Created:       2018/06/19 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds index for ZAI from sensitivity calculation ZAI list    */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindSensZAIIndex:"

/*****************************************************************************/

long FindSensZAIIndex(long ZAI)
{
  long loc0, ptr, idx, izai, myzai;

  /* Get pointer to sensitivity block or return */

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return -1;

  /* Get pointer to zai list or return -1 for no match */

  if ((ptr = (long)RDB[loc0 + SENS_PTR_ZAI_ARR]) < VALID_PTR)
    return -1;

  /* Loop over ZAI list to find ZAI */

  idx = 0;
  while ((myzai = (long)RDB[ptr + idx]) >= 10010)
    {

      if (myzai == ZAI)
        break;

      idx++;
    }

  /* Check if found */

  if (myzai < 10010)
    return -1;

  /* Was found */

  /* Increment index due to tot and sum if needed */

  izai = idx;

  /* Increment index to take in account none, total and sum bins */

  izai = izai + 1;

  if ((long)RDB[loc0 + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_TOT)
    izai = izai + 1;

  if ((long)RDB[loc0 + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_SUM)
    izai = izai + 1;

  return izai;
}
