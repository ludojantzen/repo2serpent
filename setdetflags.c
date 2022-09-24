/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setdetflags.c                                  */
/*                                                                           */
/* Created:       2016/08/01 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Sets / resets detector flags                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetDetFlags:"

/*****************************************************************************/

void SetDetFlags(long det, long part)
{
  long ptr, tst;
  
  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "part", DATA_ARRAY, part);

  /* Loop over flags */
  
  ptr = (long)RDB[det + DET_PTR_FLAGGING];
  while (ptr > VALID_PTR)
    {
      /* Get flag number */
      
      tst = (long)RDB[ptr + DET_FBIN_FLAG_NUMBER];
      CheckValue(FUNCTION_NAME, "tst", "", tst, 1, 64);
      
      /* Convert */
      
      tst = (long)(pow(2.0, (double)tst - 1.0));
      CheckValue(FUNCTION_NAME, "tst", "", tst, 1, LONG_MAX);
      
      /* Check set and reset options */
      
      if ((long)RDB[ptr + DET_FBIN_FLAG_OPTION] == DET_FLAG_OPT_SET)
        SetOption(part + PARTICLE_DET_FLAGS, tst);
      else if ((long)RDB[ptr + DET_FBIN_FLAG_OPTION] == DET_FLAG_OPT_RESET)
        ResetOption(part + PARTICLE_DET_FLAGS, tst);
      
      /* Next flag */
      
      ptr = NextItem(ptr);
    }
}

/*****************************************************************************/
