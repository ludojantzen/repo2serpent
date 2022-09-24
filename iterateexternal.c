/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iterateexternal.c                              */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2018/06/01 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calls external solver by POSIX-signalling and reads new      */
/*              interfaces                                                   */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterateExternal:"

/*****************************************************************************/

void IterateExternal()
{
  long ptr;

  /* If no external coupling */

  if (RDB[DATA_CC_SIG_MODE] == (double)SIG_MODE_NONE)
    return;

  PrintInterfaceOutput();

  /* Also write out the detector data */

  DetectorOutput();

  /* Send continue this depletion step signal */
  /* and sleep until Serpent is needed again  */

  SignalExternal(SIGUSR1);

  /* Read updated interfaces */

  ptr = (long)RDB[DATA_PTR_IFC0];

  while(ptr > VALID_PTR)
    {
      /* Update interface */

      ReadInterface(ptr, YES);

      /* Next interface */

      ptr = NextItem(ptr);

    }

  /* Process updated interfaces */

  ProcessInterface(YES);

  return;
}
