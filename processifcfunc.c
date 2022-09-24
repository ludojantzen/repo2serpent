/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcfunc.c                               */
/*                                                                           */
/* Created:       2015/02/02 (VVa)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes functional multi-physics interfaces                */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCFunc:"

/*****************************************************************************/

void ProcessIFCFunc(long loc0, long update)
{

  /***********************************************************************/

  if (update == (long)NO)
    {
      /* Link materials */

      LinkInterfaceMaterials(loc0);
    }

  /* Set or check TMS limits for linked materials */

  SetIFCTMSLimits(loc0, update);

  return;

}
