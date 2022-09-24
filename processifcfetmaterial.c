/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcfetmaterial.c                        */
/*                                                                           */
/* Created:       2018/02/23 (BWe)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes FET multi-physics interface materials              */
/*                                                                           */
/* Comments: - Copied and modified from `processifcfunc.c`                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCFETMaterial:"

/*****************************************************************************/

void ProcessIFCFETMaterial(long loc0, long update)
{

  if (update == (long)NO)
    {
      /* Link materials */

      LinkInterfaceMaterials(loc0);
    }

  /* Set or check TMS limits for linked materials */

  SetIFCTMSLimits(loc0, update);

  return;
}
