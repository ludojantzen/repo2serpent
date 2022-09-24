/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processburnupegroups.c                         */
/*                                                                           */
/* Created:       2011/11/08 (JLe)                                           */
/* Last modified: 2012/06/04 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: Sets up energy group structure used with spectrum-collapse   */
/*              method.                                                      */
/*                                                                           */
/* Comments: - Koko helvetin multi-group -hässäkkä poistettiin 4.6.2012      */
/*             (2.1.6), ja tämä rutiini pelkästään asettaa noi rajat         */
/*             vastaamaan ures-sämpläyksen rajoja.                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessBurnupEGroups:"

/*****************************************************************************/

void ProcessBurnupEGroups()
{
  /* Check burnup mode */

 if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
   return;

 /* Set limits */

 WDB[DATA_BU_URES_EMIN] = RDB[DATA_URES_EMIN];
 WDB[DATA_BU_URES_EMAX] = RDB[DATA_URES_EMAX];
}

/*****************************************************************************/
