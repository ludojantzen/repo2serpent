/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocstathistory.c                             */
/*                                                                           */
/* Created:       2011/05/18 (JLe)                                           */
/* Last modified: 2012/10/25 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Allocates memory for cycle-wise statistics                   */
/*                                                                           */
/* Comments: - Nää voisi varastoida ihan omaan taulukkoon                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocStatHistory:"

/*****************************************************************************/

void AllocStatHistory(long loc0)
{
  long ntot, ncyc, ptr;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get stat size */
  
  ntot = (long)RDB[loc0 + SCORE_STAT_SIZE];
  
  /* Get number of criticality cycles */
  
  ncyc = (long)RDB[DATA_CRIT_CYCLES] + (long)RDB[DATA_CRIT_SKIP];

  /* Compare to number of source batches */

  if (ncyc < (long)RDB[DATA_SRC_BATCHES])
    ncyc = (long)RDB[DATA_SRC_BATCHES];

  /* Allocate memory */

  ptr = ReallocMem(RES1_ARRAY, ntot*ncyc);
  
  /* Put pointer */

  WDB[loc0 + SCORE_PTR_HIS] = (double)ptr;  
}

/*****************************************************************************/
