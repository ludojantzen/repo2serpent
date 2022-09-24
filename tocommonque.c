/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tocommonque.c                                  */
/*                                                                           */
/* Created:       2019/02/21 (JLe)                                           */
/* Last modified: 2019/02/21 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Puts neutron / photon in common que                          */
/*                                                                           */
/* Comments: Used for load balancing, e.g. production of secondary photons.  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToCommonQue:"

/*****************************************************************************/

void ToCommonQue(long ptr)
{
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Put particle in que */

#ifdef OPEN_MP
#pragma omp critical (que)
#endif
  {
    AddItem(DATA_PART_PTR_COMMON_QUE, ptr);
  }
}

/*****************************************************************************/
