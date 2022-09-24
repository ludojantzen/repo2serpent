/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addddres.c                                     */
/*                                                                           */
/* Created:       2018/03/14 (JLe)                                           */
/* Last modified: 2018/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Adds value in RES3 data block                                */
/*                                                                           */
/* Comments: - Used only with domain decomposition to store transmutation    */
/*             cross sections.                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddDDRes:"

/*****************************************************************************/

void AddDDRes(long ptr, double val)
{
  /* Check pointer and value */

  CheckPointer(FUNCTION_NAME, "(ptr)", RES3_ARRAY, ptr);
  CheckValue(FUNCTION_NAME, "val", "", val, 0.0, INFTY);

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES3 array not ready for access");

  /* Check domain decomposition */
  
  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    Die(FUNCTION_NAME, "Access to RES3 withou domain decomposition");

  /* Put data */

#ifdef OPEN_MP
#pragma omp atomic
#endif
  
  RES3[ptr] += val;
}

/*****************************************************************************/
