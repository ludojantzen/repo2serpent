/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reduceprivateres.c                             */
/*                                                                           */
/* Created:       2010/11/12 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Reduces data from OpenMP private results array to thread 0   */
/*                                                                           */
/* Comments: - NOTE: Tähän liittyy nyt sellanen potentiaalinen ongelma että  */
/*             tonne ykkössegmenttiin voi kerääntyä isoja arvoja, ja sit     */
/*             numeriikka alkaa aiheuttamaan häikkää kun niihin aletaan      */
/*             lisäämään pieniä arvoja muista.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReducePrivateRes:"

/*****************************************************************************/

void ReducePrivateRes()
{
  long sz, i, n, max;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");

  /* Check and set reduced flag */

  if ((long)RDB[DATA_RES2_REDUCED] == YES)
    return;
  else
    WDB[DATA_RES2_REDUCED] = (double)YES;

  /* Return if shared */

  if ((long)RDB[DATA_OPTI_SHARED_RES2] == YES)
    return;

  /* Get buffer segment and data size */

  sz = (long)RDB[DATA_REAL_RES2_SIZE];
  max = (long)RDB[DATA_ALLOC_RES2_SIZE];

  /* Loop over OpenMP threads, copy and reset data */

  for (i = 1; i < (long)RDB[DATA_OMP_MAX_THREADS]; i++)
    for (n = 0; n < max; n++)
      {
        RES2[n] = RES2[n] + RES2[i*sz + n];
        RES2[i*sz + n] = 0.0;
      }

  /****************************************************************************/
}

/*****************************************************************************/
