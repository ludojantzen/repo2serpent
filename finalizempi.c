/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : main.c                                         */
/*                                                                           */
/* Created:       2011/11/29 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finalizes MPI before exiting program                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "FinalizeMPI:"

/*****************************************************************************/

void FinalizeMPI()
{
#ifdef MPI

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Finalize */

  if (MPI_Finalize() != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

#endif
}

/*****************************************************************************/
