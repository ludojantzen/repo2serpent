/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkfinishddsynch.c                           */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Checks for termination of the tracking loop in DD mode using */
/*              synchronous MPI communications.                              */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - Checks if all particles sent have been received.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckFinishDDSynch:"

/*****************************************************************************/

void CheckFinishDDSynch(long *done)
{

#ifdef MPI

  long dd_part_count_global;

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    {
      *done = 1;
      return;
    }

  /* Get the global particle balance */

  MPI_Allreduce(&dd_part_count, &dd_part_count_global, 1, MPI_LONG, MPI_SUM,
                my_comm);

  /* Check if all particles sent have been received */

  *done = dd_part_count_global == 0;

  /* Clean up after each transport cycle */

  if (*done)
    CleanUpDDComms();

  #endif
}

/*****************************************************************************/
