/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pause.c                                        */
/*                                                                           */
/* Created:       2020/03/06 (JLe)                                           */
/* Last modified: 2020/03/06 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Waits for any key to continue.                               */
/*                                                                           */
/* Comments: - For debugging.                                                */
/*                                                                           */
/*           - Use cautiously in parallel loops.                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Pause:"

/*****************************************************************************/

void Pause(char *routine)
{
  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Print message */

  fprintf(outp, "\nPaused in %s Press enter key to continue...\n", routine);

  /* Wait for key */

  getchar();
}

/*****************************************************************************/
