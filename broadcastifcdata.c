/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : broadcastifcdata.c                             */
/*                                                                           */
/* Created:       2019/02/11 (VVa)                                           */
/* Last modified: 2019/02/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Broadcasts updated interfaces between MPI-tasks              */
/*                                                                           */
/* Comments: -Also transfers the iteration flag                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BroadcastIFCData:"

/*****************************************************************************/

void BroadcastIFCData()
{
#ifdef MPI
  long loc0, sz, ptr;
#endif

  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == NO)
    return;

#ifdef MPI

  /* Broadcast iteration flag from task 0 to other tasks */

  MPITransfer(&WDB[DATA_ITERATE], NULL, 1, 0, MPI_METH_BC);

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];

  while (loc0 > VALID_PTR)
    {
      sz = (long)RDB[loc0 + IFC_MEM_SIZE];

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Broadcast data from task 0 to other tasks */

      MPITransfer(&WDB[loc0], NULL, sz, 0, MPI_METH_BC);

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Next interface */

      loc0 = NextItem(loc0);

    }

  /* Loop over data interfaces */

  loc0 = (long)RDB[DATA_PTR_DATAIFC0];

  while (loc0 > VALID_PTR)
    {
      sz = (long)RDB[loc0 + DATAIFC_DATA_MEM_SIZE];

      /* Get pointer to actual data */

      ptr = (long)RDB[loc0 + DATAIFC_PTR_DATA];

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Broadcast data from task 0 to other tasks */

      MPITransfer(&WDB[ptr], NULL, sz, 0, MPI_METH_BC);

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Next interface */

      loc0 = NextItem(loc0);
    }
#endif
}
