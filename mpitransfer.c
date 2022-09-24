/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : mpitransfer.c                                  */
/*                                                                           */
/* Created:       2010/11/23 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Broadcasts data block to MPI tasks                           */
/*                                                                           */
/* Comments: - From Serpent 1.1.12                                           */
/*           - Se ongelma mikä korjattiin tolla batch-jutulla saattoi johtua */
/*             siitä että luvut oli single-precision. Rutiini taisi tosin    */
/*             nopeutuakin tolla.                                            */
/*           - Serpent 2:ssa kokeillaan ilman batchaysta ja vaihdetaan       */
/*             siihen jos tulee ongelmia                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MPITransfer:"

/*****************************************************************************/

void MPITransfer(double *dat, double *buf, long sz, long root, long meth)
{

#ifdef MPI

  long tot, sz0, rc;

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Avoid compiler warning */

  rc = 0;
  sz0 = 0;

  /* Get MPI batch size */

  if (RDB == NULL)
    Die(FUNCTION_NAME, "Main data block not allocated (task %d)", mpiid);
  else if ((sz0 = (long)RDB[DATA_OPTI_MPI_BATCH_SIZE]) < 1)
    Die(FUNCTION_NAME, "Batch size not defined (task %d)", mpiid);

  /* Check data size */

  if (sz < sz0)
    {
      /***********************************************************************/

      /***** Transfer data as a single block *********************************/

      /* Broadcast or reduce block */

      if (meth == MPI_METH_BC)
        rc = MPI_Bcast(dat, sz, MPI_DOUBLE, root, my_comm);
      else if (meth == MPI_METH_RED)
        rc = MPI_Reduce(dat, buf, sz, MPI_DOUBLE, MPI_SUM, root,
                        my_comm);
      else if (meth == MPI_METH_ALL_RED)
        rc = MPI_Allreduce(dat, buf, sz, MPI_DOUBLE, MPI_SUM, my_comm);
      else
        Die(FUNCTION_NAME, "Invalid function type %ld", meth);

      /* Synchronise */

      MPI_Barrier(my_comm);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Divide transfer into batches ************************************/

      /* Reset total transfered size */

      tot = 0;

      /* Loop over batches */

      do
        {
          /* Compare to total */

          if (tot + sz0 > sz)
            sz0 = sz - tot;

          /* Broadcast or reduce batch */

          if (meth == MPI_METH_BC)
            rc = MPI_Bcast(&dat[tot], sz0, MPI_DOUBLE, root, my_comm);
          else if (meth == MPI_METH_RED)
            rc = MPI_Reduce(&dat[tot], &buf[tot], sz0, MPI_DOUBLE, MPI_SUM,
                            root, my_comm);
          else
            Die(FUNCTION_NAME, "Invalid function type %ld", meth);

          /* Break on error */

          if (rc != MPI_SUCCESS)
            break;

          /* Add to total */

          tot = tot + sz0;

          /* Synchronise */

          MPI_Barrier(my_comm);
        }
      while (tot < sz);

      /***********************************************************************/
    }

  /* Check error */

  if (rc != MPI_SUCCESS)
    Die(FUNCTION_NAME, "Data transfer failed with error condition %ld", rc);

#endif
}

/*****************************************************************************/
