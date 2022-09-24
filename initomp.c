/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initomp.c                                      */
/*                                                                           */
/* Created:       2011/11/11 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Inits OpenMP related stuff                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitOMP:"

/*****************************************************************************/

void InitOMP()
{
  int opt, ptr;

  /***************************************************************************/

  /***** Init threads and set options ****************************************/

  /* Set number of threads */

#ifdef OPEN_MP

  if ((long)RDB[DATA_OMP_MAX_THREADS] > 0)
    omp_set_num_threads((int)RDB[DATA_OMP_MAX_THREADS]);
  else
    Die(FUNCTION_NAME, "Error in number of threads");

#else

  /* Check that thread number is set to 1 */

  if ((long)RDB[DATA_OMP_MAX_THREADS] != 1)
    Die(FUNCTION_NAME, "Multiple threads without OpenMP");

#endif

  /* Use shared arrays if no OpenMP parallelization */

  if ((long)RDB[DATA_OMP_MAX_THREADS] == 1)
    {
      WDB[DATA_OPTI_SHARED_BUF] = (double)YES;
      WDB[DATA_OPTI_SHARED_RES2] = (double)YES;
    }

  /* No MPI reproducibility of no OpenMP reporucibility */

  if (((long)RDB[DATA_OMP_MAX_THREADS] == 1) &&
      ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == NO))
    WDB[DATA_OPTI_MPI_REPRODUCIBILITY] = (double)NO;

  /***************************************************************************/

  /***** Allocate memory from private array **********************************/

  /* Allocate dummy values to disable zeros from pointer check */

  AllocPrivateData(DATA_FIXED_BLOCK_SIZE, PRIVA_ARRAY);
  AllocPrivateData(DATA_FIXED_BLOCK_SIZE, BUF_ARRAY);

  /* RES2-arraylle tehdään tällainen kamala viritys, jotta myöhemmin   */
  /* readinput.c:ssä vaihdettava optio ei sotke muistinvarausta.       */
  /* BUF-arraylle vastaavaa ei tarvitse tehdä, niin kauan kuin         */
  /* initdata.c:ssä asetetaan defaultiksi NO. (JLe/22.5.2018/2.1.31)   */

  opt = (long)RDB[DATA_OPTI_SHARED_RES2];
  WDB[DATA_OPTI_SHARED_RES2] = (double)NO;

  AllocPrivateData(DATA_FIXED_BLOCK_SIZE, RES2_ARRAY);
  WDB[DATA_OPTI_SHARED_RES2] = (double)opt;

  /* Allocate memory for variables in the private data block */

  WDB[DATA_PTR_COLLISION_COUNT] = (double)AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_DBRC_COUNT] = (double)AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_DBRC_EXCEED_COUNT] = (double)AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_OMP_HISTORY_COUNT] = (double)AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_PARTICLE_INITIAL_TIME] =
    (double)AllocPrivateData(1, PRIVA_ARRAY);

  /* Allocate memory for completed flags */

  ptr = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PTR_OMP_COMPLETED] = (double)ptr;

  /***************************************************************************/
}

/*****************************************************************************/
