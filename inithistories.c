/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : inithistories.c                                */
/*                                                                           */
/* Created:       2011/04/01 (JLe)                                           */
/* Last modified: 2019/10/16 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Initializes particle stacks, ques, source and bank for       */
/*              transport simulation                                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitHistories:"

/*****************************************************************************/

void InitHistories()
{
  long ptr, loc0, sens, id, np, n;

  fprintf(outp, "Allocating memory for particle structures...\n");

  /***************************************************************************/

  /***** Create lists ********************************************************/

  /* Neutron stacks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_NSTACK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Gamma stacks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_GSTACK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Precursor stacks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_PSTACK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Particle ques */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_QUE] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Banks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_BANK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Limbos (for domain decomposition) */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
      WDB[DATA_PART_PTR_LIMBO] = (double)loc0;

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
          WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
          WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
        }
    }

  /* Separate banks for track plotting */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_TRK_BANK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Source (no division to threads) */

  ptr = NewItem(DATA_PART_PTR_SOURCE, PARTICLE_BLOCK_SIZE);
  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;

  /* Precursor source (no division to threads) */

  ptr = NewItem(DATA_PART_PTR_PSOURCE, PARTICLE_BLOCK_SIZE);
  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;

  /* Set fraction of weight to emit on next interval to 2.0 */
  /* Precursor source will be sorted based on this in a     */
  /* descending order (precursors should have this value    */
  /* in interval (0,1) */

  WDB[ptr + PARTICLE_U] = 2.0;

  /* Particle pool (no division to threads) */

  ptr = NewItem(DATA_PART_PTR_POOL, PARTICLE_BLOCK_SIZE);
  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;

  /* Common que (no division to threads) */

  ptr = NewItem(DATA_PART_PTR_COMMON_QUE, PARTICLE_BLOCK_SIZE);
  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;

  /***************************************************************************/

  /***** Allocate memory for neutrons and photons ****************************/

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
    {
      /* Number of particles */

      if (RDB[DATA_CRIT_POP] > RDB[DATA_SRC_POP])
        np = (long)(RDB[DATA_PART_NBUF_FACTOR]*RDB[DATA_CRIT_POP]);
      else if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN) ||
               ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN))
        np = (long)(RDB[DATA_PART_NBUF_FACTOR]*RDB[DATA_SRC_POP]);
      else
        /*
        np = (long)(RDB[DATA_PART_NBUF_FACTOR]*RDB[DATA_SRC_POP]*
                    (1.0 + 0.25*(RDB[DATA_OMP_MAX_THREADS] - 1.0)));
        */
        np = (long)(RDB[DATA_PART_NBUF_FACTOR]*RDB[DATA_OMP_MAX_THREADS]);

#ifdef MPI_MODE1

      /* Check number of tasks */

      if (mpitasks > 1)
        np = (long)((double)np/(double)mpitasks);

#endif

      /* Allocate memory for neutrons */

      AllocParticleStack(PARTICLE_TYPE_NEUTRON, np);
    }

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Number of particles */

      np = (long)(RDB[DATA_PART_GBUF_FACTOR]*RDB[DATA_OMP_MAX_THREADS]);

#ifdef MPI_MODE1

      /* Check number of tasks */

      if (mpitasks > 1)
        np = (long)((double)np/(double)mpitasks);

#endif

      /* Allocate memory for photons */

      AllocParticleStack(PARTICLE_TYPE_GAMMA, np);
    }


  if ((long)RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    {
      /* Number of particles */

      if (RDB[DATA_CRIT_POP] > RDB[DATA_SRC_POP])
        np = (long)(RDB[DATA_PART_PBUF_FACTOR]*RDB[DATA_CRIT_POP]);
      else
        /*
        np = (long)(RDB[DATA_PART_PBUF_FACTOR]*RDB[DATA_SRC_POP]*
                    (1.0 + 0.25*(RDB[DATA_OMP_MAX_THREADS] - 1.0)));

        np = (long)(RDB[DATA_PART_PBUF_FACTOR]*RDB[DATA_OMP_MAX_THREADS]);
        */
        np = (long)(RDB[DATA_PART_PBUF_FACTOR]*RDB[DATA_SRC_POP]);


#ifdef MPI_MODE1

      /* Check number of tasks */

      if (mpitasks > 1)
        np = (long)((double)np/(double)mpitasks);

#endif

      /* Allocate memory for precursors */

      AllocParticleStack(PARTICLE_TYPE_PRECURSOR, np);
    }

  /***************************************************************************/

  /***** Allocate memory for events ******************************************/

  /* Record if any flags are set */

  if ((long)RDB[DATA_EVENT_RECORD_FLAGS] > 0)
    {
      /* Multiply number of events by (single thread) population size */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        np = (long)(RDB[DATA_EVENT_BANK_SZ]*RDB[DATA_CRIT_POP]/
                    RDB[DATA_OMP_MAX_THREADS]);
      else
        np = (long)(RDB[DATA_EVENT_BANK_SZ]*RDB[DATA_SRC_POP]/
                    RDB[DATA_OMP_MAX_THREADS]);

      /* Allocate memory for event banks */

      loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
      WDB[DATA_PTR_EVENT_BANK] = (double)loc0;

      /* Allocate memory for events */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          for (n = 0; n < np; n++)
            NewLIFOItem(loc0, EVENT_BLOCK_SIZE);

          loc0++;
        }

      /* Put new value */

      WDB[DATA_EVENT_BANK_SZ] = (double)np*RDB[DATA_OMP_MAX_THREADS];
    }

  /* Check if Sensitivity events are recorded using event blocks */

  if (((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_SENS)
      && ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_DIRECT))
    {
      /* Number of event blocks for each particle */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        np = (long)(RDB[DATA_EVENT_MAX_GEN]*RDB[DATA_EBLOCK_BANK_SZ]
                    *RDB[DATA_CRIT_POP]/RDB[DATA_OMP_MAX_THREADS]);
      else
        np = (long)(RDB[DATA_EVENT_MAX_GEN]*RDB[DATA_EBLOCK_BANK_SZ]
                    *RDB[DATA_SRC_POP]/RDB[DATA_OMP_MAX_THREADS]);

      /* Allocate memory for event block banks */

      loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
      WDB[DATA_PTR_EBLOCK_BANK] = (double)loc0;

      /* Allocate memory for events */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          for (n = 0; n < np; n++)
            NewLIFOItem(loc0, SENS_EBLOCK_BLOCK_SIZE);

          loc0++;
        }

      /* Put total bank size */

      WDB[DATA_EBLOCK_BANK_SZ] = (double)np*RDB[DATA_OMP_MAX_THREADS];

      /* Get pointer to sensitivity block */

      sens = (long)RDB[DATA_PTR_SENS0];
      CheckPointer(FUNCTION_NAME, "(sens)", DATA_ARRAY, sens);

      /* Get maximum event label */

      np = (long)RDB[sens + SENS_MAX_LABEL];

      /* Preallocate memory for data */

      PreallocMem((long)RDB[DATA_EBLOCK_BANK_SZ]*np, DATA_ARRAY);

      /* Loop over threads */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          /* Get pointer to first item in bank */

          loc0 = (long)RDB[OMPPtr(DATA_PTR_EBLOCK_BANK, id)];

          /* Loop over blocks to allocate memory */

          while (loc0 > VALID_PTR)
            {
              /* Allocate memory for scores */

              ptr = ReallocMem(DATA_ARRAY, np);

              /* Store pointer */

              WDB[loc0 + SENS_EBLOCK_PTR_DATA] = (double)ptr;

              /* Next item */

              loc0 = NextItem(loc0);
            }
        }

    }

  /* History index for debugging */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_PRIVA_HIS_IDX] = (double)ptr;

  /* Number of independent eigenvalue simulations */

  np = (long)RDB[DATA_N_POP_EIG];
  CheckValue(FUNCTION_NAME, "np", "", np, 1, 10000);

  /* Allocate memory for cycle-keffs */

  ptr = ReallocMem(DATA_ARRAY, np);
  WDB[DATA_PTR_CYCLE_EIG_KEFF] = (double)ptr;

  /* Init */

  for (n = 0; n < np; n++)
    WDB[ptr + n] = RDB[DATA_CYCLE_KEFF];

  /***************************************************************************/

  /* Everything is OK */

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
