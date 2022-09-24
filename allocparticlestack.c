/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocparticlestacs.c                           */
/*                                                                           */
/* Created:       2012/10/17 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Allocates memory for particle histories (stacks)             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocParticleStack:"

/*****************************************************************************/

void AllocParticleStack(long type, long np)
{
  long ptr, loc0, id, n, sz;

#ifdef OLD_IFP

  long m;

#endif

  /* Calculate size for preallocation (per particle) */

  sz = PARTICLE_BLOCK_SIZE;

  /* Avoid compiler warning */

  loc0 = -1;

  /* Check type */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Get pointer */

      loc0 = DATA_PART_PTR_NSTACK;

      /* Update stack size */

      WDB[DATA_PART_ALLOC_N] = RDB[DATA_PART_ALLOC_N] + (double)np;

      /* Minimum stack size */

      if ((long)RDB[DATA_PART_PTR_MIN_NSTACK] < VALID_PTR)
        {
          ptr = AllocPrivateData(1, PRIVA_ARRAY);
          WDB[DATA_PART_PTR_MIN_NSTACK] = (double)ptr;
        }

#ifdef OLD_IFP

      /* Add IFP fission progeny size */

      sz = sz + (long)RDB[DATA_IFP_CHAIN_LENGTH]*FISS_PROG_BLOCK_SIZE;

      /* Add fission progeny list common size */

      sz = sz + LIST_COMMON_DATA_SIZE;

#endif

    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      /* Get pointer */

      loc0 = DATA_PART_PTR_GSTACK;

      /* Update stack size */

      WDB[DATA_PART_ALLOC_G] = RDB[DATA_PART_ALLOC_G] + (double)np;

      /* Minimum stack size */

      if ((long)RDB[DATA_PART_PTR_MIN_GSTACK] < VALID_PTR)
        {
          ptr = AllocPrivateData(1, PRIVA_ARRAY);
          WDB[DATA_PART_PTR_MIN_GSTACK] = (double)ptr;
        }
    }
  else if (type == PARTICLE_TYPE_PRECURSOR)
    {
      /* Get pointer */

      loc0 = DATA_PART_PTR_PSTACK;

      /* Update stack size */

      WDB[DATA_PART_ALLOC_P] = RDB[DATA_PART_ALLOC_P] + (double)np;

      /* Minimum stack size */

      if ((long)RDB[DATA_PART_PTR_MIN_PSTACK] < VALID_PTR)
        {
          ptr = AllocPrivateData(1, PRIVA_ARRAY);
          WDB[DATA_PART_PTR_MIN_PSTACK] = (double)ptr;
        }
    }
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Preallocate memory (speeds things up by a lot) */

  PreallocMem(sz*np, DATA_ARRAY);

  /* Reset OpenMP index */

  id = 0;

  /* Loop over particles */

  for (n = 0; n < np; n++)
    {
      /* Allocate memory */

      ptr = NewItem(OMPPtr(loc0, id), PARTICLE_BLOCK_SIZE);

      /* Put type */

      WDB[ptr + PARTICLE_TYPE] = (double)type;

      /* Allocate memory for fission progenies */

#ifdef OLD_IFP

      if (type == PARTICLE_TYPE_NEUTRON)
        for (m = 0; m < (long)RDB[DATA_IFP_CHAIN_LENGTH]; m++)
          NewItem(ptr + PARTICLE_PTR_FISS_PROG, FISS_PROG_BLOCK_SIZE);
#endif

      /* Update OpenMP id */

      if (++id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
        id = 0;
    }

  /* Update memory size */

  WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + (double)MemCount();
}

/*****************************************************************************/
