/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : flushbank.c                                    */
/*                                                                           */
/* Created:       2012/10/10 (JLe)                                           */
/* Last modified: 2017/12/09 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Removes all particles from bank                              */
/*                                                                           */
/* Comments: - Stackien kasvattaminen ei oo ehkä ihan ideaalisinta tehdä     */
/*             täällä.                                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FlushBank:"

/*****************************************************************************/

void FlushBank()
{
  long ptr, part, id, minn, ming, n, m;

  /* Get pointer to source */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Check that source is empty */

  if (ListSize(ptr) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Loop until bank is empty */

      while ((part = FromBank(id)) > VALID_PTR)
        {
          /* Put particle back to stack */

          ToStack(part, id);
        }
    }

  /* Reset minimum sizes */

  minn = 10000000000;
  ming = 10000000000;

  /* Loop over OpenMP threads */

  for (n = 0; n < (long)RDB[DATA_OMP_MAX_THREADS]; n++)
    {
      /*
      printf("%2ld :", n);
      */
      /* Neutron stack */

      if ((ptr = (long)RDB[DATA_PART_PTR_MIN_NSTACK]) > VALID_PTR)
        {
          /* Get value */

          m = (long)GetPrivateData(ptr, n);

          /* Compare to minimum */
          
          if ((m > 0) && (m < minn))
            minn = (long)GetPrivateData(ptr, n);
          /*
          printf(" N: %ld / %ld", (long)GetPrivateData(ptr, n),
                 (long)(RDB[DATA_PART_ALLOC_N]/RDB[DATA_OMP_MAX_THREADS]));
          */
          /* Reset */

          PutPrivateData(ptr, 0, n);
        }

      /* Photon stack */

      if ((ptr = (long)RDB[DATA_PART_PTR_MIN_GSTACK]) > VALID_PTR)
        {
          /* Get value */

          m = (long)GetPrivateData(ptr, n);

          /* Compare to minimum */
          
          if ((m > 0) && (m < ming))
            ming = (long)GetPrivateData(ptr, n);
          /*
          printf(" G: %ld / %ld", (long)GetPrivateData(ptr, n),
                 (long)(RDB[DATA_PART_ALLOC_G]/RDB[DATA_OMP_MAX_THREADS]));
          */
          /* Reset */

          PutPrivateData(ptr, 0, n);
        }
      /*
      printf("\n");
      */
    }

  /* Check minimum neutron stack size */

  if (RDB[DATA_OMP_MAX_THREADS]*minn/RDB[DATA_PART_ALLOC_N] < 0.2)
    {
      /* Allow memory allocation */

      Mem(MEM_ALLOW);

      /* Allocate 20% more particles */

      AllocParticleStack(PARTICLE_TYPE_NEUTRON, 
                         (long)(0.2*RDB[DATA_PART_ALLOC_N]));

      /* Disallow memory allocation */

      Mem(MEM_DENY);
    }

  /* Check minimum photon stack size */

  if (RDB[DATA_OMP_MAX_THREADS]*ming/RDB[DATA_PART_ALLOC_G] < 0.2)
    {
      /* Allow memory allocation */
      
      Mem(MEM_ALLOW);

      /* Allocate 20% more particles */

      AllocParticleStack(PARTICLE_TYPE_GAMMA, 
                         (long)(0.2*RDB[DATA_PART_ALLOC_G]));

      /* Disallow memory allocation */

      Mem(MEM_DENY);
    }

  /* Re-distribute stacks (does now work in track plot mode) */
  
  if ((long)RDB[DATA_STOP_AFTER_PLOT] != STOP_AFTER_PLOT_TRACKS)
    ReDistributeStacks();
}

/*****************************************************************************/
