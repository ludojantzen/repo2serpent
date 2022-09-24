/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resizefissionsrc.c                             */
/*                                                                           */
/* Created:       2014/11/01 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Resizes fission source population to DATA_CRIT_POP           */
/*              if fission source passing is used                            */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResizeFissionSrc:"

/*****************************************************************************/

void ResizeFissionSrc()
{
  long ptr, part, id, n, N, m;
  long nbatch, nsrc, np, nb, mul;
  double wgt0, wgt1, wgt2, P;

  /***************************************************************************/

  /* Get wanted criticality population */

#ifdef MPI_MODE1

  /* Check number of tasks */

  if (mpitasks > 1)
    {
      /* Calculate number of particles per task */

      nbatch = (long)(RDB[DATA_CRIT_POP]/((double)mpitasks));

    }
  else
    nbatch = (long)RDB[DATA_CRIT_POP];

#else

  nbatch = (long)RDB[DATA_CRIT_POP];

#endif

  /* Reset current population */

  nsrc = 0;

  /* Reset current weight */

  wgt0 = 0;

  /* Calculate current bank counts */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while((ptr = NextItem(ptr)) > VALID_PTR)
        {

          /* Add to count */

          nsrc++;

          /* Add to weight */

          wgt0 = wgt0 + RDB[ptr + PARTICLE_WGT];
        }
    }

  /* Reset new population */

  np = nsrc;

  /* Reset new weight */

  wgt1 = wgt0;

  /* Compare population size to batch size given in input */

  if (nsrc < nbatch)
    {
      /* Calculate multiplication */

      P = ((double)nbatch)/((double)nsrc);
      mul = (long)P;
      P = P - (double)mul;

      /* Loop over bank */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {

          /* Get pointer to bank */

          ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get number of neutrons in bank */

          nb = ListSize(ptr) - 1;

          /* Get pointer to dummy */

          ptr = FirstItem(ptr);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get pointer to first item after dummy */

          if((ptr = NextItem(ptr)) < VALID_PTR)
            continue;

          /* Loop over old bank and sample multiplication */

          for (n = 0; n < nb; n++)
            {
              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Sample multiplication */

              if (RandF(0) < P)
                N = mul;
              else
                N = mul - 1;

              /* Loop over multiplication */

              for (m = 0; m < N; m++)
                {
                  /* Duplicate neutron */

                  part = DuplicateParticle(ptr, id);

                  /* Add to new weight and population size */

                  wgt1 = wgt1 + RDB[ptr + PARTICLE_WGT];
                  np = np + 1;

                  /* Put particle in source*/

                  ToBank(part, id);
                }

              /* Next particle */

              ptr = NextItem(ptr);

            }

        }

    }
  else if (nsrc > nbatch)
    {
      /* Calculate probability to remove */

      P = 1.0 - ((double)nbatch)/((double)nsrc);

      /* Loop over bank */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {

          /* Get pointer to bank */

          ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get number of neutrons in bank */

          nb = ListSize(ptr) - 1;

          /* Get pointer to dummy */

          ptr = FirstItem(ptr);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get pointer to first item after dummy */

          if((ptr = NextItem(ptr)) < VALID_PTR)
            continue;

          /* Loop over old bank and sample removal */

          for (n = 0; n < nb; n++)
            {

              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get pointer and next */

              part = ptr;
              ptr = NextItem(ptr);

              /* Sample removal */

              if (RandF(0) < P)
                {
                  /* Remove particle */

                  RemoveItem(part);

                  /* subtract from weight and population size */

                  wgt1 = wgt1 - RDB[part + PARTICLE_WGT];
                  np = np - 1;

                  /* Subtract from current bank size so that */
                  /* for loop does not overflow              */

                  nb = nb - 1;

                  /* Put particle back to stack */

                  ToStack(part, id);

                }
            }
        }

    }

  /* Reset final weight */

  wgt2 = 0.0;

  /* Re-normalize total weight to conserve weight */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while((ptr = NextItem(ptr)) > VALID_PTR)
        {

          /* Re-normalize weight */

          WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*wgt0/wgt1;

          /* Add to final weight */

          wgt2 = wgt2 + RDB[ptr + PARTICLE_WGT];
        }
    }

  /* Check final weight */

  if (wgt0 > 0)
    if(fabs(wgt0 - wgt2)/wgt0 > 1e-5)
      Die(FUNCTION_NAME,"Mismatch in weight %E %E (%E)", wgt0, wgt2, wgt0 - wgt2);

  /* Set batch size */

  WDB[DATA_SIMUL_BATCH_SIZE] = (double)np;
}

/*****************************************************************************/
