/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizedynsrc.c                              */
/*                                                                           */
/* Created:       2012/09/23 (JLe)                                           */
/* Last modified: 2019/03/16 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Normalizes source distribution for dynamic criticality       */
/*              source simulation                                            */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizeDynSrc:"

/*****************************************************************************/

long NormalizeDynSrc()
{
  long ptr, pos, part, n, m, nsrc, nbatch, mul, N, np, id, stp, loc0, min;
  long type, rr;
  double wgt0, wgt, P, wmin;
  unsigned long seed;

  /***************************************************************************/

  /***** Move banked neutrons to source **************************************/

  /* Flag for russian roulette */

  if ((wmin = RDB[DATA_DYN_WMIN]) > 0.0)
    rr = YES;
  else
    rr = NO;

  /* Get pointer to source */

  pos = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);

  /* Check that source is empty */

  if (ListSize(pos) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Reset total weight and source size */

  wgt0 = 0.0;
  nsrc = 0;

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      while ((ptr = FromBank(id)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Check actual cut-off in interval mode */

          if ((RDB[ptr + PARTICLE_T] > RDB[DATA_TIME_CUT_TMAX]) &&
              (RDB[DATA_DYN_DT] > 0.0))
            {
              /* Score population and weight balance */

              stp = (long)RDB[RES_N_BALA_LOSS];
              CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
              AddBuf(1.0, 1.0, stp, 0, -1, BALA_N_LOSS_CUT, 0);
              AddBuf(RDB[ptr + PARTICLE_WGT], 1.0, stp, 0, -1, 
                     BALA_N_LOSS_CUT, 1);

              /* Put particle back to stack */

              ToStack(ptr, id);
            }
          else
            {
              /* Add to total weight and source size */
              
              wgt0 = wgt0 + RDB[ptr + PARTICLE_WGT];
              nsrc = nsrc + 1;
              
              /* Put neutron in source */
              
              if ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == YES)
                AddSortItem(DATA_PART_PTR_SOURCE, pos, ptr, PARTICLE_RNG_IDX,
                            SORT_MODE_ASCEND);
              else
                AddItem(DATA_PART_PTR_SOURCE, ptr);
              
              /* Update position */
              
              pos = ptr;
            }
        }
    }

  /* Score initial weight for balance */

  stp = (long)RDB[RES_N_BALA_SRC];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
  AddBuf(-wgt0, 1.0, stp, 0, -1, BALA_N_SRC_VR, 1);

  /* Check if source is empty */

  if (nsrc == 0)
    {
      /* Print warning */

      if ((RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DELDYN) &&
          (RDB[DATA_DYN_DT] < ZERO))
        Warn(FUNCTION_NAME, "Source is empty");

      /* Exit subroutine */

      return -1;
    }

  /* Check extinction and divergence */
  
  if (wgt0/RDB[DATA_SRC_POP] < 1E-16)
    {
      /* Reset id */

      id = 0;

      /* Pointer to first item after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Loop over source */

      while (ptr > VALID_PTR)
        {
          /* Get poiner */

          part = ptr;

          /* Pointer to Next */

          ptr = NextItem(ptr);

          /* Remove particle */

          RemoveItem(part);

          /* Put particle back to stack */

          ToStack(part, id++);

          /* Check id */

          if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
            id = 0;
        }

      /* Exit subroutine */

      return -1;
    }
  else if (wgt0/RDB[DATA_SRC_POP] > 1E+16)
    Error(0,"Population exceeds maximum at %1.2E seconds, adjust time cut-off",
          WDB[DATA_TIME_CUT_TMAX]);

  /* Set batch size */

  if (rr == YES)
    nbatch = (long)RDB[DATA_DYN_SRC_PREV_POP];
  else
    nbatch = (long)RDB[DATA_SRC_POP];
 
 /* Try to get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
    {
      /* If precursor detector is in use, population control is done in */
      /* ResizeDynSrc and SampleDelnu */

      /* Set nbatch to equal nsrc to not do population control again */
      /* Due to the stochastic nature of the population control, nsrc might */
      /* not be exactly DATA_SRC_POP */

      nbatch = nsrc;
    }

  /***************************************************************************/

  /***** Population control **************************************************/

  /* Get next rng index and update */

  n = (long)RDB[DATA_NHIST_TOT];
  WDB[DATA_NHIST_TOT] = RDB[DATA_NHIST_TOT] + 1.0;

  /* Init random number sequence (NOTE: tämä täytyy tehdä siksi, */
  /* että jos transport-sykli menee läpi ilman että thread 0:lle */
  /* tulee yhtään historiaa, sen seedi jää alustamatta). */

  seed = ReInitRNG(n);
  SEED[0] = seed;

  /* Reset weight and number of particles */

  wgt = wgt0;
  np = nsrc;

  /* Compare population size to batch size given in input */

  if (nsrc < nbatch)
    {
      /* Calculate multiplication */

      P = ((double)nbatch)/((double)nsrc);
      mul = (long)P;
      P = P - (double)mul;

      /* Pointer to first neutron after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Loop over population and sample particles for duplication */

      for (n = 0; n < nsrc; n++)
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

              part = DuplicateParticle(ptr, 0);

              /* Add to weight and population size */

              wgt = wgt + RDB[ptr + PARTICLE_WGT];
              np = np + 1;

              /* Score population balance */
              
              stp = (long)RDB[RES_N_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
              AddBuf(1.0, 1.0, stp, 0, -1, BALA_N_SRC_VR, 0);

              /* Put particle in source */

              AddItem(DATA_PART_PTR_SOURCE, part);
            }

          /* Next particle */

          ptr = NextItem(ptr);
        }
    }
  else if (nsrc > nbatch)
    {
      /* Calculate probability */

      P = 1.0 - ((double)nbatch)/((double)nsrc);

      /* Pointer to first neutron after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Reset id */

      id = 0;

      /* Loop over population and sample particles for duplication */

      for (n = 0; n < nsrc; n++)
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

              wgt = wgt - RDB[part + PARTICLE_WGT];
              np = np - 1;

              /* Score population balance */

              stp = (long)RDB[RES_N_BALA_LOSS];
              CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
              AddBuf(1.0, 1.0, stp, 0, -1, BALA_N_LOSS_CUT, 0);

              /* Put particle back to stack */

              ToStack(part, id++);

              /* Check id */

              if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
                id = 0;
            }
        }
    }

  /* Compare population size to minimum and maximum */

  if (nsrc < (long)RDB[DATA_DYN_POP_MIN])
    WDB[DATA_DYN_POP_MIN] = (double)nsrc;
  if (nsrc > (long)RDB[DATA_DYN_POP_MAX])
    WDB[DATA_DYN_POP_MAX] = (double)nsrc;

  /* Check */

  if (np < 1)
    return -1;

  /***************************************************************************/

  /***** Normalize source ****************************************************/

  /* Pointer to first item after dummy */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  ptr = NextItem(ptr);

  /* Normalize weights */

  while(ptr > VALID_PTR)
    {
      /* Normalize */

      WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*wgt0/wgt;

      /* Next */

      ptr = NextItem(ptr);
    }

  /***************************************************************************/

  /***** Russian roulette ****************************************************/

  /* Check option */

  if (rr == YES)
    {
      /* Pointer to first item after dummy */
      
      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);
      
      /* Reset OpenMP id */
      
      id = 0;
      
      /* Normalize weights */
      
      while (ptr > VALID_PTR)
        {
          /* Particle type */
          
          if ((type = (long)RDB[ptr + PARTICLE_TYPE]) != PARTICLE_TYPE_NEUTRON)
            break;
          
          /* Get pointer and next */
          
          part = ptr;
          ptr = NextItem(ptr);
          
          /* Get particle weight */
          
          wgt = RDB[part + PARTICLE_WGT];
          CheckValue(FUNCTION_NAME, "wgt", "", wgt, ZERO, INFTY);
          
          /* Check weight */

          if (wgt < wmin)
            {
              /* Calculate survival probability */

              P = wgt/wmin;
              CheckValue(FUNCTION_NAME, "P", "", P, ZERO, 1.0);

              /* Russian roulette */
              
              if (RandF(0) < P)
                {
                  /* Increase weight */
                  
                  WDB[part + PARTICLE_WGT] = wgt/P;
                }
              else
                {
                  /* Remove particle */
                  
                  RemoveItem(part);
                  
                  /* Put particle back to stack */
                  
                  ToStack(part, id++);
                  
                  /* Check id */
                  
                  if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
                    id = 0;
                  
                  /* Score cut-off */
                  
                  stp = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
                  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
                  AddBuf1D(1.0, wgt, stp, 0, 0);              
                  
                  /* Score particle balance */

                  stp = (long)RDB[RES_N_BALA_LOSS];
                  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
                  AddBuf(1.0, 1.0, stp, 0, -1, BALA_N_LOSS_CUT, 0);
                }
            }
        }
    }

  /***************************************************************************/

  /***** Set particle indexes ************************************************/

  /* Pointer to source buffer */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get next rng index */

  n = (long)RDB[DATA_NHIST_TOT];

  /* Reset total weight and source size */

  wgt = 0.0;
  nsrc = 0;

  /* Get pointer to last particle */

  ptr = LastItem(ptr);

  /* Loop over list and set indexes */

  while (1 != 2)
    {
      /* Break if dummy */

      if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
        break;

      /* Put index */

      WDB[ptr + PARTICLE_RNG_IDX] = (double)(++n);

      /* Add to total weight and source size */

      wgt = wgt + RDB[ptr + PARTICLE_WGT];
      nsrc = nsrc + 1;

      /* Score initial source rate and source weight for */
      /* interval in case of dynamic mode                */
      /* For MODE_SRC these are scored in samplesrcpoint */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
        {
          /* Check particle type */

          if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_GAMMA)
            {
              /* Score source rate */

              stp = (long)RDB[RES_TOT_PHOTON_SRCRATE];
              CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
              AddBuf1D(1.0, wgt, stp, 0, 0);
            }
          else
            {
              /* Score source rate */

              stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
              CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
              AddBuf1D(1.0, RDB[ptr + PARTICLE_WGT], stp, 0, 0);

              /* Score initial source weight */

              stp = (long)RDB[RES_INI_SRC_WGT];
              CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
              AddBuf1D(1.0, RDB[ptr + PARTICLE_WGT], stp, 0, 0);
            }
        }

      /* Next particle */

      ptr = PrevItem(ptr);
    }
  
  /* Check weight */

  if ((rr == NO) && (fabs(wgt/wgt0 - 1.0) > 1E-6))
    Die(FUNCTION_NAME, "Mismatch in weight");

  /* Score final weight for balance */

  stp = (long)RDB[RES_N_BALA_SRC];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
  AddBuf(wgt, 1.0, stp, 0, -1, BALA_N_SRC_VR, 1);

  /* Set new batch size */

  WDB[DATA_DYN_SRC_PREV_POP] = (double)nsrc;

  /* Put number of histories */

  WDB[DATA_NHIST_TOT] = (double)n;

  /* Score mean population size */

  ptr = (long)RDB[RES_MEAN_POP_SIZE];
  AddStat((double)nsrc, ptr, 0);

  /* Reset minimum size */
  
  min = 10000000000;

  /* Pointer to stack */

  ptr = (long)RDB[DATA_PART_PTR_MIN_NSTACK];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Get minimum stack size */

  for (n = 0; n < (long)RDB[DATA_OMP_MAX_THREADS]; n++)
    if (((m = (long)GetPrivateData(ptr, n)) > 0) && (m < min))
      min = m;

  /* Check minimum stack size */

  if (RDB[DATA_OMP_MAX_THREADS]*min/RDB[DATA_PART_ALLOC_N] < 0.2)
    {
      /* Allow memory allocation */

      Mem(MEM_ALLOW);

      /* Allocate 20% more particles */

      AllocParticleStack(PARTICLE_TYPE_NEUTRON,
                         (long)(0.2*RDB[DATA_PART_ALLOC_N]));

      /* Disallow memory allocation */

      Mem(MEM_DENY);
    }

  /* Re-distribute stacks */

  ReDistributeStacks();

  /* Exit */

  return 0;

  /***************************************************************************/
}

/*****************************************************************************/
