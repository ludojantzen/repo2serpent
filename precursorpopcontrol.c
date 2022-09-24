/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : precursorpopcontrol.c                          */
/*                                                                           */
/* Created:       2015/09/15 (VVa)                                           */
/* Last modified: 2017/12/02 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Does population control for precursors in PSOURCE, all       */
/*              precursors are normalized for the same emission during       */
/*              upcoming time-interval                                       */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrecursorsPopControl:"

/*****************************************************************************/

void PrecursorPopControl()
{
  long ptr, part, n, m, nsrc, nbatch, mul, N, np, id, loc0, gbin, ng, min;
  long *groupnums0, *groupnums1;
  double wgt0, wgt2, wgt1, wgt, P, emit0, emit, emit1, emit2, t0, t1, dt;
  double lambda, aveemit;
  double *groupwgts0, *groupwgts1;

  /***************************************************************************/

  if ((long)RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_MESH)
    return;

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "precursorpopcontrol.c -->\n");
#endif

  /***** Calculate initial source weight and size*****************************/

  /* Reset total weight and source size */

  wgt0 = 0.0;
  emit0 = 0.0;

  nsrc = 0;

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  dt = t1 - t0;

  /* Get number of delayed neutron groups */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Allocate memory for temporary lists */

  groupwgts0  = (double *)Mem(MEM_ALLOC, ng, sizeof(double));
  groupwgts1  = (double *)Mem(MEM_ALLOC, ng, sizeof(double));
  groupnums0  = (long *)Mem(MEM_ALLOC, ng, sizeof(double));
  groupnums1  = (long *)Mem(MEM_ALLOC, ng, sizeof(double));

  /* Loop over source to calculate initial source size and weight */

  /* Pointer to dummy */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Pointer to first item after dummy */

  ptr = NextItem(ptr);

  /* Loop over source */

  while (ptr > VALID_PTR)
    {

      /* Get particle weight */

      wgt = RDB[ptr + PARTICLE_WGT];

      /* Get delayed neutron group bin */

      gbin = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Get decay constant */

      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Add to total weight and source size */

      wgt0 += wgt;
      nsrc += 1;

      /* Add to initial group weight */

      groupwgts0[gbin] += wgt;

      /* Add to initial group population */

      groupnums0[gbin] += 1;

      /* Add to total activity */

      emit0 += wgt*(1-exp(-lambda*dt));

      /* Next particle */

      ptr = NextItem(ptr);
    }

  /* Check, whether the precursor source is empty */

  if (nsrc == 0)
    {
      /* Put total weight to emit on next interval */

      WDB[loc0 + PRECDET_W_EMIT] = 0.0;

      /* Warn about an empty precursor source */

      Warn(FUNCTION_NAME, "Precursor source empty at population control");

      /* Free temporary lists */

      Mem(MEM_FREE, groupwgts0);
      Mem(MEM_FREE, groupwgts1);
      Mem(MEM_FREE, groupnums0);
      Mem(MEM_FREE, groupnums1);

      /* Return */

      return;

    }

  /* Normalize to average emission on upcoming time interval */
  /* Current average emission */

  aveemit = emit0/(double)nsrc;

#ifdef DNPRINT
  fprintf(outp, "Average emission was %E neutrons during interval of %E s\n",aveemit, dt);
#endif

  /***************************************************************************/
  /***** Population control **************************************************/

  /* Normalize to the wanted population */

  nbatch = (long)RDB[DATA_SRC_POP]*(long)RDB[DATA_PREC_SRC_FACT];

#ifdef MPI_MODE1

      /* Check number of tasks */

      if (mpitasks > 1)
        nbatch = (long)((double)nbatch/(double)mpitasks);

#endif

  /* Average emission for the wanted number of neutrons */

  aveemit = emit0/(double)nbatch;

  /* Store initial emission */

  emit1 = emit0;

  /* Store initial weight */

  wgt1 = wgt0;

  /* Store initial number of precursors */

  np = nsrc;

  /* Store current group weights and populations*/

  for (n = 0; n < ng; n++)
    {
      groupwgts1[n] = groupwgts0[n];
      groupnums1[n] = groupnums0[n];
    }

  /* Reset id */

  id = 0;

  /* Loop over source */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  ptr = NextItem(ptr);

  /* Loop over population and sample particles for duplication */

  for (n = 0; n < nsrc; n++)
    {
      /* Get decay constant */

      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Get DN group */

      gbin = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Calculate number of neutrons to emit */

      emit = RDB[ptr + PARTICLE_WGT]*(1-exp(-lambda*dt));

      /* remove initial weight from total weight */

      wgt1 = wgt1 - RDB[ptr + PARTICLE_WGT];

      /* remove initial weight from group weight */

      groupwgts1[gbin] -= RDB[ptr + PARTICLE_WGT];

      /* Remove from group population */

      groupnums1[gbin] -= 1;

      /* remove initial emission from total emission */

      emit1 = emit1 - emit;

      if (emit > aveemit)
        {
          /* Sample splitting */
          /* Calculate multiplication */

          P = emit/aveemit;
          mul = (long)P;
          P = P - (double)mul;

          /* Sample multiplication */

          if (drand48() < P)
            N = mul;
          else
            N = mul - 1;

          /* Put new weight of particle */

          wgt = aveemit/emit*RDB[ptr + PARTICLE_WGT];

          WDB[ptr + PARTICLE_WGT] = wgt;

          /* Add new weight to total weight */

          wgt1 += wgt;

          /* Add new weight to group weight */

          groupwgts1[gbin] += wgt;

          /* Add to group population */

          groupnums1[gbin] += 1;

          /* Add new emission to total emission */

          emit1 += wgt*(1-exp(-lambda*dt));

          /* Create N particles */

          for (m = 0; m < N; m++)
            {

              /* Duplicate neutron */

              part = DuplicateParticle(ptr, id++);

              /* Check id */

              if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
                id = 0;

              /* Add to weight and population size */

              wgt1 += wgt;

              /* Add new weight to group weight */

              groupwgts1[gbin] += wgt;

              /* Add to group population */

              groupnums1[gbin] += 1;

              /* Add new emission to total emission */

              emit1 += wgt*(1-exp(-lambda*dt));

              np = np + 1;

              /* Put particle in source */
              /* Will be sorted later   */

              AddItem(DATA_PART_PTR_PSOURCE, part);
            }

          ptr = NextItem(ptr);

        }
      /* Never remove the last precursor from a group, because the weights */
      /* are scaled separately for each group and this would lead to lost  */
      /* weight as the weight in this group could not be conserved */
      else if ((emit < aveemit) && (groupnums1[gbin] != 0))
        {
          /* Get pointer and next */

          part = ptr;
          ptr = NextItem(ptr);

          /* Sample russian roulette */

          if (drand48() < emit/aveemit)
            {
              /* Passed russian roulette */
              /* Increase particle weight */
              wgt = RDB[part + PARTICLE_WGT]*aveemit/emit;

              /* Store new particle weight */

              WDB[part + PARTICLE_WGT] = wgt;

              /* Add to weight and population size */

              wgt1 += wgt;

              /* Add new weight to group weight */

              groupwgts1[gbin] += wgt;

              /* Add to group population */

              groupnums1[gbin] += 1;

              /* Add new emission to total emission */

              emit1 += wgt*(1-exp(-lambda*dt));

            }
          else
            {
              /* Failed russian roulette, remove from source */

              RemoveItem(part);

              /* subtract from population size */

              np = np - 1;

              /* Put particle back to stack */

              ToStack(part, id++);

              /* Check id */

              if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
                id = 0;

            }

        }
      else
        {
          /* Put particle back as is */

          wgt = RDB[ptr + PARTICLE_WGT];

          /* Add to weight and population size */

          wgt1 += wgt;

          /* Add new weight to group weight */

          groupwgts1[gbin] += wgt;

          /* Add to group population */

          groupnums1[gbin] += 1;

          /* Add new emission to total emission */

          emit1 += wgt*(1-exp(-lambda*dt));

          /* Next particle */

          ptr = NextItem(ptr);
        }

    }


#ifdef DNPRINT

  fprintf(outp, "%ld precursors in source (was %ld), %ld wanted, emit is now %E was %E\n",
          np, nsrc, nbatch, emit1, emit0);

  fprintf(outp, "Wgt0 %E, wgt1 %E\n", wgt0, wgt1);

  for (n = 0; n < ng; n++)
    {
      if (groupwgts1[n] < 0)
        {
          groupwgts1[n] = 0;
          fprintf(outp, "group %ld, num0 %ld, num1 %ld, wgt0 %E, wgt1 %E, cannot scale\n", n, groupnums0[n], groupnums1[n], groupwgts0[n], groupwgts1[n]);
        }
      else
        fprintf(outp, "group %ld, num0 %ld, num1 %ld, wgt0 %E, wgt1 %E, should scale by %f\n", n, groupnums0[n], groupnums1[n], groupwgts0[n], groupwgts1[n], groupwgts0[n]/groupwgts1[n]);
    }

#endif

  /***** Normalize source ****************************************************/

  /* Pointer to first item after dummy */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  ptr = NextItem(ptr);

  emit2 = 0.0;

  /* Normalize weights */

  wgt2 = 0.0;

  while(ptr > VALID_PTR)
    {
      /* Get decay constant */

      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Get DN group */

      gbin = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Scale weight to conserve weight in each group */

      if (groupwgts1[gbin] > 0)
        {
          wgt = RDB[ptr + PARTICLE_WGT]*groupwgts0[gbin]/groupwgts1[gbin];
        }
      else if (groupwgts1[gbin] == 0)
        {
          wgt = 0;
        }
      else
        {
          Die(FUNCTION_NAME, "Invalid group weight");
        }

      /* Put new weight */

      WDB[ptr + PARTICLE_WGT] = wgt;

      /* Calculate proportion of total emission by this precursor */
      /* Total emission should be conserved also */

      WDB[ptr + PARTICLE_U] =  wgt*(1-exp(-lambda*dt))/emit0;

      /* Add to total emission */

      emit2 += wgt*(1-exp(-lambda*dt));

      /* Add to total weight */

      wgt2 += wgt;

      /* Next precursor */

      ptr = NextItem(ptr);
    }

  /* Store total weight to emit */

  WDB[loc0 + PRECDET_W_EMIT] = emit2/RDB[DATA_NORM_COEF_N];

  /* Calculate new average emission during next interval */

  aveemit = emit2/(double)np;

  /* Store new average emission during next interval */
  /* Will be used in precdet as a threshold for Russian Roulette */

  WDB[loc0 + PRECDET_AVE_EMIT] = aveemit;

  /* Free temporary lists */

  Mem(MEM_FREE, groupwgts0);
  Mem(MEM_FREE, groupwgts1);
  Mem(MEM_FREE, groupnums0);
  Mem(MEM_FREE, groupnums1);

#ifdef DNPRINT

  fprintf(outp, "Emit0 %E, emit2 %E\n", emit0, emit2);

  fprintf(outp, "Neutrons to emit %E\n", emit2);
  fprintf(outp, "Total weight to emit %E\n", emit2/RDB[DATA_NORM_COEF_N]);
#endif

  /* Check weight */

  if (fabs(wgt2/wgt0 - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Mismatch in weight %E %%", (wgt2/wgt0 - 1.0)*100.0);

  /* Check emission */

  if (fabs(emit2/emit0 - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Mismatch in emission %E %%", (emit2/emit0 - 1.0)*100.0);

  /* Reset minimum size */
  
  min = 10000000000;

  /* Pointer to stack */

  ptr = (long)RDB[DATA_PART_PTR_MIN_PSTACK];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Loop over OpenMP threads */

  for (n = 0; n < (long)RDB[DATA_OMP_MAX_THREADS]; n++)
    if ((long)GetPrivateData(ptr, n) < min)
      min = (long)GetPrivateData(ptr, n);
  
  /* Check minimum stack size */

  if (RDB[DATA_OMP_MAX_THREADS]*min/RDB[DATA_PART_ALLOC_P] < 0.2)
    {
      /* Allow memory allocation */

      Mem(MEM_ALLOW);

      /* Allocate 20% more particles */

      AllocParticleStack(PARTICLE_TYPE_PRECURSOR,
                         (long)(0.2*RDB[DATA_PART_ALLOC_P]));

      /* Disallow memory allocation */

      Mem(MEM_DENY);
    }

  /* Re-distribute stacks */

  /* Stacks will be redistributed at normalizedynsrc.c */
  /*
  ReDistributeStacks();
  */

#ifdef DNPRINT
  fprintf(outp, "<-- precursorpopcontrol.c\n\n");
#endif

  /***************************************************************************/
}

/*****************************************************************************/
