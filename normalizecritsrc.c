/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizecritsrc.c                             */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Sets up normalized fission source for criticality source     */
/*              simulation                                                   */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizeCritSrc:"

/*****************************************************************************/

void NormalizeCritSrc()
{
  long ptr, pts, pos, n, mat, stp, nsrc, nbatch, id, fmx, idx, i, neig, min;
  double wgt, w0, keff, kw, P, kp, *wgtg, *wgtg0;

  /***************************************************************************/

  /***** Move banked neutrons to source **************************************/

  /* Get pointer to source */

  pos = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);

  /* Check that source is empty */

  if (ListSize(pos) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Reset total weight and source size */

  wgt = 0.0;
  nsrc = 0;

  /* Number of eigenvalue simulations */

  neig = (long)RDB[DATA_N_POP_EIG];
  CheckValue(FUNCTION_NAME, "(neig)", "", neig, 1, 10000);

  /* Temporary array for data */

  wgtg = WorkArray(DATA_PTR_WORK_GRID1, DATA_ARRAY, neig, 0);
  memset(wgtg, 0.0, neig*sizeof(double));

  wgtg0 = WorkArray(DATA_PTR_WORK_GRID2, DATA_ARRAY, neig, 0);
  memset(wgtg0, 0.0, neig*sizeof(double));

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      while ((ptr = FromBank(id)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Get neutron weight */

          w0 = RDB[ptr + PARTICLE_WGT];

          /* Add to total weight and source size */

          wgt = wgt + w0;
          nsrc = nsrc + 1;

          /* Get eigenvalue simulation index */

          i = (long)RDB[ptr + PARTICLE_EIG_IDX];
          CheckValue(FUNCTION_NAME, "i", "", i, 0, 9999);

          /* Add to weight */

          wgtg[i] = wgtg[i] + w0;

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

  /* Add to gpop size */

  WDB[DATA_GROW_POP_NHIST] = RDB[DATA_GROW_POP_NHIST] + (double)nsrc;

  /* Check weights */

  for (i = 0; i < neig; i++)
    if (wgtg[i] == 0.0)
      Error(0, "Unable to initiate self-sustaining chain reaction");

#ifdef MPI

  /* Check for domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Reduce weights */

      MPITransfer(wgtg, wgtg0, neig, 0, MPI_METH_RED);

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Broadcast weights */

      MPITransfer(wgtg0, NULL, neig, 0, MPI_METH_BC);

      /* Copy */

      for (n = 0; n < neig; n++)
        wgtg[n] = wgtg0[n];
    }

#endif

  /* Set batch size */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    nbatch = (long)RDB[DATA_CRIT_POP];
  else if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    {
      /* Check if growing population is used */

      if ((long)RDB[DATA_GROW_POP_SIM] == YES)
        {
          /* Calculate parameter M during the initial cycle */
          /* M = (initial batch/parameter c)^2              */

          if ((long)RDB[DATA_CYCLE_IDX] == 0)
            WDB[DATA_GROW_POP_M] =
              (RDB[DATA_SIMUL_BATCH_SIZE]/RDB[DATA_GROW_POP_C])
              *(RDB[DATA_SIMUL_BATCH_SIZE]/RDB[DATA_GROW_POP_C]);

          /* Store the non-modified batch size for keff calculation */

          WDB[DATA_GROW_POP_PREV] = RDB[DATA_SIMUL_BATCH_SIZE];

          /* Modify the batch size used in normalization */

          nbatch = (long)(RDB[DATA_GROW_POP_C]
                          *round(sqrt(RDB[DATA_GROW_POP_NHIST]
                                      + RDB[DATA_GROW_POP_M])));

          /* Compare to maximum */

          if (nbatch > (long)RDB[DATA_GROW_POP_MAX_POP])
            nbatch = (long)RDB[DATA_GROW_POP_MAX_POP];

          /* Put the modified value in the data block */

          WDB[DATA_SIMUL_BATCH_SIZE] = (double)nbatch;
        }
      else
        nbatch = (long)RDB[DATA_SIMUL_BATCH_SIZE];
    }
  else
    nbatch = (long)RDB[DATA_CRIT_POP];

  /* Set cycle batch size and weight */

  WDB[DATA_CYCLE_BATCH_SIZE] = (double)nsrc;

  /* Score mean population size */

  ptr = (long)RDB[RES_MEAN_POP_SIZE];
  AddStat((double)nsrc, ptr, 0);

  /* Score mean population weight */

  ptr = (long)RDB[RES_MEAN_POP_WGT];
  AddStat(wgt, ptr, 0);

  /* keff's for independent eigenvalue simulations */

  ptr = (long)RDB[DATA_PTR_CYCLE_EIG_KEFF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  WDB[DATA_CYCLE_KEFF] = 0.0;

  for (n = 0; n < neig; n++)
    {
      /* Get previous cycle-wise k-eff */

      keff = RDB[ptr + n];

      /* Check growing population mode */

      if ((long)RDB[DATA_GROW_POP_SIM] == YES)
        {
          /* Normalize keff to batch size set for previous cycle      */
          /* NOTE: DATA_GROW_POP_PREV should be storet separately for */
          /* different eigenvalue simulations. Now the calculation is */
          /* terminated if neig > 1. */

          if (neig > 1)
            Die(FUNCTION_NAME, "multiple eigenvalue simulations with gpop");
          else
            keff = keff*wgtg[n]/(RDB[DATA_GROW_POP_PREV]);
        }
      else
        keff = keff*wgtg[n]*((double)neig)/((double)(nbatch));

      /* Put value */

      WDB[ptr + n] = keff;

      /* Add to mean */

      WDB[DATA_CYCLE_KEFF] = RDB[DATA_CYCLE_KEFF] + keff;
    }

  /* Calculate mean */

  WDB[DATA_CYCLE_KEFF] = RDB[DATA_CYCLE_KEFF]/((double)neig);
  keff = RDB[DATA_CYCLE_KEFF];

  /* Check if Wieland shift is used */

  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      /* Avoid compiler error */

      kw = -1.0;
      kp = -1.0;
      P = -1.0;

      /* Check mode */

      if ((long)RDB[DATA_WIELANDT_MODE] == WIELANDT_MODE_FIX_K)
        {
          /* Get user-specified k-eff */

          kw = RDB[DATA_WIELANDT_KEFF];
          CheckValue(FUNCTION_NAME, "kw", "", kw, ZERO, INFTY);

          /* Calculate backtransformed k-eff */

          kp = kw*keff/(kw + keff);
          CheckValue(FUNCTION_NAME, "kp", "", kp, ZERO, INFTY);

          /* Calculate probability of banking the neutron */

          P = 1.0 - kp/kw;
          CheckValue(FUNCTION_NAME, "P", "", P, ZERO, 1.0);
        }
      else if ((long)RDB[DATA_WIELANDT_MODE] == WIELANDT_MODE_FIX_P)
        {
          /* Get user-specified probability */

          P = RDB[DATA_WIELANDT_P];
          CheckValue(FUNCTION_NAME, "P", "", P, ZERO, 1.0);

          /* Calculate backtransformed k-eff */

          kp = keff*P;
          CheckValue(FUNCTION_NAME, "kp", "", kp, ZERO, INFTY);

          /* Calculate k-eff */

          kw = kp/(1.0 - P);
          CheckValue(FUNCTION_NAME, "kw", "", kw, ZERO, INFTY);
        }
      else
        Die(FUNCTION_NAME, "Invalid Wielandt mode");

      /* Store values */

      WDB[DATA_WIELANDT_KEFF] = kw;
      WDB[DATA_WIELANDT_KP] = kp;
      WDB[DATA_WIELANDT_P] = P;

      /* Store kw and P */

      ptr = (long)RDB[RES_WIELANDT_K];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(kw, 1.0, ptr, 0, 0);

      ptr = (long)RDB[RES_WIELANDT_P];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(P, 1.0, ptr, 0, 0);

      /* Store backtransformed keff */

      ptr = (long)RDB[RES_ANA_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(kp, 1.0, ptr, 0, 0);
    }
  else
    {
      /* Store keff */

      ptr = (long)RDB[RES_ANA_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(keff, 1.0, ptr, 0, 0);
    }

  /***************************************************************************/

  /***** Re-normalize source *************************************************/

  /* Reset prompt and delayed weights */

  WDB[DATA_CYCLE_PROMPT_WGT] = 0.0;
  WDB[DATA_CYCLE_DELAYED_WGT] = 0.0;

  /* Pointer to first item after dummy */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  ptr = NextItem(ptr);

  /* Loop over source */

  while (ptr > VALID_PTR)
    {
      i = (long)RDB[ptr + PARTICLE_EIG_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, 9999);

      /* Normalize weight */

      w0 = RDB[ptr + PARTICLE_WGT]*((double)nbatch)/((double)neig)/wgtg[i];

      /* Set value */

      WDB[ptr + PARTICLE_WGT] = w0;

      /* Check delayed neutron group */

      if ((long)RDB[ptr + PARTICLE_DN_GROUP] == 0)
        WDB[DATA_CYCLE_PROMPT_WGT] = RDB[DATA_CYCLE_PROMPT_WGT] + w0;
      else
        WDB[DATA_CYCLE_DELAYED_WGT] = RDB[DATA_CYCLE_DELAYED_WGT] + w0;

      /* Get material pointer (may be null for initial source) */

      mat = (long)RDB[ptr + PARTICLE_PTR_MAT];

      /* Score source rate */

      if ((mpiid == 0) || ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO))
        {
          stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
          CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
          AddBuf1D(1.0, w0, stp, 0, 0);

          /* Score source rate in fissile and non-fissile materials */

          if (mat > VALID_PTR)
            {
              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                AddBuf1D(1.0, w0, stp, 0, 1);
              else
                AddBuf1D(1.0, w0, stp, 0, 2);
            }
        }

      /* Score fission matrix source term */

      if ((idx = (long)RDB[ptr + PARTICLE_FMTX_IDX]) > -1)
        {
          /* Get pointer */

          fmx = (long)RDB[DATA_PTR_FMTX];
          CheckPointer(FUNCTION_NAME, "(fmx)", DATA_ARRAY, fmx);

          stp = (long)RDB[fmx + FMTX_PTR_SRC];
          CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

          /* Score total */

          AddBuf(1.0, w0, stp, 0, -1, 0, idx);

          /* Score prompt or delayed */

          if ((long)RDB[ptr + PARTICLE_DN_GROUP] == 0)
            AddBuf(1.0, w0, stp, 0, -1, 1, idx);
          else
            AddBuf(1.0, w0, stp, 0, -1, 2, idx);
        }

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Calculate entropies */

  CalculateEntropies();

  /***************************************************************************/

  /***** Set particle indexes ************************************************/

  /* Get next rng index */

  n = (long)RDB[DATA_NHIST_TOT];

  /* Reset MPI index counter */

  id = 0;

  /* Reset weight */

  wgt = 0.0;

  /* Pointer to source buffer */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Pointer to particle balance */

  pts = (long)RDB[RES_N_BALA_SRC];
  CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

  /* Get pointer to last particle */

  ptr = LastItem(ptr);

  /* Loop over list and set indexes */

  while(1 != 2)
    {
      /* Break if dummy */

      if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
        break;

      /* Put index */

      WDB[ptr + PARTICLE_RNG_IDX] = (double)(++n);
      WDB[ptr + PARTICLE_HISTORY_IDX] = RDB[ptr + PARTICLE_RNG_IDX];

      /* Particle weight */

      w0 = RDB[ptr + PARTICLE_WGT];

      /* Score particle balance */

      AddBuf(1.0, 1.0, pts, id, -1, BALA_N_SRC_FISS, 0);
      AddBuf(w0, 1.0, pts, id, -1, BALA_N_SRC_FISS, 1);

      /* Add to weight */

      wgt = wgt + w0;

      /* Put MPI id */

      if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
        WDB[ptr + PARTICLE_MPI_ID] = (double)mpiid;
      else
        WDB[ptr + PARTICLE_MPI_ID] = (double)(id++);

      /* Check id */

      if (id == mpitasks)
        id = 0;

      /* Update number of histories */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_NHIST_CYCLE] + 1.0;

      /* Next particle */

      ptr = PrevItem(ptr);
    }

#ifdef MPI

  /* Check for domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Reduce total weight */

      w0 = 0;
      MPITransfer(&wgt, &w0, 1, 0, MPI_METH_RED);

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Broadcast weight */

      MPITransfer(&w0, NULL, neig, 0, MPI_METH_BC);

      /* Copy */

      wgt = w0;
    }

#endif

  /* Check total weight */

  if(fabs(wgt/((double)nbatch) - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Total weight not preserved");

  /* Put number of histories */

  WDB[DATA_NHIST_TOT] = (double)n;

  /* Reset minimum size */

  min = 10000000000;

  /* Pointer to stack */

  ptr = (long)RDB[DATA_PART_PTR_MIN_NSTACK];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Get minimum stack size */

  for (n = 0; n < (long)RDB[DATA_OMP_MAX_THREADS]; n++)
    if ((i = (long)GetPrivateData(ptr, n)) > 0)
      if (i < min)
        min = i;

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

  /* Re-distribute stacks (NOTE: Tää if-lauseke korjaa sellaisen ongelman    */
  /* että toi ReDistributeStacks() antaa erroria jos track-plot animaatiota  */
  /* yrittää tehdä kriittisyyslähdemoodissa useammalla OpenMP -threadilla.   */
  /* Perimmäinen syy pitäisi kuitenkin selvittää (JLe / 2.1.28 / 29.9.2016). */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] < 1)
    ReDistributeStacks();

  /* Check that source is sorted */

#ifdef DEBUG

  /* Check reproducibility */

  if ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == YES)
    {
      /* Pointer to first after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Loop over source */

      while (ptr > VALID_PTR)
        {
          /* Compare */

          if ((pos = NextItem(ptr)) > VALID_PTR)
            if (RDB[pos + PARTICLE_RNG_IDX] >= RDB[ptr + PARTICLE_RNG_IDX])
              Die(FUNCTION_NAME, "Sorting failed");

          /* Next */

          ptr = NextItem(ptr);
        }
    }

#endif

  /* Plot source point distribution */

  GeometryPlotter(NO);

  /***************************************************************************/
}

/*****************************************************************************/
