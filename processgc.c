/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processgc.c                                    */
/*                                                                           */
/* Created:       2011/06/10 (JLe)                                           */
/* Last modified: 2019/12/18 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Process data needed for group constant generation            */
/*                                                                           */
/* Comments: - Noiden score-pointtereiden järjestys vaikuttaa silmukkaan     */
/*             ainakin coefoutput.c:ssä. Tätä aliohjelmaa pitää kutsua vasta */
/*             processtats.c:n jälkeen.                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessGC:"

/*****************************************************************************/

void ProcessGC()
{
  long gcu, ptr, adf, ppw, alb, loc0, loc1, uni, surf, lat, nfg, nmg, n, m;
  long mat, mat0;
  double *E;

  /* Check option */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    {
      /* Reset pointers */

      WDB[DATA_MICRO_PTR_EGRID] = NULLPTR;
      WDB[DATA_ERG_FG_PTR_GRID] = NULLPTR;

      /* Exit subroutine */

      return;
    }

  /* Check batch sizes */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Get sizes */

      n = (long)RDB[DATA_CRIT_CYCLES];
      m = (long)RDB[DATA_MICRO_CALC_BATCH_SIZE];

      /* Compare to number of active cycles */

      if (n < m)
        Error(0, "Number of active cycles must be %ld or more "
              "when generating group constants.\nGroup constant generation "
              "can be switched off with \"set gcu -1\".", m);
      else if (n % m)
        Error(0, "Number of active cycles must be divisible by %ld "
              "when generating group constants.\nGroup constant generation "
              "can be switched off with \"set gcu -1\".", m);
    }
  else
    {
      /* Get sizes */

      n = (long)RDB[DATA_SRC_BATCHES];
      m = (long)RDB[DATA_MICRO_CALC_BATCH_SIZE];

      /* Compare to number of active cycles */

      if (n < m)
        Error(0, "Number of source batches must be %ld or more "
              "when generating group constants.\nGroup constant generation "
              "can be switched off with \"set gcu -1\".", m);
      else if (n % m)
        Error(0, "Number of source batches must be divisible by %ld "
              "when generating group constants.\nGroup constant generation "
              "can be switched off with \"set gcu -1\".", m);
    }

  fprintf(outp, "Processing data for group constant generation:\n\n");

  /***************************************************************************/

  /***** Few-group structure *************************************************/

  /* Check pre-defined group structure */

  if ((ptr = (long)RDB[DATA_ERG_FG_PTR_PREDEF]) > VALID_PTR)
    {
      /* Number of energy groups */

      nfg = (long)RDB[ptr + ENE_NB];

      /* Pointer to data */

      ptr = (long)RDB[ptr + ENE_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);

      /* Put pointer and number of groups */

      WDB[DATA_ERG_FG_PTR_GRID] = (double)ptr;
      WDB[DATA_ERG_FG_NG] = (double)nfg;
    }

  /* Get number of groups */

  nfg = (long)RDB[DATA_ERG_FG_NG] + 1;

  /* Check */

  if (nfg < 2)
    Error(0, "Invalid group structure for group constant generation");

  /* Allcate memory for temporary array */

  E = (double *)Mem(MEM_ALLOC, nfg, sizeof(double));

  /* Check pointer to array */

  if ((ptr = (long)RDB[DATA_ERG_FG_PTR_GRID]) < VALID_PTR)
    Error(0, "Missing few-group structure for group constant generation");
  else
    {
      /* Read data in array */

      for (n = 0; n < nfg; n++)
        E[n] = RDB[ptr++];
    }

  /* Check order */

  for (n = 1; n < nfg; n++)
    if (E[n] <= E[n - 1])
      Error(0, "Values in few-group structure must be in ascending order");

  /* Make energy grid */

  ptr = MakeEnergyGrid(nfg, 0, 0, -1, E, EG_INTERP_MODE_LIN);

  /* Put pointer */

  WDB[DATA_ERG_FG_PTR_GRID] = (double)ptr;

  /* Free temporary array */

  Mem(MEM_FREE, E);

  /* Check that transport correction is defined if TRC_TRANSPXS */
  /* Is used with homogeneous diffusion flux solver */

  if ((long)RDB[DATA_HOMOFLUX_DIFFCOEF] == HOMOFLUX_DIFF_TRC)
    if ((long)RDB[DATA_PTR_TRC0] < VALID_PTR)
      Error(0,
            "Transport correction must be defined for diffusion calculation");

  /* Check that TRC calculation is on if used with leakage correction */

  if (((long)RDB[DATA_PTR_TRC0] < VALID_PTR) &&
      ((long)RDB[DATA_B1_FM_DIFF] == FM_DIFF_TRC))
    Error(0, "Transport correction must be defined for leakage correction");

  /* Check that CMM calculation is on if used with leakage correction */

  if (((long)RDB[DATA_CMM_CALC] == NO) &&
      ((long)RDB[DATA_B1_FM_DIFF] == FM_DIFF_CMM))
    Die(0, "Cumulative migration method must be on for leakage correction");

  /***************************************************************************/

  /***** Process micro-group structure ***************************************/

  /* Find energy grid */

  if ((ptr = (long)RDB[DATA_PTR_ENE0]) < VALID_PTR)
    Error(0, "No energy group structures defined for micro-group calculation");
  else if ((ptr = SeekListStr(ptr, ENE_PTR_NAME,
                              GetText(DATA_MICRO_PTR_EGRID))) < VALID_PTR)
    Error(0, "Energy grid %s used for micro-group calculation is not defined",
          GetText(DATA_MICRO_PTR_EGRID));

  /* Pointer to grid structure */

  ptr = (long)RDB[ptr + ENE_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr3)", DATA_ARRAY, ptr);

  /* Put pointer */

  WDB[DATA_MICRO_PTR_EGRID] = (double)ptr;

  /* Number of points */

  nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

  /* Get pointer to values */

  loc0 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Pointer to few-group structure */

  ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

  /* Get pointer to values */

  loc1 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  /* Make sure that all points match */

  for (n = 1; n < nfg - 1; n++)
    {
      /* Find interval in micro-group structure */

      if ((m = SearchArray(&RDB[loc0], RDB[loc1 + n], nmg + 1)) > -1)
        {
          /* Pick closest match */

          if ((RDB[loc1 + n] - RDB[loc0 + m]) >
              (RDB[loc0 + m + 1] - RDB[loc1 + n]))
            m = m + 1;

          /* Put value */

          WDB[loc0 + m] = RDB[loc1 + n];
        }
    }

  /* Check that micro-group grid is still in ascending order */

  for (n = 1; n < nfg; n++)
    if (RDB[loc1 + n] <= RDB[loc1 + n - 1])
      Error(0, "Mismatch between micro and macro-group structure");

  /* Allocate memory for indexes */

  ptr = ReallocMem(DATA_ARRAY, nmg);
  WDB[DATA_MICRO_PTR_IDX_MAP] = (double)ptr;

  /* Map indexes */

  m = 0;
  for (n = 1; n < nfg - 1; n++)
    {
      /* Loop over micro-group structure */

      while (m < nmg)
        {
          /* Check boundary */

          if (RDB[loc0 + m] == RDB[loc1 + n])
            break;

          /* Put index */

          WDB[ptr + nmg - m - 1] = (double)(nfg - 1 - n);

          /* Next */

          m++;
        }

      /* Check match */

      if (m == nmg)
        Error(0,
              "Few-group boundary %1.5E not found in micro-group energy grid",
              RDB[loc1 + n]);
    }

  /* Put remaining */

  while (m < nmg)
    WDB[ptr + nmg - m++ - 1] = (double)(nfg -1 - n);

  /***************************************************************************/

  /***** Add ADF universes to GCU list ***************************************/

  /* NOTE: Jos näitä lisätään niin ne pitää ottaa huomioon myös */
  /* invokebranch.c:ssä */

  /* Loop over ADF's */

  adf = (long)RDB[DATA_PTR_ADF0];
  while (adf > VALID_PTR)
    {
      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(gcu + GCU_PTR_UNIV, adf + ADF_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu < VALID_PTR)
        {
          /* Not defined, add to list */

          gcu = NewItem(DATA_PTR_GCU0, GCU_BLOCK_SIZE);
          WDB[gcu + GCU_PTR_UNIV] = RDB[adf + ADF_PTR_GCU];
        }

      /* Next ADF */

      adf = NextItem(adf);
    }

  /***************************************************************************/

  /***** Add PPW universes to GCU list ***************************************/

  /* Loop over pin power distributions */

  ppw = (long)RDB[DATA_PTR_PPW0];
  while (ppw > VALID_PTR)
    {
      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(gcu + GCU_PTR_UNIV, ppw + PPW_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu < VALID_PTR)
        {
          /* Not defined, add to list */

          gcu = NewItem(DATA_PTR_GCU0, GCU_BLOCK_SIZE);
          WDB[gcu + GCU_PTR_UNIV] = RDB[ppw + PPW_PTR_GCU];
        }

      /* Next PPW */

      ppw = NextItem(ppw);
    }

  /***************************************************************************/

  /***** Add MDEP universes to GCU list **************************************/

  /* Loop over micro depletion definitions */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(gcu + GCU_PTR_UNIV, loc0 + MDEP_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu < VALID_PTR)
        {
          /* Not defined, add to list */

          gcu = NewItem(DATA_PTR_GCU0, GCU_BLOCK_SIZE);
          WDB[gcu + GCU_PTR_UNIV] = RDB[loc0 + MDEP_PTR_GCU];
        }

      /* Next definition */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Add ALB universes to GCU list ***************************************/

  /* Loop over ALB's */

  alb = (long)RDB[DATA_PTR_ALB0];
  while (alb > VALID_PTR)
    {
      /* Check if multiple albedos are defined */

      if (ListSize(alb) > 1)
        {
          /* Check direction */

          if ((long)RDB[alb + ALB_DIR] == 1)
            Note(alb,
                 "Albedo calculation with multiple surfaces may not work");
        }
      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(gcu + GCU_PTR_UNIV, alb + ALB_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu < VALID_PTR)
        {
          /* Not defined, add to list */

          gcu = NewItem(DATA_PTR_GCU0, GCU_BLOCK_SIZE);
          WDB[gcu + GCU_PTR_UNIV] = RDB[alb + ALB_PTR_GCU];
        }

      /* Next albedo */

      alb = NextItem(alb);
    }

  /***************************************************************************/

  /***** Link universes ******************************************************/

  /* Check duplicates */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Loop over remaining */

      ptr = NextItem(gcu);
      while (ptr > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(ptr + GCU_PTR_UNIV, gcu + GCU_PTR_UNIV))
            Error(0, "Universe %s included multiple times in GC generation",
                  GetText(ptr + GCU_PTR_UNIV));

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next */

      gcu = NextItem(gcu);
    }

  /* Pointer is zero if list is not given */

  if ((long)RDB[DATA_PTR_GCU0] == 0)
    {
      /* Allocate memory */

      gcu = NewItem(DATA_PTR_GCU0, GCU_BLOCK_SIZE);

      /* Pointer to root universe */

      ptr = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
      CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);

      /* Set name */

      WDB[gcu + GCU_PTR_UNIV] = RDB[ptr + UNIVERSE_PTR_NAME];
    }

  /* Loop over gcu structures */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Find universe */

      uni = (long)RDB[DATA_PTR_U0];
      if ((uni = SeekListStr(uni, UNIVERSE_PTR_NAME,
                             GetText(gcu + GCU_PTR_UNIV))) < VALID_PTR)
        Error(0, "Universe %s in group constant generation does not exist",
              GetText(gcu + GCU_PTR_UNIV));

      /* Put pointers */

      WDB[gcu + GCU_PTR_UNIV] = (double)uni;
      WDB[uni + UNIVERSE_PTR_GCU] = (double)(gcu);

      /* Next */

      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Process discontinuity factos ****************************************/

  /* Loop over ADF's */

  adf = (long)RDB[DATA_PTR_ADF0];
  while (adf > VALID_PTR)
    {
      /* Stop tracks at outer boundary */

      WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Get pointer */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Compare */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, adf + ADF_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu > VALID_PTR)
        {
          /* Link ADF to universe */

          if ((long)RDB[gcu + GCU_PTR_ADF] > VALID_PTR)
            Error(0, "Universe %s is associated with multiple df's",
                  GetText(adf + ADF_PTR_GCU));
          else
            WDB[gcu + GCU_PTR_ADF] = (double)adf;

          /* Link universe to ADF */

          WDB[adf + ADF_PTR_GCU] = (double)gcu;
        }
      else
        Error(0, "Universe %s needed for df calculation does not exist",
              GetText(adf + ADF_PTR_GCU));

      /* Find surface */

      surf = (long)RDB[DATA_PTR_S0];
      if ((surf = SeekListStr(surf, SURFACE_PTR_NAME,
                              GetText(adf + ADF_PTR_SURF))) < VALID_PTR)
        Error(0, "Surface %s needed for df calculation does not exist",
              GetText(adf + ADF_PTR_SURF));

      /* Put pointer */

      WDB[adf + ADF_PTR_SURF] = (double)surf;

       /* Get pointer to surface parameters */

      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
      CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);

      /* Get number of faces and corners */

      switch ((long)RDB[surf + SURFACE_TYPE])
        {
        case SURF_PX:
        case SURF_PY:
        case SURF_PZ:
        case SURF_PLANE:
          {
            /* Planes */

            n = 1;
            m = 0;

            /* Set volume */

            WDB[adf + ADF_VOL] = 1.0;

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 1);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0] = 1.0;

            break;
          }
        case SURF_SQC:
          {
            /* Square prism */

            n = 4;
            m = 4;

            /* Set volume */

            WDB[adf + ADF_VOL] = 4.0*RDB[ptr + 2]*RDB[ptr + 2];

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*RDB[ptr + 2];
            WDB[loc0++] = 2.0*RDB[ptr + 2];
            WDB[loc0++] = 2.0*RDB[ptr + 2];
            WDB[loc0++] = 2.0*RDB[ptr + 2];

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*ADF_MID_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*ADF_MID_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*ADF_MID_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*ADF_MID_WIDTH*RDB[ptr + 2];

            /* Set corner areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*ADF_CORN_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*ADF_CORN_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*ADF_CORN_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*ADF_CORN_WIDTH*RDB[ptr + 2];

            break;
          }
        case SURF_TRIAG:
          {
            /* Check type */

            if (((long)RDB[surf + SURFACE_N_PARAMS] > 3) &&
                ((long)RDB[ptr + 3] != 1))
              Error(adf, "Triangle orientation must be 1 for ADF calculation");

            /* Triangular prism */

            if ((long)RDB[surf + SURFACE_N_PARAMS] == 6)
              {
                /* Truncated 3D prism */

                n = 5;
                m = 3;
              }
            else
              {
                /* Infinite 2D prism */

                n = 3;
                m = 3;
              }

            /* Set volume */

            if (n == 3)
              WDB[adf + ADF_VOL] = 3.0*SQRT3*RDB[ptr + 2]*RDB[ptr + 2];
            else
              WDB[adf + ADF_VOL] = 3.0*SQRT3*RDB[ptr + 2]*RDB[ptr + 2]*
                (RDB[ptr + 5] - RDB[ptr + 4]);

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, n);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*SQRT3*RDB[ptr + 2];
            WDB[loc0++] = 2.0*SQRT3*RDB[ptr + 2];
            WDB[loc0++] = 2.0*SQRT3*RDB[ptr + 2];

            if (n == 5)
              {
                WDB[loc0++] = 3.0*SQRT3*RDB[ptr + 2]*RDB[ptr + 2];
                WDB[loc0++] = 3.0*SQRT3*RDB[ptr + 2]*RDB[ptr + 2];
              }

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, n);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*SQRT3*ADF_MID_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*SQRT3*ADF_MID_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*SQRT3*ADF_MID_WIDTH*RDB[ptr + 2];

            /* Set corner areas */

            loc0 = ReallocMem(DATA_ARRAY, m);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*SQRT3*ADF_CORN_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*SQRT3*ADF_CORN_WIDTH*RDB[ptr + 2];
            WDB[loc0++] = 2.0*SQRT3*ADF_CORN_WIDTH*RDB[ptr + 2];

            break;
          }
        case SURF_RECT:
          {
            /* Rectangular prism */

            n = 4;
            m = 4;

            /* Set volume */

            WDB[adf + ADF_VOL] = (RDB[ptr + 1] - RDB[ptr])*
              (RDB[ptr + 3] - RDB[ptr + 2]);

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = (RDB[ptr + 3] - RDB[ptr + 2]);
            WDB[loc0++] = (RDB[ptr + 1] - RDB[ptr]);
            WDB[loc0++] = (RDB[ptr + 3] - RDB[ptr + 2]);
            WDB[loc0++] = (RDB[ptr + 1] - RDB[ptr]);

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2]);
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 1] - RDB[ptr]);
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2]);
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 1] - RDB[ptr]);

            /* Set corner areas (pick smallest width) */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;

            if ((RDB[ptr + 3] - RDB[ptr + 2]) < (RDB[ptr + 1] - RDB[ptr]))
              {
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2]);
              }
            else
              {
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr]);
              }

            break;
          }
        case SURF_HEXXC:
        case SURF_HEXYC:
          {
            /* Hexagonal prisms */

            n = 6;
            m = 6;

            /* Set volume */

            WDB[adf + ADF_VOL] = 6.0*TAN30*RDB[ptr + 2]*RDB[ptr + 2];

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2];
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2];
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2];
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2];
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2];
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2];

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3;

            /* Set corner areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3;

            break;
          }
        case SURF_CUBE:
          {
            /* Cube */

            n = 6;
            m = 4;

            /* Set volume */

            WDB[adf + ADF_VOL] = 8.0*RDB[ptr + 3]*RDB[ptr + 3]*RDB[ptr + 3];

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = 4.0*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*RDB[ptr + 3]*RDB[ptr + 3];

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = 4.0*ADF_MID_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*ADF_MID_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*ADF_MID_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*ADF_MID_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];

            /* Set corner areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;
            WDB[loc0++] = 4.0*ADF_CORN_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*ADF_CORN_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*ADF_CORN_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];
            WDB[loc0++] = 4.0*ADF_CORN_WIDTH*RDB[ptr + 3]*RDB[ptr + 3];

            break;
          }
        case SURF_CUBOID:
          {
            /* Cuboid */

            n = 6;
            m = 4;

            /* Set volume */

            WDB[adf + ADF_VOL] = (RDB[ptr + 1] - RDB[ptr])*
              (RDB[ptr + 3] - RDB[ptr + 2])*(RDB[ptr + 5] - RDB[ptr + 4]);

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = (RDB[ptr + 3] - RDB[ptr + 2])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = (RDB[ptr + 1] - RDB[ptr])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = (RDB[ptr + 3] - RDB[ptr + 2])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = (RDB[ptr + 1] - RDB[ptr])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = (RDB[ptr + 1] - RDB[ptr])
              *(RDB[ptr + 3] - RDB[ptr + 2]);
            WDB[loc0++] = (RDB[ptr + 1] - RDB[ptr])
              *(RDB[ptr + 3] - RDB[ptr + 2]);

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 1] - RDB[ptr])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2])
              *(RDB[ptr + 5] - RDB[ptr + 4]);
            WDB[loc0++] = ADF_MID_WIDTH*(RDB[ptr + 1] - RDB[ptr])
              *(RDB[ptr + 5] - RDB[ptr + 4]);

            /* Set corner areas (pick smallest width) */

            loc0 = ReallocMem(DATA_ARRAY, 4);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;

            if ((RDB[ptr + 3] - RDB[ptr + 2]) < (RDB[ptr + 1] - RDB[ptr]))
              {
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 3] - RDB[ptr + 2])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
              }
            else
              {
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
                WDB[loc0++] = ADF_CORN_WIDTH*(RDB[ptr + 1] - RDB[ptr])
                  *(RDB[ptr + 5] - RDB[ptr + 4]);
              }

            break;
          }
        case SURF_HEXXPRISM:
        case SURF_HEXYPRISM:
          {
            /* Hexagonal prisms */

            n = 8;
            m = 6;

            /* Set volume */

            WDB[adf + ADF_VOL] = 6.0*TAN30*RDB[ptr + 2]*RDB[ptr + 2]
              *(RDB[ptr + 4] - RDB[ptr + 3]);

            /* Set surface areas */

            loc0 = ReallocMem(DATA_ARRAY, 8);
            WDB[adf + ADF_PTR_SURF_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2]*(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2]*(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2]*(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2]*(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2]*(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*TAN30*RDB[ptr + 2]*(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 6.0*TAN30*RDB[ptr + 2]*RDB[ptr + 2];
            WDB[loc0++] = 6.0*TAN30*RDB[ptr + 2]*RDB[ptr + 2];

            /* Set mid areas */

            loc0 = ReallocMem(DATA_ARRAY, 8);
            WDB[adf + ADF_PTR_MID_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_MID_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);

            /* Set corner areas */

            loc0 = ReallocMem(DATA_ARRAY, 6);
            WDB[adf + ADF_PTR_CORN_AREA] = (double)loc0;
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);
            WDB[loc0++] = 2.0*RDB[ptr + 2]*ADF_CORN_WIDTH/SQRT3
              *(RDB[ptr + 4] - RDB[ptr + 3]);

            break;
          }
        default:
          {
            Error(0, "Surface %s is wrong type for df calculation",
                  GetText(surf + SURFACE_PTR_NAME));
          }
        }

      /* Store values */

      WDB[adf + ADF_NSURF] = (double)n;
      WDB[adf + ADF_NCORN] = (double)m;

      /* Next */

      adf = NextItem(adf);
    }

  /***************************************************************************/

  /***** Process pin-power distributions *************************************/

  /* Loop over PPW's */

  ppw = (long)RDB[DATA_PTR_PPW0];
  while (ppw > VALID_PTR)
    {
      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Get pointer */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Compare */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, ppw + PPW_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu > VALID_PTR)
        {
          /* Link PPW to universe */

          if ((long)RDB[gcu + GCU_PTR_PPW] > VALID_PTR)
            Error(0, "Universe %s is associated with multiple pin-power distributions",
                  GetText(ppw + PPW_PTR_GCU));
          else
            WDB[gcu + GCU_PTR_PPW] = (double)ppw;

          /* Link universe to PPW */

          WDB[ppw + PPW_PTR_GCU] = (double)gcu;
        }
      else
        Error(0, "Universe %s needed for pin-power distribution does not exist",
              GetText(ppw + PPW_PTR_GCU));

      /* Find lattice */

      lat = (long)RDB[DATA_PTR_L0];
      if ((lat = SeekListStr(lat, LAT_PTR_NAME,
                              GetText(ppw + PPW_PTR_LAT))) < VALID_PTR)
        Error(0, "Lattice %s needed for pin-power distribution does not exist",
              GetText(ppw + PPW_PTR_LAT));

      /* Put pointer */

      WDB[ppw + PPW_PTR_LAT] = (double)lat;

      /* Put type */

      WDB[ppw + PPW_LAT_TYPE] = RDB[lat + LAT_TYPE];

      /* Get number of pins */

      if ((n = (long)RDB[lat + LAT_NTOT]) > 0)
        WDB[ppw + PPW_NP] = (double)n;
      else
        Die(FUNCTION_NAME, "Number of pins is zero");

      /* Next distribution */

      ppw = NextItem(ppw);
    }

  /***************************************************************************/

  /***** Process albedos *****************************************************/

  /* Loop over albedos */

  alb = (long)RDB[DATA_PTR_ALB0];
  while (alb > VALID_PTR)
    {
      /* Stop tracks at outer boundary */

      WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Get pointer */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Compare */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, alb + ALB_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu > VALID_PTR)
        {
          /* Link ALB to universe */

          if ((long)RDB[gcu + GCU_PTR_ALB] > VALID_PTR)
            Error(0, "Universe %s is associated with multiple albedos",
                  GetText(alb + ALB_PTR_GCU));
          else
            WDB[gcu + GCU_PTR_ALB] = (double)alb;

          /* Link universe to ALB */

          WDB[alb + ALB_PTR_GCU] = (double)gcu;
        }
      else
        Error(0, "Universe %s needed for albedo calculation does not exist",
              GetText(alb + ALB_PTR_GCU));

      /* Find surface */

      surf = (long)RDB[DATA_PTR_S0];
      if ((surf = SeekListStr(surf, SURFACE_PTR_NAME,
                              GetText(alb + ALB_PTR_SURF))) < VALID_PTR)
        Error(0, "Surface %s needed for albedo calculation does not exist",
              GetText(alb + ALB_PTR_SURF));

      /* Put pointer */

      WDB[alb + ALB_PTR_SURF] = (double)surf;

      /* Get number of faces */

      switch ((long)RDB[surf + SURFACE_TYPE])
        {
        case SURF_PX:
        case SURF_PY:
        case SURF_PZ:
        case SURF_PLANE:
          {
            /* Planes */

            n = 1;

            break;
          }
        case SURF_SQC:
          {
            /* Square prism */

            n = 4;

            break;
          }
        case SURF_HEXXC:
        case SURF_HEXYC:
          {
            /* Hexagonal prisms */

            n = 6;

            break;
          }
        case SURF_CUBE:
        case SURF_CUBOID:
          {
            /* Cube and cuboid */

            n = 6;

            break;
          }
        case SURF_HEXXPRISM:
        case SURF_HEXYPRISM:
          {
            /* Hexagonal prisms */

            n = 8;

            break;
          }
        default:
          {
            Error(0, "Surface %s is wrong type for albedo calculation",
                  GetText(surf + SURFACE_PTR_NAME));
          }
        }

      /* Store values */

      WDB[alb + ALB_NSURF] = (double)n;

      /* Next */

      alb = NextItem(alb);
    }

  /***************************************************************************/

  /***** Process micro depletion stuff ***************************************/

  /* Loop over definitions */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Get pointer */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Compare */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, loc0 + MDEP_PTR_GCU))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check pointer */

      if (gcu > VALID_PTR)
        {
          /* Link to universe */

          if ((long)RDB[gcu + GCU_PTR_MDEP] > VALID_PTR)
            Error(0, "Universe %s is associated with multiple micro depletion definitions",
                  GetText(loc0 + MDEP_PTR_GCU));
          else
            WDB[gcu + GCU_PTR_MDEP] = (double)loc0;

          /* Link universe to MDEP */

          WDB[loc0 + MDEP_PTR_GCU] = (double)gcu;
        }
      else
        Error(0, "Universe %s needed for pin-power distribution does not exist",
              GetText(loc0 + MDEP_PTR_GCU));

      /* Link materials */

      if ((loc1 = (long)RDB[loc0 + MDEP_PTR_MAT]) > VALID_PTR)
        {
          /* Reset pointer */

          WDB[loc0 + MDEP_PTR_MAT] = NULLPTR;

          /* Loop over materials */

          while ((long)RDB[loc1] > VALID_PTR)
            {
              /* Reset count */

              n = 0;

              /* Loop over materials and find match */

              mat = (long)RDB[DATA_PTR_M0];
              while (mat > VALID_PTR)
                {
                  /* Skip parent materials */

                  if ((long)RDB[mat + MATERIAL_DIV_TYPE] ==
                      MAT_DIV_TYPE_PARENT)
                    {
                      /* Pointer to next */

                      mat = NextItem(mat);

                      /* Cycle loop */

                      continue;
                    }

                  /* Check pointer to parent */

                  if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
                      < VALID_PTR)
                    mat0 = mat;

                  /* Compare */

                  if (CompareStr(mat0 + MATERIAL_PTR_NAME, loc1))
                    {
                      /* Allocate memory */

                      ptr = NewItem(loc0 + MDEP_PTR_MAT,
                                    MDEP_MAT_BLOCK_SIZE);

                      /* Put pointer */

                      WDB[ptr + MDEP_MAT_PTR_MAT] = (double)mat;

                      /* Add count */

                      n++;
                    }

                  /* Next material */

                  mat = NextItem(mat);
                }

              /* Check pointer */

              if (n == 0)
                Error(loc0, "Material %s is not defined", GetText(loc1));

              /* Next material in list */

              loc1++;
            }
        }
      else
        {
          /* Include all burnable materials */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Check burnable flag and divisor */

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
                  {
                    /* Allocate memory */

                    ptr = NewItem(loc0 + MDEP_PTR_MAT, MDEP_MAT_BLOCK_SIZE);

                    /* Put pointer */

                    WDB[ptr + MDEP_MAT_PTR_MAT] = (double)mat;
                  }

              /* Next material */

              mat = NextItem(mat);
            }
        }

      /* Check that list was created */

      if ((loc1 = (long)RDB[loc0 + MDEP_PTR_MAT]) < VALID_PTR)
        Error(loc0, "No burnable materials");

      /* Next definition */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Allocate memory for statistics **************************************/

  /* NOTE: noiden NewStat() -kutsujen pitää tulla peräkkäin koska */
  /* CoefOutput() luuppaa ensimmäisestä viimeiseen. */

  /* Get number of universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
  n = ListSize(gcu);

  /* Preallocate memory from scoring buffer to speed up calculation */
  /* nfg*nfg ~ 34, nfg ~ 156 (approximate values) */

  PreallocMem(n*(34*34 + 156)*nfg*BUF_BLOCK_SIZE, BUF_ARRAY);

  /* Loop over gcu structures */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /***********************************************************************/

      /***** Micro-group data ************************************************/

      /* Allocate memory for micro-group data */

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_FLX] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_FISS_FLX] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_TOT] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_ABS] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_FISS] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_CHIT] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_CHIP] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_CHID] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_NSF] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_FISSE] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_INV_V] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT0] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP0] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT1] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP1] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT2] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP2] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT3] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP3] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT4] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP4] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT5] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP5] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT6] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP6] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATT7] = (double)ptr;

      ptr = AllocPrivateData(nmg*nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_SCATTP7] = (double)ptr;

      /* Check poison calculation */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
        {
          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_I135_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_XE135_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_XE135M_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM147_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM148_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM148M_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM149_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_SM149_YIELD] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_I135_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_XE135_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_XE135M_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM147_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM148_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM148M_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_PM149_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_SM149_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_XE135_MACRO_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_XE135M_MACRO_ABS] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_SM149_MACRO_ABS] = (double)ptr;
        }

      /* B1 flux spectrum and diffusion coefficient */

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_B1_FLX] = (double)ptr;

      ptr = AllocPrivateData(nmg, RES2_ARRAY);
      WDB[gcu + GCU_MICRO_B1_DIFFCOEF] = (double)ptr;

      /* Critical spectrum burnup correction */

      ptr = ReallocMem(DATA_ARRAY, nmg);
      WDB[gcu + GCU_PTR_BURNUP_SPEC_CORR] = (double)ptr;

      /* Eddington factors (Andrew Hall 2016/3/31) */

      if ((long)RDB[DATA_OPTI_EDDINGTON_CALC] == YES)
        {
          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_FLX_XX] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_FLX_XY] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_FLX_XZ] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_FLX_YY] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_FLX_YZ] = (double)ptr;

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_FLX_ZZ] = (double)ptr;
        }

      /***********************************************************************/

      /***** Few-group constants in infinite spectrum ************************/

      ptr = NewStat("INF_MICRO_FLX", 1, nmg);
      WDB[gcu + GCU_INF_MICRO_FLX] = (double)ptr;

      ptr = NewStat("INF_FLX", 1, nfg - 1);
      WDB[gcu + GCU_INF_FLX] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      /* Remember first value for coefoutput.c */

      WDB[gcu + GCU_PTR_FIRST_STAT] = (double)ptr;

      ptr = NewStat("INF_FISS_FLX", 1, nfg - 1);
      WDB[gcu + GCU_INF_FISS_FLX] = (double)ptr;

      ptr = NewStat("INF_KINF", 1, 1);
      WDB[gcu + GCU_INF_KINF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_TOT", 1, nfg - 1);
      WDB[gcu + GCU_INF_TOT] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_CAPT", 1, nfg - 1);
      WDB[gcu + GCU_INF_CAPT] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_FISS", 1, nfg - 1);
      WDB[gcu + GCU_INF_FISS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_NSF", 1, nfg - 1);
      WDB[gcu + GCU_INF_NSF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_KAPPA", 1, nfg - 1);
      WDB[gcu + GCU_INF_KAPPA] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_INVV", 1, nfg - 1);
      WDB[gcu + GCU_INF_INVV] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_NUBAR", 1, nfg - 1);
      WDB[gcu + GCU_INF_NUBAR] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_REMXS", 1, nfg - 1);
      WDB[gcu + GCU_INF_REMXS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_RABSXS", 1, nfg - 1);
      WDB[gcu + GCU_INF_RABSXS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_CHIT", 1, nfg - 1);
      WDB[gcu + GCU_INF_CHIT] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_CHIP", 1, nfg - 1);
      WDB[gcu + GCU_INF_CHIP] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_CHID", 1, nfg - 1);
      WDB[gcu + GCU_INF_CHID] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_I135_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_I135_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_XE135_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_XE135_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_XE135M_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_XE135M_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM147_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM147_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM148_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM148_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM148M_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM148M_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM149_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM149_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SM149_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_INF_SM149_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_I135_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_I135_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_XE135_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_XE135_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_XE135M_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_XE135M_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM147_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM147_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM148_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM148_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM148M_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM148M_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_PM149_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_PM149_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SM149_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_SM149_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_XE135_MACRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_XE135_MACRO_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_XE135M_MACRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_XE135M_MACRO_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SM149_MACRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_INF_SM149_MACRO_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S0", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S1", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S2", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S3", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S4", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S5", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S6", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_S7", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_S7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP0", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP1", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP2", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP3", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP4", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP5", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP6", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SP7", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_INF_SP7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT0", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT1", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT2", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT3", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT4", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT5", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT6", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATT7", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATT7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP0", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP1", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP2", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP3", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP4", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP5", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP6", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_SCATTP7", 1, nfg - 1);
      WDB[gcu + GCU_INF_SCATTP7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_TRANSPXS", 1, nfg - 1);
      WDB[gcu + GCU_INF_TRANSPXS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("INF_DIFFCOEF", 1, nfg - 1);
      WDB[gcu + GCU_INF_DIFFCOEF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      /* Eddington factors (Andrew Hall 2016/3/31) */

      if ((long)RDB[DATA_OPTI_EDDINGTON_CALC] == YES)
        {
          ptr = NewStat("INF_FLX_EDDI_XX", 1, nfg - 1);
          WDB[gcu + GCU_INF_FLX_XX] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          ptr = NewStat("INF_FLX_EDDI_XY", 1, nfg - 1);
          WDB[gcu + GCU_INF_FLX_XY] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          ptr = NewStat("INF_FLX_EDDI_XZ", 1, nfg - 1);
          WDB[gcu + GCU_INF_FLX_XZ] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          ptr = NewStat("INF_FLX_EDDI_YY", 1, nfg - 1);
          WDB[gcu + GCU_INF_FLX_YY] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          ptr = NewStat("INF_FLX_EDDI_YZ", 1, nfg - 1);
          WDB[gcu + GCU_INF_FLX_YZ] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          ptr = NewStat("INF_FLX_EDDI_ZZ", 1, nfg - 1);
          WDB[gcu + GCU_INF_FLX_ZZ] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);
        }

      /* Remember last value for coefoutput.c */

      WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

      /***********************************************************************/

      /***** Few-group constants in critical spectrum ************************/

      /* Spectrum */

      ptr = NewStat("B1_MICRO_FLX", 1, nmg);
      WDB[gcu + GCU_B1_MICRO_FLX] = (double)ptr;

      /* Constants */

      ptr = NewStat("B1_KINF", 1, 1);
      WDB[gcu + GCU_B1_KINF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_KEFF", 1, 1);
      WDB[gcu + GCU_B1_KEFF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_B2", 1, 1);
      WDB[gcu + GCU_B1_B2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_ERR", 1, 1);
      WDB[gcu + GCU_B1_ERR] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_FLX", 1, nfg - 1);
      WDB[gcu + GCU_B1_FLX] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_FISS_FLX", 1, nfg - 1);
      WDB[gcu + GCU_B1_FISS_FLX] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_TOT", 1, nfg - 1);
      WDB[gcu + GCU_B1_TOT] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_CAPT", 1, nfg - 1);
      WDB[gcu + GCU_B1_CAPT] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_FISS", 1, nfg - 1);
      WDB[gcu + GCU_B1_FISS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_NSF", 1, nfg - 1);
      WDB[gcu + GCU_B1_NSF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_KAPPA", 1, nfg - 1);
      WDB[gcu + GCU_B1_KAPPA] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_INVV", 1, nfg - 1);
      WDB[gcu + GCU_B1_INVV] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_NUBAR", 1, nfg - 1);
      WDB[gcu + GCU_B1_NUBAR] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_REMXS", 1, nfg - 1);
      WDB[gcu + GCU_B1_REMXS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_RABSXS", 1, nfg - 1);
      WDB[gcu + GCU_B1_RABSXS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_CHIT", 1, nfg - 1);
      WDB[gcu + GCU_B1_CHIT] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_CHIP", 1, nfg - 1);
      WDB[gcu + GCU_B1_CHIP] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_CHID", 1, nfg - 1);
      WDB[gcu + GCU_B1_CHID] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_I135_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_I135_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_XE135_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_XE135_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_XE135M_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_XE135M_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM147_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM147_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM148_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM148_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM148M_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM148M_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM149_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM149_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SM149_YIELD", 1, nfg - 1);
      WDB[gcu + GCU_B1_SM149_YIELD] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_I135_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_I135_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_XE135_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_XE135_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_XE135M_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_XE135M_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM147_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM147_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM148_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM148_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM148M_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM148M_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_PM149_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_PM149_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SM149_MICRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_SM149_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_XE135_MACRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_XE135_MACRO_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_XE135M_MACRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_XE135M_MACRO_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SM149_MACRO_ABS", 1, nfg - 1);
      WDB[gcu + GCU_B1_SM149_MACRO_ABS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S0", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S1", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S2", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S3", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S4", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S5", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S6", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_S7", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_S7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP0", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP1", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP2", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP3", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP4", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP5", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP6", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SP7", 2, nfg - 1, nfg - 1);
      WDB[gcu + GCU_B1_SP7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT0", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT1", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT2", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT3", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT4", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT5", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT6", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATT7", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATT7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP0", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP0] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP1", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP1] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP2", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP2] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP3", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP3] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP4", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP4] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP5", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP5] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP6", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP6] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_SCATTP7", 1, nfg - 1);
      WDB[gcu + GCU_B1_SCATTP7] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_TRANSPXS", 1, nfg - 1);
      WDB[gcu + GCU_B1_TRANSPXS] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      ptr = NewStat("B1_DIFFCOEF", 1, nfg - 1);
      WDB[gcu + GCU_B1_DIFFCOEF] = (double)ptr;
      if((long)RDB[DATA_GC_STAT_TESTS] == YES)
        AllocStatHistory(ptr);

      /* Remember last value for coefoutput.c */

      WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

      /***********************************************************************/

      /***** Time constants **************************************************/

      ptr = NewStat("BETA_EFF", 1, 9);
      WDB[gcu + GCU_MEULEKAMP_BETA_EFF] = (double)ptr;

      ptr = NewStat("LAMBDA", 1, 9);
      WDB[gcu + GCU_MEULEKAMP_LAMBDA] = (double)ptr;

      ptr = NewStat("TOT_FISS", 1, 1);
      WDB[gcu + GCU_MEULEKAMP_TOT_FISS] = (double)ptr;

      /* Remember last value for coefoutput.c */

      WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

      /***********************************************************************/

      /***** ADF's ***********************************************************/

      /* Get pointer to adf's */

      if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) > VALID_PTR)
        {
          /* Get number of surfaces and corners */

          n = (long)RDB[adf + ADF_NSURF];
          m = (long)RDB[adf + ADF_NCORN];

          /* Allocate memory for micro-group fluxes and currents */

          if (n > 0)
            {
              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_SURF_FLUX] = (double)ptr;

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_SURF_IN_CURR] = (double)ptr;

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_SURF_OUT_CURR] = (double)ptr;

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_MID_IN_CURR] = (double)ptr;

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_MID_OUT_CURR] = (double)ptr;

              /* Sign moments of discontinuity factors */

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_SGN_SURF_FLUX] = (double)ptr;

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_SGN_SURF_IN_CURR] = (double)ptr;

              ptr = AllocPrivateData(nmg*n, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_SGN_SURF_OUT_CURR] = (double)ptr;
            }

          if (m > 0)
            {
              ptr = AllocPrivateData(nmg*m, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_CORN_FLUX] = (double)ptr;

              ptr = AllocPrivateData(nmg*m, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_CORN_IN_CURR] = (double)ptr;

              ptr = AllocPrivateData(nmg*m, RES2_ARRAY);
              WDB[gcu + GCU_MICRO_ADF_CORN_OUT_CURR] = (double)ptr;
            }

          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_ADF_CELL_FLUX] = (double)ptr;

          /* Allocate memory for surface constants */

          if (n > 0)
            {
              ptr = NewStat("DF_HET_SURF_FLUX", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_HET_SURF_FLUX] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_HOM_SURF_FLUX", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_HOM_SURF_FLUX] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_SURF_DF", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SURF_DF] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_SURF_IN_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SURF_IN_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_SURF_OUT_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SURF_OUT_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_SURF_NET_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SURF_NET_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_MID_IN_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_MID_IN_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_MID_OUT_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_MID_OUT_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_MID_NET_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_MID_NET_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              /* Sign moments of discontinuity factors */

              ptr = NewStat("DF_SGN_HET_SURF_FLX", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SGN_HET_SURF_FLUX] = (double)ptr;
              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                  AllocStatHistory(ptr);

              ptr = NewStat("DF_SGN_HOM_SURF_FLUX", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SGN_HOM_SURF_FLUX] = (double)ptr;
              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                  AllocStatHistory(ptr);

              ptr = NewStat("DF_SGN_SURF_DF", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SGN_SURF_DF] = (double)ptr;
              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                  AllocStatHistory(ptr);

              ptr = NewStat("DF_SGN_SURF_IN_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SGN_SURF_IN_CURR] = (double)ptr;
              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                  AllocStatHistory(ptr);

              ptr = NewStat("DF_SGN_SURF_OUT_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SGN_SURF_OUT_CURR] = (double)ptr;
              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                  AllocStatHistory(ptr);

              ptr = NewStat("DF_SGN_SURF_NET_CURR", 2, n, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_SGN_SURF_NET_CURR] = (double)ptr;
              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                  AllocStatHistory(ptr);
            }

          /* Allocate memory for corner constants */

          if (m > 0)
            {
              ptr = NewStat("DF_HET_CORN_FLUX", 2, m, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_HET_CORN_FLUX] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_HOM_CORN_FLUX", 2, m, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_HOM_CORN_FLUX] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_CORN_DF", 2, m, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_CORN_DF] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_CORN_IN_CURR", 2, m, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_CORN_IN_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_CORN_OUT_CURR", 2, m, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_CORN_OUT_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);

              ptr = NewStat("DF_CORN_NET_CURR", 2, m, nfg - 1);
              WDB[gcu + GCU_RES_FG_DF_CORN_NET_CURR] = (double)ptr;
              if((long)RDB[DATA_GC_STAT_TESTS] == YES)
                AllocStatHistory(ptr);
            }

          /* Allocate memory for integral values */

          ptr = NewStat("DF_HET_VOL_FLUX", 1, nfg - 1);
          WDB[gcu + GCU_RES_FG_DF_HET_VOL_FLUX] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          ptr = NewStat("DF_HOM_VOL_FLUX", 1, nfg - 1);
          WDB[gcu + GCU_RES_FG_DF_HOM_VOL_FLUX] = (double)ptr;
          if((long)RDB[DATA_GC_STAT_TESTS] == YES)
            AllocStatHistory(ptr);

          /* Remember last value for coefoutput.c */

          WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;
        }

      /***********************************************************************/

      /***** Pin-power distribution ******************************************/

      if ((ppw = (long)RDB[gcu + GCU_PTR_PPW]) > VALID_PTR)
        if ((n = (long)RDB[ppw + PPW_NP]) > 0)
          {
            /* Allocate memory for macro-group distribution */

            ptr = NewStat("PPW_POW", 2, n, nfg - 1);
            WDB[gcu + GCU_RES_FG_PPW_POW] = (double)ptr;

            /* Allocate memory for macro-group flux */

            ptr = NewStat("PPW_HOM_FLUX", 2, n, nfg - 1);
            WDB[gcu + GCU_RES_FG_PPW_HOM_FLUX] = (double)ptr;

            /* Allocate memory for macro-group form factors */

            ptr = NewStat("PPW_FF", 2, n, nfg - 1);
            WDB[gcu + GCU_RES_FG_PPW_FF] = (double)ptr;

            /* Allocate memory for macro-group coordinates */
            /* (tosi pätevä termi). */

            ptr = NewStat("PPW_XYZ", 3, 3, n, nfg - 1);
            WDB[gcu + GCU_RES_FG_PPW_XYZ] = (double)ptr;

            /* Remember last value for coefoutput.c */

            WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

            /* Allocate memory for power distribution and coordinates */
            /* NOTE: tää varataan tuolta RES2 arraystä nfg-jaolla.    */
            /* yksi ylimääräinen ryhmäbini on kokonaisarvoa varten.   */

            ptr = AllocPrivateData(nfg*n, RES2_ARRAY);
            WDB[gcu + GCU_MICRO_PPW_POW] = (double)ptr;

            ptr = AllocPrivateData(3*nfg*n, RES2_ARRAY);
            WDB[gcu + GCU_MICRO_PPW_XYZ] = (double)ptr;
          }

      /***********************************************************************/

      /***** Micro-depletion *************************************************/

      if ((loc0 = (long)RDB[gcu + GCU_PTR_MDEP]) > VALID_PTR)
        {
          /* Get size */

          n = (long)RDB[loc0 + MDEP_N_REA];
          CheckValue(FUNCTION_NAME, "n", "", n, 1, 1000000);

          /* Allocate memory atomic densities */

          ptr = NewStat("MDEP_ADENS", 1, n);
          WDB[gcu + GCU_RES_FG_MICRO_DEP_ADENS] = (double)ptr;

          /* Allocate memory for xs */

          ptr = NewStat("MDEP_XS", 2, n, nfg - 1);
          WDB[gcu + GCU_RES_FG_MICRO_DEP_XS] = (double)ptr;

          /* Remember last value for coefoutput.c */

          WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

          /* NOTE: tää varataan tuolta RES2 arraystä nfg-jaolla. */

          ptr = AllocPrivateData(nfg*n, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_MICRO_DEP_XS] = (double)ptr;
        }

      /***********************************************************************/

      /***** CMM diffusion coefficients **************************************/

      /* Allocate memory for scoring cumulative square distances */

      ptr = AllocPrivateData(5*nmg, RES2_ARRAY);
      WDB[gcu + GCU_BUF_CMM_CUMU_R2] = (double)ptr;

      ptr = NewStat("CMM_CUMU_R2", 2, nfg, 5);
      WDB[gcu + GCU_CMM_CUMU_R2] = (double)ptr;

      /* Allocate memory for stats */

      ptr = NewStat("CMM_MIGR_AREA", 2, nfg - 1, 4);
      WDB[gcu + GCU_CMM_MIGR_AREA] = (double)ptr;

      ptr = NewStat("CMM_REMXS", 1, nfg - 1);
      WDB[gcu + GCU_CMM_REMXS] = (double)ptr;

      ptr = NewStat("CMM_DIFFCOEF", 1, nfg - 1);
      WDB[gcu + GCU_CMM_DIFFCOEF] = (double)ptr;

      ptr = NewStat("CMM_TRANSPXS", 1, nfg - 1);
      WDB[gcu + GCU_CMM_TRANSPXS] = (double)ptr;

      ptr = NewStat("CMM_DIFFCOEF_X", 1, nfg - 1);
      WDB[gcu + GCU_CMM_DIFFCOEF_X] = (double)ptr;

      ptr = NewStat("CMM_TRANSPXS_X", 1, nfg - 1);
      WDB[gcu + GCU_CMM_TRANSPXS_X] = (double)ptr;

      ptr = NewStat("CMM_DIFFCOEF_Y", 1, nfg - 1);
      WDB[gcu + GCU_CMM_DIFFCOEF_Y] = (double)ptr;

      ptr = NewStat("CMM_TRANSPXS_Y", 1, nfg - 1);
      WDB[gcu + GCU_CMM_TRANSPXS_Y] = (double)ptr;

      ptr = NewStat("CMM_DIFFCOEF_Z", 1, nfg - 1);
      WDB[gcu + GCU_CMM_DIFFCOEF_Z] = (double)ptr;

      ptr = NewStat("CMM_TRANSPXS_Z", 1, nfg - 1);
      WDB[gcu + GCU_CMM_TRANSPXS_Z] = (double)ptr;

      /* Remember last value for coefoutput.c */

      WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

      /***********************************************************************/

      /***** Albedos *********************************************************/

      if ((alb = (long)RDB[gcu + GCU_PTR_ALB]) > VALID_PTR)
        if ((n = (long)RDB[alb + ALB_NSURF]) > 0)
          {
            /* Allocate memory for currents */

            ptr = NewStat("ALB_IN_CURR", 2, n, nfg - 1);
            WDB[gcu + GCU_RES_FG_ALB_IN_CURR] = (double)ptr;

            ptr = NewStat("ALB_OUT_CURR", 4, n, (nfg - 1), n, (nfg - 1));
            WDB[gcu + GCU_RES_FG_ALB_OUT_CURR] = (double)ptr;

            /* Allocate memory for total and partial albedos */

            ptr = NewStat("ALB_TOT_ALB", 2, nfg - 1, nfg - 1);
            WDB[gcu + GCU_RES_FG_TOT_ALB] = (double)ptr;

            ptr = NewStat("ALB_PART_ALB", 4, n, (nfg - 1), n, (nfg - 1));
            WDB[gcu + GCU_RES_FG_PART_ALB] = (double)ptr;

            /* Remember last value for coefoutput.c */

            WDB[gcu + GCU_PTR_LAST_STAT] = (double)ptr;

            /* Allocate memory for scoring albedo currents.        */
            /* NOTE: tää varataan tuolta RES2 arraystä nfg-jaolla. */

            ptr = AllocPrivateData((nfg - 1)*n, RES2_ARRAY);
            WDB[gcu + GCU_MICRO_ALB_IN_CURR] = (double)ptr;

            ptr = AllocPrivateData((nfg - 1)*(nfg - 1)*n*n, RES2_ARRAY);
            WDB[gcu + GCU_MICRO_ALB_OUT_CURR] = (double)ptr;
          }

      /***********************************************************************/

      /* Next */

      gcu = NextItem(gcu);
    }

  /* Allocate memory for universe pointer */

  AllocValuePair(DATA_GCU_PTR_UNI);

  /***************************************************************************/

  /***** MORA stuff **********************************************************/

  loc0 = (long)RDB[DATA_PTR_MORA0];
  while (loc0 > VALID_PTR)
    {
      /* Avoid compiler warning */

      uni = -1;

      /* Loop over gcu list */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Pointer to universe */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Compare names */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, loc0 + MORA_PTR_UNIV))
            break;

          /* Next */

          gcu = NextItem(gcu);
        }

      /* Check */

      if (gcu < VALID_PTR)
        Error(0, "Universe %s is not included in gc calculation",
              GetText(loc0 + MORA_PTR_UNIV));

      /* Put pointers */

      WDB[loc0 + MORA_PTR_UNIV] = (double)uni;
      WDB[gcu + GCU_PTR_MORA] = (double)loc0;

      /* Find energy group structure */

      loc1 = (long)RDB[DATA_PTR_ENE0];
      while (loc1 > VALID_PTR)
        {
          /* Compare names */

          if (CompareStr(loc1 + ENE_PTR_NAME, loc0 + MORA_PTR_EG))
            break;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Check */

      if (loc1 < VALID_PTR)
        Error(0, "Energy group structure %s is not defined",
              GetText(loc0 + MORA_PTR_EG));

      /* Put pointer */

      WDB[loc0 + MORA_PTR_EG] = (double)RDB[loc1 + ENE_PTR_GRID];

      /* Get number of energy groups and cosine bins */

      n = (long)RDB[loc1 + ENE_NB];
      m = (long)RDB[loc0 + MORA_N_COS];

      /* Put number of energy groups */

      WDB[loc0 + MORA_N_EG] = (double)n;

      /* Allocate memory for stats */

      ptr = NewStat("TOT", 1, n);
      WDB[loc0 + MORA_PTR_TOT] = (double)ptr;

      ptr = NewStat("CAPT", 1, n);
      WDB[loc0 + MORA_PTR_CAPT] = (double)ptr;

      ptr = NewStat("FISS", 1, n);
      WDB[loc0 + MORA_PTR_FISS] = (double)ptr;

      ptr = NewStat("PNU", 1, n);
      WDB[loc0 + MORA_PTR_PNU] = (double)ptr;

      ptr = NewStat("DNU", 1, n);
      WDB[loc0 + MORA_PTR_DNU] = (double)ptr;

      ptr = NewStat("KAPPA", 1, n);
      WDB[loc0 + MORA_PTR_KAPPA] = (double)ptr;

      ptr = NewStat("CHIP", 1, n + 1);
      WDB[loc0 + MORA_PTR_CHIP] = (double)ptr;

      ptr = NewStat("CHID", 1, n + 1);
      WDB[loc0 + MORA_PTR_CHID] = (double)ptr;

      ptr = NewStat("FLX", 1, n + 1);
      WDB[loc0 + MORA_PTR_FLX] = (double)ptr;

      ptr = NewStat("SCATTP", 3, n, n, m);
      WDB[loc0 + MORA_PTR_SCATTP] = (double)ptr;

      ptr = NewStat("SCATTW", 3, n, n, m);
      WDB[loc0 + MORA_PTR_SCATTW] = (double)ptr;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check nested universes **********************************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Get pointer to universe */

      uni = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Loop over remaining */

      loc0 = NextItem(gcu);
      while (loc0 > VALID_PTR)
        {
          /* Get pointer to universe */

          ptr = (long)RDB[loc0 + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Compare levels */

          if ((long)RDB[uni + UNIVERSE_LEVEL] !=
              (long)RDB[ptr + UNIVERSE_LEVEL])
            WDB[DATA_MULTI_LEVEL_GCU] = (double)YES;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Next */

      gcu = NextItem(gcu);
    }

  /* Print */

  fprintf(outp, " - %ld energy groups in micro-group structure\n", nmg);
  fprintf(outp, " - %ld energy groups in macro-group structure\n", nfg - 1);

  if ((long)RDB[DATA_B1_CALC] == YES)
    fprintf(outp, " - B1 fundamental mode calculation is run\n");
  else
    fprintf(outp, " - B1 fundamental mode calculation is not run\n");

  ptr = (long)RDB[DATA_PTR_GCU0];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  fprintf(outp, " - Group constants generated in %ld universes\n",
          ListSize(ptr));

  if ((ptr = (long)RDB[DATA_PTR_ADF0]) > VALID_PTR)
    fprintf(outp, " - Discontinuity factors calculated in %ld universes\n",
            ListSize(ptr));
  else
    fprintf(outp, " - Discontinuity factors are not calculated\n");

  if ((ptr = (long)RDB[DATA_PTR_PPW0]) > VALID_PTR)
    fprintf(outp, " - Pin-power distributions are calculated in %ld universes\n",
            ListSize(ptr));
  else
    fprintf(outp, " - Pin-power distributions are not calculated\n");

  if ((ptr = (long)RDB[DATA_PTR_ALB0]) > VALID_PTR)
    fprintf(outp, " - Albedos are calculated in %ld universes\n",
            ListSize(ptr));
  else
    fprintf(outp, " - Albedos are not calculated\n");

  if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
    fprintf(outp, " - Poison cross sections are calculated\n");
  else
    fprintf(outp, " - Poison cross sections are not calculated\n");

  fprintf(outp, "\n");

  /* For EDo (18.1.2014) */

  DiffCoefED(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  /***************************************************************************/
}

/*****************************************************************************/
