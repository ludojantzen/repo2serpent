/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatedtmajorants.c                         */
/*                                                                           */
/* Created:       2011/11/03 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates neutron and photon majorants for delta-tracking   */
/*                                                                           */
/* Comments: - Tässä on rutiinit otettu vanhasta calculatemajorant.c:stä     */
/*                                                                           */
/*           - Noiden makroskooppisten vaikutusalojen rekonstruoinnin        */
/*             voisi laittaa omalle aliohjelmalleen, jos niitä käytetään     */
/*             muuallakin kuin täällä. (typerä idea, koska esim.             */
/*             materialtotals.c:ssä vaikutusaloja kerrotaan ennen            */
/*             summaamista)                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateDTMajorants:"

/* Use local function to simplify OpenMP implementation */

void CalculateDTMajorants0(long, double *, long);

/*****************************************************************************/

void CalculateDTMajorants()
{
  long erg, ne, loc0, loc1, loc2, loc3, mat, rea, ptr, n, i, i0, sz, nt, rls, nuc;
  double *tot, *xs, **maj, adens, Emin, Emax;

  /* Check DT flag */

  if ((long)RDB[DATA_OPT_USE_DT] == NO)
    return;

  /***************************************************************************/

  /***** Multi-group mode ****************************************************/

  if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
    {
      fprintf(outp, "Calculating DT majorant cross sections...\n");

      /* Get pointer to energy grid */

      erg = (long)RDB[DATA_COARSE_MG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Number of points */

      ne = (long)RDB[erg + ENERGY_GRID_NE];

      /* Pointer to energy array */

      loc0 = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Get pointer to majorant reaction data */

      loc1 = (long)RDB[DATA_PTR_MAJORANT];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get pointer to data */

      loc1 = (long)RDB[loc1 + REACTION_PTR_MGXS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Reset data */

      memset(&WDB[loc1], 0.0, ne*sizeof(double));

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if domain decomposition is in use */

          if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Check if material is included in majorant */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT)
            {
              /* Pointer to total xs */

              rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Pointer to data */

              loc2 = (long)RDB[rea + REACTION_PTR_MGXS];
              CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

              /* Loop over data and find maxima */

              for (n = 0; n < ne; n++)
                if (RDB[loc2 + n] > RDB[loc1 + n])
                  WDB[loc1 + n] = RDB[loc2 + n];
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Find majorant minimum */

      WDB[DATA_MIN_NMACROXS] = INFTY;

      for (n = 0; n < ne; n++)
        if ((RDB[loc1 + n] > 0.0) && (RDB[loc1 + n] < RDB[DATA_MIN_NMACROXS]))
          WDB[DATA_MIN_NMACROXS] = RDB[loc1 + n];

      /* Tää lasketaan eri tavalla kuin Serpent 1:ssä (sielläkin väärin) */
      /* etsi tolle joku järkevä arvo ja korvaa koko hässäkkä sillä */
      /*
      WDB[DATA_MIN_NMACROXS] = 1E-2;
      */
      /* NOTE: Tää ei nyt toimi gammatransportlaskujen kanssa kun niille */
      /* ei lasketa majoranttia. */

      fprintf(outp, "OK.\n\n");

      /* Combine majorants with other domains */

      DistributeDDMajorants();

      return;
    }

  /* Continuous-energy mode, reset thread numbers for OpenMP parallelization */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset thread number */

      WDB[mat + MATERIAL_OMP_ID] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Neutron majorant ****************************************************/

  /* Check number of nuclides */

  if ((long)RDB[DATA_N_TRANSPORT_NUCLIDES] > 0)
    {
      fprintf(outp, "Calculating DT neutron majorant cross section:\n\n");

      /* Get pointer to energy grid */

      erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Number of points */

      ne = (long)RDB[erg + ENERGY_GRID_NE];

      /* Number of OpenMP threads */

      nt = (long)RDB[DATA_OMP_MAX_THREADS];

      /* Allocate memory for temporary majorant data */

      maj = (double **)Mem(MEM_ALLOC, nt, sizeof(double *));

      for(i = 0; i < nt; i++)
        {
          /* Allocate memory */

          maj[i] = (double *)Mem(MEM_ALLOC, ne, sizeof(double));
        }

      /* Start parallel timer */

      StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private (mat)
#endif
      {
        /* Print */

        PrintProgress(0, 0);

        /* Loop over materials */

        mat = (long)RDB[DATA_PTR_M0];
        while (mat > VALID_PTR)
          {

            /* Check if material is included in majorant */

            if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT))
              {
                /* Next material */

                mat = NextItem(mat);

                /* Cycle loop */

                continue;
              }

#ifdef OPEN_MP
#pragma omp critical
#endif
            {
              /* Grab material */

              if ((long)RDB[mat + MATERIAL_OMP_ID] == -1)
                WDB[mat + MATERIAL_OMP_ID] = OMP_THREAD_NUM;
            }

            /* Check thread number */

            if ((long)RDB[mat + MATERIAL_OMP_ID] == OMP_THREAD_NUM)
              {
                /* Process */

                CalculateDTMajorants0(mat, maj[OMP_THREAD_NUM], erg);

                /* Print */

                PrintProgress(mat, 1);
              }

            /* Next material */

            mat = NextItem(mat);
          }
      }

      /* Avoid compiler warning */

      i = -1;

      /* Compare data from parallel calculation */

#ifdef OPEN_MP
#pragma omp parallel private(i,n)
#endif
      {
#ifdef OPEN_MP
#pragma omp for
#endif

      for (n = 0; n < ne; n++)
        {
          for (i = 1; i < nt; i++)
            if (maj[i][n] > maj[0][n])
              maj[0][n] = maj[i][n];
        }
      }

      /* Stop parallel timer */

      StopTimer(TIMER_OMP_PARA);

      /* Get pointer to majorant reaction data */

      ptr = (long)RDB[DATA_PTR_MAJORANT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to data */

      ptr = (long)RDB[ptr + REACTION_PTR_MAJORANT_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Put data */

      memcpy(&WDB[ptr], maj[0], ne*sizeof(double));

      /* Find majorant minimum */

      WDB[DATA_MIN_NMACROXS] = INFTY;

      for (n = 0; n < ne; n++)
        if ((RDB[ptr + n] > 0.0) && (RDB[ptr + n] < RDB[DATA_MIN_NMACROXS]))
          WDB[DATA_MIN_NMACROXS] = RDB[ptr + n];

      /* Tää lasketaan eri tavalla kuin Serpent 1:ssä (sielläkin väärin) */
      /* etsi tolle joku järkevä arvo ja korvaa koko hässäkkä sillä */
      /*
      WDB[DATA_MIN_NMACROXS] = 1E-2;
      */
      /* Free allocated memory */

      for(i = 0; i < nt; i++)
        Mem(MEM_FREE, maj[i]);

      Mem(MEM_FREE, maj);

      /* Print */

      PrintProgress(0, 100);
    }

  /***************************************************************************/

  /***** Photon majorant *****************************************************/

  /* Check number of nuclides */

  if ((long)RDB[DATA_N_PHOTON_NUCLIDES] > 0)
    {
      fprintf(outp, "Calculating DT photon majorant cross section...\n");

      /* Check domain decomposition */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        Die(FUNCTION_NAME, "Domain decomposition in use");

      /* Get pointer to energy grid */

      erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Number of points */

      ne = (long)RDB[erg + ENERGY_GRID_NE];

      /* Pointer to energy array */

      loc0 = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Get pointer to majorant reaction data */

      loc1 = (long)RDB[DATA_PTR_PHOTON_MAJORANT];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get pointer to data */

      loc1 = (long)RDB[loc1 + REACTION_PTR_MAJORANT_XS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Reset data */

      memset(&WDB[loc1], 0.0, ne*sizeof(double));

      /* Allocate memory for temporary arrays */

      tot = (double *)Mem(MEM_ALLOC, ne, sizeof(double));
      xs = (double *)Mem(MEM_ALLOC, ne, sizeof(double));

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Get pointer to total */

          rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];

          /* Check total pointer and that material is included in majorant */

          if ((!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT)) ||
              (rea < VALID_PTR))
            {
              /* Next material */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Check if energy grids match */

          if ((long)RDB[rea + REACTION_PTR_EGRID] == erg)
            {
              /***************************************************************/

              /***** Same energy grid, copy total ****************************/

              /* Pointer to cross section array */

              ptr = (long)RDB[rea + REACTION_PTR_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Copy data */

              memcpy(tot, &RDB[ptr], ne*sizeof(double));

              /***************************************************************/
            }
          else
            {
              /***************************************************************/

              /***** No pre-calculated total, reconstruct ********************/

              /* Reset array */

              memset(tot, 0.0, ne*sizeof(double));

              /* Get pointer to partial list */

              ptr = (long)RDB[rea + REACTION_PTR_PARTIAL_LIST];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Reset reaction pointer (rewind list) */

              rea = -1;

              /* Loop over reactions */

              while ((rls = NextReaction(ptr, &rea, &adens, &Emin, &Emax, 0))
                     > VALID_PTR)
                {
                  /* We'll need to use zero atomic density for nuclides whose adens is */
                  /* given through a data interface as they will be added to the DT-   */
                  /* majorant during tracking */

                  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Check EXT_ADENS flag and set atomic density to zero if needed */

                  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS)
                    adens = 0.0;

                  /* Compare density to list maximum */

                  if (adens > RDB[rls + RLS_DATA_MAX_ADENS])
                    Die(FUNCTION_NAME, "density exceeds maximum 1");
                  else
                    adens = RDB[rls + RLS_DATA_MAX_ADENS];

                  /* Pointer to energy grid */

                  loc2 = (long)RDB[rea + REACTION_PTR_EGRID];
                  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                  /* Pointer to energy array */

                  loc2 = (long)RDB[loc2 + ENERGY_GRID_PTR_DATA];
                  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                  /* Pointer to cross section data */

                  loc3 = (long)RDB[rea + REACTION_PTR_XS];
                  CheckPointer(FUNCTION_NAME, "(loc3)", DATA_ARRAY, loc3);

                  /* First energy point and number of points */

                  i0 = (long)RDB[rea + REACTION_XS_I0];
                  sz = (long)RDB[rea + REACTION_XS_NE];

                  /* Reconstruct data */

                  InterpolateData(&RDB[loc0], xs, ne, &RDB[loc2 + i0],
                                  &RDB[loc3], sz, 0, NULL, NULL, NO);

                  /* Add to total */

                  for (n = 0; n < ne; n++)
                    tot[n] = tot[n] + adens*xs[n];
                }

              /***************************************************************/
            }

          /* Loop over array and compare */

          for (n = 0; n < ne; n++)
            if (tot[n] > RDB[loc1 + n])
              WDB[loc1 + n] = tot[n];

          /* Pointer to next */

          mat = NextItem(mat);
        }

      /* Find majorant minimum */

      WDB[DATA_MIN_PMACROXS] = INFTY;

      for (n = 0; n < ne; n++)
        if ((RDB[loc1 + n] > 0.0) &&
            (RDB[loc1 + n] < RDB[DATA_MIN_PMACROXS]))
          WDB[DATA_MIN_PMACROXS] = RDB[loc1 + n];

      /* Tää lasketaan eri tavalla kuin Serpent 1:ssä (sielläkin väärin) */
      /* etsi tolle joku järkevä arvo ja korvaa koko hässäkkä sillä */

      WDB[DATA_MIN_PMACROXS] = 1E-2;

      /* Free temporary arrays */

      Mem(MEM_FREE, tot);
      Mem(MEM_FREE, xs);

      fprintf(outp, "OK.\n\n");
    }

  /* Combine majorants with other domains */

  DistributeDDMajorants();

  /***************************************************************************/
}

/*****************************************************************************/

void CalculateDTMajorants0(long mat, double *maj, long erg)
{
  long rea, ne, ptr, loc0, loc1, loc2, i0, sz, n, n0, id, rls, nuc, idx, m, arr;
  double *tot, *xs, adens, Emin, Emax;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check if domain decomposition is in use */

  if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
      ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
      ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
    return;

  /* Get pointer to temperature majorant or total */

  if ((rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS]) < VALID_PTR)
    rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];

  /* If divided material and no pre-calculated data, include only first    */
  /* zone (pelkän indeksin testaaminen ei toimi jos ensimmäinen materiaali */
  /* ei ole fyysisesti geometriassa. Sitä varten tohon lisättiin pointteri */
  /* (7.8.2016 / 2.1.27 / JLe) */

  if ((long)RDB[rea + REACTION_PTR_EGRID] != erg)
    if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
      if ((ptr = (long)RDB[ptr + MATERIAL_DIV_PTR_FIRST]) > VALID_PTR)
        if (ptr != mat)
          return;

  /* Check energy grid pointer */

  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Check divisor type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Divided parent material");

  /* Number of points */

  ne = (long)RDB[erg + ENERGY_GRID_NE];

  /* Pointer to energy array */

  loc0 = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get OpenMP id */

  id = OMP_THREAD_NUM;

  /* Get pointers to temporary arrays */

  tot = WorkArray(DATA_PTR_WORK_PRIVA_GRID1 , PRIVA_ARRAY, ne, id);
  xs = WorkArray(DATA_PTR_WORK_PRIVA_GRID2, PRIVA_ARRAY, ne, id);

  /* Check if energy grids match */

  if ((long)RDB[rea + REACTION_PTR_EGRID] == erg)
    {
      /***********************************************************************/

      /***** Same energy grid, copy total ************************************/

      /* Pointer to cross section array */

      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Copy data */

      memcpy(tot, &RDB[ptr], ne*sizeof(double));

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** No pre-calculated total, reconstruct ****************************/

      /* Reset array */

      memset(tot, 0.0, ne*sizeof(double));

      /* Get pointer to partial list */

      ptr = (long)RDB[rea + REACTION_PTR_PARTIAL_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while ((rls = NextReaction(ptr, &rea, &adens, &Emin, &Emax, id))
             > VALID_PTR)
        {
          /* We'll need to use zero atomic density for nuclides whose adens is */
          /* given through a data interface as they will be added to the DT-   */
          /* majorant during tracking */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check EXT_ADENS flag and set atomic density to zero if needed */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS)
            adens = 0.0;

          /* Compare density to list maximum */

          if (adens > RDB[rls + RLS_DATA_MAX_ADENS])
            Die(FUNCTION_NAME, "density exceeds maximum");
          else
            adens = RDB[rls + RLS_DATA_MAX_ADENS];

          /* Get pointer to nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check poisons */

          if (((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES) &&
              ((long)RDB[nuc + NUCLIDE_ZAI] == 541350))
            if (adens > 0.0)
              Die(FUNCTION_NAME, "Xe-135 density is not zero");

          if (((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES) &&
              ((long)RDB[nuc + NUCLIDE_ZAI] == 621490))
            if (adens > 0.0)
              Die(FUNCTION_NAME, "Sm-149 density is not zero");

          /* Don't include critical iteration nuclides in majorant */

          if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_CRIT_ITER) &&
              ((long)RDB[mat + MATERIAL_PTR_ITER_ISO_LIST] > VALID_PTR))
            adens = 0.0;

          /* Set density to zero if nuclide is not to be included in material */
          /* totals and majorant due to atomic density being supplied         */
          /* separately */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS)
            {
              /* Loop through data interfaces to check if this material is */
              /* part of the list of materials for which this nuclides     */
              /* atomic density is supplied separately */
              if ((arr = (long)RDB[mat + MATERIAL_PTR_DATAIFC_ARR]) > VALID_PTR)
                {
                  /* Loop over data interfaces for this material */

                  while ((ptr = (long)RDB[arr]) > VALID_PTR)
                    {
                      /* ptr is a DATAIFC block, check nuclide */

                      if ((long)RDB[ptr + DATAIFC_PTR_NUCLIDE] == nuc)
                        {
                          /* Exclude from material total */

                          adens = 0.0;
                          /*
                          fprintf(outp, "Excluding %s from DT majorant (material %s)\n",
                                  GetText(nuc + NUCLIDE_PTR_NAME),
                                  GetText(mat + MATERIAL_PTR_NAME));
                          */
                          /* No need to check further */

                          break;
                        }
                      /* Next data interface */

                      arr++;
                    }
                }
            }


          /* Pointer to energy grid */

          loc1 = (long)RDB[rea + REACTION_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Pointer to energy array */

          loc1 = (long)RDB[loc1 + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Pointer to cross section data */

          loc2 = (long)RDB[rea + REACTION_PTR_XS];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* First energy point and number of points */

          i0 = (long)RDB[rea + REACTION_XS_I0];
          sz = (long)RDB[rea + REACTION_XS_NE];

          /* Do not reconstruct TMS majorants (histogram type */
          /* must be handled differently) */

          if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_NONE)
            {
              /* Reconstruct data */

              InterpolateData(&RDB[loc0], xs, ne, &RDB[loc1 + i0],
                              &RDB[loc2], sz, 0, NULL, NULL, NO);


              /* Add to total */

              for (n = 0; n < ne; n++)
                tot[n] = tot[n] + adens*xs[n];
            }

          /* Add contribution of the histogram majorant */

          else
            {
              idx = 0;

              for (n = 0; n < ne; n++)
                {
                  /* Find correct interval from nuclide reaction grid */

                  while(RDB[loc1 + idx + 1] < RDB[loc0 + n + 1])
                    idx++;

                  /* Add contribution */

                  tot[n] = tot[n] + adens*RDB[loc2 + idx];
                }
            }

        }

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Ures corrections ****************************************************/

  /* Vaihdettiin nuklidigridistä unionizediin 15.2.2016 / 2.1.25 (JLe) */

  /* Get pointer to ures list */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_TOT_URES_LIST]) > VALID_PTR)
    {
      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while ((rls = NextReaction(ptr, &rea, &adens, &Emin, &Emax, id))
             > VALID_PTR)
        {
          /* We'll need to use zero atomic density for nuclides whose adens is */
          /* given through a data interface as they will be added to the DT-   */
          /* majorant during tracking */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check EXT_ADENS flag and set atomic density to zero if needed */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS)
            adens = 0.0;

          /* Compare density to list maximum */

          if (adens > RDB[rls + RLS_DATA_MAX_ADENS])
            Die(FUNCTION_NAME, "density exceeds maximum");
          else
            adens = RDB[rls + RLS_DATA_MAX_ADENS];

          /* Check reaction pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Pointer to energy grid */

          loc1 = (long)RDB[rea + REACTION_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Pointer to energy array */

          loc1 = (long)RDB[loc1 + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Pointer to cross section data */

          loc2 = (long)RDB[rea + REACTION_PTR_XS];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* First energy point and number of points in data */

          i0 = (long)RDB[rea + REACTION_XS_I0];
          sz = (long)RDB[rea + REACTION_XS_NE];

          /* Reconstruct infinite-dilute data */

          InterpolateData(&RDB[loc0], xs, ne, &RDB[loc1 + i0], &RDB[loc2 + i0],
                          sz, 0, &n0, NULL, NO);

          /* First energy point and number of points in ures region */

          i0 = (long)RDB[rea + REACTION_URES_MAX_N0];
          sz = (long)RDB[rea + REACTION_URES_MAX_NP];

          /* Pointer to ures maximum data */

          loc2 = (long)RDB[rea + REACTION_PTR_URES_MAX];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Reset count */

          m = 0;

          /* adjust total */

          for (n = 0; n < sz; n++)
            {
              /* Check */

              if (xs[i0 + n] > RDB[loc2 + n])
                m++;

              /* Adjust */

              tot[i0 + n] = tot[i0 + n] - adens*xs[i0 + n] +
                1.01*adens*RDB[loc2 + n];
            }

#ifdef DEBUG

          /* Print warning */

          if (m > 0)
            Warn(FUNCTION_NAME,
                 "Infinite dilute exceeds ures max in %ld points", m);

#endif
        }
    }

  /***************************************************************************/

  /* Loop over array and compare */

  for (n = 0; n < ne; n++)
    if (tot[n] > maj[n])
      maj[n] = tot[n];
}

/*****************************************************************************/
