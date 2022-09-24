/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcovariancedata.c                        */
/*                                                                           */
/* Created:       2018/06/04 (VVa)                                           */
/* Last modified: 2018/09/12 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes covariance data to a form from which it can be     */
/*              easily scored during neutron tracking                        */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessCovarianceData:"

void permuteZAImt(long * list, long *list2, long **COVxy, long i, long j, long N)
{
  long tmp, rowi, coli, n;

  /* Permute the ZAImt-list */

  tmp = list[i];
  list[i] = list[j];
  list[j] = tmp;

  /* Permute the second list if needed */

  if (list2 != NULL)
    {
      tmp = list2[i];
      list2[i] = list2[j];
      list2[j] = tmp;
    }

  /* Swap diagonal first */

  tmp = COVxy[i][i];
  COVxy[i][i] = COVxy[j][j];
  COVxy[j][j] = tmp;

  /* Swap the rest of the rows and columns here */

  for (n = 0; n < N; n++)
    {

      if ((n == i) || (n == j))
        continue;

      rowi = COVxy[i][n];
      coli = COVxy[n][i];

      COVxy[i][n] = COVxy[j][n];
      COVxy[n][i] = COVxy[n][j];

      COVxy[j][n] = rowi;
      COVxy[n][j] = coli;
    }

  return;
}

/* Depth first recursive algorithm for marking all elements in a connected */
/* block for the covariance matrix */

void AddToConBlock(long *conblock, long **COVxy, long i, long Nmin, long Nmax)
{
  long j;

  /* Mark this element as already connected so that it will not be processed */
  /* again */

  conblock[i] = 1;

  /* Loop over elements between Nmin and Nmax that have not been connected */
  /* to this block already and see if they should be connected */

  for (j = Nmin; j < Nmax; j++)
    if ((COVxy[i][j] != 0) && (conblock[j] == 0))
      AddToConBlock(conblock, COVxy, j, Nmin, Nmax);
}

/*****************************************************************************/

void ProcessCovarianceData()
{
  long loc0, loc1, idx, *blocksizes;
  long ncov, zai, mt, zai1, zai2, mt1, mt2, zaimt, zaimt1, zaimt2, nuniq, n, m, i, j;
  long maxsz, ng, pos, cursz, block, ptr, ptr1, ptr2, egrid;
  long matarr, matrix, sens, izai, imat, irea, iene, blockidx;
  long *ZAImts, *conblock, **COVxy, noncon, connections, con, spot, processed, nb;

  if ((loc0 = (long)RDB[DATA_PTR_COVMTX0]) < VALID_PTR)
    return;

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  fprintf(outp, "Processing covariance data...\n\n");

  /* Calculate the max number of unique reaction modes in the covariance data */

  ncov = ListSize(loc0);

  fprintf(outp, "%ld covariance matrices can include %ld unique reactions\n", ncov, 2*ncov);

  ZAImts = (long *)Mem(MEM_ALLOC, ncov+1, sizeof(long));

  /**********************************************************************/
  /* Get unique covariance matrices that are needed for the calculation */
  /**********************************************************************/

  GetRequiredCovMatrices(ZAImts);

  /* Get the pointer to the first covariance matrix or return */
  /* In some cases all of the covariance matrices might have been removed */

  if ((loc0 = (long)RDB[DATA_PTR_COVMTX0]) < VALID_PTR)
    return;

  /* Get number of groups in covariance matrix (currently all matrices need to */
  /* have the same number of energy groups */

  ng = (long)RDB[loc0 + COVMTX_NG];

  /* Check number of energy groups  */

  if (ng != (long)RDB[sens + SENS_N_ENE])
    Die(FUNCTION_NAME, "Covariance data has %ld energy groups while sensitivity "
        "calculation uses %ld. These should be equal. You'll most likely want to "
        "modify \"sens opt egrid\".", ng, (long)RDB[sens + SENS_N_ENE]);

  /* */ /* */ /* */ /* */
  /* Actually we should also check that the energy grids match here */
  /* */ /* */ /* */ /* */

  /* Get number of groups in covariance matrix (currently all matrices need to */
  /* have the same number of energy groups */

  egrid = (long)RDB[loc0 + COVMTX_PTR_EGRID];

  /* Get unique reaction modes */

  while (loc0 > VALID_PTR)
    {
      /* Check number of energy groups */

      if (ng != (long)RDB[loc0 + COVMTX_NG])
        Die(FUNCTION_NAME, "Different covariance matrices have different numbers "
            "of energy groups (%ld groups vs %ld groups). Such data is not yet "
            "supported.", ng, (long)RDB[loc0 + COVMTX_NG]);

      /* Check energy grid */

      if (egrid != (ptr = (long)RDB[loc0 + COVMTX_PTR_EGRID]))
        {
          /* The energy grid doesn't match, check that the grid structure is */
          /* still the same */

          /* Check the number of points in the grid (should be OK since the */
          /* number of groups was already checked)                          */

          if ((n = (long)RDB[egrid + ENERGY_GRID_NE]) != (long)RDB[ptr + ENERGY_GRID_NE])
            Die(FUNCTION_NAME, "Different covariance matrices have energy grids "
                "with different number of points (%ld points vs %ld points). "
                "Such data is not yet supported.", n, (long)RDB[ptr + ENERGY_GRID_NE]);
          else
            {
              /* Check grid points */

              ptr1 = (long)RDB[egrid + ENERGY_GRID_PTR_DATA];
              ptr2 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];

              for (i = 0; i < n; i++)
                {
                  if (RDB[ptr1 + i] != RDB[ptr2 + i])
                    Die(FUNCTION_NAME, "Different covariance matrices have different "
                        "energy grids. Grid point %ld is %E MeV in first grid, "
                        "%E MeV in the second.", i+1, RDB[ptr1+i], RDB[ptr2+i]);

                }
            }
        }

      /* Get unique zai-mt pairs from covariance matrices */

      for (i = 0; i < 2; i++)
        {
          /* Get ZAI and mt*/

          if (i == 0)
            {
              zai = (long)RDB[loc0 + COVMTX_ZAI1];
              mt = (long)RDB[loc0 + COVMTX_MT1];
            }
          else
            {
              zai = (long)RDB[loc0 + COVMTX_ZAI2];
              mt = (long)RDB[loc0 + COVMTX_MT2];
            }

          /* Put these into a single variable */

          zaimt = zai*10000 + mt;

          /* Find from array */

          n = 0;
          while (ZAImts[n] > 1001*10000)
            {
              /* Compare to current element of the array */

              if (ZAImts[n] == zaimt)
                break;

              /* Next element of the array */

              n++;
            }

          /* If the current element of the array is zero, we didn't find  */
          /* what we were searching for, which means that we should store */

          ZAImts[n] = zaimt;
        }

      /* Next covariance matrix */

      loc0 = NextItem(loc0);
    }

  /* Loop and print all the unique ZAImt-pairs */
  fprintf(outp, "Unique ZAI-mt pairs:\n");
  n = 0;
  while (ZAImts[n] > 1001*10000)
    {
      /* Print out the current element of the array */

      fprintf(outp, "  %ld\n", ZAImts[n]);

      /* Next element of the array */

      n++;
    }

  fprintf(outp, "%ld in total. There can be %ld*%ld = %ld cross covariances, but we know that there are only %ld.\n", n, n, n, n*n, ncov);

  nuniq = n;

  COVxy = (long **)Mem(MEM_ALLOC, nuniq, sizeof(long *));

  for (n = 0; n < nuniq; n++)
    COVxy[n] = (long *)Mem(MEM_ALLOC, nuniq, sizeof(long *));

  /* Check which of the cross covariances have data */

  for (i = 0; i < nuniq; i++)
    {
      for (j = 0; j < nuniq; j++)
        {
          /* Get zais and mts corresponding to this cross covariance */

          zai1 = (long)(ZAImts[i]/10000);
          mt1 = ZAImts[i]-zai1*10000;

          zai2 = (long)(ZAImts[j]/10000);
          mt2 = ZAImts[j]-zai2*10000;

          /* Get unique reaction modes */

          loc0 = (long)RDB[DATA_PTR_COVMTX0];

          while (loc0 > VALID_PTR)
            {
              /* Get ZAI and mt*/

              zai = (long)RDB[loc0 + COVMTX_ZAI1];
              mt = (long)RDB[loc0 + COVMTX_MT1];

              /* Put these into a single variable */

              zaimt1 = zai*10000 + mt;

              /* second pair */

              zai = (long)RDB[loc0 + COVMTX_ZAI2];
              mt = (long)RDB[loc0 + COVMTX_MT2];

              /* Put these into a single variable */

              zaimt2 = zai*10000 + mt;

              /* Check current i,j pair against covariance matrix */

              if ((zaimt1 == ZAImts[i]) && (zaimt2 == ZAImts[j]))
                break;

              /* Check current j,i pair against covariance matrix */

              if ((zaimt1 == ZAImts[j]) && (zaimt2 == ZAImts[i]))
                break;

              loc0 = NextItem(loc0);
            }

          /* Check if covariance matrix for current i,j pair was found */

          if (loc0 > VALID_PTR)
            {
              fprintf(outp, "*");
              COVxy[i][j] = loc0;
            }
          else
            {
              fprintf(outp, " ");
              COVxy[i][j] = 0;
            }
        }
      fprintf(outp, "\n");
    }

  /**************************************************************/
  /* Take the non-connected diagonals and put them in the front */
  /**************************************************************/

  noncon = 0;

  for (i = 0; i < nuniq; i++)
    {
      /* Check if this is connected to something else than itself */

      connections = 0;

      for (j = 0; j < nuniq; j++)
        if (COVxy[i][j] > 0)
          connections++;

      /* Skip processing for connected */

      if (connections > 1)
        continue;

      /* Permute non-connected ones to the front */

      permuteZAImt(ZAImts, NULL, COVxy, noncon, i, nuniq);

      /* Increment number of non-connected ones found */

      noncon++;
    }

  /********************************************************/
  /* Print out and build the modified connectivity matrix */
  /********************************************************/

  for (i = 0; i < nuniq; i++)
    {
      for (j = 0; j < nuniq; j++)
        {
          /* Get zais and mts corresponding to this cross covariance */

          zai1 = (long)(ZAImts[i]/10000);
          mt1 = ZAImts[i]-zai1*10000;

          zai2 = (long)(ZAImts[j]/10000);
          mt2 = ZAImts[j]-zai2*10000;

          /* Get unique reaction modes */

          loc0 = (long)RDB[DATA_PTR_COVMTX0];

          while (loc0 > VALID_PTR)
            {
              /* Get ZAI and mt*/

              zai = (long)RDB[loc0 + COVMTX_ZAI1];
              mt = (long)RDB[loc0 + COVMTX_MT1];

              /* Put these into a single variable */

              zaimt1 = zai*10000 + mt;

              /* second pair */

              zai = (long)RDB[loc0 + COVMTX_ZAI2];
              mt = (long)RDB[loc0 + COVMTX_MT2];

              /* Put these into a single variable */

              zaimt2 = zai*10000 + mt;

              /* Check current i,j pair against covariance matrix */

              if ((zaimt1 == ZAImts[i]) && (zaimt2 == ZAImts[j]))
                break;

              /* Check current j,i pair against covariance matrix */

              if ((zaimt1 == ZAImts[j]) && (zaimt2 == ZAImts[i]))
                break;

              loc0 = NextItem(loc0);
            }

          /* Check if covariance matrix for current i,j pair was found */

          if (loc0 > VALID_PTR)
            {
              if (COVxy[i][j]  < VALID_PTR)
                fprintf(outp,  "Ö");
              else
                fprintf(outp, "*");
            }
          else
            {
              if (COVxy[i][j] != 0)
                fprintf(outp, "Ä");
              else
                fprintf(outp, " ");

            }
        }
      fprintf(outp, "\n");
    }

  /***********************************************************************/
  /* Make a depth first search to permute everything to connected blocks */
  /***********************************************************************/

  /* Create a list of the block sizes */

  blocksizes = (long *)Mem(MEM_ALLOC, nuniq, sizeof(long));

  /* Mark block size of 1 for the non-connected elements in the beginning */

  for (nb = 0; nb < noncon; nb++)
    blocksizes[nb] = 1;

  /* Create an empty list into which collect each block before permuting */

  conblock = (long *)Mem(MEM_ALLOC, nuniq, sizeof(long));

  i = noncon;

  while (i < nuniq)
    {
      /* Add flags to current connected */

      AddToConBlock(conblock, COVxy, i, i, nuniq);

      /* Calculate number of connections */

      connections = 0;
      for (j = 0; j < nuniq; j++)
        if (conblock[j] != 0)
          connections++;

      /* Store block size */

      blocksizes[nb++] = connections;

      /************************************************************/
      /* Permute elements so that the connected block is together */
      /************************************************************/

      /* Get spot into which the first element can be permuted */
      /* (last connected block + 1) */

      spot = i;

      /* Process all the connections */

      processed = 0;
      while (processed < connections)
        {

          /* Find first connection (the previous for loop could just save */
          /* a list of these). There probably aren't that many */

          con = 0;
          for (j = i; j < nuniq; j++)
            if ((con = conblock[j]) > 0)
              break;

          if (con <= 0)
            Die(FUNCTION_NAME, "Couldn't find the connection for permuting");

          /* Process this connection now */

          con = j;
          conblock[j] = 0;
          processed++;

          /* If it is already at the spot, just take the next connection */

          if (con == spot)
            {
              spot++;
              continue;
            }

          /* Permute this connected element to the spot */

          permuteZAImt(ZAImts, conblock, COVxy, spot, con, nuniq);

          /* The next element will take the next spot */

          spot++;
        }

      /* Increment i which keeps track of number of elements at their */
      /* correct places */

      i += processed;

      /* Zero out the connections block */

      memset(conblock, 0, nuniq*sizeof(long));
    }

  /********************************************************/
  /* Print out and build the modified connectivity matrix */
  /********************************************************/

  fprintf(outp, "---\n");
  for (i = 0; i < nuniq; i++)
    {
      for (j = 0; j < nuniq; j++)
        {
          /* Get zais and mts corresponding to this cross covariance */

          zai1 = (long)(ZAImts[i]/10000);
          mt1 = ZAImts[i]-zai1*10000;

          zai2 = (long)(ZAImts[j]/10000);
          mt2 = ZAImts[j]-zai2*10000;

          /* Get unique reaction modes */

          loc0 = (long)RDB[DATA_PTR_COVMTX0];

          while (loc0 > VALID_PTR)
            {
              /* Get ZAI and mt*/

              zai = (long)RDB[loc0 + COVMTX_ZAI1];
              mt = (long)RDB[loc0 + COVMTX_MT1];

              /* Put these into a single variable */

              zaimt1 = zai*10000 + mt;

              /* second pair */

              zai = (long)RDB[loc0 + COVMTX_ZAI2];
              mt = (long)RDB[loc0 + COVMTX_MT2];

              /* Put these into a single variable */

              zaimt2 = zai*10000 + mt;

              /* Check current i,j pair against covariance matrix */

              if ((zaimt1 == ZAImts[i]) && (zaimt2 == ZAImts[j]))
                break;

              /* Check current j,i pair against covariance matrix */

              if ((zaimt1 == ZAImts[j]) && (zaimt2 == ZAImts[i]))
                break;

              loc0 = NextItem(loc0);
            }

          /* Check if covariance matrix for current i,j pair was found */

          if (loc0 > VALID_PTR)
            {
              if (COVxy[i][j] < VALID_PTR)
                fprintf(outp,  "Ö");
              else
                fprintf(outp, "*");
            }
          else
            {
              if (COVxy[i][j] != 0)
                fprintf(outp, "Ä");
              else
                fprintf(outp, " ");

            }
        }
      fprintf(outp, "\n");
    }

  /*****************************************************/
  /* Now our global matrix is in a block diagonal form */
  /* with the size of each block in some variable      */
  /*****************************************************/

  /* Get maximum block size */

  maxsz = 0;
  for (i = 0; i < nuniq; i++)
    if (maxsz < blocksizes[i])
      maxsz = blocksizes[i];

  /* Check maximum block size */

  if (maxsz < 1)
    Die(FUNCTION_NAME, "Maximum connected block size is %ld. "
        "Should be at least 1.", maxsz);

  fprintf(outp, "Maximum block size is %ld\n", maxsz);

  /* Get number of neutron groups */

  loc0 = (long)RDB[DATA_PTR_COVMTX0];
  ng = (long)RDB[loc0 + COVMTX_NG];

  /**************************************************************/
  /* Process each block and do the following procedures         */
  /*  1. Create the real matrix corresponding to the block      */
  /*  2. Make the matrix symmetric: A = 0.5*(A + A^T)           */
  /*  3. Deflate negative eigenpairs out of the matrix to make  */
  /*     it positive semidefinite.                              */
  /*  4. Eigendecompose the symmetric positive semidefinite     */
  /*     matrix.                                                */
  /*  5. Sum up the basis functions to be scored and link so    */
  /*     that the value to be scored can be easily found during */
  /*     transport                                              */
  /**************************************************************/
  blockidx = 0;
  if (1==2)
    {
      /* Reset the position in the ZAImt list */

      pos = 0;

      for (i = 0; i < nuniq; i++)
        {
          /* Check if all blocks have been processed */

          if ((cursz = blocksizes[i]) == 0)
            break;

          /* Allocate memory for the matrix that contains the data for the block */

          matrix = ReallocMem(DATA_ARRAY, (cursz*ng)*(cursz*ng));

          /* Create the (symmetric) data matrix from the block  */
          /* Maybe pass the COVBLOCK to this routine and let it */
          /* allocate, create and store the data matrix for it? */

          CovMatrixFromBlock(ZAImts, COVxy, pos, cursz, matrix, blockidx);
          blockidx++;
          /* Store data as a covariance block */

          block = NewItem(DATA_PTR_COVBLOCK0, COVBLOCK_BLOCK_SIZE);

          /* Store pointer to the data matrix*/

          WDB[block + COVBLOCK_PTR_MATRIX_DATA] = (double)matrix;

          /* Store block matrix order to the block */

          WDB[block + COVBLOCK_ORDER] = (double)cursz;

          /* Allocate memory for the ZAI-mt array */

          loc1 = ReallocMem(DATA_ARRAY, cursz);

          /* Store pointer to the ZAI-mt array */

          WDB[block + COVBLOCK_PTR_ZAIMT_ARRAY] = (double)loc1;

          /* Copy the data for the ZAI-mt array */

          for (n = 0; n < cursz; n++)
            WDB[loc1 + n] = (double)ZAImts[pos + n];

          /* Allocate memory for the covariance matrix pointer array */

          loc1 = ReallocMem(DATA_ARRAY, cursz*cursz);

          /* Store pointer to the covariance matrix pointer array */

          WDB[block + COVBLOCK_PTR_COVMTX_ARRAY] = (double)loc1;

          /* Copy the data for the covariance matrix pointer array */

          for (n = 0; n < cursz; n++)
            for (m = 0; m < cursz; m++)
              WDB[loc1 + n*cursz + m] = (double)COVxy[pos + n][pos + m];

          /* Store number of energy groups */

          WDB[block + COVBLOCK_NG] = (double)ng;

          /* Store pointer to energy grid */

          WDB[block + COVBLOCK_PTR_EGRID] = RDB[loc0 + COVMTX_PTR_EGRID];

          /* Allocate memory for the sensitivity index array */

          loc1 = ReallocMem(DATA_ARRAY, cursz*ng);

          /* Store pointer to sensitivity index array */

          WDB[block + COVBLOCK_PTR_SENS_INDICES] = (double)loc1;

          /* Get pointer to sensitivity block or don't link */
          /* Put this into a separate subroutine?           */

          if ((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
            {
              /* Get pointer to material index array */

              matarr = (long)RDB[sens + SENS_PTR_MAT_INDICES];
              CheckPointer(FUNCTION_NAME, "(matarr)", DATA_ARRAY, matarr);

              /* Fill the sensitivity index array */

              for (n = 0; n < cursz; n++)
                {
                  /* Get ZAI and mt */

                  zai = (long)(ZAImts[pos + n]/10000);
                  mt = ZAImts[pos + n]-zai*10000;

                  /******************/
                  /* Find zai index */
                  /******************/

                  izai = FindSensZAIIndex(zai);

                  /* Check ZAI index */

                  if (izai < 0)
                    Die(FUNCTION_NAME, "Could not find sensitivity index for ZAI "
                        "%ld used in uncertainty propagation", zai);

                  /***********************/
                  /* Find reaction index */
                  /***********************/

                  if (!((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT))
                    Die(FUNCTION_NAME, "MT-based XS perturbation flag was not on");

                  /* Get pointer to pre-built array */

                  ptr = (long)RDB[sens + SENS_PTR_MT_INDICES];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Read index from pre-built array */

                  irea = (long)RDB[ptr + mt];

                  /* Check reaction index */

                  if (irea < 0)
                    Die(FUNCTION_NAME, "Could not find reaction index for MT "
                        "%ld used in uncertainty propagation", mt);

                  /*****************************/
                  /* Find total material index */
                  /*****************************/

                  imat = -1;
                  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
                    imat = (long)RDB[matarr + SENS_TOT_MAT_IDX];
                  else
                    Die(FUNCTION_NAME, "Total material is not flagged for scoring "
                        "of sensitivities even though covariance data is supplied");

                  /* Loop over energy groups and store index */

                  for (m = 0; m < ng; m++)
                    {
                      /* Get the energy group index for sensitivity (0 == total) */

                      iene = m + 1;

                      /* Get the index */

                      idx = CompressSensLabel(imat, izai, irea, iene, 1);

                      /* Store direct index */

                      WDB[loc1 + n*ng + m] = (double)idx;
                    }
                }
            }

          /* Deflate the negative eigenpairs out of the matrix */
          /* (requires solving the eigenpairs first)           */

          /*
            MakeCovMatrixPSD();
          */

          /* Eigendecompose the symmetric positive semidef. covariance  */
          /* matrix and sum up the eigenvectors*eigenvalue to be scored */
          /* during transport */

          /*
            ProcessCovMatrixToScorable();
          */

          /* Increment position */

          pos += cursz;

        }
    }
  else
    {
      /* Each ZAI-mt *2 pair is its own block */

      for (i = 0; i < nuniq; i++)
        {
          /* Only handle upper triangle part of the full COV matrix*/

          for (j = i; j < nuniq; j++)
            {
              /* Check if there is covariance data for this pair */

              if (COVxy[i][j] > VALID_PTR)
                {
                  /*************************************/
                  /* Create a block based on this pair */
                  /*************************************/

                  if (i == j)
                    cursz = 1;
                  else
                    cursz = 2;

                  /* Allocate memory for the matrix that contains the data */
                  /* for the block */

                  matrix = ReallocMem(DATA_ARRAY, (cursz*ng)*(cursz*ng));

                  /* Get the data matrix */

                  CovMatrixFromSingle(ZAImts, COVxy, i, j, matrix, blockidx);
                  blockidx++;
                  /* Store data as a covariance block */

                  block = NewItem(DATA_PTR_COVBLOCK0, COVBLOCK_BLOCK_SIZE);

                  /* Store pointer to the data matrix*/

                  WDB[block + COVBLOCK_PTR_MATRIX_DATA] = (double)matrix;

                  /* Store block matrix order to the block */

                  WDB[block + COVBLOCK_ORDER] = (double)cursz;

                  /* Allocate memory for the ZAI-mt array */

                  loc1 = ReallocMem(DATA_ARRAY, cursz);

                  /* Store pointer to the ZAI-mt array */

                  WDB[block + COVBLOCK_PTR_ZAIMT_ARRAY] = (double)loc1;

                  /* Copy the data for the ZAI-mt array */

                  WDB[loc1 + 0] = (double)ZAImts[i];

                  if (cursz > 1)
                    WDB[loc1 + 1] = (double)ZAImts[j];

                  /* Allocate memory for the covariance matrix pointer array */

                  loc1 = ReallocMem(DATA_ARRAY, cursz*cursz);

                  /* Store pointer to the covariance matrix pointer array */

                  WDB[block + COVBLOCK_PTR_COVMTX_ARRAY] = (double)loc1;

                  /* Copy the data for the covariance matrix pointer array */

                  if (cursz == 1)
                    WDB[loc1] = (double)COVxy[i][i];
                  else
                    {
                      WDB[loc1 + 0*cursz + 1] = (double)COVxy[i][j];
                      WDB[loc1 + 1*cursz + 0] = (double)COVxy[j][i];
                    }

                  /* Store number of energy groups */

                  WDB[block + COVBLOCK_NG] = (double)ng;

                  /* Store pointer to energy grid */

                  WDB[block + COVBLOCK_PTR_EGRID] = RDB[loc0 + COVMTX_PTR_EGRID];

                  /* Allocate memory for the sensitivity index array */

                  loc1 = ReallocMem(DATA_ARRAY, cursz*ng);

                  /* Store pointer to sensitivity index array */

                  WDB[block + COVBLOCK_PTR_SENS_INDICES] = (double)loc1;

                  /* Get pointer to sensitivity block or don't link */
                  /* Put this into a separate subroutine?           */

                  if ((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
                    {
                      /* Get pointer to material index array */

                      matarr = (long)RDB[sens + SENS_PTR_MAT_INDICES];
                      CheckPointer(FUNCTION_NAME, "(matarr)", DATA_ARRAY, matarr);

                      /* Fill the sensitivity index array */

                      for (n = 0; n < cursz; n++)
                        {
                          /* Get ZAI and mt */

                          if (n == 0)
                            {
                              zai = (long)(ZAImts[i]/10000);
                              mt = ZAImts[i]-zai*10000;
                            }
                          else
                            {
                              zai = (long)(ZAImts[j]/10000);
                              mt = ZAImts[j]-zai*10000;
                            }

                          /******************/
                          /* Find zai index */
                          /******************/

                          izai = FindSensZAIIndex(zai);

                          /* Check ZAI index */

                          if (izai < 0)
                            Die(FUNCTION_NAME, "Could not find sensitivity index for ZAI "
                                "%ld used in uncertainty propagation", zai);

                          /***********************/
                          /* Find reaction index */
                          /***********************/

                          if (!((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT))
                            Die(FUNCTION_NAME, "MT-based XS perturbation flag was not on");

                          /* Get pointer to pre-built array */

                          ptr = (long)RDB[sens + SENS_PTR_MT_INDICES];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          /* Read index from pre-built array */

                          irea = (long)RDB[ptr + mt];

                          /* Check reaction index */

                          if (irea < 0)
                            Die(FUNCTION_NAME, "Could not find reaction index for MT "
                                "%ld used in uncertainty propagation", mt);

                          /*****************************/
                          /* Find total material index */
                          /*****************************/

                          imat = -1;
                          if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
                            imat = (long)RDB[matarr + SENS_TOT_MAT_IDX];
                          else
                            Die(FUNCTION_NAME, "Total material is not flagged for scoring "
                                "of sensitivities even though covariance data is supplied");

                          /* Loop over energy groups and store index */

                          for (m = 0; m < ng; m++)
                            {
                              /* Get the energy group index for sensitivity (0 == total) */

                              iene = m + 1;

                              /* Get the index */

                              idx = CompressSensLabel(imat, izai, irea, iene, 1);

                              /* Store direct index */

                              WDB[loc1 + n*ng + m] = (double)idx;
                            }
                        }
                    }

                }
            }
        }
    }

  /* Free temporary arrays */

  Mem(MEM_FREE, ZAImts);
  Mem(MEM_FREE, conblock);
  Mem(MEM_FREE, blocksizes);

  for (n = 0; n < nuniq; n++)
    Mem(MEM_FREE, COVxy[n]);

  Mem(MEM_FREE, COVxy);

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
