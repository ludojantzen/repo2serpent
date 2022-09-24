/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sampleifcdata.c                                */
/*                                                                           */
/* Created:       2017/12/15 (VVa)                                           */
/* Last modified: 2019/02/13 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Produces temperature and density data on user chosen         */
/*              cartesian grid.                                              */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleIFCData:"

/*****************************************************************************/

void SampleIFCData(long ini)
{
  long n, m, k, cell, mat, id, loc0;
  long nx, ny, nz, numsample;
  double ***Tmatrix, ***ADmatrix, ***MDmatrix;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w;
  double dummy, xt, yt, zt, f, T, t;
  FILE *fp;
  char fname[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check if any samples are defined */

  if ((long)RDB[DATA_PTR_SAMPLE0] < VALID_PTR)
    return;

  /* Print output */

  fprintf(outp, "Sampling temperature and density data to file.\n");

  /**************************************************************************/

  /***** Initialization *****************************************************/

  /* Reset dummy weight and source coordinates used for boundary conditions */

  dummy = -1.0;
  xt = 0.0;
  yt = 0.0;
  zt = 0.0;

  /* In dynamic mode use end time of the current time interval in sampling */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    t = RDB[DATA_TIME_CUT_TMAX];
  else
    t = -1.0;

  /* Get pointer to samples */

  loc0 = (long)RDB[DATA_PTR_SAMPLE0];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  numsample = 0;

  /* Loop over samples */

  while (loc0 > VALID_PTR)
    {
      /* Get sampling mesh-parameters */

      xmin = RDB[loc0 + SAMPLE_MESH_MIN0];
      xmax = RDB[loc0 + SAMPLE_MESH_MAX0];
      nx   = (long)RDB[loc0 + SAMPLE_MESH_N0];

      ymin = RDB[loc0 + SAMPLE_MESH_MIN1];
      ymax = RDB[loc0 + SAMPLE_MESH_MAX1];
      ny   = (long)RDB[loc0 + SAMPLE_MESH_N1];

      zmin = RDB[loc0 + SAMPLE_MESH_MIN2];
      zmax = RDB[loc0 + SAMPLE_MESH_MAX2];
      nz   = (long)RDB[loc0 + SAMPLE_MESH_N2];

      /* Allocate memory */

      Tmatrix = (double ***)Mem(MEM_ALLOC, nx, sizeof(double **));
      ADmatrix = (double ***)Mem(MEM_ALLOC, nx, sizeof(double **));
      MDmatrix = (double ***)Mem(MEM_ALLOC, nx, sizeof(double **));

      for(n = 0; n < nx; n++)
        {

          Tmatrix[n] = (double **)Mem(MEM_ALLOC, ny, sizeof(double *));
          ADmatrix[n] = (double **)Mem(MEM_ALLOC, ny, sizeof(double *));
          MDmatrix[n] = (double **)Mem(MEM_ALLOC, ny, sizeof(double *));

          for (m = 0; m < ny; m++)
            {
              Tmatrix[n][m] = (double *)Mem(MEM_ALLOC, nz, sizeof(double));
              ADmatrix[n][m] = (double *)Mem(MEM_ALLOC, nz, sizeof(double));
              MDmatrix[n][m] = (double *)Mem(MEM_ALLOC, nz, sizeof(double));
            }
        }

#ifdef OPEN_MP
#pragma omp parallel private(n, m, k, x, y, z, u, v, w, cell, mat, f, T, id)
#endif
      {
        /* Get Open MP thread id */

        id = OMP_THREAD_NUM;

        /* Avoid compiler warnings */

        x = 0;
        y = 0;
        z = 0;
        u = 0;
        v = 0;
        w = 1.0;

        /* Loop over geometry */

#ifdef OPEN_MP
#pragma omp for
#endif
        for (k = 0; k < nz; k++)
          {
            for (m = 0; m < ny; m++)
              {
                for (n = 0; n < nx; n++)
                  {
                    /* Reset pixel value */

                    Tmatrix[n][m][k] = 0;
                    ADmatrix[n][m][k] = 0;
                    MDmatrix[n][m][k] = 0;

                    /* Calculate Co-ordinates */


                    if (nx == 1)
                      x = xmin;
                    else
                      x = (n/(nx - 1.0))*(xmax - xmin) + xmin;

                    if (ny == 1)
                      y = ymin;
                    else
                      y = (m/(ny - 1.0))*(ymax - ymin) + ymin;

                    if (nz == 1)
                      z = zmin;
                    else
                      z = (k/(nz - 1.0))*(zmax - zmin) + zmin;

                    u = 0.0;
                    v = 0.0;
                    w = 1.0;

                    /* Find location */

                    if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                      BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                         &xt, &yt, &zt, &dummy, id);

                    /* Check return value */

                    if (cell < 0)
                      {
                        /* Set value */

                        Tmatrix[n][m][k] = -cell;
                        ADmatrix[n][m][k] = -cell;
                        MDmatrix[n][m][k] = -cell;
                      }
                    else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                      {
                        /* Get pointer */

                        mat = MatPtr(mat, id);
                        CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

                        /* Reset density and temperature */

                        f = 1.0;
                        T = 0.0;

                        /* Get density factor */

                        IFCPoint(mat, &f, &T, t, id);

                        /* Check for undefined density */

                        if (f < 0.0)
                          {
                            /* Calculation will currently use density factor of 1.0 */
                            /* if it cannot be obtained from the interface */

                            f = 1.0;
                          }

                        /* Check for undefined temperature */

                        if (T < ZERO)
                          {
                            /* Get Temperature */

                            if (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
                              {
                                /* Use minimum TMS temperature */

                                T = RDB[mat + MATERIAL_TMS_TMIN];
                              }
                            else
                              {

                                /* Use material doppler temperature */

                                T = RDB[mat + MATERIAL_DOPPLER_TEMP];
                              }

                          }

                        /* Store data to matrix (ADENS contains either material or mass density) */

                        Tmatrix[n][m][k] = T;
                        ADmatrix[n][m][k] = f*RDB[mat + MATERIAL_ADENS];
                        MDmatrix[n][m][k] = f*RDB[mat + MATERIAL_MDENS];
                      }

                  }
              }
          }
      }

      /* Print data to file */

      if (ini)
        sprintf(fname, "%s_sample%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                numsample + 1);
      else
        {
          if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
            {
              if ((long)RDB[DATA_RUN_CC] == YES)
                sprintf(fname, "%s_sample%ld_bstep%ld_iter%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                        numsample + 1, (long)RDB[DATA_BURN_STEP], (long)RDB[DATA_SOL_REL_ITER]);
              else
                sprintf(fname, "%s_sample%ld_bstep%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                        numsample + 1, (long)RDB[DATA_BURN_STEP]);
            }
          else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
            {
              if ((long)RDB[DATA_RUN_CC] == YES)
                sprintf(fname, "%s_sample%ld_tstep%ld_iter%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                        numsample + 1, (long)RDB[DATA_DYN_TB], (long)RDB[DATA_SOL_REL_ITER]);
              else
                sprintf(fname, "%s_sample%ld_tstep%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                        numsample + 1, (long)RDB[DATA_DYN_TB]);
            }
          else if ((long)RDB[DATA_RUN_CC] == YES)
            sprintf(fname, "%s_sample%ld_iter%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                      numsample+1, (long)RDB[DATA_SOL_REL_ITER]);
          else
            Die(FUNCTION_NAME, "Should not go here");
        }

      /* Open output file */

      if ((fp = fopen(fname, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Print coordinate vectors */

      fprintf(fp, "SAMPLE_X = [\n");

      if (nx == 1)
        fprintf(fp, "%12.5E ", xmin);
      else
        {
          for (n = 0; n < nx; n++)
            fprintf(fp, "%12.5E ", (n/(nx - 1.0))*(xmax - xmin) + xmin);
        }

      fprintf(fp, "\n];\n\n");

      fprintf(fp, "SAMPLE_Y = [\n");

      if (ny == 1)
        fprintf(fp, "%12.5E ", ymin);
      else
        {
          for (n = 0; n < ny; n++)
            fprintf(fp, "%12.5E ", (n/(ny - 1.0))*(ymax - ymin) + ymin);
        }

      fprintf(fp, "\n];\n\n");

      fprintf(fp, "SAMPLE_Z = [\n");

      if (nz == 1)
        fprintf(fp, "%12.5E ", zmin);
      else
        {
          for (n = 0; n < nz; n++)
            fprintf(fp, "%12.5E ", (n/(nz - 1.0))*(zmax - zmin) + zmin);
        }
      fprintf(fp, "\n];\n\n");

      /* Print temperature data */

      fprintf(fp, "sampled_T = [\n");

      for (k = 0; k < nz; k++)
        for (m = 0; m < ny; m++)
          for (n = 0; n < nx; n++)
            fprintf(fp, "%12.5E ", Tmatrix[n][m][k]);

      fprintf(fp, "\n];\n\n");

      fprintf(fp, "\nsampled_T = reshape(sampled_T, [%ld, %ld, %ld]);\n\n", nx, ny, nz);

      /* Print density data */

      fprintf(fp, "sampled_adens = [\n");

      for (k = 0; k < nz; k++)
        for (m = 0; m < ny; m++)
          for (n = 0; n < nx; n++)
            fprintf(fp, "%12.5E ", ADmatrix[n][m][k]);

      fprintf(fp, "\n];\n\n");

      fprintf(fp, "\nsampled_adens = reshape(sampled_adens, [%ld, %ld, %ld]);\n\n", nx, ny, nz);

      fprintf(fp, "sampled_mdens = [\n");

      for (k = 0; k < nz; k++)
        for (m = 0; m < ny; m++)
          for (n = 0; n < nx; n++)
            fprintf(fp, "%12.5E ", MDmatrix[n][m][k]);

      fprintf(fp, "\n];\n\n");

      fprintf(fp, "\nsampled_mdens = reshape(sampled_mdens, [%ld, %ld, %ld]);\n\n", nx, ny, nz);

      /* Close file */

      fclose(fp);

      /* Free T and D matrices */

      for(n = 0; n < nx; n++)
        {
          for (m = 0; m < ny; m++)
            {
              Mem(MEM_FREE, Tmatrix[n][m]);
              Mem(MEM_FREE, ADmatrix[n][m]);
              Mem(MEM_FREE, MDmatrix[n][m]);
            }

          Mem(MEM_FREE, Tmatrix[n]);
          Mem(MEM_FREE, ADmatrix[n]);
          Mem(MEM_FREE, MDmatrix[n]);
        }

      Mem(MEM_FREE, Tmatrix);
      Mem(MEM_FREE, ADmatrix);
      Mem(MEM_FREE, MDmatrix);

      /* Increment sample number */

      numsample++;

      /* Next sample */

      loc0 = NextItem(loc0);

    }

  fprintf(outp, "OK.\n\n");

  /****************************************************************************/

}

/*****************************************************************************/
