/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writesourcefile.c                              */
/*                                                                           */
/* Created:       2012/10/19 (JLe)                                           */
/* Last modified: 2020/05/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Writes source distribution to file                           */
/*                                                                           */
/* Comments: - The distribution is actually stored into a buffer, and        */
/*             written only after the buffer is full.                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteSourceFile:"

/*****************************************************************************/

void WriteSourceFile(long det, double x, double y, double z, double u,
                     double v, double w, double E, double wgt, double t,
                     double flx, long id)
{
  long loc0, ptr, idx, sz, i, eof, skip;
  double spd, P, val;
  double x0, y0, z0, v0, u0, w0, wgt0, t0, E0;
  char fname[MAX_STR];
  FILE *fp, *fin;

#ifdef MPI

  int myturn;
  MPI_Status status;

#endif

  /* Get number of inactive batches */

  if(RDB[DATA_USE_FSP] == (double)NO)
    {
      /* No fission source passing*/

      skip = (long)RDB[DATA_CRIT_SKIP];
    }
  else if ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
         && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    {
      /* First transportcycle with fission source passing */

      skip = (long)RDB[DATA_CRIT_SKIP];
    }
  else
    {
      /* Subsequent transportcycle with fission source passing */

      skip = (long)RDB[DATA_FSP_CRIT_SKIP];
    }

  /* Do not write points during inactive cycles */

  if (WDB[DATA_CYCLE_IDX] < skip)
    return;

  /* Check pointer */

  if (det > VALID_PTR)
    {
      /***********************************************************************/

      /***** Write point in buffer *******************************************/

      /* Get pointer to buffer */

      if ((loc0 = (long)RDB[det + DET_WRITE_PTR_BUF]) < VALID_PTR)
        return;

      /* Calculate probability */

      if (flx < 0.0)
        {
          /* Store all points with probability DET_WRITE_PROB  */
          /* Surface crossing or Fission points                */

          P = RDB[det + DET_WRITE_PROB];
        }
      else
        {
          /* storing porbability weighted with neutron density */
          /* weighted with 1/v*flx    (flx = 1/totxs)          */

          spd = Speed(PARTICLE_TYPE_NEUTRON, E);

          /* MinXS can depend on energy so we don't want to use it */
          /* to normalize */
          /*
          P = RDB[det + DET_WRITE_PROB]*wgt*
            Speed(PARTICLE_TYPE_NEUTRON, RDB[DATA_NEUTRON_EMIN])/spd*
            MinXS(PARTICLE_TYPE_NEUTRON, spd, id)*flx;
          */

          P = RDB[det + DET_WRITE_PROB]*wgt*
            Speed(PARTICLE_TYPE_NEUTRON, RDB[DATA_NEUTRON_EMIN])/spd*
            flx;
        }

      if (P > 1)
        Warn(FUNCTION_NAME, "P larger than 1 (%E)", P);

      /* Rejection */

      if (RandF(id) > P)
        return;

      /* Get buffer size */

      sz = (long)RDB[det + DET_WRITE_BUF_SZ];

      /* Set OpenMP barrier */

#ifdef OPEN_MP
#pragma omp critical
#endif
      {
        /* Get index */

        idx = (long)RDB[det + DET_WRITE_BUF_IDX];
        CheckValue(FUNCTION_NAME, "idx", "", idx, 0, sz - 1);

        /* Get pointer */

        ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Put data */

        WDB[ptr + SRC_BUF_X] = x;
        WDB[ptr + SRC_BUF_Y] = y;
        WDB[ptr + SRC_BUF_Z] = z;
        WDB[ptr + SRC_BUF_U] = u;
        WDB[ptr + SRC_BUF_V] = v;
        WDB[ptr + SRC_BUF_W] = w;
        WDB[ptr + SRC_BUF_E] = E;
        WDB[ptr + SRC_BUF_WGT] = wgt;
        WDB[ptr + SRC_BUF_T] = t;

        /* Update index */

        idx++;

        /* Check if buffer is full */

        if (idx == sz)
          {
            /* Print filename */

            if (mpitasks > 1)
              sprintf(fname, "%s%d",GetText(det + DET_WRITE_PTR_FILE), mpiid);
            else
              sprintf(fname, "%s",GetText(det + DET_WRITE_PTR_FILE));

            /* Open file */

            if ((fp = fopen(fname, "a")) == NULL)
              Error(det, "Unable to open file \"%s\" for writing",
                    GetText(det + DET_WRITE_PTR_FILE));

            /* Loop over buffer */

            for (idx = 0; idx < sz; idx++)
              {
                /* Get pointer */

                ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                  {
                    /* Write binary data */

                    val = RDB[ptr + SRC_BUF_X];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_Y];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_Z];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_U];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_V];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_W];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_E];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_WGT];
                    fwrite(&val, sizeof(double), 1, fp);

                    val = RDB[ptr + SRC_BUF_T];
                    fwrite(&val, sizeof(double), 1, fp);
                  }
                else
                  {
                    /* Write ASCII data */

                    fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
                    fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
                    fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
                    fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
                    fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
                    fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
                    fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
                    fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
                    fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
                    fprintf(fp, "\n");
                  }
              }

            /* Close file */

            fclose(fp);

            /* Reset index */

            idx = 0;
          }

        /* Put new index */

        WDB[det + DET_WRITE_BUF_IDX] = (double)idx;
      }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Dump all buffers ************************************************/

#ifdef DNPRINT
  fprintf(outp, "writesourcefile.c (dumping buffers) -->\n\n");
#endif

#ifdef MPI

      myturn = 0;

      /* First task has its turn first */

      if ( mpiid == 0 )
        myturn = 1;
      /* Other tasks wait for go-ahead from the previous task */
      else
        MPI_Recv(&myturn, 1, MPI_INT, mpiid-1, 1, my_comm, &status);

      /* Check the go-ahead */

      if (myturn == 1)
        {

#endif
          /* Check OpenMP thread number */

          if (OMP_THREAD_NUM != 0)
            Die(FUNCTION_NAME, "Called from an OpenMP parallel loop");

          /* Loop over detectors */

          det = (long)RDB[DATA_PTR_DET0];
          while (det > VALID_PTR)
            {
              /* Check pointer to write buffer */

              if ((loc0 = (long)RDB[det + DET_WRITE_PTR_BUF]) > VALID_PTR)
                {
                  /* Get buffer size */

                  sz = (long)RDB[det + DET_WRITE_BUF_IDX];

                  /* Print filename */

                  sprintf(fname, "%s",GetText(det + DET_WRITE_PTR_FILE));

                  /* Open file */

                  if ((fp = fopen(fname, "a")) == NULL)
                    Error(det, "Unable to open file \"%s\" for writing",
                          GetText(det + DET_WRITE_PTR_FILE));

                  /* Loop over buffer */

                  for (idx = 0; idx < sz; idx++)
                    {
                      /* Get pointer */

                      ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                        {
                          /* Write binary data */

                          val = RDB[ptr + SRC_BUF_X];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_Y];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_Z];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_U];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_V];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_W];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_E];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_WGT];
                          fwrite(&val, sizeof(double), 1, fp);

                          val = RDB[ptr + SRC_BUF_T];
                          fwrite(&val, sizeof(double), 1, fp);
                        }
                      else
                        {
                          /* Write ASCII data */

                          fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
                          fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
                          fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
                          fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
                          fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
                          fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
                          fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
                          fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
                          fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
                          fprintf(fp, "\n");
                        }

                    }

                  /* Add the points written to task-specific file during */
                  /* simulation */

                  if (mpitasks > 1)
                    {
                      /* Print filename */

                      sprintf(fname, "%s%d",GetText(det + DET_WRITE_PTR_FILE),
                              mpiid);

                      /* Open file if it has been created */

                      if ((fin = fopen(fname, "r")) != NULL)
                        {
                          while (1 == 1)
                            {
                              if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                                {
                                  /* Read and write binary data */

                                  if (fread(&val, sizeof(double), 1, fin) == 0)
                                    break;
                                  else
                                    fwrite(&val, sizeof(double), 1, fp);
                                }
                              else
                                {
                                  /* Write ASCII data */

                                  if ((eof = fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x0, &y0, &z0, &u0, &v0, &w0, &E0, &wgt0, &t0)) == EOF)
                                    break;
                                  else
                                    fprintf(fp, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %11.5E %11.5E %11.5E\n", x0, y0, z0, u0, v0, w0, E0, wgt0, t0);
                                }
                            }

                          /* Close file */

                          fclose(fin);

                          /* Remove file */

                          remove(fname);
                        }
                    }

                  /* Close file */

                  fclose(fp);

                  /* Reset index */

                  WDB[det + DET_WRITE_BUF_IDX] = 0.0;
                }
              /* Next source */

              det = NextItem(det);
            }

          /* Criticality source detector */

          if ((det = (long)RDB[DATA_PTR_CRIT_SRC_DET]) > VALID_PTR)
            {
              /* Get pointer to write buffer */

              loc0 = (long)RDB[det + DET_WRITE_PTR_BUF];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Get buffer size */

              sz = (long)RDB[det + DET_WRITE_BUF_IDX];

              /* Print filename */

              sprintf(fname, "%s",GetText(det + DET_WRITE_PTR_FILE));

              /* Open file */

              if ((fp = fopen(fname, "a")) == NULL)
                Error(det, "Unable to open file \"%s\" for writing",
                      GetText(det + DET_WRITE_PTR_FILE));

              /* Loop over buffer */

              for (idx = 0; idx < sz; idx++)
                {
                  /* Get pointer */

                  ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                    {
                      /* Write binary data */

                      val = RDB[ptr + SRC_BUF_X];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_Y];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_Z];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_U];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_V];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_W];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_E];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_WGT];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_T];
                      fwrite(&val, sizeof(double), 1, fp);
                    }
                  else
                    {
                      /* Write ASCII data */

                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
                      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
                      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
                      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
                      fprintf(fp, "\n");
                    }
                }

              /* Add the points written to task-specific file during */
              /* simulation */

              if (mpitasks > 1)
                {
                  /* Print filename */

                  sprintf(fname, "%s%d",GetText(det + DET_WRITE_PTR_FILE),
                          mpiid);

                  /* Open file if it has been created */

                  if ((fin = fopen(fname, "r")) != NULL)
                    {
                      while (1 == 1)
                        {
                          if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                            {
                              /* Read and write binary data */

                              if (fread(&val, sizeof(double), 1, fin) == 0)
                                break;
                              else
                                fwrite(&val, sizeof(double), 1, fp);
                            }
                          else
                            {
                              /* Write ASCII data */

                              if ((eof = fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x0, &y0, &z0, &u0, &v0, &w0, &E0, &wgt0, &t0)) == EOF)
                                break;
                              else
                                fprintf(fp, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %11.5E %11.5E %11.5E\n", x0, y0, z0, u0, v0, w0, E0, wgt0, t0);
                            }
                        }

                      /* Close file */

                      fclose(fin);

                      /* Remove file */

                      remove(fname);
                    }
                }

              /* Close file */

              fclose(fp);

              /* Reset index */

              WDB[det + DET_WRITE_BUF_IDX] = 0.0;
            }

          /* Secondary photon source */

          if ((det = (long)RDB[DATA_PTR_NGAMMA_SRC_DET]) > VALID_PTR)
            {
              /* Get pointer to write buffer */

              loc0 = (long)RDB[det + DET_WRITE_PTR_BUF];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Get buffer size */

              sz = (long)RDB[det + DET_WRITE_BUF_IDX];

              /* Print filename */

              sprintf(fname, "%s",GetText(det + DET_WRITE_PTR_FILE));

              /* Open file */

              if ((fp = fopen(fname, "a")) == NULL)
                Error(det, "Unable to open file \"%s\" for writing",
                      GetText(det + DET_WRITE_PTR_FILE));

              /* Loop over buffer */

              for (idx = 0; idx < sz; idx++)
                {
                  /* Get pointer */

                  ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                    {
                      /* Write binary data */

                      val = RDB[ptr + SRC_BUF_X];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_Y];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_Z];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_U];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_V];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_W];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_E];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_WGT];
                      fwrite(&val, sizeof(double), 1, fp);

                      val = RDB[ptr + SRC_BUF_T];
                      fwrite(&val, sizeof(double), 1, fp);
                    }
                  else
                    {
                      /* Write ASCII data */

                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
                      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
                      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
                      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
                      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
                      fprintf(fp, "\n");
                    }
                }

              /* Add the points written to task-specific file during */
              /* simulation */

              if (mpitasks > 1)
                {
                  /* Print filename */

                  sprintf(fname, "%s%d",GetText(det + DET_WRITE_PTR_FILE),
                          mpiid);

                  /* Open file if it has been created */

                  if ((fin = fopen(fname, "r")) != NULL)
                    {
                      while (1 == 1)
                        {
                          if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                            {
                              /* Read and write binary data */

                              if (fread(&val, sizeof(double), 1, fin) == 0)
                                break;
                              else
                                fwrite(&val, sizeof(double), 1, fp);
                            }
                          else
                            {
                              /* Write ASCII data */

                              if ((eof = fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x0, &y0, &z0, &u0, &v0, &w0, &E0, &wgt0, &t0)) == EOF)
                                break;
                              else
                                fprintf(fp, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %11.5E %11.5E %11.5E\n", x0, y0, z0, u0, v0, w0, E0, wgt0, t0);
                            }
                        }

                      /* Close file */

                      fclose(fin);

                      /* Remove file */

                      remove(fname);
                    }
                }

              /* Close file */

              fclose(fp);

              /* Reset index */

              WDB[det + DET_WRITE_BUF_IDX] = 0.0;
            }

          /* Precursor detector */

#ifdef DNPRINT
          fprintf(outp, "Checking the existence of a precursor detector\n");
#endif

          if ((long)RDB[DATA_PTR_PREC_DET] > VALID_PTR)
            for (i = 0; i < 2; i++)
              {
                if (i == 0)
                  {
                    loc0 = (long)RDB[DATA_PTR_PREC_DET];

                    if ((det = (long)RDB[loc0 + PRECDET_PTR_PREC_DET]) <
                        VALID_PTR)
                      {
#ifdef DNPRINT
                        fprintf(outp, "No precursor detector linked.\n");
#endif

                        continue;
                      }
                  }
                else if (i == 1)
                  {
                    loc0 = (long)RDB[DATA_PTR_PREC_DET];

                    if ((RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) ||
                        ((det = (long)RDB[loc0 + PRECDET_PTR_FILE_DET])
                         < VALID_PTR))
                      {
#ifdef DNPRINT
                        fprintf(outp, "No live neutron detector linked.\n");
#endif

                        continue;
                      }
                  }

                /* Get pointer to write buffer */

                loc0 = (long)RDB[det + DET_WRITE_PTR_BUF];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Get buffer size */

                sz = (long)RDB[det + DET_WRITE_BUF_IDX];

                /* Print filename */

                sprintf(fname, "%s",GetText(det + DET_WRITE_PTR_FILE));

                /* Open file */

                if ((fp = fopen(fname, "a")) == NULL)
                  Error(det, "Unable to open file \"%s\" for writing",
                        GetText(det + DET_WRITE_PTR_FILE));

#ifdef DNPRINT
                fprintf(outp, "Proceeding to store %ld particles to %s.\n", sz, fname);
#endif
                /* Loop over buffer */

                for (idx = 0; idx < sz; idx++)
                  {
                    /* Get pointer */

                    ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                      {
                        /* Write binary data */

                        val = RDB[ptr + SRC_BUF_X];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_Y];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_Z];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_U];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_V];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_W];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_E];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_WGT];
                        fwrite(&val, sizeof(double), 1, fp);

                        val = RDB[ptr + SRC_BUF_T];
                        fwrite(&val, sizeof(double), 1, fp);
                      }
                    else
                      {
                        /* Write ASCII data */

                        fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
                        fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
                        fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
                        fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
                        fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
                        fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
                        fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
                        fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
                        fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
                        fprintf(fp, "\n");
                      }
                  }

                /* Add the points written to task-specific file during */
                /* simulation */

                if (mpitasks > 1)
                  {
                    /* Print filename */

                    sprintf(fname, "%s%d",GetText(det + DET_WRITE_PTR_FILE),
                            mpiid);

                    /* Open file if it has been created */

                    if ((fin = fopen(fname, "r")) != NULL)
                      {
                        while (1 == 1)
                          {
                            if ((long)RDB[det + DET_WRITE_BINARY] == YES)
                              {
                                /* Read and write binary data */

                                if (fread(&val, sizeof(double), 1, fin) == 0)
                                  break;
                                else
                                  fwrite(&val, sizeof(double), 1, fp);
                              }
                            else
                              {
                                /* Read and write ASCII data */

                                if ((eof = fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x0, &y0, &z0, &u0, &v0, &w0, &E0, &wgt0, &t0)) == EOF)
                                  break;
                                else
                                  fprintf(fp, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %11.5E %11.5E %11.5E\n", x0, y0, z0, u0, v0, w0, E0, wgt0, t0);
                              }
                          }

                        /* Close file */

                        fclose(fin);

                        /* Remove file */

                        remove(fname);

                      }
                  }

                /* Close file */

                fclose(fp);

                /* Reset index */

                WDB[det + DET_WRITE_BUF_IDX] = 0.0;
              }
#ifdef MPI
        }

      /* All detectors written, give go-ahead to next task */

      if (mpiid != mpitasks-1)
        MPI_Send(&myturn, 1, MPI_INT, mpiid+1, 1, my_comm);
#endif

#ifdef DNPRINT
  fprintf(outp, "<-- writesourcefile.c (dumping buffers)\n\n");
#endif

      /***********************************************************************/
    }
}

/*****************************************************************************/
