/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initprecdet.c                                  */
/*                                                                           */
/* Created:       2015/05/12 (VVa)                                           */
/* Last modified: 2018/12/12 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Initializes values for precursor detectors and calculates    */
/*              weight to emit on first time interval                        */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitPrecDet:"

/*****************************************************************************/

void InitPrecDet()
{
  long loc0, ptr, pos, part, maxi;
  long i, j, k, id;
  long stp, tbin, gbin, ibin, src;
  long ng, ntf, ngf;
  long n0, n1, n2, n0f, n1f, n2f;
  double val, relerr, lambda, *filelams;
  double x, y, z, u, v, w, E, wgt, t, atot[2], t0, normf;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* If no precursor file has been defined we will return */

  if (RDB[loc0 + PRECDET_PTR_IN_FNAME] < VALID_PTR)
    return;

  /* We don't need to initialize values for criticality source mode */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    return;

#ifdef DNPRINT
  fprintf(outp, "initprecdet.c -->\n");

  /* Print precursor tracking mode */

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    fprintf(outp, "Tracking individual point-wise precursors\n");
  else
    fprintf(outp, "Tracking precursor concentrations on a regular mesh\n");

#endif

  /* Get initial simulation time */

  t0 = RDB[DATA_TIME_CUT_TMIN];

  /* Get normalization factor (default = 1.0) */

  normf = RDB[DATA_TRANS_NORM_FACT];

  /* Print main file name to tmpstr */
  /* File name */

  sprintf(tmpstr, "%s.main", GetText(loc0 + PRECDET_PTR_IN_FNAME));

  /* Open main input file */

  if ((fp = fopen(tmpstr,"r")) == NULL)
    Die(FUNCTION_NAME, "Could not open dyn source file %s for reading", tmpstr);

  /* Read live neutron number and relerr from file */

  if (fscanf(fp, "%lf %lf", &val, &relerr) != 2)
    Die(FUNCTION_NAME, "Could not read values from precursor file");

  /* Store live neutron value */

  WDB[loc0 + PRECDET_W_LIVE] = val*normf;

  /* Read binning from precursor detector */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of mesh bins */

  n0 = (long)RDB[ptr + MESH_N0];
  n1 = (long)RDB[ptr + MESH_N1];
  n2 = (long)RDB[ptr + MESH_N2];

  /* Read binning from file */

  if (fscanf(fp, "%ld %ld %ld %ld %ld\n", &ntf, &ngf, &n0f, &n1f, &n2f) != 5)
    Die(FUNCTION_NAME, "Could not read binning from precursor file");

  /* Check that binning is the same */

  if ((n0 != n0f) || (n1 != n1f) || (n0 != n1f))
    {
      Die(FUNCTION_NAME,
          "Mesh binning of precursors different (%ld %ld %ld in the file)"
          " vs (%ld %ld %ld in this simulation).",
          n0f, n1f, n2f, n0, n1, n2);
    }
  else if (ng != ngf)
    {
      Die(FUNCTION_NAME, "Number of delayed neutron groups changed since "
          "source generation (%ld in the source files) vs (%ld in this "
          "simulation). Potential cause is generating source distributions "
          " with one library and running the transient with another.",
          ngf, ng);
    }

  /* Read time from file */

  if (fscanf(fp, "%lf", &t) != 1)
    Die(FUNCTION_NAME, "Could not read current time from file");

  /* Check that the lambdas of the groups are the same */

  /* Read lambdas from file */

  filelams  = (double *)Mem(MEM_ALLOC, ngf, sizeof(double));

  for (j = 0; j < ngf; j++)
    if(fscanf(fp,"%lf", &(filelams[j])) != 1)
      Die(FUNCTION_NAME, "Could not read lambda from file");

  /* Get pointer to stored lambdas */

  ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];

  /* Check that the lambdas are equal */

  for (j = 0; j < ngf; j++)
    if (fabs(RDB[ptr + j] - filelams[j]) / RDB[ptr + j] > 1e-5)
      Die(FUNCTION_NAME,
          "Different decay constants in file and in this simulation");

  /* Free the temporary array */

  free(filelams);

  /* Close main file */

  fclose(fp);

  /**************************************************/
  /* Read the precursor populations from .prec file */
  /**************************************************/

  /* File name */

  sprintf(tmpstr, "%s.prec", GetText(loc0 + PRECDET_PTR_IN_FNAME));

  /* Read input file */

  if ((fp = fopen(tmpstr,"r")) == NULL)
    Die(FUNCTION_NAME, "Could not open dyn source file %s for reading", tmpstr);

#ifdef DNPRINT
  fprintf(outp, "Reading precursor populations from %s\n", tmpstr);
#endif

  /* Loop over the lines to get to the last timestep */

  for (i = 0; i < ntf - 1; i++)
    for (j = 0; j < ng; j++)
      for (k = 0; k < n0*n1*n2; k++)
        if (fscanf(fp, "%ld %ld %ld %lf %lf",
                   &tbin, &gbin, &ibin, &val, &relerr) != 5)
          Die(FUNCTION_NAME, "Could not read values from precursor file");

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);


  atot[0] = 0.0;

  /* Loop over the lines of the last timestep and get values */

  for (j = 0; j < ng; j++)
    for (k = 0; k < n0*n1*n2; k++)
      {
        if (fscanf(fp, "%ld %ld %ld %lf %lf",
                   &tbin, &gbin, &ibin, &val, &relerr) != 5)
          Die(FUNCTION_NAME, "Could not read values from precursor file");

        /* Modify normalization if needed */

        val = val*normf;

        /* These are steady state populations */

        /* Store stable values to first time bin of the detector */

        AddBuf(val, 1.0, stp, 0, -1, 0, gbin, ibin);

        /*
        fprintf(outp, "Stable number of precursors for group %ld is %E\n", gbin
        + 1, val);
        */

        /* Get pointer to stored lambdas */

        ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];

        /* Get lambda from precursor group */

        lambda = RDB[ptr + gbin];

        /* Add to total activity */

        atot[0] += val*lambda;

      }

  fclose(fp);

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    {
      /*************************************/
      /* Read precursor points from a file */
      /*************************************/

      /* File name */

      sprintf(tmpstr, "%s.precpoints", GetText(loc0 + PRECDET_PTR_IN_FNAME));

      /* Read input file */

      if ((fp = fopen(tmpstr,"r")) == NULL)
        Die(FUNCTION_NAME, "Could not open dyn source file %s for reading",
            tmpstr);

      /* Get pointer to precursor source definition */

      src = (long)RDB[loc0 + PRECDET_PTR_PREC_SRC];

      /* Get pointer to source */

      pos = (long)RDB[DATA_PART_PTR_PSOURCE];
      CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);

      /* Reset id */

      id = 0;

      /* Get 15*nsrc neutrons from file  */
      /* We'll have to get a large amount, because in a long time interval */
      /* the short-lived ones are not so important */
      /* TODO: more? less?               */

      maxi = (long)RDB[DATA_SRC_POP]*(long)RDB[DATA_PREC_SRC_FACT];

#ifdef MPI_MODE1

      /* Check number of tasks */

      if (mpitasks > 1)
        maxi = (long)((double)maxi/(double)mpitasks);

#endif
      /* Reset total activity for point-wise precursors */

      atot[1] = 0.0;

#ifdef DNPRINT
      fprintf(outp, "Reading precursors from source file\n");
#endif

      /* Get pointer to first precursor (dummy) */

      part = (long)RDB[DATA_PART_PTR_PSOURCE];

      /* Check that precursor source is empty */

      if (ListSize(part) != 1)
        Die(FUNCTION_NAME, "Precursor source not empty");

      /* Get precursors from source file */
      /* OMP parallelization? */

      for (i = 0 ; i < maxi ; i++)
        {
          ReadSourceFile(src, &x, &y, &z, &u, &v, &w, &E, &wgt, &t);

          /* Time coordinate is used for group bin */

          gbin = (long)t;

          /* Get new precursor from stack */

          part = FromStack(PARTICLE_TYPE_PRECURSOR, id++);

          /* Check id */

          if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
            id = 0;

          /* Get lambda from precursor group */

          /* Get pointer to stored lambdas */

          ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];

          /* Get lambda from precursor group */

          lambda = RDB[ptr + gbin];

          /* Put values */

          WDB[part + PARTICLE_X] = x;
          WDB[part + PARTICLE_Y] = y;
          WDB[part + PARTICLE_Z] = z;

          /* Note: Precursor DN group corresponds to "global" i.e. non-     */
          /* nuclide wise group and may be larger than 6/8, starting from 0 */

          WDB[part + PARTICLE_DN_GROUP] = (double)gbin;
          WDB[part + PARTICLE_DN_LAMBDA] = lambda;

          WDB[part + PARTICLE_WGT] = wgt;

          /* Put initial time of the precursor */
          /* this is the beginning of the first time interval */

          WDB[part + PARTICLE_T0] = t0;

          /* Put current time of the precursor             */
          /* Will be changed when the precursor is decayed */
          /* Over time-intervals */

          WDB[part + PARTICLE_T0] = t0;

          /* Put RNG index */
          /* Used for sorting particles to get reproducible MPI */
          /* Calculations */

          WDB[part + PARTICLE_RNG_IDX] = (double)i
            + (double)maxi*RDB[DATA_CYCLE_IDX];

          /* Put precursor to precursor source */

          if ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == YES)
            AddSortItem(DATA_PART_PTR_PSOURCE, pos, part, PARTICLE_RNG_IDX,
                        SORT_MODE_ASCEND);
          else
            AddItem(DATA_PART_PTR_PSOURCE, part);

          /* Update position */

          pos = part;
        }

      /* Close precursor point file */

      fclose(fp);

    }

#ifdef DNPRINT
  fprintf(outp, "OK.\n");
  fprintf(outp, "<-- initprecdet.c\n\n");
#endif

}
