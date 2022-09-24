/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writedynsrc.c                                  */
/*                                                                           */
/* Created:       2015/05/15 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Writes out dynamic source at the end of the dynamic          */
/*              simulation                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteDynSrc:"

/*****************************************************************************/

void WriteDynSrc()
{
  long loc0, ptr, id, nlive, mul, N, det, stp, i, j, k, ng, nm, nt, group;
  double Wlive, P, norm;
  double x, y, z, u, v, w, E, wgt, t, Wtresh, *Wlist;


  /***************************************************************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* Check if output filename is given */

  if ((long)RDB[loc0 + PRECDET_PTR_OUT_FNAME] < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "savedynsrc.c -->\n");
#endif

  /******************************************************************/
  /***** Store the live neutron banks  ******************************/
  /******************************************************************/

  /* Reset total weight and source size */

  Wlive = 0.0;
  nlive = 0;

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while ((ptr = NextItem(ptr)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Add to total weight */

          Wlive = Wlive + RDB[ptr + PARTICLE_WGT];
          nlive++;
        }
    }

#ifdef DNPRINT
  fprintf(outp, "%ld live neutrons reached the end of simulation\n", nlive);
  fprintf(outp, "Mean weight was %E\n", Wlive/nlive);
#endif

  /* Proceed to write neutrons to <filename_base>.live */
  /* Normalize their weights to 1.0 */

  Wtresh = Wlive/(double)nlive;

  /* Get pointer to live neutron detector */

  det = (long)RDB[loc0 + PRECDET_PTR_FILE_DET];
  CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

  /* Loop over banks and write particles to file */
  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while ((ptr = NextItem(ptr)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Get particle weight */

          wgt = RDB[ptr + PARTICLE_WGT];
          x = RDB[ptr + PARTICLE_X];
          y = RDB[ptr + PARTICLE_Y];
          z = RDB[ptr + PARTICLE_Z];
          u = RDB[ptr + PARTICLE_U];
          v = RDB[ptr + PARTICLE_V];
          w = RDB[ptr + PARTICLE_W];
          E = RDB[ptr + PARTICLE_E];
          t = RDB[ptr + PARTICLE_T];

          if (1 == 1)
            {

              /* Normalize weights to 1.0 */
              if (wgt < Wtresh)
                {
                  /* Russian roulette */

                  /* Sample russian roulette */
                  /* Write to file if passed */

                  if (drand48() < wgt/Wtresh)
                    WriteSourceFile(det, x, y, z, u, v, w, E, 1.0, t, -1.0,
                                    id);

                }
              else if (wgt > Wtresh)
                {
                  /* Calculate multiplication */

                  P = wgt/Wtresh;
                  mul = (long)P;
                  P = P - (double)mul;

                  /* Sample extra point */

                  if (drand48() < P)
                    N = mul + 1;
                  else
                    N = mul;

                  /* Write to file N times */

                  for (i = 0; i < N; i++)
                    WriteSourceFile(det, x, y, z, u, v, w, E, 1.0, t, -1.0,
                                    id);

                }
              else
                WriteSourceFile(det, x, y, z, u, v, w, E, 1.0, t, -1.0,
                                id);

            }
          else
            {
              /* Do not normalize weights */

                WriteSourceFile(det, x, y, z, u, v, w, E, wgt, t, -1.0,
                                id);

            }
        }
    }

  /**********************************************************/
  /************ Save total live population ******************/
  /**********************************************************/

  /* Get normalization for neutrons */

  norm = RDB[DATA_NORM_COEF_N];

  /* Get pointer to live population detector */

  det = (long)RDB[loc0 + PRECDET_PTR_LIVE_DET];
  CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

  /* Pointer to statistics */

  stp = (long)RDB[det + DET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp1)", DATA_ARRAY, stp);

  /* Add statistics from this batch */
  /* This will be written to .main later by printprecdet.c */

  AddStat(Wlive*norm, stp, 0, 0);

#ifdef DNPRINT
  fprintf(outp, "Live population was %E neutrons\n", Wlive*norm);
#endif

  /**********************************************************/
  /*********** Store total precursor population *************/
  /**********************************************************/

  /* Get number of groups */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Allocate memory for temporaray list of weights */

  Wlist  = (double *)Mem(MEM_ALLOC, ng, sizeof(double));

  /* Loop over precursor source to calculate source population */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];

  /* Loop over source to calculate total weight */

  while ((ptr = NextItem(ptr)) > VALID_PTR)
    {
      /* Get group of precursor*/

      group = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Add to weight of this group*/

      Wlist[group] += RDB[ptr + PARTICLE_WGT];
    }

  /* Get number of time bins */

  nt = (long)RDB[loc0 + PRECDET_NT];

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Calculate number of mesh bins */

  nm = 1;
  nm = nm*(long)RDB[ptr + MESH_N0];
  nm = nm*(long)RDB[ptr + MESH_N1];
  nm = nm*(long)RDB[ptr + MESH_N2];

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp2)", DATA_ARRAY, stp);

  /* Loop over bins to print stats (calculate stable concentrations by dividing with lambda) */

  for (i = nt - 1; i < nt; i++)
    for (j = 0; j < ng; j++)
      {
        for (k = 0; k < nm; k++)
          {

            /* Store pointwise mesh populations */
            /* These would be list based populations */
            /* Should the output values be based on list or mesh? */

            /*
              AddStat(Wlist[j]*norm/(double)(nm), stp, i, j, k);
            */
          }
#ifdef DNPRINT
        fprintf(outp, "Precursor population for group was %E precursors (list)\n", Wlist[j]*norm);
#endif
      }

  /* Free temporary list */

  Mem(MEM_FREE, Wlist);

  /**********************************************************/
  /************ Store the current precursor list ************/
  /**********************************************************/

  /* Get pointer to detector */

  det = (long)RDB[loc0 + PRECDET_PTR_PREC_DET];
  CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

  /* Loop over precursor source to store particles */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];

  /* First item is dummy */

  while ((ptr = NextItem(ptr)) > VALID_PTR)
    {

      /* Get particle data */

      wgt = RDB[ptr + PARTICLE_WGT]/RDB[DATA_NORM_COEF_N];
      x = RDB[ptr + PARTICLE_X];
      y = RDB[ptr + PARTICLE_Y];
      z = RDB[ptr + PARTICLE_Z];
      t = RDB[ptr + PARTICLE_T];

      /* Get group of precursor*/

      group = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Write precursor to file */
      /* Should we cycle over all thread id's? */

      WriteSourceFile(det, x, y, z, 1.0, 0.0, 0.0, 1.0,
                      wgt, (double)group, -1.0, 0);

    }

#ifdef DNPRINT
  fprintf(outp, "<-- savedynsrc.c\n\n");
#endif
}
