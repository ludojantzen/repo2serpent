/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resizedynsrc.c                                 */
/*                                                                           */
/* Created:       2015/18/05 (VVa)                                           */
/* Last modified: 2018/09/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Resizes live source population to PRECDET_N_LIVE at time     */
/*              interval boundaries                                          */
/*                                                                           */
/* Comments:   -Number to normalize is calculated based on live weight and   */
/*              weight to emit on the upcoming interval.                     */
/*             -Live weight is calculated in countdynsrc.c                   */
/*             -Weight to emit is calculated in decayprecdet.c               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResizeDynSrc:"
#define DNPRINT
/*****************************************************************************/

void ResizeDynSrc()
{
  long loc0, ptr, part, id, N, m, mat;
  long nbatch, nsrc, np;
  long Nlive, Nemit;
  double wgt, wgt0, P, mul, x, y, z, E, t;
  double Wlive, Wemit, avewgt;


  /***************************************************************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "resizedynsrc.c -->\n");
#endif

  /* Get live weight at interval boundary */

  Wlive = RDB[loc0 + PRECDET_W_LIVE];

  /* Get weight to emit over next interval */

  Wemit = RDB[loc0 + PRECDET_W_EMIT];

#ifdef DNPRINT
  fprintf(outp, "Wlive %E\n", Wlive);
  fprintf(outp, "Wemit %E\n", Wemit);
#endif

  /* Get wanted criticality population */

  nbatch = (long)RDB[DATA_SRC_POP];

  /* Calculate number of live neutrons and neutrons to emit */

  Nlive = (long)round((double)nbatch*Wlive/(Wlive + Wemit));
  Nemit = (long)round((double)nbatch*Wemit/(Wlive + Wemit));

#ifdef DNPRINT
  fprintf(outp, "Nlive %ld\n", Nlive);
  fprintf(outp, "Nemit %ld\n", Nemit);
#endif

  /* Store numbers */

  WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
  WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

  /* Calculate average weight */

  avewgt = (Wemit + Wlive)/RDB[DATA_NORM_COEF_N]/(double)nbatch;

  /* Store average weight */

  WDB[loc0 + PRECDET_W_AVE] = avewgt;

  /* Reset total weight and source size */

  wgt0 = 0.0;
  nsrc = 0;

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

      while((ptr = NextItem(ptr)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Add to total weight and source size */

          wgt0 = wgt0 + RDB[ptr + PARTICLE_WGT];
          nsrc = nsrc + 1;

        }
    }

#ifdef DNPRINT
  fprintf(outp, "%ld particles in banks\n", nsrc);
#endif

  /***** Population control **************************************************/

  /* Reset weight and number of particles */

  wgt = wgt0;
  np = nsrc;

  /* Compare population size to number of live neutrons */

  if (nsrc < Nlive)
    {
      /* Calculate multiplication */

      P = ((double)Nlive)/((double)nsrc);
      mul = (double)((long)P);
      P = P - (double)mul;

      /* Loop over threads */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          /* Get particles from bank */

          /* Get pointer to bank */

          ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get pointer to last item */
          /* We have to loop backwards since new neutrons are  */
          /* added to the end of the bank and we don't want to */
          /* multiply them */

          ptr = LastItem(ptr);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          while(ptr > VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
                {
                  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
                    break;
                  else
                    Die(FUNCTION_NAME, "Invalid particle type");
                }

              /* Sample multiplication */

              if (RandF(0) < P)
                N = (long)mul;
              else
                N = (long)mul - 1;

              /* Loop over multiplication */

              for (m = 0; m < N; m++)
                {
                  /* Duplicate neutron */

                  part = DuplicateParticle(ptr, id);

                  /* Add to weight and population size */

                  wgt = wgt + RDB[ptr + PARTICLE_WGT];
                  np = np + 1;

                  /* Put particle in source*/

                  ToBank(part, id);
                }

              /* Next particle */

              ptr = PrevItem(ptr);

            }

        }
    }
  else if (nsrc > Nlive)
    {
      /* Calculate probability */

      P = 1.0 - ((double)Nlive)/((double)nsrc);

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

          ptr = NextItem(ptr);

          while(ptr > VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
                Die(FUNCTION_NAME, "Invalid particle type");

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

                  /* Put particle back to stack */

                  ToStack(part, id);

                }
            }
        }
    }

  /* Check */

  if (np < 1)
    {
      /* All initial neutrons on next interval will be emitted */

      /* Calculate number of live neutrons and neutrons to emit */

      Nlive = 0;
      Nemit = nbatch;

#ifdef DNPRINT
      fprintf(outp, "Nlive %ld\n", Nlive);
      fprintf(outp, "Nemit %ld\n", Nemit);
#endif

      /* Store numbers */

      WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
      WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

    }
  else if (np > nbatch)
    {
      /* All initial neutrons on next interval will be live */

      /* Calculate number of live neutrons and neutrons to emit */

      Nlive = np;
      Nemit = 0;

#ifdef DNPRINT
      fprintf(outp, "Nlive %ld\n", Nlive);
      fprintf(outp, "Nemit %ld\n", Nemit);
#endif

      /* Store numbers */

      WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
      WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

    }
  else
    {

      /* Calculate eventual number of live neutrons and neutrons to emit */

      Nlive = np;
      Nemit = nbatch - np;

#ifdef DNPRINT
      fprintf(outp, "Nlive %ld\n", Nlive);
      fprintf(outp, "Nemit %ld\n", Nemit);
#endif

      /* Store numbers */

      WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
      WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

    }

  /***************************************************************************/

  /***** Normalize source ****************************************************/

  /* Normalize weights */

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

      while((ptr = NextItem(ptr)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Normalize */

          WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*wgt0/wgt;
        }
    }

  /***************************************************************************/

  /***** Score live neutron primary source ***********************************/

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

      while((ptr = NextItem(ptr)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Get required information from particle */

          x = RDB[ptr + PARTICLE_X];
          y = RDB[ptr + PARTICLE_Y];
          z = RDB[ptr + PARTICLE_Z];

          E = RDB[ptr + PARTICLE_E];

          t = RDB[ptr + PARTICLE_T];

          wgt = RDB[ptr + PARTICLE_WGT];

          /* We could find this out if required */

          mat = NULLPTR;

          /* Score time-dependent source rate detectors */

          ScoreTimeSource(ptr, mat, MT_PRIMARY_LIVE_SOURCE,
                          x, y, z, E, t+ZERO, 1, wgt, id);

        }
    }

#ifdef DNPRINT
  fprintf(outp, "<-- resizedynsrc.c\n\n");
#endif

}

/*****************************************************************************/
