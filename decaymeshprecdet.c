/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : decaymeshprecdet.c                             */
/*                                                                           */
/* Created:       2015/05/15 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates decay of precursors and emission of delayed       */
/*              neutrons for the next interval for the mesh based mode       */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DecayMeshPrecDet:"

/*****************************************************************************/

void DecayMeshPrecDet()
{
  long loc0, ptr, vplist;
  long stp, j, k; 
  long ng, tb;
  long n0, n1, n2, zeroes, nonzeroes;
  double val, lambda, t0, t1, dt, Wemit, Wtal, Wall, Wmax;
  

  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* We don't need to do this for criticality source mode */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    return;

  /* If we are not tracking any precursors, we can return */

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_NONE)
    return;

#ifdef DNPRINT
  fprintf(outp, "decaymeshprecdet.c -->\n");  
#endif

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  /* Calculate length of time interval */

  dt = t1 - t0;

  /* Get current time bin */

  tb = (long)RDB[DATA_DYN_TB];
  
  /* Read group binning from precursor detector */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Get pointer to precursor spatial mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of mesh bins */

  n0 = (long)RDB[ptr + MESH_N0];
  n1 = (long)RDB[ptr + MESH_N1];
  n2 = (long)RDB[ptr + MESH_N2];    

  /* Get pointer to decay constant array */

  ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];

  /* Reset total weight to emit */

  Wall = 0.0;
  Wmax = 0.0;
  Wemit = 0.0;

  zeroes = 0;
  nonzeroes = 0;

  /* Loop over spatial and group bins to calculate */
  /* values to emit and tally during interval      */

  for (j = 0; j < ng; j++)
    for (k = 0; k < n0*n1*n2; k++)
      {

        /* Close buffer for reading */
                  
        WDB[DATA_BUF_REDUCED] = (double)YES;

        /* Get value of the buffer */

        val = BufVal(stp, tb, j, k);

        /* Get decay constant for this group */

        lambda = RDB[ptr + j];

        /* Calculate amount that survives next interval */

        Wtal = val*exp(-lambda*dt);

        /* Calculate amount to emit on next interval */

        Wemit = val*(1-exp(-lambda*dt));        

        if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_MESH)
          {
            /* Get pointer to value pair list */

            vplist = (long)RDB[loc0 + PRECDET_PTR_MESH_LIST];

            /* We can store this value pair anywhere on the list    */
            /* as long as we store each value pair into a different */
            /* position */

            /* Let's store zero-valued bins to the end of the list */
            /* so that when we sample bins based on this list, the */
            /* nonzero bins are looped over first */
            /* This will save time in samplemeshdelnu.c */

            if (Wemit == 0)
              {

                /* Get pointer from the end of list */

                vplist = ListPtr(vplist, ng*n0*n1*n2 - 1 - zeroes);

                /* Increment number of stored zero values */

                zeroes++;

              }
            else
              {

                /* Get pointer from the beginning of list */

                vplist = ListPtr(vplist, nonzeroes);

                /* Increment number of stored nonzero values */

                nonzeroes++;

              }

            /* Put index to bin */

            WDB[vplist + VALUE_PAIR_VAL1] = (double)(j*n0*n1*n2 + k);

            /* Put value to bin */

            WDB[vplist + VALUE_PAIR_VAL2] = Wemit;
          }        

        /* Check maximum */

        if (Wmax < Wemit)         
          Wmax = Wemit;

        /* Re-open buffer for writing */
                  
        WDB[DATA_BUF_REDUCED] = (double)NO;

        /* Tally remaining value to next time step */

        AddBuf(Wtal, 1.0, stp, 0, -1, tb + 1, j, k);

        /* Add weight to emit to total counter */

        Wall = Wall + Wemit;
      }

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_MESH)
    {
      /* Store maximum bin value */
      /* Used as a majorant in sampling mesh bin to emit */

      WDB[loc0 + PRECDET_MAX_WGT] = Wmax;

      /* Store weight emitted from spatial mesh */

      WDB[loc0 + PRECDET_W_EMIT] = Wall/RDB[DATA_NORM_COEF_N];

      /* Check that we stored every bin emission value to vp-list */

      if (nonzeroes + zeroes != ng*n0*n1*n2)
        Die(FUNCTION_NAME, "Stored %ld emission values, should have stored %ld", 
            nonzeroes + zeroes, ng*n0*n1*n2);

      /* Get pointer to value pair list */

      ptr = (long)RDB[loc0 + PRECDET_PTR_MESH_LIST];

      /* Sort value pair list for sampling */
      /* This might take very long with big meshes */
      /* So it is commented for now */

      /*SortList(ptr, VALUE_PAIR_VAL2, SORT_MODE_DESCEND);*/

#ifdef DNPRINT            
      /* Print out weight to emit over next time interval */

      fprintf(outp, "Weight to emit: %f (%E neutrons)\n", Wall/RDB[DATA_NORM_COEF_N], 
              Wall*RDB[DATA_NORM_COEF_N]);

      /* Print some things */

      fprintf(outp, "Weight of live neutrons    %f\n", 1.0);
#endif
    }

#ifdef DNPRINT
  fprintf(outp, "<-- decaymeshprecdet.c\n\n");  
#endif

}
