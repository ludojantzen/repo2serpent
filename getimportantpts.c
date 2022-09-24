/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getimportantpts.c                              */
/*                                                                           */
/* Created:       2010/12/13 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Reads local minima and maxima from criss section data to an  */
/*              garray                                                       */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetImportantPts:"

/*****************************************************************************/

/* NOTE: kynnysenergiat menee pieleen ellei datassa ole nollaa */

double *GetImportantPts(long ace, long *NXS, long *JXS, long *ni)
{
  long NES, n0, i, itot, ptr, xss, n, nr, mt, ne, loc0, loc1, add;
  double *tmp, *Ei;

  /* Get maximum number of energy points */
                  
  NES = NXS[2];

  /* Allocate memory for data */

  tmp = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

  /* Reset values */

  for (n = 0; n < NES; n++)
    tmp[n] = 0.0;

  /* Reset total number of important points */

  itot = 0;

  /* Pointer to XSS array */
          
  xss = (long)ACE[ace + ACE_PTR_XSS];
        
  /* Loop over reaction channels */
          
  for (nr = -1; nr < NXS[3]; nr++)
    {
      /* Reset pointers to energy array and xs data */

      loc0 = -1;
      loc1 = -1;

      /* Reset indexes */

      n0 = -1;
      ne = -1;

      if (nr == -1)
        {
          /* Elastic scattering. Set mt */

          mt = 2;

          /* Set number of energy points */
          
          ne = NES;
          
          /* Set index to first point */

          n0 = 0;

          /* Pointer to energy array */

          loc0 = xss + JXS[0] - 1;

          /* Get pointer to xs data */

          loc1 = xss + JXS[0] - 1 + 3*NES;
        }
      else
        {
          /* Other reactions, get mt */

          ptr = xss + JXS[2] - 1 + nr;
          mt = (long)ACE[ptr];

          /* Check mt and specials flag (tää on vähän tarpeeton) */
          
          if (((mt > 15) && (mt < 100)) || ((mt > 101) && (mt < 200)) ||
              ((long)RDB[DATA_OPTI_INCLUDE_SPECIALS] == YES))
            {
              /* Get pointer to SIG-block (Table F-10, page F-17) */
              
              ptr = xss + JXS[5] - 1 + nr;
              ptr = xss + (long)ACE[ptr] + JXS[6] - 1;

              /* Get number of energy points */
                      
              ne = (long)ACE[ptr];              
              
              /* Get index to first point */
                      
              n0 = NXS[2] - (long)ACE[ptr];

              /* Pointer to energy array */

              loc0 = xss + JXS[0] - 1 + n0;
              
              /* Get pointer to xs data */
              
              loc1 = ptr + 1;
            }
        }

      /* Find important points */

      if (loc0 > 0)
        {
          /* Check indexes */

          CheckValue(FUNCTION_NAME, "n0", "", n0, 0, NES);
          CheckValue(FUNCTION_NAME, "ne", "", n0, 0, NES);
          CheckValue(FUNCTION_NAME, "n0 + ne", "", n0 + ne, 0, NES);

          /* Reset number of important points */

          i = 0;

          /* Loop over grid */

          for (n = 0; n < ne; n++)
            {
              /* Reset add flag */

              add = NO;

              /* First two points */

              if (n < 2)
                add = YES;

              /* Last two points */

              else if (n > ne - 3)
                add = YES;

              /* Local minimum */

              else if ((ACE[loc1 + n] < ACE[loc1 + n - 1]) &&
                       (ACE[loc1 + n] < ACE[loc1 + n + 1]))
                add = YES;

              /* Local maximum */

              else if ((ACE[loc1 + n] > ACE[loc1 + n - 1]) &&
                       (ACE[loc1 + n] > ACE[loc1 + n + 1]))
                add = YES;

              /* Check add-flag */

              if (add == YES)
                {
                  /* Add new or compare to existing value */

                  if (tmp[n0 + n] == 0.0)
                    {
                      tmp[n0 + n] = ACE[loc0 + n];
                      itot++;
                    }
                  else if (tmp[n0 + n] != ACE[loc0 + n])
                    {
                      fprintf(errp, "%s Mismatch in energy data.\n", 
                              FUNCTION_NAME);
                      exit(-1);
                    }

                  i++;
                }
            }
        }
    }
  
  /* Allocate memory for final important data */

  Ei = (double *)Mem(MEM_ALLOC, itot, sizeof(double));

  /* Clean data (remove zeros) */
  
  i = 0;

  for (n = 0; n < NES; n++)
    if (tmp[n] > 0.0)
      Ei[i++] = tmp[n];

  /* Check pointer */

  if (i != itot)
    {
      fprintf(errp, "%s i != itot.\n", FUNCTION_NAME);
      exit(-1);
    }

  /* Free temporary array */

  if (tmp != NULL)
    Mem(MEM_FREE, tmp);

  /* Set array size */

  *ni = itot;

  /* Return pointer */

  return Ei;
}

/*****************************************************************************/
