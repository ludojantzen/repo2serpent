/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : expodecfit.c                                   */
/*                                                                           */
/* Created:       2016/09/08 (JLe)                                           */
/* Last modified: 2016/09/08 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Fits exponential functions to decay heat curve.              */
/*                                                                           */
/* Comments: - Tää voi feilata jos decay data librarya ei lueta.             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ExpoDecFit:"

/*****************************************************************************/

void ExpoDecFit()
{
  long nf, nt, n, mat, nuc, iso;
  double tmin, tmax, vol, adens, lambda, A, H, *t, *tgt, *param;

  /* Get number of exponential functions and check if fit is defined */

  if ((nf = (long)RDB[DATA_EXPO_DEC_FIT_NF]) < 1)
    return;

  /***************************************************************************/
  
  /***** Form target function ************************************************/

  /* Number of time points and limits */

  nt = (long)RDB[DATA_EXPO_DEC_FIT_NT];
  CheckValue(FUNCTION_NAME, "nt", "", nt, 2, 100000);

  tmin = RDB[DATA_EXPO_DEC_FIT_TMIN];
  CheckValue(FUNCTION_NAME, "tmin", "", tmin, 0, INFTY);

  tmax = RDB[DATA_EXPO_DEC_FIT_TMAX];
  CheckValue(FUNCTION_NAME, "tmax", "", tmax, tmin, INFTY);

  /* Create time array (linear) */

  t = MakeArray(tmin, tmax, nt, 1);
  
  /* Allocate memory for target array */

  tgt = Mem(MEM_ALLOC, nt, sizeof(double));

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get volume */
      
      vol = RDB[mat + MATERIAL_VOLUME];

      /* Check physical flag and volume */
      
      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT) ||
          (vol <= 0.0) || (vol >= 1E+18))
        {
          /* Pointer to next */
          
          mat = NextItem(mat);
          
          /* Cycle loop */
          
          continue;
        }
      
      /* Check divisor type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        Die(FUNCTION_NAME, "Divided parent material");
      
      /* Loop over composition */
          
      if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) < VALID_PTR)
        iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */
          
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
          
          /* Get atomic density and decay constant */
          
          adens = RDB[iso + COMPOSITION_ADENS];
          lambda = RDB[nuc + NUCLIDE_LAMBDA];

          /* Calculate activity and decay heat */
          
          A = vol*adens*lambda*1E+24;
          H = A*RDB[nuc + NUCLIDE_DECAY_E]*MEV;

          /* Check values */
              
          CheckValue(FUNCTION_NAME, "A", "", A, 0.0, 1E+30);
          CheckValue(FUNCTION_NAME, "H", "", H, 0.0, 1E+30);
         
          /* Loop over time points and add to target array*/

          for (n = 0; n < nt; n++)
            tgt[n] = tgt[n] + H*exp(-lambda*t[n]);
 
          /* Next isotope */
              
          iso = NextItem(iso);
        }              

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Make exponential fits ***********************************************/

  /* Allocate memory for fitting parameters */

  param = Mem(MEM_ALLOC, 2*nf, sizeof(double));

  /* Printataan kohdefunktio */

  for (n = 0; n < nt; n++)
    printf("%E %E\n", t[n], tgt[n]);

  /***************************************************************************/

  /* Free temporary arrays */

  Mem(MEM_FREE, t);
  Mem(MEM_FREE, tgt);
  Mem(MEM_FREE, param);


  /***************************************************************************/
}

/*****************************************************************************/
