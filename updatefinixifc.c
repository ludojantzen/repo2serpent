/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writefinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Writes new values from FINIX to IFC                          */
/*                                                                           */
/* Comments:  - Also calculates convergence criterion                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "UpdateFinixIFC:"

/*****************************************************************************/

void UpdateFinixIFC()
{
  long axi, fib, fpe, ifc, i, j, k, ang, tptr, crptr, hrptr, dfptr;
  double f, rc, rc0, rh, rh0, T;
  Rod *rod;
  Results *results;
  Options *options;

  fib = (long)RDB[DATA_PTR_FIN0];

  /* Get pointer to interface */

  ifc = (long)RDB[fib + FINIX_PTR_IFC];

  /* Pointer to FUEP BLOCK */

  fpe = (long)RDB[ifc + IFC_PTR_FUEP];

  /* Loop over pins */

  while (fpe > VALID_PTR)
    {

      /* Get pointer to FINIX block based on fpe */

      fib = FinixPtrFromFpe(fpe);

      /* Get Finix pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);

      /* Get axial bin pointer */

      axi = (long)RDB[fpe + IFC_FUEP_PTR_AX];
      CheckPointer(FUNCTION_NAME, "(axi)", DATA_ARRAY, axi);

      /* Reset axial index */

      i = 0;

      /* Loop over axial zones to update temperatures and hot radii */

      while (axi > VALID_PTR)
        {

          if (RDB[axi + IFC_FUEP_AX_N_ANG] > 1)
            Die(FUNCTION_NAME, "Angular zones not supported with FINIX.");

          ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];

          /* Tän seuraavan vois periaatteessa korvata memcpylla jos ei ois */
          /* noita R2:ia. Vielä simppelimpi systeemi olisi, että IFC ja FIN*/
          /* pointterit osoittaisi samaan muistialueeseen (yksiköt)*/

          tptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP];
          crptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
          hrptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
          dfptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF];

          /* Reset radial index*/

          k = 0;

          /* Handle central hole */

          if (rod->pellet_inner_radius != 0.0)
            {
              /* Get temperature */

              T = results->temperature[i][0];

              /* Get radial node positions */

              rc = results->radial_node_position_cold[i][0];
              rh = results->radial_node_position[i][0];

              /* Radius is zero */

              WDB[crptr + 0] = 0.0;
              WDB[hrptr + 0] = 0.0;

              /* Calculate density factor of central hole       */
              /* (Should be dependent on the other gas volumes) */

              f = (rc*rc)/(rh*rh);

              /* Cutoff for larger than unity density factors */

              if (f > 1.0)
                f = 1.0;

              /* Store density factor */

              WDB[dfptr + 0] =  f;

              /* Update centerline temperatures */

              WDB[tptr  + 0] = T;

              /* Increment radial node count */

              k = 1;
            }

          /* Loop over remaining radial nodes */

          for (j = 0; j < options->pellet_radial_nodes +
                 options->clad_radial_nodes ; j++)
            {

              /* Get temperature */

              T = results->temperature[i][j];

              /* Get radial node positions */

              rc = results->radial_node_position_cold[i][j];
              rh = results->radial_node_position[i][j];

              /* Get radial node positions for prev node */

              if (j > 0)
                {
                  rc0 = results->radial_node_position_cold[i][j-1];
                  rh0 = results->radial_node_position[i][j-1];
                }
              else
                {
                  rc0 = 0.0;
                  rh0 = 0.0;
                }

              /* Write new cold radius */

              WDB[crptr + k] = (rc*100)*(rc*100);

              /* Write new hot radius */

              WDB[hrptr + k] = (rh*100)*(rh*100);

              /* Calculate density factor */
              /* Handle pellet centerline densityfactor */

              if (rh > 0)
                f = (rc*rc - rc0*rc0)/(rh*rh - rh0*rh0);
              else
                f = 1.0;

              /* Cutoff for larger than unity density factors   */
              /* TODO: Majorant density can be now increased in */
              /* mat-card -> Cut-off should only come to play   */
              /* if majorant density would be exceeded          */

              if (f > 1.0)
                f = 1.0;

              /* Store density factor */

              WDB[dfptr + k] = f;

              /* Update temperature */

              WDB[tptr  + k] = T;

              /* Increment radial node count */

              k++;

            }

          /* Next axial segment */

          axi = NextItem(axi);
          i++;

        }

      /* Next rod */

      fpe = NextItem(fpe);
    }
}

#endif

/*****************************************************************************/
