/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processphotonatt.c                             */
/*                                                                           */
/* Created:       2015/05/26 (JLe)                                           */
/* Last modified: 2017/03/04 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Process photon attenuation data                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "photon_attenuation.h"

#define FUNCTION_NAME "ProcessPhotonAtt:"

/*****************************************************************************/

void ProcessPhotonAtt()
{
  long mat, iso, nuc, Z, ne, n, N, i0, i, ptr;
  double Ei[100], fi[100], *E, *f, mdens, mfrac, val;

  /* Avoid compiler warning (photon_attenuation.h sisältää dataa jota */
  /* käytetään useammassa aliohjelmassa) */

  n = idx1[0][0];
  val = dat1[0][0];

  /* Adjust coincident points */

  ne = (long)idx0[91][0] + (long)idx0[91][1];

  for (n = 1; n < ne; n++)
    if (dat0[n][1] == dat0[n - 1][1])
      dat0[n][1] = dat0[n][1] + 1E-11;

  fprintf(outp, "Processing response functions for photon dose rates...\n");

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /***********************************************************************/

      /***** Unionize energy grid ********************************************/

      /* Reset pointer */

      E = NULL;
      N = 0;

      /* Get mass density */

      if ((mdens = RDB[mat + MATERIAL_MDENS]) == 0.0)
        {
          /* Skip zero-density material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
       { 
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get Z */

          Z = (long)RDB[nuc + NUCLIDE_Z];

          /* Data only goes up to 92 (ZAI = -1 is lost) */

          if ((Z > 92) || ((long)RDB[nuc + NUCLIDE_ZAI] == -1))
            {
              /* Skip */

              iso = NextItem(iso);

              /* Cycle loop */

              continue;
            }

          /* Get index to data array and number of points */

          i0 = (long)idx0[Z - 1][0];
          ne = (long)idx0[Z - 1][1];

          /* Check first and last energy point */

          if (dat0[i0][1] != 1E-3)
            Die(FUNCTION_NAME, "Mismatch in first energy point");
          else if (dat0[i0 + ne - 1][1] != 20.0)
            Die(FUNCTION_NAME, "Mismatch in last energy point");

          /* Read data */

          for (n = 0; n < ne; n++)
            {
              /* Check Z and read value */

              if (dat0[i0 + n][0] != Z)
                Die(FUNCTION_NAME, "Mismatch in Z");
              else
                Ei[n] = dat0[i0 + n][1];
            }

          /* Add points to main array */

          E = AddPts(E, &N, Ei, ne);

          /* Next */

          iso = NextItem(iso);
        }

      /* Check order */

      for (n = 1; n < N; n++)
        if (E[n] <= E[n - 1])
          Die(FUNCTION_NAME, "Error in order");

      /* Skip the material if no data was found (TKa 9.2.2017 / 2.1.29) */

      if (N == 0) 
        {
          /* Pointer to next */
          
          mat = NextItem(mat);

          /* Cycle loop */
          
          continue;
        }

      /***********************************************************************/

      /***** Reconstruct data ************************************************/

      /* Put number of points */

      WDB[mat + MATERIAL_PHOTON_ATT_NE] = (double)N;
      
      /* Put energy array in material structure */

      ptr = ReallocMem(DATA_ARRAY, N);
      WDB[mat + MATERIAL_PTR_PHOTON_ATT_E] = (double)ptr;

      for (n = 0; n < N; n++)
        WDB[ptr++] = E[n];

      /* Allocate memory for data array */

      ptr = ReallocMem(DATA_ARRAY, N);
      WDB[mat + MATERIAL_PTR_PHOTON_ATT_F] = (double)ptr;
      f = &WDB[ptr];

      /* Reset vector (just to be safe) */

      memset(f, 0.0, N*sizeof(double));

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Calculate mass fraction */

          mfrac = RDB[nuc + NUCLIDE_AW]*RDB[iso + COMPOSITION_ADENS]/
            mdens/N_AVOGADRO;

          /* Get Z */

          Z = (long)RDB[nuc + NUCLIDE_Z];

          /* Data only goes up to 92 (ZAI = -1 is lost) */

          if ((Z > 92) || ((long)RDB[nuc + NUCLIDE_ZAI] == -1))
            {
              /* Skip */

              iso = NextItem(iso);

              /* Cycle loop */

              continue;
            }

          /* Get index to data array and number of points */

          i0 = (long)idx0[Z - 1][0];
          ne = (long)idx0[Z - 1][1];

          /* Read data */

          for (n = 0; n < ne; n++)
            {
              Ei[n] = dat0[i0 + n][1];
              fi[n] = dat0[i0 + n][3];
            }

          /* Loop over unionized grid */

          for (n = 0; n < N - 1; n++)
            {
              /* Find index */

              if ((i = SearchArray(Ei, E[n], ne)) < 0)
                Die(FUNCTION_NAME, "Point not in grid: %ld %E", n, E[n]);

              /* Interpolate */

              val = ENDFInterp(5, E[n], Ei[i], Ei[i + 1], fi[i], fi[i + 1]);

              /* Add to data */

              f[n] = f[n] + val*mfrac;
            }

          /* Add last point */

          f[N - 1] = f[N -1] + fi[ne - 1]*mfrac;

          /* Next */

          iso = NextItem(iso);
        }

      /***********************************************************************/

      /* Free memory */

      if (E != NULL)
        Mem(MEM_FREE, E);

      /* Next material */

      mat = NextItem(mat);
    }

  /* Exit OK */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
