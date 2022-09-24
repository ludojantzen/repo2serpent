/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printsrcimp.c                                  */
/*                                                                           */
/* Created:       2017/01/28 (JLe)                                           */
/* Last modified: 2017/06/07 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Prints maximum source importances for weight window mesh     */
/*              normalization                                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintSrcImp:"

/*****************************************************************************/

void PrintSrcImp()
{
  long wwd, i, n, imax, nmax, ne, nb, msh, ptr, loc0, erg;
  double max, val, x, y, z, E;

  /* Check calculation */

  if ((long)RDB[DATA_SRC_IMP_CALC] == NO)
    return;

  /* Pointer to weight window structure */
  
  if ((wwd = (long)RDB[DATA_PTR_WWD0]) < VALID_PTR)
    return;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Loop over weight window structures */

  while (wwd > VALID_PTR)
    {
      /* Reset maximum importance and indexex */

      max = -INFTY;
      imax = -1;
      nmax = -1;

      /* Get number of energy groups */

      if ((long)RDB[wwd + WWD_PTR_ERG] > VALID_PTR)
        ne = (long)RDB[wwd + WWD_NE];
      else 
        ne = 1;
      
      /* Pointer to mesh */
      
      msh = (long)RDB[wwd + WWD_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get number of spatial bins */
      
      nb = (long)(RDB[msh + MESH_N0]*RDB[msh + MESH_N1]*RDB[msh + MESH_N2]);

      /* Pointer to statistics */

      ptr = (long)RDB[wwd + WWD_PTR_SRC_IMP_DIS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over data */

      for (i = 0; i < nb; i++)
        for (n = 0; n < ne; n++)
          if ((val = BufVal(ptr, n, i)) > max)
            {
              /* Store */

              max = val;
              imax = i;
              nmax = n;
            }

      /* Check indexes */

      if ((imax < 0) || (nmax < 0))
        Die(FUNCTION_NAME, "fail");

      /* Check mesh type */

      if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
        Die(FUNCTION_NAME, "Error in mesh type");

      /* Get pointer to data */

      ptr = (long)RDB[msh + MESH_PTR_PTR] + imax;
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Avoid compiler warning */

      x = 0.0;
      y = 0.0;
      z = 0.0;

      /* Find point in cell (tää on tosi typerä tapa) */

      for (n = 0; n < 1000000; n++)
        {
          /* Sample coordinates */

          x = RandF(0)*(RDB[msh + MESH_MAX0] - RDB[msh + MESH_MIN0])
            + RDB[msh + MESH_MIN0]; 
          y = RandF(0)*(RDB[msh + MESH_MAX1] - RDB[msh + MESH_MIN1])
            + RDB[msh + MESH_MIN1]; 
          z = RandF(0)*(RDB[msh + MESH_MAX2] - RDB[msh + MESH_MIN2])
            + RDB[msh + MESH_MIN2]; 

          /* Check */

          if (ptr == MeshPtr(msh, x, y, z))
            break;
        }

      if (n == 1000000)
        {
          fprintf(outp, "\nFailed. Try again.\n\n");
          exit(-1);
        }

      /* Pointer to structure */

      loc0 = (long)RDB[ptr];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Check type */

      if ((long)RDB[wwd + WWD_TYPE] == WWD_MESH_TYPE_MCNP)
        {
          /* Pointer to importance vector */

          ptr = (long)RDB[loc0 + WWD_MESH_PTR_IMP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Get importance */
          
          val = RDB[ptr + nmax]*RDB[wwd + WWD_NORM_FACT];
          
          /* Pointer to energy array */

          erg = (long)RDB[wwd + WWD_PTR_ERG];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
          
          /* Calculate mean */

          E = 0.5*(RDB[erg + nmax] + RDB[erg + nmax + 1]);
        }
      else
        {
          val = RDB[loc0 + WWD_MESH_IMP]*RDB[wwd + WWD_NORM_FACT];
          E = 1.0;
        }

      fprintf(outp, "\nTry:\n\n");

      fprintf(outp, "wn 1.0 %12.5E %12.5E %12.5E %12.5E\n", x, y, z, E);

      fprintf(outp, "\nfor weight window mesh normalization\n");

      /* Next */

      wwd = NextItem(wwd);
    }

  /* Exit OK */

  fprintf(outp, "\nOK.\n\n");

  /* Terminate run */

  exit(-1);
}

/*****************************************************************************/

#ifdef mmmmmmmmmmmmmmmm


  /* Check neutron importance */

  if (RDB[DATA_NEUTRON_MAX_SRC_IMP] > 0.0)
    {
      fprintf(outp, "\nMaximum neutron importance %1.5E at:\n\n",
              RDB[DATA_NEUTRON_MAX_SRC_IMP]);
      fprintf(outp, "x: %12.5E cm\n", RDB[DATA_NEUTRON_MAX_SRC_X]);
      fprintf(outp, "y: %12.5E cm\n", RDB[DATA_NEUTRON_MAX_SRC_Y]);
      fprintf(outp, "z: %12.5E cm\n", RDB[DATA_NEUTRON_MAX_SRC_Z]);
      fprintf(outp, "E: %12.5E MeV\n\n", RDB[DATA_NEUTRON_MAX_SRC_E]);

      fprintf(outp, "For normalization use:\n\n");

      fprintf(outp, "wn 1.0 %12.5E %12.5E %12.5E %12.5E\n",
              RDB[DATA_NEUTRON_MAX_SRC_X],
              RDB[DATA_NEUTRON_MAX_SRC_Y],
              RDB[DATA_NEUTRON_MAX_SRC_Z],
              RDB[DATA_NEUTRON_MAX_SRC_E]);
    }

  /* Check photon importance */

  if (RDB[DATA_PHOTON_MAX_SRC_IMP] > 0.0)
    {
      fprintf(outp, "\nMaximum photon importance %1.5E at:\n\n",
              RDB[DATA_PHOTON_MAX_SRC_IMP]);
      fprintf(outp, "x: %12.5E cm\n", RDB[DATA_PHOTON_MAX_SRC_X]);
      fprintf(outp, "y: %12.5E cm\n", RDB[DATA_PHOTON_MAX_SRC_Y]);
      fprintf(outp, "z: %12.5E cm\n", RDB[DATA_PHOTON_MAX_SRC_Z]);
      fprintf(outp, "E: %12.5E MeV\n\n", RDB[DATA_PHOTON_MAX_SRC_E]);

      fprintf(outp, "For normalization use:\n\n");

      fprintf(outp, "wn 1.0 %12.5E %12.5E %12.5E %12.5E\n",
              RDB[DATA_PHOTON_MAX_SRC_X],
              RDB[DATA_PHOTON_MAX_SRC_Y],
              RDB[DATA_PHOTON_MAX_SRC_Z],
              RDB[DATA_PHOTON_MAX_SRC_E]);

#endif
