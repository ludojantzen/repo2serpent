/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectvrmeshdata.c                            */
/*                                                                           */
/* Created:       2012/04/19 (JLe)                                           */
/* Last modified: 2019/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects mesh data for variance reduction routines           */
/*                                                                           */
/* Comments: Tätä tarvitaan koska RES2-array, johon data kerätään,           */
/*           resetoidaan clearprivateres.c:ssä.                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectVRMeshData:"

/*****************************************************************************/

void CollectVRMeshData()
{
  long ptr, msh, n0, n1, n2, i, j, k;
  double f, avg;

  /***************************************************************************/

  /***** Source mesh for UFS *************************************************/

  if (((long)RDB[DATA_UFS_MODE] != UFS_MODE_NONE) &&
      ((long)RDB[DATA_UFS_MODE] != UFS_MODE_RMX))
    {
      /* Get pointer to mesh */

      msh = (long)RDB[DATA_UFS_PTR_SRC_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get mesh size */

      n0 = (long)RDB[msh + MESH_N0];
      n1 = (long)RDB[msh + MESH_N1];
      n2 = (long)RDB[msh + MESH_N2];

      /* Check values */

      CheckValue(FUNCTION_NAME, "n0", "", n0, 1, 10000);
      CheckValue(FUNCTION_NAME, "n1", "", n1, 1, 10000);
      CheckValue(FUNCTION_NAME, "n2", "", n2, 1, 10000);
      
      /* Loop over mesh and calculate number of non-zero cells */

      avg = 0.0;
      for (i = 0; i < n0; i++)
        for (j = 0; j < n1; j++)
          for (k = 0; k < n2; k++)
            if (ReadMesh(msh, i, j, k) > 0.0)
              avg = avg + 1.0;

      /* Calculate average (tätä kutsutaan turhaan ihan simulaation alussa */
      /* clearstats.c:stä, kun mitään dataa ei vielä ole kerätty */

      if (avg == 0.0)
        return;
      else
        avg = MeshTot(msh)/avg;
      
      /* Check value */

      CheckValue(FUNCTION_NAME, "avg", "", avg, ZERO, INFTY);

      /* Get pointer to data */

      ptr = (long)RDB[DATA_UFS_PTR_FACTORS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Reset data */

      memset(&WDB[ptr], 0.0, n0*n1*n2*sizeof(double));

      /* Loop over mesh */

      for (i = 0; i < n0; i++)
        for (j = 0; j < n1; j++)
          for (k = 0; k < n2; k++)
            {
              /* Get factor */
              
              f = ReadMesh(msh, i, j, k)/avg;
              
              /* Adjust */
          
              if (f > 0.0)
                f = pow(1.0/f, fabs(RDB[DATA_UFS_ORDER]));
              else
                f = 1.0;
              
              /* Store value */

              WDB[ptr + i + j*n0 + k*n0*n1] = f;
            }
    }
  
  /***************************************************************************/
}

/*****************************************************************************/
