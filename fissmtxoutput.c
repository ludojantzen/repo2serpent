/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fissmtxoutput.c                                */
/*                                                                           */
/* Created:       2012/09/02 (JLe)                                           */
/* Last modified: 2017/05/04 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Writes (analog) fission matrixes into file                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FissMtxOutput:"

/*****************************************************************************/

void FissMtxOutput()
{
  long loc0, mat, uni, sz, ptr, n, m;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check pointer */

  if ((loc0 = (long)RDB[DATA_PTR_FMTX]) < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Open file for writing */
  
  sprintf(tmpstr, "%s_fmtx%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);
  
  if ((fp = fopen(tmpstr, "w")) == NULL)
    Warn(FUNCTION_NAME, "Unable to open file for writing");

  fprintf(fp, "\n%% ----- Analog fission matrixes\n\n");

  /* Check type */

  if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_MAT)
    {
      /* Print material indexes */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Print name and index */
          
          if ((n = (long)RDB[mat + MATERIAL_FMTX_IDX]) > -1)
            fprintf(fp, "fmtx_mat (%3ld, [1 : %2ld]) = '%s' ;\n", 
                    n + 1, strlen(GetText(mat + MATERIAL_PTR_NAME)),
                    GetText(mat + MATERIAL_PTR_NAME));
          
          /* Next material */
          
          mat = NextItem(mat);
        }
      
      /* Newline */
      
      fprintf(fp, "\n");
    }
  else if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_UNI)
    {
      /* Print universe indexes */

      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Print name and index */
          
          if ((n = (long)RDB[uni + UNIVERSE_FMTX_IDX]) > -1)
            fprintf(fp, "fmtx_uni (%3ld, [1 : %2ld]) = '%s' ;\n", 
                    n + 1, strlen(GetText(uni + UNIVERSE_PTR_NAME)),
                    GetText(uni + UNIVERSE_PTR_NAME));
          
          /* Next universe */
          
          uni = NextItem(uni);
        }
      
      /* Newline */
      
      fprintf(fp, "\n");
    }
  else if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_XYZ)
    {
      /* Print mesh data */

      ptr = (long)RDB[loc0 + FMTX_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "fmtx_xmin = %1.5E ;\n", RDB[ptr + MESH_MIN0]);
      fprintf(fp, "fmtx_xmax = %1.5E ;\n", RDB[ptr + MESH_MAX0]);
      fprintf(fp, "fmtx_nx = %ld ;\n", (long)RDB[ptr + MESH_N0]);

      fprintf(fp, "fmtx_ymin = %1.5E ;\n", RDB[ptr + MESH_MIN1]);
      fprintf(fp, "fmtx_ymax = %1.5E ;\n", RDB[ptr + MESH_MAX1]);
      fprintf(fp, "fmtx_ny = %ld ;\n", (long)RDB[ptr + MESH_N1]);

      fprintf(fp, "fmtx_zmin = %1.5E ;\n", RDB[ptr + MESH_MIN2]);
      fprintf(fp, "fmtx_zmax = %1.5E ;\n", RDB[ptr + MESH_MAX2]);
      fprintf(fp, "fmtx_nz = %ld ;\n", (long)RDB[ptr + MESH_N2]);
      
      /* Newline */
      
      fprintf(fp, "\n");
    }
  else
    Die(FUNCTION_NAME, "Invalid matrix type");
     
  /* Get size */
      
  sz = (long)RDB[loc0 + FMTX_SIZE];
  
  /* Pointer to stats */
  
  ptr = (long)RDB[loc0 + FMTX_PTR_MTX];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
  /* Print total */

  fprintf(fp, "fmtx_t = zeros(%ld,%ld);\n", sz, sz);      
  fprintf(fp, "fmtx_t_err = zeros(%ld,%ld);\n\n", sz, sz);      

  for (n = 0; n < sz; n++)
    for (m = 0; m < sz; m++)
      /*if ((Mean(ptr, 0, n, m) != 0.0) && (RelErr(ptr, 0, n, m) < 0.1))*/
      if ((Mean(ptr, 0, n, m) != 0.0))
        {
          fprintf(fp, "fmtx_t (%3ld, %3ld) = %1.5E ;  ", m + 1, 
                  n + 1, Mean(ptr, 0, n, m));
          fprintf(fp, "fmtx_t_err (%3ld, %3ld) = %1.5E ;\n", 
                  m + 1, n + 1, RelErr(ptr, 0, n, m));
        }
  
  /* Newline */
  
  fprintf(fp, "\n");                      
  
#ifdef mmmmmmmmmmmmmmmmm

  /* Print prompt */

  fprintf(fp, "fmtx_p = zeros(%ld,%ld);\n", sz, sz);
  fprintf(fp, "fmtx_p_err = zeros(%ld,%ld);\n\n", sz, sz);      
  
  for (n = 0; n < sz; n++)
    for (m = 0; m < sz; m++)
      if (Mean(ptr, 1, n, m) != 0.0)
        {
          fprintf(fp, "fmtx_p (%3ld, %3ld) = %1.5E ;  ", m + 1, 
                  n + 1, Mean(ptr, 1, n, m));
          fprintf(fp, "fmtx_p_err (%3ld, %3ld) = %1.5E ;\n", 
                  m + 1, n + 1, RelErr(ptr, 1, n, m));
        }
  
  /* Newline */
  
  fprintf(fp, "\n");                      
  
  /* Print delayed */
  
  fprintf(fp, "fmtx_d = zeros(%ld,%ld);\n", sz, sz);      
  fprintf(fp, "fmtx_d_err = zeros(%ld,%ld);\n\n", sz, sz);      

  for (n = 0; n < sz; n++)
    for (m = 0; m < sz; m++)
      if (Mean(ptr, 2, n, m) != 0.0)
        {
          fprintf(fp, "fmtx_d (%3ld, %3ld) = %1.5E ;  ", m + 1, 
                  n + 1, Mean(ptr, 2, n, m));
          fprintf(fp, "fmtx_d_err (%3ld, %3ld) = %1.5E ;\n", 
                  m + 1, n + 1, RelErr(ptr, 2, n, m));
        }
  
  /* Newline */
  
  fprintf(fp, "\n");                      

#endif
  
  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
