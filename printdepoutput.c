/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printdepoutput.c                               */
/*                                                                           */
/* Created:       2011/05/24 (JLe)                                           */
/* Last modified: 2017/07/30 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Prints output for depletion calculation                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintDepOutput:"

/*****************************************************************************/

void PrintDepOutput()
{
  long nbu, nnuc, nmat, n, m, i, j, k, l, sz, ptr;
  char name[MAX_STR];
  double **val, **tot, *bu, *days, *vol, *mbu, *totbu, *V, adens;
  FILE *fp, *fo;

  struct depnuc *nuc;

  /* Check mpi task */

  if (mpiid > 0)
    return;
  
  /***************************************************************************/

  /***** Read nuclide data ***************************************************/

  /* Open file for reading */

  sprintf(name, "%s.dep", GetText(DATA_PTR_INPUT_FNAME));
      
  if ((fp = fopen(name, "r")) == NULL)
    {
      /* Check re-deplete mode */

      if ((long)RDB[DATA_PARTICLE_REDEPLETE_MODE] == YES)
        Error(0, "Depletion file \"%s\" does not exist", name);
      else 
        Die(FUNCTION_NAME, "Unable to open depletion file for reading");
    }

  /* Open file for writing */
      
  sprintf(name, "%s_dep.m", GetText(DATA_PTR_INPUT_FNAME));

  if ((fo = fopen(name, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open depletion output for writing");

  /* Seek to end of file */

  fseek(fp, -sizeof(long), SEEK_END);
  
  /* Get number of burnup steps */

  if (fread(&nbu, sizeof(long), 1, fp) == 0)
    Die(FUNCTION_NAME, "fread error");

  nbu++;

  /* Check */

  CheckValue(FUNCTION_NAME, "nbu", "", nbu, 1, 1E+5);

  /* Rewind */

  rewind(fp);

  /* Read number of materials */
          
  if (fread(&nmat, sizeof(long), 1, fp) == 0)
    Die(FUNCTION_NAME, "fread error");

  CheckValue(FUNCTION_NAME, "nmat", "", nmat, 1, 1E+7);

  /* Read number of nuclides */
          
  if (fread(&nnuc, sizeof(long), 1, fp)  == 0)
    Die(FUNCTION_NAME, "fread error");

  CheckValue(FUNCTION_NAME, "nnuc", "", nnuc, 1, 1E+4);

  /* Allocate memory for nuclide data */

  nuc = (struct depnuc *)Mem(MEM_ALLOC, nnuc + 1, sizeof(struct depnuc));

  /* Read data */

  for (n = 0; n < nnuc; n++)
    {
      /* Read data */

      if (fread(&nuc[n + 1].ZAI, sizeof(long), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].Z, sizeof(long), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].AW, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].lambda, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].dh, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].sf, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].gI, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].ingtox, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");

      if (fread(&nuc[n + 1].inhtox, sizeof(double), 1, fp) == 0)
        Die(FUNCTION_NAME, "fread error");
    }

  /***************************************************************************/

  /***** Allocate memory for arrays ******************************************/

  /* Allocate memory for data */

  bu = (double *)Mem(MEM_ALLOC, nbu, sizeof(double));
  days = (double *)Mem(MEM_ALLOC, nbu, sizeof(double));
  vol = (double *)Mem(MEM_ALLOC, nbu, sizeof(double));
  mbu = (double *)Mem(MEM_ALLOC, nbu, sizeof(double));
  totbu = (double *)Mem(MEM_ALLOC, nbu, sizeof(double));
  V = (double *)Mem(MEM_ALLOC, nbu, sizeof(double));
  val = (double **)Mem(MEM_ALLOC, nbu, sizeof(double *));
  tot = (double **)Mem(MEM_ALLOC, nbu, sizeof(double *));
 
  for (n = 0; n < nbu; n++)
    {
      val[n] = (double *)Mem(MEM_ALLOC, nnuc + 1, sizeof(double));
      tot[n] = (double *)Mem(MEM_ALLOC, nnuc + 1, sizeof(double));
    }

  /**************************************************************************/

  /***** Print common data **************************************************/
  
  /* Print ZAI */

  fprintf(fo, "\nZAI = [\n");

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      fprintf(fo, "%ld\n", (long)RDB[ptr + INVENTORY_ZAI]);

      /* Next nuclide in list */
      
      ptr = NextItem(ptr);
    }

  fprintf(fo, "666\n");
  fprintf(fo, "0\n");
  fprintf(fo, "];\n");

  /* Print names */

  fprintf(fo, "\nNAMES = [\n");

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      fprintf(fo, "'%-16s'\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */
      
      ptr = NextItem(ptr);
    }

  fprintf(fo, "'%-16s'\n", "lost");
  fprintf(fo, "'%-16s'\n", "total");
  fprintf(fo, "];\n\n");

  /* Indexes */

  n = 1;
  
  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      fprintf(fo, "i%-6ld = %ld; i%-16s = %ld;\n", 
              (long)RDB[ptr + INVENTORY_ZAI], n,
              GetText(ptr + INVENTORY_PTR_NAME), n);
      
      /* Update index */

      n++;

      /* Next nuclide in list */
      
      ptr = NextItem(ptr);
    }

  fprintf(fo, "iLOST    = %ld;\n", n);
  fprintf(fo, "iTOT     = %ld;\n", n + 1);

  fprintf(fp, "\n");

  /**************************************************************************/

  /***** Material-wise values ***********************************************/

  /* Loop over materials */

  for (m = 0; m < nmat; m++)
    {
      /* Reset data */

      for (n = 0; n < nbu; n++)
        memset(val[n], 0.0, (nnuc + 1)*sizeof(double));
      
      memset(bu, 0.0, nbu*sizeof(double));
      memset(days, 0.0, nbu*sizeof(double));
      memset(vol, 0.0, nbu*sizeof(double));
      memset(mbu, 0.0, nbu*sizeof(double));

      /* Seek to beginning of data */
      
      sz = 2*sizeof(long) + nnuc*sizeof(struct depnuc);
      fseek(fp, sz, SEEK_SET);

      /* Loop over burnup steps */
      
      for (i = 0; i < nbu; i++)
        {
          /* Loop over materials */

          while (1 != 2)
            {
              /* Read burnup and burn time */
          
              if (fread(&bu[i], sizeof(double), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              if (fread(&totbu[i], sizeof(double), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              if (fread(&days[i], sizeof(double), 1, fp) == 0)         
                Die(FUNCTION_NAME, "fread error");
              
              /* Read number of nuclides */

              if (fread(&sz, sizeof(long), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              /* Read volume */

              if (fread(&vol[i], sizeof(double), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              /* Read burnup */

              if (fread(&mbu[i], sizeof(double), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              /* Read name */

              if (fread(&n, sizeof(long), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              if (fread(&name, sizeof(char), n, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              /* Put eof */

              name[n] = '\0';

              /* Read index and burnup step */

              if (fread(&n, sizeof(long), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              if (fread(&l, sizeof(long), 1, fp) == 0)
                Die(FUNCTION_NAME, "fread error");

              /* Compare */

              if ((n == m) && (l == i))
                {
                  for (k = 0; k < sz; k++)
                    {
                      /* Read index and atomic density */
                      
                      if (fread(&j, sizeof(long), 1, fp) == 0)
                        Die(FUNCTION_NAME, "fread error");

                      if (fread(&adens, sizeof(double), 1, fp) == 0)
                        Die(FUNCTION_NAME, "fread error");

                      /* Check values */
                      
                      CheckValue(FUNCTION_NAME, "j", "", j, 0, nnuc);
                      CheckValue(FUNCTION_NAME, "adens", "", adens, ZERO, 
                                 INFTY);
                      
                      /* Add to composition */
                      
                      val[i][j] = val[i][j] + adens;
                      tot[i][j] = tot[i][j] + adens*vol[i];
                    }
                  
                  /* Add to total volume */
                  
                  V[i] = V[i] + vol[i];
                  
                  /* Read burnup step (not used) */

                  if (fread(&k, sizeof(long), 1, fp) == 0)
                    Die(FUNCTION_NAME, "fread error");

                  /* Break loop */
                  
                  break;
                }          
              else
                {
                  /* Seek to the end of material */
                  
                  sz = sz*(sizeof(long) + sizeof(double)) + sizeof(long);
                  fseek(fp, sz, SEEK_CUR);
                }
            }
        }

      /* Print material-wise data */

      PrintDepVals(fo, name, nuc, nnuc + 1, nbu, val, vol, mbu);
    }

  /***************************************************************************/

  /***** More printing *******************************************************/

  /* Print totals */

  PrintDepVals(fo, NULL, nuc, nnuc + 1, nbu, tot, V, totbu);

  /* Print burnup vector */
  
  fprintf(fo, "\nBU = [ ");
  
  for (i = 0; i < nbu; i++)
    fprintf(fo, "%1.5E ", bu[i]);
  
  fprintf(fo, "];\n");
  
  /* Print day vector */
  
  fprintf(fo, "\nDAYS = [ ");
  
  for (i = 0; i < nbu; i++)
    fprintf(fo, "%1.5E ", days[i]);
  
  fprintf(fo, "];\n");

  /***************************************************************************/
  
  /* Free memory */
  
  for (n = 0; n < nbu; n++)
    {
      Mem(MEM_FREE, val[n]);
      Mem(MEM_FREE, tot[n]);
    }
  
  Mem(MEM_FREE, nuc);
  Mem(MEM_FREE, val);
  Mem(MEM_FREE, vol);
  Mem(MEM_FREE, mbu);
  Mem(MEM_FREE, totbu);
  Mem(MEM_FREE, V);
  Mem(MEM_FREE, tot);
  Mem(MEM_FREE, bu);
  Mem(MEM_FREE, days);
  
  /* Close files */
  
  fclose(fp);
  fclose(fo);
}
  
/*****************************************************************************/
