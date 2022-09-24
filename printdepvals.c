/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printdepvals.c                                 */
/*                                                                           */
/* Created:       2011/05/24 (JLe)                                           */
/* Last modified: 2015/06/17 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Formatted output for PrintDepOutput()                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintDepVals:"

/*****************************************************************************/

void PrintDepVals(FILE *fp, char *mat, struct depnuc *nuc, long nnuc, 
                  long nbu, double **val, double *vol, double *mbu)
{
  long i, j,ptr, ZAI;
  double sum;
  char name[MAX_STR];

  /* Name for output */
  
  if (mat != NULL)
    sprintf(name, "MAT_%s", mat);
  else
    sprintf(name, "TOT");
  
  /* Volume */

  fprintf(fp, "\n%s_VOLUME = [ ", name);

  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", vol[i]);

  fprintf(fp, "];\n");

  /* Burnup */

  if (mbu != NULL)
    {
      fprintf(fp, "\n%s_BURNUP = [ ", name);
      
      for (i = 0; i < nbu; i++)
        fprintf(fp, "%1.5E ", mbu[i]);
      
      fprintf(fp, "];\n");
    }

  /***************************************************************************/

  /***** Atomic density ******************************************************/

  fprintf(fp, "\n%s_ADENS = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) || (nuc[j].Z == ZAI))
              sum = sum + val[i][j];
             
          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", sum);
          else if (vol[i] > 0.0)
            fprintf(fp, "%1.5E ", sum/vol[i]);
          else
            fprintf(fp, "%1.5E ", 0.0);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    {
      if (mat != NULL)
        fprintf(fp, "%1.5E ", val[i][0]);
      else if (vol[i] > 0.0)
        fprintf(fp, "%1.5E ", val[i][0]/vol[i]);
      else
        fprintf(fp, "%1.5E ", 0.0);
    }

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j];

      /* Print value */

      if (mat != NULL)
        fprintf(fp, "%1.5E ", sum);
      else if (vol[i] > 0.0)
        fprintf(fp, "%1.5E ", sum/vol[i]);
      else
        fprintf(fp, "%1.5E ", 0.0);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/

  /***** Mass density ********************************************************/

  if (mat != NULL)
    fprintf(fp, "\n%s_MDENS = [\n", name);
  else
    fprintf(fp, "\n%s_MASS = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */
      
      ZAI = (long)RDB[ptr + INVENTORY_ZAI];
          
      /* Loop over burnup steps */
      
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;
          
          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].AW/N_AVOGADRO;
          
          /* Print value */
        
          fprintf(fp, "%1.5E ", sum);
        }
      
      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));
      
      /* Next nuclide in list */
      
      ptr = NextItem(ptr);
    }
  
  /* Lost data */
  
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);
  
  fprintf(fp, "%% lost data\n");
      
  /* Total */
  
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].AW/N_AVOGADRO;
      
      /* Print value */
        
      fprintf(fp, "%1.5E ", sum);  
    }
      
  fprintf(fp, "%% total\n");
  
  fprintf(fp, "];\n");
  
  /***************************************************************************/

  /***** Activity ************************************************************/

  fprintf(fp, "\n%s_A = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].lambda;

          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
          else
            fprintf(fp, "%1.5E ", sum/BARN);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].lambda;

      /* Print value */
      
      if (mat != NULL)
        fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
      else
        fprintf(fp, "%1.5E ", sum/BARN);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/

  /***** Decay heat **********************************************************/

  fprintf(fp, "\n%s_H = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].dh;

          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", vol[i]*sum/BARN*MEV);
          else
            fprintf(fp, "%1.5E ", sum/BARN*MEV);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].dh;

      /* Print value */

      if (mat != NULL)
        fprintf(fp, "%1.5E ", vol[i]*sum/BARN*MEV);
      else
        fprintf(fp, "%1.5E ", sum/BARN*MEV);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/

  /***** Spontaneous fission rate ********************************************/

  fprintf(fp, "\n%s_SF = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].sf;

          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
          else
            fprintf(fp, "%1.5E ", sum/BARN);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].sf;

      /* Print value */

      if (mat != NULL)
        fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
      else
        fprintf(fp, "%1.5E ", sum/BARN);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/

  /***** Photon source rate **************************************************/

  fprintf(fp, "\n%s_GSRC = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].gI;

          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
          else
            fprintf(fp, "%1.5E ", sum/BARN);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].gI;

      /* Print value */

      if (mat != NULL)
        fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
      else
        fprintf(fp, "%1.5E ", sum/BARN);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/

  /***** Ingestion toxicity **************************************************/

  fprintf(fp, "\n%s_ING_TOX = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].lambda*nuc[j].ingtox;

          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
          else
            fprintf(fp, "%1.5E ", sum/BARN);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].lambda*nuc[j].ingtox;

      /* Print value */

      if (mat != NULL)
        fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
      else
        fprintf(fp, "%1.5E ", sum/BARN);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/

  /***** Inhalation toxicity *************************************************/

  fprintf(fp, "\n%s_INH_TOX = [\n", name);

  /* Nuclides in inventory list */

  ptr = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (ptr > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[ptr + INVENTORY_ZAI];

      /* Loop over burnup steps */
          
      for (i = 0; i < nbu; i++)
        {
          /* Reset sum */
          
          sum = 0.0;

          /* Loop over nuclides and add to sum */
          
          for (j = 1; j < nnuc; j++)
            if ((nuc[j].ZAI == ZAI) ||(nuc[j].Z == ZAI))
              sum = sum + val[i][j]*nuc[j].lambda*nuc[j].inhtox;

          /* Print value */

          if (mat != NULL)
            fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
          else
            fprintf(fp, "%1.5E ", sum/BARN);
        }

      fprintf(fp, "%% %s\n", GetText(ptr + INVENTORY_PTR_NAME));

      /* Next nuclide in list */

      ptr = NextItem(ptr);
    }

  /* Lost data */
          
  for (i = 0; i < nbu; i++)
    fprintf(fp, "%1.5E ", 0.0);

  fprintf(fp, "%% lost data\n");

  /* Total */
        
  for (i = 0; i < nbu; i++)
    {
      /* Reset sum */
      
      sum = 0.0;
      
      /* Loop over nuclides and add to sum */
      
      for (j = 1; j < nnuc; j++)
        sum = sum + val[i][j]*nuc[j].lambda*nuc[j].inhtox;

      /* Print value */

      if (mat != NULL)
        fprintf(fp, "%1.5E ", vol[i]*sum/BARN);
      else
        fprintf(fp, "%1.5E ", sum/BARN);
    }

  fprintf(fp, "%% total\n");

  fprintf(fp, "];\n");

  /***************************************************************************/
}

/*****************************************************************************/
