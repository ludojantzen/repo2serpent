/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readdirectoryfile.c                            */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2016/10/04 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description:   Reads cross section directory files into ACE data          */
/*                structure                                                  */
/*                                                                           */
/* Comments: - Converted from Serpent 1.1.14 readdirectoryfiles.c            */
/*           - List functions are not used because data is read in ACE       */
/*             block.                                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadDirectoryFile:"

/*****************************************************************************/

void ReadDirectoryFile()
{
  long n, ptr, ZA, Z, A, I, type;
  long pta, ace;
  double AW, T;
  char tmpstr[MAX_STR], chk[MAX_STR], alias[MAX_STR], name[MAX_STR];
  char path[MAX_STR];
  FILE *fp;

  /* Check pointer */
  
  if ((pta = (long)RDB[DATA_PTR_ACEDATA_FNAME_LIST]) < 1)
    return;

  fprintf(outp, "Reading ACE directory files...\n");
  
  /***************************************************************************/

  /***** Read transport, dosimetry and thermal scattering data ***************/

  /* Reset data pointer */

  ace = -1;

  /* Loop over files */

  while ((long)RDB[pta] > 0)
    {
      /* Test format */

      WDB[DATA_DUMMY] = RDB[pta];
      TestDOSFile(GetText(DATA_DUMMY));

      /* Open file for reading */

      fp = OpenDataFile(pta, "cross section directory");
      
      /* Loop over directory file */
      
      while ((n = fscanf(fp, "%s", alias)) > 0)
        {          
          /* Check for MCNP xsdir file format */

          sprintf(chk, "%s", alias);
          chk[8] = '\0';
          if (!strcasecmp(chk, "datapath"))
            Error(0, "acelib entry points to MCNP format xsdir file");

          /* Reset EOF flag */
          
          n = 0;
          
          /* Name */
          
          if (fscanf(fp, "%s", name) < 1)
            break;
          
          /* Read type */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            type = AtoI(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
          
          /* Get ZA */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            ZA = AtoI(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
          
          /* Check value */
          
          if ((type != XS_TYPE_SAB) && ((ZA < 1000) || (ZA > 120000)))
            Die(FUNCTION_NAME, "Invalid ZA %ld", ZA);
          
          /* Separate Z and A */
          
          Z = (long)((double)ZA/1000.0);
          A = ZA - 1000*Z;

          /* Check A */

          if (A > 300)
            Note(0, "Nuclide %s has incorrect ZA (%ld) in directory file", 
                 name, ZA);

          /* Get isomeric state */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            I = AtoI(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
  
          /* Check value */

          CheckValue(FUNCTION_NAME, "I", "", I, 0, 2);
              
          /* Atomic weight */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            AW = AtoF(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
          
          /* Temperature */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            T = AtoF(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;

          /* Check value */

          CheckValue(FUNCTION_NAME, "T", "", T, 0.0, 10000.0);
          
          /* Binary format flag */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            {
              /* Check that type is 0 */
              
              if (AtoI(tmpstr, FUNCTION_NAME, NULL, -1) != 0)
                Die(FUNCTION_NAME, "Data format %s not supported by %s %s",
                    tmpstr, CODE_NAME, CODE_VERSION);
            }
          else
            break;
          
          /* File path */
          
          if (fscanf(fp, "%s", path) < 1)
            break;

          /* Check type */
          
          if ((type < 1) || (type > 7) || (type == 4) || (type == 6))
            Die(FUNCTION_NAME, "Invalid data type %ld", type);
          else if (type != 7)
            {
              /* Allocate memory for data */

              ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);

              /* Set null pointer to next */
              
              ACE[ace + ACE_PTR_NEXT] = NULLPTR;
              
              /* Check if previous exists (ei voi k‰ytt‰‰ VALID_PTR) */
              
              if ((ptr = (long)RDB[DATA_PTR_ACE0]) < 0)
                {
                  /* First definition, set pointer */
                  
                  WDB[DATA_PTR_ACE0] = (double)ace;
                }
              else
                {
                  /* Find last block  (tohon ei VALID_PTR) */
                  
                  while ((long)ACE[ptr + ACE_PTR_NEXT] > 0)
                    ptr = (long)ACE[ptr + ACE_PTR_NEXT];
                  
                  /* Set pointer to new */
                  
                  ACE[ptr + ACE_PTR_NEXT] = (double)ace;
                }

              /* Put alias, name and library id */
          
              ACE[ace + ACE_PTR_ALIAS] = (double)PutText(alias);
              ACE[ace + ACE_PTR_NAME] = (double)PutText(name);
              ACE[ace + ACE_PTR_LIB_ID] = 
                (double)PutText(&name[strlen(name) - 3]);
              
              /* Put type */
              
              ACE[ace + ACE_TYPE] = (double)type;
          
              /* Put ZA, I and ZAI */

              ACE[ace + ACE_ZA] = (double)ZA;
              ACE[ace + ACE_I] = (double)I;
              ACE[ace + ACE_ZAI] = (double)I + (double)ZA*10.0;

              /* Put atomic weight and temperature */
          
              ACE[ace + ACE_AW] = AW;
              ACE[ace + ACE_TEMP] = T;


              /* Put file path */
          
              ACE[ace + ACE_PTR_FILE] = (double)PutText(path);
            }
          
          /*******************************************************************/
        }
      
      /* Check error flag */
      
      if (n == 0)
        Die(FUNCTION_NAME, "Error in ace directory file %s", GetText(pta));
      
      /* Close file */
      
      fclose(fp);

      /* Next file */

      pta++;
    }

  /***************************************************************************/

  /***** Transmutation reaction data *****************************************/

  /* Loop over files */
  
  pta = (long)RDB[DATA_PTR_ACEDATA_FNAME_LIST];
  while ((long)RDB[pta] > 0)
    {
      /* Open file for reading */

      fp = OpenDataFile(pta, "cross section directory");
      
      /* Loop over directory file */
      
      while ((n = fscanf(fp, "%s", alias)) > 0)
        {          
          /* Check for MCNP xsdir file format */

          sprintf(chk, "%s", alias);
          chk[8] = '\0';
          if (!strcasecmp(chk, "datapath"))
            Error(0, "acelib entry points to MCNP format xsdir file");

          /* Reset EOF flag */
          
          n = 0;
          
          /* Name */
          
          if (fscanf(fp, "%s", name) < 1)
            break;
          
          /* Read type */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            type = AtoI(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
          
          /* Get ZA */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            ZA = AtoI(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
          
          /* Check value */
          
          if ((type != XS_TYPE_SAB) && ((ZA < 1000) || (ZA > 120000)))
            Die(FUNCTION_NAME, "Invalid ZA %ld", ZA);
          
          /* Get isomeric state */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            I = AtoI(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
  
          /* Check value */

          CheckValue(FUNCTION_NAME, "I", "", I, 0, 2);
              
          /* Atomic weight */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            AW = AtoF(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;
          
          /* Temperature */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            T = AtoF(tmpstr, FUNCTION_NAME, NULL, -1);
          else
            break;

          /* Check value */

          CheckValue(FUNCTION_NAME, "T", "", T, 0.0, 10000.0);
          
          /* Binary format flag */
          
          if (fscanf(fp, "%s", tmpstr) > 0)
            {
              /* Check that type is 0 */
              
              if (AtoI(tmpstr, FUNCTION_NAME, NULL, -1) != 0)
                Die(FUNCTION_NAME, "Data format %s not supported by %s %s",
                    tmpstr, CODE_NAME, CODE_VERSION);
            }
          else
            break;
          
          /* File path */
          
          if (fscanf(fp, "%s", path) < 1)
            break;

          /* Check type */
          
          if ((type < 1) || (type > 7) || (type == 4) || (type == 6))
            Die(FUNCTION_NAME, "Invalid data type %ld", type);
          else if (type == 7)
            {
              /* Allocate memory for data */

              ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);

              /* Set null pointer to next */
              
              ACE[ace + ACE_PTR_NEXT] = NULLPTR;
              
              /* Check if previous exists (ei voi k‰ytt‰‰ VALID_PTR) */
              
              if ((ptr = (long)RDB[DATA_PTR_ACE0]) < 1)
                Die(FUNCTION_NAME, "No transport data");
              else
                {
                  /* Find last block  (tohon ei VALID_PTR) */
                  
                  while ((long)ACE[ptr + ACE_PTR_NEXT] > 0)
                    ptr = (long)ACE[ptr + ACE_PTR_NEXT];
                  
                  /* Set pointer to new */
                  
                  ACE[ptr + ACE_PTR_NEXT] = (double)ace;
                }

              /* Put pointer to first */

              if ((long)RDB[DATA_PTR_TRANSMU_ACE0] < 1)
                WDB[DATA_PTR_TRANSMU_ACE0] = (double)ace;

              /* Put alias, name and library id */
          
              ACE[ace + ACE_PTR_ALIAS] = (double)PutText(alias);
              ACE[ace + ACE_PTR_NAME] = (double)PutText(name);
              ACE[ace + ACE_PTR_LIB_ID] = 
                (double)PutText(&name[strlen(name) - 3]);
              
              /* Put type */
              
              ACE[ace + ACE_TYPE] = (double)type;
          
              /* Put ZA, I and ZAI */

              ACE[ace + ACE_ZA] = (double)ZA;
              ACE[ace + ACE_I] = (double)I;
              ACE[ace + ACE_ZAI] = (double)I + (double)ZA*10.0;

              /* Put atomic weight and temperature */
          
              ACE[ace + ACE_AW] = AW;
              ACE[ace + ACE_TEMP] = T;


              /* Put file path */
          
              ACE[ace + ACE_PTR_FILE] = (double)PutText(path);
            }
          
          /*******************************************************************/
        }
      
      /* Check error flag */
      
      if (n == 0)
        Die(FUNCTION_NAME, "Error in ace directory file %s", GetText(pta));
      
      /* Close file */
      
      fclose(fp);

      /* Next file */

      pta++;
    }

  /***************************************************************************/

  /* Exit */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
