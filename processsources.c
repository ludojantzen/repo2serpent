/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsources.c                               */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2018/11/27 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes source definitions                                 */
/*                                                                           */
/* Comments: - Source cells and universes are first processed in             */
/*             creategeometry.c. The pointer search is repeated here to get  */
/*             STL geometries working.                                       */
/*                                                                           */
/*           - Ton radioactive decay sourcen kanssa ei pitäisi sallia muita  */
/*             fotonilähteitä, muuten normeeraus menee pieleen.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSources:"

/*****************************************************************************/

void ProcessSources()
{
  long src, surf, mat, loc0, loc1, loc2, loc3, n, ptr, cell, ne, i;
  char *name, *str, tmpstr[MAX_STR];
  double tot, max;
  FILE *fp;

  /***************************************************************************/

  /***** Source for RIA calculation ******************************************/

  if ((long)RDB[DATA_PTR_RIA0] > VALID_PTR)
    {
      /* Allocate memory for source structure (not included in list) */

      src = ReallocMem(DATA_ARRAY, SRC_BLOCK_SIZE);

      /* Put source weight and type */

      WDB[src + SRC_WGT] = 1.0;
      WDB[src + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;

      /* Reset boundaries */
      
      WDB[src + SRC_XMIN] = -INFTY;
      WDB[src + SRC_XMAX] =  INFTY;
      WDB[src + SRC_YMIN] = -INFTY;
      WDB[src + SRC_YMAX] =  INFTY;
      WDB[src + SRC_ZMIN] = -INFTY;
      WDB[src + SRC_ZMAX] =  INFTY;
      
      /* Reset point */
      
      WDB[src + SRC_X0] = -INFTY;
      WDB[src + SRC_Y0] = -INFTY;
      WDB[src + SRC_Z0] = -INFTY;
      
      /* Reset energy */

      WDB[src + SRC_E] = -INFTY;
      
      /* Put file name */

      sprintf(tmpstr, "%s.src", GetText(DATA_PTR_INPUT_FNAME));
      WDB[src + SRC_READ_PTR_FILE] = (double)PutText(tmpstr);

      /* Allocate memory */

      ptr = ReallocMem(DATA_ARRAY, 
                       (long)RDB[DATA_SRC_FILE_BUF_SIZE]*SRC_BUF_BLOCK_SIZE);
      
      /* Put pointer and buffer size */
      
      WDB[src + SRC_READ_PTR_BUF] = (double)ptr;
      WDB[src + SRC_READ_BUF_SZ] = RDB[DATA_SRC_FILE_BUF_SIZE];
      
      /* Reset file and buffer indexes */
      
      WDB[src + SRC_READ_FILE_POS] = 0.0;
      WDB[src + SRC_READ_BUF_IDX] = RDB[DATA_SRC_FILE_BUF_SIZE];
      
      /* Put type */

      WDB[src + SRC_READ_FILE_TYPE] = (double)SRC_FILE_TYPE_S1_RENORM;

      /* Put pointer */

      WDB[DATA_PTR_RIA_SRC] = (double)src;
    }

  /***************************************************************************/

  /***** Source for calculation with precursors ******************************/

  /* Check if there is a source input */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
    if ((long)RDB[loc0 + PRECDET_PTR_IN_FNAME] > VALID_PTR)
      {
        /*******************************************/
        /* Create initial precursor neutron source */
        /*******************************************/

        /* Allocate memory for source structure (not included in list) */

        src = ReallocMem(DATA_ARRAY, SRC_BLOCK_SIZE);

        /* Reset weight */

        WDB[src + SRC_WGT] = 1.0;
              
        /* Reset boundaries */

        WDB[src + SRC_XMIN] = -INFTY;
        WDB[src + SRC_XMAX] =  INFTY;
        WDB[src + SRC_YMIN] = -INFTY;
        WDB[src + SRC_YMAX] =  INFTY;
        WDB[src + SRC_ZMIN] = -INFTY;
        WDB[src + SRC_ZMAX] =  INFTY;

        /* Reset point */

        WDB[src + SRC_X0] = -INFTY;
        WDB[src + SRC_Y0] = -INFTY;
        WDB[src + SRC_Z0] = -INFTY;

        /* Set default type to neutron */

        WDB[src + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;

        /* Reset energy */

        WDB[src + SRC_E] = -INFTY;

        /* Source name */
          
        WDB[src + SRC_PTR_NAME] = (double)PutText("DynsrcPrec");

        /* Print precursor filename to tmpstr */

        sprintf(tmpstr, "%s.precpoints", GetText(loc0 + PRECDET_PTR_IN_FNAME));

        /* Set file name */
              
        WDB[src + SRC_READ_PTR_FILE] = (double)PutText(tmpstr);
              
        /* Allocate memory */

        ptr = ReallocMem(DATA_ARRAY, (long)RDB[DATA_SRC_FILE_BUF_SIZE]
                         *SRC_BUF_BLOCK_SIZE);
      
        /* Put pointer and buffer size */
      
        WDB[src + SRC_READ_PTR_BUF] = (double)ptr;
        WDB[src + SRC_READ_BUF_SZ] = RDB[DATA_SRC_FILE_BUF_SIZE];
      
        /* Reset file and buffer indexes */
      
        WDB[src + SRC_READ_FILE_POS] = 0.0;
        WDB[src + SRC_READ_BUF_IDX] = RDB[DATA_SRC_FILE_BUF_SIZE];
      
        /* Put type */

        WDB[src + SRC_READ_FILE_TYPE] = (double)SRC_FILE_TYPE_SERPENT1;

        /* Set binary input on */

        WDB[src + SRC_READ_BINARY] = (double)YES;

        /* Put pointer */

        WDB[loc0 + PRECDET_PTR_PREC_SRC] = (double)src;
      }

  /***************************************************************************/

  /***** Process sources *****************************************************/

  /* Check that source definition exists in external source mode */
  
  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) &&
      ((long)RDB[DATA_PTR_SRC0] < VALID_PTR))
    Error(0, "External source mode without source definition");

  /* Check that source definitions exist */

  if ((long)RDB[DATA_PTR_SRC0] < VALID_PTR)
    return;
  
  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];  
  while (src > VALID_PTR)
    {
      /* Check type and reaction mode */

      if (((long)RDB[src + SRC_TYPE] == PARTICLE_TYPE_GAMMA) &&
          ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT))
        Error(src, "Gamma source %s in criticality souce simulation",
              GetText(src + SRC_PTR_NAME));

      /* Check if type not set */

      if ((long)RDB[src + SRC_TYPE] == 0)
        {
          /* Check combined mode */

          if ((long)RDB[DATA_MULTI_PARTICLE_TRANSPORT] == YES)
            Error(src, "Source type must be set in multi-particle mode");
          else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
            WDB[src + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;
          else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
            WDB[src + SRC_TYPE] = (double)PARTICLE_TYPE_GAMMA;
        }
      else if (((long)RDB[DATA_PHOTON_PRODUCTION] > 0) && 
               ((long)RDB[src + SRC_TYPE] != PARTICLE_TYPE_NEUTRON))
        {
          /* Type must be neutron */

          Error(src, "Only neutron sources allowed in multi-particle mode");
        }

      /* Check for multiple types (NOTE: useampi tyyppi tarvitaan jos      */
      /* esim. halutaan yhdistää radioaktiivisuuslähde kytkettyyn moodiin  */
      /* Tämä kuitenkin edellyttää esim. koko normeerauksen kirjoittamista */
      /* uusiksi. */

      if ((ptr = PrevItem(src)) > VALID_PTR)
        if ((long)RDB[src + SRC_TYPE] != (long)RDB[ptr + SRC_TYPE])
          Error(src, "Multiple source particle types not allowed");

      /***********************************************************************/
      
      /***** Link materials to sources ***************************************/
      
      /* Check pointer */
      
      if ((long)RDB[src + SRC_PTR_MAT] > VALID_PTR)
        {
          /* Find material */

          mat = (long)RDB[DATA_PTR_M0];
          if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
                                 GetText(src + SRC_PTR_MAT))) > VALID_PTR)
            {
              /* Set pointer */
                      
              WDB[src + SRC_PTR_MAT] = (double)mat;

              /* Set used-flag */

              SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
            }
          else
            Error(src, "Material %s in source %s not defined", 
                  GetText(src + SRC_PTR_MAT), 
                  GetText(src + SRC_PTR_NAME));
        }
      
      /***********************************************************************/

      /***** Link cells to sources *******************************************/

      /* Check pointer */
      
      if ((long)RDB[src + SRC_PTR_CELL] > VALID_PTR)
        {
          /* Find cell */

          cell = (long)RDB[DATA_PTR_C0];
          if ((cell = SeekListStr(cell, CELL_PTR_NAME, 
                                 GetText(src + SRC_PTR_CELL))) > VALID_PTR)
            {
              /* Set pointer */
                      
              WDB[src + SRC_PTR_CELL] = (double)cell;

              /* Set used-flag */

              SetOption(cell + CELL_OPTIONS, OPT_USED);
            }
          else
            Error(src, "Cell %s in source %s not defined", 
                  GetText(src + SRC_PTR_CELL), 
                  GetText(src + SRC_PTR_NAME));
        }

      /***********************************************************************/

      /***** Link radioactive materials to sources ***************************/
      
      /* Check pointer */
      
      if ((long)RDB[src + SRC_PTR_RAD_SRC_MAT] > VALID_PTR)
        {
          /* Check that decay data library is set (NOTE: this does not */
          /* guaranee that the data is read?) */

          if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] < VALID_PTR)
            Error(src, "Decay library must be defined for radiation source");
            
          /* Check special */

          if (!strcmp(GetText(src + SRC_PTR_RAD_SRC_MAT), "-1"))
            WDB[src + SRC_PTR_RAD_SRC_MAT] = -1.0;
          else
            {
              /* Find material */

              mat = (long)RDB[DATA_PTR_M0];
              if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
                                     GetText(src + SRC_PTR_RAD_SRC_MAT))) 
                  > VALID_PTR)
                {
                  /* Set pointer */
                  
                  WDB[src + SRC_PTR_RAD_SRC_MAT] = (double)mat;

                  /* Set used-flag */

                  SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
                }
              else
                Error(src, "Material %s in source %s not defined", 
                      GetText(src + SRC_PTR_RAD_SRC_MAT), 
                      GetText(src + SRC_PTR_NAME));
            }

          /* Put global pointer for normalization (tarkista ettei oo eri */
          /* materiaali, tms? */

          if ((long)RDB[DATA_NORM_PTR_RAD_SRC_MAT] != 0)
            Warn(FUNCTION_NAME, 
                 "Multiple radioactive decay sources not allowed");
          else
            WDB[DATA_NORM_PTR_RAD_SRC_MAT] = RDB[src + SRC_PTR_RAD_SRC_MAT];  
        }
      
      /***********************************************************************/
      
      /***** Link surfaces to sources ****************************************/

      /* Check if surface is given */

      if ((long)RDB[src + SRC_PTR_SURF] > VALID_PTR)
        {
          /* Get name */
          
          str = GetText(src + SRC_PTR_SURF);
          
          if (*str == '-')
            name = &str[1];
          else
            name = str;
          
          /* Find surface */

          surf = (long)RDB[DATA_PTR_S0];
          if ((surf = SeekListStr(surf, SURFACE_PTR_NAME, name)) > VALID_PTR)
            {
              /* Set pointer */
                      
              WDB[src + SRC_PTR_SURF] = (double)surf;

              /* Set used-flag */

              SetOption(surf + SURFACE_OPTIONS, OPT_USED);
            }
          else
            Error(src, "Surface %s in source %s not defined", name,
                  GetText(src + SRC_PTR_NAME));

          /* Set direction */
          
          if (*str == '-')
            WDB[src + SRC_SURF_SIDE] = -1.0;
          else
            WDB[src + SRC_SURF_SIDE] = 1.0;
        }

      /***********************************************************************/

      /***** Source files ****************************************************/

      if ((long)RDB[src + SRC_READ_PTR_FILE] > VALID_PTR)
        {
          /* Open file */

          if ((fp = fopen(GetText(src + SRC_READ_PTR_FILE), "r")) == NULL)
            Error(src, "Unable to open source file \"%s\"", 
                  GetText(src + SRC_READ_PTR_FILE));
          else
            fclose(fp);

          /* Test file format */

          if (!((long)RDB[src + SRC_READ_BINARY]))
            TestDOSFile(GetText(src + SRC_READ_PTR_FILE));

          /* Allocate memory */

          ptr = ReallocMem(DATA_ARRAY, (long)RDB[DATA_SRC_FILE_BUF_SIZE]
                           *SRC_BUF_BLOCK_SIZE);

          /* Put pointer and buffer size */

          WDB[src + SRC_READ_PTR_BUF] = (double)ptr;
          WDB[src + SRC_READ_BUF_SZ] = RDB[DATA_SRC_FILE_BUF_SIZE];

          /* Reset file and buffer indexes */

          WDB[src + SRC_READ_FILE_POS] = 0.0;
          WDB[src + SRC_READ_BUF_IDX] = RDB[DATA_SRC_FILE_BUF_SIZE];
        }

      /***********************************************************************/

      /***** Fusion plasma sources *******************************************/

      /* Check type */

      if ((long)RDB[src + SRC_READ_FILE_TYPE] == SRC_FILE_TYPE_FUSION_PLASMA)
        {
          /* Read and process plasma source */

          ReadPlasmaSrc(src);
        }

      /***********************************************************************/

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Renormalize directions **********************************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Calculate total */

      tot = sqrt(RDB[src + SRC_U0]*RDB[src + SRC_U0] + 
                 RDB[src + SRC_V0]*RDB[src + SRC_V0] + 
                 RDB[src + SRC_W0]*RDB[src + SRC_W0]);

      /* Normalize vectors */

      if (tot > 0.0)
        {
          WDB[src + SRC_U0] = RDB[src + SRC_U0]/tot;
          WDB[src + SRC_V0] = RDB[src + SRC_V0]/tot;
          WDB[src + SRC_W0] = RDB[src + SRC_W0]/tot;
        }
      else
        {
          /* Set x-component to infinity to indicate isotropic */

          WDB[src + SRC_U0] = INFTY;
        }

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Renormalize weights *************************************************/

  /* Reset total */

  tot = 0.0;

  /* Loop over sources to calculate total */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Add to sum */

      tot = tot + RDB[src + SRC_WGT];

      /* Next source */

      src = NextItem(src);
    }

  /* Check sum */

  if (tot == 0.0)
    Error(0, "Sum of source weights is zero");

  /* Loop over sources to normalize */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Divide by total */

      WDB[src + SRC_WGT] = RDB[src + SRC_WGT]/tot;

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Check energy bin order and normalize ********************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Pointer to energies bins */

      if ((loc0 = (long)RDB[src + SRC_PTR_EBINS]) > VALID_PTR)
        {
          /* Get Number of points */

          ne = (long)RDB[loc0 + SRC_EBIN_NE];
          CheckValue(FUNCTION_NAME, "(ne)", "", ne, 2, 10000000);

          /* Get interpolation */

          i = (long)RDB[loc0 + SRC_EBIN_INTT];
          CheckValue(FUNCTION_NAME, "i", "", i, 0, 4);

          /* Pointer to energy grid and pdf */

          loc1 = RDB[loc0 + SRC_EBIN_PTR_E];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          loc2 = RDB[loc0 + SRC_EBIN_PTR_PDF];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Check first point */

          if (RDB[loc2] != 0.0)
            Error(src, "Probability of first energy point must be zero");

          /* Allocate memory for cdf */

          loc3 = ReallocMem(DATA_ARRAY, ne);
          WDB[loc0 + SRC_EBIN_PTR_CDF] = (double)loc3;

          /* Integrate */

          for (n = 1; n < ne; n++)
            {
              /* Check energies */

              if (RDB[loc1 + n] < RDB[loc1 + n - 1])
                Error(src, "Source energies must be in ascending order");

              /* Check interpolation */

              if (i == 0)
                WDB[loc3 + n] = RDB[loc3 + n - 1] + RDB[loc2 + n];
              else if (i == 1)
                WDB[loc3 + n] = RDB[loc3 + n - 1] + 
                  (RDB[loc1 + n] - RDB[loc1 + n - 1])*RDB[loc2 + n];
              else if (i == 2)
                WDB[loc3 + n] = RDB[loc3 + n - 1] + 
                  (RDB[loc1 + n] - RDB[loc1 + n - 1])*
                  (RDB[loc2 + n] + RDB[loc2 + n - 1])*0.5;
              else
                Error(src, "Unsupported interpolation type");
            }

          /* Normalize */

          if ((max = RDB[loc3 + ne - 1]) <= 0.0)
            Error(src, "Error in energy distribution");

          CheckValue(FUNCTION_NAME, "max", "", max, ZERO, INFTY);

          for (n = 0; n < ne; n++)
            {
              WDB[loc2 + n] = RDB[loc2 + n]/max;
              WDB[loc3 + n] = RDB[loc3 + n]/max;
            }
        }

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/
}  

/*****************************************************************************/
