/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processicm.c                                   */
/*                                                                           */
/* Created:       2013/10/01 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Process data needed for ICM                                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessICM:"

/*****************************************************************************/

void ProcessICM()
{
  long ptr, loc0, surf, lat, n, m, sz, ene, npmax, ng0, ng1, nseg, nsub;
  long nmu0, nmu1, nmu2, nmua, nmus;
  char tmpstr[MAX_STR];

  /* Check option */

  if((long)RDB[DATA_ICM_CALC] == NO)
    return;
  
  fprintf(outp, "Processing data for interface current method...\n");

  /***************************************************************************/

  /***** Link energy group data **********************************************/

  /* Check pointer to main grid */

  if ((long)RDB[DATA_ICM_PTR_ENE0] < VALID_PTR)
    Error(0, "Missing energy grid in ICM calculation");

  /* Get name */
      
  sprintf(tmpstr, "%s", GetText(DATA_ICM_PTR_ENE0));

  /* Find grid */

  ene = (long)RDB[DATA_PTR_ENE0];
  if ((ene = SeekListStr(ene, ENE_PTR_NAME, tmpstr)) > VALID_PTR)
    {
      /* Pointer to grid data */

      ene = (long)RDB[ene + ENE_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ene)", DATA_ARRAY, ene);
      
      /* Put pointer */
      
      WDB[DATA_ICM_PTR_ENE0] = (double)ene;

      /* Set grid size */

      WDB[DATA_ICM_NG0] = RDB[ene + ENERGY_GRID_NE] - 1.0;
    }
  else
    Error(0, "ICM energy grid %s is not defined", tmpstr);

  /* Check pointer to reconstruction grid */

  if ((long)RDB[DATA_ICM_PTR_ENE1] < VALID_PTR)
    Error(0, "Missing energy grid in ICM calculation");
  
  /* Get name */
      
  sprintf(tmpstr, "%s", GetText(DATA_ICM_PTR_ENE1));

  /* Find grid */

  ene = (long)RDB[DATA_PTR_ENE0];
  if ((ene = SeekListStr(ene, ENE_PTR_NAME, tmpstr)) > VALID_PTR)
    {  
      /* Pointer to grid data */

      ene = (long)RDB[ene + ENE_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ene)", DATA_ARRAY, ene);

      /* Put pointer */
      
      WDB[DATA_ICM_PTR_ENE1] = (double)ene;

      /* Set grid size */

      WDB[DATA_ICM_NG1] = RDB[ene + ENERGY_GRID_NE] - 1.0;
    }
  else if (!strcmp(tmpstr, "-1"))
    WDB[DATA_ICM_PTR_ENE1] = NULLPTR;
  else
    Error(0, "ICM energy grid %s is not defined", tmpstr);

  /***************************************************************************/
  
  /***** Link surfaces *******************************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_ICM0];
  while (loc0 > VALID_PTR)
    {
      /* Find surface */

      surf = (long)RDB[DATA_PTR_S0];
      if ((surf = SeekListStr(surf, SURFACE_PTR_NAME, 
                             GetText(loc0 + ICM_PTR_SURF))) < VALID_PTR)
        Error(loc0, "Surface %s is not defined", 
              GetText(loc0 + ICM_PTR_SURF));

      /* Put pointer */

      WDB[loc0 + ICM_PTR_SURF] = (double)surf;

      /* Get pointer to surface parameters */

      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get number of faces */

      switch ((long)RDB[surf + SURFACE_TYPE])
        {
        case SURF_SQC:
          {
            /* Set number of surfaces */
            
            n = 4;

            /* Set origin */

            WDB[loc0 + ICM_X0] = RDB[ptr];
            WDB[loc0 + ICM_Y0] = RDB[ptr + 1];
            WDB[loc0 + ICM_Z0] = 0.0;

            break;
          }
        default:
          Error(loc0, "Surface %s is wrong type for ICM calculation",
                GetText(surf + SURFACE_PTR_NAME));
        }
      
      /* Get number of sub-segments */

      m = (long)RDB[DATA_ICM_NSUB];
      CheckValue(FUNCTION_NAME, "m", "", m, 1, 50);
    
      /* Store values */
      
      if (((long)RDB[DATA_ICM_NSEG] > 0) && 
          ((long)RDB[DATA_ICM_NSEG] != n*m))
        Error(loc0, "Mismatch in surface type");
      else
        WDB[DATA_ICM_NSEG] = (double)(n*m);

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Stop tracks at outer boundary */

  WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

  /***************************************************************************/
  
  /***** Link lattices *******************************************************/

  /* Reset maximum number of pins */

  npmax = 0;

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_ICM0];
  while (loc0 > VALID_PTR)
    {
      /* Check pointer */

      if ((long)RDB[loc0 + ICM_PTR_LAT] < VALID_PTR)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Find lattice */

      lat = (long)RDB[DATA_PTR_L0];
      if ((lat = SeekListStr(lat, LAT_PTR_NAME, 
                             GetText(loc0 + ICM_PTR_LAT))) < VALID_PTR)
        Error(loc0, "Lattice %s is not defined", 
              GetText(loc0 + ICM_PTR_LAT));

      /* Put pointer */
      
      WDB[loc0 + ICM_PTR_LAT] = (double)lat;

      /* Get number of pins */
          
      if ((n = (long)RDB[lat + LAT_NTOT]) > 0)
        WDB[loc0 + ICM_NP] = (double)n;
      else
        Die(FUNCTION_NAME, "Number of pins is zero");
      
      /* Compare to maximum */
      
      if (n > npmax)
        npmax = n;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check discretization ************************************************/

  /* Get number of groups, segments and angular bins */
  
  if ((ng0 = (long)RDB[DATA_ICM_NG0]) < 1)
    Die(FUNCTION_NAME, "Invalid number of energy groups");

  if ((ng1 = (long)RDB[DATA_ICM_NG1]) < 0)
    Die(FUNCTION_NAME, "Invalid number of energy groups");
  
  if ((nseg = (long)RDB[DATA_ICM_NSEG]) < 4)
    Die(FUNCTION_NAME, "Invalid number of segments");

  nsub = (long)RDB[DATA_ICM_NSUB];
  CheckValue(FUNCTION_NAME, "nsub", "", nsub, 1, 50);

  nmu0 = (long)RDB[DATA_ICM_NMU0];
  CheckValue(FUNCTION_NAME, "nmu0", "", nmu0, 1, 50);

  nmu1 = (long)RDB[DATA_ICM_NMU1];
  CheckValue(FUNCTION_NAME, "nmu1", "", nmu1, 1, 50);

  nmu2 = (long)RDB[DATA_ICM_NMU2];
  CheckValue(FUNCTION_NAME, "nmu2", "", nmu2, 1, 50);

  /* Put sub-segments if not defined */
  
  if ((ptr = (long)RDB[DATA_ICM_PTR_SUB]) < VALID_PTR)
    {
      ptr = ReallocMem(DATA_ARRAY, nsub + 1);
      WDB[DATA_ICM_PTR_SUB] = (double)ptr;
      
      for (n = 0; n < nsub + 1; n++)
        WDB[ptr++] = ((double)n)/((double)nsub);

      WDB[DATA_ICM_NSUB] = (double)nsub;
    }
  else
    {
      /* Check order */

      for (n = 1; n < nsub + 1; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Error(0, "ICM sub-segmentation must be in ascending order");
    }

  /* Put angular bins if not defined */
  
  if ((ptr = (long)RDB[DATA_ICM_PTR_MU0]) < VALID_PTR)
    {
      ptr = ReallocMem(DATA_ARRAY, nmu0 + 1);
      WDB[DATA_ICM_PTR_MU0] = (double)ptr;
      
      for (n = 0; n < nmu0 + 1; n++)
        WDB[ptr++] = ((double)n)/((double)nmu0);

      WDB[DATA_ICM_NMU0] = (double)nmu0;
    }
  else
    {
      /* Check order */

      for (n = 1; n < nmu0 + 1; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Error(0, "ICM angular binning must be in ascending order");
    }

  if ((ptr = (long)RDB[DATA_ICM_PTR_MU1]) < VALID_PTR)
    {
      ptr = ReallocMem(DATA_ARRAY, nmu1 + 1);
      WDB[DATA_ICM_PTR_MU1] = (double)ptr;
      
      for (n = 0; n < nmu1 + 1; n++)
        WDB[ptr++] = 2.0*((double)n)/((double)nmu1) - 1.0;

      WDB[DATA_ICM_NMU1] = (double)nmu1;
    }
  else
    {
      /* Check order */

      for (n = 1; n < nmu1 + 1; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Error(0, "ICM angular binning must be in ascending order");
    }

  if ((ptr = (long)RDB[DATA_ICM_PTR_MU2]) < VALID_PTR)
    {
      ptr = ReallocMem(DATA_ARRAY, nmu2 + 1);
      WDB[DATA_ICM_PTR_MU2] = (double)ptr;
      
      for (n = 0; n < nmu2 + 1; n++)
        WDB[ptr++] = 2.0*((double)n)/((double)nmu2) - 1.0;

      WDB[DATA_ICM_NMU2] = (double)nmu2;
    }
  else
    {
      /* Check order */

      for (n = 1; n < nmu2 + 1; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Error(0, "ICM angular binning must be in ascending order");
    }

  /***************************************************************************/

  /***** Allocate memory for results *****************************************/

  /* Asymmetric and symmetric bins */

  nmua = nmu1;
  nmus = nmu0*nmu2;
  
  /* Calculate number of values */
  
  sz = 2*nseg*nmua*nmus*ng0*nseg*nmua*nmus*ng0;
  sz = sz + 2*nseg*nmua*nmus*ng0*ng1;
  sz = sz + 9*nseg*nmua*nmus*ng0;
  sz = sz + 2*nseg*nmua*nmus*ng0*npmax;
  sz = sz + 2*nseg*nmua*nmus*ng0*ng1*npmax;

  /* Multiply by number of structures */
  
  loc0 = (long)RDB[DATA_PTR_ICM0];
  sz = sz*ListSize(loc0);
  
  /* Preallocate memory from buffer array */
  
  PreallocMem(sz*BUF_BLOCK_SIZE, BUF_ARRAY);
  
  /* Loop over structures and allocate memory for results */

  loc0 = (long)RDB[DATA_PTR_ICM0];
  while (loc0 > VALID_PTR)
    {
      ptr = NewStat("ICM_CURR0", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_CURR0] = (double)ptr;

      ptr = NewStat("CC1", 8, nseg, nmua, nmus, ng0, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_CC1] = (double)ptr;

      ptr = NewStat("CC2", 8, nseg, nmua, nmus, ng0, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_CC2] = (double)ptr;

      if (ng1 > 0)
        {
          ptr = NewStat("AFLX1", 5, nseg, nmua, nmus, ng0, ng1);
          WDB[loc0 + ICM_RES_AFLX1] = (double)ptr;
          
          ptr = NewStat("AFLX2", 5, nseg, nmua, nmus, ng0, ng1);
          WDB[loc0 + ICM_RES_AFLX2] = (double)ptr;
        }

      ptr = NewStat("ASRC1", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_ASRC1] = (double)ptr;

      ptr = NewStat("ASRC2", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_ASRC2] = (double)ptr;

      ptr = NewStat("AFISS1", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_AFISS1] = (double)ptr;

      ptr = NewStat("AFISS2", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_AFISS2] = (double)ptr;

      ptr = NewStat("AABS1", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_AABS1] = (double)ptr;

      ptr = NewStat("AABS2", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_AABS2] = (double)ptr;

      ptr = NewStat("LEAK1", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_LEAK1] = (double)ptr;

      ptr = NewStat("LEAK2", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_LEAK2] = (double)ptr;

      ptr = NewStat("APOW1", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_APOW1] = (double)ptr;

      ptr = NewStat("APOW2", 4, nseg, nmua, nmus, ng0);
      WDB[loc0 + ICM_RES_APOW2] = (double)ptr;

      /* Get number of pins */

      if ((m = (long)RDB[loc0 + ICM_NP]) > 0)
        {
          ptr = NewStat("PPOW1", 5, nseg, nmua, nmus, ng0, m);
          WDB[loc0 + ICM_RES_PPOW1] = (double)ptr;

          ptr = NewStat("PPOW2", 5, nseg, nmua, nmus, ng0, m);
          WDB[loc0 + ICM_RES_PPOW2] = (double)ptr;

          if (ng1 > 0)
            {
              ptr = NewStat("PFLX1", 6, nseg, nmua, nmus, ng0, ng1, m);
              WDB[loc0 + ICM_RES_PFLX1] = (double)ptr;

              ptr = NewStat("PFLX2", 6, nseg, nmua, nmus, ng0, ng1, m);
              WDB[loc0 + ICM_RES_PFLX2] = (double)ptr;
            }
        }

      /* Allocate memory for break counter */
      
      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[loc0 + ICM_BREAK_PTR_COUNT] = (double)ptr;

      /***********************************************************************/

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
