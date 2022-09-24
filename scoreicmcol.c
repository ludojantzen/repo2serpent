/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreicmcol.c                                  */
/*                                                                           */
/* Created:       2013/10/01 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Scores reaction rates needed for ICM reconstruction factors  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreICMCol:"

/*****************************************************************************/

void ScoreICMCol(double flx, double capt, double fiss, double nsf, 
                 double fissE, double sprod, long mat, long part, double x, 
                 double y, double z, double E, double wgt, double g, long id)
{
  long ptr, lat, icm, ng1, ntot, ncol, n0, ma0, ms0, ng0, m;
  double wgt0;

  /* Check mode */

  if ((long)RDB[DATA_ICM_CALC] == NO)
    return;
  
  /* Check pointer */
      
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);
        
  /* Get data */
  
  n0 = (long)RDB[part + PARTICLE_ICM_IDX];
  ma0 = (long)RDB[part + PARTICLE_ICM_MUA];
  ms0 = (long)RDB[part + PARTICLE_ICM_MUS];
  ng0 = (long)RDB[part + PARTICLE_ICM_G];
  wgt0 = RDB[part + PARTICLE_ICM_WGT];
  
  /* Check */
  
  if ((n0 < 0) || (ma0 < 0) || (ms0 < 0) || (ng0 < 0))
    return;

  /* Get pointer */

  if ((icm = (long)RDB[part + PARTICLE_ICM_PTR_ICM]) < VALID_PTR)
    return;
    
  /* Get few-group index */

  ng1 = -1;
  
  if ((ptr = (long)RDB[DATA_ICM_PTR_ENE1]) > VALID_PTR)
    if ((ng1 = GridSearch(ptr, E)) > -1)
      {
        /* Total number of groups */
        
        ntot = (long)RDB[DATA_ICM_NG1];
        
        /* Convert index */
        
        ng1 = ntot - ng1 - 1;
        CheckValue(FUNCTION_NAME, "ng1", "", ng1, 0, ntot - 1);
      }

  /***************************************************************************/

  /***** Assembly-wise data **************************************************/

  /* Score */
  
  ptr = (long)RDB[icm + ICM_RES_AFISS1];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(fiss, wgt, ptr, id, -1, n0, ma0, ms0, ng0);

  ptr = (long)RDB[icm + ICM_RES_AFISS2];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(wgt0*fiss, wgt, ptr, id, -1, n0, ma0, ms0, ng0);

  ptr = (long)RDB[icm + ICM_RES_AABS1];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  /*
  AddBuf(fiss + capt - sprod, wgt, ptr, id, -1, n0, ma0, ms0, ng0);
  */
  AddBuf(fiss + capt, wgt, ptr, id, -1, n0, ma0, ms0, ng0);

  ptr = (long)RDB[icm + ICM_RES_AABS2];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  /*
  AddBuf(wgt0*(fiss + capt - sprod), wgt, ptr, id, -1, n0, ma0, ms0, ng0);
  */
  AddBuf(wgt0*(fiss + capt), wgt, ptr, id, -1, n0, ma0, ms0, ng0);

  ptr = (long)RDB[icm + ICM_RES_ASRC1];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(nsf, wgt0, ptr, id, -1, n0, ma0, ms0, ng0);
  
  ptr = (long)RDB[icm + ICM_RES_ASRC2];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(wgt*nsf, wgt0, ptr, id, -1, n0, ma0, ms0, ng0);
      
  ptr = (long)RDB[icm + ICM_RES_APOW1];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(fissE, wgt, ptr, id, -1, n0, ma0, ms0, ng0);

  ptr = (long)RDB[icm + ICM_RES_APOW2];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(wgt0*fissE, wgt, ptr, id, -1, n0, ma0, ms0, ng0);
  
  /* Check energy group */

  if (ng1 > -1)
    {
      /* Score flux */
      
      ptr = (long)RDB[icm + ICM_RES_AFLX1];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(flx, wgt, ptr, id, -1, n0, ma0, ms0, ng0, ng1);

      ptr = (long)RDB[icm + ICM_RES_AFLX2];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(wgt0*flx, wgt, ptr, id, -1, n0, ma0, ms0, ng0, ng1);
    }
  
  /***************************************************************************/

  /***** Pin-wise data *******************************************************/

  /* Check lattice pointer */
  
  if ((lat = (long)RDB[icm + ICM_PTR_LAT]) > VALID_PTR)
    /*
    if ((m = (long)TestValuePair(lat + LAT_COL_COUNT, (double)ncol, id)) > -1)
    */
    {
      /* Get relative coordinates */

      x = x - RDB[icm + ICM_X0];
      y = y - RDB[icm + ICM_Y0];
      z = z - RDB[icm + ICM_Z0];

      /* Get region */
      
      if ((m = FindLatticeRegion(lat, -1, &x, &y, &z, NULL, id)) > -1)
        {
          /* Check maximum */

          if (m > (long)RDB[icm + ICM_NP] - 1)
            Die(FUNCTION_NAME, "Error in index");

          /* Score power */
      
          ptr = (long)RDB[icm + ICM_RES_PPOW1];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(fissE, wgt, ptr, id, -1, n0, ma0, ms0, ng0, m);
          
          ptr = (long)RDB[icm + ICM_RES_PPOW2];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(wgt0*fissE, wgt, ptr, id, -1, n0, ma0, ms0, ng0, m);

          /* Score flux */
          
          if (ng1 > -1)
            {        
              ptr = (long)RDB[icm + ICM_RES_PFLX1];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(flx, wgt, ptr, id, -1, n0, ma0, ms0, ng0, ng1, m);

              ptr = (long)RDB[icm + ICM_RES_PFLX2];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(wgt0*flx, wgt, ptr, id, -1, n0, ma0, ms0, ng0, ng1, m);
            }
        }
    }

  /***************************************************************************/
}

/*****************************************************************************/
