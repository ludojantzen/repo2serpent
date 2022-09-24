/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initprecdetsource.c                            */
/*                                                                           */
/* Created:       2015/11/03 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Sets up the live neutron source for dynamic simulations      */
/*              with delayed neutrons                                        */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitPrecDetSource:"

/*****************************************************************************/

void InitPrecDetSource()
{
  long loc0, loc1, tme;
  double t0;
  char tmpstr[MAX_STR];

  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  if ((long)RDB[loc0 + PRECDET_PTR_IN_FNAME] < VALID_PTR)
    return;

  /* Get pointer to time bin structure */

  tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
  CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);
      
  /* Get pointer to bins */

  tme = (long)RDB[tme + TME_PTR_BINS];
  CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    

  /* Get value of beginning of first time bin */

  t0 = RDB[tme];
  CheckValue(FUNCTION_NAME, "t0", "", t0, 0.0, INFTY);
    
#ifdef DNPRINT
  fprintf(outp, "initprecdetsrc.c-->\n");  
#endif

  /**************************************/
  /* Create initial live neutron source */
  /**************************************/

  loc1 = NewItem(DATA_PTR_SRC0, SRC_BLOCK_SIZE);

  /* Reset weight */

  WDB[loc1 + SRC_WGT] = 1.0;
              
  /* Reset boundaries */

  WDB[loc1 + SRC_XMIN] = -INFTY;
  WDB[loc1 + SRC_XMAX] =  INFTY;
  WDB[loc1 + SRC_YMIN] = -INFTY;
  WDB[loc1 + SRC_YMAX] =  INFTY;
  WDB[loc1 + SRC_ZMIN] = -INFTY;
  WDB[loc1 + SRC_ZMAX] =  INFTY;

  /* Reset point */

  WDB[loc1 + SRC_X0] = -INFTY;
  WDB[loc1 + SRC_Y0] = -INFTY;
  WDB[loc1 + SRC_Z0] = -INFTY;

  /* Set type to neutron */

  WDB[loc1 + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;

  /* Reset energy */

  WDB[loc1 + SRC_E] = -INFTY;

  /* Source name */
          
  WDB[loc1 + SRC_PTR_NAME] = (double)PutText("DynsrcLive");

  /* Print live neutron filename to tmpstr */

  sprintf(tmpstr, "%s.live", GetText(loc0 + PRECDET_PTR_IN_FNAME));

  /* Set file name */
              
  WDB[loc1 + SRC_READ_PTR_FILE] = (double)PutText(tmpstr);
              
  /* Set type: Weights will be set to zero */
  /* If we are starting our simulation from time 0.0, our live source will start from 0.0 */

  if (t0 == 0.0)
    WDB[loc1 + SRC_READ_FILE_TYPE] = (double)SRC_FILE_TYPE_S1_RENORM;
  else if (t0 > 0.0)
    WDB[loc1 + SRC_READ_FILE_TYPE] = (double)SRC_FILE_TYPE_WGT_RENORM;
  else
    Die(FUNCTION_NAME, 
        "Simulation begins at negative time or could not read initial time: %E", t0);    

  /* Set binary input on */

  WDB[loc1 + SRC_READ_BINARY] = (double)YES;

#ifdef DNPRINT
  fprintf(outp, "<-- initprecdetsrc.c\n\n");  
#endif

}
