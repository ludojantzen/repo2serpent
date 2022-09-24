/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printtmsdiagnostics.c                          */
/*                                                                           */
/* Created:       2016/04/12 (Tvi)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Checks if the proportion of erroneous TMS samples is too     */
/*              high and advices the user to decrease qparams                */
/*                                                                           */
/* Comments:  - Serpent 1:ssä kutsutaan lopussa rutiinia joka printtaa       */
/*              kootusti kaikki variotukset. Jotain vastaavaa voisi harkita  */
/*              tähänkin versioon.                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h" 
#include "locations.h"

#define FUNCTION_NAME "PrintTMSDiagnostics:"

/*****************************************************************************/

void PrintTMSDiagnostics()
{
  double mean, std, div;
  long ptr;

  /* TMS */

  if ((long)RDB[DATA_TMS_MODE] != TMS_MODE_NONE)
    { 
      ptr = (long)RDB[RES_TMS_FAIL_STAT];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
      
      mean = Mean(ptr, 1);
      std = RelErr(ptr, 1);
      
      if ((std < 0.2) && (mean > 1E-6))
        {
          fprintf(outp, "WARNING:\n");
          fprintf(outp, "The TMS majorant cross section was exceeded frequently (proportion %.2E > 1E-6)\n", mean);
          fprintf(outp, "This may have some effect on the results.\n");
          fprintf(outp, "Please re-run the calculation with \"set qparam_tms 1E-6\".\n\n");
        }
    }

  /* DBRC */
    
  if ((long)RDB[DATA_USE_DBRC] == YES)
    {
      /* DBRC majorant exceed fraction */
      
      ptr = (long)RDB[DATA_PTR_DBRC_COUNT];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      div = SumPrivateData(ptr);
      
      ptr = (long)RDB[DATA_PTR_DBRC_EXCEED_COUNT];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);      

      if (div > 0.0)
        mean = SumPrivateData(ptr)/div; 
      else
        mean = 0.0;
          
      if(mean > 5E-6)
        {
          fprintf(outp, "WARNING:\n");
          fprintf(outp, "The DBRC majorant cross section was exceeded frequently (proportion %.2E > 5E-6)\n", mean);
          fprintf(outp, "This may have some effect on the results.\n");
          fprintf(outp, "Please re-run the calculation with \"set qparam_dbrc 1E-6\".\n\n");      
        }
    }
}

/*****************************************************************************/
