/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : diffcoefed.c                                   */
/*                                                                           */
/* Created:       2014/01/18 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Routines for calculating new types of diffusion coefficients */
/*                                                                           */
/* Comments: - Called from several subroutines.                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DiffcoefED:"

/*****************************************************************************/

void DiffCoefED(long call, double x0, double y0, double z0, double u, double v, 
                double w, double d, double E, double wgt, long id)
{
  long ptr, in0, in1, nfg, nmg, gcu;
  long ncol, ng, n, i;
  double norm, val, x, y, z, *f;
  const double *imap;

  /* Remove this to activate */

  return;

  /* Stop tracks at outer boundary (this is new Aug 15 2014 / 2.1.22) */

  WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

  /***************************************************************************/

  /***** Part I -- Allocate memory for data **********************************/

  /* This part of the subroutine allocates memory for two structures: */
  /*                                                                  */
  /* 1) Micro-group array for storing temporary data (batch-wise      */
  /*    results) using the group structure defined in the             */
  /*    "set micro" card)                                             */
  /*                                                                  */
  /* 1) Few-group structure for storing the statistics using the      */
  /*    common few-group structure ("set nfg")                        */
  /*                                                                  */
  /* The allocation is done for every universe listed in the          */
  /* "set gcu" card.                                                  */

  if (call == 1)
    {
      /* Get pointer to microgroup energy grid */
      
      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);        
      
      /* Number of groups */
      
      nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
      
      /* Get number of groups in few-group structure */
      
      nfg = (long)RDB[DATA_ERG_FG_NG];
      
      /* Loop over gcu structures */
      
      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Allocate memory for micro-group data */
          
          ptr = AllocPrivateData(nmg, RES2_ARRAY);
          WDB[gcu + GCU_MICRO_DIFFCOEF_ED] = (double)ptr;
          
          /* Allocate memory for few-group data */
          
          ptr = NewStat("DIFFCOEF_ED", 1, nfg); 
          WDB[gcu + GCU_RES_DIFFCOEF_ED] = (double)ptr;
          
          /* Next */
          
          gcu = NextItem(gcu);
        }
    }
     
  /***************************************************************************/
  
  /***** Part II -- Score surface crossings **********************************/

  /* This part checks if the surface is crossed and stores neutron       */
  /* weight in the temporary micro-group array. The surface is currently */
  /* hard-coded (plane at x = 0), but it is possible to link any surface */
  /* from the geometry definition to this subroutine. This is done,      */
  /* for example, in superdet.c, which scores surface current for a      */
  /* super-imposed detector. */

  else if (call == 2)
    {
      /* Check that group constants are calculated */

      if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
        return;

      /* Check if active cycle */

      if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
        return;

      /* Get collision number (the following is done to get pointer to the */
      /* group constant structure. To work correctly here, the the gc      */
      /* universe should cover the entire geometry. */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);

      /* Avoid compiler warning */

      gcu = -1;

      /* Check for multiple levels */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
        {
          /* Single level, get pointer */
          
          if ((gcu = (long)TestValuePair(DATA_GCU_PTR_UNI, (double)ncol, id)) 
              < VALID_PTR)
            return;
        }
      else
        Die(FUNCTION_NAME, "This doesn't work with multiple GCU's");

      /* Start at previous position (the following takes the starting  */
      /* point of the previous neutron track, moves neutron to the end */
      /* and checks if the surface was crossed) */

      x = x0;
      y = y0;
      z = z0;

      /* Check initial position (x < 0 --> in, x >= 0 --> out) */
              
      if (x < 0.0)
        in0 = YES;
      else
        in0 = NO;

      /* Move neutron to final position (EXTRAP_L is needed in surface-   */
      /* tracking mode to make sure the neutron is moved over the surface */

      x = x + (d + EXTRAP_L)*u;
      y = y + (d + EXTRAP_L)*v;
      z = z + (d + EXTRAP_L)*w;

      /* Check final position (x < 0 --> in, x >= 0 --> out) */
              
      if (x < 0.0)
        in1 = YES;
      else
        in1 = NO;

      /* Check if surface was crossed from negative to positive side    */
      /* (to get result corresponding to the detector in the test case) */

      if ((in0 == YES) && (in1 == NO))
        {
          /* Get pointer to microgroup energy grid */
      
          ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);        
      
          /* Number of groups */
      
          nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
      
          /* Get group index (exit if energy is out of bounds) */
      
          if ((ng = GridSearch(ptr, E)) < 0)
            return;
          
          /* Convert index (to invert order) */
          
          ng = nmg - ng - 1;
          CheckValue(FUNCTION_NAME, "ng2", "", ng, 0, nmg - 1);

          /* Score surface crossing: the stored value is the second    */
          /* argument in AddPrivateRes(), currently the neutron weight */

          ptr = (long)RDB[gcu + GCU_MICRO_DIFFCOEF_ED];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng, wgt, id);
        }
    }

  /***************************************************************************/
  
  /***** Part III -- Handle the batch statistics *****************************/

  /* This part retrieves the data stored in the micro-group array during */
  /* the batch and stores it in the statistics structure. The result can */
  /* be combined here to any other data collected during the batch.      */
  /* Similar procedure for other group constants is done in microcalc.c  */

  /* NOTE: This is done in two parts, but that's just because of the way */
  /*       the data structures are handled in Serpent. */

  else if (call == 3)
    {
      /* Get pointer to micro-group structure */
  
      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Number of micro- and macro-groups (few-group structure) */

      nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
      nfg = (long)RDB[DATA_ERG_FG_NG];

      /* Index map (to map from micro- to macro-groups) */

      ptr = (long)RDB[DATA_MICRO_PTR_IDX_MAP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      imap = &RDB[ptr];

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Pointer to micro group array */

          ptr = (long)RDB[gcu + GCU_MICRO_DIFFCOEF_ED];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          f = &RES2[ptr];

          /* Loop over micro-group structure */
          
          for (n = 0; n < nmg; n++)
            {
              /* Get macro-group index */
              
              i = (long)imap[n];
              CheckValue(FUNCTION_NAME, "i", "", i, 0, nfg - 1);

              /* Add value to macro-group statistics (the value is added */
              /* to bin i in a buffer) */

              ptr = (long)RDB[gcu + GCU_RES_DIFFCOEF_ED];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);          
              AddBuf(f[n], 1.0, ptr, 0, i);
            }
          
          /* Clear micro-group array */
          
          memset(f, 0.0, nmg*sizeof(double));

          /* Next */
          
          gcu = NextItem(gcu);
        }
    }
  else if (call == 4)
    {
      /* Get normalization coefficient */

      norm = RDB[DATA_MICRO_CALC_NORM]/RDB[DATA_MICRO_CALC_BATCH_SIZE];
      CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);
      
      /* Check criticality source mode */
      
      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        norm = norm/RDB[DATA_MICRO_CALC_BATCH_SIZE];

      /* Number of macro-groups (few-group structure) */

      nfg = (long)RDB[DATA_ERG_FG_NG];

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Loop over macro-group structure and add buffered values */
          /* to final statistics */

          for (i = 0; i < nfg; i++)
            {
              /* Pointer to statistics */

              ptr = (long)RDB[gcu + GCU_RES_DIFFCOEF_ED];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);          

              /* Get value: Function BufVal() returns sum stored in buffer, */
              /* mean of stored values can be calculated using BufMean(). */

              val = BufVal(ptr, i);
              AddStat(norm*val, ptr, i);
            }

          /* Next */
          
          gcu = NextItem(gcu);
        }
    }

  /***************************************************************************/
  
  /***** Part IV -- Print results ********************************************/

  /* This part prints the results in the standard output after every */
  /* active cycle. */

  else if (call == 5)
    {
      /* Get number of groups */
      
      nfg = (long)RDB[DATA_ERG_FG_NG];

      /* Loop over universes */

      gcu = (long)RDB[DATA_PTR_GCU0];
      while (gcu > VALID_PTR)
        {
          /* Pointer to statistics */
          
          ptr = (long)RDB[gcu + GCU_RES_DIFFCOEF_ED];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);          
          
          /* Loop over groups */

          for (i = 0; i < nfg; i++)
            {
              /* Print mean value and relative statistical error */

              printf("D(%ld) = %1.5E (%1.5f)\n", i + 1, Mean(ptr, i), 
                     RelErr(ptr, i));
            }

          /* Next */
          
          gcu = NextItem(gcu);
        }
    }
  
  /***************************************************************************/
}

/*****************************************************************************/
