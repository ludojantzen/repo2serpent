/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : linkgcumaterials.c                             */
/*                                                                           */
/* Created:       2012/03/19 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Links universes involved with group constant generation to   */
/*              burnable materials when b1 spectrum correction is applied.   */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LinkGCUMaterials:"

/*****************************************************************************/

void LinkGCUMaterials(long uni, long gcu)
{
  long nst, n, reg, ptr, mat, cell, lat, pbd, loc0;

  /* Check universe pointer */

  if (uni == -1)
    {
      /* Pointer to root universe */

      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];

      /* Reset gcu pointer */

      gcu = -1;
    }

  /* Check universe pointer */

  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check and set gcu pointer */

  if ((long)RDB[uni + UNIVERSE_PTR_GCU] > VALID_PTR)
    {
      /* Check pointer */

      if (gcu > VALID_PTR)
        Error(0, "Nested universes involved with group constant generation");
      else
        gcu = (long)RDB[uni + UNIVERSE_PTR_GCU];
    }

  /* Check universe type */
  
  switch((long)RDB[uni + UNIVERSE_TYPE])
    {
    case UNIVERSE_TYPE_NEST:
      {
        /*********************************************************************/

        /***** Nest universe *************************************************/

        /* Pointer to nest */
        
        nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
        CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);
        
        /* Get pointer to regions */
        
        loc0 = (long)RDB[nst + NEST_PTR_REGIONS];
        CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
        
        /* Loop over regions */ 
        
        n = 0;
        while ((reg = ListPtr(loc0, n++)) > VALID_PTR)
          {
            /* Put region index */

            WDB[reg + NEST_REG_IDX] = (double)(n - 1);

            /* Pointer to cell */

            cell = (long)RDB[reg + NEST_REG_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Check fill and material pointers */
            
            if ((ptr = (long)RDB[reg + NEST_REG_PTR_FILL]) > VALID_PTR)
              {
                /* Filled region, call recursively */
                
                LinkGCUMaterials(ptr, gcu);
              }
            else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
              if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
                {
                  /* Check previous definition and set pointer */
          
                  if (((long)RDB[mat + MATERIAL_PTR_GCU] > VALID_PTR) &&
                      ((long)RDB[mat + MATERIAL_PTR_GCU] != gcu))
                    Error(mat, 
                          "Material used in multiple universes in B1 calculation");
                  else
                    WDB[mat + MATERIAL_PTR_GCU] = (double)gcu;
                }
          }
        
        /* Break case */
            
        break;
        
        /*********************************************************************/
      }     
    case UNIVERSE_TYPE_CELL:
      {
        /***** Cell universe *************************************************/
            
        /* Pointer to cell list */
        
        loc0 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
        CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);  
        
        /* Loop over cell list */
        
        n = 0;
        while ((cell = ListPtr(loc0, n++)) > VALID_PTR)
          {
            /* Pointer to cell */
            
            cell = (long)RDB[cell + CELL_LIST_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Check fill and material pointers */
            
            if ((ptr = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
              {
                /* Filled region, call recursively */

                LinkGCUMaterials(ptr, gcu);
              }
            else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
              if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
                {
                  /* Check previous definition and set pointer */
                  
                  if (((long)RDB[mat + MATERIAL_PTR_GCU] > VALID_PTR) &&
                      ((long)RDB[mat + MATERIAL_PTR_GCU] != gcu))
                    Error(mat, 
                          "Material used in multiple universes in B1 calculation");
                  else
                    WDB[mat + MATERIAL_PTR_GCU] = (double)gcu;
                }
          }
        
        /* Break case */
        
        break;
        
        /*********************************************************************/
      }      
    case UNIVERSE_TYPE_LATTICE:
      {
        /***** Lattice universe **********************************************/
        
        /* Pointer to lattice */
        
        lat = (long)RDB[uni + UNIVERSE_PTR_LAT];
        CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);
        
        /* Check type */
        
        if ((long)RDB[lat + LAT_TYPE] == LAT_TYPE_CLU)
          {
            /***** Circular array ********************************************/
            
            /* Get pointer to rings */
            
            loc0 = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
            
            /* Loop over rings */ 
            
            n = 0;
            while ((reg = ListPtr(loc0, n++)) > VALID_PTR)
              {
                /* Pointer to items */ 
                
                ptr = (long)RDB[reg + RING_PTR_FILL];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                
                /* Loop over items */
                
                while ((long)RDB[ptr] > VALID_PTR)
                  {
                    /* Call recursively */

                    LinkGCUMaterials((long)RDB[ptr], gcu);
                    
                    /* Next */
                    
                    ptr++;
                  }
              }
            
            /*****************************************************************/
          }
        else
          {
            /***** Simple types **********************************************/
            
            /* Pointer to items */ 
            
            ptr = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

            /* Loop over items and call recursively (tuolla voi olla */
            /* se NULLPTR välissä) */
            
            for (n = 0; n < (long)RDB[lat + LAT_NTOT]; n++)
              if ((uni = (long)RDB[ptr + n]) > VALID_PTR)
                LinkGCUMaterials(uni, gcu);
                        
            /*****************************************************************/
          }

        /* Break case */
        
        break;
        
        /*********************************************************************/
      }
    case UNIVERSE_TYPE_PBED:
      {
        /***** Explicit stochastic geometry **********************************/
        
        /* Pointer to geometry */
        
        pbd = (long)RDB[uni + UNIVERSE_PTR_PBED];
        CheckPointer(FUNCTION_NAME, "(pbd)", DATA_ARRAY, pbd);        
        
        /* Pointer to background universe */

        loc0 = (long)RDB[pbd + PBED_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);        

        /* Call recursively */
        
        LinkGCUMaterials(loc0, gcu);

        /* Loop over pebble types */
      
        loc0 = (long)RDB[pbd + PBED_PTR_PEBBLE_TYPES];
        while (loc0 > VALID_PTR)
          {
            /* Pointer to universe */

            ptr = (long)RDB[loc0 + PEBTYPE_PTR_UNIV];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);        

            /* Counter */

            if ((n = (long)RDB[loc0 + PEBTYPE_COUNT]) < 1)
              Die(FUNCTION_NAME, "Zero pebble counter");

            /* Call recursively */
                
            LinkGCUMaterials(loc0, gcu);

            /* Next type */
            
            loc0 = NextItem(loc0);
          }

        /* Break case */
        
        break;
        
        /*********************************************************************/
      }
    default:
      {
        /* Invalid type */
        
        Die(FUNCTION_NAME, "Invalid universe type");
      }
    }
}

/*****************************************************************************/
