/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : zonecount.c                                    */
/*                                                                           */
/* Created:       2012/05/10 (JLe)                                           */
/* Last modified: 2018/10/16 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Counts the number of zones at geometry levels                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ZoneCount:"

/*****************************************************************************/

void ZoneCount(long uni, long lvl, long recu)
{
  long nst, ptr, reg, n, cell, lat, loc0, pbd, umsh, stl, max, lst, i, prod;

  /* Update level pointer or get pointer to first level and universe */

  if (lvl < VALID_PTR)
    {
      lvl = (long)RDB[DATA_PTR_LVL0];
      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
    }
  else
    lvl = NextItem(lvl);

  /* Print */

  if (recu == 0)
    fprintf(outp, "Counting geometry zones...\n\n");

  /* Check universe pointer */

  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check level (30.7.2013 / 2.1.15) */

  if (lvl < VALID_PTR)
    {
      /* Check recursion and geometry levels */

      if (recu == (long)RDB[DATA_GEOM_LEVELS])
        Error(0, "Major geometry error in universe structure");
      else
        Die(FUNCTION_NAME, "WTF?");
    }
  else if ((long)RDB[lvl + LVL_NUMBER] != (long)RDB[uni + UNIVERSE_LEVEL])
    {
      /* Print warning (first occasion only) */

      if (RDB[uni + UNIVERSE_WARN_MULTI_LVL] == 0)
        Note(0, "Universe %s is used at multiple levels (%ld and %ld)",
             GetText(uni + UNIVERSE_PTR_NAME), (long)RDB[lvl + LVL_NUMBER],
             (long)RDB[uni + UNIVERSE_LEVEL]);

      /* Add counter */

      WDB[uni + UNIVERSE_WARN_MULTI_LVL] =
        RDB[uni + UNIVERSE_WARN_MULTI_LVL] + 1.0;
    }

  /* Check infinite loop */

  if (recu++ > 1000)
    Die(FUNCTION_NAME, "Infinite geometry loop involving universe %s",
        GetText(uni + UNIVERSE_PTR_NAME));

  /* Reset maximum number of zones */

  max = 0;

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

        reg = (long)RDB[nst + NEST_PTR_REGIONS];
        CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

        /* Reset region index */

        i = 0;

        /* Loop over regions */

        while (reg > VALID_PTR)
          {
            /* Put region index */

            WDB[reg + NEST_REG_IDX] = (double)(i++);

            /* Update counter */

            max++;

            /* Check fill pointer and call recursively */

            if ((ptr = (long)RDB[reg + NEST_REG_PTR_FILL]) > VALID_PTR)
              ZoneCount(ptr, lvl, recu);

            /* Next region */

            reg = NextItem(reg);
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_CELL:
      {
        /*********************************************************************/

        /***** Cell universe *************************************************/

        /* Pointer to cell list */

        lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
        CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

        /* Reset region index */

        i = 0;

        /* Loop over cell list */

        while (lst > VALID_PTR)
          {
            /* Put region index */

            WDB[lst + CELL_LIST_REG_IDX] = (double)(i++);

            /* Pointer to cell */

            cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Update counter */

            max++;

            /* Check fill pointer and call recursively */

            if ((ptr = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
              ZoneCount(ptr, lvl, recu);

            /* Next */

            lst = NextItem(lst);
          }

        /* Pointer to source cell list */

        lst = (long)RDB[uni + UNIVERSE_PTR_SRC_CELL_LIST];
        CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

        /* Reset region index */

        i = 0;

        /* Loop over cell list */

        while (lst > VALID_PTR)
          {
            /* Put region index */

            WDB[lst + CELL_LIST_REG_IDX] = (double)(i++);

            /* Next */

            lst = NextItem(lst);
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_LATTICE:
      {
        /*********************************************************************/

        /***** Lattice universe **********************************************/

        /* Pointer to lattice */

        lat = (long)RDB[uni + UNIVERSE_PTR_LAT];
        CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

        /* Get maximum */

        max = (long)RDB[lat + LAT_NTOT];
        CheckValue(FUNCTION_NAME, "max", "", max, 1, 1000000000);

        /* Check type */

        if ((long)RDB[lat + LAT_TYPE] == LAT_TYPE_CLU)
          {
            /*****************************************************************/

            /***** Circular array ********************************************/

            /* Get pointer to rings */

            reg  = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

            /* Loop over rings */

            while (reg > VALID_PTR)
              {
                /* Pointer to items */

                ptr = (long)RDB[reg + RING_PTR_FILL];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                /* Loop over items */

                while ((long)RDB[ptr] > VALID_PTR)
                  {
                    /* Call recursively */

                    ZoneCount((long)RDB[ptr], lvl, recu);

                    /* Next */

                    ptr++;
                  }

                /* Next ring */

                reg = NextItem(reg);
              }

            /*****************************************************************/
          }
        else
          {
            /*****************************************************************/

            /***** Simple types **********************************************/

            /* Pointer to items */

            ptr = (long)RDB[lat + LAT_PTR_FILL];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

            /* Loop over items and call recursively (tuolla voi olla */
            /* se NULLPTR välissä) */

            for (n = 0; n < (long)RDB[lat + LAT_NTOT]; n++)
              if ((uni = (long)RDB[ptr + n]) > VALID_PTR)
                ZoneCount(uni, lvl, recu);

            /*****************************************************************/
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_PBED:
      {
        /*********************************************************************/

        /***** Explicit stochastic geometry **********************************/

        /* Pointer to geometry */

        pbd = (long)RDB[uni + UNIVERSE_PTR_PBED];
        CheckPointer(FUNCTION_NAME, "(pbd)", DATA_ARRAY, pbd);

        /* Pointer to background universe */

        ptr = (long)RDB[pbd + PBED_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get count */

        max = (long)RDB[pbd + PBED_N_PEBBLES] + 1;

        /* Call recursively */

        ZoneCount(ptr, lvl, recu);

        /* Loop over pebble types */

        loc0 = (long)RDB[pbd + PBED_PTR_PEBBLE_TYPES];
        while (loc0 > VALID_PTR)
          {
            /* Pointer to universe */

            ptr = (long)RDB[loc0 + PEBTYPE_PTR_UNIV];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

            /* Counter */

            if ((long)RDB[loc0 + PEBTYPE_COUNT] < 1)
              Die(FUNCTION_NAME, "Zero pebble counter");

            /* Call recursively */

            ZoneCount(ptr, lvl, recu);

            /* Next type */

            loc0 = NextItem(loc0);
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_UMSH:
      {
        /*********************************************************************/

        /***** Unstructured mesh based geometry ******************************/

        /* Pointer to geometry */

        umsh = (long)RDB[uni + UNIVERSE_PTR_UMSH];
        CheckPointer(FUNCTION_NAME, "(umsh)", DATA_ARRAY, umsh);

        /* Pointer to cells */

        ptr = (long)RDB[umsh + UMSH_PTR_IFC];
        CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

        max = (long)RDB[ptr + IFC_NC];

        /* Pointer to background universe */

        ptr = (long)RDB[umsh + UMSH_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(ptr3)", DATA_ARRAY, ptr);

        /* Call recursively */

        ZoneCount(ptr, lvl, recu);

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_STL:
      {
        /*********************************************************************/

        /***** STL based geometry ********************************************/

        /* Pointer to geometry */

        stl = (long)RDB[uni + UNIVERSE_PTR_STL];
        CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);

        /* Pointer to solids */

        ptr = (long)RDB[stl + STL_PTR_SOLIDS];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get count */

        max = ListSize(ptr) + 1;

        /* Pointer to background universe */

        ptr = (long)RDB[stl + STL_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Call recursively */

        ZoneCount(ptr, lvl, recu);

        /* Reset region index (index 0 is for background) */

        i = 1;

        /* Loop over solids */

        loc0 = (long)RDB[stl + STL_PTR_SOLIDS];
        while (loc0 > VALID_PTR)
          {
            /* Add to count */

            WDB[loc0 + STL_SOLID_REG_IDX] = (double)(i++);

            /* NOTE: Noita cellejä ei oo vielä prosessoitu tässä */
            /* vaiheessa */

            /* Pointer to  cell */
            /*
            cell = (long)RDB[loc0 + STL_SOLID_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
            */
            /* Check fill pointer and call recursively */
            /*
            if ((ptr = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
              ZoneCount(ptr, lvl, recu);
            */
            /* Next solid */

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

  /* Compare number of regions to maximum */

  if (max > (long)RDB[lvl + LVL_MAX_REGIONS])
    WDB[lvl + LVL_MAX_REGIONS] = (double)max;

  /* Reset product */

  prod = 1;

  /* Calculate cumulative */

  if ((long)RDB[lvl + LVL_NUMBER] == 0)
    {
      /* Loop over levels */

      while (lvl > VALID_PTR)
        {
          /* Put multiplier */

          WDB[lvl + LVL_ZONE_IDX_MULT] = (double)prod;

          /* Check numerical limit */

          if (RDB[lvl + LVL_MAX_REGIONS]*RDB[lvl + LVL_ZONE_IDX_MULT] > 9E+18)
            Error(0, "Indexing overflow, reduce number of geometry levels");

          /* Update multiplier */

          prod = prod*((long)WDB[lvl + LVL_MAX_REGIONS]);

          /* Put maximum cumulative */

          WDB[lvl + LVL_CUM_MAX_REGIONS] = (double)prod;

          /* Next level */

          lvl = NextItem(lvl);
        }

      /* Print */

      if ((long)RDB[DATA_GEOM_LEVELS] > 1)
        {
          fprintf(outp,
                  "The geometry consists of %ld levels:\n\n",
                  (long)RDB[DATA_GEOM_LEVELS]);

          /* Loop over levels and print */

          lvl = (long)RDB[DATA_PTR_LVL0];
          while (lvl > VALID_PTR)
            {
              /* Print */

              fprintf(outp, "Level %ld size: max %ld zones\n",
                      (long)RDB[lvl + LVL_NUMBER],
                      (long)WDB[lvl + LVL_MAX_REGIONS]);

              /* Next level */

              lvl = NextItem(lvl);
            }

          fprintf(outp, "\n");
        }

      /* Allocate memory for regional collision index */

      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[DATA_PTR_ZONE_IDX] = (double)ptr;
    }
}

/*****************************************************************************/
