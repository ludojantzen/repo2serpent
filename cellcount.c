/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cellcount.c                                    */
/*                                                                           */
/* Created:       2011/07/03 (JLe)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Counts the number of cells                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CellCount:"

/*****************************************************************************/

void CellCount(long uni, long lvl, long recu, long add)
{
  long nst, ptr, reg, n, cell, lat, loc0, loc1, pbd, umsh, stl, ifc, prnt;
  long tetlist, ntet, tet, i;

  /* Print */

  if (recu == 0)
    fprintf(outp, "Counting cells...\n");

  /* Check if first */

  if (lvl < VALID_PTR)
    {
      /* Reset cell counters */

      cell = (long)RDB[DATA_PTR_C0];
      while (cell > VALID_PTR)
        {
          /* Reset counter */

          WDB[cell + CELL_VOL_COUNT] = 0.0;

          /* Next */

          cell = NextItem(cell);
        }

      /* Reset nest counters */

      nst = (long)RDB[DATA_PTR_NST0];
      while (nst > VALID_PTR)
        {
          /* Reset counter */

          WDB[nst + NEST_COUNT] = 0.0;

          /* Next */

          nst = NextItem(nst);
        }

      /* Get pointer to first level and universe */

      lvl = (long)RDB[DATA_PTR_LVL0];
      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
    }
  else
    {
      /* Update level pointer */

      lvl = NextItem(lvl);
    }

  /* Check level and universe pointers */

  CheckPointer(FUNCTION_NAME, "(lvl)", DATA_ARRAY, lvl);
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check infinite loop */

  if (recu++ > 1000)
    Die(FUNCTION_NAME, "Infinite geometry loop involving universe %s",
        GetText(uni + UNIVERSE_PTR_NAME));

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
            /* Get pointer to cell */

            cell = (long)RDB[reg + NEST_REG_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Update counter */

            WDB[cell + CELL_VOL_COUNT] = RDB[cell + CELL_VOL_COUNT]
              + (double)add;

            /* Check fill pointer */

            if ((ptr = (long)RDB[reg + NEST_REG_PTR_FILL]) > VALID_PTR)
              {
                /* Filled region, call recursively */

                CellCount(ptr, lvl, recu, add);
              }
          }

        /* Add nest counter */

        WDB[nst + NEST_COUNT] = RDB[nst + NEST_COUNT] + 1.0;

        /* Break case */

        break;

        /*********************************************************************/
      }
    case UNIVERSE_TYPE_CELL:
      {
        /*********************************************************************/

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

            /* Update counter */

            WDB[cell + CELL_VOL_COUNT] = RDB[cell + CELL_VOL_COUNT] + 1.0;

            /* Check fill pointer */

            if ((ptr = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
              {
                /* Filled region, call recursively */

                CellCount(ptr, lvl, recu, add);
              }
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

                    CellCount((long)RDB[ptr], lvl, recu, add);

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
                CellCount(uni, lvl, recu, add);

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

        loc0 = (long)RDB[pbd + PBED_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

        /* Call recursively */

        CellCount(loc0, lvl, recu, add);

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

            CellCount(ptr, lvl, recu, add*n);

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

        /* Pointer to background universe */

        loc0 = (long)RDB[umsh + UMSH_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

        /* Call recursively */

        CellCount(loc0, lvl, recu, add);

        /* Pointer to interface structure */

        ifc = (long)RDB[umsh + UMSH_PTR_IFC];
        CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

        tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];

        ntet = (long)RDB[ifc + IFC_NC];

        for (i = 0; i < ntet; i++)
          {
            /* Get pointer to tet */

            tet = (long)RDB[tetlist + i];

            /* Get pointer to parent */

            prnt = (long)RDB[tet + TET_PTR_PARENT];

            /* Pointer to geometry cell */

            cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Update counter */

            WDB[cell + CELL_VOL_COUNT] = RDB[cell + CELL_VOL_COUNT] + 1.0;

            /* Check fill pointer */

            if ((long)RDB[cell + CELL_PTR_FILL] > VALID_PTR)
              Die(FUNCTION_NAME, "Fill pointer not null");
          }

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

        /* Pointer to background universe */

        loc0 = (long)RDB[stl + STL_PTR_BG_UNIV];
        CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

        /* Call recursively */

        CellCount(loc0, lvl, recu, add);

        /* Pointer to solids */

        loc1 = (long)RDB[stl + STL_PTR_SOLIDS];
        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

        /* Loop over solids */

        while (loc1 > VALID_PTR)
          {
            /* Pointer to geometry cell */

            cell = (long)RDB[loc1 + STL_SOLID_PTR_CELL];
            CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

            /* Update counter */

            WDB[cell + CELL_VOL_COUNT] = RDB[cell + CELL_VOL_COUNT] + 1.0;

            /* Check fill pointer */

            if ((ptr = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
              {
                /* Filled region, call recursively */

                CellCount(ptr, lvl, recu, add);
              }

            /* Next solid */

            loc1 = NextItem(loc1);
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

  /* Print */

  if (recu == 1)
    fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
