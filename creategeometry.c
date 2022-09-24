/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : creategeometry.c                               */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2018/09/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Creates physical and super-imposed universes               */
/*                                                                           */
/* Comments: - Toi "level" on huono termi, kun oikeammin kyse on             */
/*             väliaikaisen datan varastoinnista tracking-rutiinin aikana    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CreateGeometry:"

/*****************************************************************************/

void CreateGeometry()
{
  long n, src, det, loc0, ptr, uni, cell, stl, icm, lat, lvl, max, gcu, adf;
  long alb;
  char name[MAX_STR];

  fprintf(outp, "Creating geometry...\n");

  /***************************************************************************/

  /***** Physical universes **************************************************/

  /* Create physical universes starting from root */

  sprintf(name, "%s", GetText(DATA_PTR_ROOT_UNIVERSE));
  ptr = CreateUniverse(0, name, 0);

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Put pointer */

  WDB[DATA_PTR_ROOT_UNIVERSE] = (double)ptr;

  /* Create levels */

  for (n = 0; n < (long)RDB[DATA_GEOM_LEVELS]; n++)
    {
      /* Create data structure */

      lvl = NewItem(DATA_PTR_LVL0, LVL_BLOCK_SIZE);

      /* Allocate memory for private data */

      ptr = AllocPrivateData(LVL_PRIV_BLOCK_SIZE, PRIVA_ARRAY);
      WDB[lvl + LVL_PTR_PRIVATE_DATA] = (double)ptr;

      /* Put index */

      WDB[lvl + LVL_NUMBER] = (double)n;
    }

  /* Remember maximum level */

  max = (long)RDB[DATA_GEOM_LEVELS];

  /* Allocate memory for collision universe */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_COLLISION_UNI] = (double)ptr;

  /***************************************************************************/

  /***** Super-imposed universes for sources *********************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Check universe pointer */

      if ((long)RDB[src + SRC_PTR_UNIV] > VALID_PTR)
        {
          /* Loop over existing */

          uni = (long)RDB[DATA_PTR_U0];
          while (uni > VALID_PTR)
            {
              /* Compare names */

              if (!strcmp(GetText(src + SRC_PTR_UNIV),
                          GetText(uni + UNIVERSE_PTR_NAME)))
                break;

                /* Next universe */

              uni = NextItem(uni);
            }

          /* Check if found */

          if (uni < VALID_PTR)
            {
              /* Create new universe */

              sprintf(name, "%s", GetText(src + SRC_PTR_UNIV));
              uni = CreateUniverse(src, name, 0);
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Check type */

              if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_CELL)
                Error(src, "Super-imposed universe must consist of cells");

              /* Put type */

              if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
                Error(src, "Sub-levels not allowed with universe %s",
                      name);
              else
                WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
            }

          /* Put pointer */

          WDB[src + SRC_PTR_UNIV] = (double)uni;
        }

      /* Check cell pointer */

      if ((long)RDB[src + SRC_PTR_CELL] > VALID_PTR)
        {
          /* Loop over cells to find match */

          cell = (long)RDB[DATA_PTR_C0];
          while (cell > VALID_PTR)
            {
              /* Compare names */

              if (!strcmp(GetText(src + SRC_PTR_CELL),
                          GetText(cell + CELL_PTR_NAME)))
                {
                  /* Set pointer (ei toimi STL-geometrioiden kanssa, */
                  /* asetetaan processources.c:ssä). */

                  /*
                  WDB[src + SRC_PTR_CELL] = (double)cell;
                  */
                  /* Check used-flag */

                  if (!((long)RDB[cell + CELL_OPTIONS] & OPT_USED))
                    {
                      /* Create universe */

                      sprintf(name, "%s", GetText(cell + CELL_PTR_UNI));
                      uni = CreateUniverse(src, name, 0);
                      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

                      /* Put pointer */

                      if (RDB[src + SRC_PTR_UNIV] < VALID_PTR)
                        WDB[src + SRC_PTR_UNIV] = (double)uni;

                      /* Put universe type */

                      if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
                        Error(src, "Sub-levels not allowed with universe %s",
                              name);
                      else
                        WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
                    }

                  break;
                }

              /* Next cell */

              cell = NextItem(cell);
            }
        }

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Super-imposed universes for group constant generation ***************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Loop over existing */

      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Compare names */

          if (!strcmp(GetText(gcu + GCU_PTR_UNIV),
                      GetText(uni + UNIVERSE_PTR_NAME)))
            break;

          /* Next universe */

          uni = NextItem(uni);
        }

      /* Check if found */

      if (uni < VALID_PTR)
        {
          /* Create new universe */

          sprintf(name, "%s", GetText(gcu + GCU_PTR_UNIV));
          uni = CreateUniverse(0, name, 0);
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Check type */

          if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_CELL)
            Error(src, "Super-imposed universe must consist of cells");

          /* Put type */

          if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
            Error(gcu, "Sub-levels not allowed with universe %s", name);
          else
            WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
        }

      /* Next */

      gcu = NextItem(gcu);
    }

  /* Loop over ADFs (GCUs not created yet) */

  adf = (long)RDB[DATA_PTR_ADF0];
  while (adf > VALID_PTR)
    {
      /* Loop over existing */

      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Compare names */

          if (!strcmp(GetText(adf + ADF_PTR_GCU),
                      GetText(uni + UNIVERSE_PTR_NAME)))
            break;

          /* Next universe */

          uni = NextItem(uni);
        }

      /* Check if found */

      if (uni < VALID_PTR)
        {
          /* Create new universe */

          sprintf(name, "%s", GetText(adf + ADF_PTR_GCU));
          uni = CreateUniverse(0, name, 0);
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Check type */

          if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_CELL)
            Error(src, "Super-imposed universe must consist of cells");

          /* Put type */

          if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
            Error(adf, "Sub-levels not allowed with universe %s", name);
          else
            WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
        }

      /* Next */

      adf = NextItem(adf);
    }

  /* Loop over ALBs (GCUs not created yet) */

  alb = (long)RDB[DATA_PTR_ALB0];
  while (alb > VALID_PTR)
    {
      /* Loop over existing */

      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Compare names */

          if (!strcmp(GetText(alb + ALB_PTR_GCU),
                      GetText(uni + UNIVERSE_PTR_NAME)))
            break;

          /* Next universe */

          uni = NextItem(uni);
        }

      /* Check if found */

      if (uni < VALID_PTR)
        {
          /* Create new universe */

          sprintf(name, "%s", GetText(alb + ALB_PTR_GCU));
          uni = CreateUniverse(0, name, 0);
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Check type */

          if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_CELL)
            Error(src, "Super-imposed universe must consist of cells");

          /* Put type */

          if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
            Error(alb, "Sub-levels not allowed with universe %s", name);
          else
            WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
        }

      /* Next */

      alb = NextItem(alb);
    }

  /***************************************************************************/

  /***** Super-imposed universes for detectors *******************************/

  /* NOTE: Pointers are set in ProcessDetectors() */

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Loop over universe bins */

      loc0 = (long)RDB[det + DET_PTR_UBINS];
      while (loc0 > VALID_PTR)
        {
          /* Loop over existing */

          uni = (long)RDB[DATA_PTR_U0];
          while (uni > VALID_PTR)
            {
              /* Compare names */

              if (!strcmp(GetText(loc0 + DET_UBIN_PTR_UNI),
                          GetText(uni + UNIVERSE_PTR_NAME)))
                break;

              /* Next universe */

              uni = NextItem(uni);
            }

          /* Check if found */

          if (uni < VALID_PTR)
            {
              /* Create new universe */

              sprintf(name, "%s", GetText(loc0 + DET_UBIN_PTR_UNI));
              uni = CreateUniverse(det, name, 0);
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Check type */

              if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_CELL)
                Error(det, "Super-imposed universe must consist of cells");

              /* Put type */

              if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
                Error(det, "Sub-levels not allowed with universe %s",
                      name);
              else
                WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
            }

          /* Next universe bin */

          loc0 = NextItem(loc0);
        }

      /* Loop over cell bins */

      loc0 = (long)RDB[det + DET_PTR_CBINS];
      while (loc0 > VALID_PTR)
        {
          /* Check for UMSH type */

          if ((long)RDB[loc0 + DET_CBIN_UMSH_PTR_UMSH] > VALID_PTR)
            {
              /* Next bin */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Loop over cells to find match */

          cell = (long)RDB[DATA_PTR_C0];
          while (cell > VALID_PTR)
            {
              /* Compare names */

              if (!strcmp(GetText(loc0 + DET_CBIN_PTR_CELL),
                          GetText(cell + CELL_PTR_NAME)))
                break;

              /* Next cell */

              cell = NextItem(cell);
            }

          /* Check pointer */

          if (cell < VALID_PTR)
            {
              /* Cells in STL solids are not yet created, but they should */
              /* not be super-imposed either */

              stl = (long)RDB[DATA_PTR_STL0];
              while (stl > VALID_PTR)
                {
                  /* Loop over bodies */

                  ptr = (long)RDB[stl + STL_PTR_BODIES];
                  while(ptr > VALID_PTR)
                    {
                      /* Compare */

                      if (CompareStr(ptr + STL_BODY_PTR_CNAME,
                                     loc0 + DET_CBIN_PTR_CELL))
                        break;

                      /* Next */

                      ptr = NextItem(ptr);
                    }

                  /* Check pointer */

                  if (ptr > VALID_PTR)
                    break;

                  /* Next */

                  stl = NextItem(stl);
                }

              /* Check pointer */

              if (stl > VALID_PTR)
                {
                  /* Next cell bin */

                  loc0 = NextItem(loc0);

                  /* Cycle loop */

                  continue;
                }
              else
                Error(det, "Cell %s in detector %s not defined",
                      GetText(loc0 + DET_CBIN_PTR_CELL),
                      GetText(det + DET_PTR_NAME));
            }

          /* Check used-flag */

          if (!((long)RDB[cell + CELL_OPTIONS] & OPT_USED))
            {
              /* Create universe */

              /* TKa -> */

              sprintf(name, "%s", GetText(cell + CELL_PTR_UNI));
              uni = CreateUniverse(det, name, 0);
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Put universe type */

              if ((long)RDB[uni + UNIVERSE_OPTIONS] & OPT_UNIV_SUB)
                Error(det, "Sub-levels not allowed with universe %s",
                      name);
              else
                WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_SUPER;
            }

          /* Next cell bin */

          loc0 = NextItem(loc0);
        }

      /* Next detector */

      det = NextItem(det);
    }

  /***************************************************************************/

  /***** Super-imposed lattices for ICM **************************************/

  /* Loop over structures */

  icm = (long)RDB[DATA_PTR_ICM0];
  while (icm > VALID_PTR)
    {
      /* Check pointer */

      if ((long)RDB[icm + ICM_PTR_LAT] < VALID_PTR)
        {
          /* Pointer to next */

          icm = NextItem(icm);

          /* Cycle loop */

          continue;
        }

      /* Find lattice */

      lat = (long)RDB[DATA_PTR_L0];
      if ((lat = SeekListStr(lat, LAT_PTR_NAME,
                             GetText(icm + ICM_PTR_LAT))) < VALID_PTR)
        Error(icm, "Lattice %s is not defined",
              GetText(icm + ICM_PTR_LAT));

      /* Set used-flag */

      SetOption(lat + LAT_OPTIONS, OPT_USED);

      /* Next */

      icm = NextItem(icm);
    }

  /***************************************************************************/

  /***** Misc. stuff *********************************************************/

  /* Check levels */

  if ((long)RDB[DATA_GEOM_LEVELS] > max)
    {
      /* Check if detectors or sources are defined */

      if (((long)RDB[DATA_PTR_DET0] > VALID_PTR) ||
          ((long)RDB[DATA_PTR_SRC0] > VALID_PTR))
        Die(FUNCTION_NAME,
            "Geometry is disjoint, check source and detector definitions");

      Die(FUNCTION_NAME, "Geometry is disjoint");
    }

  /* Close list */

  lvl = (long)RDB[DATA_PTR_LVL0];
  CloseList(lvl);

  /* Delta-tracking enforce flag */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_DT_ENFORCE_NEXT_TRACK] = (double)ptr;

  /* Surface-tracking mode flag */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_ST_USE_STL_MODE] = (double)ptr;

  /* Allocate memory for undefined cell */

  cell = NewItem(DATA_PTR_UNDEFINED_CELL, CELL_BLOCK_SIZE);

  /* Put name and material */

  WDB[cell + CELL_PTR_NAME] = (double)PutText("undefined");
  WDB[cell + CELL_PTR_MAT] = NULLPTR;
  WDB[cell + CELL_TYPE] = CELL_TYPE_VOID;

  /* Allocate memory for search list option */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_CELL_SEARCH_LIST] = (double)ptr;

  /* Exit OK */

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
