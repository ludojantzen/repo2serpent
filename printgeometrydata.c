/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printgeometrydata.c                            */
/*                                                                           */
/* Created:       2012/09/06 (JLe)                                           */
/* Last modified: 2012/09/06 (JLe)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: Prints universes to output                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintGeometryData:"

/*****************************************************************************/

void PrintGeometryData()
{
  long table, uni, lvl, lst, cell, mat, ptr;
  char outfile[MAX_STR];
  FILE *fp;

  /*********************************************************************/

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Reset table counter (TODO: lue to DATA:sta) */

  table = 0;

  /* Set output file */

  sprintf(outfile, "%s.out", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for writing */

  if ((fp = fopen(outfile, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Skip this for now */

  fclose(fp);
  return;

  /***************************************************************************/

  /***** Universe data *******************************************************/

  table++;

  fprintf(fp, "\n --- Table %2ld: Summary of Geometry and universes: \n\n", 
          table);

  /* Print level data */

  fprintf(fp, "The geometry consists of %ld levels:\n\n",
          (long)RDB[DATA_GEOM_LEVELS]);
  
  lvl = (long)RDB[DATA_PTR_LVL0];
  while (lvl > VALID_PTR)
    {
      /* Print */
      
      fprintf(fp, "Level %ld size: max %ld zones\n",
              (long)RDB[lvl + LVL_NUMBER], 
              (long)WDB[lvl + LVL_MAX_REGIONS]);
      
      /* Next level */
      
      lvl = NextItem(lvl);
    }

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      fprintf(fp, "\nUniverse \"%s\":\n\n", GetText(uni + UNIVERSE_PTR_NAME));

      fprintf(fp, " - Universe is at level %ld\n", 
              (long)RDB[uni + UNIVERSE_LEVEL]);

      /* Check type */

      if (((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_CELL) ||
          ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_SUPER))
        {
          if ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_CELL)         
            {
              fprintf(fp, " - %ldD Universe with maximum dimensions:\n\n", 
                      (long)RDB[uni + UNIVERSE_DIM]);
              
              fprintf (fp, "   x : [%12.5E %12.5E]\n", RDB[uni + UNIVERSE_MINX],
                       RDB[uni + UNIVERSE_MAXX]);
              fprintf (fp, "   y : [%12.5E %12.5E]\n", RDB[uni + UNIVERSE_MINY],
                       RDB[uni + UNIVERSE_MAXY]);
              
              if ((long)RDB[uni + UNIVERSE_DIM] == 3)
                fprintf (fp, "   z : [%12.5E %12.5E]\n", RDB[uni + UNIVERSE_MINZ],
                         RDB[uni + UNIVERSE_MAXZ]);
            }
          else
            fprintf(fp, " - Super-imposed %ldD Universe\n", 
                    (long)RDB[uni + UNIVERSE_DIM]);

          /* Pointer to cell list */
        
          lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
          CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);  
          
          fprintf(fp, "\nUniverse consists %ld cells:\n\n",
                  ListSize(lst));

          fprintf(fp, "---------------------------------------------------------\n");
          fprintf(fp, "Cell       Content         Volume (calculated by Serpent)\n");
          fprintf(fp, "---------------------------------------------------------\n");
          
          /* Loop over list */
          
          while (lst > VALID_PTR)
            {
              /* Pointer to cell */
              
              cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
              
              fprintf (fp, "%-10s ", GetText(cell + CELL_PTR_NAME));
              
              /* Check type */
              
              if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_MAT)
                {
                  /* Pointer to material */
                  
                  mat = (long)RDB[cell + CELL_PTR_MAT];
                  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
                  
                  fprintf(fp, "%-15s ", GetText(mat + MATERIAL_PTR_NAME));
                }
              else if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_VOID)
                fprintf(fp, "%-15s ", "void");
              else if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
                fprintf(fp, "%-15s ", "outside");
              else if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_FILL)
                {
                  ptr = (long)RDB[cell + CELL_PTR_FILL];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  
                  fprintf(fp, "fill %-10s ", GetText(ptr + UNIVERSE_PTR_NAME));
                }

              fprintf(fp, "%11.5E ", RDB[cell + CELL_VOLUME]);
              
              
              fprintf(fp, "\n");
              
              /* Next cell */
              
              lst = NextItem(lst);
            }

          fprintf(fp, "---------------------------------------------------------\n");
        }

      /* Next universe */

      uni = NextItem(uni);
    }
  
  /***************************************************************************/

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
