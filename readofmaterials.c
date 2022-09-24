/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofmaterials.c                              */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/11/08 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface materials file        */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFMaterials:"

/*****************************************************************************/

void ReadOFMaterials(long ifc)
{
  char *line, tmpstr[MAX_STR], fname[MAX_STR];
  long n, dim[7], ptr, i, cgns, type, mat, marr;
  long nc, nmat, cell, *matarr, disttype;
  double value;
  FILE *fp;

  /* Get interface type */

  type = (long)RDB[ifc + IFC_TYPE];

  /* Check interface type */

  if (type == IFC_TYPE_OPENFOAM)
    {
      /* The whole geometry consists of a single material */

      /* Create only one geometry cell */

      cell = NewItem(ifc + IFC_PTR_GCELL_LIST, CELL_BLOCK_SIZE);

      /* Allocate memory for cell collision counter */

      AllocValuePair(cell + CELL_COL_COUNT);

      /* Allocate memory for collision tet-cell */

      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[cell + CELL_PTR_PREV_TET] = (double)ptr;

      /* Put name for cell */

      WDB[cell + CELL_PTR_NAME] = (double)PutText("ifc cell");

      /* Get pointer to interface material array */

      marr = (long)RDB[ifc + IFC_PTR_MAT_ARR];

      /* Get pointer to the first interface material name in */
      /* the ASCII-block */

      mat = (long)RDB[marr + 0];

      /* Put negative pointer to cell material (indicates not linked) */
      /* Will be linked in processifctetmesh.c */

      WDB[cell + CELL_PTR_MAT] = -mat;

      /* Close the cell list */

      CloseList(cell);

      /* Link tet-cells to this geometry cell */

      /* Loop over cells */

      cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];

      while (cgns > VALID_PTR)
        {

          /* Put cell pointer to tet */

          WDB[cgns + IFC_TET_PRNT_PTR_CELL] = (double)cell;

          /* Next */

          cgns = NextItem(cgns);
        }

      return;
    }

  /********************************************/

  /* Different cells have different materials */

  /********************************************/

  /* Get materials file filename */

  sprintf(fname, "%s", GetText(ifc + IFC_PTR_OF_MFILE));

  /* Check file format */

  TestDOSFile(fname);

  /* Open materials file */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(ifc, "Materials file \"%s\" does not exist",
          fname);

  /* Read header data */

  ReadOFHeader(fp, &n, &i, (long *)dim, &disttype, &value);

  /* Get number of cells */

  nc = (long)RDB[ifc + IFC_NC];

  /* Check size */

  if (i != nc)
    Error(ifc, "Invalid number of entries in material file");

  /* Create an array for unique material names */

  matarr = (long *)Mem(MEM_ALLOC, 1, sizeof(long));

  nmat = 0;

  /* Loop over cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
  while (cgns > VALID_PTR)
    {
      /* Read material name */

      line = ReadOFData(fp, OF_FILE_MATERIAL);

      if (sscanf(line, "%s", tmpstr) == EOF)
        Error(ifc, "Not enough entries in materials file");

      /* Check for filled cell */

      if (!strcasecmp(tmpstr, "fill"))
        Error(ifc, "Filled cells not allowed with unstructured meshes");

      /* Check if the material name has been stored */

      for (i = 0; i < nmat; i++)
        {
          if (strcmp(&ASCII[matarr[i]], tmpstr) == 0)
            break;
        }

      if (i >= nmat)
        {

          /* New material */

          /* Reallocate material name array */

          matarr = (long*)Mem(MEM_REALLOC, matarr, (nmat+1)*sizeof(long));

          /* Put pointer to material name */

          matarr[nmat] = PutText(tmpstr);

          nmat++;
        }

      WDB[cgns + IFC_TET_PRNT_PTR_CELL] = (double)i;

      /* Next */

      cgns = NextItem(cgns);
    }

  /* Close file */

  fclose(fp);

  /* Reset cell pointer */

  cell = -1;

  /* Create nmat cells */

  for (i = 0; i < nmat; i++)
    {
      /* Create cell */

      cell = NewItem(ifc + IFC_PTR_GCELL_LIST, CELL_BLOCK_SIZE);

      /* Allocate memory for cell collision counter */

      AllocValuePair(cell + CELL_COL_COUNT);

      /* Allocate memory for collision tet-cell */

      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[cell + CELL_PTR_PREV_TET] = (double)ptr;

      /* Put name for cell */
      /* NB: GetText can be used only with ascii pointers */
      /* in RDB[], that's why we have to first write the pointer */
      /* from matarr[] to CELL_PTR_NAME */

      WDB[cell + CELL_PTR_NAME] = (double)matarr[i];
      sprintf(tmpstr, "ifc cell %s", GetText(cell + CELL_PTR_NAME));
      WDB[cell + CELL_PTR_NAME] = (double)PutText(tmpstr);

      /* Put material name */
      /* Linked in processifctetmesh.c */

      WDB[cell + CELL_PTR_MAT] = -(double)matarr[i];

    }

  CloseList(cell);

  /* Link geometry cells */

  /* Loop over cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
  while (cgns > VALID_PTR)
    {

      /* Get material number from tet */

      i = (long)RDB[cgns + IFC_TET_PRNT_PTR_CELL];

      /* Get cell pointer from cell list */

      cell = ListPtr(cell, i);

      /* Put cell pointer to tet */

      WDB[cgns + IFC_TET_PRNT_PTR_CELL] = (double)cell;

      /* Next */

      cgns = NextItem(cgns);
    }


  /* One of the materials is already known */

  nmat--;

  /* Store the number of "extra" materials */
  /* This is needed for geometry plotter   */

  WDB[DATA_OF_N_EXTRAMAT] += (double)nmat;

  /* Free the material array*/

  Mem(MEM_FREE,matarr);

}
