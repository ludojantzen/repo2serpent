/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : surfacesrc.c                                   */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Samples neutron on a surface                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SurfaceSrc:"

/*****************************************************************************/

void SurfaceSrc(long src, long surf, double *x, double *y, double *z, 
		double *u, double *v, double *w, long id)
{
  long ptr, type;
  double dir, x0, y0, z0, r, l;

  /* Check surface pointer */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get direction */

  dir = RDB[src + SRC_SURF_SIDE];

  /* Pointer to parameter list */
	      
  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];

  /* Sample coordinates and direction accordin to type  */
      
  switch((type = (long)RDB[surf + SURFACE_TYPE]))
    {
      /***********************************************************************/

    case SURF_CYL:
      {
	x0 = RDB[ptr];
	y0 = RDB[ptr + 1];
	r = RDB[ptr + 2];

	/* Sample direction vector */
	
	do
	  {
	    *u = 1.0 - 2*RandF(id);
	    *v = 1.0 - 2*RandF(id);
	  }
	while ((l = sqrt(*u**u + *v**v)) > 1.0);

	/* Normalize */

	*u = *u/l;
	*v = *v/l;
	*w = 0.0;

	/* Put coordinates */

	*x = x0 + r*(*u);
	*y = y0 + r*(*v);

	/* Put direction */

	*u = dir*(*u);
	*v = dir*(*v);

	break;
      }

      /***********************************************************************/

    case SURF_SPH:
      {
	x0 = RDB[ptr];
	y0 = RDB[ptr + 1];
	z0 = RDB[ptr + 2];
	r = RDB[ptr + 3];

	/* Sample direction vector */

	IsotropicDirection(u, v, w, id);

	/* Put coordinates */

	*x = x0 + r*(*u);
	*y = y0 + r*(*v);
	*z = z0 + r*(*w);

	/* Put direction */

	*u = dir*(*u);
	*v = dir*(*v);
	*w = dir*(*w);

	break;
      }

      /***********************************************************************/

    default:
      {
	Error(src, "Invalid surface type for source definition (surface %s)",
	      GetText(surf + SURFACE_PTR_NAME));
      }
    }
}

/*****************************************************************************/
