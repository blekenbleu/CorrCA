/*----------------------------------------------------------------------------

  Image data type and basic functions.

  Copyright 2010-2011 rafael grompone von gioi (grompone@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------- Image Data Types -----------------------------*/
/*----------------------------------------------------------------------------*/
/** @file image.c
	Image data types.
	@author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include "image.h"
#include "spline.h"

/*----------------------------------------------------------------------------*/
/** Free memory used in image_char 'i'.
 */
void free_image_char(image_char i)
{
	if (i != NULL && i->data != NULL)
	{
		free((void*)i->data);
		free((void*)i);
	} else error("free_image_char: invalid input image.");
}

/*----------------------------------------------------------------------------*/
/** Create a new image_char of size 'xsize' times 'ysize'.
 */
void new_image_char(image_char &image, unsigned int xsize, unsigned int ysize)
{
  /* check parameters */
  if( xsize == 0 || ysize == 0 ) error("new_image_char: invalid image size.");

  /* get memory */
  if(image = (image_char)malloc(sizeof(struct image_char_s)))
  {
	image->data = (unsigned char *)calloc((size_t)(xsize*ysize),
										  sizeof(unsigned char));
  	/* set image size */
	image->xsize = xsize;
	image->ysize = ysize;
  }
  else error("not enough memory.");
}

image_char new_image_char(unsigned int xsize, unsigned int ysize)
{
	image_char image;
	new_image_char(image, xsize, ysize);
	return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_char of size 'xsize' times 'ysize',
	initialized to the value 'fill_value'.
 */
image_char new_image_char_ini(unsigned int xsize, unsigned int ysize,
							  unsigned char fill_value)
{
	image_char image;
	new_image_char_ini(image, xsize, ysize, fill_value);
	return image;
}

void new_image_char_ini(image_char &image, unsigned int xsize, unsigned int ysize,
						unsigned char fill_value)
{

  if (image = new_image_char(xsize, ysize)) {	/* create image */
	unsigned char *d = image->data;		/* initialize */
	for(unsigned char *dx = d + xsize*ysize; d < dx; d++)
		*d = fill_value;
  }
}

/*----------------------------------------------------------------------------*/
/** Create a new image_char, copy of image 'in'.
 */
image_char new_image_char_copy(image_char in)
{
  image_char image = new_image_char(in->xsize,in->ysize); /* create image */
  unsigned int N = in->xsize * in->ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = in->data[i];

  return image;
}

unsigned char image_char_pixel(image_char &image, unsigned int row, unsigned int col)
{
	return image->data[col + row * image->xsize];
}

/*----------------------------------------------------------------------------*/
/** Free memory used in image_int 'i'.
 */
void free_image_int(image_int i)
{
  if(i != NULL && i->data != NULL )
  {
	free((void *) i->data);
	free((void *) i);
  } else error("free_image_int: invalid input image.");
}

/*----------------------------------------------------------------------------*/
/** Create a new image_int of size 'xsize' times 'ysize'.
 */
image_int new_image_int(unsigned int xsize, unsigned int ysize)
{
  /* check parameters */
  if(xsize == 0 || ysize == 0) error("new_image_int: invalid image size.");

  /* get memory */
  image_int image = (image_int) malloc(sizeof(struct image_int_s));
  if(image != NULL )
  {
	image->data = (int *)calloc( (size_t)(xsize*ysize), sizeof(int) );

	/* set image size */
	image->xsize = xsize;
	image->ysize = ysize;
  } else error("not enough memory.");
  return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_int of size 'xsize' times 'ysize',
	initialized to the value 'fill_value'.
 */
image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
							 int fill_value )
{
  image_int image = new_image_int(xsize,ysize); /* create image */
  unsigned int N = xsize*ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = fill_value;

  return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_int, copy of image 'in'.
 */
image_int new_image_int_copy(image_int in)
{
  image_int image = new_image_int(in->xsize,in->ysize); /* create image */
  unsigned int N = in->xsize * in->ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = in->data[i];

  return image;
}

/*----------------------------------------------------------------------------*/
/** Free memory used in image_double 'i'.
 */
void free_image_double(image_double i)
{
  if(i != NULL && i->data != NULL )
  {
	free((void *) i->data);
	free((void *) i);
  } else error("free_image_double: invalid input image.");
}

/*----------------------------------------------------------------------------*/
/** Create a new image_double of size 'xsize' times 'ysize'.
 */
void new_image_double(image_double &image, unsigned int xsize, unsigned int ysize)
{
  image = {};

  /* check parameters */
  if( xsize == 0 || ysize == 0 ) error("new_image_double: invalid image size.");

  /* get memory */
  image = (image_double) malloc( sizeof(struct image_double_s) );
  if(image != NULL)
  {
	/* set image size */
	image->xsize = xsize;
	image->ysize = ysize;
	image->data = (double *) calloc( (size_t) (xsize*ysize), sizeof(double) );
  	if( image->data == NULL )
	  error("not enough memory.");
  }
  else error("not enough memory.");
}

/*----------------------------------------------------------------------------*/
/** Create a new image_double of size 'xsize' times 'ysize',
	initialized to the value 'fill_value'.
 */
void new_image_double_ini(image_double &image, unsigned int xsize,
							unsigned int ysize, double fill_value)
{
  new_image_double(image, xsize,ysize); /* create image */

  /* initialize */
  double *d = image->data;
  for(double *dx = d + xsize*ysize; d < dx; d++) *d = fill_value;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_double, copy of image 'in'.
 */
void new_image_double_copy(image_double &image, image_double in)
{
  new_image_double(image, in->xsize,in->ysize); /* create image */
  unsigned int N = in->xsize * in->ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = in->data[i];
}

void new_image_double_copy(image_double &image, image_char in)
{
  new_image_double(image, in->xsize,in->ysize); /* create image */
  unsigned int N = in->xsize * in->ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = in->data[i];
}

/*----------------------------------------------------------------------------*/
/** Linear Interpolation.
 */
double interpolate_image_double(image_double& in, int order, double u, double v)
{
	double color;
	if (!interpolate_spline(in, order, u, v, color)) 
		color = 0;
	return color;
}

/*----------------------------------------------------------------------------*/
/** Is pixel inside image?
 */
bool valid_image_double(image_double& in, int x, int y)
{
	return (0 <= x && x < (int)in->xsize && 0 <= y && y < (int)in->ysize);
}

/*----------------------------------------------------------------------------*/
