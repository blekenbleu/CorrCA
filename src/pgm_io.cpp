/*----------------------------------------------------------------------------

  PGM image file format I/O functions.

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
/*------------------------------ PGM image I/O -------------------------------*/
/*----------------------------------------------------------------------------*/
/** @file pgm_io.c
	PGM image format read and write functions.
	@author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "misc.h"
#include "image.h"
#include "pgm_io.h"

/*----------------------------------------------------------------------------*/
/** Skip white characters and comments in a PGM file.
 */
static void skip_whites_and_comments(FILE * f)
{
  int c;
  do
	{
	  while(isspace(c=getc(f))); /* skip spaces */
	  if(c=='#') /* skip comments */
		while( c!='\n' && c!='\r' && c!=EOF )
		  c=getc(f);
	}
  while( c == '#' || isspace(c) );
  if( c != EOF && ungetc(c,f) == EOF )
	error("unable to 'ungetc' while reading PGM file.");
}

/*----------------------------------------------------------------------------*/
/** Read a ASCII number from a PGM file.
 */
static unsigned int get_num(FILE * f)
{
  unsigned int num;
  int c;

  while(isspace(c=getc(f)));
  if(!isdigit(c))
	error("corrupted PGM file.");
  num = (unsigned int) (c - '0');
  while( isdigit(c=getc(f)) ) num = 10 * num + c - '0';
  if( c != EOF && ungetc(c,f) == EOF )
	error("unable to 'ungetc' while reading PGM file.");

  return num;
}

static unsigned int next_value(FILE *f)
{
	skip_whites_and_comments(f);
	return get_num(f);
}

/*----------------------------------------------------------------------------*/
/** Parse a PNM header
 * eliminate duplicate code in file readers
 */
FILE *read_pnm_header(char * name, unsigned int &x, unsigned int &y, int &bin, char &type)
{
  FILE *fin = (0 == strcmp(name, "-")) ? stdin : fopen(name,"rb");

  if (NULL == fin)
  {
	  error("Input image file %s open failed", name);
	  return NULL;
  }

  if(getc(fin) != 'P')		// first char in PNM
	error("not a PNM file!");

  if ((type = getc(fin)) == '2' || '3' == type)
	bin = FALSE;
  else if(type == '5' || '6' == type)
	bin = TRUE;
  else error("not a supported PNM file!");

  x = next_value(fin);			// columns
  y = next_value(fin);			// rows
  if (255 < next_value(fin))	// pixel range
	fprintf(stderr, "Warning: some values out of char range\n");

  if(!isspace(getc(fin)))		// '\n' expected
	error("missing newline after PNM file header.");

  return fin;
}

/*----------------------------------------------------------------------------*/
/** Read a PPM file into an "image_char".
	If the name is "-" the file is read from standard input.
 */
image_char read_ppm_image_char(char *fin)
{
  image_char chimg;
  int bin;
  char type;
  unsigned int columns, rows, width, i;

  FILE *f = read_pnm_header(fin, columns, rows, bin, type);

  if (NULL == fin)
    return NULL;
  if('3' != type && '6' != type)
    error("not a PPM file!");

  width = 3 * columns;	// 3 bytes per pixel
  chimg = new_image_char(width, rows);
  chimg->xsize = columns;	// 3 bytes per pixel

  for(i=0; i < rows * width; i++)
	  chimg->data[i] = bin ? (unsigned char)getc(f) : (unsigned char)get_num(f);

  if(f != stdin && fclose(f) == EOF)
	  error("unable to close PPM file %s after reading.", fin);

  return chimg;
}

/*----------------------------------------------------------------------------*/
/** Read a PPM file into an "image_double".
	If the name is "-" the file is read from standard input.
 */
void read_ppm_image_double(image_double& imageR, image_double& imageG, image_double& imageB,
							FILE *f, int bin, unsigned int xsize, unsigned int ysize)
{
  int c, g = 0;
  unsigned int x, y;

  imageR = new_image_double(xsize, ysize);
  imageG = new_image_double(2 * xsize, 2 *ysize);
  imageB = new_image_double(xsize, ysize);

  /* read data */
  if (bin)
	for(y = 0; y < ysize; y++)
	{
	  c = y * xsize;
	  g = 4 * c;
	  for (int xend = c + xsize; c < xend; c++)
	  {
	  	imageR->data[c] = (double) getc(f);
	  	imageG->data[g] = (double) getc(f);
	  	imageB->data[c] = (double) getc(f);
		g += 2;
	  }
	  // duplicate last column;
	  imageG->data[g - 1] = imageG->data[g - 2];
	}
  else for(y = 0; y < ysize; y++)
  {
	  c = y * xsize;
	  g = 4 * c;
	  for (int xend = c + xsize; c < xend; c++)
	  {
	  	imageR->data[c] = (double) get_num(f);
	  	imageG->data[g] = (double) get_num(f);	// even columns in even rows
	  	imageB->data[c] = (double) get_num(f);
		g += 2;
	  }
	  imageG->data[g - 1] = imageG->data[g - 2];
  }

  if (stdin != f && EOF == fclose(f))
	  error("failed to close PPM file after reading.");

  // fill in G by interpolation, first horizontally
  for (y = 0; y < ysize; y++)
  {
	  g = 1 + 4 * y * xsize;	// odd columns in even rows
	  for (x = 1; x < xsize; x++)
	  {
		imageG->data[g] = 0.5 * (imageG->data[g - 1] + imageG->data[g + 1]);
		g += 2;
	  }
  }
  // duplicate last row
  c = g - 2 * xsize;
  for (x = 0; x < xsize; x++)
	imageG->data[g++] = imageG->data[c++];

  // interpolate vertically
  int b = 0;		// column 0, row 0
  g = 2 * xsize;	// column 0, row 1 (unpopulated)
  c = 2 * g;		// column 0, row 2
  for (y = 1; y < ysize; y++)
  {
  	for (x = 0; x < 2 * xsize; x++)
		imageG->data[g++] = 0.5 * (imageG->data[b++] + imageG->data[c++]);

	b = g;
	g = c;
	c += 2 * xsize;
  }
}

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_char".
    If the name is "-" the file is read from standard input.
 */
image_char read_pgm_image_char(char * name)
{
  FILE * f; 
  int bin = FALSE;
  char type;
  unsigned int xsize,ysize,x,y;
  image_char image;

  /* open file */
  if (NULL == (f = read_pnm_header(name, xsize, ysize, bin, type))) return NULL;

  if('2' != type && '5' != type)
	error("not a PGM file!");

  /* get memory */
  image = new_image_char(xsize,ysize);

  /* read data */
  for(y=0;y<ysize;y++)
    for(x=0;x<xsize;x++)
      image->data[ x + y * xsize ] = bin ? (unsigned char) getc(f)
                                         : (unsigned char) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
      error("unable to close file %s while reading PGM file.", name);

  return image;
}

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_int".
	If the name is "-" the file is read from standard input.
 */
image_int read_pgm_image_int(char * name)
{
  FILE * f;
  char c;
  int bin = FALSE;
  unsigned int xsize,ysize,x,y;
  image_int image;

  /* open file */
  if (NULL == (f = read_pnm_header(name, xsize, ysize, bin, c))) return NULL;

  if('2' != c && '5' != c)
  {
	error("not a PGM file!");
	return NULL;
  }

  /* get memory */
  image = new_image_int(xsize,ysize);

  /* read data */
  for(y=0;y<ysize;y++)
	for(x=0;x<xsize;x++)
	  image->data[ x + y * xsize ] = bin ? getc(f) : (int) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
	  error("unable to close file %s while reading PGM file.", name);

  return image;
}

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_double".
	If the name is "-" the file is read from standard input.
 */
image_double read_pgm_image_double(char * name)
{
  FILE * f;
  int bin=FALSE;
  char c;
  unsigned int xsize,ysize,x,y;
  image_double image;

  /* open file */
  if (NULL == (f = read_pnm_header(name, xsize, ysize, bin, c))) return NULL;

  if('2' != c && '5' != c)
	error("not a PGM file!");

  /* get memory */
  image = new_image_double(xsize,ysize);

  /* read data */
  for(x = ysize * xsize, y=0; y < x; y++)
	image->data[y] = bin ? (double) getc(f) : (double) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
	  error("unable to close file %s while reading PGM file.", name);

  return image;
}

static FILE *pnm_open(char *name, char type, unsigned int x, unsigned int y, unsigned int max)
{
	FILE *f = (0 == strcmp(name, "-")) ? stdout : fopen(name, ('5' == type || '6' == type) ? "wb" : "w");
	
  	if (NULL == f)
		error("unable to open output image file %s", name);
	else
	{	/* write PNM header */
		fprintf(f, "P%c\n", type);
  		fprintf(f, "# Generated by chromaberrat\n");
		fprintf(f, "%u %u\n%u\n", x, y, max);
	}

  	return f;
}

/*----------------------------------------------------------------------------*/
/** Write an "image_char" into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_char(image_char image, char * name)
{
  FILE *f = pnm_open(name, '2', image->xsize, image->ysize, 255);

  if (NULL == f)
	return;

  unsigned int x,y,n;

  /* write data */
  for(n=0,y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	{
		fprintf(f,"%d ",image->data[ x + y * image->xsize ]);
		if(++n==8)  /* lines should not be longer than 70 characters  */
		{
			fprintf(f,"\n");
			n = 0;
		}
	}

  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("unable to close PGM file %s after writing.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_int" into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int(image_int image, char * name)
{
  unsigned int x,y,n;
  int v,max,min;

  /* check min and max values */
  max = min = 0;
  for(y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	{
		v = image->data[ x + y * image->xsize ];
		if( v > max ) max = v;
		if( v < min ) min = v;
	}
  if (min < 0)
	  fprintf(stderr, "Warning: write_pgm_image_int: negative values in '%s'.\n", name);
  if (max > 65535)
	  fprintf(stderr, "Warning: write_pgm_image_int: values exceeding 65535 in '%s'.\n", name);

  /* open file */
  FILE *f = pnm_open(name, '2', image->xsize, image->ysize, max);

  if (NULL == f)
    return;

  /* write data */
  for(n=0,y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	{
		fprintf(f,"%d ",image->data[ x + y * image->xsize ]);
		if(++n==8)  /* lines should not be longer than 70 characters  */
		{
			fprintf(f,"\n");
			n = 0;
		}
	}

  /* close file if needed */
  if(f != stdout && fclose(f) == EOF)
	  error("unable to close file %s while writing PGM file.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_int" normalized to [0,255] into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int_normalized(image_int image, char * name)
{
  unsigned int x,y,n;
  int v,max,min;
  double factor;

  /* check min and max values */
  max = min = 0;
  for(y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	{
		v = image->data[ x + y * image->xsize ];
		if( v > max ) max = v;
		if( v < min ) min = v;
	}
  factor = (0 == (max-min)) ?  1.0 : 255.0 / (double) (max-min);

  /* open file */
  FILE *f = pnm_open(name, '2', image->xsize, image->ysize, 255);

  if (NULL == f)
    return;

  /* write data */
  for(n = y = 0; y<image->ysize; y++)
	for(x = 0; x < image->xsize; x++)
	{
		v = (int)(factor * (image->data[ x + y * image->xsize ] - min));
		fprintf(f,"%d ",v);
		if(++n == 8)  /* lines should not be longer than 70 characters  */
		{
			fprintf(f,"\n");
			n = 0;
		}
	}

  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("unable to close PGM file %s after writing", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_double" into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double(image_double image, char * name)
{
  char *buffer = (char*)calloc(image->xsize, sizeof(char)), *cp = buffer;
  if (NULL == buffer)
	error("not enough memory.");

  /* open file */
  FILE *f = pnm_open(name, '5', image->xsize, image->ysize, 255);

  if (NULL == f)
    return;

  /* write data */
  double* id = image->data;
  for (size_t y = 0; y < image->ysize; y++)
  {
	cp = buffer;
	for (double *x = id + image->xsize; id < x; id++)
		*cp++ = (char)*id;
	fwrite(buffer, sizeof(char), image->xsize, f);
  }

  free(buffer);
  if( f != stdout && fclose(f) == EOF )	// close file if needed
	  error("unable to close PGM file %s after writing.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_double" normalized to [0,255] into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double_normalized(image_double image, char * name)
{
  double v,max,min,factor;
  unsigned int x,y,n;

  /* check min and max values */
  max = min = 0;
  for(y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	  {
		v = image->data[ x + y * image->xsize ];
		if( v > max ) max = v;
		if( v < min ) min = v;
	  }
  factor = ((max-min) == 0.0) ? 1.0 : 255.0 / (max-min);

  /* open file */
  FILE *f = pnm_open(name, '2', image->xsize, image->ysize, 255);

  if (NULL == f)
    return;

  /* write data */
  for(n=0,y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	  {
		v = factor * (image->data[ x + y * image->xsize ] - min);
		fprintf(f,"%d ",(unsigned char) v);
		if(++n==8)  /* lines should not be longer than 70 characters  */
		  {
			fprintf(f,"\n");
			n = 0;
		  }
	  }

  /* close file if needed */
  if(f != stdout && fclose(f) == EOF)
	  error("unable to close file %s while writing PGM file.", name);
}

/*----------------------------------------------------------------------------*/
/** Write R,G,B "image_double" into a PPM file.
	If the name is "-" the file is written to standard output.
	imageG width is 3x imageR, imageB
 */
void write_ppm_image_double(image_double imageR, image_double imageG, image_double imageB, char * name)
{
  size_t length = 3 * imageR->xsize;
  char *buffer = (char*)calloc(length, sizeof(char));
  // green plane may be same size or 2x red
  size_t bump = imageG->xsize / imageR->xsize, glen = imageG->xsize * (bump - 1);

  if (NULL == buffer)
  {
	error("not enough memory.");
	return;
  }

  /* open file */
  FILE *f = pnm_open(name, '6', imageR->xsize, imageR->ysize, 255);

  if (NULL == f)
    return;

  /* write data */
  double* r = imageR->data, *g = imageG->data, *b = imageB->data;
  char *s = buffer + length;
  for(size_t y = 0; y < imageR->ysize; y++)
  {
	for(char *cp = buffer; cp < s;)
	{
	  *cp++ = (char)(0.5 + *r++);
	  *cp++ = (char)(0.5 + *g);
	  *cp++ = (char)(0.5 + *b++);
	  g += bump;
	}
	g += glen;
	fwrite(buffer, sizeof(char), length, f);
  }

  free(buffer);
  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("failed to close PPM file %s after writing.", name);
}

/*----------------------------------------------------------------------------*/
