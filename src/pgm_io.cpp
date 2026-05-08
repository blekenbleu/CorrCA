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

typedef unsigned int uint;

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
static FILE *read_pnm_header(char * name, unsigned int &x, unsigned int &y, int &bin, char &type)
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
/** Read a PPM file into an "image_double".
	Files named "-" are read from standard input.
	For compatibility with deBayer(),
	green planes must be twice red or blue width and height.
	Since deBayer() expects pixels in this matrix:
   rgrgr <- row 0;  unread pixels in lower case
   gbgbg
   rgRGR
   gbGBG
   rgRGR
   ... note position of "real" green pixels
   and duplicate scaling used in deBayer(),
   while special-case upper left and lower right corner green pixel components
   by averaging e.g. row 0 column 1 with row 1 column 0,
   so 1/4 of green pixels should have original Bayer values,
   while averaging other green values from 4 vertically and horizontally adjacent pixels.
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
/** Write R,G,B "image_double" into a PPM file.
	If the name is "-" the file is written to standard output.
	imageG width is 3x imageR, imageB
 */
void write_ppm_image_double(image_double imageR, image_double imageG, image_double imageB, char * name)
{
  size_t length = imageR->xsize; length *= 3;
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
  double *r = imageR->data, *g = 1 + imageG->data, *b = imageB->data;
  double *g2 = g + glen - 1;
  char *s = buffer + length;
  for(size_t y = 0; y < imageR->ysize; y++)
  {
	for(char *cp = buffer; cp < s;)
	{
	  *cp++ = (char)(0.5 + *r++);
	  *cp++ = (char)(0.5 + 0.5 * (*g + *g2));
	  *cp++ = (char)(0.5 + *b++);
	  g += bump;
	  g2 += bump;
	}
	g += glen;
	g2 += glen;
	fwrite(buffer, sizeof(char), length, f);
  }

  free(buffer);
  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("failed to close PPM file %s after writing.", name);

  printf(" imageR %dx%d, imageG %dx%d imageB %dx%d", imageR->xsize, imageR->ysize,
			imageG->xsize, imageG->ysize, imageB->xsize, imageB->ysize);
}

void read_pnm_double(image_double &imageR, image_double &imageG, image_double &imageB, char *fnameRGB)
{
	int bin;
	char type;
	uint wiG, heG;

	FILE* f = read_pnm_header(fnameRGB, wiG, heG, bin, type);

	if (NULL == f)
		error("NULL pnm FILE*");

	if ('6' == type || '3' == type)
	{
		read_ppm_image_double(imageR, imageG, imageB, f, bin, wiG, heG);
		wiG *= 2; heG *= 2;
	}
	else if('5' == type || '2' == type)
	{
		image_double image_bayer = read_pgm_image_double(fnameRGB);
		int wi = image_bayer->xsize, he = image_bayer->ysize;
		int wiRB = wi/2, heRB = he/2;
		wiG = wiRB*2; heG = heRB*2;
		imageR = new_image_double_ini(wiRB, heRB, 255);
		imageG = new_image_double_ini(wiG, heG, 255);
		imageB = new_image_double_ini(wiRB, heRB, 255);
		deBayer<double>(image_bayer, imageR, imageG, imageB);
		free_image_double(image_bayer);
	} else error("not a PNM file!\n");
  
	printf("imageR %dx%d imageG %dx%d imageB %dx%d\n", imageR->xsize,
			imageR->ysize, imageG->xsize, imageG->ysize, imageB->xsize, imageB->ysize);
	printf("\nSaving uncorrected.ppm\n");
	write_ppm_image_double(imageR, imageG, imageB, "R:/Temp/uncorrected.ppm");
/*
	exit(0);
 */
	printf("\nSaving uncorrected PGMs\n");
	write_pgm_image_double(imageR, "R:/Temp/uncorrectedR.pgm");
	write_pgm_image_double(imageG, "R:/Temp/uncorrectedG.pgm");
	write_pgm_image_double(imageB, "R:/Temp/uncorrectedB.pgm");
//	exit(0);
}
/*----------------------------------------------------------------------------*/
