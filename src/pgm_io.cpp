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
  if(!isdigit(c)) error("corrupted PGM file.");
  num = (unsigned int) (c - '0');
  while( isdigit(c=getc(f)) ) num = 10 * num + c - '0';
  if( c != EOF && ungetc(c,f) == EOF )
	error("unable to 'ungetc' while reading PGM file.");

  return num;
}

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_char".
	If the name is "-" the file is read from standard input.
 */
image_char read_pgm_image_char(char * name)
{
  FILE * f;
  int c,bin=FALSE;
  unsigned int xsize,ysize,depth,x,y;
  image_char image;

  /* open file */
  if(strcmp(name,"-") == 0 ) f = stdin;
  else f = fopen(name,"rb");
  if (f == NULL)
  {
	  error("unable to open input image file %s", name);
	  return NULL;
  }
  /* read header */
  if( getc(f) != 'P' ) error("not a PGM file!");
  if( (c=getc(f)) == '2' ) bin = FALSE;
  else if( c == '5' ) bin = TRUE;
  else error("not a PGM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);			/* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);			/* Y size */
  skip_whites_and_comments(f);
  depth = get_num(f);			/* depth */
  if( depth > 255 ) fprintf(stderr,"Warning: some values out of char range\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("corrupted PGM file.");

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
/** Read a PPM file into an "image_char".
	If the name is "-" the file is read from standard input.
 */
image_char read_ppm_image_char(char * name)
{
  FILE * f;
  int c, bin=FALSE;
  unsigned int xsize, ysize, depth, x, y, width, i;
  image_char image;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdin;
  else f = fopen(name,"rb");
  if(f == NULL)
  {
	error("unable to open input PPM file %s", name);
	return NULL;
  }
  /* read header */
  if( getc(f) != 'P' )
	error("not a PNM file!");
  if((c=getc(f)) == '3') bin = FALSE;
  else if(c == '6') bin = TRUE;
  else error("not a PPM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);			/* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);			/* Y size */
  skip_whites_and_comments(f);
  depth = get_num(f);			/* depth */
  if( depth > 255 ) fprintf(stderr,"Warning: some values out of char range\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("corrupted PGM file.");

  /* get memory */
  width = 3 * xsize;	// 3 bytes per pixel
  image = new_image_char(width, ysize);
  image->xsize = xsize;	// 3 bytes per pixel

  /* read data */
  for(i=y=0; y < ysize; y++)
	for(x=0; x < width; x++)
	  image->data[i++] = bin ? (unsigned char) getc(f)
							 : (unsigned char) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
	  error("unable to close PPM file %s after reading.", name);

  return image;
}

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_int".
	If the name is "-" the file is read from standard input.
 */
image_int read_pgm_image_int(char * name)
{
  FILE * f;
  int c,bin=FALSE;
  unsigned int xsize,ysize,x,y;
  image_int image;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdin;
  else f = fopen(name,"rb");
  if( f == NULL ) error("unable to open input image file %s", name);

  /* read header */
  if( getc(f) != 'P' ) error("not a PGM file!");
  if( (c=getc(f)) == '2' ) bin = FALSE;
  else if( c == '5' ) bin = TRUE;
  else error("not a PGM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);			/* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);			/* Y size */
  skip_whites_and_comments(f);
  /*unsigned int depth =*/ get_num(f);			/* depth */
  /* white before data */
  if(!isspace(c=getc(f))) error("corrupted PGM file.");

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
  int c,bin=FALSE;
  unsigned int xsize,ysize,depth,x,y;
  image_double image;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdin;
  else
	  f = fopen(name,"rb");
  if( f == NULL )
	  error("unable to open input image file %s", name);

  /* read header */
  if( getc(f) != 'P' ) error("not a PGM file!");
  if( (c=getc(f)) == '2' ) bin = FALSE;
  else if( c == '5' ) bin = TRUE;
  else error("not a PGM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);			/* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);			/* Y size */
  skip_whites_and_comments(f);
  depth = get_num(f);			/* depth */
  if(depth==0) fprintf(stderr,"Warning: depth=0, probably invalid PGM file\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("corrupted PGM file.");

  /* get memory */
  image = new_image_double(xsize,ysize);

  /* read data */
  for(y=0;y<ysize;y++)
	for(x=0;x<xsize;x++)
	  image->data[ x + y * xsize ] = bin ? (double) getc(f)
										 : (double) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
	  error("unable to close file %s while reading PGM file.", name);

  return image;
}

/*----------------------------------------------------------------------------*/
/** Read a PPM file into an "image_double".
	If the name is "-" the file is read from standard input.
 */
void read_ppm_image_double(image_double& imageR, image_double& imageG, image_double& imageB, char * name)
{
  FILE * f;
  int c, g = 0, bin=FALSE;
  unsigned int xsize, ysize, depth, x, y;

  /* open file */
  if( strcmp(name, "-") == 0 ) f = stdin;
  else
	  f = fopen(name, "rb");
  if( f == NULL )
	  error("unable to open input image file %s", name);

  /* read header */
  if( getc(f) != 'P' ) error("not a PNM file!");
  if( (c=getc(f)) == '3' ) bin = FALSE;
  else if( c == '6' ) bin = TRUE;
  else error("not a PPM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);			/* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);			/* Y size */
  skip_whites_and_comments(f);
  depth = get_num(f);			/* depth */
  if(depth==0) fprintf(stderr, "Warning: depth=0, probably invalid PPM file\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("corrupted PPM file.");

  /* get memory */
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

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
	  error("failed to close PPM file %s after reading.", name);

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
/** Write an "image_char" into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_char(image_char image, char * name)
{
  FILE * f;
  unsigned int x,y,n;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdout;
  else f = fopen(name,"w");
  if (f == NULL)
  {
	  error("unable to open output image file %s", name);
	  return;
  }
  /* write header */
  fprintf(f,"P2\n");
  fprintf(f,"%u %u\n",image->xsize,image->ysize);
  fprintf(f,"255\n");

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
	  error("unable to close file %s while writing PGM file.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_int" into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int(image_int image, char * name)
{
  FILE * f;
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
  if(strcmp(name,"-") == 0)
	  f = stdout;
  else f = fopen(name,"w");
  if (f == NULL)
  {
	  error("unable to open output image file %s", name);
	  return;
  }
  /* write header */
  fprintf(f,"P2\n");
  fprintf(f,"%u %u\n",image->xsize,image->ysize);
  fprintf(f,"%d\n",max);

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
	  error("unable to close file %s while writing PGM file.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_int" normalized to [0,255] into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int_normalized(image_int image, char * name)
{
  FILE * f;
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
  if( (max-min) == 0 ) factor = 1.0;
  else factor = 255.0 / (double) (max-min);

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdout;
  else f = fopen(name,"w");
  if( f == NULL )
  {
	error("unable to open output image file %s", name);
	return;
  }
  /* write header */
  fprintf(f,"P2\n");
  fprintf(f,"%u %u\n",image->xsize,image->ysize);
  fprintf(f,"255\n");

  /* write data */
  for(n = y = 0; y<image->ysize; y++)
	for(x = 0; x < image->xsize; x++)
	{
		v = (int)(factor * (image->data[ x + y * image->xsize ] - min));
		fprintf(f,"%d ",v);
		if(++n==8)  /* lines should not be longer than 70 characters  */
		{
			fprintf(f,"\n");
			n = 0;
		}
	}

  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("unable to close file %s while writing PGM file", name);
}

/*----------------------------------------------------------------------------*/
/** Write R,G,B "image_double" into a PPM file.
	If the name is "-" the file is written to standard output.
	imageG width is 3x imageR, imageB
 */
void write_ppm_image_double(image_double imageR, image_double imageG, image_double imageB, char * name)
{
  FILE * f;
  unsigned int x,y,n,d,g;

  /* open file */
  if (strcmp(name,"-") == 0 ) f = stdout;
  else f = fopen(name,"wb");	// binary mode to prevent x0D0A newlines
  if( f == NULL )
  {
	error("unable to open output image file %s", name);
	return;
  }

  size_t length = (size_t)(3 * imageR->xsize);
  char *data = (char*)calloc(length, sizeof(char));
  char *dd = data;
  if (NULL == data)
  {
	error("not enough memory.");
	return;
  }

  /* write header */
  fprintf(f, "P6\n");	// binary ppm
  fprintf(f, "# Generated by chromaberrat\n");
  fprintf(f, "%u %u\n", imageR->xsize, imageR->ysize);
  fprintf(f, "255"); fputc(10, f);

  /* write data */
  for(g = d = n = y = 0; y < imageR->ysize; y++)
  {
	dd = data;
	for(x = 0; x < imageR->xsize; x++)
	{
	  *dd++ = (char)imageR->data[d];
	  *dd++ = (char)imageG->data[g++];
	  *dd++ = (char)imageB->data[d++];
/*
	  fprintf(f,"%d %d %d ",
				(unsigned int) imageR->data[d],
				(unsigned int) imageG->data[g++],
				(unsigned int) imageB->data[d++]);
	  if(++n == 5)  // lines should be < 70 characters
	  {
		fprintf(f,"\n");
		n = 0;
	  }
 */
	  g++;
	}
	g += imageG->xsize;
	fwrite(data, sizeof(char), length, f);
  }

  free(data);
  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("failed to close PPM file %s after writing.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_double" into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double(image_double image, char * name)
{
  FILE * f;
  unsigned int x,y,n;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdout;
  else f = fopen(name,"w");
  if( f == NULL )
  {
	error("unable to open output image file %s", name);
	return;
  }
  /* write header */
  fprintf(f,"P2\n");
  fprintf(f,"%u %u\n",image->xsize,image->ysize);
  fprintf(f,"255\n");

  /* write data */
  for(n=0,y=0; y<image->ysize; y++)
	for(x=0; x<image->xsize; x++)
	  {
		fprintf(f,"%d ",(unsigned int) image->data[ x + y * image->xsize ]);
		if(++n==8)  /* lines should not be longer than 70 characters  */
		  {
			fprintf(f,"\n");
			n = 0;
		  }
	  }

  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
	  error("unable to close file %s while writing PGM file.", name);
}

/*----------------------------------------------------------------------------*/
/** Write an "image_double" normalized to [0,255] into a PGM file.
	If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double_normalized(image_double image, char * name)
{
  FILE * f;
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
  if( (max-min) == 0.0 ) factor = 1.0;
  else factor = 255.0 / (max-min);

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdout;
  else f = fopen(name,"w");
  if (f == NULL)
  {
	  error("unable to open output image file %s", name);
	  return;
  }
  /* write header */
  fprintf(f,"P2\n");
  fprintf(f,"%u %u\n",image->xsize,image->ysize);
  fprintf(f,"255\n");

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
  if( f != stdout && fclose(f) == EOF )
	  error("unable to close file %s while writing PGM file.", name);
}

/*----------------------------------------------------------------------------*/
