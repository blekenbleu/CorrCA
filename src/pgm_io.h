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
/** @file pgm_io.h
    PGM image format read and write functions.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef PGMIO_HEADER
#define PGMIO_HEADER

#include "image.h"

/*----------------------------------------------------------------------------*/
/** Read a Bayer-interleaved PGM file into an "image_double".
    If the name is "-" the file is read from standard input.
 */

template <typename T>
void deBayer(image_double &img_bayer, image_double &imgR, image_double &imgG, image_double &imgB);
image_double read_pgm_image_double(char * name);

/*----------------------------------------------------------------------------*/
/** Read a PPM file into 3 "image_double".
    If the name is "-" the file is read from standard input.
 */
void read_ppm_image_double(image_double& imageR, image_double& imageG, image_double& imageB,
							FILE *f, int bin, unsigned int xsize, unsigned int ysize);

/*----------------------------------------------------------------------------*/
/** Read a PPM or Bayer PGM file into 3 "image_double";
    file name "-" is read from standard input.
 */
void read_pnm_double(image_double imageR, image_double imageG, image_double imageB, char *fnameRGB);

/*----------------------------------------------------------------------------*/
/** Write an "image_double" into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double(image_double image, char * name);
void write_ppm_image_double(image_double imageR, image_double imageG, image_double imageB, char * name);

#endif /* !PGMIO_HEADER */
/*----------------------------------------------------------------------------*/
