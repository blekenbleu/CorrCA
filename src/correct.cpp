#include <cmath>
#include <fstream>
#include "pgm_io.h"
typedef unsigned int uint;
/* Shift Red and Blue pixel components
 ; according to polymomials previously calculated,
 ; with unique coefficient sets for each of Red dx, Red dy, Blue dx, Blue dy
 ; Reasonable surface fits over a test image wants double C[4][10] coefficients
 ; e.g. C[0][] for Red dx, using polynomial
 ; C[0][0] + C[0][1]*x + C[0][2]*y + C[0][3]*x*x + C[0][4]*y*y + C[0][5]*x*x*x
 ; + C[0][6]*y*y*y + C[0][5]*x*y + C[0][5]*x*x*y + C[0][5]*x*y*y
 ; ... where x and y are pixel coordinates, rescaled to range 0:1
 ; - calculate floating point row and column offsets for some input pixel[row, column]
 ; - resample Red and Blue values using 4x4 pixel neighborhoods around offsets
 */

// https://danceswithcode.net/engineeringnotes/interpolation/interpolation.html
static double CatmullRom(double ya, double yb, double yc, double yd, double x)
{
	double a0 = yb, a1 = 0.5*((double)yc - ya);
	double a2 = 0.5*((double)2*ya - 5*yb + 4*yc - yd);
	double a3 = 0.5*((double)yd - ya + 3*yb - 3*yc);

	return a0 + a1*x + a2*x*x + a3*x*x*x;
}

// interpolate a pixel at floating point row and column in plane
static unsigned char get_CR(double row, double column, image_char plane)
{
	// handle borders
	uint Rmax = plane->ysize - 1, Cmax = plane->xsize - 1;
	int Rfloor = (int)floor(row), Cfloor = (int)floor(column);
	if (row <= 0)
		return (Cfloor <= 0) ? image_char_pixel(plane, 0,0)
				: column >= Cmax ? image_char_pixel(plane, 0, Cmax)
				: (unsigned char)(image_char_pixel(plane, 0, Cfloor)
					+ (column - Cfloor)*(image_char_pixel(plane, 0, 1 + Cfloor)
					- image_char_pixel(plane, 0, Cfloor)));
	else if (row >= Rmax)
		return (Cfloor <= 0) ? image_char_pixel(plane, Rmax,0)
				: column >= Cmax ? image_char_pixel(plane, Rmax, Cmax)
				: (unsigned char)(0.5 + image_char_pixel(plane, Rmax, Cfloor)
					+ (column - Cfloor)*(image_char_pixel(plane, Rmax, 1 + Cfloor)
					- image_char_pixel(plane, Rmax, Cfloor)));
	else if (column >= Cmax)
		return (unsigned char)(0.5+image_char_pixel(plane, Rfloor, Cmax)
				+ (row - Rfloor)*(image_char_pixel(plane, 1 + Rfloor, Cmax)
				- image_char_pixel(plane, Rfloor, Cmax)));
	else if (column <= 0)
		return (unsigned char)(0.5 + image_char_pixel(plane, Rfloor, 0)
				+ (row - Rfloor)*(image_char_pixel(plane, 1 + Rfloor, 0)
				- image_char_pixel(plane, Rfloor, 0)));

	/* border pixels should now be handled
	 ; Catmull-Rom wants 4x4 pixel neighborhoods;
	 ; duplicate bottom or top rows, left- or right-most columns
	 ; for pixels in row 1 or Rmax - 1, column 1 or Cmax -1
     */
	unsigned char buffer[16]{}, *l[4]{}, i = 0;

	 if (1 == Rfloor)
	 {
		l[0] = plane->data + Cfloor;
		l[1] = l[0];
		l[2] = l[1] + plane->xsize;
		l[3] = l[2] + plane->xsize;
	} else if (Rmax == Rfloor - 1) {
		l[3] = plane->data + Cfloor + Rmax * plane->xsize;
		l[2] = l[3];
		l[1] = l[2] - plane->xsize;
		l[0] = l[1] - plane->xsize;
	}
	if (1 == Cfloor) {
		buffer[0] = l[0][0];
		buffer[4] = l[1][0];
		buffer[8] = l[2][0];
		buffer[12] = l[3][0];
		for (int i = 0; i < 3; i++)
		{
			buffer[1 + i] = l[0][i];
			buffer[5 + i] = l[1][i];
			buffer[9 + i] = l[2][i];
			buffer[13+ i] = l[3][i];
		}
		l[0] = buffer; l[1] = buffer + 4;
		l[2] = buffer + 8; l[3] = buffer + 12;
	}
	else if (Cmax == Cfloor + 1) {
		buffer[3] = l[0][3];
		buffer[7] = l[1][3];
		buffer[11] = l[2][3];
		buffer[15] = l[3][3];
		for (int i = 0; i < 3; i++)
		{
			buffer[i] = l[0][i];
			buffer[4 + i] = l[1][i];
			buffer[8 + i] = l[2][i];
			buffer[12 + i] = l[3][i];
		}
		l[0] = buffer; l[1] = buffer + 4;
		l[2] = buffer + 8; l[3] = buffer + 12;
	}
	// l[][] is now populated
	double x = row - Rfloor, y[4]{};
    for (int i = 0; i < 4; i++)
		y[i] = CatmullRom(l[i][0], l[i][1], l[i][2], l[i][3], x);
	x = CatmullRom(y[0], y[1], y[2], y[3], row - Cfloor);
	if (0 > x)
		x = 0;
	else if (255 < x)
		x = 255;
	return (unsigned char)(0.5 + x);
}
