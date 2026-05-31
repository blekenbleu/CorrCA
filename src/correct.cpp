/*
#include <cmath>
#include <fstream>
#include "pgm_io.h"
#include "array.h"
using namespace libNumerics;

typedef unsigned int uint;
typedef unsigned char uchar; */

// included in gnuplot2file.h for visibility to array.h

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
static double a1(double ya, double yc) { return 0.5*(yc - ya); }
static double a2(double ya, double yb, double yc, double yd)
{ return 0.5*(2*ya - 5*yb + 4*yc - yd); }
static double a3(double ya, double yb, double yc, double yd)
{ return 0.5*(yd - ya + 3*yb - 3*yc); }
static double CatmullRom(double ya, double yb, double yc, double yd, double x)
{ return yb + a1(ya, yc)*x + a2(ya, yb, yc, yd)*x*x + a3(ya, yb, yc, yd)*x*x*x; }
static double CatmullRom(double *y, double x)
{ return y[1] + a1(y[0], y[2])*x + a2(y[0], y[1], y[2], y[3])*x*x + a3(y[0], y[1], y[2], y[3])*x*x*x; }

// interpolate a pixel at floating point row and column in plane
template <typename T>
static unsigned char get_CR(matrix<T> plane, int spline_order, double row, double column)
{
	// handle borders
	uint Rmax = plane.nrow() - 1, Cmax = plane.ncol() - 1;
	int Rfloor = (int)floor(row), Cfloor = (int)floor(column);
	if (row <= 0)
		return (Cfloor <= 0) ? plane(0)
				: column >= Cmax ? plane(0, Cmax)
				: (uchar)(0.5 + plane(0, Cfloor)
					+ (column - Cfloor)*plane(0, 1 + Cfloor)
					- plane(0, Cfloor));
	else if (row >= Rmax)
		return (Cfloor <= 0) ? plane(Rmax,0)
				: column >= Cmax ? plane(Rmax, Cmax)
				: (unsigned char)(0.5 + plane(Rmax, Cfloor)
					+ (column - Cfloor)*plane(Rmax, 1 + Cfloor)
					- plane(Rmax, Cfloor));
	else if (column >= Cmax)
		return (unsigned char)(0.5+plane(Rfloor, Cmax)
				+ (row - Rfloor)*(plane(1 + Rfloor, Cmax)
				- plane(Rfloor, Cmax)));
	else if (column <= 0)
		return (unsigned char)(0.5 + plane(Rfloor, 0)
				+ (row - Rfloor)*plane(1 + Rfloor, 0)
				- plane(Rfloor, 0));

	/* border pixels should now be handled
	 ; Catmull-Rom wants 4x4 pixel neighborhoods;
	 ; duplicate bottom or top rows, left- or right-most columns
	 ; for pixels in row 1 or Rmax - 1, column 1 or Cmax -1
     */
	unsigned char buffer[16] = {0}, *l[4] = { buffer, 4 + buffer, 8 + buffer, 12 + buffer }, i = 0;

	 if (1 == Rfloor)
	 {
		l[0] = plane.data(Cfloor);
		l[1] = l[0];
		l[2] = l[1] + plane.ncol();
		l[3] = l[2] + plane.ncol();
	} else if (Rmax == Rfloor - 1) {
		l[3] = plane.data(Cfloor + Rmax * plane.ncol());
		l[2] = l[3];
		l[1] = l[2] - plane.ncol();
		l[0] = l[1] - plane.ncol();
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

#if 0
/* interpolate a pixel at floating point row and column in plane
static unsigned char get_CR(double row, double column, matrix<double> plane)
{
	// handle borders
	uint Rmax = plane.nrow() - 1, Cmax = plane.ncol() - 1;
	int Rfloor = (int)floor(row), Cfloor = (int)floor(column);
	if (row <= 0)
		return (unsigned char)(0.5 + (Cfloor <= 0) ? plane(0,0)
				: column >= Cmax ? plane(0, Cmax)
				: plane(0, Cfloor)
					+ (column - Cfloor)*(plane(0, 1 + Cfloor)
					- plane(0, Cfloor)));
	else if (row >= Rmax)
		return (unsigned char)(0.5 + (Cfloor <= 0) ? plane(Rmax,0)
				: column >= Cmax ? plane(Rmax, Cmax)
				: plane(Rmax, Cfloor)
					+ (column - Cfloor)*(plane(Rmax, 1 + Cfloor)
					- plane(Rmax, Cfloor)));
	else if (column >= Cmax)
		return (unsigned char)(0.5+plane(Rfloor, Cmax)
				+ (row - Rfloor)*(plane(1 + Rfloor, Cmax)
				- plane(Rfloor, Cmax)));
	else if (column <= 0)
		return (unsigned char)(0.5 + plane(Rfloor, 0)
				+ (row - Rfloor)*(plane(1 + Rfloor, 0)
				- plane(Rfloor, 0)));

	/* border pixels should now be handled
	 ; Catmull-Rom wants 4x4 pixel neighborhoods;
	 ; duplicate bottom or top rows, left- or right-most columns
	 ; for pixels in row 1 or Rmax - 1, column 1 or Cmax -1
	 :
	 :	[Rmax,0] . . . . . . . . . . . . . . . . . . . . . . . . . . .[Rmax, Cmax]
	 :			. . . . . . . . . . . . . . . . . . . . . . . . . . .
	 :			. . . . . . . . . . . . . . . . . . . . . . . . . . .
	 :
	 :			[Rfloor + 1, Cfloor]		[Rfloor + 1, Cfloor + 1]
								[row, column]
	 :			[Rfloor, Cfloor]			[Rfloor, Cfloor + 1]
	 :			. . . . . . . . . . . . . . . . . . . . . . . . . . .
	 :			. . . . . . . . . . . . . . . . . . . . . . . . . . .
	 :
	 :	[0,0]   . . . . . . . . . . . . . . . . . . . . . . . . . . . [0, Cmax]
     */
	double d, x = row - Rfloor, y[4]{};
	int i, j, l[4]{};

	if (0 == Rfloor)						// duplicate first row	
	{
	  if (0 == Cfloor) {					// duplicate first row and column
		l[0] = 0; l[1] = plane.ncol(); l[2] = l[1] + plane.ncol();
		for (i = 0; i < 3; i++)		
			y[i + 1] = CatmullRom(d = plane(l[0]+i), d, plane(l[1]+i), plane(l[2]+i), x); 
	  } else if (Cmax == Cfloor + 1) {		// duplicate first row, last column
		l[0] = Cfloor - 2; l[1] = l[0] + plane.ncol(); l[2] = l[1] + plane.ncol();
		for (int i = 0; i < 3; i++)
			y[1 + i] = CatmullRom(plane(l[0]+i), plane(l[1]+i), d = plane(l[2]+i), d, x);
	  } else for (int i = 0; i < 3; i++)
			y[1 + i] = CatmullRom(plane.data(i * plane.ncol()), x);

	  y[0] = y[1];
	} else if (Rmax == Rfloor - 1) {	// duplicate top row
	  j = (Rmax - 2) * plane.ncol();
	  if (0 == Cfloor) {				// duplicate first column and top row
		for (int i = 0; i < 3; i++) {
			y[i] = CatmullRom(d = plane(j), d, plane(1 + j), plane(2 + j), x);
			j += plane.ncol();
		}
	  }
	  else if (Cmax == Cfloor + 1) {	// duplicate last column and top row
		j += Cmax - 2;
		for (int i = 0; i < 3; i++)
		{
			y[i] = CatmullRom(plane(j), plane(1 + j), d = plane(2 + j), d, x);
			j += plane.ncol();
		} 
	  } else {
		j = Cfloor + (Rfloor - 1) * plane.ncol() - 1;
		for (i = 0; i < 3; i++) {
	 		y[i] = CatmullRom(plane.data(j), x);
			j += plane.ncol();
		}
	  }		
	  y[3] = y[2];
	} else {							// 4x4 fully inside plane
		j = Cfloor + (Rfloor - 1) * plane.ncol() - 1;
		for (i = 0; i < 4; i++) {
			y[i] = CatmullRom(plane.data(j), x);
			j += plane.ncol();
		}
	}
	// y[] is now populated
	x = CatmullRom(y[0], y[1], y[2], y[3], row - Cfloor);
	if (0 > x)
		x = 0;
	else if (255 < x)
		x = 255;
	return (unsigned char)(0.5 + x);
} */
#endif
