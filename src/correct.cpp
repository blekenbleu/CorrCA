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
template <typename T>
void matrix_shift(double &dy, double &dx, matrix<T> &coef, int color, double y, double x)
{
	// y and x are in range [0:1]
	int cc = coef.ncol();
	int c = color * cc;
	T x2 = x*x, y2 = y*y;
	dx = coef(c++) + x*coef(c++) + y*coef(c++) + x2*coef(c++) + y2*coef(c++)
	   + x*x2*coef(c++) + y*y2*coef(c++) + x*y*coef(c++) + y*x2*coef(c++) + x*y2*coef(c++);
	dy = coef(c++) + x*coef(c++) + y*coef(c++) + x2*coef(c++) + y2*coef(c++)
	   + x*x2*coef(c++) + y*y2*coef(c++) + x*y*coef(c++) + y*x2*coef(c++) + x*y2*coef(c++);
}

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
// !!! do NOT play with plane.data() elsewhere
#if 1
template <typename T>
static unsigned char get_CR(matrix<T> &plane, double row, double column)
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
	int p[4]{}, t = 0;		// for plane(p[i])

	 if (1 == Rfloor)
	 {
		p[0] = Cfloor;
		p[1] = p[0];
		p[2] = p[1] + plane.ncol();
		p[3] = p[2] + plane.ncol();
	} else if (Rmax == Rfloor - 1) {
		p[3] = Cfloor + Rmax * plane.ncol();
		p[2] = p[3];
		p[1] = p[2] - plane.ncol();
		p[0] = p[1] - plane.ncol();
	}
	if (1 == Cfloor) {
		t = 1;
		buffer[0] = plane(p[0]);
		buffer[4] = plane(p[1]);
		buffer[8] = plane(p[2]);
		buffer[12] = plane(p[3]);
		for (int i = 0; i < 3; i++)
		{
			buffer[1 + i] = plane(p[0]+i);
			buffer[5 + i] = plane(p[1]+i);
			buffer[9 + i] = plane(p[2]+i);
			buffer[13+ i] = plane(p[3]+i);
		}
		l[0] = buffer; l[1] = buffer + 4;
		l[2] = buffer + 8; l[3] = buffer + 12;
	}
	else if (Cmax == Cfloor + 1) {
		t = 1;
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
	if (t)
    	for (int i = 0; i < 4; i++)
			y[i] = CatmullRom(l[i][0], l[i][1], l[i][2], l[i][3], x);
	else for (int i = 0; i < 4; i++)
		y[i] = CatmullRom(plane(p[i]), plane(p[i]+1), plane(p[i]+2), plane(p[i]+3), x);
	x = CatmullRom(y[0], y[1], y[2], y[3], column - Cfloor);
	if (0 > x)
		x = 0;
	else if (255 < x)
		x = 255;
	return (unsigned char)(0.5 + x);
}

#else

// interpolate a pixel at floating point row and column in plane
// 4 interpolations along successive rows, then a columnar interpolation
static unsigned char get_CR(double row, double column, matrix<double> &plane)
{
	// handle borders
	int ncol = plane.ncol()
	uint Rmax = plane.nrow() - 1, Cmax = ncol - 1;
	int Rfloor = (int)floor(row), Cfloor = (int)floor(column);
	if (row <= 0)
		return (unsigned char)(0.5 + (Cfloor <= 0) ? plane(0)
				: column >= Cmax ? plane(Cmax)
				: plane(Cfloor)
				+ (column - Cfloor)
				* (plane(1 + Cfloor) - plane(Cfloor)));
	else if (row >= Rmax)
		return (unsigned char)(0.5 + (Cfloor <= 0) ? plane(Rmax,0)
				: column >= Cmax ? plane(Rmax, Cmax)
				: plane(Rmax, Cfloor)
				+ (column - Cfloor)
				* (plane(Rmax, 1 + Cfloor) - plane(Rmax, Cfloor)));
	else if (column >= Cmax)
		return (unsigned char)(0.5 + plane(Rfloor, Cmax)
				+ (row - Rfloor)
				* (plane(1 + Rfloor, Cmax) - plane(Rfloor, Cmax)));
	else if (column <= 0)
		return (unsigned char)(0.5 + plane(Rfloor, 0)
				+ (row - Rfloor)
				* (plane(1 + Rfloor, 0) - plane(Rfloor, 0)));

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
	int i, j, l0, l1, l2, l3;

	if (0 == Rfloor)						// duplicate first row	
	{
	  if (0 == Cfloor) {					// duplicate first row and column
		l0 = 0; l1 = ncol; l2 = l1 + ncol;
		for (i = 0; i < 3; i++)		
			y[i + 1] = CatmullRom(d = plane(l0+i), d, plane(l1+i), plane(l2+i), x); 
	  }
	  else if (Cmax == Cfloor + 1) {		// duplicate first row, last column
		l0 = Cfloor - 1; l1 = l0 + ncol; l2 = l1 + ncol;
		for (int i = 0; i < 3; i++)
			y[1 + i] = CatmullRom(plane(l0+i), plane(l1+i), d = plane(l2+i), d, x);
	  }
	  else for (int i = 0; i < 3; i++)
		y[1 + i] = CatmullRom(plane.data(Cfloor + i * ncol - 1), x);

	  y[0] = y[1];
	}
	else if (Rmax == Rfloor - 1) {	// duplicate top row
	  j = (Rmax - 2) * ncol;
	  if (0 == Cfloor) {				// duplicate first column and top row
		for (int i = 0; i < 3; i++) {
			y[i] = CatmullRom(d = plane(j), d, plane(1 + j), plane(2 + j), x);
			j += ncol;
		}
	  }
	  else if (Cmax == Cfloor + 1) {	// duplicate last column and top row
		j += Cmax - 2;
		for (int i = 0; i < 3; i++)
		{
			y[i] = CatmullRom(plane(j), plane(1 + j), d = plane(2 + j), d, x);
			j += ncol;
		} 
	  } else {
		j = Cfloor + (Rfloor - 1) * ncol - 1;
		for (i = 0; i < 3; i++) {
	 		y[i] = CatmullRom(plane.data(j), x);
			j += ncol;
		}
	  }		
	  y[3] = y[2];
	}
	else {							// 4x4 fully inside plane
		j = Cfloor + (Rfloor - 1) * ncol - 1;
		for (i = 0; i < 4; i++) {
			y[i] = CatmullRom(plane.data(j), x);
			j += ncol;
		}
	}

	// y[] is now populated
	x = CatmullRom(y[0], y[1], y[2], y[3], column - Cfloor);
	if (0 > x)
		return 0;
	if (255 < x)
		return 255;
	return (unsigned char)(0.5 + x);
}
#endif

static int read_coef(matrix<double> &coef, const char *fname)
{
	int rc = -1;
	if (4 > coef.nrow()) {
		printf("invalid matrix<double> coef\n");
		return rc;
	}
	std::ifstream f(fname);
	if(f) {
		char buf[151] { '\0' };
		f.getline(buf, 150);
		if ('#' == buf[0] && 0 == (rc = strncmp(2 + buf, "polyXR(x,y)", 11))) {
			int cs = coef.nrow() * coef.ncol();
			double d;
			char *more;
			for (int c = 0; c < cs && f; c++)
			{
				f.getline(buf, 150);
				if ('#' == *buf)
				{
					--c;
					continue;
				}
				coef(c) = d = strtod(buf, &more);
			}
		} else printf("invalid coef file %s\n", fname);
		f.close();
	} else printf("unable to read %s\n", fname);
	return rc;
}
