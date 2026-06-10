#include "regress.h"

/* Fractionally shift Red and Blue pixel components using polymomials
 ; with polynomial coefficient sets unique to Red sxR, Red syR, Blue sxB, Blue syB
 ; Reasonable surface fits over a test image want polynomial matrix<double> coef(4, 11)
 ; with a row for each of sxR, syR, sxB, syB shifts
 ; - matrix_shift() calculates floating point sy and sx shifts for some input pixel[y, x]
 ;		with pixel coordinates y and x are rescaled to range 0:1
 ; - getCR() resamples Red and Blue values using 4x4 pixel neighborhoods around shifted pixels
 */
void matrix_shift(double &sy, double &sx, matrix<double> &coef, int color, double y, double x)
{
	// y and x are in range [0:1]
	int cc = coef.ncol();
	int c = color * cc;
	double x2 = x*x, y2 = y*y;

	// should match gnuplot equations
	sx = coef.get(c) + x*coef.get(c+1) + y*coef.get(c+2)
		+ x2*coef.get(c+3) + y2*coef.get(c+4) + x*x2*coef.get(c+5)
		+ y*y2*coef.get(c+6) + x*y*coef.get(c+7) + y*x2*coef.get(c+8)
		+ x*y2*coef.get(c+9) + x2*y2*coef.get(c+10);
	c += cc;
	sy = coef.get(c) + x*coef.get(c+1) + y*coef.get(c+2)
		+ x2*coef.get(c+3) + y2*coef.get(c+4) + x*x2*coef.get(c+5)
		+ y*y2*coef.get(c+6) + x*y*coef.get(c+7) + y*x2*coef.get(c+8)
		+ x*y2*coef.get(c+9) + x2*y2*coef.get(c+10);
	return;
}

// https://danceswithcode.net/engineeringnotes/interpolation/interpolation.html
static double a1(double ya, double yc) { return 0.5*(yc - ya); }
static double a2(double ya, double yb, double yc, double yd)
{ return 0.5*(2*ya - 5*yb + 4*yc - yd); }
static double a3(double ya, double yb, double yc, double yd)
{ return 0.5*(yd - ya + 3*yb - 3*yc); }

// 0 <= d <= 1; usually column - Cfloor
static double CatmullRom(double ya, double yb, double yc, double yd, double d)
{ return yb + a1(ya, yc)*d + a2(ya, yb, yc, yd)*d*d + a3(ya, yb, yc, yd)*d*d*d; }
static double CatmullRom(double *y, double d)
{ return y[1] + a1(y[0], y[2])*d + a2(y[0], y[1], y[2], y[3])*d*d + a3(y[0], y[1], y[2], y[3])*d*d*d; }

// interpolate a pixel at floating point row and column in plane
// !!! do NOT play with plane.data() elsewhere

#if 1

unsigned char get_CR(matrix<uchar> &plane, double row, double column)
{
	// handle borders
	int ncol = plane.ncol(), Rmax = plane.nrow() - 1, Cmax = ncol - 1;	// last addressable
	int Rfloor = (int)floor(row), Cfloor = (int)floor(column);
	// final interpolation is usually vertical
	double d = row - Rfloor, x = column - Cfloor, y[4] = {0};
	int y0 = 0, y1 = 0, y2 = 0, y3 = 0, l0 = 0, l1 = 0, l2 = 0, l3 = 0, l[4] = {0}, i = 0;

	if (row <= 0 || row >= Rmax) {	// no vertical interpolation
		if (column <= 0)			// no horizontal interpolation
			return plane(row <= 0 ? 0 : Rmax * ncol);

		if (column >= Cmax)
			return plane(row <= 0 ? Cmax : Cmax + Rmax * ncol);

		l0 = Cfloor;				// one horizontal interpolation
		if (0 == Cfloor)
		{ l1 = l0; l2 = l1 + 1; l3 = 1 + l2; }
		else { l1 = 1 + l0; l2 = 1 + l1;
				l3 = (Cfloor == Cmax - 1) ? l2 : 1 + l2; }
		x = CatmullRom(plane(l0), plane(l1), plane(l2), plane(l3), x);
		return (0 >= x) ? 0 : (255 <= x) ? 255 : (unsigned char)(0.5 + x); 
	}

	// border columns;  border rows already handled
	if (column <= 0 || column >= Cmax) {	// one vertical interpolation
		y0 = (0 == Rfloor) ? 0 : ncol;	// possibly duplicate a row
		l0 = Rfloor - y0 + (column >= Cmax) ? Cmax - 2 : 0;
		l1 = (0 == Rfloor) ? l0 : ncol + l0;
		l2 = l1 + ncol; l3 = (Rmax == Rfloor - 1) ? l2 : l2 + ncol;
		x = CatmullRom(plane(l0), plane(l1), plane(l2), plane(l3), d);
		return (0 >= x) ? 0 : (255 <= x) ? 255 : (unsigned char)(0.5 + x); 
	}

	// vertical interpolation, duplicating one column
	if (0 == Cfloor || Cfloor == Cmax - 1) {
		y0 = ncol * (Rfloor - 1) - (0 == Cfloor) ? 0 : 1;
		y1 = y0 + ncol; y2 = y1 + ncol; y3 = y2 + ncol;
		l0 = (0 == Cfloor) ? 1 : 0;
		for (i = 0; i < 3; i++)
			y[l0 + i] = CatmullRom(l0 + i, l1 + i, l2 + i, l3 + i, d);
		if (0 == Cfloor)
			y[0] = y[1];
		else y[3] = y[2];
		d = x;				// horizontal last interpolation
	}
	else { // horizontal interpolation, perhaps duplicating a row
		y0 = (0 == Rfloor) ? 0 : ncol;
		y1 = (0 == Cfloor) ? 1 : 0;		// duplicate a column?
		y2 = (0 == Cfloor || Cfloor == Cmax - 1) ? 3 : 4;
		l0 = Rfloor - y0 - y1;
		l1 = l0 + 1 - y1; l2 = 1 + l1; l3 = l2 + y1;
		for (i = 0; i < y2; i++) {
			y[y1 + i] = CatmullRom(l0, l1, l2, l3, x);
			l0 += ncol; l1 += ncol; l2 += ncol; l3 += ncol;
		}
		if (0 == Rfloor)
			y[0] = y[1];
		else if (Rfloor == Rmax - 1)
			y[3] = y[2];
	}
	x = CatmullRom(y[0], y[1], y[2], y[3], d);
	return (0 >= x) ? 0 : (255 <= x) ? 255 : (unsigned char)(0.5 + x); 
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

int read_coef(matrix<double> &coef, const char *fname)
{
	int rc = -1;
	if (4 != coef.nrow()) {
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

// called in gnuplot2file() after plane populates coef
void test_matrix_shift(matrix<double> &shifted, matrix<double> &x, matrix<double> &coef)
{
    int i, r, Grows = x.nrow();
	printf("test_matrix_shift()\n");
	double sy, sx;
	for (i = r = 0; r < Grows; r++) {
		matrix_shift(sy, sx, coef, 0, x(r, 2), x(r, 1));
		shifted.set(sx, i++); shifted.set(sy, i++);
		matrix_shift(sy, sx, coef, 2, x(r, 2), x(r, 1));
		shifted.set(sx, i++); shifted.set(sy, i++);
	}
}
