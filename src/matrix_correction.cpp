/*
#include "array.h"
*/
typedef unsigned int uint;
typedef unsigned char uchar;

template <typename T>
void read_poly(const char *fname, matrix<T> &coef)
{
	coef.init(4, 10, 1);
}

/* Gcols can be 2x imgF.nchar();
 ; since polynomial coefficients are based on image dimensions rescaled tp [0:1],
 ; interpolation for pixel shift also works for resampling...
 */
static void correct_channel(matrix<uchar> &imgIn, matrix<uchar> &imgOut,
					matrix<double> params, int color,
					int Gcols, int Grows)
{
	double dy = 0, scaleY = 1.0 / Grows, dx = 0, xs = 1.0, scaleX = 1.0 / Gcols;
	printf("calculating channel correction... ");
	for (int i = 0, y = 0, yc = 0; y < Grows; y++) {
		double ys = y * scaleY, ry = ys * imgIn.nrow();
		for (int x = 0; x < Gcols; x++, i++) {
			// polynomials are for normalized pixel locations
			matrix_shift(dy, dx, params, color, ys, xs = x * scaleX);
			// +0.5 to compensate -0.5 in interpolation function?
			// x, y rescaled to imgIn dimensions, then add dy, dx
			imgOut(i) = get_CR(imgIn, dy + ry, dx + xs * imgIn.ncol());
		}
		double percent = (100.0 * y) / Grows;
		if (!(y % (int)(0.2*Grows))) printf("%i%c", (int)percent+1, '%');
		else if (!(y % (int)(0.04*Grows))) printf(".");
	}
	printf("done.\n");
}

template <typename T>
void matrix_correction(int argc, const char ** argv, bool clr)
{
	printf("\nAberration correction... \n");
	const char* fnameRGB = argv[1];
	const char* fnamePoly = argv[2];
	const char *fnameR = argv[4], *fnameG, *fnameB;
	if (7 == argc)
	{
		fnameG = argv[5];
		fnameB = argv[6];
	}
	int degX = 11, degY = 11;
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	matrix<uchar> Rin{}, Gin{}, Bin{};

	/* matrix<> data are 1D;
	 ; matrix<> plane(row, column) == plane(column + row * plane.ncol())
	 ; consequently, process an entire matrix<> old to matrix<> new by:
	 ; for (int i = 0, x = plane.nrow()*plane.ncol(); i < ix; i++)
	 ;	  new(i) = process(i, old(i));
	 */
	read_matrix(Rin, Gin, Bin, fnameRGB);
	uint Gcols = Gin.ncol(), Grows = Gin.nrow();

	matrix<T> params{}; read_poly<T>(fnamePoly, params);

//	printf("Gcols = %d;  Grows = %d, Gin.ncol() = %d, Gin.nrow() = %d for %s\n",
//			Gcols, Grows, Gin.ncol(), Gin.nrow(), fnameRGB);
//	Gcols = 5634;  Grows = 3752, Gin.ncol() = 5634, Gin.nrow() = 3752 for ../../../../data/_MG_7626.pgm

	int spline_order = 3;
	matrix<uchar> Rout; Rout.init(Grows, Gcols);
	matrix<uchar> Bout; Bout.init(Grows, Gcols);
	T xp = 0.2 + 0.5 * Gin.ncol(), yp = 0.2 + 0.5 * Gin.nrow();

	printf("Red  ");
	correct_channel(Rin, Rout, params, 0, Gcols, Grows);
	printf("Blue ");
	correct_channel(Bin, Bout, params, 2, Gcols, Grows);

	printf("\nSaving images to file... \n");
	if (7 == argc)
	{
		write_pgm_matrix(argv[4], Rout);
		write_pgm_matrix(argv[5], Gin);
		write_pgm_matrix(argv[6], Bout);
	}
	else write_matrix(argv[7], Rout, Gin, Bout);
}
