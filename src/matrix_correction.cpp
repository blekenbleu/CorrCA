#include "regress.h"
#include "matrix_pnm.h"

/* Gcols can be 2x imgF.nchar();
 ; since polynomial coefficients are based on image dimensions rescaled tp [0:1],
 ; interpolation for pixel shift also works for resampling...
 */
static void correct_channel(matrix<uchar> &imgIn, matrix<uchar> &imgOut,
					matrix<double> &coef, int color)	// red color = 0; blue color = 2
{
	int Grows = imgIn.nrow();
	double dy = 0, scaleY = 1.0 / Grows, dx = 0, xs = 1.0, scaleX = 1.0 / imgIn.ncol();
	int wIn = imgIn.ncol(), hIn = imgIn.nrow();
	double hs = scaleY * hIn;
	printf("channel correction... ");
	for (int i = 0, y = 0, yc = 0; y < Grows; y++) {
		double ys = y * scaleY, ry = ys * hIn;
		for (double xs = 0, rx = 0, x = 0; x < 1; x += scaleX, i++) {
			// polynomials are based on pixel locations ys, xs scaled [0:1]
			matrix_shift(dy, dx, coef, color, ys, xs);
			// +0.5 to compensate -0.5 in interpolation function?
			// dy, dx are direction shifts in imgIn pixel dimensions
			imgOut(i) = get_CR(imgIn, dy + ry, dx + rx);
			rx += hs;
		}
		double percent = (100.0 * y) / Grows;
		if (!(y % (int)(0.2 * Grows))) printf("%i%c", (int)percent+1, '%');
		else if (!(y % (int)(0.04 * Grows))) printf(".");
	}
	printf("done.\n");
}

int matrix_correction(int argc, const char **argv, matrix<double> &coef, bool clr)
{
	int rc = -1;
	printf("\nCA reduction... \n");
	const char* fnameRGB = argv[1];

	if (6 < argc && 0 != (rc = read_coef(coef, argv[6])))
        return rc;

	if(4 != coef.nrow()) {
		printf("\t>>>> invalid matrix<> coef <<<<\n");
		return rc;
	}
	/* matrix<> data are 1D;
	 ; matrix<> plane(row, column) == plane(column + row * plane.ncol())
	 ; consequently, process an entire matrix<> old to matrix<> new by:
	 ; for (int i = 0, x = plane.nrow()*plane.ncol(); i < ix; i++)
	 ;	  new(i) = process(i, old(i));
	 */
	matrix<uchar> Rin{}, Gin{}, Bin{};
	read_matrix(Rin, Gin, Bin, fnameRGB);
	uint Gcols = Gin.ncol(), Grows = Gin.nrow();

	matrix<uchar> Rout(Grows, Gcols);
	matrix<uchar> Bout(Grows, Gcols);

	printf("Red  "); correct_channel(Rin, Rout, coef, 0);
	printf("Blue "); correct_channel(Bin, Bout, coef, 2);

	printf("\nSaving images to file... \n");
	if (9 < argc)
	{
		write_pgm_matrix(argv[7], Rout);
		write_pgm_matrix(argv[8], Gin);
		write_pgm_matrix(argv[9], Bout);
	}
	else write_matrix(argv[5], Rout, Gin, Bout);
	return 0;
}
