#include "regress.h"
#include "matrix_pnm.h"

/* Gcols can be 2x imgF.nchar();
 ; Polynomial coefficients are for image dimensions rescaled to [0:1];
 ; interpolation works for both pixel shift and resampling...
 */
static void poly_color(matrix<uchar> &imgIn, matrix<uchar> &imgOut,
					matrix<double> &coef, int color)	// red color = 0; blue color = 2
{
	int wIn = imgIn.ncol(), hIn = imgIn.nrow(), io = 0, cy = 0;
	double scaleY = 1.0 / hIn, scaleX = 1.0 / wIn;
	printf("poly_color(%d)... ", color);
	for (double sy = 0; sy < 1.0; sy += scaleY) {
		int rx = 0;
		for (double dx, dy, sx = 0; sx < 1; sx += scaleX) {
			// polynomials are based on pixel locations sy, sx scaled [0:1]
			matrix_shift(dy, dx, coef, color, sy, sx);
			// +0.5 to compensate -0.5 in interpolation function?
			// dy, dx are direction shifts in imgIn pixel dimensions
			imgOut(io++) = get_CR(imgIn, dy + cy, dx + rx++);
		}
		double percent = (100.0 * cy) / hIn;
		if (!(cy % (int)(0.2 * hIn))) printf("%i%c", (int)percent+1, '%');
		else if (!(cy % (int)(0.04 * hIn))) printf(".");
		cy++;
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
	uint Gcols = Gin.ncol(), hIn = Gin.nrow();

	matrix<uchar> Rout(hIn, Gcols);
	matrix<uchar> Bout(hIn, Gcols);

	printf("Red  "); poly_color(Rin, Rout, coef, 0);
	printf("Blue "); poly_color(Bin, Bout, coef, 2);

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
