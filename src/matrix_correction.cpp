/*
#include "array.h"
*/
typedef unsigned int uint;
typedef unsigned char uchar;

// Gcols can be 2x imgF.nchar()
void correct_channel(matrix<uchar> &imgF, matrix<uchar> &imgFz,
					vector<double> paramsXF, vector<double> paramsYF,
					int spline_order, int degX, int degY, double xp, double yp,
					int Gcols, int Grows, double scale)
{
	printf("calculating channel correction... ");
	for (int i = 0; i < Gcols; i++) {
		for (int j = 0; j < Grows; j++) {
			double p1=0, p2=0;
			// apply polynomials to green pixel locations
			undistortPixel(p1, p2, paramsXF, paramsYF, i, j, xp, yp, degX, degY);
			// +0.5 to compensate -0.5 in interpolation function
			imgFz(i+j*imgFz.ncol()) = get_CR(imgF, spline_order, p1 / scale + 0.5, p2 / scale + 0.5);
		}
		double percent = ((double)i / (double)Gcols)*100;
		if (!(i % (int)(0.2*Gcols))) printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*Gcols))) printf(".");
	}
	printf("done.\n");
}

template <typename T>
void matrix_correction(int argc, const char ** argv, bool clr)
{
	printf("\nAberration correction... \n");
	const char* fnameRGB = argv[1];
	const char* fnamePolyR = argv[2];
	const char* fnamePolyB = argv[3];
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

	read_matrix(Rin, Gin, Bin, fnameRGB);
	uint Gcols = Gin.ncol(), Grows = Gin.nrow();

	vector<T> paramsR = read_poly<T>(fnamePolyR, degX, degY);
	vector<T> paramsB = read_poly<T>(fnamePolyB, degX, degY);
	vector<T> paramsXR = paramsR.copyRef(0, sizex-1);
	vector<T> paramsYR = paramsR.copyRef(sizex, sizex+sizey-1);
	vector<T> paramsXB = paramsB.copyRef(0, sizex-1);
	vector<T> paramsYB = paramsB.copyRef(sizex, sizex+sizey-1);

//	printf("Gcols = %d;  Grows = %d, Gin.ncol() = %d, Gin.nrow() = %d for %s\n",
//			Gcols, Grows, Gin.ncol(), Gin.nrow(), fnameRGB);
//	Gcols = 5634;  Grows = 3752, Gin.ncol() = 5634, Gin.nrow() = 3752 for ../../../../data/_MG_7626.pgm

	int spline_order = 3;
	matrix<uchar> Rout; Rout.init(Grows, Gcols);
	matrix<uchar> Bout; Bout.init(Grows, Gcols);
	T xp = 0.2 + 0.5 * Gin.ncol(), yp = 0.2 + 0.5 * Gin.nrow();

	printf("Red  ");
	correct_channel(Rin, Rout, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, Gcols, Grows, 2);
	printf("Blue ");
//	correct_channel<T>(Bin, Bout, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, Gcols, Grows, 2);

	printf("\nSaving images to file... \n");
	if (7 == argc)
	{
		write_pgm_matrix(argv[4], Rout);
		write_pgm_matrix(argv[5], Gin);
		write_pgm_matrix(argv[6], Bout);
	}
	else write_matrix(argv[7], Rout, Gin, Bout);
}
