/*
#include "array.h"
*/
typedef unsigned int uint;
typedef unsigned char uchar;

template <typename T>
void matrix_correction(int argc, char ** argv, bool clr)
{
	printf("\nAberration correction... \n");
	char* fnameRGB = argv[1];
	char* fnamePolyR = argv[2];
	char* fnamePolyB = argv[3];
	char *fnameR = argv[4], *fnameG, *fnameB;
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
	T xp = 0.2, yp = 0.2;
	
		xp += Gin.ncol()/2;
		yp += Gin.nrow()/2;
 

	printf("Red  ");
//	correct_channel<T>(Rin, Rout, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, Gcols, Grows, 2);
	printf("Blue ");
//	correct_channel<T>(Bin, Bout, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, Gcols, Grows, 2);

	printf("\nSaving images to file... \n");
	if (7 == argc)
	{
		write_pgm_matrix(argv[4], Rout);
		write_pgm_matrix(argv[5], Gin);
		write_pgm_matrix(argv[6], Bout);
	}
	else write_matrix(argv[4], Rout, Gin, Bout);
}
