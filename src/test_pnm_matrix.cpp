#include "regress.h"
#include "matrix_pnm.h"

int test_pnm_matrix(const char **argv, matrix<double> &coef)
{
	const char *fname = argv[1];
	matrix<uchar> red {}, green {}, blue {};
//	int rc = read_pgm_matrix(red, fname);
//	if (0 <= rc)
//		return write_pgm_matrix("R:/Temp/write_matrix.pgm", red);
	int rc = read_matrix(red, green, blue, fname);
	if (0 <= rc)
	{
		write_pgm_matrix("R:/Temp/red_matrix.pgm", red);
		write_pgm_matrix("R:/Temp/green_matrix.pgm", green);
		write_pgm_matrix("R:/Temp/blue_matrix.pgm", blue);
		write_Bayer_matrix("R:/Temp/write_Bayer_matrix.pgm", red, green, blue);
		// 6 disables read_coef() in matrix_correction(),
		// then depends directly on matrix<> coef from polyEstimation()
		return matrix_correction(6, argv, coef, true);
	}
	return rc;
}
