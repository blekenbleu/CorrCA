typedef unsigned char uchar;

int test_read_matrix(const char *fname)
{
	matrix<uchar> red {}, green {}, blue {};
	int rc = read_pgm_matrix(red, fname);
	if (0 <= rc)
		return write_pgm_matrix("R:/Temp/write_matrix.pgm", red);
//	if (0 < read_matrix(red, green, blue, fname))
//		write_matrix(red, green, blue, "R:/Temp/write_matrix.ppm");
	return rc;
}

