typedef unsigned char uchar;

int test_read_matrix(const char *fname)
{
	matrix<uchar> red {}, green {}, blue {};
	if (0 < read_matrix(red, fname))
		write_matrix(red, "R:/Temp/write_matrix.pgm");
//	if (0 < read_matrix(red, green, blue, fname))
//		write_matrix(red, green, blue, "R:/Temp/write_matrix.ppm");
}

