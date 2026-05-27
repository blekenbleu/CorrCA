typedef unsigned char uchar;

int test_read_matrix(const char *fname)
{
	matrix<uchar> red {}, green {}, blue {};
	return read_matrix(red, green, blue, fname);
}

