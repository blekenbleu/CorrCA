int read_pgm_matrix(matrix<T> &img, const char *filename);

int read_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB,
				 const char *filename);

int write_pgm_matrix(const char *filename, matrix<T> &img);

int write_Bayer_matrix(const char *fname, matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB);

int write_matrix(const char *fname, matrix<double> &imgR, matrix<unsigned char> &imgG,
				  matrix<double> &imgB);
