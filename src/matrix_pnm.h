template <class T>
int read_pgm_matrix(matrix<T> &img, const char *filename);

int read_matrix(matrix<uchar> &imgR, matrix<uchar> &imgG, matrix<uchar> &imgB, const char *filename);

int write_pgm_matrix(const char *filename, matrix<uchar> &img);

int write_Bayer_matrix(const char *fname, matrix<uchar> &imgR, matrix<uchar> &imgG, matrix<uchar> &imgB);

int write_matrix(const char *fname, matrix<uchar> &imgR, matrix<uchar> &imgG, matrix<uchar> &imgB);

int matrix_correction(int argc, const char **argv, matrix<double> &coef, bool clr);
