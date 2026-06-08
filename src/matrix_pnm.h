template <class T>
int read_pgm_matrix(matrix<T> &img, const char *filename);

int read_matrix(matrix<uchar> &imgR, matrix<uchar> &imgG, matrix<uchar> &imgB,
				 const char *filename);

int write_pgm_matrix(const char *filename, matrix<uchar> &img);

template <class T>
int write_Bayer_matrix(const char *fname, matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB);

int write_matrix(const char *fname, matrix<uchar> &imgR, matrix<uchar> &imgG, matrix<uchar> &imgB);
