int read_pgm_matrix(matrix<T> &img, const char *filename);

int read_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB,
				 const char *filename);

int write_pgm_matrix(const char *filename, matrix<T> &img);

void write_matrix(const matrix<T> &imgR, const matrix<T> &imgG,
				  const matrix<T> &imgB, const char *filename);
