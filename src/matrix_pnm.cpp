#ifdef ARRAY_H // Do nothing if not included from array.h

namespace libNumerics {
#define _CRT_SECURE_NO_DEPRECATE    // fopen(), fscanf() warnings; must be macro;  unavoidable VCR101

template <typename T> int matrix<T>::read_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB,
                 char *fname)
{
	std::ifstream f(fname);
	if ( f.good() ) {
		char buf[200];
		f.getline(buf, 19, '\n');
		int len = f.gcount();
		buf[19] = '\0'
		if ('P' != buf[0]) {
			error("not a PNM file!");
			return -1;
		}

		if ((type = buf[1]) == '2' || '3' == type)
			bin = false;
	  	else if(type == '5' || '6' == type)
    		bin = true;
	 	else {
			error("not a supported PNM file!");
			return -2;
		}
		printf("read_matrix() %s P%c len %d strlen %d\n", fname, type, len, strlen(buf));
	}
}
} // namespace libNumerics

#endif
