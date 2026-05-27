#ifdef ARRAY_H /* Only included from array.h
 ; class template member functions' declarations and definitions
 ; must all be in the same header file, in this case array.h.
 ; As of 27 May 2026, array.h includes matrix_pnm.h, vector.cpp,
 ; matrix.cpp (which includes matrix_pnm.cpp)
 */
#include <stdlib.h>

// https://en.wikipedia.org/wiki/Netpbm#Description

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB, const char *fname)
{
	char whitespace[] = " \t";
	int bin = true, len = 0, rows, columns, max;
	char buf[200], type, *first, *str_end;
	double stuff;

	std::ifstream f(fname);
	if ( f.good() ) {
		f.getline(buf, 19, '\n');
		len = (int)f.gcount();
		buf[19] = '\0';
		if ('P' != buf[0]) {
			error("not a PNM file!");
			return -1;
		}

		if ((type = buf[1]) == '2' || '3' == type)
			bin = false;
	  	else if(type != '5' && '6' != type) {
			error("not a supported PPM or PGM file!");
			return -2;
		}
	} else return -3;

	if ( f.good() )
	  do {
        f.getline(buf, 19, '\n');
        len = (int)f.gcount();
        buf[19] = '\0';
		first = buf + strspn(buf, whitespace);
	  } while ('#' == *first);
	columns = strtol(first, &str_end, 10);
	if (str_end) {
		first = str_end;
		rows = strtol(first, &str_end, 10);
	}
	if (f.good() && '\0' == *str_end)
	  do {
        f.getline(buf, 19, '\n');
        len = (int)f.gcount();
        buf[19] = '\0';
		first = buf + strspn(buf, whitespace);
	  } while ('#' == *first);
	max = strtol(first, &str_end, 10);
	if (f.good() && '\0' == *str_end) {
        f.getline(buf, 198, '\n');
        len = (int)f.gcount();
        buf[199] = '\0';
		if (2 < len)
			stuff = strtod(buf, &str_end);
	}
		

		// printf("read_matrix() %s P%c len %d strlen %d\n", fname, type, len, (int)strlen(buf));
		// read_matrix() ../../../../data/_MG_7626.pgm P5 len 3 strlen 2
		if (2 == strlen(buf))
		{
			f.getline(buf, 19, '\n');
		}

	return len;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_matrix(matrix<T> &imgR, const char *fname)
{
	char whitespace[] = " \t";
	int bin = true, len = 0, rows, columns, max;
	char buf[200], type, *first, *str_end;
	double stuff;

	std::ifstream f(fname);
	if ( f.good() ) {
		f.getline(buf, 19, '\n');
		len = (int)f.gcount();
		buf[19] = '\0';
		if ('P' != buf[0]) {
			error("not a PNM file!");
			return -1;
		}

		if ((type = buf[1]) == '2')
			bin = false;
	  	else if(type != '5') {
			error("not a supported PGM file!");
			return -2;
		}
	} else return -3;

	if ( f.good() )
	  do {
        f.getline(buf, 19, '\n');
        len = (int)f.gcount();
        buf[19] = '\0';
		first = buf + strspn(buf, whitespace);
	  } while ('#' == *first);
	columns = strtol(first, &str_end, 10);
	if (str_end) {
		first = str_end;
		rows = strtol(first, &str_end, 10);
	}
	if (f.good() && '\0' == *str_end)
	  do {
        f.getline(buf, 19, '\n');
        len = (int)f.gcount();
        buf[19] = '\0';
		first = buf + strspn(buf, whitespace);
	  } while ('#' == *first);
	max = strtol(first, &str_end, 10);
	if (f.good() && '\0' == *str_end) {
        f.getline(buf, 198, '\n');
        len = (int)f.gcount();
        buf[199] = '\0';
		if (2 < len)
			stuff = strtod(buf, &str_end);
	}
		

		// printf("read_matrix() %s P%c len %d strlen %d\n", fname, type, len, (int)strlen(buf));
		// read_matrix() ../../../../data/_MG_7626.pgm P5 len 3 strlen 2
		if (2 == strlen(buf))
		{
			f.getline(buf, 19, '\n');
		}

	return len;
}

#endif // ARRAY_H
