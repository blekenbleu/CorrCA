#ifdef ARRAY_H /* Only included from array.h
 ; class template member functions' declarations and definitions
 ; must all be in the same header file, in this case array.h.
 ; As of 27 May 2026, array.h includes matrix_pnm.h, vector.cpp,
 ; matrix.cpp (which includes matrix_pnm.cpp)
 */

// https://en.wikipedia.org/wiki/Netpbm#Description

// seemingly cannot pass ifstream to read_pnm_header() unless it is also template <class T>, requiring a T argument
template <class T>
int read_pnm_header(int &bin, char &type, unsigned int &rows, unsigned int &columns, unsigned int &max, T *buf, std::ifstream &f)
{
	int len = 0;
	bin = true;
	char *first, *str_end;
	char whitespace[] = " \t";
	char b[200]{};

 
	f.getline(b, 19, '\n');
	len = (int)f.gcount();
	b[19] = '\0';
	if ('P' != b[0]) {
		error("not a PNM file!");
		return -2;
	}

	if ((type = b[1]) == '2' || '3' == type)
		bin = false;
  	else if(type != '5' && '6' != type) {
		error("not a supported PPM or PGM file!");
		return -3;
	}

	do {
		if ( !f.good() )
			return -4;
        f.getline(b, 19, '\n');
        len = (int)f.gcount();
        b[19] = '\0';
		first = b + strspn(b, whitespace);
	} while ('#' == *first);

	columns = strtol(first, &str_end, 10);
	if (str_end) {
		first = str_end;
		rows = strtol(first, &str_end, 10);
	} else return -5;

	do {
		if (!f.good() || '\0' != *str_end)
			return -6;
        f.getline(b, 19, '\n');
        len = (int)f.gcount();
        b[19] = '\0';
		first = b + strspn(b, whitespace);
	} while ('#' == *first);

	max = strtol(first, &str_end, 10);
	if (f.good() && '\0' == *str_end)
		return len;

	return -7;
}

/* separate img_bayer interleaved Bayer matrix into red, blue, and green pixel planes:
   rgrgr <- row 0;  unread pixels in lower case
   gbgbg
   rgRGR
   gbGBG
   rgRGR
   ... where each input pixel has only one color component.
   Note that first and last columns and rows are ignored...???
   .. with the first red value  from row 2, column 2,
      the first blue pixel from row 3, column 3
      and first green pixels from r2, column 3 and row 3, column 2.
   .. then red is written to imgR row 1 column 1, leaving a black pixel border.
   Green pixel plane has a 2-pixel wide black border.
 */
template <class T>
void deBayer_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB,
	int &bin, char &type, unsigned int &rows, unsigned int &columns, unsigned int &max, std::ifstream &f)
{
	printf("de-Bayer into separate red, green, blue planes... ");
	unsigned int c, c2 = 2 * columns;
	T *cR = imgR.data(), *g2 = imgG.data(), *cB = imgB.data(),
	  *g0 = g2, *g1 = 1 + columns + g0; 
	
	*cR++ = f.get();
	*g2 = f.get(); g2 += 2;
	for (c = 2; c < columns; c++)
	{
		*cR++ = f.get();
		*g2 = f.get();
		*(g2 - 1) = (T)(*(g2 - 2) + *g2 + 1) >> 1;	// first row 1D missing green averaging
		g2 += 2;
	}

	g2++;
	for (c = 0; c < columns; c++)
	{
		*g2 = f.get();
		g2 += 2;			// *g1 will fill in
		*cB++ = f.get();
	}
	
	for (unsigned int r = 2; r < rows; r++)
	{
		if (1 & r)
		{
			g1 = g2 - columns;
			g2--;
			g0 = g2 - c2;
			*g2 = f.get();
			*g1 = (T)(1 + *g2 + *g0);	// vertically interpolate first green column 
			g2 += 2;
			*cB++ = f.get();
			for (c = 2; c < columns; c++)
			{
				*g2 = f.get();
				*g1 = (T)(2 + *g0 + *(g0 + 2) + *g2 + *(g2 - 1))>>2;
				g0 += 2;
				g1 += 2;
				g2 += 2;
				*cB = f.get();
			}
		}
		else {
			g1 = g2 - columns;
			g2++;
			g0 = g2 - c2;
			*cR++ = f.get();
			*g2 = f.get(); g2 += 2;
			
			for (g2 = g1, c = 0; c < columns; c++)
			{
				*cR++ = f.get();
				*g2 = f.get();
				*g1 = (T)(2 + *g0 + *(g0 + 2) + *g2 + *(g2 - 2))>>2;
				g0 += 2;
				g1 += 2;
				g2 += 2;
			}
		}
	}
/*
			blue = img_bayer->data + i * 2 + 1 + (j * 2 + 1) * img_bayer->xsize;
			imgB->data[i + j * columns] = *blue;

			green = img_bayer->data + i*2+1+j*2*img_bayer->xsize;
			imgG->data[i*2+1+j*2*imgG->xsize] = *green;
			green = img_bayer->data + i * 2 + (j * 2 + 1) * img_bayer->xsize;
			imgG->data[i * 2 + (j * 2 + 1) * imgG->xsize] = *green;

			imgG->data[i * 2 + j * 2 * imgG->xsize] = 
			(
			 2 +  img_bayer->data[i * 2 + 1 + j * 2 * img_bayer->xsize]
			 + img_bayer->data[i * 2 - 1 + j * 2 * img_bayer->xsize]
			 + img_bayer->data[i * 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + (j * 2 - 1) * img_bayer->xsize]
			) >> 2;
			imgG->data[i * 2 + 1 + (j * 2 + 1) * imgG->xsize] = 
			(
			 2 + img_bayer->data[i * 2 + 1 + j * 2 * img_bayer->xsize]
		 	 + img_bayer->data[i * 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + 1 + (j * 2 + 2) * img_bayer->xsize]
			) >> 2;
 */
	printf("done\n");
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB, const char *fname)
{
	int bin = true, rc = 0;
	unsigned int len = 0, rows = 0, columns = 0, max = 0, r;
	if (sizeof(T) != sizeof(unsigned char))
	{
		printf("read_matrix() supports only unsigned char\n");
		return -1;
	}
	T buf[200]{};
	unsigned char *b = (unsigned char *)buf, * first = b, * str_end = b;
	char type = '0';
	double stuff = 0;
	char whitespace[] = " \t";

	std::ifstream f(fname, std::ios::binary);
	rc = read_pnm_header(bin, type, rows, columns, max, buf, f);
	if (0 > rc)
	{
		printf("read_pnm_header(%s):  failed\n", fname);
		return rc;
	}
	
	imgR.init(columns, rows);
	imgB.init(columns, rows);
	T *cR = imgR.data(), * cG = imgG.data(), * cB = imgB.data(), * cx = cR;

	if ('6' == type)
	{
		imgG.init(columns, rows);
		for (r = 0; r < rows; r++)
			for (cx = cR + columns; cR < cx && f.good(); cR++) {
				*cR = f.get(); *cG++ = f.get(); *cB++ = f.get();
			}
	}
	else {
		imgG.init(columns*2, 2*rows);
		deBayer_matrix(imgR, imgG, imgB, bin, type, rows, columns, max, f);
		return len;
	}

	if (rows > r || cx > cR)
	{
		printf("read_matrix(%s): stopped at row %d/%d, column %d/%d\n",
				fname, r, rows, columns - (int)(cx - cR), columns);
		return -9;
	}

	return len;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_matrix(matrix<T> &img, const char *fname)
{
	int rc = 0, bin = true;
	unsigned int rows = 0, columns = 0, max = 0, len = 0;
	T buf[200]{};
	char type = '0';
	char *b = (char *)buf, *first = b, *str_end = b;

	std::ifstream f(fname, std::ios::binary);
	rc = read_pnm_header(bin, type, rows, columns, max, buf, f);
    if (0 > rc)
        return rc;

	if ('2' != type && '5' != type)
	{
		error("not a supported PGM file!");
        return -2;
	}

	unsigned int r, c = columns, d;
	img.init(rows, columns);
	char *dest = (char* )img.data(), *end = dest;

	if (bin && 255 == max)
	{
		int g = true;
		for (r = 0; r < rows && g; r++)
			for (c = columns; c > 0 && g; )
			{
				f.read(dest, c);
				d = (int)f.gcount();
				dest += d;
				c -= d;
				if (!f.good())
					g = false;
			}
		f.close();
		if (c != 0)
			printf("read_matrix(%s):  %d/%d columns unread at row %d\n",
					fname, c, columns, r);
		return len;
		
	} else {
		T *dT = img.data(), *xT;

		for (r = 0; r < rows && f.good(); r++)
			for (xT = dT + columns; dT < xT && f.good(); dT++)
			{
				f.getline(b, 19, ' ');
       			b[19] = '\0';
				*dT = (T)strtod(b, &str_end);
			}

		if (r < rows || dest < end)
			printf("read_matrix(%s): stopped at row %d/%d, column %d/%d\n",
				fname, r, rows, columns - (int)(end - dest), columns);
	}
	return len;
}

#endif // ARRAY_H
