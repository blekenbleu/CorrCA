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
void deBayer_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB, std::ifstream &f)
{
	int columns = imgG.ncol(), rows = imgG.nrow(), c2 = 2 * columns;
	printf("de-Bayer into separate red, green, blue planes... ");
	int cR = 0, g2 = 0, cB = 0, g0 = g2, g1 = 1 + columns + g0, row = 0,
		xR = cR, xB = cB, x = g2 + columns, col = 0, stop = rows * columns;
	
	imgR.set(f.get(), cR++);
	imgG.set(f.get(), g2); g2 += 2;
	while (g2 < x)
	{
		imgR.set(f.get(), cR++);
		imgG.set(f.get(), g2);
		imgG.set((T)(1 + imgG(g2 - 2) + imgG(g2)) >> 1, g2 - 1);	// first row 1D missing green averaging
		g2 += 2;
	}

	row++;
	col += columns;
	x = g2 + columns;
	g2++;
	while (g2 < x)
	{
		imgG.set(f.get(), g2);
		g2 += 2;			// *g1 will fill in
		imgB.set(f.get(), cB++);
	}

	row++;
	col += columns;
	while (g2 < stop)
	{
		if (1 & row)
		{
			g1 = g2 - columns;
			x = g2 + columns;
			g2++;
			imgG.set(f.get(), g2);
			imgG.set((T)(1 + imgG(g2) + imgG(g0))/2, g1);	// vertically interpolate first green column 
			g2 += 2;
			g0 = g2 - c2;
			imgB.set(f.get(), cB++);
			while (g2 < x)
			{
				imgG.set(f.get(), g2);
				imgG.set((T)(2 + imgG(g0) + imgG(g0 - 2) + imgG(g2) + imgG(g2 - 1))>>2, g1);
				g0 += 2;
				g1 += 2;
				g2 += 2;
				imgB.set(f.get(), cB++);
			}
		}
		else {
			imgR.set(f.get(), cR++);
			g1 = g2 - columns;
			g2--;
			x = g2 + columns;
			imgG.set(f.get(), g2); g2 += 2;
			g0 = g2 - c2;
			
			while (g2 < x)
			{
				imgR.set(f.get(), cR++);
				imgG.set(f.get(), g2);
				imgG.set((T)(2 + imgG(g0) + imgG(g0 - 2) + imgG(g2) + imgG(g2 - 2))>>2, g1);
				g0 += 2;
				g1 += 2;
				g2 += 2;
			}
		}
		col += columns;
		row++;
	}
/*
			blue = i * 2 + 1 + (j * 2 + 1) * img_bayer->xsize;
			imgB.set(img_bayer(blue), i + j * columns);

			green = i*2+1+j*2*img_bayer->xsize;
			imgG.set(img_bayer(green), i*2+1+j*2*imgG->xsize);
			green = i * 2 + (j * 2 + 1) * img_bayer->xsize;
			imgG.set(img_bayer(green), i * 2 + (j * 2 + 1) * imgG->xsize);

			imgG.set(
			 2 + img_bayer(i * 2 + 1 + j * 2 * img_bayer->xsize)
			 + img_bayer(i * 2 - 1 + j * 2 * img_bayer->xsize)
			 + img_bayer(i * 2 + (j * 2 + 1) * img_bayer->xsize)
			 + img_bayer(i * 2 + (j * 2 - 1) * img_bayer->xsize)
			) >> 2, i * 2 + j * 2 * imgG->xsize);
			imgG.set(
			 2 + img_bayer(i * 2 + 1 + j * 2 * img_bayer->xsize)
		 	 + img_bayer(i * 2 + (j * 2 + 1) * img_bayer->xsize)
			 + img_bayer(i * 2 + 2 + (j * 2 + 1) * img_bayer->xsize)
			 + img_bayer(i * 2 + 1 + (j * 2 + 2) * img_bayer->xsize)
			) >> 2, i * 2 + 1 + (j * 2 + 1) * imgG->xsize);
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
	
	imgG.init(rows, columns);
	int cG = 0;

	if ('6' == type)
	{
		imgR.init(rows, columns);
		imgB.init(rows, columns);
		int cR = 0, cB = 0, cx = cR;
		for (r = 0; r < rows; r++)
			for (cx = cR + columns; cR < cx && f.good(); cR++) {
				imgR.set(f.get(), cR); imgG.set(f.get(), cG++); imgB.set(f.get(), cB++);
			}

		if (rows > r || cx > cR)
		{
			printf("read_matrix(%s): stopped at row %d/%d, column %d/%d\n",
				fname, r, rows, columns - (int)(cx - cR), columns);
			return -9;
		}
	}
	else {
		imgR.init(r = (1 + rows)/2, columns/2);
		imgB.init(rows - r, columns/2);		// perhaps 1 more red than blue row
		deBayer_matrix(imgR, imgG, imgB, f);
		return len;
	}


	return len;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_pgm_matrix(matrix<T> &img, const char *fname)
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
	int dest = 0, end = dest;

	if (bin && 255 == max)
	{
		int g = true;
		for (r = 0; r < rows && g; r++)
			for (c = columns; c > 0 && g; )
			{
				f.read(img.data(dest), c);
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
		int dT = 0, xT;

		for (r = 0; r < rows && f.good(); r++)
			for (xT = dT + columns; dT < xT && f.good(); dT++)
			{
				f.getline(b, 19, ' ');
       			b[19] = '\0';
				img.set((T)strtod(b, &str_end), dt);
			}

		if (r < rows || dest < end)
			printf("read_matrix(%s): stopped at row %d/%d, column %d/%d\n",
				fname, r, rows, columns - (int)(end - dest), columns);
	}
	return len;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int write_pgm_matrix(const char *fname, matrix<T> &img)
{
	if (sizeof(T) != sizeof(unsigned char))
	{
		printf("write_matrix() supports only bytes\n");
		return -1;
	}
	char header[35]; sprintf(header, "P5\n%d %d\n# write_matrix\n255\n",
		img.ncol(), img.nrow());
	std::ofstream f(fname, std::ios::binary);
	if (f)
	{
		int len = (int)strlen(header), buf = 0;
		f.write(header, len);
		len = img.ncol();
		for (int rx = buf + img.nrow() * len; buf < rx; buf += len)
			f.write((char *)img.data(buf), len);
		f.close();
		return len;
	} else return -1;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int write_Bayer_matrix(const char *fname, matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB)
{
	int w = imgG.ncol(), rowsB = imgB.nrow(), colsR = imgR.ncol(), rowsG = imgG.nrow();
	if (sizeof(T) != sizeof(unsigned char))
	{
		printf("write_pgm_matrix() supports only bytes\n");
		return -1;
	}
	if (2 * colsR != w) {
		printf("write_pgm_matrix() expects green 2x red or blue\n");
		return -2;
	}
	int b = 0, g = 0, r = 0, g2 = 1 + g;
	char header[35]; sprintf(header, "P5\n%d %d\n# write_matrix\n255\n",
		imgG.ncol(), imgG.nrow());
	int len = (int)strlen(header);
	std::ofstream f(fname, std::ios::binary);
	if (f)
	{
	  int gx;

	  f.write(header, len);
	  for (int c = 0; c < rowsB; c++)
	  {
		for (gx = g + w; g < gx; g += 2) {
			f.put(imgR(r++)); f.put(imgG(g)); }
		g++;
		for (gx = g + w; g < gx; g += 2) {
			f.put(imgG(g)); f.put(imgB(b++)); }
		g--;
	  }
	  if(imgR.nrow() > rowsB)
		for (gx = g + w; g < gx; g += 2) {
            f.put(imgR(r++)); f.put(imgG(g)); }
	  f.close();
	} else return -1;
	return w;
}
#endif // ARRAY_H
