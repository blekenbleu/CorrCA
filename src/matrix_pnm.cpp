#ifdef ARRAY_H /* Only included from array.h
 ; class template member functions' declarations and definitions
 ; must all be in the same header file, in this case array.h.
 ; As of 27 May 2026, array.h includes matrix_pnm.h, vector.cpp,
 ; matrix.cpp (which includes matrix_pnm.cpp)
 */

// https://en.wikipedia.org/wiki/Netpbm#Description

// seemingly cannot pass ifstream to read_pnm_header()
// unless it is also template <class T>, requiring a T argument
template <class T>
int read_pnm_header(char &type, matrix<T> &img, std::ifstream &f)
{
	int len = 0, max = -7, rows = 0, columns = 0;
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

	if ((b[1] != '2' && '5' != b[1])
  		&& ((b[1] != '3' && '6' != b[1]) || '3' != type)) {
		error("not a supported PPM or PGM file!");
		return -3;
	}
	type = b[1];

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
	{
		img.init(rows, columns);
		return max;
	}

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
int deBayer_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB, std::ifstream &f)
{
	int columns = imgG.ncol(), rows = imgG.nrow(), c2 = 2 * columns;
	printf("de-Bayer into separate red, green, blue planes... ");
	int cR = 0, cB = 0, cG = 0, g1 = 0, row = 0, stop = cG + columns;
	
	imgR.set(f.get(), cR++);		// row 0
	cG++;
	imgG.set(f.get(), cG++); cG++;
	while (cG < stop)
	{
		imgR.set(f.get(), cR++);
		imgG.set(f.get(), cG);
		imgG.set((T)(1 + imgG(cG - 2) + imgG(cG)) >> 1, cG - 1);	// first row 1D missing green averaging
		cG += 2;
	}

	row++;							// row 1
	cG--;
	stop = cG + columns;
	while (cG < stop)
	{
		imgG.set(f.get(), cG);
		cG += 2;			// will fill in
		imgB.set(f.get(), cB++);
	}

	row++;							// row 2
	stop = rows * columns;
	while (cG < stop)
	{
		int x;

		if (1 & row)
		{
			cG--;
			x = cG + columns;
			imgG.set(f.get(), cG);
			imgG.set((T)(1 + imgG(cG) + imgG(cG - c2))/2, cG - columns);	// vertically interpolate first green column 
			cG += 2;
			imgB.set(f.get(), cB++);
			for (; cG < x; cG += 2)
			{
				imgG.set(f.get(), cG);
				g1 = cG - columns;
				imgG.set((2 + imgG(cG - c2) + imgG(g1 - 1) + imgG(cG) + imgG(g1 + 1))>>2, g1);
				imgB.set(f.get(), cB++);
			}
		}
		else {
			imgR.set(f.get(), cR++);
			x = cG + columns;
			cG++;				// vertically interpolated in next row
			imgG.set(f.get(), cG); cG += 2;
			for (; cG < x; cG += 2)
			{
				imgR.set(f.get(), cR++);
				imgG.set(f.get(), cG);
				g1 = cG - columns;
				imgG.set((2 + imgG(cG - c2) + imgG(g1 - 1) + imgG(cG) + imgG(g1 + 1))>>2, g1);
			}
		}
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
	return 0;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_matrix(matrix<T> &imgR, matrix<T> &imgG, matrix<T> &imgB, const char *fname)
{
	int max = 0, r;
	unsigned int len = 0;
	if (sizeof(T) != sizeof(unsigned char))
	{
		printf("read_matrix() supports only unsigned char\n");
		return -1;
	}
	T buf[200]{};
	unsigned char *b = (unsigned char *)buf, * first = b, * str_end = b;
	char type = '3';	// color plane count
	double stuff = 0;
	char whitespace[] = " \t";

	std::ifstream f(fname, std::ios::binary);
	max = read_pnm_header(type, imgG, f);
	if (0 > max)
	{
		printf("read_pnm_header(%s):  failed\n", fname);
		return max;
	}
	
	if ('6' == type)
	{
		imgR.init(imgG.nrow(), imgG.ncol());
		imgB.init(imgG.nrow(), imgG.ncol());
		int cR = 0, cG = 0, cB = 0, cx = cR;
		for (r = 0; r < imgG.nrow(); r++)
			for (cx = cR + imgG.ncol(); cR < cx && f.good(); cR++) {
				imgR.set(f.get(), cR); imgG.set(f.get(), cG++); imgB.set(f.get(), cB++);
			}

		if (imgG.nrow() > r || cx > cR)
		{
			printf("read_matrix(%s): stopped at row %d/%d, column %d/%d\n",
				fname, r, imgG.nrow(), imgG.ncol() - (int)(cx - cR), imgG.ncol());
			return -9;
		}

		return 0;
	}
	else if ('5' == type) {
		imgR.init(r = (1 + imgG.nrow())/2, imgG.ncol()/2);
		imgB.init(imgG.nrow() - r, imgG.ncol()/2);		// perhaps 1 more red than blue row
		return deBayer_matrix(imgR, imgG, imgB, f);
	}
	else printf("read_matrix(%s): type P%c not yet supported\n", fname, type);

	return -10;
}

template <class T>	// https://users.cis.fiu.edu/~weiss/Deltoid/vcstl/templates
int read_pgm_matrix(matrix<T> &img, const char *fname)
{
	int rc = 0, max = 0;
	char buf[200]{};
	char type = '1';

	std::ifstream f(fname, std::ios::binary);
	max = read_pnm_header(type, img, f);
	if (0 > max)
		return max;

	unsigned int r, c = img.ncol(), d;
	int dest = 0, end = dest;

	if ('5' == type && 255 == max)
	{
		int g = true;
		for (r = 0; r < img.nrow() && g; r++)
			for (c = img.ncol(); c > 0 && g; )
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
					fname, c, img.ncol(), r);
		return c;
		
	} else {
		int dT = 0, xT;

		for (r = 0; r < img.nrow() && f.good(); r++)
			for (xT = dT + img.ncol(); dT < xT && f.good(); dT++)
			{
				f.getline(b, 19, ' ');
	   			b[19] = '\0';
				img.set((T)strtod(b, &str_end), dt);
			}

		if (r < img.nrow() || dest < end) {
			printf("read_matrix(%s): stopped at row %d/%d, column %d/%d\n",
				fname, r, img.nrow(), img.ncol() - (int)(end - dest), img.ncol());
			return r;
		}
	}
	return 0;
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
	int b = 0, g = 0, r = 0, cG = 1 + g;
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
