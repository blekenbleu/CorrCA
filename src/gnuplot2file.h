#include "regress.h"
#include "report.h"
#include "correct.cpp"
#include "test_read_matrix.h"

// included in main_centering.cpp for visibility to array.h templates
// about class templates:  https://isocpp.org/wiki/faq/templates

char *gph =
{
	"set title '%s differences from green centers'\n"
	"set grid\n unset xrange\n unset yrange\n unset zrange\n"
	"unset zeroaxis\n"
	"set xrange [0:1]\n"
	"set yrange [0:1]\n"
	"set xlabel 'green center pixel column' rotate parallel\n"
	"set ylabel 'green center pixel row' rotate parallel\n"
	"set zlabel '%s center differences' rotate parallel \n\n"
	"set datafile separator ' ,'\n\n"
};

void plane(char *data, char *plotfile, matrix<double> &x, matrix<double> &y, matrix<double> &coef)
{
  char *factor[] = {"intercept", "x", "y", "x*x", "y*y", "x*x*x", "y*y*y", "x*y", "x*x*y", "x*y*y"};
  char *colors[] = { "red", "blue" }, fsn[100] = { '\0' };

  for (int col = 0; col < 4; col++)
  {
	int i = 1 & (col >> 1);
	char *color = colors[i];
    char axis = "xy"[1 & col];
	char dep[] = "dxR";
	dep[1] = axis; dep[2] = (1 == i) ? 'B' : 'R';
	char cx[10] = { '\0' };
	sprintf(cx, "%s d%c", color, axis);
	printf("\nfit %s coefficients for column %d of y\n", cx, col);
	sprintf(fsn, FOLDER "%s.gp", dep);
	if (FILE *gnuplot = fopen(fsn, "wt"))
	{
		int c10 = col*10;
		fprintf(gnuplot, gph, cx, cx);
		report(coef.data(c10), x, y, col, dep, factor);
		Metrics m; regress(m, x, y, col);
		mprint(m, factor, dep);
		for (int c = 0; c < x.ncol(); c++)
			coef(c + c10) = m.B(c,0);		// set poly coefficients
		fprintf(gnuplot,
				"splot '%s' using 1:2:%d with points"
				" pt 7 ps 0.5 lc rgb '%s' title '%s',\\\n",
				data, 3 + col, color, cx);
		fprintf(gnuplot, "%.3f", coef(c10));
		for (i = 1; i < x.ncol(); i++)
			fprintf(gnuplot, " + %.3f*%s", coef(c10 + i), factor[i]);
 		fprintf(gnuplot, "\n");
		fclose(gnuplot);
	} else printf("cannot open file %s\n", fsn);
  }
  // while we're at it, write the polynomial coefficients
  char* p = strrchr(plotfile, '/'), copy[50] = { '\0' };
  if (NULL == p)
	p = plotfile;
  char *dot = strrchr(p, '.');
  if (NULL != dot) {
	  strncpy(copy, p, dot - p);
	  p = copy;
  }
  sprintf(fsn, FOLDER "Poly_%s.txt", p);
  std::ofstream f(fsn, std::ios::binary);
  if (f) {
  	printf("writing coefficients to %s\n", p);
	char *hdr[] = { "XR", "YR", "XB", "YB" };
	double cf;
	for (int col = 0, c = 0; col < 4; col++) {
		sprintf(fsn, "# poly%s(x,y) :\n", hdr[col]);
		f.write(fsn, strlen(fsn));
		sprintf(fsn, "%.16g\n", coef(c++));
		f.write(fsn, strlen(fsn));
		for(int r = 1; r < 10; r++) {
			if (0 > (cf = coef(c++)))
				sprintf(fsn, "%.16g * %s\n", cf, factor[r]);
			else sprintf(fsn, "+ %.16g * %s\n", cf, factor[r]);
			f.write(fsn, strlen(fsn));
		}
	}
	f.close();
  }
}

template <typename T>
void gnuplot2file(char *plotfile,	// red, green, blue centers
	vector<T> &xR, vector<T> &yR, vector<T> &xG, vector<T> &yG,
	vector<T> &xB, vector<T> &yB,
	image_char &imgR, image_char &imgG, matrix<T>coef)
{
	uint len = 16 + (uint)strlen(plotfile);
	char *fsn = (char *)calloc(len, sizeof(char));
	if (0 == fsn)
	{
		printf("gnuplot2file(): calloc(%d) failed\n", len);
		return;
	}

	len = xR.size();
	// create and populate regress() input
	matrix<double> x = matrix<T>(len, 10), y = matrix<T>(len, 4);

	sprintf(fsn, FOLDER "%sG.txt", plotfile);
	if (FILE *txtplot = fopen(fsn, "wt"))
	{
		char *gfmt = "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f \n";
		// scale range of green pixel centers [0:1]
		double xm = imgG->xsize, ym = imgG->ysize;
		
		printf("\nSaving uncorrected centers to gnuplot file... ");
//		fprintf(txtplot, "# rows %d\n", (uint)len);
		// gnuplot: green x, y centers; red center diffs x, y; blue diffs x,y
		fprintf(txtplot, "xG,yG,dxR,dyR,dxB,dyB,xG2,yG2,xG3,yG3,xyG,xxyG,xyyG\n");
		double scale = xG[0];
		scale /= xR[0];
		scale = (1.5 < scale) ? 2.0 : 1.0;

		for (uint i = 0; i < len; i++)
		{
			T xGi = xG[i], yGi = yG[i];
			T xg1 = xGi / xm, yg1 = yGi / ym;	// rescaled [0:1]
			T sxR = scale * xR[i];
			T syR = scale * yR[i];
			T sxB = scale * xB[i], syB = scale * yB[i];
			// x, y matrices for regress(), called in plane();
			x(i, 0) = 1.0; // x(i, 0) are intercepts
			fprintf(txtplot, gfmt, x(i, 1) = xg1, x(i, 2) = yg1,
					y(i, 0) = sxR - xGi,	y(i, 1) = syR - yGi,
					y(i, 2) = sxB - xGi,	y(i, 3) = syB - yGi,
					x(i, 3) = xg1*xg1,		x(i, 4) = yg1*yg1,
					x(i, 5) = xg1*xg1*xg1,	x(i, 6) = yg1*yg1*yg1,
					x(i, 7) = xg1*yg1,    	x(i, 8) = xg1*xg1*yg1,
					x(i, 9) = xg1*yg1*yg1);
		}
		fclose(txtplot);
		sprintf(fsn, "%sG.txt", plotfile);
		plane(fsn, plotfile, x, y, coef);
	} else printf("gnuplot2file():  cannot open file %s\n", fsn);
	free(fsn);

	printf(" done.");
}

