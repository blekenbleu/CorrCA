#define _CRT_SECURE_NO_DEPRECATE
#include <fstream>
#include "regress.h"
#include "pgm_io.h"	// FOLDER

void pad9(char *parm)
{
	int i, j;
	for (i = 0, j = (int)strlen(parm); i < j; i++)
		putc(parm[i], stdout);
	for (; i < 9; i++)
		putc(' ', stdout);
}

void mprint(Metrics &m, matrix<double> &coef, int c, char **factor)
{
	int ncoef = coef.ncol();
	for (int i = 1; i < ncoef; i++)
		printf(" + %.3f %s", coef(c++), factor[i]);
	printf("\n%.3f Residuals Sum of Squares, %.3f T critical value\n"
			"Estimate T-value\n", m.RSS, m.critical_value);
	for (int i = 0; i < ncoef; i++)
	{
		pad9(factor[i]);
		printf(" %.3f\n", m.t_value[i]);
	}
}

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

void report(char *data, char *plotfile, matrix<double> &x, matrix<double> &y, matrix<double> &coef)
{
  char* factor[] = { "intercept", "x", "y", "x*x", "y*y", "x*x*x", "y*y*y", "x*y", "x*x*y", "x*y*y" }; //, "x*x*y*y" };
  char *colors[] = { "red", "blue", "orange", "purple" }, fsn[100] = { '\0' };
  int ncoef = x.ncol(), np = y.ncol();
  matrix<double> xinv(x.ncol(), x.ncol());

//inversematrix(xinv, x);
  TestTimer timer;
  timer.tic();
  mpinvert(xinv, x);
  timer.toc("mpinvert(xinv, x)");
  printf("done.\nfit d[xy][RB] coefficients to corresponding y columns\n");
  // np polynomials for CA reduction
  for (int c = 0, poly = 0; poly < np; poly++)
  {
	int i = 1 & (poly >> 1);
	char *color = colors[i], *shift = colors[2+i];
    char axis = "xy"[1 & poly];
	char dep[] = "dxR";
	dep[1] = axis; dep[2] = (1 == i) ? 'B' : 'R';
	char cx[11] = { '\0' };
	sprintf(cx, "%s d%c", color, axis);
	// solve poly coefficients given i independent x(i,) and dependent y(i,poly)
	Metrics m; regress(m, xinv, coef, c, x, y, poly);
//	printf("%s = %.3f", dep, coef(c));
//	mprint(m, coef, c, factor);
	sprintf(fsn, FOLDER "%s.gp", dep);
	if (FILE *gnuplot = fopen(fsn, "wt"))
	{
		fprintf(gnuplot, gph, cx, cx);
		fprintf(gnuplot,
				"splot '%s' using 1:2:%d with points"
				" pt 7 ps 0.5 lc rgb '%s' title '%s',\\\n"
				"\t '%s' using 1:2:%d with points"
				" pt 7 ps 0.5 lc rgb '%s' title 'shifted %s',\\\n",
				data, 3 + poly, color, cx, data, 7 + poly, shift, cx);

		fprintf(gnuplot, "%.3f", coef(c++));
		for (i = 1; i < ncoef; i++)
			fprintf(gnuplot, " + %.3f*%s", coef(c++), factor[i]);

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
	for (int poly = 0, c = 0; poly < np; poly++) {
		sprintf(fsn, "# poly%s(x,y) :\n", hdr[poly]);
		f.write(fsn, strlen(fsn));
		sprintf(fsn, "%.16g\n", coef(c++));
		f.write(fsn, strlen(fsn));
		for(int n = 1; n < ncoef; n++) {
			if (0 > (cf = coef(c++)))
				sprintf(fsn, "%.16g * %s\n", cf, factor[n]);
			else sprintf(fsn, "+%.16g * %s\n", cf, factor[n]);
			f.write(fsn, strlen(fsn));
		}
	}
	f.close();
  } else printf("could not write coefficients to %s\n", p);
}

/* solved in reduceCA.cpp matrix_shift()
 ; dxR = coef(0,0) + coef(0,1)*x + coef(0,2)*y + coef(0,3)*x*x + coef(0,4)*y*y
 ;		+ coef(0,5)*x*x*x + coef(0,6)*y*y*y + coef(0,7)*x*y + coef(0,8)*x*x*y
 ;		+ coef(0,9)*x*y*y // + coef(0,10)*x*x*y*y;
 */
void gnuplot2file(char *plotfile,	// red, green, blue centers
	vector<double> &xR, vector<double> &yR, vector<double> &xG, vector<double> &yG,
	vector<double> &xB, vector<double> &yB,
	image_char &imgR, image_char &imgG, matrix<double> &coef)
{
	// create and populate regress() input
	uint len = xR.size();
	matrix<double> x(len, coef.ncol()), y(len, coef.nrow());
	double scale = xG[0];
	scale /= xR[0];
	scale = (1.5 < scale) ? 2.0 : 1.0;
	// scale range of green pixel centers [0:1]
	double xm = imgG->xsize, ym = imgG->ysize;
	for (uint i = 0; i < len; i++)
	{
		double xGi = xG[i], yGi = yG[i];
		double xg1 = xGi / xm, yg1 = yGi / ym;	// rescaled [0:1]
		double x2 = xg1 * xg1, y2 = yg1 * yg1;
		double sxR = scale * xR[i], syR = scale * yR[i];
		double sxB = scale * xB[i], syB = scale * yB[i];
		// x, y matrices for regress(), called in report();
		x(i, 0) = 1.0; // x(i, 0) are intercepts
		x(i, 1) = xg1; x(i, 2) = yg1; x(i, 3) = x2;	x(i, 4) = y2; x(i, 5) = xg1*x2;
		x(i, 6) = yg1*y2; x(i, 7) = xg1*yg1; x(i, 8) = x2*yg1; x(i, 9) = xg1*y2;
//		x(i, 10) = x2*y2;
		y(i, 0) = sxR - xGi; y(i, 1) = syR - yGi;
		y(i, 2) = sxB - xGi; y(i, 3) = syB - yGi;
	}

	char fsn[180]{};	// fully qualified filename buffer
	sprintf(fsn, "%sG.txt", plotfile);
	// calculate coef for independent variables x (pixel coordinates)
	// that best fit measured y (red, blue misregistrations)
	// write solutions to fsn derived from plotfile
	report(fsn, plotfile, x, y, coef);
	matrix<double> shifted(len, 4);	// red, blue corrections by matrix_shift()
	// shifted are results for x(,1) x(,2) test spot pixel pixel coordinates using coef
	test_matrix_shift(shifted, x, coef);

	sprintf(fsn, FOLDER "%sG.txt", plotfile);
	if (FILE *txtplot = fopen(fsn, "wt"))
	{
		char *gfmt = "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f \n";
		
		printf("\nSaving spot centers and scaled pixel polynomial values to gnuplot file... ");
		fprintf(txtplot, "# rows %d\n", (uint)len);
		// gnuplot: green x, y centers; red center diffs x, y; blue diffs x,y, shifted xR, yR, xB, yB
		fprintf(txtplot, "xG,yG,dxR,dyR,dxB,dyB,sxR,syR,sxB,syB\n");
		double scale = xG[0];
		scale /= xR[0];
		scale = (1.5 < scale) ? 2.0 : 1.0;
		for (uint i = 0; i < len; i++)
			fprintf(txtplot, gfmt, x(i, 1), x(i, 2), y(i, 0), y(i, 1), y(i, 2), y(i, 3),
					shifted(i, 0), shifted(i, 1), shifted(i, 2), shifted(i, 3));
		fclose(txtplot);
	} else printf("gnuplot2file():  cannot open file %s\n", fsn);

	printf(" done.\n");
}
