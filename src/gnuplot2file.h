#include "regress.h"
#include "report.h"

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

void plane(char *data, char *fn, matrix<double> x, matrix<double> y, char *color, char axis, int col)
{
	char cx[10] = { '\0' };
	sprintf(cx, "%s %c", color, axis);
	printf("\nfit %s coefficients for column %d of y\n", cx, col);
	char fsn[100] = { '\0' };
	vector<int> ix = vector<int>::index(7);;
	sprintf(fsn, FOLDER "%s.gp", fn);
	if (FILE *gnuplot = fopen(fsn, "wt"))
	{
		fprintf(gnuplot, gph, cx, cx);
		matrix<double> B = report(regress(x, y, col), fn, ix);
		fprintf(gnuplot,
					"splot '%s' using 1:2:%d with points pt 7 ps 0.5 lc rgb '%s' title '%s',\\\n",
					data, 3 + col, color, cx);
		fprintf(gnuplot, "%.3f + %.3f*x + %.3f*y + %.3f*x*x + %.3f*y*y + %.3f*x*x*x + %.3f*y*y*y\n",
				B(0,0), B(1,0), B(2,0), B(3,0), B(4,0), B(5,0), B(6,0));
		fclose(gnuplot);
	} else printf("cannot open file %s\n", fsn);
}

template <typename T>
void gnuplot2file(char *plotfile,	// red, green, blue centers
	vector<T> &xR, vector<T> &yR, vector<T> &xG, vector<T> &yG, vector<T> &xB, vector<T> &yB,
	image_double &imgR, image_double &imgG)
{
	uint len = 6 + (uint)strlen(plotfile);
	char *fsn = (char *)calloc(len, sizeof(char));
	if (0 == fsn)
	{
		printf("gnuplot2file(): calloc(%d) failed\n", len);
		return;
	}

	len = xR.size();
	// create and populate regress() input
	matrix<double> x = matrix<T>(len, 7), y = matrix<T>(len, 4);

	sprintf(fsn, "%sG.txt", plotfile);
	if (FILE *txtplot = fopen(fsn, "wt"))
	{
		char *gfmt = "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.4f,%.6f,%.4f\n";
		// scale range of green pixel centers [0:1]
		double xm = imgG->xsize, ym = imgG->ysize;
		double scale = imgG->ysize;
		
		scale /= imgR->ysize;	// green plane may be 2x red, blue
		printf("Saving uncorrected centers to gnuplot file... ");
		fprintf(txtplot, "# rows %d\n", (uint)len);
		// gnuplot: green x, y centers; red center diffs x, y; blue diffs x,y
		fprintf(txtplot, "xG,yG,dxR,dyR,dxB,dyB,xG2,yG2,xG3,yG3\n");

		for (uint i = 0; i < len; i++)
		{
			T xg1 = xG[i] /xm, yg1 = yG[i] /ym, sxR = scale * xR[i];
			T syR = scale * yR[i], sxB = scale * xB[i], syB = scale * yB[i];
			x(i, 0) = 1.0;	// x, y matrix for regress; first x column is intercept
			fprintf(txtplot, gfmt, x(i, 1) = xg1, x(i, 2) = yg1,
					y(i, 0) = sxR - xG[i],	y(i, 1) = syR - yG[i],
					y(i, 2) = sxB - xG[i],	y(i, 3) = syB - yG[i],
					x(i, 3) = xg1*xg1,		x(i, 4) = yg1*yg1,
					x(i, 5) = xg1*xg1*xg1,	x(i, 6) = yg1*yg1*yg1);
		}
		fclose(txtplot);
		plane(fsn, "dxR", x, y, "red", 'x', 0);
		plane(fsn, "dyR", x, y, "red", 'y', 1);
		plane(fsn, "dxB", x, y, "blue", 'x', 2);
		plane(fsn, "dyB", x, y, "blue", 'y', 3);
	} else printf("gnuplot2file():  cannot open file %s\n", fsn);
	free(fsn);

	printf(" done.");
}

