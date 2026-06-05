#include "regress.h"
#include "report.h"
#include "correct.cpp"
#include "test_read_matrix.h"

// included in main_centering.cpp for visibility to array.h templates
// about class templates:  https://isocpp.org/wiki/faq/templates

template <typename T>
void gnuplot2file(char *plotfile,	// red, green, blue centers
	vector<T> &xR, vector<T> &yR, vector<T> &xG, vector<T> &yG,
	vector<T> &xB, vector<T> &yB,
	image_char &imgR, image_char &imgG, matrix<T> &coef)
{
	// create and populate regress() input
	uint len = xR.size();
	// set polynomial coefficient count 10
	// and polynomial count 4 { XR, YR, XB, YB }
	matrix<T> x(len, 10), y(len, 4);
	double scale = xG[0];
	scale /= xR[0];
	scale = (1.5 < scale) ? 2.0 : 1.0;
	// scale range of green pixel centers [0:1]
	double xm = imgG->xsize, ym = imgG->ysize;
	for (uint i = 0; i < len; i++)
	{
		T xGi = xG[i], yGi = yG[i];
		T xg1 = xGi / xm, yg1 = yGi / ym;	// rescaled [0:1]
		double x2 = xg1 * xg1, y2 = yg1 * yg1;
		T sxR = scale * xR[i];
		T syR = scale * yR[i];
		T sxB = scale * xB[i], syB = scale * yB[i];
		// x, y matrices for regress(), called in report();
		x(i, 0) = 1.0; // x(i, 0) are intercepts
		x(i, 1) = xg1; x(i, 2) = yg1; x(i, 3) = x2;	x(i, 4) = y2; x(i, 5) = xg1*x2;
		x(i, 6) = yg1*y2; x(i, 7) = xg1*yg1; x(i, 8) = x2*yg1; x(i, 9) = xg1*y2;
		y(i, 0) = sxR - xGi; y(i, 1) = syR - yGi;
		y(i, 2) = sxB - xGi; y(i, 3) = syB - yGi;
	}

	char fsn[180]{};	// fully qualified filename buffer
	sprintf(fsn, "%sG.txt", plotfile);
	// calculate coef for independent variables x (pixel coordinates)
	// that best fit measured y (red, blue misregistrations)
	// write solutions to fsn derived from plotfile
	report(fsn, plotfile, x, y, coef);
	matrix<T> shifted(len, 4);	// red, blue corrections by matrix_shift()
	// shifted are results for x(,1) x(,2) test spot pixel pixel coordinates using coef
	test_matrix_shift(shifted, x, coef);

	sprintf(fsn, FOLDER "%sG.txt", plotfile);
	if (FILE *txtplot = fopen(fsn, "wt"))
	{
		char *gfmt = "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f \n";
		
		printf("\nSaving spot centers and scaled pixel polynomial values to gnuplot file... ");
//		fprintf(txtplot, "# rows %d\n", (uint)len);
		// gnuplot: green x, y centers; red center diffs x, y; blue diffs x,y
		fprintf(txtplot, "xG,yG,dxR,dyR,dxB,dyB,xG2,yG2,xG3,yG3,xyG,xxyG,xyyG,-sxR,-syR,-sxB,-syB\n");
		double scale = xG[0], dxR, dyR, dxB, dyB;
		scale /= xR[0];
		scale = (1.5 < scale) ? 2.0 : 1.0;
		for (uint i = 0; i < len; i++)
			fprintf(txtplot, gfmt, x(i, 1), x(i, 2), y(i, 0), y(i, 1),
					y(i, 2), y(i, 3), x(i, 3), x(i, 4), x(i, 5), x(i, 6),
					x(i, 7), x(i, 8), x(i, 9), dxR = shifted(i, 0), dyR = shifted(i, 1), dxB = shifted(i, 2), dyB = shifted(i, 3));
		fclose(txtplot);
	} else printf("gnuplot2file():  cannot open file %s\n", fsn);

	printf(" done.\n");
}

