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

	char fsn[180]{};
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
		for (uint i = 0, j = 0, k = 0; i < len; i++)
		{
			T xGi = xG[i], yGi = yG[i];
			T xg1 = xGi / xm, yg1 = yGi / ym;	// rescaled [0:1]
			double x2 = xg1 * xg1, y2 = yg1 * yg1;
			T sxR = scale * xR[i];
			T syR = scale * yR[i];
			T sxB = scale * xB[i], syB = scale * yB[i];
			// x, y matrices for regress(), called in plane();
			x(j++) = 1.0; // x(i, 0) are intercepts
			fprintf(txtplot, gfmt, x(j++) = xg1, x(j++) = yg1,
					y(k++) = sxR - xGi,		y(k++) = syR - yGi,
					y(k++) = sxB - xGi,		y(k++) = syB - yGi,
					x(j++) = x2,			x(j++) = y2,
					x(j++) = xg1*x2,		x(j++) = yg1*y2,
					x(j++) = xg1*yg1,    	x(j++) = x2*yg1,
					x(j++) = xg1*y2);
		}
		fclose(txtplot);
		sprintf(fsn, "%sG.txt", plotfile);
		plane(fsn, plotfile, x, y, coef);
	} else printf("gnuplot2file():  cannot open file %s\n", fsn);

	printf(" done.");
}

