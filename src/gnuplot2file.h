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

