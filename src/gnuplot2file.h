#include "regress.h"
#include "report.h"

template <typename T>
void gnuplot2file(char *plotfile,	// red, green, blue centers
	vector<T> &xR, vector<T> &yR, vector<T> &xG, vector<T> &yG, vector<T> &xB, vector<T> &yB,
	image_double &imgR, image_double &imgG)
{
	uint len = (int)strlen(plotfile);
	len += 6;
	char *fsn = (char *)calloc(len, sizeof(char));
	if (0 != fsn)
	{
		sprintf(fsn, "%s.gp", plotfile);
		if (FILE *gnuplot = fopen(fsn, "wt"))
		{
				fprintf(gnuplot, "set title 'red and blue differences from green centers'\n"
						"set grid\n unset xrange\n unset yrange\n unset zrange\n"
						"unset zeroaxis\n"
						"set xlabel 'green center pixel column' rotate parallel\n"
						"set ylabel 'green center pixel row' rotate parallel\n"
						"set zlabel 'red, blue center differences' rotate parallel \n\n"
						"set datafile separator ' ,'\n\n");

			sprintf(fsn, "%sG.txt", plotfile);
   		 	if (FILE *txtplot = fopen(fsn, "wt"))
			{
				len = xR.size();
				// create and populate regress() input
				matrix<double> x = matrix<T>(len, 7), y = matrix<T>(len, 4);
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

				fprintf(gnuplot, 
						"splot '%s' using 1:2:3 with points pt 7 ps 0.5 lc rgb 'orange' title 'red x',\\\n", fsn);
/*
				fprintf(gnuplot,
					  	"'%s' using 1:2:4 with points pt 6 ps 0.7 lc rgb 'magenta' title 'red y', \\\n"
					  	"'%s' using 1:2:5 with points pt 7 ps 0.5 lc rgb 'blue' title 'blue x', \\\n"
					  	"'%s' using 1:2:6 with points pt 6 ps 0.7 lc rgb 'cyan' title 'blue y'", fsn, fsn, fsn);
 */

				printf("\nfit coefficients for column 0 of y\n");
				matrix<double> B = report(regress(x, y, 0), "dxR");
				fprintf(gnuplot, "%.3f + %.3f*x + %.3f*y + %.3f*x*x + %.3f*y*y + %.3f*x*x*x + %.3f*y*y*y\n",
						B(0,0), B(1,0), B(2,0), B(3,0), B(4,0), B(5,0), B(6,0));
			} else printf("cannot open file %s\n", fsn);
			fclose(gnuplot);
		} else printf("cannot open file %s\n", fsn);
		free(fsn);
	}

	printf(" done.");
}

