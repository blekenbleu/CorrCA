
template <typename T>
void gnuplot2file(char *plotfile, double scale,	// red, green, blue centers
	vector<T> &xR, vector<T> &yR, vector<T> &xG, vector<T> &yG, vector<T> &xB, vector<T> &yB)
{
	uint len = (int)strlen(plotfile);
	len += 6;
	char *fsn = (char *)calloc(len, sizeof(char));
	if (0 != fsn)
	{
		sprintf(fsn, "%s.gp", plotfile);
		if (FILE *gnuplot = fopen(fsn, "wt"))
		{
			sprintf(fsn, "%sG.txt", plotfile);
   		 	if (FILE *txtplot = fopen(fsn, "wt"))
			{
				printf("Saving uncorrected centers to gnuplot file... ");
				fprintf(txtplot, "# rows %d\n", (uint)(len = xR.size()));
				// gnuplot: green x, y centers; red center diffs x, y; blue diffs x,y
				fprintf(txtplot, "xG,yG,dxR,dyR,dxB,dyB,xG2,yG2,xG3,yG3,r,r2,r3\n");
    			char *gfmt = "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.4f,%.6f,%.4f,%.4f,%.4f,%.4f\n";
				double xm = xG[len - 1], ym = yG[len - 1];
				for (uint i = 0; i < len; i++)
				{
					T xg1 = xG[i] /xm, yg1 = yG[i] /ym;
					T r3, r2 = xg1 - 0.505, r = yg1 - 0.505;
					r2 = r*r + r2*r2;
					r = (T)sqrt(r2);
					r3 = r2 * r;
					T sxR = scale * xR[i], syR = scale * yR[i], sxB = scale * xB[i], syB = scale * yB[i];
					fprintf(txtplot, gfmt, xg1, yg1, sxR - xG[i], syR - yG[i], sxB - xG[i], syB - yG[i],
							xg1*xg1, yg1*yg1, xg1*xg1*xg1, yg1*yg1*yg1, r, r2, r3);
				}
				fclose(txtplot);

				fprintf(gnuplot, "set title 'red and blue differences from green centers'\n"
						"set grid\n unset xrange\n unset yrange\n unset zrange\n"
						"unset zeroaxis\n"
						"set xlabel 'green center pixel column' rotate parallel\n"
						"set ylabel 'green center pixel row' rotate parallel\n"
						"set zlabel 'red, blue center differences' rotate parallel \n\n"
						"set datafile separator ' ,'\n\n");

				fprintf(gnuplot, 
						"splot '%s' using 1:2:3 with points pt 7 ps 0.5 lc rgb 'orange' title 'red x',\\\n", fsn);
				fprintf(gnuplot,
					  	"'%s' using 1:2:4 with points pt 6 ps 0.7 lc rgb 'magenta' title 'red y', \\\n"
					  	"'%s' using 1:2:5 with points pt 7 ps 0.5 lc rgb 'blue' title 'blue x', \\\n"
					  	"'%s' using 1:2:6 with points pt 6 ps 0.7 lc rgb 'cyan' title 'blue y'", fsn, fsn, fsn);
			} else printf("cannot open file %s\n", fsn);
			fclose(gnuplot);
		} else printf("cannot open file %s\n", fsn);
		free(fsn);
	}

	printf(" done.");
}

