int suspect(Metrics m)
{
	int i, j = 0;
	double low = m.critical_value;
	for (i = m.B.nrow() - 1; 0 <= i; i--)
	if (low > abs(m.t_value[i]))
		low = abs(m.t_value[j = i]);
	return j;
}

void pad9(char *parm)
{
	int i, j;
	for (i = 0, j = (int)strlen(parm); i < j; i++)
		putc(parm[i], stdout);
	for (; i < 9; i++)
		putc(' ', stdout);
}

void mprint(Metrics &m, char **factor, char *dep)
{
	printf("%s = %.3f", dep, m.B(0,0));
	for (int i = 1; i < 10; i++)
		printf(" + %.3f %s", m.B(i, 0), factor[i]);
	printf("\n");
	printf("%.3f Residuals Sum of Squares, %.3f T critical value\n"
			"Estimate T-value\n", m.RSS, m.critical_value);
	for (int i = 0; i < 10; i++)
		if (0 != m.B(i, 0))
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
  } else printf("could not write coefficients to %s\n", p);
}
