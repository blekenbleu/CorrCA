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

void mprint(Metrics m, char **factor, vector<int> v, char *dep)
{
	printf("%s = %.3f", dep, m.B(0,0));
	for (int i = 1; i < v.size(); i++)
		printf(" + %.3f %s", m.B(i, 0), factor[v(i)]);
	printf("\n");
	printf("%.3f Residuals Sum of Squares, %.3f T critical value\n"
			"Estimate T-value\n", m.RSS, m.critical_value);
	int count = v.nrow();
	for (int i = 0; i < count; i++)
		if (0 != m.B(i, 0))
		{
			pad9(factor[v(i)]);
			printf(" %.3f\n", m.t_value[i]);
		}
}

void report(double *coef, matrix<double> x, matrix<double> y,
			 int col, char *dep, vector<int> v, char *factor[])
{
	Metrics m;	//, mj;
	regress(m, x, y, col);
	mprint(m, factor, v, dep);
	for (int i = 0; i < x.ncol(); i++)
		coef[i] = m.B(i,0);
}
