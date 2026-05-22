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

void mprint(Metrics m, char **l, vector<int> v, char *dep)
{
	printf("%s = %.3f", dep, m.B(0,0));
	for (int i = 1; i < v.size(); i++)
		printf(" + %.3f %s", m.B(i, 0), l[v(i)]);
	printf("\n");
	printf("%.3f Residuals Sum of Squares, %.3f T critical value\nEstimate T-value\n", m.RSS, m.critical_value);
	int count = v.nrow();
	for (int i = 0; i < count; i++)
		if (0 != m.B(i, 0))
		{
			pad9(l[v(i)]);
			printf(" %.3f\n", m.t_value[i]);
		}
}

matrix<double> report(matrix<double> x, matrix<double> y, int col, char *dep, vector<int> v)
{
	char *l[] = { "intercept","xG","yG","xG2","yG2","xG3","yG3" };
	Metrics m = regress(x, y, col), mj;
	matrix<double> xj, xk;
	vector <int> vj, vk;
	int j, k;

	mprint(m, l, v, dep);	
	if (0 < (j = suspect(m)))
	{
		printf("\nsuspect %s %s t-value %f\n", dep, l[j], m.t_value[j]);
		int row = x.nrow();
		xj.without(row, j, x);
	 	vj.without(j, v);
		mj = regress(xj, y, col);
		if (abs(mj.t_value[0]) > abs(m.t_value[0]))
			mprint(mj, l, vj, dep);	
		else {
			printf("\t no improvement:  old:new intercept T-value %f:%f\n", m.t_value[0], mj.t_value[0]);
			return m.B;
		}
	} else
		return m.B;

	if (0 < (k = suspect(mj)))
	{
		printf("\nsuspect %s %s t-value %f\n", dep, l[k], mj.t_value[k]);
		int row = xj.nrow();
		xk.without(row, k, xj);
	 	vk.without(k, vj);
		Metrics mk = regress(xk, y, col);
		if (abs(mk.t_value[0]) > abs(mj.t_value[0]))
		{
			mprint(mk, l, vk, dep);
			return mk.B;
		}
		else printf("\t no improvement:  old:new intercept T-value %f:%f\n", mj.t_value[0], mk.t_value[0]);
	}

	return mj.B;
}
