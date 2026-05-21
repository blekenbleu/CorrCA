int suspect(Metrics m)
{
	int i, j = 0;
	double low = m.critical_value;
	for (i = m.B.nrow() - 1; 0 <= i; i--)
	if (low > abs(m.t_value[i]))
		low = abs(m.t_value[j = i]);
	return j;
}

matrix<double> report(matrix<double> x, matrix<double> y, int col, char *dep, vector<int> v)
{
	char *l[] = { "intercept","xG","yG","xG2","yG2","xG3","yG3" };
	Metrics m = regress(x, y, col);
	int j;
	if (0 < (j = suspect(m)))
		printf("suspect %s %s t-value %f\n", dep, l[j], m.t_value[j]);

	printf("%s = %.3f", dep, m.B(0,0));
	for (int i = 1; i < v.size(); i++)
		printf(" + %.3f %s", m.B(i, 0), l[(i)]);
	printf("\n");
	printf("%.3f Residuals Sum of Squares, %.3f T critical value\nEstimate T-value\n", m.RSS, m.critical_value);
	int count = v.nrow();
	for (int i = 0; i < count; i++)
		if (0 != m.B(i, 0))
			printf("%s %.3f\n", l[v(i)], m.t_value[i]);

	return m.B;
}
