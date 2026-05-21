matrix<double> report(Metrics m, char *dep, vector<int> v)
{
	char *l[] = { "intercept","xG","yG","xG2","yG2","xG3","yG3" };

	printf("%s = %.3f", dep, m.B(0,0));
	for (int i = 1; i < v.size(); i++)
		printf(" + %.3f %s", m.B(i, 0), l[(i)]);
	printf("\n");
	printf("%.3f Residuals Sum of Squares, %.3f T critical value\n\t T-value\n", m.RSS, m.critical_value);
	int count = v.nrow();
	for (int i = 0; i < count; i++)
		if (0 != m.B(i, 0))
			printf("%s %.3f\n", l[v(i)], m.t_value[i]);

	return m.B;
}
