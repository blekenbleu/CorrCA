void report(Metrics m, char *dep)
{
	char *l[] = { "intercept","xG","yG","xG2","yG2","xG3","yG3" };

	printf("%s = %.3f + %.3f %s +  %.3f %s +  %.3f %s +  %.3f %s +  %.3f %s +  %.3f  %s\n",
			dep, m.B(0, 0), m.B(1, 0), l[1], m.B(2, 0), l[2], m.B(3, 0), l[3],
			m.B(4, 0), l[4], m.B(5, 0), l[5], m.B(6, 0), l[6]);
	printf("%.3f Residuals Sum of Squares, %.3f T critical value\n\t T-value\n", m.RSS, m.critical_value);
	for (int i = 0; i < 7; i++)
		if (0 != m.B(i, 0))
			printf("%s %.3f\n", l[i], m.t_value[i]);
	return;
}
