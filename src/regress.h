// https://statisticsbyjim.com/hypothesis-testing/how-to-find-p-value/

	typedef unsigned char uchar;

	struct MetricsStruct {};
	typedef struct NamedType : MetricsStruct
	{
		double RSS = 0, critical_value = 0, t_value[10] = { 0 };	// solution statistics
	} Metrics;

/*
	double T_value(double sample_mean, double null_value,
					double stddev, int sample_size)
	{
		return (sample_mean - null_value) / (stddev / sqrt(sample_size));
	}

	static double mean(matrix<double> m, unsigned int column)
	{
		double sum = 0;
		column %= m.ncol();

		for (int i = 0; i < m.nrow(); i++)
			sum += m(i, column);
		return sum / m.nrow();
	}
 */
