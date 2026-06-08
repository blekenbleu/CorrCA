// https://statisticsbyjim.com/hypothesis-testing/how-to-find-p-value/

typedef unsigned char uchar;

struct MetricsStruct {};
typedef struct NamedType : MetricsStruct
{
	double RSS = 0, critical_value = 0, t_value[10] = { 0 };	// solution statistics
} Metrics;

#include "array.h"

	using namespace libNumerics;
	void inversematrix(matrix<double>& xinv, matrix<double>& x);
    void regress(Metrics& mm, matrix<double>& xinv, matrix<double>& coef, int c, matrix<double>& x, matrix<double>& y, uint ycol);
	void report(char *data, char *plotfile, matrix<double> &x, matrix<double> &y, matrix<double> &coef);
	void matrix_shift(double &dy, double &dx, matrix<double> &coef, int color, double y, double x);
	unsigned char get_CR(matrix<uchar> &plane, double row, double column);
	void test_matrix_shift(matrix<double> &shifted, matrix<double> &x, matrix<double> &coef);

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
