#include "regress.h"
	
	float critical_value(unsigned int df)
	{

	// two-tailed significance 0.05
	// https://statisticsbyjim.com/hypothesis-testing/t-distribution-table/
	float t_table[] = { 12.71f, 4.303f, 3.182f, 2.776f, 2.571f, 2.447f,
		2.365f, 2.306f, 2.262f, 2.228f, 2.201f, 2.179f, 2.160f, 2.145f,
		2.131f, 2.120f, 2.110f, 2.101f, 2.093f, 2.086f, 2.080f, 2.074f,
		2.069f, 2.064f, 2.060f, 2.056f, 2.052f, 2.048f, 2.045f, 2.042f,
		2.021f, 2.000f, 1.990f, 1.984f, 1.962f, 1.960f };

		if (30 >= df)
			return t_table[df - 1];
		else if (df >= 1000)
			return t_table[34];
		else if (df >= 100)
			return t_table[33] + (df - 100) * (t_table[34] - t_table[33]) / 900;
		else if (df >= 80)
			return t_table[32] + (df - 80) * (t_table[33] - t_table[32]) / 20;
		else if (df >= 60)
			return t_table[31] + (df - 60) * (t_table[32] - t_table[31]) / 20;
		else if (df >= 40)
			return t_table[30] + (df - 40) * (t_table[31] - t_table[30]) / 20;
		else
			return t_table[29] + (df - 30) * (t_table[30] - t_table[29]) / 10;
	}

	// sum rows of xinv columns multiplied by result into coef
	void vector_matrix_mult(matrix<double>& coef, int c, matrix<double>& inv, matrix<double>& result)
	{
		double sum = 0;
		int ic = inv.ncol(), ir = inv.nrow();
		for (int i, r = 0; r < ir; r++)
		{
			for (sum = i = 0; i < ic; i++)
				sum += inv(r, i) * result(i);
			coef(c++) = sum;
		}
	}

	void vector_matrix_mult(matrix<double>& Yhat, matrix<double>& x, matrix<double>& coef, int c)
	{
		double sum = 0;
		int xc = x.ncol(), xr = x.nrow();
		for (int i, r = 0; r < xr; r++)
		{
			for (sum = i = 0; i < xc; i++)
				sum += x(r, i) * coef(c++);
			Yhat(r) = sum;
		}
	}

	void multitransmatrix(matrix<double>& m, matrix<double>& b)
	{
		double sum = 0;
		int nc = b.ncol(), nr = b.nrow();

		for (int k, j, i = 0; i < nc; i++)
			for (j = 0; j < nc; j++)
			{
				for (sum = k = 0; k < nr; k++)
					sum += b(k, i) * b(k, j);
				m(i, j) = sum;
			}
	}

	// return a single column matrix for b(, bindex) 
	matrix<double> multitranscolumn(matrix<double>& result, matrix<double>& a, matrix<double>& b, uint column)
	{
		if (a.nrow() != b.nrow()) {
			printf("Can't multiply matrices");
			return a;
		}

		for (int i = 0; i < a.ncol(); i++) {
			double sum = 0;
			for (int k = 0; k < a.nrow(); k++)
				sum = sum + a(k, i) * b(k, column);
			result(i) = sum;
		}

		return result;
	}

	int squarematrix(matrix<double>& m, double (*square)[12])
	{
		if (m.nrow() != m.ncol() || m.nrow() > 24)
		{
			printf("Matrix isn't square!");
			return 0;
		}

		for (int i = 0; i < m.nrow(); i++)
			for (int j = 0; j < m.nrow(); j++)
				square[i][j] = m(i, j);

		return m.nrow();
	}

	double determinant(double (*a)[12], double k) {
		double s = 1, det = 0, b[12][12] = { 0 };
		if (k == 1) {
			return (a[0][0]);
		}
		else {
			int i, j, m, n, c;
			det = 0;
			for (c = 0; c < k; c++) {
				for (m = n = i = 0; i < k; i++) {
					for (j = 0; j < k; j++) {
						b[i][j] = 0;
						if (i != 0 && j != c) {
							b[m][n] = a[i][j];
							if (n < (k - 2))
								n++;
							else {
								n = 0;
								m++;
							}
						}
					}
				}
				det = det + s * (a[0][c] * determinant(b, k - 1));
				s = -1 * s;
			}
		}
		return det;
	}

	double mdeterminant(matrix<double>& a, double k) {
		double s = 1, det = 0, b[12][12] = { 0 };
		if (k == 1) {
			return a(0);
		}
		else {
			int i, j, m, n, c;
			det = 0;
			for (c = 0; c < k; c++) {
				for (m = n = i = 0; i < k; i++) {
					for (j = 0; j < k; j++) {
						b[i][j] = 0;
						if (i != 0 && j != c) {
							b[m][n] = a(i, j);
							if (n < (k - 2))
								n++;
							else {
								n = 0;
								m++;
							}
						}
					}
				}
				det = det + s * (a(c) * determinant(b, k - 1));
				s = -1 * s;
			}
		}
		return det;
	}

	void trans(matrix<double>& num, double (*fac)[12], double r) {
		double d = mdeterminant(num, r);

		for (int i = 0; i < r; i++)
			for (int j = 0; j < r; j++)
				num(i, j) = fac[j][i] / d;
	}

	// https://www.cuemath.com/algebra/cofactor-matrix/
	void cofactors(matrix<double>& num, double f) {
		double b[12][12] = { 0 }, fac[12][12] = { 0 };
		int p, q, m, n, i, j;
		for (q = 0; q < f; q++) {
			for (p = 0; p < f; p++) {
				for (m = n = i = 0; i < f; i++)
					for (j = 0; j < f; j++) {
						b[i][j] = 0;
						if (i != q && j != p) {
							b[m][n] = num(i, j);
							if (n < (f - 2))
								n++;
							else {
								n = 0;
								m++;
							}
						}
					}
				fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
			}
		}
		trans(num, fac, f);
	}

	void inversematrix(matrix<double>& xinv, matrix<double>& x)
	{
		int ncol = x.ncol();
		multitransmatrix(xinv, x);
		if (0 == mdeterminant(xinv, ncol))
			printf("\ninversematrix(): MATRIX IS NOT INVERSIBLE\n");
		else cofactors(xinv, ncol);
	}

	// fit coef vector to x for column ycol of y
	void regress(Metrics& mm, matrix<double>& xinv, matrix<double>& coef, int c, matrix<double>& x, matrix<double>& y, uint ycol)
	{
		int nrow = x.nrow(), ncol = x.ncol();
		if (nrow != y.nrow()) {
			printf("Can't multiply matrices");
			return;
		}

		matrix<double> result(ncol, 1);
		multitranscolumn(result, x, y, ycol);
		// a column of coefficients for y column ycol
		vector_matrix_mult(coef, c, xinv, result);

		/* generate statistics (or not...)
		matrix<double> Yhat(nrow, 1);		//  y estimates
		vector_matrix_mult(Yhat, x, coef, c);

		// Residuals Sum of Squares (RSS):  Unexplained Variance
		mm.RSS = 0;
		for (int i = 0; i < nrow; i++)
		{
		  double diff = y(i, ycol) - Yhat(i, 0);
		  mm.RSS += diff * diff;
		}
		mm.critical_value = critical_value(nrow - 1);

	  // doomed:  pruning covariants with p-value > 0.1 (t-value < ~1.65 for hundreds of samples)
	  // https://www.statology.org/how-to-calculate-a-p-value-from-a-t-test-by-hand/
	  // t-test statistic = Model coefficient / regression coefficient standard error
		int dof = nrow - x.ncol();
		mm.t_value[0] = coef(c) / sqrt(xinv(0, 0) * mm.RSS / dof--);
		double md = mm.RSS / dof;
		for(int j = 1; j < x.ncol(); j++)	// independent variables
		  mm.t_value[j] = coef(c++) / sqrt(xinv(j, j) * md);
	   */
	}
