#include "t_test.h"

double mean(matrix<double> m, unsigned int column)
{
	double sum = 0;
	column %= m.ncol();

	for (int i = 0; i < m.nrow(); i++)
		sum += m(i, column);
	return sum / m.nrow();
}

// sum rows of xinv columns multiplied by result into coef
void vector_matrix_mult(matrix<double> &coef, int c, matrix<double> &inv, matrix<double> &result)
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

void vector_matrix_mult(matrix<double> &Yhat, matrix<double> &x, matrix<double> &coef, int c)
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

void multitransmatrix(matrix<double> &m, matrix<double> &b)
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
matrix<double> multitranscolumn(matrix<double> &result, matrix<double> &a, matrix<double> &b, uint column)
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

int squarematrix(matrix<double> &m, double (*square)[12])
{
  if (m.nrow() != m.ncol() || m.nrow() > 24)
  {
	printf("Matrix isn't square!");
	return 0;
  }

  for(int i = 0; i < m.nrow(); i++)
	  for(int j = 0; j < m.nrow(); j++)
		square[i][j] = m(i, j);

  return m.nrow();
}

double determinant(double (*a)[12], double k) {
	double s = 1, det = 0, b[12][12] = { 0 };
	if (k == 1) {
		return (a[0][0]);
	} else {
		int i, j, m, n, c;
		det = 0;
		for (c = 0; c < k; c++) {
			for (m =  n = i = 0; i < k; i++) {
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

double mdeterminant(matrix<double> &a, double k) {
	double s = 1, det = 0, b[12][12] = { 0 };
	if (k == 1) {
		return a(0);
	} else {
		int i, j, m, n, c;
		det = 0;
		for (c = 0; c < k; c++) {
			for (m =  n = i = 0; i < k; i++) {
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

void trans(matrix<double> &num, double (*fac)[12], double r) {
	double d = mdeterminant(num, r);

	for (int i = 0; i < r; i++)
		for (int j = 0; j < r; j++)
			num(i, j) = fac[j][i] / d;
}

// https://www.cuemath.com/algebra/cofactor-matrix/
void cofactors(matrix<double> &num, double f) {
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

void inversematrix(matrix<double> &xinv, matrix<double> &x)
{
  int ncol = x.ncol();
  multitransmatrix(xinv, x);
  if (0 == mdeterminant(xinv, ncol))
	  printf("\ninversematrix(): MATRIX IS NOT INVERSIBLE\n");
  else cofactors(xinv, ncol);
}

// fit coef vector to x for column ycol of y
void regress(Metrics &mm, matrix<double> &xinv, matrix<double> &coef, int c, matrix<double> &x, matrix<double> &y, uint ycol)
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
