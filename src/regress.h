#include "t_test.h"

double mean(matrix<double> m, unsigned int column)
{
    double sum = 0;
    column %= m.ncol();

    for (int i = 0; i < m.nrow(); i++)
        sum += m(i, column);
    return sum / m.nrow();
}

matrix<double> multimatrix(matrix<double> a, matrix<double> b)
{
  if(a.ncol() != b.nrow()) {
    printf("Can't multiply matrices");
    return a;
  }

  else {
    matrix<double> result = matrix<double>(a.nrow(), b.ncol());

    for (int i = 0; i < a.nrow(); i++) {
        for (int j = 0; j < b.ncol(); j++) {
            double sum = 0;
            for (int k = 0; k < a.ncol(); k++)
                sum = sum + a(i, k) * b(k, j);
            result(i, j) = sum;
        }
    }

    return result;
  }
}

matrix<double> multitransmatrix(matrix<double> b)
{
    matrix<double> result = matrix<double>(b.ncol(), b.ncol());

    for (int i = 0; i < b.ncol(); i++)
	{
        for (int j = 0; j < b.ncol(); j++)
		{
            double sum = 0;
            for (int k = 0; k < b.nrow(); k++)
                sum += b(k, i) * b(k, j);
            result(i, j) = sum;
        }
    }

    return result;
}

// return a single column matrix for b(, bindex) 
matrix<double> multitransmatrix(matrix<double> a, matrix<double> b, uint bindex)
{
  if(a.nrow() != b.nrow()) {
    printf("Can't multiply matrices");
    return a;
  }

  matrix<double> result = matrix<double>(a.ncol(), 1);

  for (int i = 0; i < a.ncol(); i++) {
    double sum = 0;
    for (int k = 0; k < a.nrow(); k++)
      sum = sum + a(k, i) * b(k, bindex);
    result(i, 0) = sum;
  }

  return result;
}

int squarematrix(matrix<double> m, double square[25][25])
{
  if(m.nrow() != m.ncol() || m.nrow() > 24)
  {
    printf("Matrix isn't square!");
    return 0;
  }

  for(int i = 0; i < m.nrow(); i++)
      for(int j = 0; j < m.nrow(); j++)
        square[i][j] = m(i, j);

  return m.nrow();
}

double determinant(double a[25][25], double k) {
    double s = 1, det = 0, b[25][25] = { 0 };
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

void trans(double num[25][25], double fac[25][25], double r) {
    int i, j = 0;
    double b[25][25] = { 0 }, inv[25][25] = { 0 }, d;
    for (i = 0; i < r; i++) {
        for (j = 0; j < r; j++) {
            b[i][j] = fac[j][i];
        }
    }
    d = determinant(num, r);
    inv[i][j] = 0;
    for (i = 0; i < r; i++) {
        for (j = 0; j < r; j++) {
            inv[i][j] = b[i][j] / d;
        }
    }

  for (i = 0; i < r; i++) {
        for (j = 0; j < r; j++) {
            num[i][j] = inv[i][j];
        }
    }
}

// https://www.cuemath.com/algebra/cofactor-matrix/
void cofactors(double num[25][25], double f) {
    double b[25][25] = { 0 }, fac[25][25] = { 0 };
    int p, q, m, n, i, j;
    for (q = 0; q < f; q++) {
        for (p = 0; p < f; p++) {
            for (m = n = i = 0; i < f; i++) {
                for (j = 0; j < f; j++) {
                    b[i][j] = 0;
                    if (i != q && j != p) {
                        b[m][n] = num[i][j];
                        if (n < (f - 2))
                               n++; else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    }
    trans(num, fac, f);
}

matrix<double> inversematrix(matrix<double> m)
{
  double squareTemp[25][25];
  memset(squareTemp, 0, 625 * sizeof(double)); // 25 * 25
  matrix<double> result = matrix<double>(m.nrow(), m.nrow());
  int n = squarematrix(m, squareTemp);
  double d = determinant(squareTemp, n);

  // printf("\nThe determinant is: %.0f", d);

  if (d == 0)
      printf("\ninverseMatrix(): MATRIX IS NOT INVERSIBLE\n");
  else cofactors(squareTemp, n);

  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      result(i, j) = squareTemp[i][j];

  return result;
}

// fit x to column yindex of y
Metrics regress(matrix<double> x, matrix<double> xinv,
				matrix<double> y, uint yindex)
{
  int dof = x.nrow() - x.ncol();
  double yMean = mean(y, yindex);
  Metrics modelMetrics = { 0 };

  // a column of (up to 11) coefficients for y column yindex
  modelMetrics.B = multimatrix(xinv, multitransmatrix(x, y, yindex));
  matrix<double> Yhat = multimatrix(x, modelMetrics.B);  //  y estimates

  // Residuals Sum of Squares (RSS):  Unexplained Variance
  modelMetrics.RSS = 0;
  for (int i = 0; i < y.nrow(); i++)
  {
      double diff = y(i, yindex) - Yhat(i, 0);
      modelMetrics.RSS += diff * diff;
  }
  modelMetrics.critical_value = critical_value(x.nrow() - 1);

// prune covariants with p-value > 0.1 (t-value < ~1.65 for hundreds of samples)
// https://www.statology.org/how-to-calculate-a-p-value-from-a-t-test-by-hand/
// t-test statistic = Model coefficient / regression coefficient standard error
  modelMetrics.t_value[0] = modelMetrics.B(0, 0) / sqrt(xinv(0, 0) * modelMetrics.RSS / dof--);
  double md = modelMetrics.RSS / dof;
  for(int j = 1; j < x.ncol(); j++)	// independent variables
    modelMetrics.t_value[j] = modelMetrics.B(j, 0) / sqrt(xinv(j, j) * md);

  return modelMetrics;
}
