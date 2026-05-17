#include "t_test.h"

double mean(matrix<double> m, unsigned int column)
{
    double sum = 0;
    column %= m.ncol();

    for (int i = 0; i < m.nrow(); i++)
        sum += m(i, column);
    return sum / m.nrow();
}

matrix<double> transmatrix(matrix<double> m)
{
  matrix<double> result = matrix<double>::zeros(m.ncol(), m.nrow());
  for(int j = 0; j < m.nrow(); j++)
    for(int i = 0; i < m.ncol(); i++)
        result(j, i) = m(i, j);
  return result;
}

matrix<double> multimatrix(matrix<double> a, matrix<double> b)
{
  if(a.ncol() != b.nrow()) {
    printf("Can't multiply matrices");
    return a;
  }

  else {
    matrix<double> result = matrix<double>::zeros(a.nrow(), b.ncol());

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
  matrix<double> result = matrix<double>::zeros(m.nrow(), m.nrow());
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

matrix<double> genCoefficients(matrix<double> x, matrix<double> y, matrix<double> xtrans)
{
  matrix<double> xmm = multimatrix(xtrans, x);
  matrix<double> xinv = inversematrix(xmm);
//free(xmm.data);
  matrix<double> xty = multimatrix(xtrans,y);
  matrix<double> B = multimatrix(xinv, xty);
//free(xinv.data);
//free(xty.data);
  return B;
}

matrix<double> stdErr(matrix<double> x, double errStdDev, matrix<double> xtrans)
{
  matrix<double> result = inversematrix(multimatrix(xtrans, x));

  double est = errStdDev * errStdDev;
  for(int i = 0; i < x.ncol(); i++)
	for (int j = 0; j < x.ncol(); j++)
      result(i, j) *= est;

  return result;
}

Metrics regress(matrix<double> x, matrix<double> y)
{
  int i, j;
  double errStdDev = 0;
  double yMean = mean(y, 0);
  matrix<double> xtrans = transmatrix(x);
  // B: a 1-D vector of (up to 11) coefficients
  matrix<double> B = genCoefficients(x, y, xtrans);
  matrix<double> Yhat = multimatrix(x, B);  //  y estimates
  // errStdDev = sqrt(sum((y-Yhat)**2)) / (result.nrow() - x.ncol()));
  // Matrix residuals = calcResiduals(x, y, Yhat, &errStdDev);
  matrix<double> stdErrmatrix = stdErr(x, errStdDev, xtrans);
  // predicted values:  x.nrow()
  Metrics modelMetrics = { 0 };

  // Residuals Sum of Squares (RSS):  Unexplained Variance
  modelMetrics.RSS = 0;
  for (i = 0; i < y.nrow(); i++)
  {
	double d = y(i, 0) - Yhat(i, 0);
    modelMetrics.RSS += d * d;
  }

  // Generate Coefficient Metrics
  modelMetrics.critical_value = critical_value(x.nrow() - 1);

  for(i = 0, j = 1; j < x.ncol(); j++)	// independent variables
  {
    // Standard error for regression coefficients
    double StdErr = sqrt(stdErrmatrix(i, j));
    // t-test statistic = Model coefficient / regression coefficient standard error
    modelMetrics.t_test[i++] = B(j, 0) / StdErr;
  }
  return modelMetrics;
}
