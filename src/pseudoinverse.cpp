// Reference:
// Fast Computation of Moore-Penrose Inverse Matrices
// Pierre Courrieu, https://arxiv.org/abs/0804.4809

#include <vector>
#include "regress.h"

void showMatrix(const matrix<double> &matG, const char *describe)
{
    int nrows{matG.nrow()}, ncols{matG.ncol()};

    if (describe)
        std::cout << describe << " : " << nrows << " x " << ncols << " Matrix:\n";

    ncols--;
    for (int row{0}, col(0); row < nrows; ++row)
    {
        std::cout << "  row[" << row << "]: ";
        for (col = 0; col < ncols; ++col)
            std::cout << matG(row, col) << ",  ";
        std::cout << matG(row, col) << ";\n";
    }
    std::cout << "\n";
}

// Matrix multiplication O(n^3) naive implementation
void matrix_mult(matrix<double> &resMat,
		const matrix<double> &matA, const matrix<double> &matB)
{
	const int nrowsA{matA.nrow()}, ncolsA{matA.ncol()},
					nrowsB{matB.nrow()}, ncolsB{matB.ncol()};
	if (ncolsA != nrowsB)
	{
		std::cout << "matrix_mult() error: dimension mismatch.\n";
		return;
	}
	resMat.init(0, nrowsA, ncolsB);
	int i{0}, j{0}, k{0};
	for (i = 0; i < nrowsA; ++i)
		for (j = 0; j < ncolsB; ++j)
			for (k = 0; k < ncolsA; ++k)
				resMat(i, j) += matA(i, k) * matB(k, j);
}

// LU decomposition-based matrix inversion [*3][*4]
void inv(matrix<double> &matLU, const matrix<double> &matG)
{
	const int nrows{matG.nrow()}, ncols{matG.ncol()};
	if (nrows != ncols) {
		std::cout << "Error:  inv(matrix) is not square.\n";
		return;
	}
	if (matG(0) == 0.0)
	{
		std::cout << "Warning:  inv(matrix) is singular.\n";
		return;
	}

	const int nSize{nrows};
	int i{0}, j{0}, k{0};

	// ------------- Step 1: row permutation (swap diagonal zeros) -----------
	std::vector<int> permuteLU; // Permute vector
	for (i = 0; i < nSize; ++i)
		permuteLU.push_back(i); // Push back row index

	matLU = matrix<double>(matG); // Simply duplicate matrix
//	showMatrix(matLU, "inv matLU before");

	// ----------- Step 2: LU decomposition (save both L & U in matLU) -------
	for (i = 1; i < nSize; ++i)
		matLU(i, 0) /= matLU(0); // Initialize first column of L matrix

	for (i = 1; i < nSize; ++i)
	{
		for (j = i; j < nSize; ++j)
			for (k = 0; k < i; ++k)
				matLU(i, j) -= matLU(i, k) * matLU(k, j); // Calculate U matrix
		if (matLU(i, i) == 0.0)
		{
			std::cout << "inv warning: singular matrix.\n";
			return;
		}
		for (k = i + 1; k < nSize; ++k)
		{
			for (j = 0; j < i; ++j)
				matLU(k, i) -= matLU(k, j) * matLU(j, i); // Calculate L matrix
			matLU(k, i) /= matLU(i, i);
		}
	}
//	showMatrix(matLU, "inv matLU after");

	// ----- Step 3: L & U inversion (save both L^-1 & U^-1 in matLU_inv) -----
	matrix<double> matLU_inv;  matLU_inv.init(0, nSize, nSize);

	// matL inverse & matU inverse
	for (i = 0; i < nSize; ++i)
	{
		// L matrix inverse, omit diagonal ones
		matLU_inv(i, i) = 1.0;
		for (k = i + 1; k < nSize; ++k)
			for (j = i; j <= k - 1; ++j)
				matLU_inv(k, i) -= matLU(k, j) * matLU_inv(j, i);
		// U matrix inverse
		matLU_inv(i, i) = 1.0 / matLU(i, i);
		for (k = i; k > 0; --k)
		{
			for (j = k; j <= i; ++j)
				matLU_inv(k - 1, i) -= matLU(k - 1, j) * matLU_inv(j, i);
			matLU_inv(k - 1, i) /= matLU(k - 1, k - 1);
		}
	}
//	showMatrix(matLU_inv, "inv matLU_inv");

	// ******************** Step 4: Calculate G^-1 = U^-1 * L^-1 ********************
	// Lower part product
	for (i = 1; i < nSize; ++i)
		for (j = 0; j < i; ++j)
		{
			const int jp{permuteLU[j]}; // Permute column back
			matLU(i, jp) = 0.0;
			for (k = i; k < nSize; ++k)
				matLU(i, jp) += matLU_inv(i, k) * matLU_inv(k, j);
		}
	// Upper part product
	for (i = 0; i < nSize; ++i)
		for (j = i; j < nSize; ++j)
		{
			const int jp{permuteLU[j]}; // Permute column back
			matLU(i, jp) = matLU_inv(i, j);
			for (k = j + 1; k < nSize; ++k)
				matLU(i, jp) += matLU_inv(i, k) * matLU_inv(k, j);
		}
//	showMatrix(matLU, "inv matLU final");
}

// Moore-Penrose pseudoinversion (same as pinv(G) in MATLAB) [*1]
void mpinvert(matrix<double> &mpi, const matrix<double> &matG,
			const double tolerance)
{
	bool useTranspose{false};
	const int nrows{matG.nrow()}, ncols{matG.ncol()};
	int nSize{ncols};

	matrix<double> matA, matGt; matG.t(matGt);
	if (nrows < nSize)
	{
		useTranspose = true;
		nSize = nrows;
		matrix_mult(matA, matG, matGt); // A = G * G'
	}
	else matrix_mult(matA, matGt, matG); // A = G' * G
//	showMatrix(matA, "mpinvert matA = GGt");

	// Full rank Cholesky decomposition of A
	int i{0}, j{0}, k{0};

	double tol{abs(matA(0))};
	for (i = 0; i < nSize; ++i)
		if (matA(i, i) > 0)
		{
			const double temp{matA(i, i)};
			if (temp < tol)
				tol = temp;
		}
	tol *= tolerance;

	matrix<double> matL; matL.init(0, nSize, nSize);
	int rankA{0};
	for (k = 0; k < nSize; ++k)
	{
		for (i = k; i < nSize; ++i)
		{
			matL(i, rankA) = matA(i, k);
			for (j = 0; j < rankA; ++j)
				matL(i, rankA) -= matL(i, j) * matL(k, j);
		}
		if (matL(k, rankA) > tol)
		{
			matL(k, rankA) = sqrt(matL(k, rankA));
			if (k < nSize)
				for (j = k + 1; j < nSize; ++j)
					matL(j, rankA) /= matL(k, rankA);
			++rankA;
		}
	}
//	showMatrix(matL, "mpinvert matL = Cholesky decomposition of A");

	if (rankA == 0) {
		mpi = matrix<double>(matGt); // All-zero matrix's transpose
		return;
	}

	// Slice L = L(:, 0:r);
	matL.keepCols(rankA);

	// Generalized inverse
	matrix<double> matLt;	matL.t(matLt);
	matrix<double> matLtL;	matrix_mult(matLtL, matLt, matL);
	matrix<double> matM;	inv(matM, matLtL);				// M = inv(L' * L)	
	matrix<double> matLM;	matrix_mult(matLM, matL, matM);	// L*M
//	showMatrix(matLM, "mpinvert matLM = L*(M=inv(LtL))");
	matrix<double> matLMM;	matrix_mult(matLMM, matLM, matM); // L*M*M
//	showMatrix(matLMM, "mpinvert matLMM");
	matrix_mult(matA, matLMM, matLt); // A = L * M * M * L'
//	showMatrix(matA, "mpinvert matA = LMMLt");

	if (useTranspose)
		matrix_mult(mpi, matGt, matA);	// pinv(G) = G' * (L * M * M * L')
	else matrix_mult(mpi, matA, matGt);	// pinv(G) = (L * M * M * L') * G'
}
