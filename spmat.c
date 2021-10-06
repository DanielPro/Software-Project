/*
 * spmat.c
 */

#include <stdlib.h>
#include <math.h>
#include "error.h"
#include "spmat.h"

#define IS_POSITIVE(X) ((X) > 0.00001)

/* Allocates a NEW empty sparse matrix given its size n and nnz non-zero elements */
spmat* allocate_spmat(int size, int nnz) {

	spmat *matrix = malloc(sizeof(spmat));

	matrix->n = size;
	matrix->columnIndices = calloc(nnz, sizeof(int));
	matrix->rowIndices = calloc(size + 1, sizeof(int));
	matrix->values = calloc(nnz, sizeof(double));

	return matrix;
}

void free_spmat(spmat *A) {

	free(A->columnIndices);
	free(A->rowIndices);
	free(A->values);
	free(A);
}

/* Adds row i to the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1)
 * Uses special diagonal format built by buildModularSpmat/BuildHatModularSpmat */
void add_row(spmat *mat, void *row, int i) {

	double *rows = row;
	int j, ind = 0, start = mat->rowIndices[i], counter = start;

	for (j = start; j < start + mat->n - i; ++j) {

		if (rows[ind] != 0) {
			mat->values[counter] = rows[ind];
			mat->columnIndices[counter] = ind + i;
			counter++;
		}
		ind++;
	}
	mat->rowIndices[i + 1] = counter;
}

/* Multiplys a given matrix and a vector */
double* mult(spmat *mat, double *v) {

	int i, j = 0, nnz, size = mat->n;
	double *result = calloc(size, sizeof(double));

	for (i = 0; i < size; ++i) {

		nnz = mat->rowIndices[i + 1] - mat->rowIndices[i];
		j = mat->rowIndices[i];

		while (nnz > 0) {
			result[i] += mat->values[j] * v[mat->columnIndices[j]];

			if (i != mat->columnIndices[j])
				result[mat->columnIndices[j]] += mat->values[j] * v[i];

			j++;
			nnz--;
		}
	}

	return result;
}

/* Returns a value stored in the sparse matrix by row and column efficiently.
 * Implements binary search. */
double get_value(spmat *mat, int row, int column) {

	int temp, lowBound, highBound, mid;

	/* Check for valid values */
	if (row >= mat->n || column >= mat->n || row < 0 || column < 0)
		process_err(
				"Attempted to read a value from a matrix with an out of bounds index.");

	/* We keep half the matrix. Switch row & column*/
	if (column < row) {
		temp = column;
		column = row;
		row = temp;
	}

	lowBound = mat->rowIndices[row];
	highBound = mat->rowIndices[row + 1];

	/* Start binary search */
	if (lowBound == highBound)
		return 0;

	while (lowBound <= highBound) {
		mid = (lowBound + highBound) / 2;

		if (column < mat->columnIndices[mid])
			highBound = mid - 1;
		else if (column > mat->columnIndices[mid])
			lowBound = mid + 1;
		else
			return mat->values[mid];
	}

	return 0;
}

/* Used by buildHatModular. Returns the sum of a row according to a division vector */
double get_partial_row_sum(spmat *mat, double *division, int row) {

	double sum = 0;
	int idx = mat->rowIndices[row], rowNnz;

	/* Check for valid values */
	if (row >= mat->n || row < 0)
		process_err(
				"Attempted to read a value from a matrix with an out of bounds index.");

	/* Sum values for column < row using get_value */
	for (rowNnz = 0; rowNnz < row; ++rowNnz)
		if (division[rowNnz] == 1)
			sum += get_value(mat, row, rowNnz);

	/* Sum the rest of the values directly. column >= row */
	for (rowNnz = mat->rowIndices[row + 1] - idx; rowNnz > 0; --rowNnz) {
		if (division[mat->columnIndices[idx]] == 1)
			sum += mat->values[idx];
		idx++;
	}

	return sum;
}

/* Returns the 1-norm of a matrix. Used for matrix shifting operation. */
double get_l_norm(spmat *mat) {

	double *colSums, lNorm = 0;
	int size = mat->n, i, j, nnz, start;

	colSums = calloc(size, sizeof(double));

	/* Main loop for summing each column, keeping the result in colSums */
	for (i = 0; i < size; i++) {

		start = mat->rowIndices[i];
		nnz = mat->rowIndices[i + 1] - mat->rowIndices[i];

		for (j = start; j < nnz + start; j++) {

			colSums[mat->columnIndices[j]] += fabs(mat->values[j]);

			if (i != mat->columnIndices[j])
				colSums[i] += fabs(mat->values[j]);
		}
	}

	/* Checking for maximun value */
	for (i = 0; i < size; i++)
		if (colSums[i] > lNorm)
			lNorm = colSums[i];

	free(colSums);
	return lNorm;

}

/* Builds & returns a NEW shifted matrix for the eigen pair computation. */
spmat* shift_matrix(spmat *mat) {

	spmat *shiftedMat;
	double lNorm = get_l_norm(mat), **fauxMatrix;
	int size = mat->n, i, j, nnz = 0;

	fauxMatrix = calloc(size, sizeof(double*));

	/* Building the fauxMatrix, as well as counting the nnz required for the new shifted matrix allocation */
	for (i = 0; i < size; ++i) {
		fauxMatrix[i] = calloc(size - i, sizeof(double));

		for (j = i; j < size; ++j) {

			if (i == j) /* Adding the 1-norm for the matrix diagonal line */
				fauxMatrix[i][j - i] = get_value(mat, i, j) + lNorm;
			else
				/* Copy the rest of the values as is */
				fauxMatrix[i][j - i] = get_value(mat, i, j);
			/* Counting non-zero elements */
			if (fauxMatrix[i][j - i] != 0)
				nnz++;
		}
	}

	shiftedMat = allocate_spmat(size, nnz);

	/* Adding all rows from fauxMatrix after allocation */
	for (i = 0; i < size; ++i) {
		add_row(shiftedMat, fauxMatrix[i], i);
		free(fauxMatrix[i]);
	}

	free(fauxMatrix);
	return shiftedMat;
}
/* Building the first modularity matrix from the graph adjacency matrix.
 * Called only once at the beginning. */
spmat* build_modular_spmat(int **mat, int size, int matNnz) {

	int i, j, k_i, nnzTotal = 0;
	double B_i_j, M = (double) matNnz, **fauxMatrix = calloc(size,
			sizeof(double*));
	spmat *modularSpmat;

	/* Check for a valid value */
	if (M == 0)
		process_err(
				"Matrix is all zeroes. This will lead to attempting division by zero when building modular matrix.");

	/* Main loop for fauxMatrix creation */
	for (i = 0; i < size; ++i) {
		fauxMatrix[i] = calloc(size - i, sizeof(double));
		k_i = mat[i][size]; /* Index degree kept in last extra column */

		/* Calculating values for each row */
		for (j = i; j < size; ++j) {
			B_i_j = (double) mat[i][j] - (double) (k_i * mat[j][size]) / M;

			/* Chceck for zero's and count nnz for allocation */
			if (IS_POSITIVE(fabs(B_i_j))) {
				fauxMatrix[i][j - i] = B_i_j;
				nnzTotal++;
			}
		}
		free(mat[i]);
	}
	free(mat);

	modularSpmat = allocate_spmat(size, nnzTotal);

	/* Adding all rows from fauxMatrix after allocation */
	for (i = 0; i < size; ++i) {
		add_row(modularSpmat, fauxMatrix[i], i);
		free(fauxMatrix[i]);
	}

	free(fauxMatrix);
	return modularSpmat;
}

spmat* build_hat_modular_spmat(spmat *modularMat, double *division) {

	int i, j, nnz = 0, newSize = 0, row = 0, column = 1, size = modularMat->n;
	double f_i_g, result, **modularMatrixSimple;
	spmat *hatModular;

	/* Counting new size from division */
	for (i = 0; i < size; ++i) {
		if (division[i] == 1) {
			newSize++;
		}
	}

	modularMatrixSimple = calloc(newSize, sizeof(double*));

	/* Main loop for building the simple matrix */
	for (i = 0; i < size; ++i) {

		if (division[i] == 1) {
			f_i_g = get_partial_row_sum(modularMat, division, i);
			modularMatrixSimple[row] = calloc(newSize - row, sizeof(double));
			result = get_value(modularMat, i, i) - f_i_g; /* Calculation for diagonal line values*/

			/* Check for zero's and count nnz for allocation */
			if (IS_POSITIVE(fabs(result))) {
				modularMatrixSimple[row][0] = result;
				nnz++;
			}

			/* Loop for non diagonal line values*/
			for (j = i + 1; j < size; ++j) {
				if (division[j] == 1) {
					result = get_value(modularMat, i, j);

					/* Check for zero's and count nnz for allocation */
					if (IS_POSITIVE(fabs(result))) {
						modularMatrixSimple[row][column] = result;
						nnz++;
					}
					column++;
				}
			}
			column = 1;
			row++;
		}
	}

	hatModular = allocate_spmat(newSize, nnz);

	/* Adding all rows from modularMatrixSimple after allocation */
	for (i = 0; i < newSize; ++i) {
		add_row(hatModular, modularMatrixSimple[i], i);
		free(modularMatrixSimple[i]);
	}

	free(modularMatrixSimple);
	return hatModular;
}
