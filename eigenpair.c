/*
 * eigenpair.c
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "error.h"
#include "vectormath.h"
#include "eigenpair.h"

EigenPair* compute_eigenpair(spmat *matrix) {

	int i, size = matrix->n;
	double *eigenVec, *tempVec, product;
	spmat *shiftedSpmat;
	EigenPair *eigenPair = malloc(sizeof(EigenPair));

	shiftedSpmat = shift_matrix(matrix);

	eigenVec = calloc(size, sizeof(double));

	for (i = 0; i < size; i++)
		eigenVec[i] = rand();

	tempVec = eigenVec;
	eigenVec = mult(shiftedSpmat, eigenVec);
	normalize(eigenVec, size);

	/* Implementation of the power iteration */
	while (!check(eigenVec, tempVec, size)) {
		free(tempVec);
		tempVec = eigenVec;
		eigenVec = mult(shiftedSpmat, tempVec);
		normalize(eigenVec, size);
	}
	free(tempVec);

	product = dot_product(eigenVec, eigenVec, size);
	if (product == 0)
		process_err(
				"Eigenvector found is the zero vector! This will lead to division by zero.");

	tempVec = mult(shiftedSpmat, eigenVec);
	eigenPair->eigenvalue = dot_product(eigenVec, tempVec, size) / product
			- get_l_norm(matrix);
	eigenPair->eigenvector = eigenVec;

	free(tempVec);
	free_spmat(shiftedSpmat);

	return eigenPair;
}

void free_eigenpair(EigenPair *A) {

	free(A->eigenvector);
	free(A);
}
