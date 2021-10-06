/*
 * spmat.h
 * Module for represnting a symethric matrix as a sparse matrix. (keeps half of it)
 * Used for representing all of our matrixes in an efficient way, since they are all symethric.
 */

#ifndef SPMAT_H
#define SPMAT_H

typedef struct _spmat {
	/* Matrix size */
	int n;

	int *rowIndices;

	int *columnIndices;

	double *values;

} spmat;

/* Frees all resources */
void free_spmat(spmat *A);

/* Allocates & builds the modularity sparse matrix given the adjacency matrix of the graph.
 * Only being called once at the beginning of the program in cluster.
 * Frees all allocated resources (int** mat) */
spmat* build_modular_spmat(int **mat, int size, int matNnz);

/* Allocates & builds the "hat" modularity matrix given a modularity matrix and a +-1 vector of a division */
spmat* build_hat_modular_spmat(spmat *modularMat, double *division);

/* Multiplyes a given matrix with a vector. returns a NEW vector containing the result */
double* mult(spmat *mat, double *v);

/* Returns the 1-norm of a matrix required for the matrix shifting operation */
double get_l_norm(spmat *mat);

/* Returns a NEW shifted matrix */
spmat* shift_matrix(spmat *mat);

#endif
