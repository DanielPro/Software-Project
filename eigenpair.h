/*
 * eigenpair.h
 * Module for eigen pair computation.
 */

#ifndef EIGENPAIR_H_
#define EIGENPAIR_H_

#include "spmat.h"

/*Structure for holding the eigen pair values*/
typedef struct _eigenpair {

	double eigenvalue, *eigenvector;

} EigenPair;

/*Creates, computes and returns a NEW eigen pair given an spmat matrix*/
EigenPair* compute_eigenpair(spmat *matrix);

/*Frees all resources*/
void free_eigenpair(EigenPair *A);

#endif
