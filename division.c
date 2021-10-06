/*
 * division.c
 */

#include <stdlib.h>
#include "eigenpair.h"
#include "graph.h"
#include "vectormath.h"
#include "spmat.h"
#include "division.h"

/* Responsible for dividing into two groups using the matrix eigen pair.
 * Returns a +-1 vector represnting the division.
 * If no division is of bigger modularity, returns "0" vector. */
double* divide(spmat *mat);

/* Responsible for improving the vector from the divide function.
 * Using it to find the optimal division into two groups. */
void improve(double *division, spmat *mat);

void recursive_division(spmat *mat, Graph *groups, int *group) {

	int *firstGroup, *secondGroup, plusCount = 0, minusCount = 0, i, size =
			mat->n;
	double *division = divide(mat);
	firstGroup = calloc(size, sizeof(int));
	secondGroup = calloc(size, sizeof(int));

	/* Splitting to groups by division. For "0" vector, all indices go to second group. (no division) */
	for (i = 0; i < size; ++i) {
		if (division[i] == 1) {
			firstGroup[plusCount] = group[i];
			plusCount++;
		} else {
			secondGroup[minusCount] = group[i];
			minusCount++;
		}
	}

	/* First case: ont of the groups is of size zero or both of size 1. Adding both groups*/
	if (plusCount == 0 || minusCount == 0
			|| (plusCount == 1 && minusCount == 1)) {
		add_group(groups, firstGroup, plusCount);
		add_group(groups, secondGroup, minusCount);
		free(firstGroup);
		free(secondGroup);
	} else {
		/* Second case: first group of size 1, second is larger. Add first, recursive call on second. */
		if (plusCount == 1) {
			add_group(groups, firstGroup, plusCount);
			free(firstGroup);
			vector_minus(division, size);
			recursive_division(build_hat_modular_spmat(mat, division), groups,
					secondGroup);
		} else {
			/* Third case: second group of size 1, first is larger. Add second, recursive call on first. */
			if (minusCount == 1) {
				add_group(groups, secondGroup, minusCount);
				free(secondGroup);
				recursive_division(build_hat_modular_spmat(mat, division),
						groups, firstGroup);
			} else { /*All other cases. Both groups of size bigger then 1 */
				recursive_division(build_hat_modular_spmat(mat, division),
						groups, firstGroup);
				vector_minus(division, size);
				recursive_division(build_hat_modular_spmat(mat, division),
						groups, secondGroup);
			}
		}
	}

	free(group);
	free_spmat(mat);
	free(division);
}

double* divide(spmat *mat) {

	int i, size = mat->n;
	double *tempVec, *division = calloc(size, sizeof(double));
	EigenPair *pair;

	pair = compute_eigenpair(mat);

	/* First condition check */
	if (pair->eigenvalue <= 0)
		return division;

	/* Vector build */
	for (i = 0; i < size; ++i) {
		if (pair->eigenvector[i] > 0)
			division[i] = 1;
		else
			division[i] = -1;
	}
	free_eigenpair(pair);

	tempVec = mult(mat, division);

	/* Second condition check */
	if (dot_product(tempVec, division, size) <= 0) {
		free(division);
		division = calloc(size, sizeof(double));
		free(tempVec);
		return division;
	}

	free(tempVec);
	improve(division, mat); /* Check for improvement */
	return division;
}

void improve(double *division, spmat *mat) {

	int i, j, max, *indices, *moved, size = mat->n;
	double Qs, *delta, *improve, deltaQ, *tempVec;

	delta = calloc(size, sizeof(double));
	improve = calloc(size, sizeof(double));
	indices = calloc(size, sizeof(int));

	do {
		moved = calloc(size, sizeof(int));

		for (i = 0; i < size; ++i) {

			tempVec = mult(mat, division);
			Qs = dot_product(tempVec, division, size);
			free(tempVec);

			/* Calculating modularity change for each indice */
			for (j = 0; j < size; ++j) {
				if (moved[j] != 1) {
					division[j] *= -1;
					tempVec = mult(mat, division);
					delta[j] = dot_product(tempVec, division, size) - Qs;
					free(tempVec);
					division[j] *= -1;
					max = j;
				}
			}

			/* Finding maximun change */
			for (j = 0; j < size; ++j)
				if (moved[j] != 1)
					if (delta[j] > delta[max])
						max = j;

			/* Moving indice accordingly */
			division[max] *= -1;
			moved[max] = 1;
			indices[i] = max;

			/* Saving current change */
			if (i == 0)
				improve[i] = delta[max];
			else
				improve[i] = improve[i - 1] + delta[max];

		}

		max = 0;

		/* Finding point of maximun improvement */
		for (i = 1; i < size; ++i)
			if (improve[i] > improve[max])
				max = i;

		/* Moving all indices to this point */
		for (i = size - 1; i > max; --i)
			division[indices[i]] *= -1;

		if (max == size - 1) /* Moved all indices -> same division, no change */
			deltaQ = 0;
		else
			deltaQ = improve[max];

		free(moved);

	} while (deltaQ > 0); /* Continue as long as there is an improvement */

	free(delta);
	free(improve);
	free(indices);
}
