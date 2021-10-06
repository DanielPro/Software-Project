/*
 * vectormath.c
 */

#include <math.h>
#include "error.h"
#include "vectormath.h"

#define IS_POSITIVE(X) ((X) > 0.00001)

double dot_product(double *vec1, double *vec2, int size) {

	int i;
	double product = 0;
	for (i = 0; i < size; ++i) {
		product += vec1[i] * vec2[i];
	}
	return product;
}

void normalize(double *vec, int size) {

	double norm = 0;
	int i;

	for (i = 0; i < size; ++i) {
		norm += vec[i] * vec[i];
	}
	norm = sqrt(norm);

	if (norm == 0)
		process_err("Attempted division by 0 when normalizing vector.");

	for (i = 0; i < size; ++i) {
		vec[i] /= norm;
	}
}

int check(double *vec1, double *vec2, int size) {

	int i;

	for (i = 0; i < size; ++i)
		if (IS_POSITIVE(fabs(vec1[i] - vec2[i])))
			return 0;

	return 1;
}

void vector_minus(double *vec, int size) {

	int i;

	for (i = 0; i < size; i++)
		vec[i] *= -1;

}
