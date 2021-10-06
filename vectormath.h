/*
 * vectormath.h
 * Module for vector arithmetics.
 */

#ifndef VECTORMATH_H_
#define VECTORMATH_H_

/* Calculates 2 vectors' dot product */
double dot_product(double *vec1, double *vec2, int size);

/* Normalizes the vector in place */
void normalize(double *vec, int size);

/* Returns 1 when no coordinate between the 2 vectors differ by more then 0.00001 */
int check(double *vec1, double *vec2, int size);

/* Multiplies a given vector by -1 in place */
void vector_minus(double *vec, int size);

#endif
