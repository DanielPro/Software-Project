/*
 * division.h
 * Module for the division operation.
 * given the modularity matrix *mat and its indices group *group
 * will generate the optimal modularity inside *groups graph.
 */

#ifndef DIVISION_H_
#define DIVISION_H_

/* Recursively divides an indices group according to the modularity matrix, keeping the result in *groups.
 * Frees all resources as well. -Final calculation- */
void recursive_division(spmat *mat, Graph *groups, int *group);

#endif
