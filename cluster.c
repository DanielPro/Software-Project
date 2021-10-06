/*
 * cluster.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"
#include "eigenpair.h"
#include "graph.h"
#include "vectormath.h"
#include "division.h"
#include "spmat.h"

int main(int argc, char *argv[]) {

	FILE *input, *output;
	int size, i, vertex, nnzTotal = 0, nnz = 0, **graphMat, *group;
	spmat *modularMat;
	Graph *groups;

	argc > 2 ? (void) 1 : process_err("Too few program arguments.");
	input = fopen(argv[1], "r");
	input != NULL ? (void) 1 : process_err("Couldn't open input file.");

	/* getting matrix size*/
	fread(&size, sizeof(int), 1, input);
	graphMat = calloc(size, sizeof(int*));

	/* initializing graph matrix, last column for vertex's degree */
	for (i = 0; i < size; i++) {

		fread(&nnz, sizeof(int), 1, input);
		nnzTotal += nnz;
		graphMat[i] = calloc(size + 1, sizeof(int));
		graphMat[i][size] = nnz;

		while (nnz > 0) {
			fread(&vertex, sizeof(int), 1, input);
			graphMat[i][vertex] = 1;
			nnz--;
		}
	}

	fclose(input);

	/* Frees GraphMat during calculation. Only use of graphMat. */
	modularMat = build_modular_spmat(graphMat, size, nnzTotal);

	groups = allocate_graph(size);
	group = calloc(size, sizeof(int));

	for (i = 0; i < size; ++i)
		group[i] = i;

	/* Frees modularMat & group after final use */
	recursive_division(modularMat, groups, group);

	/* Initializing & writing to output file */
	output = fopen("groups.out", "w");
	fwrite(&groups->groupCount, sizeof(int), 1, output) == 1 ?
			(void) 1 : process_err("Writing to output file failed.");

	for (i = 1; i <= groups->groupCount; ++i) {
		nnzTotal = groups->groupSizes[i] - groups->groupSizes[i - 1];
		fwrite(&nnzTotal, sizeof(int), 1, output) == 1 ?
				(void) 1 : process_err("Writing to output file failed.");

		for (nnz = groups->groupSizes[i - 1]; nnz < groups->groupSizes[i];
				++nnz)
			fwrite(&groups->groups[nnz], sizeof(int), 1, output) == 1 ?
					(void) 1 : process_err("Writing to output file failed.");
	}

	free_graph(groups);

	printf("----Finished----");
	return 0;
}
