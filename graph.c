/*
 * graph.c
 */

#include <stdlib.h>
#include "graph.h"

Graph* allocate_graph(int size) {

	Graph *graph = malloc(sizeof(Graph));

	graph->size = size;
	graph->groupCount = 0;
	graph->groups = calloc(size, sizeof(int));
	graph->groupSizes = calloc(size + 1, sizeof(int));

	return graph;
}

void free_graph(Graph *A) {

	free(A->groups);
	free(A->groupSizes);
	free(A);
}

void add_group(Graph *graph, int *group, int size) {

	int i, prev = graph->groupSizes[graph->groupCount];

	/* Checks the group is not empty and adds it in the given order. (which is always increasing order from division.c) */
	if (size > 0) {
		graph->groupCount += 1;
		graph->groupSizes[graph->groupCount] = prev + size;

		for (i = prev; i < prev + size; ++i)
			graph->groups[i] = group[i - prev];
	}
}
