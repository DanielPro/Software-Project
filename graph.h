/*
 * graph.h
 * Module for representing the indices division into groups.
 */

#ifndef GRAPH_H_
#define GRAPH_H_

/* Representing a division of indices into groups. Each group will have it's own section in *groups
 * based on the values inside groupSizes */
typedef struct _graph {

	/* Total size */
	int size;

	/* number of groups */
	int groupCount;

	int *groups;

	int *groupSizes;

} Graph;

/*Creates and returns a new empty graph of a given size*/
Graph* allocate_graph(int size);

/*Frees all resources*/
void free_graph(Graph *A);

/* Adds new group to the graph given their indices and the group size*/
void add_group(Graph *graph, int *group, int size);

#endif
