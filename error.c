/*
 * error.c
 */

#include "error.h"
#include <stdlib.h>
#include <stdio.h>

void process_err(char *error_text) {

	printf("Error! Reason: %s", error_text);
	exit(-1);
}
