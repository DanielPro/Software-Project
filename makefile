FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: cluster
clean:
	rm -rf *.o cluster

cluster: cluster.o division.o eigenpair.o error.o graph.o spmat.o vectormath.o
	gcc cluster.o division.o eigenpair.o error.o graph.o spmat.o vectormath.o -o cluster $(LIBS)
cluster.o: cluster.c division.h eigenpair.h error.h graph.h spmat.h vectormath.h
	gcc $(FLAGS) -c cluster.c

division.o: division.h division.c eigenpair.h graph.h spmat.h vectormath.h
	gcc $(FLAGS) -c division.c
eigenpair.o: eigenpair.h eigenpair.c error.h vectormath.h
	gcc $(FLAGS) -c eigenpair.c
error.o: error.h error.c
	gcc $(FLAGS) -c error.c
graph.o: graph.h graph.c
	gcc $(FLAGS) -c graph.c
spmat.o: spmat.h spmat.c error.h
	gcc $(FLAGS) -c spmat.c
vectormath.o: vectormath.h vectormath.c error.h
	gcc $(FLAGS) -c vectormath.c