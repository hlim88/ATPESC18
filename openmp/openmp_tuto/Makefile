#CC = gcc
CC = gcc-6 # For Mac

SRC1 = hello_world.c
SRC2 = omp_task_example.c
all:
	$(CC) -fopenmp $(SRC1) -o hello_world
	$(CC) -fopenmp $(SRC2) -o omp_task_example
clean:
	rm -rf *.o hello_world omp_task_example 
