#CC = gcc
CC = gcc-6 # For mac gnu

SRC1 = linked_serial.c
SRC2 = linked_parallel.c
SRC3 = linked_parallel_task.c

all:
	$(CC) -fopenmp $(SRC1) -o linked_serial
	$(CC) -fopenmp $(SRC2) -o linked_parallel
	$(CC) -fopenmp $(SRC3) -o linked_parallel_task
clean:
	rm -rf *.o linked_serial linked_parallel linked_parallel_task
