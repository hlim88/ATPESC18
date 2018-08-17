CC = gcc-6

SRC1 = pi_serial.c
SRC2 = pi_parallel.c
SRC3 = pi_parallel_pad.c
SRC4 = pi_parallel_fs.c
SRC5 = pi_parallel_loop.c
SRC6 = pi_parallel_task.c

all:
	$(CC) -fopenmp $(SRC1) -o pi_serial
	$(CC) -fopenmp $(SRC2) -o pi_parallel
	$(CC) -fopenmp $(SRC3) -o pi_parallel_pad
	$(CC) -fopenmp $(SRC4) -o pi_parallel_fs
	$(CC) -fopenmp $(SRC5) -o pi_parallel_loop
	$(CC) -fopenmp $(SRC6) -o pi_parallel_task
clean:
	rm -rf *.o pi_serial pi_parallel pi_parallel_pad pi_parallel_fs pi_parallel_loop pi_parallel_task
