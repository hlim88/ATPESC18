CC = gcc-6

SRC1 = mandel_prob.c
SRC2 = mandel_sol.c

all:
	$(CC) -fopenmp $(SRC1) -o mandel_prob
	$(CC) -fopenmp $(SRC2) -o mandel_sol
clean:
	rm -rf *.o mandel_prob mandel_sol 
