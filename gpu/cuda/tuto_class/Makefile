CC = nvcc

SRC1 = addVecs.cu
SRC2 = rev-array.cu
SRC3 = inplace_rev-array.cu
all:
	$(CC) $(SRC1) -o addVecs
	$(CC) $(SRC2) -o rev-array
	$(CC) $(SRC3) -o inplace_rev-array
clean:
	rm -rf *.o addVecs rev-array inplace_rev-array
