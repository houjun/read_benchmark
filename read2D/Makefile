all:
	CC read2D.c -o read2D -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich

debug:
	CC read2D.c -o read2D -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich -g

clean:
	rm -f *.o read2D
