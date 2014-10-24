all:
#CC read2D.c -o read2D -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich
	mpicc read2D.c space_filling/space_filling.c -o read2D -g

debug:
#CC read2D.c -o read2D -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich -g
	mpicc read2D.c space_filling/space_filling.c -o read2D -g

clean:
	rm -f *.o read2D
