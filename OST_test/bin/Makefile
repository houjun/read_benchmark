all:
	CC procNode.c -o procNode -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich
	CC OST_test.c -o OST_test -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich
	CC new_OST_test.c -o new_OST_test -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich
	CC new_OST_multitest.c -o new_OST_multitest -I$(MPICH_ROOT)/include -L$(MPICH_ROOT)/lib -lmpich


clean:
	rm -f *.o procNode OST_test new_OST_test new_OST_multitest
