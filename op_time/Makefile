TARGET=op_time
CFLAGS+=-I/home/tanghoujun/build/include/ -O2
LDFLAGS+=-L/home/tanghoujun/build/lib -O2
CC=mpicc
all: $(TARGET)

.c.o:
	mpicc $(CFLAGS) -c $< -o $@

$(TARGET): $(TARGET).o
	mpicc $(TARGET).o  $(LDFLAGS) -lgsl -lgslcblas -lm -o $@

clean:
	rm -f $(TARGET) *.o
