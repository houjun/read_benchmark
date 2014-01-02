CC=/home/tanghoujun/build/mpich/bin/mpicc
#CFLAGS=-g
#CFLAGS=-O2
RM=rm -f

SRCS=$(wildcard *.c)
OUTS=$(patsubst %.c,%,$(SRCS))

all:$(OUTS)

$(OUTS):%:%.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(OUTS)


