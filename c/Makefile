CC=gcc
CFLAGS=-O3 -flto -march=native
LDFLAGS=-lm
OBJS=main.o bvp.o grid.o
TARGET=star

all: $(TARGET)

clean:
	rm -f *.o
	rm -f $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

main.o: main.c grid.h bvp.h
bvp.o: bvp.c bvp.h grid.h
grid.o: grid.c grid.h