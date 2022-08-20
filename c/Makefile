CC=gcc
CFLAGS=
LDFLAGS=-lm
OBJS=main.o multigrid.o
TARGET=star

all: $(TARGET)

clean:
	rm -f *.o
	rm -f $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

main.o: multigrid.h main.c
multigrid.o: multigrid.h multigrid.c