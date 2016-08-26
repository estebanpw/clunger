CC=gcc
CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE
LIBS=-lm
DEPS=common.h structs.h NWscore2rows.h

all: buildcluster
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

buildcluster: buildcluster.o common.o NWscore2rows.o

.PHONY: clean
clean:
	rm -f *.o buildcluster