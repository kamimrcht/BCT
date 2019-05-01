CC=g++
CFLAGS= -Wall -Ofast -std=c++11  -flto -pipe -funit-at-a-time -fopenmp -lz -Isparsepp -flto
LDFLAGS=-flto -lpthread -fopenmp -lz  -Isparsepp  -flto



EXEC=clean_homopoly recover_tail

all: $(EXEC)



recover_tail:   recover_tail.o
	$(CC) -o $@ $^ $(LDFLAGS)

recover_tail.o: RecoverTail.cpp zstr.hpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean_homopoly:   clean_homopoly.o
	$(CC) -o $@ $^ $(LDFLAGS)

clean_homopoly.o: clean_homopoly.cpp zstr.hpp
	$(CC) -o $@ -c $< $(CFLAGS)
