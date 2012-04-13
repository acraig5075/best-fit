CC = g++
DEP = /usr/lib/libboost_program_options.so

CFLAGS = -g -Wall

bf: main.o Shapes.o BestFit.o
	$(CC) $(CFLAGS) -o bf main.o Shapes.o BestFit.o $(DEP)
main.o: main.cpp Shapes.h
	$(CC) $(CFLAGS) -c main.cpp
Shapes.o: Shapes.cpp Shapes.h
	$(CC) $(CFLAGS) -c Shapes.cpp
BestFit.o: BestFit.cpp BestFit.h
	$(CC) $(CFLAGS) -c BestFit.cpp

clean:
	rm main.o Shapes.o BestFit.o
