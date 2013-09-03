CC = g++
RM = rm
CFLAGS = -g -Wall
TARGET = bf
LIBS = -L/usr/lib -lboost_program_options

$(TARGET): main.o Shapes.o Double.o BestFit.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o Shapes.o Double.o BestFit.o $(LIBS)
main.o: main.cpp BestFit.h
	$(CC) $(CFLAGS) -c main.cpp
Shapes.o: Shapes.cpp BestFit.h
	$(CC) $(CFLAGS) -c Shapes.cpp
Double.o: Double.cpp BestFit.h
	$(CC) $(CFLAGS) -c Double.cpp
BestFit.o: BestFit.cpp BestFit.h
	$(CC) $(CFLAGS) -c BestFit.cpp

clean:
	$(RM) main.o Shapes.o BestFit.o Double.o
