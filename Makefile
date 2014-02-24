# Makefile for best-fit
# Alasdair Craig
#
# Type "make" to build best-fit as an executable.
#
# Type "make library" to build best-fit as a static library.
#
# Type "make clean" to delete all object files.
#
# Type "make distclean" to delete all object files and the executable, and the library file.

CC = g++
RM = rm -f
AR = ar
RANLIB = ranlib
#CFLAGS = -g -Wall -DNDEBUG
CFLAGS = -DNDEBUG
ARFLAGS = crv
TARGET = best-fit
TARGETLIB = libbest-fit.a
LIBS = -L/usr/lib -lboost_program_options

$(TARGET): main.o Shapes.o Double.o BestFit.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o Shapes.o Double.o BestFit.o $(LIBS)
main.o: main.cpp BestFit.h
	$(CC) $(CFLAGS) -c main.cpp

library: Shapes.o Double.o BestFit.o
	$(AR) $(ARFLAGS) $(TARGETLIB) Shapes.o Double.o BestFit.o
	$(RANLIB) $(TARGETLIB)

Shapes.o: Shapes.cpp BestFit.h
	$(CC) $(CFLAGS) -c Shapes.cpp
Double.o: Double.cpp BestFit.h
	$(CC) $(CFLAGS) -c Double.cpp
BestFit.o: BestFit.cpp BestFit.h
	$(CC) $(CFLAGS) -c BestFit.cpp

clean:
	$(RM) main.o Shapes.o BestFit.o Double.o
distclean:
	$(RM) $(TARGET) $(TARGETLIB) main.o Shapes.o BestFit.o Double.o
