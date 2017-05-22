#compilar con "make -f Makefile F=-g" para depurar
CC = g++
# CC = icc
CFLAGS = $(F) -ansi -pedantic -Wall
LIBGEN = ar

all: 
	make libreria F="$(F) -O3"

libreria: data20.o node20.o util20.o variable20.o fuzzyset20.o tree20.o definitions20.o

	$(LIBGEN) -r libCART.a data20.o node20.o util20.o variable20.o fuzzyset20.o tree20.o definitions20.o

data20.o: data20.cpp
	@echo Compiling... data20.cpp
	$(CC) $(CFLAGS) -c data20.cpp -I ../../nucleo/

node20.o: node20.cpp
	@echo Compiling... node20.cpp
	$(CC) $(CFLAGS) -c node20.cpp -I ../../nucleo/

util20.o: util20.cpp
	@echo Compiling... util20.cpp
	$(CC) $(CFLAGS) -c util20.cpp

variable20.o: variable20.cpp
	@echo Compiling... variable20.cpp
	$(CC) $(CFLAGS) -c variable20.cpp

#policy20.o: policy20.cpp
#	@echo Compiling... policy20.cpp
#	$(CC) $(CFLAGS) -c policy20.cpp
#
fuzzyset20.o: fuzzyset20.cpp
	@echo Compiling... fuzzyset20.cpp
	$(CC) $(CFLAGS) -c fuzzyset20.cpp

tree20.o: tree20.cpp
	@echo Compiling... tree20.cpp
	$(CC) $(CFLAGS) -c tree20.cpp -I ../../nucleo/

definitions20.o: definitions20.cpp
	@echo Compiling... definitions20.cpp
	$(CC) $(CFLAGS) -c definitions20.cpp

clean:
	rm -f *.o *.a
	
debug:
	make libreria F="$(F) -g"

