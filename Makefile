CC = g++

all:
	make debug
	@echo ------------------------------------------------------ OK -----
	@echo ---------------------------------------------------------------

runrunrun:
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------- FES/Extracto -----
	make -C ./FES/Extracto/ F="$(F)"
	@echo ---------------------------------------------------------------
	@echo ---------------------------------------------------- C4.5 -----
	make -C ./R8/SrcAdapted/ F="$(F)"
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------------- nucleo -----
	make -C ./nucleo/ F="$(F) -static"
	@echo ------------------------------------------------------ OK -----
	@echo ---------------------------------------------------------------

runrunrun_sdp:
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------- FES/Extracto -----
	make -C ./FES/Extracto/ F="$(F)"
	@echo ---------------------------------------------------------------
	@echo ---------------------------------------------------- C4.5 -----
	make -C ./R8/SrcAdapted/ F="$(F)"
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------------- nucleo -----
	make -C ./nucleo/ F="$(F) -static -DCOMP_SDP" all_sdp
	@echo ------------------------------------------------------ OK -----
	@echo ---------------------------------------------------------------

debug_ef:
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------- FES/Extracto -----
	make -C ./FES/Extracto/ debug
	@echo ---------------------------------------------------------------
	@echo ---------------------------------------------------- C4.5 -----
	make -C ./R8/SrcAdapted/ debug
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------------- nucleo -----
	make -C ./nucleo/ debug_ef

debug:
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------- FES/Extracto -----
	make -C ./FES/Extracto/ debug
	@echo ---------------------------------------------------------------
	@echo ---------------------------------------------------- C4.5 -----
	make -C ./R8/SrcAdapted/ debug
	@echo ---------------------------------------------------------------
	@echo -------------------------------------------------- nucleo -----
	make -C ./nucleo/ F="$(F) -DCOMP_SDP" debug

clean:
	make -C ./FES/Extracto/ clean
	make -C ./R8/SrcAdapted/ clean
	make -C ./nucleo/ clean

backup:
	backup.sh

utils:
	$(CC) -O3 mean.c -pedantic -ansi -Wall -o mean nucleo/Matriz.cpp \
	-I NR++/ -I nucleo/ NR++/nrr.cpp -static
	$(CC) -O3 traspose.c -pedantic -ansi -Wall -o traspose nucleo/Matriz.cpp \
	-I NR++/ -I nucleo/ NR++/nrr.cpp -static
	$(CC) -O3 matop.c -pedantic -ansi -Wall -o matop nucleo/Matriz.cpp \
	-I NR++/ -I nucleo/ NR++/nrr.cpp -static
	$(CC) -O3 even.cpp -pedantic -ansi -Wall -o even
	$(CC) -O3 odd.cpp -pedantic -ansi -Wall -o odd
	$(CC) -O3 every.c -pedantic -ansi -Wall -o every
	$(CC) -O3 angloma.c -pedantic -ansi -Wall -o angloma -I nucleo/ -I NR++/ \
	NR++/nrr.cpp nucleo/Matriz.cpp -static
	@echo ---------------------------------------------------------------
	@echo --------------------------------------- Libreria dinamica -----
	make -C ./viewer/int/

analisis:
	$(CC) -g analisis.c -pedantic -ansi -o run nucleo/Matriz.cpp \
	-I NR++/ -I nucleo/ NR++/nrr.cpp

