# This makefile is on purpose not made with configure, to show how to use the library
# The make file requires that the fann library is installed (see ../README)

GCC=gcc

TARGETS = xor_train xor_test xor_test_fixed simple_train steepness_train simple_test robot mushroom
DEBUG_TARGETS = xor_train_debug xor_test_debug xor_test_fixed_debug

all: $(TARGETS)

%: %.c Makefile
	$(GCC) -O3 -lm -lfann $< -o $@

%_fixed: %.c Makefile
	$(GCC) -O3 -lm -lfixedfann -DFIXEDFANN $< -o $@

clean:
	rm -f $(TARGETS) $(DEBUG_TARGETS) xor_fixed.data *.net *~

runtest: $(TARGETS)
	@echo
	@echo Training network
	./xor_train

	@echo
	@echo Testing network with floats
	./xor_test

	@echo
	@echo Testing network with fixed points
	./xor_test_fixed

#below this line is only for debugging the fann library

rundebugtest: $(DEBUG_TARGETS)
	@echo
	@echo Training network
	./xor_train_debug

	@echo
	@echo Testing network with floats
	./xor_test_debug

	@echo
	@echo Testing network with fixed points
	./xor_test_fixed_debug

#compiletest is used to test whether the library will compile easily in other compilers
compiletest:
	gcc -O3 -ggdb -lm -DDEBUG -Wall -Wformat-security -Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wsign-compare -pedantic -ansi -I../src/ -I../src/include/ ../src/floatfann.c xor_train.c -o xor_train
	gcc -O3 -ggdb -lm -DDEBUG -Wall -Wformat-security -Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wsign-compare -pedantic -ansi -DFIXEDFANN -I../src/ -I../src/include/ ../src/fixedfann.c xor_test.c -o xor_test
	g++ -O3 -ggdb -lm -DDEBUG -Wall -Wformat-security -Wfloat-equal -Wpointer-arith -Wcast-qual -Wsign-compare -pedantic -ansi -I../src/ -I../src/include/ ../src/floatfann.c xor_train.c -o xor_train

debug: $(DEBUG_TARGETS)

%_debug: %.c Makefile ../src/*c ../src/include/*h
	$(GCC) -O3 -ggdb -lm -DDEBUG -Wall -ansi -I../src/ -I../src/include/ ../src/floatfann.c $< -o $@

%_fixed_debug: %.c Makefile
	$(GCC) -O3 -ggdb -lm -DDEBUG -Wall -ansi -DFIXEDFANN -I../src/ -I../src/include/ ../src/fixedfann.c $< -o $@

rundebug: $(DEBUG_TARGETS)
	@echo
	@echo Training network
	valgrind --leak-check=yes --show-reachable=yes --leak-resolution=high ./xor_train_debug

	@echo
	@echo Testing network with floats
	valgrind --leak-check=yes --show-reachable=yes --leak-resolution=high ./xor_test_debug

	@echo
	@echo Testing network with fixed points
	valgrind --leak-check=yes --show-reachable=yes --leak-resolution=high ./xor_test_fixed_debug
