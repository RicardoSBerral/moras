# This makefile was written to compile a distribution of pyfann for
# GNU platforms (cygwin included.)

TARGETS = _libfann.dll

PYTHON=python2.3
LIBS=-L. -L/usr/local/lib -L/usr/bin -l$(PYTHON) ../src/fann*.o

all: $(TARGETS)

_%.dll: %_wrap.o fann_helper.o
	gcc $(LIBS) -shared -dll $^ -o $@

%.o: %.c 
	gcc -c $< -I/usr/include/$(PYTHON)/

%_wrap.c: %.i 
	swig -python $<

clean:
	rm -f $(TARGETS) *_wrap.* fann_helper.o fann.pyc *.so
