# to install executable into $(BINDIR),
# and function library into $(LIBDIR),
#	type "make".

CC	= cc
AR	= ar
CFLAGS	= $(F)
OBJS	= hull.o ch.o io.o rand.o pointops.o fg.o geom_utils.o
HDRS	= hull.h points.h pointsites.h stormacs.h geom_utils.h
SRC	= hull.c ch.c io.c rand.c pointops.c fg.c geom_utils.c
LIB	= libhull.a


#all	: $(PROG) rsites
#	cp $(PROG) $(BINDIR)/.
#	cp rsites $(BINDIR)/.

all	: lib 


$(OBJS) : $(HDRS)

lib	: $(OBJS) 
	$(AR) rcv $(LIB) $(OBJS)

clean	:
	-rm -f $(OBJS) hullmain.o core a.out $(PROG)

