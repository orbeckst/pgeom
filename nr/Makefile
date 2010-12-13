# Makefile for numerical recipes lib

#DOUBLE_PREC := yes
ifeq ($(strip $(DOUBLE)),yes)
DOUBLE_PREC := yes
endif


LIBFILES := nrutil.c minimisation.c integrate3d.c
LIB_O := $(subst .c,.o,$(LIBFILES))
LIBHEADER := nrutil.h

ifdef DOUBLE_PREC
LIB_NR   := libnr_d.a
else
LIB_NR   := libnr.a
endif

# can use -DDOUBLE to compile in double precisison

CFLAGS  := -Wall  -g -O3 -ffast-math -funroll-loops -finline-functions 
LDFLAGS := -lm

ifdef DOUBLE_PREC
CFLAGS += -DDOUBLE
endif



$(LIB_NR): $(LIB_O)
	ar rv $@ $^
	ranlib $@

.phony: distclean clean
clean: 
	-rm -f *.o core *~ 
distclean: clean
	-rm -f *.a	
