# $Id: Makefile,v 1.19 2009/01/11 23:01:42 oliver Exp $

PROG   := pgeom
DPROG  := -DPGEOM
MAIN := $(PROG).c xg.c
MAIN_H := xg.h 

DOUBLE := yes

BIN_DIR := $(HOME)/bin

#------------------------------------

COMMON_SRC := xg.c
COMMON_H   := xg.h 

PGEOM   := pgeom
D_PGEOM := -DPGEOM
PGEOM_SRC := $(PGEOM).c

LIBFILES = mol_data.c util.c cylint.c
LIB_H := $(subst .c,.h,$(LIBFILES))
LIB_O := $(subst .c,.o,$(LIBFILES))
ifdef DOUBLE
LIBPG := libpg_d.a
else
LIBPG := libpg.a
endif

SOURCE := $(MAIN) 
HEADERS := $(MAIN_H) $(LIB_H)

# numerical recipes routines

NR_D := ./nr

ifdef DOUBLE
LIBS_AUX := -L. -lpg_d -L$(NR_D) -lnr_d
else
LIBS_AUX := -L. -lpg -L$(NR_D) -lnr
endif
# can use -DDOUBLE
#CFLAGS := -I. -I../nr -Wall -Wno-unused  -g -DDEBUG 
CFLAGS := -I. -I$(NR_D) -Wall  -Wno-unused -g -O3 -ffast-math -funroll-loops -finline-functions   -DTESTCASE
ifdef DOUBLE
CFLAGS += -DDOUBLE
endif
LDFLAGS := -lm

CC = gcc

.PHONY: all install NR
all: $(PGEOM) 

NR:
	cd $(NR_D) && $(MAKE) DOUBLE=$(DOUBLE)


$(PGEOM): NR $(PGEOM_SRC) $(PGEOM_H) $(COMMON_SRC) $(COMMON_H) $(LIBPG) $(LIB_H)
	gcc $(CFLAGS) $(D_PGEOM) $(LDFLAGS) -o $@ $(SOURCE) $(LIBS_AUX) 


$(LIBPG): $(LIB_O)
	ar rv $@ $^
	ranlib $@

install: $(PGEOM)
	test -d $(BIN_DIR) || install -v -d $(BIN_DIR)
	install --mode=755 $< $(BIN_DIR)

.PHONY: test clean dist-clean
TEST_PDB := smallpore.pdb
TEST_ITP := smallpore.itp
TEST_FILES := $(TEST_ITP) $(TEST_PDB) pore.pdb pore.itp
test:
	$(PROG) -debug 80 -R 12 -P 2 4 -M 4 4 -o $(TEST_PDB) -s $(TEST_ITP)

clean:
	-rm -f $(TEST_FILES) *~ \#* *.o core
	cd $(NR_D) && $(MAKE) clean

distclean: clean
	-rm -f $(PROG) $(LIBPG)
	cd $(NR_D) && $(MAKE) distclean

proto.h: $(SOURCE)
	egrep '^(int|struct|void|double|char) .*\(.*\)' $< | sort -k 2 > $@.tmp

.PHONY: mtags ptags
mtags: $(MGEOM_SRC) $(MGEOM_H) $(COMMON_SRC) $(COMMON_H) $(LIBFILES) $(LIB_H)
	etags $^

ptags: $(PGEOM_SRC) $(PGEOM_H) $(COMMON_SRC) $(COMMON_H) $(LIBFILES) $(LIB_H)
	etags $^

TAGS: $(SOURCE) $(LIBFILES) $(HEADERS)
	etags $^
