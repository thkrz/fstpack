.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f .f90

VERSION = 1.0
SONUM := $(shell echo $(VERSION) | cut -d '.' -f 1)

PREFIX = /usr/local
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
MANDIR = $(PREFIX)/share/man

AR = ar
FC = gfortran

LEGACYFLAGS = -std=legacy -ffixed-form -w -O3
FFLAGS = -std=f2008 -ffree-form -fmax-errors=1 \
	-pedantic -Wall -O3
LDFLAGS = -s -L./ -static -lfstpack

FFTSRC := $(wildcard ./src/fftpack/*.f)
SRC = src/fftpack.f90 src/hilbrt.f90 src/fstpack.f90
OBJ = $(FFTSRC:.f=.o) $(SRC:.f90=.o)

%.o: %.f
	@echo FC $<
	@$(FC) -o $@ -c $(LEGACYFLAGS) $<

%.o: %.f90
	@echo FC $<
	@$(FC) -o $@ -c $(FFLAGS) $<

all: libfstpack

libfstpack: $(OBJ)
	@echo AR $(@).a
	@$(AR) rcs $(@).a $^
	@echo LD $(@).so
	@$(FC) -fPIC -shared -o $(@).so.$(VERSION) $^
	@[ -s $(@).so.$(SONUM) ] || ln -s $(@).so.$(VERSION) $(@).so.$(SONUM)
	@[ -s $(@).so ] || ln -s $(@).so.$(SONUM) $(@).so

tests: libfstpack tfft1 #tfst1 tfst2

tfft1: test/tfft1.o
	@echo LD $<
	@$(FC) -o $@ $< $(LDFLAGS)

tfst1: test/tfst1.o
	@echo LD $<
	@$(FC) -o $@ $< $(LDFLAGS)

tfst2: test/tfst2.o
	@echo LD $<
	@$(FC) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(OBJ) test/*.o libfstpack.a libfstpack.so* *.mod tfst*

install:
	install -m644 fstpack.mod $(DESTDIR)$(INCDIR)/fstpack.mod
	install -m644 libfstpack.a $(DESTDIR)$(LIBDIR)/libfstpack.a
	install -m644 libfstpack.so.$(VERSION) $(DESTDIR)$(LIBDIR)/libfstpack.so.$(VERSION)
	cp -P libfstpack.so.$(SONUM) $(DESTDIR)$(LIBDIR)/libfstpack.so.$(SONUM)
	cp -P libfstpack.so $(DESTDIR)$(LIBDIR)/libfstpack.so

.PHONY: all clean install tests
