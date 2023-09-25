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

SPHINXBUILD = sphinx-build
SPHINXOPTS =

LEGACYFLAGS = -std=legacy -ffixed-form -w -O3
FFLAGS = -std=f2018 -ffree-form -fmax-errors=1 \
	-pedantic -Wall
LDFLAGS = -s -L./ -static -lfstpack

FFTSRC := $(wildcard ./src/fftpack/*.f)
SRC = src/mutl.f90 src/fftpack.f90 src/hilbrt.f90 src/fstpack.f90
OBJ = $(FFTSRC:.f=.o) $(SRC:.f90=.o)

%.o: %.f
	@echo FC $<
	@$(FC) -o $@ -c $(LEGACYFLAGS) $<

%.o: %.f90
	@echo FC $<
	@$(FC) -o $@ -c $(FFLAGS) $<

all: libfstpack pyfstpack tests

libfstpack: $(OBJ)
	@echo AR $(@).a
	@$(AR) rcs $(@).a $^
	@echo LD $(@).so
	@$(FC) -fPIC -shared -o $(@).so.$(VERSION) $^
	@[ -s $(@).so.$(SONUM) ] || ln -s $(@).so.$(VERSION) $(@).so.$(SONUM)
	@[ -s $(@).so ] || ln -s $(@).so.$(SONUM) $(@).so

help:
	$(SPHINXBUILD) -b html $(SPHINXOPTS) doc docs

tests:
	python3 -m unittest test.tfst

pyfstpack: python/st.pyf
	@echo F2PY $<
	@f2py3 --fcompiler=gnu95 --f90exec=$(FC) \
		-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION \
		libfstpack.a -c $< ${<:.pyf=.f90} --quiet

clean:
	rm -f $(OBJ) test/*.o libfstpack.a libfstpack.so* *.mod tfst* *.cpython*

install:
	install -m644 fstpack.mod $(DESTDIR)$(INCDIR)/fstpack.mod
	install -m644 libfstpack.a $(DESTDIR)$(LIBDIR)/libfstpack.a
	install -m644 libfstpack.so.$(VERSION) $(DESTDIR)$(LIBDIR)/libfstpack.so.$(VERSION)
	cp -P libfstpack.so.$(SONUM) $(DESTDIR)$(LIBDIR)/libfstpack.so.$(SONUM)
	cp -P libfstpack.so $(DESTDIR)$(LIBDIR)/libfstpack.so

.PHONY: all clean help install tests
