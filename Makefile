FC = gfortran
AR = ar
PREFIX = /usr/local
LIBDIR = $(PREFIX)/lib
PYTHON := $(shell python3 -c 'import sys; i=sys.version_info; print(f"python{i.major}.{i.minor}")')
PYDIST = $(LIBDIR)/$(PYTHON)/dist-packages

FFLAGS = -std=f2018 -Wall -pedantic -O2 -fPIC -fmax-errors=1 -I. -J.
LEGACYFLAGS = -std=legacy -ffixed-form -w -O2 -fPIC

FFTSRC = $(wildcard src/fftpack/*.f)
SRC = src/mutl.f90 src/fftpack.f90 src/hilbrt.f90 src/fstpack.f90
OBJ = $(FFTSRC:.f=.o) $(SRC:.f90=.o)

all: pyfstpack

%.o: %.f
	$(FC) $(LEGACYFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

src/fftpack.o: src/mutl.o
src/hilbrt.o: src/mutl.o src/fftpack.o
src/fstpack.o: src/mutl.o src/fftpack.o src/hilbrt.o

libfstpack.a: $(OBJ)
	$(AR) rcs $@ $(OBJ)

pyfstpack: libfstpack.a python/st.pyf python/st.f90
	FC=$(FC) python3 -m numpy.f2py -c python/st.pyf python/st.f90 --f90flags="-I$(CURDIR)" -L$(CURDIR) -lfstpack
	touch $@

tests: pyfstpack
	python3 -m unittest test.tfst

clean:
	rm -f $(OBJ) *.mod *.o libfstpack.a libfstpack.so* fstpack*.so pyfstpack tfst* test/*.o
	rm -rf build test/__pycache__

install: pyfstpack
	mkdir -p $(DESTDIR)$(PYDIST)
	install -m644 fstpack*.so $(DESTDIR)$(PYDIST)

uninstall:
	rm -f $(DESTDIR)$(PYDIST)/fstpack*.so

.DELETE_ON_ERROR:
.PHONY: all tests clean install uninstall
