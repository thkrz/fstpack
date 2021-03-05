include config.mk

SRC = src/fftpack.f90 src/hilbrt.f90 src/fstpack.f90
FFTSRC = $(shell find ./src/fftpack -name '*.f')
OBJ = ${FFTSRC:.f=.o} ${SRC:.f90=.o}
LIB = -L./ -lfstpack

%.o: %.f
	@echo FC $<
	@${FC} -o $@ -c ${LEGACYFLAGS} $<

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: libfstpack

libfstpack: ${OBJ}
	@echo AR ${@}.a
	@${AR} rcs ${@}.a $^

fst: ${OBJ}
	@echo F2PY $@
	@python3 -m numpy.f2py -c -m $@ \
		--no-wrap-functions --fcompiler=gnu95 --f90exec=${FC} \
		--f90flags="${FFLAGS}" \
		-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION \
		--quiet ${OBJ} $@.f90

tests: libfstpack tfst1

tfst1: test/tfst1.o
	@echo LD $<
	@${FC} -o $@ $< ${LDFLAGS} ${LIB}

clean:
	find . \( -name '*.o' -or -name '*.mod' \) -exec rm {} \;
	rm -f *.a *.so

install:
	@echo installing

.PHONY: all clean install tests
.SUFFIXES: .o .f .f90
