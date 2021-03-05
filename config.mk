PREFIX=/usr/local

AR = ar
#FC = nvfortran
FC = gfortran

LEGACYFLAGS = -std=legacy -ffixed-form -w -O3
FFLAGS = -std=f2008 \
	-ffree-form -fmax-errors=1 \
	-pedantic -Wall \
	-O3
LDFLAGS = -s
