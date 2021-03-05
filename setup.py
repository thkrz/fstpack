from numpy.distutils.core import Extension, setup

setup(
    name="fstpack",
    description="A Fortran library for Fast Stockwell transforms",
    author="Thomas Kreuzer",
    author_email="thomas.kreuzer@uni-wuerzburg.de",
    ext_modules=[Extension(name="fst", sources=["fst.f90"])],
)
