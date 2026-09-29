project('fstpack', 'c', 'fortran', version: '0.1.0',
        meson_version: '>=1.4.0')

py = import('python').find_installation(pure: false)
numpy = dependency('numpy')

f2py_include = run_command(
  py, '-c',
  'import numpy.f2py; print(numpy.f2py.get_include())',
  check: true
).stdout().strip()

wrappers = custom_target(
  'f2py-wrappers',
  input: 'python/st.pyf',
  output: [
    'fstpackmodule.c',
    'fstpack-f2pywrappers.f',
    'fstpack-f2pywrappers2.f90',
  ],
  command: [
    'sh', '-c',
    '"$1" -m numpy.f2py "$2" && touch "$3" "$4"',
    'sh', py, '@INPUT@', '@OUTPUT1@', '@OUTPUT2@',
  ],
)

legacy = static_library(
  'fftpack',
  run_command(
    'sh', '-c', 'printf "%s\n" src/fftpack/*.f',
    check: true
  ).stdout().strip().split('\n'),
  fortran_args: [
    '-std=legacy',
    '-ffixed-form',
    '-w',
  ],
  pic: true,
)

core = static_library(
  'fstpack-core',
  files(
    'src/mutl.f90',
    'src/fftpack.f90',
    'src/hilbrt.f90',
    'src/fstpack.f90',
    'python/st.f90',
  ),
  fortran_args: [
    '-std=f2018',
    '-Wall',
    '-pedantic',
  ],
  link_with: legacy,
  pic: true,
)

py.extension_module(
  'fstpack',
  [
    wrappers,
    f2py_include / 'fortranobject.c',
  ],
  dependencies: numpy,
  include_directories: include_directories(f2py_include),
  link_with: core,
  install: true,
)
