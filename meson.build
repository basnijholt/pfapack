project('pfapack',
  ['c', 'fortran'],
  version: '0.3.1',
  default_options: [
    'buildtype=release',
    'c_std=c99',
    'fortran_std=legacy',
  ],
)

# Dependencies
if host_machine.system() == 'darwin'
  add_project_link_arguments('-framework', 'Accelerate', language: ['c', 'fortran'])
  lapack_dep = declare_dependency()
  blas_dep = declare_dependency()
else
  lapack_dep = dependency('lapack')
  blas_dep = dependency('blas')
endif

pfapack_dir = 'original_source'

# Fortran sources in correct order matching the Makefile.FORTRAN
fortran_sources = files(
  # Single precision real
  join_paths(pfapack_dir, 'fortran', 'slasktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'slasktrf.f'),
  join_paths(pfapack_dir, 'fortran', 'sskbpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'sskbpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'sskbtrd.f'),
  join_paths(pfapack_dir, 'fortran', 'sskmv.f'),
  join_paths(pfapack_dir, 'fortran', 'sskpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'sskpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'sskr2.f'),
  join_paths(pfapack_dir, 'fortran', 'sskr2k.f'),
  join_paths(pfapack_dir, 'fortran', 'ssktd2.f'),
  join_paths(pfapack_dir, 'fortran', 'ssktf2.f'),
  join_paths(pfapack_dir, 'fortran', 'ssktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'ssktrf.f'),
  # Double precision real
  join_paths(pfapack_dir, 'fortran', 'dlasktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'dlasktrf.f'),
  join_paths(pfapack_dir, 'fortran', 'dskbpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'dskbpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'dskbtrd.f'),
  join_paths(pfapack_dir, 'fortran', 'dskmv.f'),
  join_paths(pfapack_dir, 'fortran', 'dskpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'dskpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'dskr2.f'),
  join_paths(pfapack_dir, 'fortran', 'dskr2k.f'),
  join_paths(pfapack_dir, 'fortran', 'dsktd2.f'),
  join_paths(pfapack_dir, 'fortran', 'dsktf2.f'),
  join_paths(pfapack_dir, 'fortran', 'dsktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'dsktrf.f'),
  # Single precision complex
  join_paths(pfapack_dir, 'fortran', 'clasktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'clasktrf.f'),
  join_paths(pfapack_dir, 'fortran', 'cskbpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'cskbpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'cskbtrd.f'),
  join_paths(pfapack_dir, 'fortran', 'cskmv.f'),
  join_paths(pfapack_dir, 'fortran', 'cskpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'cskpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'cskr2.f'),
  join_paths(pfapack_dir, 'fortran', 'cskr2k.f'),
  join_paths(pfapack_dir, 'fortran', 'csktd2.f'),
  join_paths(pfapack_dir, 'fortran', 'csktf2.f'),
  join_paths(pfapack_dir, 'fortran', 'csktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'csktrf.f'),
  # Double precision complex
  join_paths(pfapack_dir, 'fortran', 'zlasktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'zlasktrf.f'),
  join_paths(pfapack_dir, 'fortran', 'zskbpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'zskbpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'zskbtrd.f'),
  join_paths(pfapack_dir, 'fortran', 'zskmv.f'),
  join_paths(pfapack_dir, 'fortran', 'zskpf10.f'),
  join_paths(pfapack_dir, 'fortran', 'zskpfa.f'),
  join_paths(pfapack_dir, 'fortran', 'zskr2.f'),
  join_paths(pfapack_dir, 'fortran', 'zskr2k.f'),
  join_paths(pfapack_dir, 'fortran', 'zsktd2.f'),
  join_paths(pfapack_dir, 'fortran', 'zsktf2.f'),
  join_paths(pfapack_dir, 'fortran', 'zsktrd.f'),
  join_paths(pfapack_dir, 'fortran', 'zsktrf.f'),
  # Other files
  join_paths(pfapack_dir, 'fortran', 'mul10.f'),
  # Fortran 90 interface
  join_paths(pfapack_dir, 'fortran', 'precision.f90'),
  join_paths(pfapack_dir, 'fortran', 'f77_interface.f90'),
  join_paths(pfapack_dir, 'fortran', 'f95_interface.f90'),
  join_paths(pfapack_dir, 'fortran', 'message.f90'),
  join_paths(pfapack_dir, 'fortran', 'skpfa.f90'),
  join_paths(pfapack_dir, 'fortran', 'skpf10.f90'),
  join_paths(pfapack_dir, 'fortran', 'skbpfa.f90'),
  join_paths(pfapack_dir, 'fortran', 'skbpf10.f90'),
  join_paths(pfapack_dir, 'fortran', 'sktrd.f90'),
  join_paths(pfapack_dir, 'fortran', 'sktd2.f90'),
  join_paths(pfapack_dir, 'fortran', 'sktrf.f90'),
  join_paths(pfapack_dir, 'fortran', 'sktf2.f90'),
  join_paths(pfapack_dir, 'fortran', 'skbtrd.f90'),
)

# C sources
c_sources = files(
  join_paths(pfapack_dir, 'c_interface', 'skpfa.c'),
  join_paths(pfapack_dir, 'c_interface', 'skpf10.c'),
  join_paths(pfapack_dir, 'c_interface', 'skbpfa.c'),
  join_paths(pfapack_dir, 'c_interface', 'skbpf10.c'),
  join_paths(pfapack_dir, 'c_interface', 'sktrf.c'),
  join_paths(pfapack_dir, 'c_interface', 'sktrd.c'),
  join_paths(pfapack_dir, 'c_interface', 'skbtrd.c'),
)

# Build Fortran library
libpfapack = shared_library('pfapack',
  fortran_sources,
  dependencies: [lapack_dep, blas_dep],
  install: true,
)

# Build C library
libcpfapack = shared_library('cpfapack',
  c_sources,
  link_with: libpfapack,
  dependencies: [lapack_dep, blas_dep],
  install: true,
)

# Python module
py = import('python').find_installation()
py.install_sources(
  ['pfapack/__init__.py',
   'pfapack/_version.py',
   'pfapack/ctypes.py',
   'pfapack/pfaffian.py'],
  pure: false,
  subdir: 'pfapack'
)
