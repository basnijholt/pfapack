project('pfapack',
  ['c', 'fortran'],
  meson_version: '>= 1.1.0',
  default_options: [
    'buildtype=release',
    'c_std=c99',
    'fortran_std=legacy',
  ],
  version: run_command('pfapack/_version.py', check: true).stdout().strip(),
)

# Get compilers
c_compiler = meson.get_compiler('c')
fortran_compiler = meson.get_compiler('fortran')

# Dependencies
thread_dep = dependency('threads')
m_dep = c_compiler.find_library('m', required: false)  # math library

if host_machine.system() == 'darwin'
  add_project_link_arguments('-framework', 'Accelerate', language: ['c', 'fortran'])
  lapack_dep = declare_dependency()
  blas_dep = declare_dependency()
elif host_machine.system() == 'windows'
  # On Windows, find OpenBLAS which includes LAPACK
  openblas_lib = dependency('openblas', required: false)
  if not openblas_lib.found()
    # Fallback to manual detection
    compiler = meson.get_compiler('c')
    openblas_lib = compiler.find_library('libopenblas',
      dirs: ['C:/msys64/mingw64/lib'],
      required: true)
  endif
  # OpenBLAS includes both BLAS and LAPACK
  lapack_dep = declare_dependency(dependencies: [openblas_lib])
  blas_dep = declare_dependency(dependencies: [openblas_lib])

else
  # Add rt library for clock_gettime
  rt_dep = c_compiler.find_library('rt', required: false)
  if rt_dep.found()
    add_project_link_arguments('-lrt', language: ['c', 'fortran'])
  endif

    # Add rt library for clock_gettime
    rt_dep = c_compiler.find_library('rt', required: false)
    if rt_dep.found()
        add_project_link_arguments('-lrt', language: ['c', 'fortran'])
    endif

    # Try multiple methods to find OpenBLAS
    openblas_dep = dependency('openblas', required: false)
    if not openblas_dep.found()
        # Fallback to manual detection
        openblas_dep = c_compiler.find_library('openblas', required: false)
    endif

    if openblas_dep.found()
        # OpenBLAS includes both BLAS and LAPACK
        blas_dep = declare_dependency(dependencies: [openblas_dep])
        lapack_dep = declare_dependency(dependencies: [openblas_dep])
    else
        # Fall back to separate BLAS and LAPACK
        blas_dep = c_compiler.find_library('blas', required: true)
        lapack_dep = c_compiler.find_library('lapack', required: true)
    endif
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

# Import Python module
py = import('python').find_installation(pure: false)

# Build Fortran library
libpfapack = shared_library('pfapack',
  fortran_sources,
  dependencies: [lapack_dep, blas_dep, thread_dep, m_dep],
  install: true,
  install_dir: py.get_install_dir() / 'pfapack',  # Install in package directory
)

# Include directories
c_include_dir = include_directories(
  join_paths(pfapack_dir, 'c_interface'),
  join_paths(pfapack_dir, 'c_interface', 'TESTING'),
)

# Build C library
libcpfapack = shared_library('cpfapack',
  c_sources,
  include_directories: c_include_dir,
  link_with: libpfapack,
  dependencies: [lapack_dep, blas_dep, thread_dep, m_dep],
  install: true,
  install_dir: py.get_install_dir() / 'pfapack',  # Install in package directory
)

# Python module
py.install_sources(
  ['pfapack/__init__.py',
   'pfapack/_version.py',
   'pfapack/ctypes.py',
   'pfapack/pfaffian.py'],
  pure: false,
  subdir: 'pfapack'
)


# Build test executables
if get_option('build_tests').allowed()
  # C interface tests
  c_test_sources = files(
    join_paths(pfapack_dir, 'c_interface', 'TESTING', 'test.c'),
    join_paths(pfapack_dir, 'c_interface', 'TESTING', 'check.f90'),
  )

  c_test_exe = executable('test_c_interface',
    c_test_sources,
    include_directories: c_include_dir,
    link_with: [libpfapack, libcpfapack],
    dependencies: [lapack_dep, blas_dep, thread_dep, m_dep],
  )

  # Common Fortran test sources
  fortran_test_common = files(
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_decomp.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_decomp_c.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_decomp_d.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_decomp_s.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_decomp_z.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_pfaffian.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_pfaffian_c.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_pfaffian_d.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_pfaffian_s.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'check_pfaffian_z.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'error.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'matrix_tools.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'matrix_tools_c.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'matrix_tools_d.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'matrix_tools_s.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'matrix_tools_z.f90'),
    join_paths(pfapack_dir, 'fortran', 'TESTING', 'test_ilaenv.f'),
  )

  # Test pfaff executable
  test_pfaff_exe = executable('test_pfaff',
    [fortran_test_common,
     join_paths(pfapack_dir, 'fortran', 'TESTING', 'test_pfaff.f90')],
    link_with: [libpfapack],
    dependencies: [lapack_dep, blas_dep, thread_dep, m_dep],
  )

  # Test skew executable
  test_skew_exe = executable('test_skew',
    [fortran_test_common,
     join_paths(pfapack_dir, 'fortran', 'TESTING', 'test_skew.f90')],
    link_with: [libpfapack],
    dependencies: [lapack_dep, blas_dep, thread_dep, m_dep],
  )

  # Register tests
  test('C interface tests', c_test_exe)
  test('Fortran pfaff tests', test_pfaff_exe)
  test('Fortran skew tests', test_skew_exe)
endif
