project('fluidity', 'fortran', 'c', 'cpp')

cxx = meson.get_compiler('cpp')
spatialindex_dep = cxx.find_library('spatialindex')
judy_dep = cxx.find_library('Judy')
scotch_deps = [
    cxx.find_library('parmetis'),
    cxx.find_library('ptscotchparmetisv3'),
    cxx.find_library('ptscotch'),
]
hypre_dep = cxx.find_library('HYPRE')
if hypre_dep.found()
    add_global_arguments('-DHAVE_HYPRE', language: ['fortran', 'cpp'])
endif

zoltan_dep = cxx.find_library('zoltan')
zoltan_deps = [zoltan_dep]
if zoltan_dep.found()
    add_global_arguments('-DHAVE_ZOLTAN', language: ['fortran', 'cpp'])
    zoltan_deps += cxx.find_library('metis')
endif

add_global_arguments('-DDOUBLEP', '-DHAVE_ADAPTIVITY', '-DHAVE_MBA_2D', language: ['fortran', 'cpp'])
add_global_arguments('-DSIGNAL', '-DSTREAM_IO', language: ['fortran', 'cpp'])
add_global_arguments('-fdefault-real-8', '-fdefault-double-8', '-ffree-line-length-none', '-ffixed-line-length-none', language: 'fortran')

fc = meson.get_compiler('fortran')
spud_f90_dep = fc.find_library('spud', has_headers: ['spud.h'])

inc = include_directories('include', 'libwm')
sys_inc = include_directories('/usr/local/include')

cgal = dependency('CGAL', required: false)
if cgal.found()
    add_global_arguments('-DHAVE_LIBCGAL', language: ['fortran', 'cpp'])
endif

hdf5_deps = [
    dependency('hdf5', language: 'c'),
    dependency('hdf5', language: 'fortran'),
]

mpi = dependency('mpi', language: 'fortran', required: false, method: 'config-tool')
if mpi.found()
    add_global_arguments('-DHAVE_MPI', language: ['fortran', 'cpp'])
endif

netcdf = dependency('netcdf', language: 'c', required: false)
if netcdf.found()
    add_global_arguments('-DHAVE_LIBNETCDF', language: ['fortran', 'cpp'])
endif

numa = dependency('numa', required: false)
if numa.found()
    add_global_arguments('-DHAVE_LIBNUMA', language: ['fortran', 'cpp'])
endif

openblas = dependency('openblas')

# for petsc in /usr/local, use its pkgconf name
#petsc = dependency('PETSc', required: false)
# for system petsc
petsc = dependency('petsc', required: false)
if petsc.found()
    add_global_arguments('-DHAVE_PETSC', language: ['fortran', 'cpp'])
endif

python_deps = []
python = dependency('python3-embed', required: false)
if python.found()
    python_deps += python
    python_deps += dependency('numpy')
    add_global_arguments('-DHAVE_PYTHON', '-DHAVE_NUMPY', language: ['fortran', 'c', 'cpp'])
endif

udunits = dependency('udunits', required: false)
if udunits.found()
    add_global_arguments('-DHAVE_LIBUDUNITS', language: ['fortran', 'c', 'cpp'])
endif

vtk = dependency(
    'VTK',
    modules: [
        'VTK::vtksys', 'VTK::CommonCore', 'VTK::CommonDataModel', 'VTK::FiltersVerdict',
	'VTK::IOXML', 'VTK::IOParallelXML', 'VTK::ParallelMPI', 'VTK::zlib'
    ]
)
if vtk.found()
    add_global_arguments('-DHAVE_VTK', language: ['cpp', 'fortran'])
endif

lib_debug = static_library(
    'debug',
    ['debug/Debug.F90', 'debug/Debug_Parameters.F90', 'debug/C++_Debug.cpp', 'debug/Diagnostic_Output.F90'],
    include_directories: inc,
    dependencies: [mpi],
)

lib_wm = static_library(
    'libwm',
    [
	'libwm/Wm4Intersector.cpp',
	'libwm/Wm4Intersector1.cpp',
	'libwm/Wm4IntrQuad2Quad2.cpp',
	'libwm/Wm4IntrTetrahedron3Tetrahedron3.cpp',
        'libwm/Wm4IntrTriangle2Triangle2.cpp',
	'libwm/Wm4Math.cpp',
	'libwm/Wm4System.cpp',
    ],
)

subdir('h5hut')
subdir('libadaptivity')
subdir('libmba2d')
subdir('libvtkfortran')
subdir('femtools')

lib_fluidity = static_library(
    'fluidity',
    ['main/mainfl.F90', 'main/Usage.cpp', 'main/Fluids.F90'],
    include_directories: [inc, sys_inc],
    dependencies: [hdf5_deps, mpi, numa, openblas, petsc, python_deps, udunits, spud_f90_dep, judy_dep, scotch_deps, zoltan_deps],
    link_with: [lib_adaptivity, lib_debug, lib_femtools, lib_wm],
)

executable(
    'fluidity',
    ['main.cpp'],
    include_directories: [inc],
    dependencies: [hdf5_deps, mpi, numa, openblas, petsc, python_deps, udunits, spud_f90_dep, judy_dep, scotch_deps, hypre_dep],
    link_with: lib_fluidity,
)

executable(
    'flredecomp',
    ['tools/Flredecomp.F90', 'tools/Flredecomp_main.cpp'],
    include_directories: [inc, sys_inc],
    dependencies: [mpi, spud_f90_dep],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)

executable(
    'fladapt',
    ['tools/Fladapt.F90', 'tools/Fladapt_main.cpp'],
    include_directories: [inc, sys_inc],
    dependencies: [mpi, spud_f90_dep, petsc, python_deps],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)

executable(
    'fldecomp',
    ['fldecomp/fldmain.cpp', 'fldecomp/fldgmsh.cpp', 'fldecomp/fldtriangle.cpp', 'fldecomp/partition.cpp'],
    include_directories: [inc, sys_inc],
    dependencies: [mpi, vtk],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)

executable(
    'gmsh2vtu',
    ['tools/gmsh2vtu.F90', 'tools/gmsh2vtu_main.cpp'],
    include_directories: [inc, sys_inc],
    dependencies: [mpi],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)

executable(
    'periodise',
    ['tools/periodise.F90'],
    include_directories: [inc, sys_inc],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)

executable(
    'project_vtu',
    ['tools/Project_Vtu.F90', 'tools/Project_Vtu_main.cpp'],
    include_directories: [inc, sys_inc],
    dependencies: [mpi, petsc],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)

executable(
    'test_pressure_solve',
    ['tools/test_pressure_solve.F90', 'tools/test_pressure_solve_main.cpp'],
    include_directories: [inc, sys_inc],
    dependencies: [mpi, petsc, python_deps],
    link_with: [lib_debug, lib_femtools, lib_fluidity],
)