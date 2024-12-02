h5hut_inc = include_directories('src/include', 'src/h5core')

write_interface = find_program('write_interface.sh')

h5hut_f90 = files('src/include/H5hut.f90')

lib_h5hut = static_library(
    'h5hut',
    [
        'src/h5core/h5_attachments.c',
	'src/h5core/h5_attribs.c',
	'src/h5core/h5_err.c',
	'src/h5core/h5_file.c',
	'src/h5core/h5_init.c',
	'src/h5core/h5_log.c',
	'src/h5core/h5_model.c',
	'src/h5core/h5_syscall.c',
	'src/h5core/h5b_io.c',
	'src/h5core/h5b_model.c',
	'src/h5core/h5u_io.c',
	'src/h5core/h5u_model.c',
	'src/h5core/private/h5_hdf5.c',
        'src/Fortran/H5.c',
	'src/Fortran/H5_attribs.c',
	'src/Fortran/H5_model.c',
	'src/Fortran/H5Part_io.c',
	'src/Fortran/H5Part_model.c',
    ],
    include_directories: h5hut_inc,
    dependencies: hdf5_deps,
    c_args: [
        '-DFC_MANGLING_LOWERCASE_SINGLE_UNDERSCORE',
	'-DPACKAGE_VERSION=""',
    ],
)