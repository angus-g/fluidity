inc_vtkfortran = include_directories('include')

lib_vtkfortran = static_library(
    'vtkfortran',
    ['vtkfortran.cpp', 'fvtkfortran.F90', 'vtkmeshio.cpp'],
    include_directories: inc_vtkfortran,
    dependencies: [vtk],
)