gfortran -I/home/jihenghu/rttov13/mod -I/home/jihenghu/rttov13/include -fPIC -O3 -fopenmp -ffree-line-length-none  -c example_k.f90 -o ./example_k.o
gfortran -o ./example_k.exe ./example_k.o \
-L/home/jihenghu/rttov13/lib -lrttov13_brdf_atlas -lrttov13_emis_atlas -lrttov13_mw_scatt -lrttov13_other -lrttov13_coef_io -lrttov13_hdf -lrttov13_parallel -lrttov13_main  \
-L/home/jihenghu/netcdf/lib -lnetcdff -L/home/jihenghu/hdf5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -fopenmp -L/home/jihenghu/eccodes/lib64 -leccodes_f90 -leccodes
