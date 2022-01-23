gfortran -c -I/home/wang/tool/netcdf480/include REASv1.F mo_util.F reastype.F
gfortran -o REASrun REASv1.o mo_util.o reastype.o -L/home/wang/tool/netcdf480/lib -lnetcdff -lnetcdf
