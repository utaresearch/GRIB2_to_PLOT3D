echo "Compiling fortran codes."

# COMPILE FORTRAN CODES
gfortran -c liutex_mod.f08
gfortran -c grib_data_read.f08
gfortran grib_data_read.o liutex_mod.o -o grib2liutex.o


echo "Running fortran programs."

# RUN FORTRAN CODE
./grib2liutex.o $1

read -p "Press enter to continue"