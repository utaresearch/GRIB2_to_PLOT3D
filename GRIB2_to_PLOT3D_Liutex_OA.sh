#!/bin/bash

# Author: Oscar Alvarez
# Email: oscar.alvarez@uta.edu
# University of Texas at Arlington
# Department of Mathematics
# Arlington, Texas U.S.

######### Instructions ############
# Requires files:
# liutex_mod.f08
# grib_data_read.f08

# Grib2 files need to be in a data folder.

# INPUT NEEDED from command line (command line arguments):
# $1 = grib file name (without file name extension)
# $2 = record number of timestamp (usually 1)

# An EXAMPLE on how to run this file:
# $> ./GRIB2_to_PLOT3D_Liutex_OA.sh FABIAN10L_200309021830_Postr_Mean 1

#### Program Start ####

# FILE NAMES and RECORD NUMBER:
filename="data/$1"
ctl_filename="data/$1.ctl"
record_number="data/$2"
grib_grid_file="data/$1_grid.dat"
grib_velocity_file="data/$1_velocity.dat"

# GET LATITUDE GRID INFO:
wgrib -V -d $record_number $filename | grep "lat" | grep -oP '(?<= lat ).*?(?=to)|(?<=to).*?(?=by)' > $grib_grid_file

# GET LONGITUDE GRID INFO:
wgrib -V -d $record_number $filename | grep "long" | grep -oP '(?<=long).*?(?=to)|(?<=to).*?(?=by)' >> $grib_grid_file

# GET THE NUMBER OF LEVELS (HEIGHTS):
cat $ctl_filename | grep -oP "(?<=zdef).*?(?=levels)" >> $grib_grid_file

# GET HEIGHTS IN FORM OF mb = millibars GRID INFO:
cat $ctl_filename | grep "zdef" -A 1 | grep -v "zdef" >> $grib_grid_file

# GET THE ORDER OF THE VELOCITY COMPONENTS
wgrib $filename | grep -v "above gnd" | egrep -m 3 "(:UGRD:|:VGRD:|:DZDT:)" >> $grib_grid_file

# GET STARTING TIMESTAMP:
wgrib -d $record_number $filename | grep -oP "(?<=TimeU=).*?(?=:)" >> $grib_grid_file

# GET VELOCITY DATA:
wgrib $filename | egrep "(:UGRD:|:VGRD:|:DZDT:)" | grep -v "above gnd" | wgrib -i -grib $filename -text -o $grib_velocity_file

echo "GRIB file data extraction complete."

# COMPILE FORTRAN CODES
gfortran -c liutex_mod.f08
gfortran -c grib_data_read.f08
gfortran grib_data_read.o liutex_mod.o -o grib2liutex.o

# RUN FORTRAN CODE
./grib2liutex.o $1
