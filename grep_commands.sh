#!/bin/bash

# INPUT NEEDED:
# + grib file name
# + record number of timestamp

# FILE NAMES and RECORD NUMBER:
filename="data/$1"
ctl_filename="data/${filename}.ctl"
record_number="data/$2"

# GET LATITUDE GRID INFO:
wgrib -V -d $record_number $filename | grep "lat" | grep -oP '(?<= lat ).*?(?=to)|(?<=to).*?(?=by)' > grib_grid_data.dat

# GET LONGITUDE GRID INFO:
wgrib -V -d $record_number $filename | grep "long" | grep -oP '(?<=long).*?(?=to)|(?<=to).*?(?=by)' >> grib_grid_data.dat

# GET THE NUMBER OF LEVELS (HEIGHTS):
cat $ctl_filename | grep -oP "(?<=zdef).*?(?=levels)" >> grib_grid_data.dat

# GET HEIGHTS IN FORM OF mb = millibars GRID INFO:
cat $ctl_filename | grep "zdef" -A 1 | grep -v "zdef" >> grib_grid_data.dat

# GET THE ORDER OF THE VELOCITY COMPONENTS
wgrib $filename | grep -v "above gnd" | egrep -m 3 "(:UGRD:|:VGRD:|:DZDT:)" >> grib_grid_data.dat

# GET STARTING TIMESTAMP:
wgrib -d $record_number $filename | grep -oP "(?<=TimeU=).*?(?=:)" >> grib_grid_data.dat

# GET VELOCITY DATA:
wgrib $filename | egrep "(:UGRD:|:VGRD:|:DZDT:)" | grep -v "above gnd" | wgrib -i -grib $filename -text -o grib_velocity_data

