#!/bin/bash

# FILE NAMES and RECORD NUMBER:
filename="data/$1"
ctl_filename="data/$1.ctl"
record_number="$2"
grib_grid_file="$1_grid.dat"
grib_velocity_file="$1_velocity.dat"

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
echo " "
echo "$1"