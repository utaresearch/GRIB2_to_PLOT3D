program main
    !!!======================================================================
    !!! Using the grib data files produced by grep to calculate Liutex.
    !!!
    !!! Author: Oscar Alvarez
    !!! Email: oscar.alvarez@uta.edu
    !!! University of Texas at Arlington
    !!! Arlington, Texas U.S.
    !!!======================================================================
    use liutex_mod

    implicit none

    real, parameter :: pi = 4.0 * atan(1.0)

    !! External functions
    real, external :: latlong_dist

    !! file handling variables
    integer, parameter :: file1=10, file2=20, file3=30, file4=40
    character(100) :: filename, grib_grid_filename, grib_velocity_filename
    character(500) :: first_velocity_component, second_velocity_component, third_velocity_component
    character(10) :: u_string, v_string, w_string
    character(100) :: grid_output_filename, fun_output_filename, data_output_filename
    real, dimension(:,:,:,:), allocatable :: f, f2
    real, dimension(:,:,:), allocatable :: lat, long, pressure
    integer, dimension(:), allocatable :: height
    integer :: n_vars, n_data_variables

    !! other variables
    real, dimension(:,:,:,:,:), allocatable :: velocity_gradient

    real, dimension(:,:,:,:), allocatable :: liutex_vector, mod_omega_liutex_vec
    real, dimension(:,:,:,:), allocatable :: liutex_magnitude_gradient, mod_omega_liutex_mag_gradient
    
    real, dimension(:,:,:), allocatable :: first, second, third
    real, dimension(:,:,:), allocatable :: u, v, w, x, y, z
    real, dimension(:,:,:), allocatable :: liutex_magnitude, mod_omega_liutex_mag
    
    real, dimension(3,3) :: a

    real :: lat_start, lat_end, lat_step, long_start, long_end, long_step
    real :: dist_x, dist_y
    real :: height_in_meters
    real :: velocity_tolerance
    real :: liutex_mag
    real :: undefined_value

    integer :: n_levels, time_start
    integer :: nx, ny
    integer :: i, j, k, r


    !! Input file names
    call get_command_argument(1, filename)
    grib_grid_filename = "output_data/" // trim(filename) //"_grid.dat"
    grib_velocity_filename = "output_data/" // trim(filename) //"_velocity.dat"


    !! Read grib grid data file.
    write(*,*)
    write(*,*) "Reading grib grid data file."
    open(file1, file=trim(grib_grid_filename), form="formatted", action="read")
    
    !! Reading latitude and longitude and their step sizes.
    read(file1, *) lat_start
    read(file1, *) lat_end
    read(file1, *) long_start
    read(file1, *) long_end

    !! Getting heights pressures from grid data file.
    read(file1, *) n_levels
    
    allocate(height(n_levels))
    
    read(file1, *) (height(i), i = 1, n_levels)

    !! Determine the order of the velocity components based on the first three records in which they appear.
    read(file1, *) first_velocity_component
    read(file1, *) second_velocity_component
    read(file1, *) third_velocity_component

    !! Initial time of beginning record.
    read(file1, *) time_start
    
    close(file1)


    !! Read the velocity values from the grib velocity data file.
    write(*,*) "Reading grib velocity data file."
    open(file2, file=trim(grib_velocity_filename), form="formatted", action="read")
    
    !! Dimensions of x and y grid directions.
    read(file2, *) nx, ny

    allocate(first(nx, ny, n_levels))
    allocate(second(nx, ny, n_levels))
    allocate(third(nx, ny, n_levels))

    !! Reading data from file.
    do k = 1, n_levels - 1

        do j = 1, ny
            do i = 1, nx
                read(file2, *) first(i,j,k)
            end do
        end do

        read(file2, *) nx, ny

        do j = 1, ny
            do i = 1, nx
                read(file2, *) second(i,j,k)
            end do
        end do

        read(file2, *) nx, ny

        do j = 1, ny
            do i = 1, nx
                read(file2, *) third(i,j,k)
            end do
        end do

        read(file2, *) nx, ny

    end do

    do j = 1, ny
        do i = 1, nx
            read(file2, *) first(i,j,n_levels)
        end do
    end do

    read(file2, *) nx, ny

    do j = 1, ny
        do i = 1, nx
            read(file2, *) second(i,j,n_levels)
        end do
    end do

    read(file2, *) nx, ny

    do j = 1, ny
        do i = 1, nx
            read(file2, *) third(i,j,n_levels)
        end do
    end do

    close(file2)

    write(*,*) "Finished file reads."
    write(*,*) "Data Dimensions (i,j,k): (", nx, ny, n_levels, ")"
    write(*,*) "Setting velocity components and calcuating grid."

    !! Set velocity components based on grib file order.
    allocate(u(nx, ny, n_levels))
    allocate(v(nx, ny, n_levels))
    allocate(w(nx, ny, n_levels))

    write(*,*)
    write(*,*) 'Velocity Order From Grib Output: '
    write(*,*) trim(first_velocity_component)
    write(*,*) trim(second_velocity_component)
    write(*,*) trim(third_velocity_component)
    
    write(*,*)
    write(*,*) 'Program Determined Velocity Order: '
    
    u_string = ':UGRD:'
    v_string = ':VGRD:'
    w_string = ':DZDT:'

    if (index(first_velocity_component, trim(u_string)) > 0) then
        u = first
        write(*,*) '1st: u'
    else if (index(second_velocity_component, trim(u_string)) > 0) then
        u = second
        write(*,*) '2nd: u'
    else if (index(third_velocity_component, trim(u_string)) > 0) then
        u = third
        write(*,*) '3rd: u'
    end if

    if (index(first_velocity_component, trim(v_string)) > 0) then
        v = first
        write(*,*) '1st: v'
    else if (index(second_velocity_component, trim(v_string)) > 0) then
        v = second
        write(*,*) '2nd: v'
    else if (index(third_velocity_component, trim(v_string)) > 0) then
        v = third
        write(*,*) '3rd: v'
    end if

    if (index(first_velocity_component, trim(w_string)) > 0) then
        w = first
        write(*,*) '1st: w'
    else if (index(second_velocity_component, trim(w_string)) > 0) then
        w = second
        write(*,*) '2nd: w'
    else if (index(third_velocity_component, trim(w_string)) > 0) then
        w = third
        write(*,*) '3rd: w'
    end if

    deallocate(first, second, third)

    write(*,*)

    !! Filtering out undefined data
    velocity_tolerance = 1.e20

    do k = 1, n_levels
        do j = 1, ny
            do i = 1, nx

                if (u(i,j,k) > velocity_tolerance) then
                    if ( isnan(u(i,j,k)) ) then
                        write(*,*) "ERROR in fortran main code: Velocity (u) is NaN value after data file read."
                        stop
                    end if
                    u(i,j,k) = 0.0
                end if
                
                if (v(i,j,k) > velocity_tolerance) then
                    if ( isnan(v(i,j,k)) ) then
                        write(*,*) "ERROR in fortran main code: Velocity (v) is NaN value after data file read."
                        stop
                    end if
                    v(i,j,k) = 0.0
                end if
                
                if (w(i,j,k) > velocity_tolerance) then
                    if ( isnan(w(i,j,k)) ) then
                        write(*,*) "ERROR in fortran main code: Velocity (w) is NaN value after data file read."
                        stop
                    end if
                    w(i,j,k) = 0.0
                end if

            end do
        end do
    end do
    
    !! Make latitute and longitude values.
    allocate(lat(nx, ny, n_levels))
    allocate(long(nx, ny, n_levels))

    lat_step = (lat_end - lat_start) / nx
    long_step = (long_end - long_start) / ny

    lat(1,:,:) = lat_start
    do i = 2, nx
        lat(i,:,:) = lat(i-1,:,:) + lat_step
    end do
    lat(nx,:,:) = lat_end

    long(:,1,:) = long_start
    do j = 2, ny
        long(:,j,:) = long(:,j-1,:) + long_step
    end do
    long(:,ny,:) = long_end

    ! !! Create 3D grid using meters instead of longitude and latitude
    ! allocate(x(nx, ny, n_levels))
    ! allocate(y(nx, ny, n_levels))
    ! allocate(z(nx, ny, n_levels))

    ! x(1,1,:) = 0.0
    ! y(1,1,:) = 0.0

    ! do j = 2, ny
    !     do i = 2, nx

    !         !! Return distance between two latitude and longitude points in meters
    !         dist_x = latlong_dist(lat(i-1, j), long(i, j), lat(i,j), long(i,j))
    !         dist_y = latlong_dist(lat(i, j), long(i, j-1), lat(i,j), long(i,j))

    !         x(i,j,:) = x(i-1, j, 1) + dist_x
    !         y(i,j,:) = y(i, j-1, 1) + dist_y
            
    !     end do
    ! end do

    !! Re-shaping height/pressure array in order to be added to function file.
    allocate(pressure(nx,ny,n_levels))

    do k = 1, n_levels
        pressure(:,:,k) = height(k)
    end do

    deallocate(height)

    !! Calculate the velocity gradient tensor for each point.
    write(*,*) "Calculating Velocity Gradient."
    allocate(velocity_gradient(nx, ny, n_levels, 3, 3))

    velocity_gradient = velocity_gradient_tensor(u, v, w, lat, long, pressure, nx, ny, n_levels)

    write(*,*) 'vel grad in main: ', velocity_gradient(nx/2, ny/2, n_levels/2, 1, 1)

    !! Calculate Liutex for each point.
    write(*,*) "Calculating Liutex Vector."
    allocate(liutex_vector(nx, ny, n_levels, 3))
    allocate(liutex_magnitude(nx, ny, n_levels))
    
    liutex_vector    = 0.0
    liutex_magnitude = 0.0

    liutex_vector = liutex(velocity_gradient, nx, ny, n_levels)

    write(*,*) "Calculating Modified Omega Liutex Vector."
    allocate(mod_omega_liutex_vec(nx, ny, n_levels, 3))
    allocate(mod_omega_liutex_mag(nx, ny, n_levels))

    ! mod_omega_liutex_vec = 0.0
    ! mod_omega_liutex_mag = 0.0

    mod_omega_liutex_vec = modified_omega_liutex(velocity_gradient, nx, ny, n_levels)

    write(*,*) 'mod omega in main: ', mod_omega_liutex_mag(nx/2, ny/2, n_levels/2), mod_omega_liutex_vec(nx/2, ny/2, n_levels/2, 1)

    deallocate(velocity_gradient)

    !! Calculating vector magnitude values.
    write(*,*) "Calculating Vector Magnitude Values."

    do k = 1, n_levels
        do j = 1, ny
            do i = 1, nx

                liutex_magnitude(i,j,k) = sqrt(liutex_vector(i,j,k,1)**2 + liutex_vector(i,j,k,2)**2   &
                                               + liutex_vector(i,j,k,3)**2)
                
                mod_omega_liutex_mag(i,j,k) = sqrt(mod_omega_liutex_vec(i,j,k,1)**2 + mod_omega_liutex_vec(i,j,k,2)**2  &
                                                   + mod_omega_liutex_vec(i,j,k,3)**2)
            
            end do
        end do
    end do


    !! Calculating magnitude gradients
    write(*,*) "Calculating Magnitude Gradient Vectors"
    allocate(liutex_magnitude_gradient(nx, ny, n_levels, 3))
    allocate(mod_omega_liutex_mag_gradient(nx, ny, n_levels, 3))

    write(*,*) "liutex mag gradient"
    liutex_magnitude_gradient = gradient(liutex_magnitude, lat, long, pressure, nx, ny, n_levels)

    write(*,*) "modified omega liutex mag gradient"
    mod_omega_liutex_mag_gradient = gradient(mod_omega_liutex_mag, lat, long, pressure, nx, ny, n_levels)


    !! Writing PLOT3D output files
    write(*,*) "Creating and writing PLOT3D grid (.xyz) and function (.fun) files."

    grid_output_filename = 'output_data/' // trim(filename) // '_liutex.xyz'
    fun_output_filename  = 'output_data/' // trim(filename) // '_liutex.fun'

    !! Writing grid file (.xyz)
    open(file2, file=trim(grid_output_filename), form='unformatted', action='write')
    write(file2) nx, ny, n_levels
    write(file2) (((lat(i, j, k),      i=1,nx), j=1,ny), k=1,n_levels),    &
                 (((long(i, j, k),     i=1,nx), j=1,ny), k=1,n_levels),    &
                 (((pressure(i, j, k), i=1,nx), j=1,ny), k=1,n_levels)

    close(file2)


    !! Writing function file (.fun)
    n_vars = 17

    allocate(f(nx, ny, n_levels, n_vars))

    f(:,:,:,1)  = u
    f(:,:,:,2)  = v
    f(:,:,:,3)  = w
    
    f(:,:,:,4)  = liutex_vector(:,:,:,1)
    f(:,:,:,5)  = liutex_vector(:,:,:,2)
    f(:,:,:,6)  = liutex_vector(:,:,:,3)
    
    f(:,:,:,7)  = liutex_magnitude

    f(:,:,:,8)  = liutex_magnitude_gradient(:,:,:,1)
    f(:,:,:,9)  = liutex_magnitude_gradient(:,:,:,2)
    f(:,:,:,10) = liutex_magnitude_gradient(:,:,:,3)
    
    f(:,:,:,11) = mod_omega_liutex_vec(:,:,:,1)
    f(:,:,:,12) = mod_omega_liutex_vec(:,:,:,2)
    f(:,:,:,13) = mod_omega_liutex_vec(:,:,:,3)
    
    f(:,:,:,14) = mod_omega_liutex_mag
    
    f(:,:,:,15) = mod_omega_liutex_mag_gradient(:,:,:,1)
    f(:,:,:,16) = mod_omega_liutex_mag_gradient(:,:,:,2)
    f(:,:,:,17) = mod_omega_liutex_mag_gradient(:,:,:,3)
    

    open(file3, file=trim(fun_output_filename), form='unformatted', action='write')
    write(file3) nx, ny, n_levels, n_vars
    
    write(file3) ((((f(i, j, k, r), i=1,nx), j=1,ny), k=1,n_levels), r=1,n_vars)
    
    close(file3)

    deallocate(f)


    !! Creating data (.dat) output file
    data_output_filename  = 'output_data/' // trim(filename) // '_liutex_UTA.dat'

    n_data_variables = 20

    write(*,*) " "
    write(*,*) "Creating and writing data (.dat) file."
    write(*,*) " "
    write(*,*) "(.dat) data file format/structure."
    write(*,*) "Line 1 - Dimensions (imax, jmax, kmax, n_data_variables): "
    write(*,*) nx, ny, n_levels, n_data_variables
    write(*,*) " "
    write(*,*) "Line 2 - All the data:"
    write(*,*) "Data Structure: "
    write(*,*) "================================================="
    write(*,*) "for variable_index = 1 to n_data_variables"
    write(*,*) "  for k = 1 to k_max"
    write(*,*) "    for j = 1 to j_max"
    write(*,*) "      for i = 1 to i_max"
    write(*,*) " "
    write(*,*) "        data(i, j, k, variable_index) = data_value"
    write(*,*) " "
    write(*,*) "end all for loops"
    write(*,*) "================================================="
    write(*,*) " "
    write(*,*) "Variables List:"
    write(*,*) " 1: latitute"
    write(*,*) " 2: longitude"
    write(*,*) " 3: pressure level"
    write(*,*) " 4: u velocity"
    write(*,*) " 5: v velocity"
    write(*,*) " 6: w velocity"
    write(*,*) " 7: liutex vector x-direction"
    write(*,*) " 8: liutex vector y-direction"
    write(*,*) " 9: liutex vector z-direction"
    write(*,*) "10: liutex magnitude"
    write(*,*) "11: liutex magnitude vector x-direction"
    write(*,*) "12: liutex magnitude vector y-direction"
    write(*,*) "13: liutex magnitude vector z-direction"
    write(*,*) "14: modified omega liutex vector x-direction"
    write(*,*) "15: modified omega liutex vector y-direction"
    write(*,*) "16: modified omega liutex vector z-direction"
    write(*,*) "17: modified omega liutex magnitude"
    write(*,*) "18: modified omega liutex magnitude gradient vector x-direction"
    write(*,*) "19: modified omega liutex magnitude gradient vector y-direction"
    write(*,*) "20: modified omega liutex magnitude gradient vector z-direction"
    write(*,*) " "

    allocate(f2(nx, ny, n_levels, n_data_variables))

    f2(:,:,:,1)  = lat
    f2(:,:,:,2)  = long
    f2(:,:,:,3)  = pressure
    f2(:,:,:,4)  = u
    f2(:,:,:,5)  = v
    f2(:,:,:,6)  = w
    f2(:,:,:,7)  = liutex_vector(:,:,:,1)
    f2(:,:,:,8)  = liutex_vector(:,:,:,2)
    f2(:,:,:,9)  = liutex_vector(:,:,:,3)
    f2(:,:,:,10) = liutex_magnitude
    f2(:,:,:,11) = liutex_magnitude_gradient(:,:,:,1)
    f2(:,:,:,12) = liutex_magnitude_gradient(:,:,:,2)
    f2(:,:,:,13) = liutex_magnitude_gradient(:,:,:,3)
    f2(:,:,:,14) = mod_omega_liutex_vec(:,:,:,1)
    f2(:,:,:,15) = mod_omega_liutex_vec(:,:,:,2)
    f2(:,:,:,16) = mod_omega_liutex_vec(:,:,:,3)
    f2(:,:,:,17) = mod_omega_liutex_mag
    f2(:,:,:,18) = mod_omega_liutex_mag_gradient(:,:,:,1)
    f2(:,:,:,19) = mod_omega_liutex_mag_gradient(:,:,:,2)
    f2(:,:,:,20) = mod_omega_liutex_mag_gradient(:,:,:,3)

    !! Replacing all zeros with a more OpenGrads friendly number.

    undefined_value = 9.999e20
    
    do r = 4, n_data_variables
        do k = 1, n_levels
            do j = 1, ny
                do i = 1, nx

                    if (f2(i,j,k,r) == 0.0) then
                        f2(i,j,k,r) = undefined_value
                    end if

                end do 
            end do
        end do
    end do

    !! Writing .dat file.

    open(file4, file=trim(data_output_filename), form='unformatted', action='write')
    write(file4) nx, ny, n_levels, n_data_variables

    write(file4) ((((f2(i, j, k, r), i=1,nx), j=1,ny), k=1,n_levels), r=1,n_data_variables)

    close(file4)
    
    write(*,*) "File writing successful."
    
    deallocate(f2)
    deallocate(lat, long, pressure)
    deallocate(u, v, w)
    deallocate(liutex_vector, liutex_magnitude)
    deallocate(mod_omega_liutex_vec, mod_omega_liutex_mag)
    write(*,*) "Program finished."

end program main


function latlong_dist(lat1, long1, lat2, long2) result(distance)
    !!! Finds the distance between two points (latitude, longitude) on Earth in meters.
    !!! Oscar Alvarez
    implicit none

    real, parameter :: pi = 4.0 * atan(1.0)
    real, parameter :: earth_radius_m = 6378137.0    !! Radius of the Earth in meters
    
    real, intent(in) :: lat1, long1, lat2, long2
    real :: distance

    real :: lat1_rad, lat2_rad, long1_rad, long2_rad
    real :: d_lat, d_long, a, c

    !! Convert latitude and longitude from degrees to radians
    lat1_rad     = lat1 * pi / 180.0
    lat2_rad     = lat2 * pi / 180.0
    long1_rad    = long1 * pi / 180.0
    long2_rad    = long2 * pi / 180.0

    !! Distance between each latitude and longtitude
    d_lat   = lat2_rad - lat1_rad
    d_long  = long2_rad - long1_rad
    
    !! Equations that will calculate the distance between the points in meters.
    a = sin(d_lat/2.0) * sin(d_lat/2.0) + cos(lat1_rad) * cos(lat2_rad) * sin(d_long/2.0) * sin(d_long/2.0)
    c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))

    !! Finish equation
    distance = earth_radius_m * c

end function latlong_dist