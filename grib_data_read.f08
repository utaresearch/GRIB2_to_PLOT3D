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

    real(8), parameter :: pi = 4.d0  * atan(1.d0)

    !! External functions
    real(8), external :: latlong_dist

    !! file handling variables
    integer, parameter :: file1=10, file2=20, file3=30, file4=40
    character(100) :: filename, grib_grid_filename, grib_velocity_filename
    character(500) :: first_velocity_component, second_velocity_component, third_velocity_component
    character(10) :: u_string, v_string, w_string
    character(100) :: grid_output_filename, fun_output_filename, data_output_filename
    real(8), dimension(:,:,:,:), allocatable :: f
    real(8), dimension(:,:,:), allocatable :: f_lat, f_long, f_height
    integer :: n_vars

    !! other variables
    real(8), dimension(:,:,:,:,:), allocatable :: velocity_gradient

    real(8), dimension(:,:,:,:), allocatable :: liutex_vector, mod_omega_liutex_vec
    
    real(8), dimension(:,:,:), allocatable :: first, second, third
    real(8), dimension(:,:,:), allocatable :: u, v, w, x, y, z
    real(8), dimension(:,:,:), allocatable ::liutex_magnitude, mod_omega_liutex_mag
    
    real(8), dimension(:,:), allocatable :: lat, long
    
    real(8), dimension(3,3) :: a

    real(8) :: lat_start, lat_end, lat_step, long_start, long_end, long_step
    real(8) :: dist_x, dist_y
    real(8) :: height_in_meters
    real(8) :: velocity_tolerance
    real(8) :: liutex_mag
    real(8) :: undefined_value

    integer, dimension(:), allocatable :: height
    integer :: n_levels, time_start
    integer :: nx, ny
    integer :: i, j, k, r


    !! Input file names
    call get_command_argument(1, filename)
    grib_grid_filename = "data/"// trim(filename) //"_grid.dat"
    grib_velocity_filename = "data/"// trim(filename) //"_velocity.dat"


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
    velocity_tolerance = 1.d20

    do k = 1, n_levels
        do j = 1, ny
            do i = 1, nx

                if (u(i,j,k) > velocity_tolerance) then
                    u(i,j,k) = 0.d0
                end if
                
                if (v(i,j,k) > velocity_tolerance) then
                    v(i,j,k) = 0.d0
                end if
                
                if (w(i,j,k) > velocity_tolerance) then
                    w(i,j,k) = 0.d0
                end if

            end do
        end do
    end do
    
    !! Make latitute and longitude values.
    allocate(lat(nx, ny))
    allocate(long(nx, ny))

    lat_step = (lat_end - lat_start) / nx
    long_step = (long_end - long_start) / ny

    lat(1,:) = lat_start
    do i = 2, nx
        lat(i,:) = lat(i-1,:) + lat_step
    end do
    lat(nx,:) = lat_end

    long(:,1) = long_start
    do j = 2, ny
        long(:,j) = long(:,j-1) + long_step
    end do
    long(:,ny) = long_end

    !! Create 3D grid using meters instead of longitude and latitude
    allocate(x(nx, ny, n_levels))
    allocate(y(nx, ny, n_levels))
    allocate(z(nx, ny, n_levels))

    x(1,1,:) = 0.d0
    y(1,1,:) = 0.d0

    do j = 2, ny
        do i = 2, nx

            !! Return distance between two latitude and longitude points in meters
            dist_x = latlong_dist(lat(i-1, j), long(i, j), lat(i,j), long(i,j))
            dist_y = latlong_dist(lat(i, j), long(i, j-1), lat(i,j), long(i,j))

            x(i,j,:) = x(i-1, j, 1) + dist_x
            y(i,j,:) = y(i, j-1, 1) + dist_y
            
        end do
    end do


    !! Calculate the velocity gradient tensor for each point.
    write(*,*) "Calculating Velocity Gradient."
    allocate(velocity_gradient(nx, ny, n_levels, 3, 3))

    velocity_gradient = velocity_gradient_tensor(u, v, w, x, y, z, nx, ny, n_levels)

    !! Calculate Liutex for each point.
    write(*,*) "Calculating Liutex."
    allocate(liutex_vector(nx, ny, n_levels, 3))
    allocate(liutex_magnitude(nx, ny, n_levels))
    
    call liutex(velocity_gradient, liutex_vector, liutex_magnitude, nx, ny, n_levels)

    write(*,*) "Calculating Modified Omega Liutex."
    allocate(mod_omega_liutex_vec(nx, ny, n_levels,3))
    allocate(mod_omega_liutex_mag(nx, ny, n_levels))

    call modified_omega_liutex(velocity_gradient, mod_omega_liutex_vec, mod_omega_liutex_mag, nx, ny, n_levels)

    deallocate(velocity_gradient)


    !! Writing PLOT3D output files
    write(*,*) "Creating and writing PLOT3D grid (.xyz) and function (.fun) files."

    grid_output_filename = 'output_data/' // trim(filename) // '.xyz'
    fun_output_filename  = 'output_data/' // trim(filename) // '.fun'

    !! Writing grid file (.xyz)
    open(file2, file=trim(grid_output_filename), form='unformatted', action='write')
    write(file2) nx, ny, n_levels
    write(file2)    (((x(i, j, k), i=1,nx), j=1,ny), k=1,n_levels),  &
                    (((y(i, j, k), i=1,nx), j=1,ny), k=1,n_levels),  &
                    (((z(i, j, k), i=1,nx), j=1,ny), k=1,n_levels)

    close(file2)

    deallocate(x, y, z)

    !! Re-shaping latitude, longitude, and height arrays in order to be added to function file.
    allocate(f_lat(nx,ny,n_levels))
    allocate(f_long(nx,ny,n_levels))
    allocate(f_height(nx,ny,n_levels))
    
    do i = 1, nx
        f_lat(i,:,:) = lat(i,1)    
    end do

    do j = 1, ny
        f_long(:,j,:) = long(1,j)
    end do
    
    do k = 1, n_levels
        f_height(:,:,k) = height(k)
    end do

    deallocate(lat, long, height)

    !! Writing function file (.fun)
    n_vars = 14

    allocate(f(nx, ny, n_levels, n_vars))

    f(:,:,:,1)  = u
    f(:,:,:,2)  = v
    f(:,:,:,3)  = w
    f(:,:,:,4)  = f_lat
    f(:,:,:,5)  = f_long
    f(:,:,:,6)  = f_height
    f(:,:,:,7)  = liutex_vector(:,:,:,1)
    f(:,:,:,8)  = liutex_vector(:,:,:,2)
    f(:,:,:,9)  = liutex_vector(:,:,:,3)
    f(:,:,:,10) = liutex_magnitude
    f(:,:,:,11) = mod_omega_liutex_vec(:,:,:,1)
    f(:,:,:,12) = mod_omega_liutex_vec(:,:,:,2)
    f(:,:,:,13) = mod_omega_liutex_vec(:,:,:,3)
    f(:,:,:,14) = mod_omega_liutex_mag
    

    open(file3, file=trim(fun_output_filename), form='unformatted', action='write')
    write(file3) nx, ny, n_levels, n_vars
    
    write(file3) ((((f(i, j, k, r), i=1,nx), j=1,ny), k=1,n_levels), r=1,n_vars)
    
    close(file3)

    !! Creating data (.dat) output file
    write(*,*) " "
    write(*,*) "Creating and writing data (.dat) file."
    write(*,*) " "
    write(*,*) "(.dat) data file format/structure."
    write(*,*) "Line 1 - Dimensions (i_max, j_max, k_max, n_data_variables): "
    write(*,*) nx, ny, n_levels, n_vars
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
    write(*,*) "11: modified omega liutex vector x-direction"
    write(*,*) "12: modified omega liutex vector y-direction"
    write(*,*) "13: modified omega liutex vector z-direction"
    write(*,*) "14: modified omega liutex magnitude"
    write(*,*) " "

    data_output_filename  = 'output_data/' // trim(filename) // '_liutex_UTA.dat'

    f(:,:,:,1)  = f_lat
    f(:,:,:,2)  = f_long
    f(:,:,:,3)  = f_height
    f(:,:,:,4)  = u
    f(:,:,:,5)  = v
    f(:,:,:,6)  = w                

    !! Replacing all zeros with a more OpenGrads friendly number.

    undefined_value = 9.999d20
    
    do r = 4, n_vars
        do k = 1, n_levels
            do j = 1, ny
                do i = 1, nx

                    if (f(i,j,k,r) == 0.d0) then
                        f(i,j,k,r) = undefined_value
                    end if

                end do 
            end do
        end do
    end do

    !! Writing .dat file.

    open(file4, file=trim(data_output_filename), form='unformatted', action='write')
    write(file4) nx, ny, n_levels, n_vars

    write(file4) ((((f(i, j, k, r), i=1,nx), j=1,ny), k=1,n_levels), r=1,n_vars)

    close(file4)
    
    write(*,*) "File writing successful."
    
    deallocate(f)
    deallocate(f_lat, f_long, f_height)
    deallocate(u, v, w)
    deallocate(liutex_vector, liutex_magnitude)
    deallocate(mod_omega_liutex_vec, mod_omega_liutex_mag)
    write(*,*) "Program finished."

end program main


function latlong_dist(lat1, long1, lat2, long2) result(distance)
    !!! Finds the distance between two points (latitude, longitude) on Earth in meters.
    !!! Oscar Alvarez
    implicit none
    real(8), parameter :: pi = 4.d0 * atan(1.d0)
    real(8), parameter :: earth_radius_m = 6378137.d0    !! Radius of the Earth in meters
    
    real(8), intent(in) :: lat1, long1, lat2, long2
    real(8) :: distance

    real(8) :: lat1_rad, lat2_rad, long1_rad, long2_rad
    real(8) :: d_lat, d_long, a, c

    !! Convert latitude and longitude from degrees to radians
    lat1_rad     = lat1 * pi / 180.d0
    lat2_rad     = lat2 * pi / 180.d0
    long1_rad    = long1 * pi / 180.d0
    long2_rad    = long2 * pi / 180.d0

    !! Distance between each latitude and longtitude
    d_lat   = lat2_rad - lat1_rad
    d_long  = long2_rad - long1_rad
    
    !! Equations that will calculate the distance between the points in meters.
    a = sin(d_lat/2.d0) * sin(d_lat/2.d0) + cos(lat1_rad) * cos(lat2_rad) * sin(d_long/2.d0) * sin(d_long/2.d0)
    c = 2.d0 * atan2(sqrt(a), sqrt(1.d0 - a))

    !! Finish equation
    distance = earth_radius_m * c

end function latlong_dist