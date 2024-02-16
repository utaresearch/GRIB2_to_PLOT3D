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
    integer, parameter :: file1=10, file2=20, file3=30
    character(100) :: grib_grid_filename, grib_velocity_filename
    character(500) :: first_velocity_component, second_velocity_component, third_velocity_component
    character(10) :: u_string, v_string, w_string
    character(100) :: grid_output_filename, fun_output_filename
    real(8), dimension(:,:,:,:), allocatable :: f
    real(8), dimension(:,:,:), allocatable :: f_lat, f_long, f_height
    integer :: n_vars

    !! other variables
    real(8), dimension(:,:,:,:,:), allocatable :: velocity_gradient
    
    real(8), dimension(:,:,:), allocatable :: first, second, third
    real(8), dimension(:,:,:), allocatable :: u, v, w, x, y, z
    real(8), dimension(:,:,:), allocatable :: liutex_x, liutex_y, liutex_z, liutex_magnitude
    
    real(8), dimension(:,:), allocatable :: lat, long
    
    real(8), dimension(3,3) :: a

    real(8), dimension(3) :: vor, liutex_vec

    real(8) :: lat_start, lat_end, lat_step, long_start, long_end, long_step
    real(8) :: dist_x, dist_y
    real(8) :: height_in_meters
    real(8) :: velocity_tolerance
    real(8) :: liutex_mag

    integer, dimension(:), allocatable :: height
    integer :: n_levels, time_start
    integer :: nx, ny
    integer :: i, j, k, r


    !! Input file names
    grib_grid_filename = "grib_grid_data.dat"
    grib_velocity_filename = "grib_velocity_data"


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

    !! Input heights are also in hPa, we convert them to meters.
    do i = 1, n_levels
        !! Convert hPa to meters
        height_in_meters = 288150.d0*(1.d0 - (height(i) / 1013.25d0)**(1.d0/5.255d0)) / 6.5d0

        z(:,:,i) = height_in_meters
    end do

    !! Calculate the velocity gradient tensor for each point.
    write(*,*) "Calculating Liutex."

    allocate(velocity_gradient(nx, ny, n_levels, 3, 3))

    velocity_gradient = velocity_gradient_tensor(u, v, w, x, y, z, nx, ny, n_levels)

    !! Calculate Liutex for each point.
    allocate(liutex_x(nx, ny, n_levels))
    allocate(liutex_y(nx, ny, n_levels))
    allocate(liutex_z(nx, ny, n_levels))
    allocate(liutex_magnitude(nx, ny, n_levels))
    
    do k = 1, n_levels
        do j = 1, ny
            do i = 1, nx

                a = velocity_gradient(i,j,k,:,:)
                vor = vorticity(a)
                
                call liutex(a, vor, liutex_vec, liutex_mag)

                liutex_x(i,j,k) = liutex_vec(1)
                liutex_y(i,j,k) = liutex_vec(2)
                liutex_z(i,j,k) = liutex_vec(3)

                liutex_magnitude(i,j,k) = liutex_mag

            end do
        end do
    end do

    deallocate(velocity_gradient)


    !! Writing PLOT3D output files
    write(*,*) "Creating and writing PLOT3D grid (.xyz) and function (.fun) files."

    grid_output_filename = 'output_data/grib_plot3d_grid.xyz'
    fun_output_filename  = 'output_data/grib_plot3d_function.fun'

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
    n_vars = 10

    allocate(f(nx, ny, n_levels, n_vars))

    f(:,:,:,1)  = u
    f(:,:,:,2)  = v
    f(:,:,:,3)  = w
    f(:,:,:,4)  = f_lat
    f(:,:,:,5)  = f_long
    f(:,:,:,6)  = f_height
    f(:,:,:,7)  = liutex_x
    f(:,:,:,8)  = liutex_y
    f(:,:,:,9)  = liutex_z
    f(:,:,:,10) = liutex_magnitude

    open(file3, file=trim(fun_output_filename), form='unformatted', action='write')
    write(file3) nx, ny, n_levels, n_vars
    
    write(file3) ((((f(i, j, k, r), i=1,nx), j=1,ny), k=1,n_levels), r=1,n_vars)
    
    close(file3)
    
    write(*,*) "File writing successful."
    
    deallocate(f)
    deallocate(f_lat, f_long, f_height)
    deallocate(u, v, w)
    deallocate(liutex_x, liutex_y, liutex_z, liutex_magnitude)
    
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