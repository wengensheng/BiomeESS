! 07/18/2025
! Ensheng Weng, read netcdf files for BiomeE global run

!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is a very simple example which writes a 2D array of
!     sample data. To handle this in netCDF we create two shared
!     dimensions, "x" and "y", and a netCDF variable, called "data".

!     This example demonstrates the netCDF Fortran 90 API. This is part
!     of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_xy_wr.f90,v 1.7 2006/12/09 18:44:58 russ Exp $

!==============================================================
!==============================================================
module netcdf_io
  use netcdf
  use datatypes

  ! CRU NetCDF file dimensions
  integer, parameter :: NDIMS = 3
  integer, parameter :: Nlon = 720, Nlat = 360, Ntime = 1460
  real, pointer :: tswrfH(:) ! Hours of tswrf (hours since 1850-01-01)
  real, pointer :: CRUData(:,:,:,:) ! N_yr*Ntime, N_vars, Nlon, Nlat
  real, pointer :: GridClimateData(:,:,:)  ! N_yr*Ntime, N_vars, N_VegGrids, for land grids
  integer, pointer :: GridLonLat(:) ! iLon, iLat
  logical :: WriteSample = .True.
  integer :: N_VegGrids

  public ReadNCfiles
  public unzip_gzip_file

contains
!===================================================
  subroutine ReadNCfiles (fpath,fields,yr_start, yr_end)
    implicit none
    character(len=*), intent(in) :: fpath
    character(len=*), intent(in) :: fields(:)
    integer, intent(in) :: yr_start, yr_end

    !-------- local vars -----------------
    character(len=256) :: command
    character (len =256) :: fname,fnc,fgz
    character (len =20) :: field_idx
    character (len =3) :: PFTID(9)
    character (len =4) :: yr_str

    integer, pointer :: GridMask(:,:)
    integer :: N_yrs,totL,N_vars
    integer :: iostat,i,j,k,m,iLon,iLat
    real :: dataarray(Nlon,Nlat,Ntime),timearray(Ntime)
    real :: PFTdata(144,90,9),VegFraction(144,90)

    ! Read in PFT map
    PFTID = [character(len=3) :: 'TEB','EGN','CDB','TDB','CDN','CAS','AAS','C3G','C4G']
    fnc = trim(fpath)//'BiomeE-PFTs.nc'
    write(*,*)trim(fnc)
    do i=1, 9
      write(*,*)"read PFT: ", PFTID(i)
      call nc_read_2D(fnc,PFTID(i),144,90,PFTdata(:,:,i))
      VegFraction(:,:) = VegFraction(:,:) + PFTdata(:,:,i)
    enddo
    do i =1, 144
      do j=1, 90
        VegFraction(i,j) = max(0.0,min(1.0, VegFraction(i,j)))
      enddo
    enddo
    !write(*,'(360(E10.4,","))')  VegFraction(:,70)
    !write(*,'(360(f22.2,","))')  PFTdata(:,70,3)

    !fields = [character(len=5) :: 'tmp', 'pre', 'tswrf', 'spfh']
    N_vars = size(fields)
    N_yrs = yr_end - yr_start + 1
    totL = N_yrs*Ntime

    ! Open a file for identifying valid grids
    fgz = trim(fpath)//'tswrf_v12_2010.nc.gz'
    fnc = trim(fpath)//'tswrf_v12_2010.nc'
    write(*,*)trim(fgz)

    call unzip_gzip_file(trim(fgz))
    call nc_read_3D(fnc,'tswrf',Nlon,Nlat,Ntime,dataarray)
    ! Remove unziped file
    command = 'rm ' // trim(fnc)
    call execute_command_line(command)

    ! Calculate number of vegetated grids
    allocate(GridMask(LowerLon:UpperLon, LowerLat:UpperLat))
    GridMask = 0
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        i = iLon/5 + 1
        j = iLat/4 + 1
        if(dataarray(ilon,ilat,1) < 9999.0 .and. VegFraction(i,j)>0.1)then
          m = m + 1
          GridMask(ilon,ilat) = 1
        endif
      enddo
    enddo
    N_VegGrids = m
    write(*,*)"Valid grids: ",N_VegGrids

    ! Allocate data arrays
    allocate(tswrfH(totL))
    !allocate(CRUData(totL, 4, LowerLon:UpperLon, LowerLat:UpperLat))
    allocate(GridClimateData(totL, N_vars, N_VegGrids))
    allocate(GridLonLat(N_VegGrids))

    ! Record LonLat in GridLonLat
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0)then
          m = m + 1
          GridLonLat(m) = iLon * 1000 + iLat
        endif
      enddo
    enddo
    write(*,*)"Re-counted: ", sum(GridMask), m

    ! Read in all data
    do j= 1, 4
      field_idx = fields(j) ! 'tmp'
      do i =1, N_yrs
        write(yr_str, '(I4)') yr_start + i - 1
        if(field_idx == 'tswrf')then
          fnc = trim(field_idx)//'_v12_'//trim(yr_str)//'.nc'
        else
          fnc = 'crujra.v2.4.5d.'//trim(field_idx)//'.'//trim(yr_str)//'.365d.noc.nc'
        endif
        fgz = trim(fnc)//'.gz'

        ! Unzip a file
        fgz = trim(fpath)//trim(fgz)
        call unzip_gzip_file(trim(fgz))

        fnc = trim(fpath)//trim(fnc)

        write(*,*)'read nc file:', fnc
        call nc_read_3D(fnc,trim(field_idx),Nlon,Nlat,Ntime,dataarray)
        m = 0
        do iLon = LowerLon, UpperLon
          do iLat = LowerLat, UpperLat
            !CRUData((i-1)*Ntime+1:i*Ntime,j,iLon,iLat) = dataarray(iLon,iLat,:)
            if(GridMask(iLon,iLat) > 0)then
              m = m + 1
              GridClimateData((i-1)*Ntime+1:i*Ntime,j,m)  = dataarray(iLon,iLat,:)
            endif
          enddo
        enddo

        if(field_idx == 'tswrf')then
          call nc_read_1D(fnc,'time',Ntime,timearray)
          tswrfH((i-1)*Ntime+1: i*Ntime) = timearray
          !write(*,*)"Time Data"
          !write(*,'(360(E15.4,","))')  timearray(1:360)
        endif

        ! Remove unziped file
        command = 'rm ' // trim(fnc)
        call execute_command_line(command)
      enddo
    enddo
    deallocate(GridMask)
  end subroutine ReadNCfiles

!=============================================================================
subroutine CRU_Interpolation(GridData,iLon,iLat,steps_per_hour,forcingData)
  implicit none
  real, intent(in) :: GridData(:,:), steps_per_hour
  integer, intent(in) :: iLon,iLat ! Column and Lines (started from -179.75 and -89.75)
  type(climate_data_type), pointer :: forcingData(:) ! output

  ! Define local variables
  type(climate_data_type), pointer :: climateData(:)
  character(len=120) :: fname3
  real, allocatable :: fdSW(:)
  real, allocatable :: timecols(:,:)
  real, allocatable :: hourly_data(:,:)
  real    :: ShiftData(12,10) ! To make the first hour is hour 0
  real    :: LT ! Local time
  real    :: Lati, Longi
  real    :: steps_in_6H ! Temporary variable, steps interpolated
  real    :: td,cosz,solarelev,solarzen,r_light
  integer :: year0, year1 ! Start and end year
  integer :: yr,doy,iday,ihour,iyr
  integer :: ndays,nyear,totalL
  integer :: m,n,i,j,k,steps_per_day
  integer :: totyr, Nlines

  ! Latitude and Longitude of this grid
  Longi = (iLon-1)*0.5 - 179.75
  Lati = (iLat-1) *0.5 - 89.75
  ! Data lines
  Nlines = SIZE(tswrfH)
  steps_per_day = int(steps_per_hour * 24)
  totalL = Nlines * int(6 * steps_per_hour) ! for hourly data
  totyr = int(Nlines/(365*4))
  year0 = int(tswrfH(1)/(365*24)) + 1850
  year1 = year0 + totyr - 1

  ! Allocate allocatable variables
  allocate(fdSW(Nlines))
  allocate(timecols(Nlines,3)) ! year, doy, hour
  allocate(hourly_data(totalL,10)) ! Interpolated data (hourly)
  timecols(:,:) = 0
  do m=1, Nlines
    timecols(m,1) = int(tswrfH(m)/(365*24)) + 1850.0       ! Year
    timecols(m,2) = MODULO(int(tswrfH(m)/24),365) + 1.0    ! Day of the year
    timecols(m,3) = MODULO(tswrfH(m), 24.0) + int(Longi/15.0) ! Local time
    timecols(m,3) = MODULO(timecols(m,3), 24.0) ! 0~23
  enddo

  ! Fraction of solar radiation
  do m=1, Nlines
    td = timecols(m,2) + timecols(m,3)/24.0
    call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
    if(cosz>0.001)then
      fdSW(m) = max(min(1.0,GridData(m,3)/(cosz*1380.0)),0.0)
    else
      fdSW(m) = 0.0
    endif
  enddo

  ! ------------- Data interpolation ------------------
  m = 0
  do i = 1, Nlines -1
    do j=0, 5
      m=m+1
      steps_in_6H = 6.0*steps_per_hour
      hourly_data(m,4) = GridData(i,1) + (GridData(i+1,1)-GridData(i,1))*j/steps_in_6H ! Tmp, K
      hourly_data(m,5) = GridData(i,2)/(6.0 * 3600.0)                  ! Precipitation, mm/second
      hourly_data(m,6) = GridData(i,3) + (GridData(i+1,3)-GridData(i,3))*j/steps_in_6H ! Shortwave radiation,need to calculate zenith angle
      hourly_data(m,7) = GridData(i,4) + (GridData(i+1,4)-GridData(i,4))*j/steps_in_6H ! spfh
      hourly_data(m,8) = 101325.0 ! Pa                                       ! pres
      hourly_data(m,9) = 1.2        ! Wind speed m/s
      hourly_data(m,10)= fdSW(i) + (fdSW(i+1)-fdSW(i))*j/steps_in_6H          ! fraction of Shortwave radiation
    enddo
  enddo
  ! Last 6 hours' data
  do j=0, 5
    hourly_data(totalL-j,4) = GridData(Nlines,1) + (GridData(1,1)-GridData(Nlines,1))*j/steps_in_6H ! Tmp, K
    hourly_data(totalL-j,5) = GridData(Nlines,2)/(6.0 * 3600.0)                    ! Precipitation, mm/second
    hourly_data(totalL-j,6) = GridData(Nlines,3) + (GridData(1,3)-GridData(Nlines,3))*j/steps_in_6H ! Shortwave radiation,need to calculate zenith angle
    hourly_data(totalL-j,7) = GridData(Nlines,4) + (GridData(1,4)-GridData(Nlines,4))*j/steps_in_6H ! spfh
    hourly_data(totalL-j,8) = 101325.0 ! Pa
    hourly_data(totalL-j,9)   = 1.2 ! m/s, wind spped
    hourly_data(totalL-j,10)  = fdSW(Nlines) + (fdSW(1)-fdSW(Nlines))*j/steps_in_6H
  enddo

  ! Shift data to fit a whole day
  if(timecols(1,3)>0.9999 .and. timecols(1,3)<11.9999) then ! Move down
    m = max(1,int(timecols(1,3)))
    ShiftData(1:m,:) = hourly_data(totalL-m+1:totalL,:)
    hourly_data(m+1:totalL,:) = hourly_data(1:totalL-m,:)
    hourly_data(1:m,:) = ShiftData(1:m,:)
  elseif (timecols(1,3)>=11.9999) then ! Move up
    m = 24 - max(12,int(timecols(1,3)))
    ShiftData(1:m,:) = hourly_data(1:m,:)
    hourly_data(1:totalL-m,:) = hourly_data(m+1:totalL,:)
    hourly_data(totalL-m+1:totalL,:) = ShiftData(1:m,:)
  endif

  !Time columns
  m=0
  do iyr = year0, year1
    do iday=1, 365
      do ihour = 1, steps_per_day ! 24
        m =  m + 1
        hourly_data(m,1) = iyr
        hourly_data(m,2) = iday
        hourly_data(m,3) = 24.0 * (ihour - 1)/steps_per_day
      enddo
    enddo
  enddo

  ! -------------- Put the data into forcing -------------------
  allocate(climateData(totalL))
  do i=1,totalL
     td = hourly_data(i,2) + hourly_data(i,3)/24.0
     call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
     if(cosz>0.001)then
       climateData(i)%radiation = 1380.0 * cosZ * hourly_data(i,10)      ! W/m2
       climateData(i)%PAR       = climateData(i)%radiation * 2.0 ! umol/m2/s
     else
       climateData(i)%radiation = 0.0 ! W/m2
       climateData(i)%PAR       = 0.0 ! umol/m2/s
     endif
     !'tmp','pre','tswrf','spfh','pres','windU'
     climateData(i)%year  = int(hourly_data(i,1))     ! Year
     climateData(i)%doy   = int(hourly_data(i,2))     ! day of the year
     climateData(i)%hod   = hourly_data(i,3)          ! hour of the day
     climateData(i)%Tair  = hourly_data(i,4)          ! air temperature, K
     climateData(i)%Tsoil = (climateData(i)%Tair - 273.16) * 0.8 + 273.16 ! soil temperature, C
     climateData(i)%rain  = hourly_data(i,5)          ! kgH2O m-2 second-1
     climateData(i)%P_air = hourly_data(i,8)          ! pa
     climateData(i)%windU = hourly_data(i,9)          ! wind velocity (m s-1)
     climateData(i)%RH    = hourly_data(i,7) / mol_h2o * mol_air * &
                            climateData(i)%P_air/  &
                            esat(climateData(i)%Tair-273.16) ! relative humidity (0.xx)
     climateData(i)%CO2   = CO2_Hist(climateData(i)%year-1700)       ! ppm
     climateData(i)%eCO2  = climateData(i)%CO2 + 200.       ! ppm
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  datalines = totalL
  days_data = totyr * 365
  yr_data   = totyr
  steps_per_day = int(24.0/step_hour)
  dt_fast_yr    = step_hour/(365.0 * 24.0)
  step_seconds  = step_hour*3600.0

  ! Write out a sample file
  if(WriteSample)then
    fname3 = 'CRU_interpo_0724_forcing.csv'
    open(15,file=trim(fname3))
    write(15,*)"YEAR,DOY,HOUR,PAR,Swdown,Tair,Tsoil,RH,RAIN,WIND,PRESSURE,aCO2,eCO2"
    do i=1,totalL
        write(15,'(2(I4,","),1(f8.2,","),30(E15.6,","))') &
          climateData(i)%year, climateData(i)%doy,      climateData(i)%hod,  &
          climateData(i)%PAR,  climateData(i)%radiation,climateData(i)%Tair, &
          climateData(i)%Tsoil,climateData(i)%RH,       climateData(i)%rain, &
          climateData(i)%windU,climateData(i)%P_air, &
          climateData(i)%CO2,  climateData(i)%eCO2
    enddo
    close(15)
    WriteSample = .False.
  endif

  !Release memory
  deallocate(hourly_data)
  deallocate(fdSW)
  deallocate(timecols)

end subroutine CRU_Interpolation

!==============================================================
  subroutine unzip_gzip_file(filename_gz)
    character (len = *), intent(in) :: filename_gz
    !----------local vars
    character(len=256) :: command
    integer :: iostat

    ! Construct the gunzip command. The -k option keeps the original .gz file.
    command = 'gunzip -k ' // trim(filename_gz)

    ! Execute the command
    call execute_command_line(command, exitstat=iostat)
    if (iostat == 0) then
      print *, 'Successfully unzipped ', trim(filename_gz)
    else
      print *, 'Error unzipping ', trim(filename_gz), ' (Exit status: ', iostat, ')'
    end if

  end subroutine unzip_gzip_file

!==============================================
  subroutine nc_read_3D(FILE_NAME,field_idx,NX,NY,Ntime,DA)

  ! This is the name of the data file we will create.
  character (len = *), intent(in) :: FILE_NAME,field_idx
  integer, intent(in) :: NX, NY, Ntime
  real, intent(inout) :: DA(:,:,:)

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  integer :: ncid, varid

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
  print *, 'ncid=',ncid
  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
  print *, 'varid=',varid
  ! Read the data.
  call check( nf90_get_var(ncid, varid, DA) )

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS reading ncfile ", FILE_NAME, "! "

  end subroutine nc_read_3D

  !==============================================
  subroutine nc_read_2D(FILE_NAME,field_idx,NX,NY,DA)
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: NX, NY
    real, intent(inout) :: DA(:,:)

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    integer :: ncid, varid

    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
    ! the file.
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    print *, 'ncid=',ncid
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
    print *, 'varid=',varid
    ! Read the data.
    call check( nf90_get_var(ncid, varid, DA) )

    ! Close the file, freeing all resources.
    call check( nf90_close(ncid) )

    print *,"SUCCESS reading ncfile ", FILE_NAME, ". "



  end subroutine nc_read_2D
  !==============================================
    subroutine nc_read_1D(FILE_NAME,field_idx,Ntime,DA)

    ! This is the name of the data file we will read
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: Ntime
    real, intent(inout) :: DA(:)

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    integer :: ncid, varid

    ! Open NF90_NOWRITE with read-only access to the file.
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    print *, 'ncid=',ncid
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
    print *, 'varid=',varid
    ! Read the data.
    call check( nf90_get_var(ncid, varid, DA) )

    ! Close the file, freeing all resources.
    call check( nf90_close(ncid) )

    print *,"*** SUCCESS reading ncfile ", FILE_NAME, "! "

  end subroutine nc_read_1D

!===================================================
  subroutine nc_write(FILE_NAME,NDIMS,NX,NY)
    ! This is the name of the data file we will create.
    character (len = *), intent(in) :: FILE_NAME
    integer, intent(in) :: NDIMS, NX, NY ! , Ntime

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    integer :: ncid, varid
    integer :: x_dimid, y_dimid

    ! This is the data array we will write. It will just be filled with
    ! a progression of integers for this example.
    integer,allocatable :: data_out(:,:) !(NY, NX)
    integer,allocatable :: dimids(:)

    ! Loop indexes, and error handling.
    integer :: x, y

    ! Create some pretend data. If this wasn't an example program, we
    ! would have some real data to write, for example, model output.
    allocate(dimids(NDIMS))
    allocate(data_out(NY,NX))
    do x = 1, NX
       do y = 1, NY
          data_out(y, x) = (x - 1) * NY + (y - 1)
       end do
    end do

    ! Always check the return code of every netCDF function call. In
    ! this example program, wrapping netCDF calls with "call check()"
    ! makes sure that any return which is not equal to nf90_noerr (0)
    ! will print a netCDF error message and exit.

    ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
    ! overwrite this file, if it already exists.
    call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )

    ! Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, "x", NX, x_dimid) )
    call check( nf90_def_dim(ncid, "y", NY, y_dimid) )

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids =  (/ y_dimid, x_dimid /)

    ! Define the variable. The type of the variable in this case is
    ! NF90_INT (4-byte integer).
    call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )

    ! End define mode. This tells netCDF we are done defining metadata.
    call check( nf90_enddef(ncid) )

    ! Write the pretend data to the file. Although netCDF supports
    ! reading and writing subsets of data, in this case we write all the
    ! data in one operation.
    call check( nf90_put_var(ncid, varid, data_out) )

    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    call check( nf90_close(ncid) )
    write(*,*)data_out
    print *, "*** SUCCESS writing ncfile simple_xy.nc! "

    deallocate(data_out)
    deallocate(dimids)

  end subroutine nc_write

!===================================================
  subroutine nc_read(FILE_NAME,NDIMS,NX,NY)
    ! This is the name of the data file we will create.
    character (len = *), intent(in) :: FILE_NAME
    integer, intent(in) :: NDIMS, NX, NY

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    integer :: ncid, varid
    ! We are reading 2D data, a 6 x 12 grid.
    integer,allocatable :: data_in(:,:) !(NY, NX)

    ! Loop indexes, and error handling.
    integer :: x, y

    allocate(data_in(NY,NX))
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
    ! the file.
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    print *, 'ncid=',ncid
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, "data", varid) )
    print *, 'varid=',varid
    ! Read the data.
    call check( nf90_get_var(ncid, varid, data_in) )

    ! Check the data.
    do x = 1, NX
       do y = 1, NY
          print *, "data_in(", y, ", ", x, ") = ", data_in(y, x)
       end do
    end do

    ! Close the file, freeing all resources.
    call check( nf90_close(ncid) )

    print *,"*** SUCCESS reading ncfile ", FILE_NAME, "! "

    deallocate(data_in)
  end subroutine nc_read

!=================================================================
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end module
!=================================================================
