! 07/18/2025
! Ensheng Weng, read netcdf files for BiomeE global run
!     Adapted from NCAR netCDF package.
!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!==============================================================
!==============================================================
module netcdf_io
  use netcdf
  use datatypes
  implicit none

  private

  public ReadNCfiles, CRU_Interpolation, CRU_end
  public read_GridLonLat, read_interpolatedCRU
  public unzip_gzip_file

contains

!===================================================
subroutine ReadNCfiles (fpath,fields,yr_start, yr_end)
    implicit none
    character(len=*), intent(in) :: fpath
    character(len=*), intent(in) :: fields(:)
    integer, intent(in) :: yr_start, yr_end

    !-------- local vars -----------------
    character(len=256) :: command,fname,fnc
    character(len=3)   :: PFTID(9)
    character(len=4)   :: yr_str
    character(len=6)   :: GridStr
    integer, pointer   :: GridMask(:,:) => null()
    integer :: N_yrs,totL,N_vars
    integer :: istat1,i,j,k,m,iLon,iLat
    real :: dataarray(Nlon,Nlat,Ntime),timearray(Ntime)
    real :: PFTdata(144,90,9),VegFraction(144,90)

    ! Read in PFT map
    !PFTID = [character(len=3) :: 'TEB','EGN','CDB','TDB','CDN','CAS','AAS','C3G','C4G']
    PFTID = [character(len=3) :: 'C4G','C3G','TEB','TDB','EGN','CDB','CDN','CAS','AAS']
    fnc = trim(fpath)//'BiomeE-PFTs.nc'
    write(*,*)trim(fnc)
    do i=1, 9
      call nc_read_2D(fnc,PFTID(i),144,90,PFTdata(:,:,i))
      write(*,*)"Map PFT: ", PFTID(i)
    enddo
    do i =1, 144
      do j=1, 90
        VegFraction(i,j) = max(0.0,min(1.0, sum(PFTdata(i,j,:))))
      enddo
    enddo

    ! -------------- Select land grids for model run ------------------!
    write(yr_str, '(I4)') yr_start
    fnc = trim(fpath)//trim(fields(1))//'/'//trim(ncversion)//trim(fields(1))//'.'//trim(yr_str)//'.365d.noc.nc'

#ifdef ZippedNCfiles
    call unzip_gzip_file(trim(fnc)//'.gz')
#endif

    call nc_read_3D(fnc,trim(fields(1)),Nlon,Nlat,Ntime,dataarray)

#ifdef ZippedNCfiles
    command = 'rm '//trim(fnc) ! Remove unziped file
    call execute_command_line(command)
#endif

    ! Tag the vegetated grids for model run
    allocate(GridMask(LowerLon:UpperLon, LowerLat:UpperLat))
    GridMask(:,:) = 0
    m = 0
    do iLon = LowerLon, UpperLon, StepLatLon
      do iLat = LowerLat, UpperLat, StepLatLon
        i = iLon/5 + 1
        j = iLat/4 + 1
        if(dataarray(ilon,ilat,1) < 9999.0 .and. VegFraction(i,j)>0.1)then
          GridMask(ilon,ilat) = 1  ! Select the grids for model run
          m = m + 1
        endif
      enddo
    enddo
    N_VegGrids = m
    write(*,*)"Valid grids: ", N_VegGrids
    allocate(GridLonLat(N_VegGrids))
    grid_No1 = 1
    grid_No2 = N_VegGrids ! ! Run all the grids in GridLonLat

    ! Put GridID into an array
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0)then
          m = m + 1
          GridLonLat(m) = iLon * 1000 + iLat
        endif
      enddo
    enddo

#ifndef Use_InterpolatedData
    !fields = [character(len=5) :: 'tmp', 'pre', 'dswrf', 'spfh']
    N_vars = size(fields)
    N_yrs  = yr_end - yr_start + 1
    totL   = N_yrs * Ntime
    ! Allocate data arrays
    !allocate(CRUData(totL, 4, LowerLon:UpperLon, LowerLat:UpperLat))
    allocate(CRUtime(totL))
    allocate(ClimData(totL, N_vars, N_VegGrids))
    allocate(LandGrid(N_VegGrids))

    ! Record GridStr and PFT coverage for each grid
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0)then
          m = m + 1
          LandGrid(m)%iLon = iLon
          LandGrid(m)%iLat = iLat
          LandGrid(m)%climate => ClimData(:,:,m)
          ! Assign PFT coverage for each grid
          i = (iLon-1)/5 + 1
          j = (iLat-1)/4 + 1
          do k = 1, 9
            LandGrid(m)%fPFT(k) = Min(1.0, Max(0.0, PFTdata(i,j,k)))
          enddo
        endif
      enddo
    enddo

    ! Write GridLonLat and forcing file names to a file
    if(WriteForcing)then
      fname = trim(filepath_out)//trim(GridListFile) ! List file name
      open(NEWUNIT=Grids_Unit,file=trim(fname),ACTION='write', IOSTAT=istat1)
    endif

    ! ----------------- Read in all data ----------------------!
    do j= 1, N_vars ! 4 or 7 (including pres, ugrd, vgrd)
      do i =1, N_yrs
        write(yr_str, '(I4)') yr_start + i - 1
        fnc = trim(fpath)//trim(fields(j))//'/'//trim(ncversion)//trim(fields(j))//'.'//trim(yr_str)//'.365d.noc.nc'

#ifdef ZippedNCfiles
        call unzip_gzip_file(trim(fnc)//'.gz') ! Unzip the nc data file
#endif

        write(*,*)'Reading: ', trim(fnc)
        call nc_read_3D(fnc,trim(fields(j)),Nlon,Nlat,Ntime,dataarray)
        m = 0
        do iLon = LowerLon, UpperLon
          do iLat = LowerLat, UpperLat
            !CRUData((i-1)*Ntime+1:i*Ntime,j,iLon,iLat) = dataarray(iLon,iLat,:)
            if(GridMask(iLon,iLat) > 0)then
              m = m + 1
              ClimData((i-1)*Ntime+1:i*Ntime,j,m) = dataarray(iLon,iLat,:)
            endif
          enddo
        enddo

        ! Read in the time array of the first variable
        if(j == 1)then
          call nc_read_1D(fnc,'time',Ntime,timearray)
          CRUtime((i-1)*Ntime+1: i*Ntime) = timearray
        endif

#ifdef ZippedNCfiles
        command = 'rm '//trim(fnc) ! Remove unziped nc data file
        call execute_command_line(command)
#endif
      enddo ! N_yrs
    enddo   ! four variables
#endif
    ! Release allocatable arrays
    deallocate(GridMask)
end subroutine ReadNCfiles

!=============================================================================
subroutine CRU_Interpolation(LandGrid,forcingData)
  implicit none
  type(grid_initial_type), pointer, intent(in) :: LandGrid
  type(climate_data_type), pointer :: forcingData(:) ! output

  !---------- local variables ------------------
  type(climate_data_type), pointer :: climateData(:) ! will be pointed by forcingData
  character(len=256):: command, fname   ! For testing output
  character(len=6)  :: GridStr  ! Used in output file name
  real, pointer     :: GridData(:,:)
  real, allocatable :: fdSW(:)
  real, allocatable :: timecols(:,:)
  real, allocatable :: hourly_data(:,:)
  real    :: Lati, Longi
  real    :: steps_in_6H ! Temporary variable, steps interpolated
  real    :: td,cosz,solarelev,solarzen,r_light
  real    :: WindS1, WindS2, SWdaily, SWmax
  real    :: tmp1(12,10), tmp2(SHshift,10)  ! Shift hourly data
  integer :: forcing_unit ! for interpolated grid forcing file writting
  integer :: iLon, iLat   ! Column and Lines (started from -179.75 and -89.75)
  integer :: year0, year1 ! Start and end year
  integer :: yr,doy,iday,ihour,iyr
  integer :: ndays,nyear,totalL,Nsteps
  integer :: iostat,m,i,j,k
  integer :: totyr, totDays, Nlines

  ! Assigne LandGrid data to local variables
  GridData => LandGrid%climate
  iLon = LandGrid%iLon
  iLat = LandGrid%iLat

  ! Latitude and Longitude of this grid
  Longi = Lon0 + (iLon - 0.5) * Wlon
  Lati  = Lat0 + (iLat - 0.5) * Wlat

  ! Data lines
  Nlines = SIZE(CRUtime)
  steps_per_hour= TargetSteps_per_hour
  steps_per_day = steps_per_hour * 24
  totalL = Nlines * Hours_NCstep * steps_per_hour ! for hourly data
  totyr = totalL/hours_per_year
  totDays = totyr * 365
  year0 = int(CRUtime(1)/365) + 1901 ! Started from 1901/1/1
  year1 = year0 + totyr - 1

  ! Allocate allocatable variables
  allocate(fdSW(Nlines))
  allocate(timecols(Nlines,3)) ! year, doy, hour
  allocate(hourly_data(totalL,10)) ! Interpolated data (hourly)
  allocate(climateData(totalL))
  timecols(:,:) = 0
  do m=1, Nlines
    timecols(m,1) = int(CRUtime(m)/365.0) + 1901       ! Year
    timecols(m,2) = MOD(int(CRUtime(m)),365) + 1.0    ! Day of the year
    timecols(m,3) = MOD(CRUtime(m)*24.0,24.0)+int(Longi/15.0) ! Local time
    timecols(m,3) = MOD(timecols(m,3), 24.0) ! Converted to 0~23
  enddo

  ! Calculate fdSW
  Nsteps = 24/Hours_NCstep ! netCDF Data daily steps, 4 for CRU
  do i = 1, totDays
    ! Actual daily radiation
    SWdaily = sum(GridData((i-1)*Nsteps+1 : i*Nsteps, 3))

    ! Calculate daily max SW (SWmax)
    SWmax = 0.0
    do j=1,96 ! 15 minutes
      td = timecols(m,2) + (j-1)/96.0
      call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
      SWmax = SWmax + cosz * solarC * seconds_per_day/96.0
    enddo

    ! Fraction of solar radiation
    fdSW((i-1)*Nsteps+1 : i*Nsteps) = Max(0.0,Min(1.0,SWdaily/(SWmax+0.0001))) *0.8 ! too high!
  enddo

  ! ------------- Data interpolation ------------------
  ! Linear interpolation (and update radiation and PAR later)
  steps_in_6H = Hours_NCstep * steps_per_hour
  m = 0
  do i = 1, Nlines -1
    do j=0, int(steps_in_6H-1.0)
      m=m+1
      hourly_data(m,4) = GridData(i,1) + (GridData(i+1,1)-GridData(i,1))*j/steps_in_6H ! Tmp, K
      hourly_data(m,5) = GridData(i,2)/(6.0 * 3600.0)                  ! Precipitation, mm/second
      hourly_data(m,7) = GridData(i,4) + (GridData(i+1,4)-GridData(i,4))*j/steps_in_6H ! spfh
      hourly_data(m,10)= fdSW(i) + (fdSW(i+1)-fdSW(i))*j/steps_in_6H          ! fraction of Shortwave radiation
      WindS1 = SQRT(GridData(i,  6)**2 + GridData(i,  7)**2)
      WindS2 = SQRT(GridData(i+1,6)**2 + GridData(i+1,7)**2)
      hourly_data(m,8) = GridData(i,5) + (GridData(i+1,5)-GridData(i,5))*j/steps_in_6H ! air presssure, Pa
      hourly_data(m,9) = WindS1 + (WindS2 - WindS1)*j/steps_in_6H
    enddo
  enddo
  ! Last 6 hours' data
  do j=0, int(steps_in_6H-1.0)
    hourly_data(totalL-j,4) = GridData(Nlines,1) + (GridData(1,1)-GridData(Nlines,1))*j/steps_in_6H ! Tmp, K
    hourly_data(totalL-j,5) = GridData(Nlines,2)/(6.0 * 3600.0)                    ! Precipitation, mm/second
    hourly_data(totalL-j,7) = GridData(Nlines,4) + (GridData(1,4)-GridData(Nlines,4))*j/steps_in_6H ! spfh
    hourly_data(totalL-j,10)= fdSW(Nlines) + (fdSW(1)-fdSW(Nlines))*j/steps_in_6H
    WindS1 = SQRT(GridData(Nlines,6)**2 + GridData(Nlines,7)**2)
    WindS2 = SQRT(GridData(1,     6)**2 + GridData(1,     7)**2)
    hourly_data(totalL-j,8) = GridData(Nlines,5) + (GridData(1,5)-GridData(Nlines,5))*j/steps_in_6H ! air presssure, Pa
    hourly_data(totalL-j,9) = WindS1 + (WindS2 - WindS1)*j/steps_in_6H
  enddo
  ! Assign air pressure and wind speed since they are not read in from NC files
  !hourly_data(:,8) = 101325.0 ! air presssure, Pa
  !hourly_data(:,9) = 1.2      ! Wind speed m/s

  ! Shift data to fit a whole day
  m = int(timecols(1,3))
  if(m >= 1) then ! Move down
    tmp1(1:m,:) = hourly_data(totalL-m+1:totalL,:)
    hourly_data(m+1:totalL,:) = hourly_data(1:totalL-m,:)
    hourly_data(1:m,:) = tmp1(1:m,:)
  endif

  ! Add time columns
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

  ! Update radiation according to solar constant, zenith angle and SW fraction
  do i=1,totalL
    td = hourly_data(i,2) + hourly_data(i,3)/24.0
    call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
    hourly_data(i,6) = solarC * cosZ * hourly_data(i,10)      ! W/m2
  enddo

  ! Shift up 182 days for Southern Hemisphere, Weng 2025-09-22
  if(ShiftSHdata .and. Lati < 0.0) then
    tmp2(1:SHshift,4:10) = hourly_data(1:SHshift,4:10)
    hourly_data(1:totalL-SHshift,4:10) = hourly_data(SHshift+1:totalL,4:10)
    hourly_data(totalL-SHshift+1:totalL,4:10) = tmp2(1:SHshift,4:10)
  endif

  ! -------------- Put the data into forcing -------------------
  do i=1,totalL
     !'tmp','pre','dswrf','spfh','pres','windU'
     climateData(i)%year  = int(hourly_data(i,1))     ! Year
     climateData(i)%doy   = int(hourly_data(i,2))     ! day of the year
     climateData(i)%hod   = hourly_data(i,3)          ! hour of the day
     climateData(i)%Tair  = hourly_data(i,4)          ! air temperature, K
     climateData(i)%Tsoil = (climateData(i)%Tair - 273.16) * 0.8 + 273.16 ! soil temperature, C
     climateData(i)%rain  = hourly_data(i,5)          ! kgH2O m-2 second-1
     climateData(i)%radiation = hourly_data(i,6)      ! W/m2
     climateData(i)%PAR       = hourly_data(i,6) * 2.0 ! umol/m2/s
     climateData(i)%RH    = Max(0.01, min(0.99, hourly_data(i,7) / mol_h2o * mol_air * &
                            hourly_data(i,8) / esat(hourly_data(i,4) - 273.16))) ! relative humidity (0.xx)
     climateData(i)%P_air = hourly_data(i,8)          ! pa
     climateData(i)%windU = hourly_data(i,9)          ! wind velocity (m s-1)
     climateData(i)%CO2   = CO2_c ! CO2_Hist(climateData(i)%year-1700)       ! ppm
     climateData(i)%eCO2  = climateData(i)%CO2 + 200.       ! ppm
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  datalines = totalL
  days_data = totyr * 365
  yr_data   = totyr
  dt_fast_yr= 1.0/(365.0 * 24.0 * steps_per_hour)
  step_hour = 1.0/steps_per_hour
  step_seconds  = 3600.0 * step_hour

  ! Write out a sample file
  if(WriteForcing)then
    write(GridStr, GridIDFMT) GridID
    fname = trim(ncversion)//trim(GridStr)//'_forcing.csv' ! Data file name
    write(Grids_Unit, '(a6,"," a35)')trim(GridStr),trim(fname)
    fname = trim(filepath_out)//trim(fname)
    open(NEWUNIT=forcing_unit,file=trim(fname))
    !write(15,*)"YEAR,DOY,HOUR,PAR,Swdown,Tair,Tsoil,RH,RAIN,WIND,PRESSURE,aCO2,eCO2"
    write(forcing_unit,*)'Swdown,Tair,RH,RAIN,WIND,PRESSURE'
    do i=1,totalL
        write(forcing_unit,'(6(E15.8,","))') &
          climateData(i)%radiation,climateData(i)%Tair,climateData(i)%RH,  &
          climateData(i)%rain,climateData(i)%windU,climateData(i)%P_air
    enddo
    close(forcing_unit)
    command = 'gzip -f ' // trim(fname)
    call execute_command_line(command, exitstat=iostat)
  endif

  !Release memory
  deallocate(hourly_data, fdSW, timecols)
end subroutine CRU_Interpolation

!==============================================
subroutine CRU_end()
  close(Grids_Unit)
  deallocate(GridLonLat)
#ifndef Use_InterpolatedData
  !deallocate(CRUData)
  deallocate(CRUtime)
  deallocate(ClimData)
  deallocate(LandGrid)
#endif
end subroutine CRU_end

!=============================================================================
! Read the interpolated data file list
subroutine read_GridLonLat(fname,file_exists)
  implicit none
  character(len=*),intent(in) :: fname
  logical, intent(inout) :: file_exists

  ! ------- Local vars ---------------
  integer :: GridNo(Nlon*Nlat),tmpNo(Nlon*Nlat) ! maximum grids, 720*360
  integer :: i,j,k,m,n,istat1
  character(len=300) :: listfile

  listfile=trim(int_fpath)//trim(fname)
  INQUIRE (file=trim(listfile), EXIST=file_exists)
  if (.not. file_exists) then
    write (*, '("read_GridLonLat: ", a, " does not exist")') trim(listfile)
    return
  endif

  ! Read the file
  open(11,file=listfile,status='old',ACTION='read',IOSTAT=istat1)
  m = 0
  do
    read(11,*,IOSTAT=istat1) GridNo(m+1)
    if(istat1 < 0)exit
    m = m + 1
  enddo

  ! Update GridNo and m with StepLatLon
  if(StepLatLon > 1 .or. UpperLon < 720)then ! Regional or partial run
    N_VegGrids = m
    tmpNo = GridNo
    m = 0
    n = 1
    do i = LowerLat, UpperLon, StepLatLon
      do j = LowerLon, UpperLat, StepLatLon
        k = i*1000 + j
        if (k > tmpNo(n)) then
          do while(k > tmpNo(n) .and. n < N_VegGrids)
            n = n + 1
          enddo
        elseif (k < tmpNo(n)) then
          cycle
        endif
        if (k == tmpNo(n))then
            m = m + 1
            GridNo(m) = tmpNo(n)
            n = n + 1
        endif
        if(n>=N_VegGrids)exit
      enddo
    enddo
  endif

  ! Update N_VegGrids and GridLonLat
  N_VegGrids = m
  allocate(GridLonLat(N_VegGrids))
  GridLonLat(:) = GridNo(1:N_VegGrids)
  grid_No1 = min(grid_No1,N_VegGrids)
  grid_No2 = min(grid_No2,N_VegGrids)

  write(*,*)"Read GridLonLat", N_VegGrids, grid_No1, grid_No2
end subroutine read_GridLonLat

!=============================================================================
subroutine read_interpolatedCRU(fpath,fprefix,GridID,year0,year1,forcingData,file_exists)
  character(len=*),intent(in) :: fpath,fprefix
  integer, intent(in) :: GridID, year0, year1
  type(climate_data_type), pointer :: forcingData(:)
  logical, intent(out) :: file_exists

  !------------local var -------------------
  type(climate_data_type), pointer :: climateData(:)
  character(len=250) :: command, climfile, fname
  character(len=6)  :: commts, GridStr
  integer, parameter :: niterms = 6 ! 6 columns in the interpolated files

  real,allocatable :: timecols(:,:), input_data(:,:)
  real    :: temp(niterms)
  integer :: istat1,istat2
  integer :: m,n,i,iyr,iday,ihour

  ! Findout the data file
  write(GridStr,GridIDFMT) GridID
  fname = trim(fprefix)//trim(GridStr)//'_forcing.csv'
  climfile=trim(fpath)//trim(fname)
  INQUIRE (file=trim(climfile)//'.gz', EXIST=file_exists)
  if (file_exists) then
    call unzip_gzip_file(trim(climfile)//'.gz')
  else
    forcingData => NULL()
    write (*, '("read_interpolatedCRU: ", a, " does not exist")') trim(climfile)
    return
  endif

  ! Total days and lines
  steps_per_day = 24 ! Make sure to be consistent with interpolation
  yr_data   = year1 - year0 + 1
  days_data = yr_data * 365
  datalines = days_data * steps_per_day

  ! Allocate arrays for reading in data
  allocate(input_data(niterms,datalines))
  allocate(timecols(3,datalines))

  ! Read in forcing data
  open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
  read(11,'(a160)',IOSTAT=istat2) commts ! One line comments
  m = 0
  do while (m < datalines .and. istat2 == 0) ! Only read in maximum of datalines
    read(11,*,IOSTAT=istat2)(temp(n), n = 1,niterms)
    m = m + 1
    input_data(:,m) = temp(:)
  enddo ! end of reading the forcing file

  ! Close the file and delete it
  close(11)    ! close forcing file
  command = 'rm '//trim(climfile) ! Remove unziped file
  call execute_command_line(command)

  ! Check the consistency between the file data lines and required
  if(m /= datalines) then
    write (*, '("In read_interpolatedCRU, File ",a," is shorter than needed: lines: ",I12)') trim(fname),m
    file_exists = .False.
    deallocate(input_data,timecols)
    return
  endif

  ! Setup the time table
  m=0
  do iyr = year0, year1
    do iday=1, 365
      do ihour = 1, steps_per_day ! 24
        m =  m + 1
        timecols(1,m) = iyr
        timecols(2,m) = iday
        timecols(3,m) = 24.0 * (ihour - 1)/steps_per_day
      enddo
    enddo
  enddo

  ! Put the data into forcing
  allocate(climateData(datalines))
  do i=1,datalines
     climateData(i)%year      = int(timecols(1,i))         ! Year
     climateData(i)%doy       = int(timecols(2,i))         ! day of the year
     climateData(i)%PAR       = input_data(1,i)*2.0        ! umol/m2/s
     climateData(i)%radiation = input_data(1,i)            ! W/m2
     climateData(i)%Tair      = input_data(2,i)  ! air temperature, K
     climateData(i)%Tsoil     = input_data(2,i)*0.8 + 273.16*0.2  ! soil temperature, K
     climateData(i)%RH        = Max(0.01, min(0.99, input_data(3,i) ))       ! relative humidity (0.xx)
     climateData(i)%rain      = input_data(4,i)        ! kgH2O m-2 s-1
     climateData(i)%windU     = input_data(5,i)        ! wind velocity (m s-1)
     climateData(i)%P_air     = input_data(6,i)        ! pa
     climateData(i)%CO2       = CO2_c        !ppm
     climateData(i)%eCO2      = CO2_c + 200. !ppm
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  write(*,*)"forcing from interpolated: hours,days,years", datalines,days_data,yr_data

  !Close opened file and release memory
  deallocate(input_data,timecols)
end subroutine read_interpolatedCRU

!==============================================================
  subroutine unzip_gzip_file(filename_gz)
    character (len = *), intent(in) :: filename_gz
    !----------local vars
    character(len=256) :: command
    integer :: iostat

    ! Construct the gunzip command. The -k option keeps the original .gz file.
    command = 'gunzip -k ' // trim(filename_gz)
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

    !----- Local vars ----------------
    integer :: ncid, varid  ! IDs were created with netCDF files

    ! Open the file with NF90_NOWRITE as read-only access
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
    call check( nf90_get_var(ncid, varid, DA) )
    call check( nf90_close(ncid) )
    print *, 'Read file: ncid=',ncid, 'varid=',varid

  end subroutine nc_read_3D

  !==============================================
  subroutine nc_read_2D(FILE_NAME,field_idx,NX,NY,DA)
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: NX, NY
    real, intent(inout) :: DA(:,:)

    !----- Local vars ----------------
    integer :: ncid, varid ! IDs were created with netCDF files

    ! Open the file with NF90_NOWRITE as read-only access
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) ) ! Get the varid of the data variable
    call check( nf90_get_var(ncid, varid, DA) )  ! Read the data.
    call check( nf90_close(ncid) ) ! Close the file, freeing all resources.
    print *, 'Read file: ncid=',ncid, 'varid=',varid
  end subroutine nc_read_2D
  !==============================================
  subroutine nc_read_1D(FILE_NAME,field_idx,Ntime,DA)
    ! This is the name of the data file we will read
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: Ntime
    real, intent(inout) :: DA(:)

    !----- Local vars ----------------
    integer :: ncid, varid ! IDs were created with netCDF files

    ! Open the file with NF90_NOWRITE as read-only access
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
    call check( nf90_get_var(ncid, varid, DA) )
    call check( nf90_close(ncid) )
    print *, 'Read file: ncid=',ncid, 'varid=',varid
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
