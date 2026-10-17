! Stand-alone CRU-JRA forcing interpolation program for BiomeE global runs.
! Reads the netCDF CRU climate, vegetation, and N deposition files, interpolates
! the 6-hourly CRU data to the model time step for each vegetated grid, and writes
! <ncversion><GridID>_forcing.csv.gz plus the grid list file to filepath_out.
! No BiomeE simulation is run. The subroutines are copied from netcdf_io.F90
! (ReadNCfiles, CRU_Interpolation, and the netCDF helpers) and behave as in
! GlobalDataRun.x with WriteForcing = .True. and the -DGlobalRun build.
! With Interpolate=F, no interpolation is done: the 6-hourly variables are written
! grid by grid to <ncversion><GridID>_forcing6h.csv.gz (no time shifts) so that BiomeE
! can interpolate them at run time.
!
! Self-contained (needs only a Fortran compiler and the netCDF-Fortran library):
!   gfortran CRU_Interp.F90 -o cru_interp -I<netcdf include> -L<netcdf lib> -lnetcdff
! Usage: ./cru_interp key=value ...   (see the program at the end of this file)
!
module cru_common
  ! Constants, types, settings, and small functions needed for the CRU interpolation
  ! (extracted from datatypes.F90 and model_utils.F90 so this file is self-contained)
  implicit none

  integer, parameter :: hours_per_year = 365 * 24  ! 8760
  real,    parameter :: seconds_per_day = 24. * 3600.
  real,    parameter :: PI      = 3.1415926
  real,    parameter :: mol_air = 28.96e-3 ! molar mass of air, kg
  real,    parameter :: mol_h2o = 18.0e-3  ! molar mass of water, kg
  real,    parameter :: solarC  = 1361.0   ! Solar constant, W/m2

  ! CRU netCDF file dimensions and grid definition
  integer, parameter :: N_Vegs = 10 ! pft_2011 vegetation types
  integer, parameter :: Nlon = 720, Nlat = 360, Ntime = 1460
  integer, parameter :: Hours_NCstep = hours_per_year/Ntime ! 6
  real,    parameter :: Lon0 = -180.0, Lat0 = -90.0 ! The first grid's coordinates
  real,    parameter :: Wlon = 360.0/Nlon, Wlat = 180.0/Nlat

  ! Data interpolation (hourly)
  integer, parameter :: TargetSteps_per_hour = 1
  integer, parameter :: SHshift = 182*24*TargetSteps_per_hour ! Southern Hemisphere data shift (to July 2nd)
  logical, parameter :: ShiftSHdata = .True.

  integer, parameter :: CO2Yrs     = 326  ! 326 = the years from 1700 to 2025
  real,    parameter :: CO2_Hist(CO2Yrs) = & ! CO2 concentration 1700~2025, ppm
  ! 1,     2,     3,     4,     5,     6,     7,     8,     9,     10
  (/276.59,276.62,276.65,276.67,276.70,276.72,276.75,276.78,276.80,276.83,  &
  276.86,276.89,276.92,276.95,276.98,277.02,277.05,277.09,277.13,277.17,  &
  277.21,277.25,277.30,277.35,277.40,277.45,277.50,277.54,277.58,277.62,  &
  277.64,277.66,277.66,277.65,277.63,277.59,277.54,277.49,277.42,277.35,  &
  277.28,277.22,277.15,277.09,277.03,276.98,276.93,276.89,276.85,276.81,  &
  276.78,276.74,276.71,276.68,276.65,276.63,276.61,276.59,276.58,276.58,  &
  276.58,276.60,276.62,276.66,276.70,276.76,276.83,276.91,276.99,277.09,  &
  277.19,277.29,277.40,277.51,277.63,277.73,277.86,277.98,278.11,278.24,  &
  278.39,278.54,278.72,278.90,279.10,279.30,279.52,279.74,279.97,280.20,  &
  280.44,280.67,280.90,281.13,281.35,281.57,281.78,281.98,282.18,282.37,  &
  282.54,282.71,282.86,283.01,283.14,283.27,283.38,283.49,283.57,283.65,  &
  283.71,283.76,283.79,283.81,283.81,283.80,283.77,283.73,283.69,283.64,  &
  283.60,283.55,283.51,283.48,283.46,283.46,283.47,283.50,283.53,283.58,  &
  283.64,283.70,283.76,283.84,283.92,284.01,284.11,284.23,284.37,284.53,  &
  284.71,284.92,285.14,285.36,285.58,285.79,285.98,286.14,286.28,286.39,  &
  286.46,286.51,286.52,286.52,286.49,286.46,286.44,286.41,286.40,286.40,  &
  286.42,286.45,286.51,286.59,286.71,286.85,287.02,287.21,287.42,287.64,  &
  287.87,288.11,288.34,288.56,288.77,288.97,289.16,289.35,289.56,289.77,  &
  290.00,290.25,290.51,290.80,291.13,291.48,291.86,292.25,292.65,293.05,  &
  293.46,293.86,294.25,294.62,294.98,295.31,295.60,295.86,296.11,296.34,  &
  296.57,296.80,297.04,297.29,297.57,297.87,298.17,298.49,298.82,299.15,  &
  299.49,299.84,300.20,300.58,300.97,301.36,301.76,302.16,302.55,302.92,  &
  303.28,303.63,303.97,304.32,304.66,305.01,305.37,305.73,306.10,306.47,  &
  306.84,307.21,307.58,307.96,308.32,308.69,309.04,309.38,309.70,310.01,  &
  310.27,310.50,310.68,310.82,310.93,311.03,311.12,311.24,311.36,311.52,  &
  311.71,311.92,312.16,312.44,312.75,313.10,313.48,313.93,314.69,315.47,  &
  316.38,317.09,317.88,318.40,319.00,319.40,320.71,321.50,322.35,323.90,  &
  324.94,325.55,326.67,328.87,329.36,330.25,331.13,332.73,334.58,335.81,  &
  337.80,339.31,340.54,342.11,343.64,345.15,346.43,348.25,350.45,352.08,  &
  353.27,354.56,355.50,356.27,357.78,359.71,361.39,362.70,365.44,367.26,  &
  368.42,370.02,372.15,374.62,376.33,378.45,380.51,382.53,384.42,386.05,  &
  388.23,390.03,392.13,394.88,396.94,399.24,402.52,404.71,406.94,409.63,  &
  414.21,416.41,418.53,421.08,424.61,427.35/)

  type :: climate_data_type
    integer :: year       ! Year
    integer :: doy        ! day of the year
    real :: hod           ! hour of the day
    real :: PAR           ! umol m-2 s-1
    real :: radiation     ! W/m2
    real :: Tair          ! air temperature,  K
    real :: Tsoil         ! soil temperature, K
    real :: RH            ! relative humidity, 0~1
    real :: rain          ! kgH2O m-2 s-1
    real :: windU         ! wind velocity (m s-1)
    real :: P_air         ! pa
    real :: CO2           ! mol/mol
    real :: eCO2          ! mol/mol
    real :: soilwater     ! soil moisture, vol/vol
    real :: N_input       ! N deposition (kgN/m2/yr)
  end type climate_data_type

  type :: grid_initial_type
    integer :: iLon ! grid number along Longitude (from -180 to 180)
    integer :: iLat ! grid number along Latitude (from -90 to 90)
    real    :: N_input      ! kgN m-2 yr-1
    real, pointer :: VegCover(:)        ! for pft2011_0.5x0.5.nc only
    real, pointer :: climate(:,:)       ! Ntimes, N_vars
  end type grid_initial_type

  ! Settings (defaults; override with key=value command-line arguments)
  character (len = 256) :: ncfilepath   = './'
  character (len = 256) :: veg_path     = './'
  character (len = 256) :: Ndp_path     = './'
  character (len = 256) :: filepath_out = './'
  character (len = 20)  :: ncversion    = 'crujra.v2.4.5d.'
  character (len = 20)  :: veg_file     = 'pft2011_0.5x0.5.nc' ! Vegetation coverage
  character (len = 50)  :: GridListFile = 'GlobalVegGridList.csv' ! in filepath_out
  !character (len = 5)   :: ncfields(7) = [character(len=5):: 'tmp','pre','dswrf','spfh','pres','ugrd','vgrd'] ! CRU-JRA
  character (len = 5)   :: ncfields(7)  = [character(len=5):: 'tmp','pre','tswrf','spfh','pres','ugrd','vgrd'] ! Radiation is from tswrf_v15
  character (len = 6)   :: GridIDFMT    = '(I6.6)' ! For the file name string (GridID)
  character (len=30)    :: Vegstr       = 'TOTAL_VEG'
  character (len=9)     :: VegID(N_Vegs)= [character(len=9) :: &
                           'SHRUBS_BD','SHRUBS_BE','SHRUBS_ND','SHRUBS_NE', &
                           'TREES_BD ','TREES_BE ','TREES_ND ','TREES_NE ', &
                           'GRASS_MAN','GRASS_NAT'] ! pft2011_0.5x0.5.nc
  character (len=6)     :: NdpID(4) = [character(len=6) :: &
                           'wetnoy','wetnhx','drynoy','drynhx']
  integer :: LowerLon   = 1
  integer :: UpperLon   = 720
  integer :: LowerLat   = 61
  integer :: UpperLat   = 320
  integer :: yr_start   = 2010
  integer :: yr_end     = 2019
  integer :: grid_No1   = 1
  integer :: grid_No2   = 56395
  integer :: N_VegGrids = 1
  integer :: StepLatLon = 1 ! Skip grids. 1: all; 2: one per 2x2 grids
  integer :: GridID     = 999999 ! = iLon*1000 + iLat
  real    :: N_input    = 0.0    ! kgN m-2 yr-1, uniform N input to all grids (closed N: 0)
  logical :: WriteForcing = .True. ! Write interpolated forcing data
  logical :: Interpolate  = .True. ! .False.: skip the hourly interpolation and write the 6-hourly data
  integer :: Grids_UN1  = 99     ! Vegetated grids list file unit

  ! Interpolation results and grid data
  integer :: datalines, data_days, data_yrs, steps_per_hour = 1, steps_per_day = 24
  real    :: step_hour = 1.0, step_seconds = 3600.0, dt_fast_yr = 1.0/(365.0*24.0)
  type(grid_initial_type), pointer :: LandGrid(:) => null()
  type(climate_data_type), pointer :: forcingData(:) => null()
  integer,allocatable     :: GridLonLat(:)     ! LonLat
  real,allocatable,target :: GridVegT(:,:)     ! Grid vegetation composition
  real,allocatable,target :: ClimData(:,:,:)   ! N_yr*Ntime, N_vars, N_VegGrids
  real,allocatable        :: CRUtime(:)        ! Days since 1901-01-01 in CRU data

contains
  ! Weng, 05/24/2018
  subroutine calc_solarzen(td,latdegrees,cosz,solarelev,solarzen)
    !* Calculate solar zenith angle **in radians**
    !* From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
    implicit none
    real,intent(in) :: td            ! day(to minute fraction)
    real,intent(in) :: latdegrees    ! latitude in degrees
    real :: hour,latrad
    real :: delta    ! declination angle
    real :: pi, rad
    real,intent(out) :: cosz         ! cosz=cos(zen angle)=sin(elev angle)
    real,intent(out) :: solarelev    ! solar elevation angle (rad)
    real,intent(out) :: solarzen     ! solar zenith angle (rad)
    pi  = 3.1415926
    rad = pi / 180.0 ! Conversion from degrees to radians.
    hour = (td-floor(td))*24.0
    latrad = latdegrees*rad
    delta  = asin(-sin(rad*23.450)*cos(2.0*pi*(td+10.0)/365.0))
    cosz = sin(latrad)*sin(delta) + &
    cos(latrad)*cos(delta)*cos(rad* 15.0*(hour-12.0))
    cosz = max (cosz, 0.0)  ! Sun's angular is 0.01

    ! compute the solar elevation and zenth angles below
    solarelev = asin(cosz)/pi*180.0  !since asin(cos(zen))=pi/2-zen=elev
    solarzen = 90.0 - solarelev ! pi/2.d0 - solarelev
  end subroutine calc_solarzen

  !===========================
  FUNCTION esat(T) ! pressure, Pa
    IMPLICIT NONE
    REAL :: esat
    REAL, INTENT(IN) :: T ! degC
    esat=610.78*exp(17.27*T/(T+237.3))
  END FUNCTION esat
end module cru_common

module cru_interp_mod
  use netcdf
  use cru_common
  implicit none
  private
  public ReadNCfiles, CRU_Interpolation, Write_6hourly, Release_grid_memory

contains
!===================================================
subroutine ReadNCfiles (cru_path, veg_path, ndp_path,file_exists)
  ! Read in all netcdf files for global run, including;
  ! climate, vegetation, and nitrogen deposition data
  ! Updated 04/04/2026
   implicit none
   character(len=*),intent(in) :: cru_path,veg_path,ndp_path
   logical, intent(out) :: file_exists

   !-------- local vars -----------------
   character(len=256) :: fcru, fveg, fndp, fgz
   character(len=256) :: fout, command
   character(len=4)   :: yr_str
   character(len=6)   :: GridStr
   logical :: gz_exists, unzip_ok ! for zipped netcdf files
   integer :: N_yrs,totL,N_vars
   integer :: istat1,istat2
   integer :: i,j,k,m,iLon,iLat
   integer :: nlon_sub, nlat_sub
   integer :: start3(3), count3(3),ndpct3(3) ! for climate and Ndpst
   integer :: start2(2), count2(2)
   real, allocatable :: GridMask(:,:)  ! Nlon, Nlat
   real, allocatable :: dataarray(:,:,:), timearray(:)
   real, allocatable :: VegCover(:,:,:),  Vegetated(:,:)
   real, allocatable :: NdpData(:,:,:),   TotalNdp(:,:) ! N deposition


   ! Read in a vegetation map
   allocate(GridMask(LowerLon:UpperLon, LowerLat:UpperLat))

   ! Sub-domain sizes and hyperslab indices for netCDF reads
   nlon_sub = UpperLon - LowerLon + 1
   nlat_sub = UpperLat - LowerLat + 1
   start2 = [LowerLon, LowerLat]
   count2 = [nlon_sub, nlat_sub]
   start3 = [LowerLon, LowerLat, 1]
   count3 = [nlon_sub, nlat_sub, Ntime] ! CRUJRE climate data
   ndpct3 = [nlon_sub, nlat_sub, 12]    ! 12 months of N deposition (kgN m-2 s-1)

   ! Allocate (sub-domain) temporary arrays for netCDF reading
   allocate(dataarray(LowerLon:UpperLon, LowerLat:UpperLat, Ntime))
   allocate(timearray(Ntime))
   allocate(Vegetated(LowerLon:UpperLon, LowerLat:UpperLat))
   allocate(VegCover(LowerLon:UpperLon, LowerLat:UpperLat, N_Vegs))
   allocate(NdpData(LowerLon:UpperLon, LowerLat:UpperLat, 12)) ! kgN m-2 s-1
   allocate(TotalNdp(LowerLon:UpperLon, LowerLat:UpperLat))    ! kgN m-2 yr-1

   ! Read in vegetation data (vegtated grids and vegetation components)
   fveg  = trim(veg_path)//trim(veg_file)
   inquire(file=trim(fveg), exist=file_exists)
   if(.not. file_exists)then
      write(*,*)trim(fveg), ' file does not exist.'
      return
   endif
   call nc_read_2D(fveg, trim(Vegstr), Vegetated, start2, count2)
   do i=1, N_Vegs
     call nc_read_2D(fveg, trim(VegID(i)), VegCover(:,:,i), start2, count2)
   enddo
   write(*,*)'Successfully Read ',trim(fveg)

   TotalNdp = N_input ! from an assumed N input for all grids in all years
    ! -------------- Select land grids for model run ------------------!
    write(yr_str, '(I4)') yr_start
    !fcru = trim(cru_path)//trim(ncfields(1))//'/'//trim(ncversion)//trim(ncfields(1))//'.'//trim(yr_str)//'.365d.noc.nc'
    fcru = trim(cru_path)//'crujra3/'//trim(ncversion)//trim(ncfields(1))//'.'//trim(yr_str)//'.365d.noc.nc'
    fgz  = trim(fcru)//'.gz'

    ! If fcru exists, use fcru. Otherwise unzip fgz
    call unzip_gzip_file(trim(fgz), ok=unzip_ok)
    if (.not. unzip_ok) then
      file_exists = .false.
      forcingData => null()
      write (*, '("ReadNCfiles: unzip failed for ", a)') trim(fcru)//'.gz'
      return
    endif

    call nc_read_3D(fcru, trim(ncfields(1)), dataarray, start3, count3)

    ! Remove fcru if fgz exists
    inquire(file=trim(fgz), exist=gz_exists)
    if(gz_exists)then
       command = 'rm '//trim(fcru) ! Remove unziped file
       call execute_command_line(command)
    endif

    ! Tag the vegetated grids for model run
    GridMask(:,:) = 0
    m = 0
    do iLon = LowerLon, UpperLon, StepLatLon
      do iLat = LowerLat, UpperLat, StepLatLon
        if(dataarray(ilon,ilat,1) < 9999.0 .and. Vegetated(ilon,ilat)> 5.)then
          GridMask(ilon,ilat) = 1  ! Select the grids for model run
          m = m + 1
        endif
      enddo
    enddo
    N_VegGrids = m
    N_vars   = size(ncfields)
    N_yrs    = yr_end - yr_start + 1
    totL     = N_yrs * Ntime
    grid_No1 = 1
    grid_No2 = N_VegGrids ! ! Run all the grids in GridLonLat
    write(*,*)'LowerLon,UpperLon,LowerLat,UpperLat', LowerLon, UpperLon,LowerLat,UpperLat
    write(*,*)'Valid grids: grid_No1, grid_No2', grid_No1, grid_No2
    if(grid_No1 > grid_No2)then
      write(*,*)'No legible grids', N_VegGrids
      error stop
    endif

    ! Allocate data arrays
    !allocate(CRUData(totL, 4, LowerLon:UpperLon, LowerLat:UpperLat))
    allocate(CRUtime(totL))
    allocate(ClimData(totL, N_vars, N_VegGrids))
    allocate(LandGrid(N_VegGrids))
    allocate(GridLonLat(N_VegGrids))
    allocate(GridVegT(N_Vegs,N_VegGrids))
    write(*,*)'Grid arrays allocated.'

    ! Set GridLonLat array and Sort grid lon-lat and climate data
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0) then
          m = m + 1
          GridLonLat(m) = iLon * 1000 + iLat
          GridVegT(:,m) = VegCover(iLon,iLat,:)
          LandGrid(m)%iLon = iLon
          LandGrid(m)%iLat = iLat
          LandGrid(m)%N_input = TotalNdp(iLon,iLat)
          LandGrid(m)%VegCover => GridVegT(:,m) ! VegCover(iLon,iLat,:)
          LandGrid(m)%climate  => ClimData(:,:,m)
        endif
      enddo
    enddo
  
    ! Write GridLonLat and forcing file names to a file
    if(WriteForcing)then
      fout = trim(filepath_out)//trim(GridListFile) ! Grid ID, VegCover, and N_input
      open(NEWUNIT=Grids_UN1,file=trim(fout),ACTION='write', IOSTAT=istat1)
      !do m=1, N_VegGrids
      !  write(Grids_UN1, '(I8,11(",",E9.4))')GridLonLat(m),(GridVegT(i,m),i=1,N_vegs),LandGrid(m)%N_input
      !enddo
    endif


    ! ----------------- Read in all climate data ----------------------!
    do j= 1, N_vars ! 7 ('tmp','pre','dswrf'/'tswrf','spfh','pres','ugrd','vgrd')
      do i =1, N_yrs
        write(yr_str, '(I4)') yr_start + i - 1
        if(index(ncfields(j),'tswrf') > 0)then
          fcru = trim(cru_path)//'Radiation/Processed/'//trim(ncfields(j))//'_v15_'//trim(yr_str)//'.nc'
        else
          !fcru = trim(cru_path)//trim(ncfields(j))//'/'//trim(ncversion)//trim(ncfields(j))//'.'//trim(yr_str)//'.365d.noc.nc'
          fcru = trim(cru_path)//'crujra3/'//trim(ncversion)//trim(ncfields(j))//'.'//trim(yr_str)//'.365d.noc.nc'
        endif
        fgz  = trim(fcru)//'.gz'

        ! If the climate file is gzipped, unzip it and returen unzip_ok. 
        ! Otherwise, return unzip_ok (file exists)
        call unzip_gzip_file(trim(fgz), ok=unzip_ok)
        if (.not. unzip_ok) then
          file_exists = .false.
          forcingData => null()
          write (*, '("read_interpolatedCRU: unzip failed for ", a)') trim(fgz)
          return
        endif

        write(*,*)'Reading: ', trim(fcru)
        call nc_read_3D(fcru, trim(ncfields(j)), dataarray, start3, count3)
        m = 0
        do iLon = LowerLon, UpperLon
          do iLat = LowerLat, UpperLat
            !CRUData((i-1)*Ntime+1:i*Ntime,j,iLon,iLat) = dataarray(iLon,iLat,:)
            if(GridMask(iLon,iLat) > 0)then
              m = m + 1
              if(index(ncfields(j),'tswrf') > 0)then ! mean shortwave flux over the 6 hours after the timestamp
                ClimData((i-1)*Ntime+1:i*Ntime,j,m) = dataarray(iLon,iLat,:) * 21600.0 ! convert to J/m2 for 6 hours
              else
                ClimData((i-1)*Ntime+1:i*Ntime,j,m) = dataarray(iLon,iLat,:)
              endif
            endif
          enddo
        enddo

        ! Read in the time array of the first variable
        if(j == 1)then
          call nc_read_1D(fcru, 'time', timearray)
          CRUtime((i-1)*Ntime+1: i*Ntime) = timearray
        endif

        ! Remove unzipped file if fgz exists
        inquire(file=trim(fgz), exist=gz_exists)
        if(gz_exists)then
           command = 'rm '//trim(fcru) ! Remove unziped file
           call execute_command_line(command)
        endif
      enddo ! N_yrs
    enddo   ! All variables

    ! ---- Release temporary allocatable arrays ----
    if(allocated(GridMask))  deallocate(GridMask)
    if(allocated(dataarray)) deallocate(dataarray)
    if(allocated(timearray)) deallocate(timearray)
    if(allocated(Vegetated)) deallocate(Vegetated)
    if(allocated(VegCover))  deallocate(VegCover)
    if(allocated(NdpData))   deallocate(NdpData)
    if(allocated(TotalNdp))  deallocate(TotalNdp)

end subroutine ReadNCfiles

!=============================================================================
subroutine CRU_Interpolation(LandGrid,forcingData)
  implicit none
  type(grid_initial_type), pointer, intent(in) :: LandGrid
  type(climate_data_type), pointer :: forcingData(:) ! output

  !---------- local variables ------------------
  type(climate_data_type), pointer :: climateData(:) ! will be pointed by forcingData
  character(len=256):: command, fout   ! For testing output
  character(len=6)  :: GridStr  ! Used in output file name
  real, pointer     :: GridData(:,:)
  real, allocatable :: fdSW(:)
  real, allocatable :: timecols(:,:)
  real, allocatable :: hourly_data(:,:)
  real    :: Lati, Longi
  real    :: steps_in_6H ! Temporary variable, steps interpolated
  real    :: td,cosz,solarelev,solarzen,r_light
  real    :: WindS1, WindS2, SWdaily, SWmax
  real    :: tmLocal ! Local time, may be negative values as the previous day 
  real    :: tmp1(23,10), tmp2(SHshift,10)  ! Shift hourly data
  integer :: forcing_unit ! for interpolated grid forcing file writting
  integer :: year0, year1 ! Start and end year
  integer :: yr,doy,iday,ihour,iyr
  integer :: ndays,nyear,totalL,Nsteps
  integer :: iostat,m,i,j,k
  integer :: totyr, totDays, Nlines
  real    :: SWmax_doy(365)               ! Daily maximum SW (J m-2) for each day of the year
  real, allocatable :: cosz_doy(:,:)      ! cos(zenith) for each day of the year and time step of the day

  ! Assigne LandGrid climate data to local variables
  GridData => LandGrid%climate

  ! Latitude and Longitude of this grid, (started from -179.75 and -89.75)
  Longi = Lon0 + (LandGrid%iLon - 0.5) * Wlon
  Lati  = Lat0 + (LandGrid%iLat - 0.5) * Wlat

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
    tmLocal       = MODULO(CRUtime(m)*24.0,24.0) + int(Longi/15.0) ! Local time
    timecols(m,3) = MODULO(tmLocal, 24.0)          ! Hour, in the range of 0~23
    timecols(m,2) = MODULO(int(CRUtime(m)),365) + 1.0 ! Day of the year, 1~365
    timecols(m,1) = int(CRUtime(m)/365.0) + 1901   ! Year, CRUtime is days since 1901-01-01
  enddo

  ! Daily maximum SW and cos(zenith) depend only on the latitude and the day of the year,
  ! so calculate them once for the 365 days (not for every day of every year)
  allocate(cosz_doy(365,steps_per_day))
  do k = 1, 365
    SWmax_doy(k) = 0.0
    do j=1,96 ! a day with 15 minutes interval
      td = real(k) + (j-1)/96.0
      call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
      SWmax_doy(k) = SWmax_doy(k) + cosz * solarC * seconds_per_day/96.0
    enddo
    do j=1,steps_per_day
      td = real(k) + (24.0 * (j - 1)/steps_per_day)/24.0
      call calc_solarzen(td,Lati,cosz_doy(k,j),solarelev,solarzen)
    enddo
  enddo

  ! Calculate fdSW
  Nsteps = 24/Hours_NCstep ! netCDF Data daily steps, 4 for CRU
  do i = 1, totDays
    ! Actual daily radiation
    SWdaily = sum(GridData((i-1)*Nsteps+1 : i*Nsteps, 3))
    ! Daily max SW (SWmax) of this day of the year
    SWmax = SWmax_doy(nint(timecols((i-1)*Nsteps+1,2))) ! Corrected by Paul, 09/08/2026
    ! Fraction of solar radiation
    fdSW((i-1)*Nsteps+1 : i*Nsteps) = Max(0.0,Min(1.0,SWdaily/(SWmax+0.0001))) ! *0.8 ! too high!
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
    j = nint(hourly_data(i,3) * steps_per_day / 24.0) + 1 ! time step of the day
    hourly_data(i,6) = solarC * cosz_doy(nint(hourly_data(i,2)),j) * hourly_data(i,10)      ! W/m2
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
     climateData(i)%CO2   = CO2_Hist(Min(CO2Yrs,Max(1,climateData(i)%year-1700+1))) ! CO2_c   ! ppm
     climateData(i)%eCO2  = climateData(i)%CO2 + 200.       ! ppm
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
     climateData(i)%N_input = LandGrid%N_input
  enddo
  forcingData => climateData
  datalines = totalL
  data_days = totyr * 365
  data_yrs   = totyr
  dt_fast_yr= 1.0/(365.0 * 24.0 * steps_per_hour)
  step_hour = 1.0/steps_per_hour
  step_seconds  = 3600.0 * step_hour

  ! Write out a sample file
  if(WriteForcing)then
    write(GridStr, GridIDFMT) GridID

    ! Write interpolated data
    fout = trim(filepath_out)//trim(ncversion)//trim(GridStr)//'_forcing.csv' ! Data file name
    open(NEWUNIT=forcing_unit,file=trim(fout))
    write(forcing_unit,*)'Swdown,Tair,RH,RAIN,WIND,PRESSURE,Ndpst'
    do i=1,totalL
        write(forcing_unit,'(7(E15.8,","))') &
          climateData(i)%radiation,climateData(i)%Tair,climateData(i)%RH,  &
          climateData(i)%rain,climateData(i)%windU,climateData(i)%P_air,   &
          climateData(i)%N_input
    enddo
    close(forcing_unit)
    command = 'gzip -f ' // trim(fout)
    call execute_command_line(command, exitstat=iostat)

    ! ------- Write Grid list files -------
    ! Write GridID, VegCover, and N_input for each grid
    write(Grids_UN1,'(I8,11(",",E9.4))') &
          GridID,(LandGrid%VegCover(i),i=1,N_vegs),LandGrid%N_input
  endif

  !Release memory
  deallocate(hourly_data, fdSW, timecols, cosz_doy)
end subroutine CRU_Interpolation


!=============================================================================
subroutine Write_6hourly(LandGrid)
  ! Write the 6-hourly CRU data of one grid without interpolation or any time shift
  ! (neither the local-time nor the Southern Hemisphere shift; BiomeE applies them
  ! when it interpolates). The columns are the values stored in ClimData: tmp (K), pre (mm per 6 hours), tswrf (J m-2
  ! per 6 hours), spfh, pres (Pa), ugrd, vgrd (m s-1), and N input (kgN m-2 yr-1).
  implicit none
  type(grid_initial_type), pointer, intent(in) :: LandGrid
  real, pointer     :: GridData(:,:)
  character(len=256):: command, fout
  character(len=6)  :: GridStr
  integer :: forcing_unit, iostat, i, Nlines

  GridData => LandGrid%climate
  Nlines = SIZE(CRUtime)

  write(GridStr, GridIDFMT) GridID
  fout = trim(filepath_out)//trim(ncversion)//trim(GridStr)//'_forcing6h.csv'
  open(NEWUNIT=forcing_unit,file=trim(fout))
  write(forcing_unit,*)'Tair,RAIN,SWdown6h,SPFH,PRESSURE,WINDU,WINDV,Ndpst'
  do i=1,Nlines
    write(forcing_unit,'(8(E15.8,","))') GridData(i,1:7), LandGrid%N_input
  enddo
  close(forcing_unit)
  command = 'gzip -f ' // trim(fout)
  call execute_command_line(command, exitstat=iostat)

  write(Grids_UN1,'(I8,11(",",E9.4))') &
        GridID,(LandGrid%VegCover(i),i=1,N_vegs),LandGrid%N_input
end subroutine Write_6hourly

!==============================================
subroutine Release_grid_memory()
  if(WriteForcing) close(Grids_UN1)
  if(allocated(GridLonLat)) deallocate(GridLonLat)
  if(allocated(GridVegT))   deallocate(GridVegT)
  !deallocate(CRUData)
  if(allocated(CRUtime))  deallocate(CRUtime)
  if(allocated(ClimData)) deallocate(ClimData)
  deallocate(LandGrid)
end subroutine Release_grid_memory

!=============================================================================

!==============================================================
subroutine unzip_gzip_file(filename_gz, ok, exitstat, out_file)
    character(len=*), intent(in)            :: filename_gz
    logical,          intent(out), optional :: ok
    integer,          intent(out), optional :: exitstat
    character(len=*), intent(out), optional :: out_file

    character(len=600) :: command
    character(len=600) :: filename_out
    logical :: gz_exists, out_exists
    integer :: est, Lfn

    if (present(ok)) ok = .false.
    if (present(exitstat)) exitstat = -999

    ! --------- Derive output filename by stripping trailing ".gz" if present
    filename_out = trim(filename_gz)
    Lfn = len_trim(filename_out)
    if (Lfn <= 3) error stop 'Specify file types, .gz, .csv, or .txt'
    if (filename_out(Lfn-2: Lfn) == '.gz') filename_out = filename_out(1: Lfn-3)
    if (present(out_file)) out_file = trim(filename_out)

    ! If output already exists, treat as success (idempotent)
    inquire(file=trim(filename_out), exist=out_exists)
    if (out_exists) then
      est = 0
      if (present(exitstat)) exitstat = est
      if (present(ok)) ok = .true.
      return
    endif

    ! --------- Check input exists
    ! Check the netcdf file so that to skip unzipping if the file exists
    inquire(file=trim(filename_gz), exist=gz_exists)
    if (.not. gz_exists) then
      write(*,'("unzip_gzip_file: missing ",a)') trim(filename_gz)
      est = 2
      if (present(exitstat)) exitstat = est
      return
    endif

    ! Unzip: keep .gz (-k), force overwrite (-f)
    command = 'gunzip -kf ' // trim(filename_gz)
    call execute_command_line(command, exitstat=est)
    if (present(exitstat)) exitstat = est

    ! --------- Verify output exists
    inquire(file=trim(filename_out), exist=out_exists)
    if (out_exists .and. est == 0) then
      if (present(ok)) ok = .true.
    else
      write(*,'("unzip_gzip_file: expected output missing: ",a)') trim(filename_out)
      if (present(ok)) ok = .false.
      return
    endif
end subroutine unzip_gzip_file

!==============================================
  subroutine nc_read_3D(file_name, var_name, da, start, count)
    ! Read a 3-D netCDF variable. Optionally read a hyperslab using start/count.
    character(len=*), intent(in) :: file_name, var_name
    real, intent(out) :: da(:,:,:)
    integer, intent(in), optional :: start(3), count(3)

    integer :: ncid, varid, ndims, xtype
    integer :: dimids(3)
    integer :: dlen(3)
    integer :: s(3), c(3)

    call check( nf90_open(trim(file_name), NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(var_name), varid) )
    call check( nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, dimids=dimids) )
    if (ndims /= 3) then
      write(*,*) 'ERROR: variable is not 3-D: ', trim(var_name), ' in ', trim(file_name)
      error stop
    endif

    call check( nf90_inquire_dimension(ncid, dimids(1), len=dlen(1)) )
    call check( nf90_inquire_dimension(ncid, dimids(2), len=dlen(2)) )
    call check( nf90_inquire_dimension(ncid, dimids(3), len=dlen(3)) )

    if (present(start) .and. present(count)) then
      s = start
      c = count
      if (size(da,1) /= c(1) .or. size(da,2) /= c(2) .or. size(da,3) /= c(3)) then
        write(*,*) 'ERROR: hyperslab shape mismatch reading ', trim(var_name)
        write(*,*) ' start/count=', s, c, ' array=', size(da,1), size(da,2), size(da,3)
        error stop
      endif
      if (any(s < 1)) then
        write(*,*) 'ERROR: invalid start (<1) reading ', trim(var_name), ' start=', s
        error stop
      endif
      if (any(s + c - 1 > dlen)) then
        write(*,*) 'ERROR: hyperslab out of bounds reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' start/count=', s, c
        error stop
      endif
      call check( nf90_get_var(ncid, varid, da, start=s, count=c) )
    else
      if (size(da,1) /= dlen(1) .or. size(da,2) /= dlen(2) .or. size(da,3) /= dlen(3)) then
        write(*,*) 'ERROR: full-field shape mismatch reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' array=', size(da,1), size(da,2), size(da,3)
        error stop
      endif
      call check( nf90_get_var(ncid, varid, da) )
    endif

    call check( nf90_close(ncid) )
  end subroutine nc_read_3D

  !==============================================
  subroutine nc_read_2D(file_name, var_name, da, start, count)
    ! Read a 2-D netCDF variable. Optionally read a hyperslab using start/count.
    character(len=*), intent(in) :: file_name, var_name
    real, intent(out) :: da(:,:)
    integer, intent(in), optional :: start(2), count(2)

    integer :: ncid, varid, ndims, xtype
    integer :: dimids(2)
    integer :: dlen(2)
    integer :: s(2), c(2)

    call check( nf90_open(trim(file_name), NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(var_name), varid) )
    call check( nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims) )
    if (ndims /= 2) then
      write(*,*) 'ERROR: variable is not 2-D: ', trim(var_name), ' in ', trim(file_name)
      error stop
    endif
    call check( nf90_inquire_variable(ncid, varid, dimids=dimids) )
    call check( nf90_inquire_dimension(ncid, dimids(1), len=dlen(1)) )
    call check( nf90_inquire_dimension(ncid, dimids(2), len=dlen(2)) )

    if (present(start) .and. present(count)) then
      s = start
      c = count
      if (size(da,1) /= c(1) .or. size(da,2) /= c(2)) then
        write(*,*) 'ERROR: hyperslab shape mismatch reading ', trim(var_name)
        write(*,*) ' start/count=', s, c, ' array=', size(da,1), size(da,2)
        error stop
      endif
      if (any(s < 1)) then
        write(*,*) 'ERROR: invalid start (<1) reading ', trim(var_name), ' start=', s
        error stop
      endif
      if (any(s + c - 1 > dlen)) then
        write(*,*) 'ERROR: hyperslab out of bounds reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' start/count=', s, c
        error stop
      endif
      call check( nf90_get_var(ncid, varid, da, start=s, count=c) )
    else
      if (size(da,1) /= dlen(1) .or. size(da,2) /= dlen(2)) then
        write(*,*) 'ERROR: full-field shape mismatch reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' array=', size(da,1), size(da,2)
        error stop
      endif
      call check( nf90_get_var(ncid, varid, da) )
    endif

    call check( nf90_close(ncid) )
  end subroutine nc_read_2D

  !==============================================
  subroutine nc_read_1D(file_name, var_name, da)
    ! Read a 1-D netCDF variable.
    character(len=*), intent(in) :: file_name, var_name
    real, intent(out) :: da(:)

    integer :: ncid, varid, ndims, xtype
    integer :: dimids(1)
    integer :: dlen(1)

    call check( nf90_open(trim(file_name), NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(var_name), varid) )
    call check( nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, dimids=dimids) )
    if (ndims /= 1) then
      write(*,*) 'ERROR: variable is not 1-D: ', trim(var_name), ' in ', trim(file_name)
      error stop
    endif
    call check( nf90_inquire_dimension(ncid, dimids(1), len=dlen(1)) )

    if (size(da,1) /= dlen(1)) then
      write(*,*) 'ERROR: shape mismatch reading ', trim(var_name), ' dims=', dlen(1), ' array=', size(da,1)
      error stop
    endif

    call check( nf90_get_var(ncid, varid, da) )
    call check( nf90_close(ncid) )
  end subroutine nc_read_1D

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop "Stopped"
    end if
  end subroutine check

end module cru_interp_mod

!==============================================================================
program CRU_Interp
  ! Settings are taken from the defaults in datatypes.F90 and can be overridden on
  ! the command line with key=value arguments (no namelist file is read), e.g.
  !   ./cru_interp ncfilepath=/data/CRU/ veg_path=/data/Veg/ filepath_out=./out/ \
  !                LowerLon=1 UpperLon=180 yr_start=1990 yr_end=2019
  ! Keys: ncfilepath ncversion veg_path veg_file ndp_path filepath_out GridListFile
  !       LowerLon UpperLon LowerLat UpperLat yr_start yr_end StepLatLon
  !       N_input Interpolate (T: hourly interpolation, F: write the 6-hourly data only)
  use cru_common, only: ncfilepath, ncversion, veg_path, veg_file, Ndp_path,  &
                       filepath_out, GridListFile, LowerLon, UpperLon,       &
                       LowerLat, UpperLat, yr_start, yr_end, StepLatLon,     &
                       N_input, WriteForcing, Interpolate, Nlon, Nlat, GridID, GridLonLat,&
                       LandGrid, forcingData, grid_No1, grid_No2
  use cru_interp_mod
  implicit none

  real    :: start_time, end_time, last_time
  integer :: m, nWritten, i, ieq
  character(len=512) :: arg
  logical :: file_exists

  call cpu_time(start_time)

  ! Closed N cycle default of the original global run (N input = 0)
  N_input = 0.0
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    ieq = index(arg, '=')
    if (ieq < 2) then
      write(*,*) 'Ignoring argument (expected key=value): ', trim(arg)
      cycle
    endif
    call set_option(trim(adjustl(arg(1:ieq-1))), trim(adjustl(arg(ieq+1:))))
  enddo
  ! Keep the grid ranges within the global grid
  LowerLon = max(LowerLon, 1);  UpperLon = min(UpperLon, Nlon)
  LowerLat = max(LowerLat, 1);  UpperLat = min(UpperLat, Nlat)

  WriteForcing = .True.   ! This program exists to write interpolated forcing files

  call ReadNCfiles(ncfilepath, veg_path, ndp_path, file_exists)
  if (.not. file_exists) then
    write(*,*) 'NetCDF file missing. Stopped'
    error stop
  endif
  call cpu_time(end_time)
  write(*,'(A,f8.2)') 'Data reading CPU time (minutes): ', (end_time-start_time)/60.

  nWritten = 0
  do m = grid_No1, grid_No2
    GridID = GridLonLat(m)
    call cpu_time(last_time)
    print '(A, I6, A, I6)', 'Interpolating grid: ', GridID, '. Grid No. ', m
    if(Interpolate)then
      call CRU_Interpolation(LandGrid(m), forcingData)
      deallocate(forcingData)
    else
      call Write_6hourly(LandGrid(m))
    endif
    nWritten = nWritten + 1
    call cpu_time(end_time)
    write(*,'(a6,f6.1,a12,I6,a17,f7.1,a9)') 'Used: ', end_time-last_time, ' seconds at ', &
         GridID, '. Total CPU time: ', (end_time-start_time)/60.0, ' minutes.'
  enddo

  call Release_grid_memory()
  write(*,'(A,I8,A)') 'Interpolated forcing written for ', nWritten, ' grids.'

contains
  subroutine set_option(key, val)
    character(len=*), intent(in) :: key, val
    select case (key)
    case ('ncfilepath');   ncfilepath   = val
    case ('ncversion');    ncversion    = val
    case ('veg_path');     veg_path     = val
    case ('veg_file');     veg_file     = val
    case ('ndp_path');     Ndp_path     = val
    case ('filepath_out'); filepath_out = val
    case ('GridListFile'); GridListFile = val
    case ('LowerLon');     read(val,*) LowerLon
    case ('UpperLon');     read(val,*) UpperLon
    case ('LowerLat');     read(val,*) LowerLat
    case ('UpperLat');     read(val,*) UpperLat
    case ('yr_start');     read(val,*) yr_start
    case ('yr_end');       read(val,*) yr_end
    case ('StepLatLon');   read(val,*) StepLatLon
    case ('N_input');      read(val,*) N_input
    case ('Interpolate');  read(val,*) Interpolate
    case default
      write(*,*) 'Unknown option: ', key
      error stop
    end select
  end subroutine set_option
end program CRU_Interp
