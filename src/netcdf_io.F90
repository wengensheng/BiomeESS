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
  use model_utils
  implicit none

  private

  public ReadNCfiles, CRU_Interpolation, CRU_end
  public read_GridLonLat, read_interpolatedCRU
  public unzip_gzip_file

contains

!===================================================
subroutine ReadNCfiles (cru_path, veg_path, ndp_path)
  ! Read in all netcdf files for global run, including;
  ! climate, vegetation, and nitrogen deposition data
  ! Updated 04/04/2026
   use, intrinsic :: ieee_arithmetic
   implicit none
   character(len=*),intent(in) :: cru_path,veg_path,ndp_path

   !-------- local vars -----------------
   character(len=256) :: fcru, fveg, fndp
   character(len=256) :: fout, command
   character(len=4)   :: yr_str
   character(len=6)   :: GridStr
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

   ! For Hurtt's farmland data
   integer :: crpct3(3)
   real, allocatable :: CropData(:,:,:),  TotCrops(:,:,:)

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
   call nc_read_2D(fveg, trim(Vegstr), Vegetated, start2, count2)
   do i=1, N_Vegs
     call nc_read_2D(fveg, trim(VegID(i)), VegCover(:,:,i), start2, count2)
   enddo
   write(*,*)'Successfully Read ',trim(fveg)

#ifdef WIEMIP_setting
   ! For Hurtt's farmland data
   crpct3 = [nlon_sub, nlat_sub, FM_Yrs]
   allocate(CropData(LowerLon:UpperLon, LowerLat:UpperLat, FM_Yrs))
   allocate(TotCrops(LowerLon:UpperLon, LowerLat:UpperLat, FM_Yrs))

   ! Read in farmland fractions from Hurtt's data
   fveg  = trim(veg_path)//'CropLand_states4_05deg.nc'
   TotCrops = 0.0
   do i = 1, N_crop
     call nc_read_3D(fveg, trim(CropID(i)), CropData, start3, crpct3)
     !write(*,*) 'file = ', trim(fveg)
     !write(*,*) 'var  = ', trim(CropID(i))
     !write(*,*) 'valid count = ', count(.not. ieee_is_nan(CropData))
     where (ieee_is_nan(CropData)) CropData = 0.0
     TotCrops = TotCrops + CropData
   enddo
   write(*,*) 'Grids with farmland>0.01 = ', count(TotCrops(:,:,1170)>0.01)
   write(*,*)'Successfully Read ',trim(fveg)

   ! Read in N deposition and calculate yearly N deposition rate (kgN/m2/yr)
   TotalNdp = 0.0
   do i=1, 4 ! Files for 4 types of N deposition
     fndp = trim(Ndp_path)//trim(NdpID(i))//'_185001-185012.nc'
     call nc_read_3D(fndp, trim(NdpID(i)), NdpData, start3, ndpct3)
     TotalNdp = TotalNdp + SUM(NdpData, DIM=3)
     write(*,*)'Successfully Read ',trim(fndp)
   enddo
   TotalNdp = TotalNdp * (3600.*24.*365.) /12.
#else
   TotalNdp = N_input ! from an assumed N input for all grids in all years
#endif
    ! -------------- Select land grids for model run ------------------!
    write(yr_str, '(I4)') yr_start
    fcru = trim(cru_path)//trim(ncfields(1))//'/'//trim(ncversion)//trim(ncfields(1))//'.'//trim(yr_str)//'.365d.noc.nc'

#ifdef ZippedNCfiles
    call unzip_gzip_file(trim(fcru)//'.gz')
#endif

    call nc_read_3D(fcru, trim(ncfields(1)), dataarray, start3, count3)

#ifdef ZippedNCfiles
    command = 'rm '//trim(fcru) ! Remove unziped file
    call execute_command_line(command)
#endif

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

    ! Allocate data arrays
    !allocate(CRUData(totL, 4, LowerLon:UpperLon, LowerLat:UpperLat))
    allocate(CRUtime(totL))
    allocate(ClimData(totL, N_vars, N_VegGrids))
    allocate(LandGrid(N_VegGrids))
    allocate(GridLonLat(N_VegGrids))
    allocate(GridVegCov(N_Vegs,N_VegGrids))
    write(*,*)'Grid arrays allocated.'

    ! Set GridLonLat array and Sort grid lon-lat and climate data
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0) then
          m = m + 1
          GridLonLat(m) = iLon * 1000 + iLat
          GridVegCov(:,m) = VegCover(iLon,iLat,:)
          LandGrid(m)%iLon = iLon
          LandGrid(m)%iLat = iLat
          LandGrid(m)%N_input = TotalNdp(iLon,iLat)
          LandGrid(m)%VegCover => GridVegCov(:,m) ! VegCover(iLon,iLat,:)
          LandGrid(m)%climate  => ClimData(:,:,m)
        endif
      enddo
    enddo

#ifdef WIEMIP_setting
    allocate(GridFarm(FM_Yrs,N_VegGrids))
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0) then
          m = m + 1
          GridFarm(:,m)    = TotCrops(iLon,iLat,:)
          LandGrid(m)%Farm => GridFarm(:,m)
        endif
      enddo
    enddo
#endif

    ! Write GridLonLat and forcing file names to a file
    if(WriteForcing)then
      fout = trim(filepath_out)//trim(GridListFile) ! Grid ID, VegCover, and N_input
      open(NEWUNIT=Grids_UN1,file=trim(fout),ACTION='write', IOSTAT=istat1)
      !do m=1, N_VegGrids
      !  write(Grids_UN1, '(I8,11(",",E9.4))')GridLonLat(m),(GridVegCov(i,m),i=1,N_vegs),LandGrid(m)%N_input
      !enddo
#ifdef WIEMIP_setting
      fout = trim(filepath_out)//'GridFarmRatio.csv' ! Data file name
      open(NEWUNIT=Grids_UN2,file=trim(fout),ACTION='write', IOSTAT=istat2)
      !do m=1, N_VegGrids
      !  write(Grids_UN2, '(I8, 1176(",",f12.4))')GridLonLat(m),(GridFarm(i,m),i=1,FM_Yrs)
      !enddo
#endif
    endif

    ! ----------------- Read in all climate data ----------------------!
    do j= 1, N_vars ! 7 ('tmp','pre','dswrf','spfh','pres','ugrd','vgrd')
      do i =1, N_yrs
        write(yr_str, '(I4)') yr_start + i - 1
        fcru = trim(cru_path)//trim(ncfields(j))//'/'//trim(ncversion)//trim(ncfields(j))//'.'//trim(yr_str)//'.365d.noc.nc'

#ifdef ZippedNCfiles
        call unzip_gzip_file(trim(fcru)//'.gz') ! Unzip the nc data file
#endif

        write(*,*)'Reading: ', trim(fcru)
        call nc_read_3D(fcru, trim(ncfields(j)), dataarray, start3, count3)
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
          call nc_read_1D(fcru, 'time', timearray)
          CRUtime((i-1)*Ntime+1: i*Ntime) = timearray
        endif

#ifdef ZippedNCfiles
        command = 'rm '//trim(fcru) ! Remove unziped nc data file
        call execute_command_line(command)
#endif
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
    if(allocated(CropData))  deallocate(CropData)
    if(allocated(TotCrops))  deallocate(TotCrops)

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
     climateData(i)%CO2   = CO2_Hist(Min(CO2Yrs,Max(1,climateData(i)%year-1700+1))) ! CO2_c   ! ppm
     climateData(i)%eCO2  = climateData(i)%CO2 + 200.       ! ppm
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
     climateData(i)%N_input = LandGrid%N_input
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
#ifdef WIEMIP_setting
    write(Grids_UN2,'(I8,1176(",",E9.4))')GridID,(LandGrid%Farm(i),i=1,FM_Yrs)
#endif
  endif

  !Release memory
  deallocate(hourly_data, fdSW, timecols)
end subroutine CRU_Interpolation

!==============================================
subroutine CRU_end()
  close(Grids_UN1)
  close(Grids_UN2)
  if(allocated(GridLonLat)) deallocate(GridLonLat)
  if(allocated(GridVegCov)) deallocate(GridVegCov)
  if(allocated(GridFarm))   deallocate(GridFarm)
#ifndef Use_InterpolatedData
  !deallocate(CRUData)
  if(allocated(CRUtime))  deallocate(CRUtime)
  if(allocated(ClimData)) deallocate(ClimData)
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
  integer, parameter :: maxGrids = Nlon*Nlat/3
  real    :: GridVF(N_Vegs, maxGrids),tmpVF(N_Vegs, maxGrids) ! Veg cover for all grids
#ifdef WIEMIP_setting
  real    :: GridFM(FM_Yrs, maxGrids),tmpFM(FM_Yrs, maxGrids)
#endif
  integer :: GridNo(maxGrids),tmpNo(maxGrids) ! maximum grids, 720*360
  integer :: GridCode, istat1, istat2
  integer :: i,j,k,m,n
  integer :: TotalFiles ! Each file represents a grid
  character(len=300) :: listfile,FarmFile

  listfile=trim(int_fpath)//trim(fname)
  INQUIRE (file=trim(listfile), EXIST=file_exists)
  if (.not. file_exists) then
    write (*, '("read_GridLonLat: ", a, " does not exist")') trim(listfile)
    return
  endif
  open(11,file=listfile,status='old',ACTION='read',IOSTAT=istat1)

#ifdef WIEMIP_setting
  FarmFile=trim(int_fpath)//'GridFarmRatio.csv'
  INQUIRE (file=trim(FarmFile), EXIST=file_exists)
  if (.not. file_exists) then
    write (*, '("read_GridLonLat: ", a, " does not exist")') trim(Farmfile)
    return
  endif
  open(12,file=Farmfile,status='old',ACTION='read',IOSTAT=istat2)
#endif

  m = 0
  do
    read(11,*,IOSTAT=istat1) GridNo(m+1), (tmpVF(i,m+1),i=1,10)
    if(istat1<0 ) exit
#ifdef WIEMIP_setting
    read(12,*,IOSTAT=istat2) GridCode, (tmpFM(i,m+1),i=1,FM_Yrs)
    if(istat2 < 0 .or. GridNo(m+1) /= GridCode) exit
#endif
    ! Next grid
    m = m + 1
  enddo
  write(*,*)'Total grids in the list fie:',m

  ! Update GridNo and m with StepLatLon
  if(StepLatLon > 1 .or. UpperLon < 720)then ! Regional or partial run
    TotalFiles = m
    tmpNo = GridNo
    m = 0
    n = 1
    do i = LowerLon, UpperLon, StepLatLon
      do j = LowerLat, UpperLat, StepLatLon
        k = i*1000 + j
        if (k > tmpNo(n)) then
          do while(k > tmpNo(n) .and. n < TotalFiles)
            n = n + 1
          enddo
        elseif (k < tmpNo(n)) then
          cycle
        endif
        if (k == tmpNo(n))then
            m = m + 1
            GridNo(m) = tmpNo(n)
            GridVF(:,m) = tmpVF(:,n)
#ifdef WIEMIP_setting
            GridFM(:,m) = tmpFM(:,n)
#endif
            n = n + 1
        endif
        if(n > TotalFiles)exit ! Exit LowerLat-UpperLat loop
      enddo
      if(n > TotalFiles)exit   ! Exit LowerLon-UpperLon loop
    enddo
  endif

  ! Update N_VegGrids and GridLonLat
  N_VegGrids = m
  grid_No1 = min(grid_No1,N_VegGrids)
  grid_No2 = min(grid_No2,N_VegGrids)

  ! Global variables
  allocate(GridLonLat(N_VegGrids))
  allocate(GridVegCov(N_Vegs,N_VegGrids))
  GridLonLat(:)   = GridNo(1:N_VegGrids)
  GridVegCov(:,:) = GridVF(:,1:N_VegGrids)
  close(11)

#ifdef WIEMIP_setting
  allocate(GridFarm(FM_Yrs,N_VegGrids))
  GridFarm(:,:)   = GridFM(:,1:N_VegGrids)
  close(12)
#endif

  write(*,*)"StepLatLon:", StepLatLon
  write(*,*)'Selected grids:',N_VegGrids
  write(*,*)"LowerLon, UpperLon, LowerLat, UpperLat:", LowerLon, UpperLon, LowerLat, UpperLat
  write(*,*)"In read_GridLonLat, N_VegGrids, grid_No1, grid_No2:", N_VegGrids, grid_No1, grid_No2
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
#ifdef WIEMIP_setting
  integer, parameter :: niterms = 7 !6 ! 7 columns (including N_input)
#else
  integer, parameter :: niterms = 6 ! No column of N deposition
#endif

  real,allocatable :: timecols(:,:), input_data(:,:)
  real    :: temp(niterms)
  integer :: istat2
  integer :: m,n,i,iyr,iday,ihour
  integer :: u, exitstat
  logical :: csv_exists, unzip_ok

  ! Find out the data file
  write(GridStr,GridIDFMT) GridID
  fname    = trim(fprefix)//trim(GridStr)//'_forcing.csv'
  climfile = trim(fpath)//trim(fname)

  ! Check gz exists
  inquire(file=trim(climfile)//'.gz', exist=file_exists)
  if (.not. file_exists) then
    forcingData => null()
    write (*, '("read_interpolatedCRU: ", a, " does not exist")') trim(climfile)//'.gz'
    return
  endif

  ! Unzip (keep original .gz) with explicit status
  call unzip_gzip_file(trim(climfile)//'.gz', ok=unzip_ok, quiet=.true.)
  if (.not. unzip_ok) then
    file_exists = .false.
    forcingData => null()
    write (*, '("read_interpolatedCRU: unzip failed for ", a)') trim(climfile)//'.gz'
    return
  endif

  ! Confirm unzipped CSV exists
  inquire(file=trim(climfile), exist=csv_exists)
  if (.not. csv_exists) then
    file_exists = .false.
    forcingData => null()
    write (*, '("read_interpolatedCRU: missing unzipped file ", a)') trim(climfile)
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

  ! Read in forcing data (use newunit to avoid unit collisions)
  open(newunit=u, file=climfile, status='old', action='read', iostat=istat2)
  if (istat2 /= 0) then
    file_exists = .false.
    forcingData => null()
    write (*, '("read_interpolatedCRU: cannot open ", a, ", iostat=", I0)') trim(climfile), istat2
    deallocate(input_data, timecols)
    return
  endif

  ! Header line
  read(u,'(a160)', iostat=istat2) commts
  if (istat2 /= 0) then
    file_exists = .false.
    forcingData => null()
    write (*, '("read_interpolatedCRU: failed reading header from ", a, ", iostat=", I0)') trim(climfile), istat2
    close(u)
    deallocate(input_data, timecols)
    return
  endif

  m = 0
  do while (m < datalines .and. istat2 == 0)
    read(u,*, iostat=istat2) (temp(n), n = 1,niterms)
    if (istat2 == 0) then
      m = m + 1
      input_data(:,m) = temp(:)
    endif
  enddo

  close(u)

  ! Remove unzipped CSV (best effort; keep .gz)
  command = 'rm ' // trim(climfile)
  call execute_command_line(command, exitstat=exitstat)
  ! Not fatal if rm fails

  ! Check the consistency between the file data lines and required
  if (m /= datalines) then
    write (*, '("In read_interpolatedCRU, File ",a," is shorter than needed: lines: ",I12)') trim(fname), m
    file_exists = .false.
    deallocate(input_data, timecols)
    forcingData => null()
    return
  endif

  ! Setup the time table
  m = 0
  do iyr = year0, year1
    do iday = 1, 365
      do ihour = 1, steps_per_day
        m = m + 1
        timecols(1,m) = iyr
        timecols(2,m) = iday
        timecols(3,m) = 24.0 * (ihour - 1)/steps_per_day
      enddo
    enddo
  enddo

  ! Put the data into forcing
  allocate(climateData(datalines))
  do i = 1, datalines
     climateData(i)%year      = int(timecols(1,i))
     climateData(i)%doy       = int(timecols(2,i))
     climateData(i)%PAR       = input_data(1,i)*2.0
     climateData(i)%radiation = input_data(1,i)
     climateData(i)%Tair      = input_data(2,i)
     climateData(i)%Tsoil     = input_data(2,i)*0.8 + 273.16*0.2
     climateData(i)%RH        = max(0.01, min(0.99, input_data(3,i) ))
     climateData(i)%rain      = input_data(4,i)
     climateData(i)%windU     = input_data(5,i)
     climateData(i)%P_air     = input_data(6,i)
     climateData(i)%CO2       = CO2_Hist(min(CO2Yrs, max(1, climateData(i)%year-1700+1)))
     climateData(i)%eCO2      = climateData(i)%CO2 + 200.
     climateData(i)%soilwater = 0.8
#ifdef WIEMIP_setting
     climateData(i)%N_input   = input_data(7,i) ! for WIEMIP N deposition, kgN m-2 yr-1
#else
     climateData(i)%N_input   = N_input
#endif
  enddo

  forcingData => climateData
  write(*,*)"forcing from interpolated: hours,days,years", datalines,days_data,yr_data

  ! Release memory
  deallocate(input_data, timecols)
end subroutine read_interpolatedCRU

!==============================================================
    subroutine unzip_gzip_file(filename_gz, ok, exitstat, out_file, quiet)
    character(len=*), intent(in)            :: filename_gz
    logical,          intent(out), optional :: ok
    integer,          intent(out), optional :: exitstat
    character(len=*), intent(out), optional :: out_file
    logical,          intent(in),  optional :: quiet

    character(len=600) :: command
    character(len=600) :: filename_out
    integer :: est
    logical :: q, gz_exists, out_exists
    integer :: L

    q = .false.
    if (present(quiet)) q = quiet

    if (present(ok)) ok = .false.
    if (present(exitstat)) exitstat = -999

    ! Check input exists
    inquire(file=trim(filename_gz), exist=gz_exists)
    if (.not. gz_exists) then
      if (.not. q) write(*,'("unzip_gzip_file: missing ",a)') trim(filename_gz)
      est = 2
      if (present(exitstat)) exitstat = est
      return
    endif

    ! Derive output filename by stripping trailing ".gz" if present
    filename_out = trim(filename_gz)
    L = len_trim(filename_out)
    if (L >= 3) then
      if (filename_out(L-2:L) == '.gz') filename_out = filename_out(1:L-3)
    endif
    if (present(out_file)) out_file = trim(filename_out)

    ! If output already exists, treat as success (idempotent)
    inquire(file=trim(filename_out), exist=out_exists)
    if (out_exists) then
      est = 0
      if (present(exitstat)) exitstat = est
      if (present(ok)) ok = .true.
      return
    endif

    ! Unzip: keep .gz (-k), force overwrite (-f)
    command = 'gunzip -kf ' // trim(filename_gz)
    call execute_command_line(command, exitstat=est)

    ! Verify output exists
    inquire(file=trim(filename_out), exist=out_exists)

    if (est /= 0 .or. .not. out_exists) then
      if (.not. q) then
        write(*,'("unzip_gzip_file: failed for ",a," exitstat=",I0)') trim(filename_gz), est
        if (.not. out_exists) write(*,'("unzip_gzip_file: expected output missing: ",a)') trim(filename_out)
      endif
      if (present(exitstat)) exitstat = est
      if (present(ok)) ok = .false.
      return
    endif

    if (present(exitstat)) exitstat = est
    if (present(ok)) ok = .true.
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
      stop
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
        stop
      endif
      if (any(s < 1)) then
        write(*,*) 'ERROR: invalid start (<1) reading ', trim(var_name), ' start=', s
        stop
      endif
      if (any(s + c - 1 > dlen)) then
        write(*,*) 'ERROR: hyperslab out of bounds reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' start/count=', s, c
        stop
      endif
      call check( nf90_get_var(ncid, varid, da, start=s, count=c) )
    else
      if (size(da,1) /= dlen(1) .or. size(da,2) /= dlen(2) .or. size(da,3) /= dlen(3)) then
        write(*,*) 'ERROR: full-field shape mismatch reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' array=', size(da,1), size(da,2), size(da,3)
        stop
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
      stop
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
        stop
      endif
      if (any(s < 1)) then
        write(*,*) 'ERROR: invalid start (<1) reading ', trim(var_name), ' start=', s
        stop
      endif
      if (any(s + c - 1 > dlen)) then
        write(*,*) 'ERROR: hyperslab out of bounds reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' start/count=', s, c
        stop
      endif
      call check( nf90_get_var(ncid, varid, da, start=s, count=c) )
    else
      if (size(da,1) /= dlen(1) .or. size(da,2) /= dlen(2)) then
        write(*,*) 'ERROR: full-field shape mismatch reading ', trim(var_name)
        write(*,*) ' dims=', dlen, ' array=', size(da,1), size(da,2)
        stop
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
      stop
    endif
    call check( nf90_inquire_dimension(ncid, dimids(1), len=dlen(1)) )

    if (size(da,1) /= dlen(1)) then
      write(*,*) 'ERROR: shape mismatch reading ', trim(var_name), ' dims=', dlen(1), ' array=', size(da,1)
      stop
    endif

    call check( nf90_get_var(ncid, varid, da) )
    call check( nf90_close(ncid) )
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
