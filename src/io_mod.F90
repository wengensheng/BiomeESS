module io_mod
! Data input and output
  use datatypes
  use model_utils

  implicit none
  private

! ------ public subroutines ---------
  public :: setup_forcingdata, setup_output_files
  public :: zip_output_files
  public :: hourly_diagnostics, daily_diagnostics, annual_diagnostics
  public :: flush_annual_diagnostics_buffers
#ifdef DO_ANIMAL
  public :: ani_annual_diagnostics
#endif

  integer, parameter :: annual_buf_max = 2048
  integer, parameter :: annual_line_len = 8192
  character(len=annual_line_len), allocatable :: annual_cohort_buf(:)
  character(len=annual_line_len), allocatable :: annual_cohort2_buf(:)
  character(len=annual_line_len), allocatable :: annual_patch_buf(:)
#ifdef DO_ANIMAL
  character(len=annual_line_len), allocatable :: annual_ani_buf(:)
#endif
  integer :: n_annual_cohort_buf  = 0
  integer :: n_annual_cohort2_buf = 0
  integer :: n_annual_patch_buf   = 0
#ifdef DO_ANIMAL
  integer :: n_annual_ani_buf     = 0
#endif

!---------------------------------
  contains

!====================== Buffered annual output ==========================
  subroutine init_annual_diagnostics_buffers()
    implicit none
    if (.not. allocated(annual_cohort_buf))  allocate(annual_cohort_buf(annual_buf_max))
    if (.not. allocated(annual_cohort2_buf)) allocate(annual_cohort2_buf(annual_buf_max))
    if (.not. allocated(annual_patch_buf))   allocate(annual_patch_buf(annual_buf_max))
#ifdef DO_ANIMAL
    if (.not. allocated(annual_ani_buf))     allocate(annual_ani_buf(annual_buf_max))
#endif

    n_annual_cohort_buf  = 0
    n_annual_cohort2_buf = 0
    n_annual_patch_buf   = 0
#ifdef DO_ANIMAL
    n_annual_ani_buf     = 0
#endif
  end subroutine init_annual_diagnostics_buffers

  subroutine append_annual_line(unit_no, line)
    implicit none
    integer, intent(in) :: unit_no
    character(len=*), intent(in) :: line

    if (unit_no == fno5) then
      call push_annual_line(annual_cohort_buf, n_annual_cohort_buf, fno5, line)
    elseif (unit_no == fno6) then
      call push_annual_line(annual_patch_buf, n_annual_patch_buf, fno6, line)
#ifdef DO_ANIMAL
    elseif (unit_no == fno7) then
      call push_annual_line(annual_ani_buf, n_annual_ani_buf, fno7, line)
#endif
#ifdef DroughtMIP
    elseif (unit_no == fno5+10) then
      call push_annual_line(annual_cohort2_buf, n_annual_cohort2_buf, fno5+10, line)
#endif
    else
      write(unit_no,'(A)') trim(line)
    endif
  end subroutine append_annual_line

  subroutine push_annual_line(buffer, nbuf, unit_no, line)
    implicit none
    character(len=annual_line_len), intent(inout) :: buffer(:)
    integer, intent(inout) :: nbuf
    integer, intent(in) :: unit_no
    character(len=*), intent(in) :: line

    if (nbuf >= size(buffer)) call flush_one_annual_buffer(buffer, nbuf, unit_no)
    nbuf = nbuf + 1
    buffer(nbuf) = line
  end subroutine push_annual_line

  subroutine flush_one_annual_buffer(buffer, nbuf, unit_no)
    implicit none
    character(len=annual_line_len), intent(inout) :: buffer(:)
    integer, intent(inout) :: nbuf
    integer, intent(in) :: unit_no
    integer :: i

    do i = 1, nbuf
      write(unit_no,'(A)') trim(buffer(i))
    enddo
    nbuf = 0
  end subroutine flush_one_annual_buffer

  subroutine flush_annual_diagnostics_buffers()
    implicit none
    logical :: is_open

    if (allocated(annual_cohort_buf)) then
      inquire(unit=fno5, opened=is_open)
      if (is_open) call flush_one_annual_buffer(annual_cohort_buf, n_annual_cohort_buf, fno5)
    endif

    if (allocated(annual_patch_buf)) then
      inquire(unit=fno6, opened=is_open)
      if (is_open) call flush_one_annual_buffer(annual_patch_buf, n_annual_patch_buf, fno6)
    endif

#ifdef DO_ANIMAL
    if (allocated(annual_ani_buf)) then
      inquire(unit=fno7, opened=is_open)
      if (is_open) call flush_one_annual_buffer(annual_ani_buf, n_annual_ani_buf, fno7)
    endif
#endif

#ifdef DroughtMIP
    if (allocated(annual_cohort2_buf)) then
      inquire(unit=fno5+10, opened=is_open)
      if (is_open) call flush_one_annual_buffer(annual_cohort2_buf, n_annual_cohort2_buf, fno5+10)
    endif
#endif
  end subroutine flush_annual_diagnostics_buffers

!============================= Subroutines ====================================
!====================== Read Forcing Data ===============================
! --------- Setup forcing data and step lenght ----------------------
  subroutine setup_forcingdata(climfile)
    implicit none
    character(len=*),intent(in) :: climfile

#ifdef DroughtPaleo
    call set_PaleoForcing(climfile,PaleoPfile,PaleoTfile,iDraw, &
    forcingData,datalines,data_days,data_yrs,step_hour)
#else
    call Read_ForcingData(climfile,forcingData,datalines,data_days,data_yrs,step_hour)
    !call read_NACPforcing(forcingData,datalines,data_days,data_yrs,step_hour)
#endif

    ! ------ Setup steps for model run ------
    steps_per_day = int(24.0/step_hour)
    dt_fast_yr    = step_hour/(365.0 * 24.0)
    step_seconds  = step_hour*3600.0
    write(*,*)'steps/day,dt_fast,s/step',steps_per_day,dt_fast_yr,step_seconds
    write(*,*)'Datalines,data_days,data_yrs,step_hour',datalines,data_days,data_yrs,step_hour
  end subroutine setup_forcingdata

!=============================================================================
! read in forcing data (Users need to write their own data input procedure)
  subroutine Read_ForcingData(fdata,forcingData,datalines,data_days,data_yrs,step_hour)
    implicit none
    character(len=*),intent(in) :: fdata
    type(climate_data_type),pointer,intent(inout) :: forcingData(:)
    integer,intent(inout) :: datalines,data_days,data_yrs
    real, intent(inout)   :: step_hour
    !------------local var -------------------
    type(climate_data_type), pointer :: climateData(:)
    character(len=80)  commts

#ifdef FACE_run
    integer, parameter :: niterms=26 ! 30 columns in FACEMDS-2
#else
    integer, parameter :: niterms=9 ! 9 columns in FACEMDS-1
#endif

    integer, allocatable :: HRMIN(:),doy_data(:),year_data(:)
    real,    allocatable :: DTIME(:),hour_data(:),input_data(:,:)
    integer :: yr,dy
    integer :: istat1,istat2,istat3
    integer :: ndays,nyear,totlines
    integer :: m,n,i

    ! Open forcing data
    climfile=trim(filepath_in)//trim(fdata)
    ! Check whether file exists
    inquire (file=climfile, iostat=istat1)
    if (istat1 /= 0) then
      write (*, '("Error: input file ", a, " does not exist")') climfile
      stop
    end if
    open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
    ! Skip 1 line of input met data file
    read(11,'(a160)') commts ! MDK data only has one line comments
#ifdef FACE_run
    read(11,'(a160)') commts ! Two lines of head in FACDMDS-2
#endif
    ! Count total lines
    totlines = 0  ! to record the lines in a file
    do
      read(11,*,IOSTAT=istat3)yr 
      if(istat3 < 0)exit
      totlines = totlines + 1
    enddo ! end of reading the forcing file
    write (*, '("Forcing file ", a, " total lines: ",I12)') trim(climfile),totlines

    ! Read forcing files
    rewind 11
    read(11,'(a160)') commts

    ! Allocate arrays for reading in data
    allocate(doy_data(totlines),year_data(totlines),hour_data(totlines))
    allocate(input_data(niterms,totlines))
#ifdef FACE_run
    allocate(DTIME(totlines),HRMIN(totlines))
    read(11,'(a160)') commts ! Two lines of head in FACDMDS-2
#endif
    ndays = 0 ! the total days in this data file
    nyear = 0 ! the total years of this data file
    dy    = -1  ! Initial value
    yr    = -1
    do m = 1, totlines
#ifdef FACE_run
      read(11,*,IOSTAT=istat3)year_data(m),DTIME(m),doy_data(m),HRMIN(m),   &
      (input_data(n,m),n=1,niterms)
#else
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      (input_data(n,m),n=1,niterms)
#endif
      ! Count days
      if(m > 1) then
        dy = doy_data(m-1)
        yr = year_data(m-1)
      endif
      if(dy /= doy_data(m)) ndays = ndays + 1
      if(yr /= year_data(m))nyear = nyear + 1
      !Remove -9999
      do n=1,niterms
        if(input_data(n,m)<-900.0 .and. m > 1)then
          input_data(n,m) = input_data(n,m-1) ! remove -9999
        endif
      enddo
    enddo
    ! Check fast time step
#ifdef FACE_run
    step_hour = (HRMIN(2) - HRMIN(1))/60.0
#else
    step_hour = hour_data(2) - hour_data(1)
#endif
    if (step_hour==1.0)then
      write(*,*)"the data freqency is hourly"
    elseif(step_hour==0.5)then
      write(*,*)"the data freqency is half hourly"
    else
      write(*,*)'hour data:',hour_data(1),hour_data(2),hour_data(3)
      write(*,*)"Please check time step!"
      stop
    endif

    ! Put the data into forcing
    allocate(climateData(totlines))
    do i=1,totlines
      climateData(i)%year      = year_data(i)          ! Year
      climateData(i)%doy       = doy_data(i)           ! day of the year
#ifdef FACE_run
      climateData(i)%PAR       = input_data(15,i)      ! umol/m2/s
      climateData(i)%radiation = input_data(15,i)/2.   ! W/m2, input_data(13,i), factor
      climateData(i)%Tair      = input_data(3,i)       ! air temperature, K
      climateData(i)%Tsoil     = input_data(3,i)       ! soil temperature, K
      climateData(i)%RH        = input_data(5,i)*0.01  ! relative humidity (0.xx)
      climateData(i)%rain      = input_data(1,i)       ! kgH2O m-2 s-1
      climateData(i)%windU     = input_data(11,i)      ! wind velocity (m s-1)
      climateData(i)%P_air     = input_data(19,i)      ! pa
      climateData(i)%CO2       = input_data(21,i) !ppm
      climateData(i)%eCO2      = input_data(22,i) !ppm
#else
      climateData(i)%PAR       = input_data(1,i)       ! umol/m2/s
      climateData(i)%radiation = input_data(2,i)       ! W/m2
      climateData(i)%Tair      = input_data(3,i) + 273.16  ! air temperature, K
      climateData(i)%Tsoil     = input_data(4,i) + 273.16  ! soil temperature, K
      climateData(i)%RH        = min(input_data(5,i),99.9) * 0.01    ! relative humidity (0.xx)
      climateData(i)%rain      = input_data(6,i)/(step_hour * 3600)! ! kgH2O m-2 s-1
      climateData(i)%windU     = input_data(7,i)        ! wind velocity (m s-1)
      climateData(i)%P_air     = input_data(8,i)        ! pa
      climateData(i)%CO2       = input_data(9,i)        !ppm
#endif
      climateData(i)%soilwater = 0.8                    ! soil moisture, vol/vol
      climateData(i)%N_input   = N_input                ! kgN m-2 yr-1

    enddo
    forcingData => climateData
    datalines = totlines
    data_days = ndays
    data_yrs  = nyear
    write(*,*)"runID:  ", runID
    write(*,*)"forcing: hours,days,years", datalines,data_days,data_yrs

    !Close opened file and release memory
    close(11)    ! close forcing file
    deallocate(doy_data,year_data,hour_data,input_data)
#ifdef FACE_run
    deallocate(DTIME,HRMIN)
#endif
  end subroutine Read_ForcingData

  !=============================================================
#ifdef Use_NACP_site_data
  ! for reading in NACP site synthesis forcing
  subroutine read_NACPforcing(forcingData,datalines,data_days,data_yrs,step_hour)
    implicit none
    type(climate_data_type),pointer,intent(inout) :: forcingData(:)
    integer,intent(inout) :: datalines,data_days,data_yrs
    real, intent(inout)   :: step_hour
    !------------local var -------------------
    type(climate_data_type), pointer :: climateData(:)
    character(len=80)  commts
    integer, parameter :: niterms=15       ! NACP site forcing
    integer, parameter :: ilines=22*366*48 ! the maxmum records
    integer,dimension(ilines) :: year_data, doy_data
    real,   dimension(ilines) :: hour_data
    real input_data(niterms,ilines)
    real inputstep
    integer :: istat1,istat2,istat3
    integer :: doy,idays
    integer :: i,j,k
    integer :: m,n

    climfile=trim(filepath_in)//trim(climfile)
    write(*,*)'inputfile: ',climfile
    ! open forcing data
    open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
    write(*,*)istat2
    ! skip 2 lines of input met data file
    read(11,'(a160)') commts
    read(11,'(a160)') commts
    m       = 0  ! to record the lines in a file
    idays   = 1  ! the total days in a data file
    data_yrs = 0 ! to record years of a dataset
    do    ! read forcing files
      m=m+1
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      (input_data(n,m),n=1,niterms)
      if(istat3<0)exit
      if(m == 1) then
        doy = doy_data(m)
      else
        doy = doy_data(m-1)
      endif
      if(doy /= doy_data(m)) idays = idays + 1
    enddo ! end of reading the forcing file

    step_hour = hour_data(2) - hour_data(1)
    write(*,*)"forcing",datalines,data_yrs,step_hour,dt_fast_yr
    if (step_hour==1.0)then
      write(*,*)"the data freqency is hourly"
    elseif(step_hour==0.5)then
      write(*,*)"the data freqency is half hourly"
    else
      write(*,*)"Please check time step!"
      stop
    endif
    close(11)    ! close forcing file
    ! Put the data into forcing
    datalines = m - 1
    data_days = idays
    data_yrs  = year_data(datalines-1) - year_data(1) + 1

    allocate(climateData(datalines))
    do i=1,datalines
      climateData(i)%year      = year_data(i)          ! Year
      climateData(i)%doy       = doy_data(i)           ! day of the year
      climateData(i)%PAR       = input_data(11,i)*2.0  ! umol/m2/s
      climateData(i)%radiation = input_data(11,i)      ! W/m2
      climateData(i)%Tair      = input_data(1,i)       ! air temperature, K
      climateData(i)%Tsoil     = input_data(1,i)       ! soil temperature, K
      climateData(i)%rain      = input_data(7,i)       ! kgH2O m-2 s-1
      climateData(i)%windU     = input_data(5,i)        ! wind velocity (m s-1)
      climateData(i)%P_air     = input_data(9,i)        ! pa
      climateData(i)%RH        = input_data(3,i)/mol_h2o*mol_air* & ! relative humidity (0.xx)
      climateData(i)%P_air/esat(climateData(i)%Tair-273.16)
      climateData(i)%CO2       = input_data(15,i) ! ppm
      climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
      climateData(i)%N_input   = N_input ! kgN m-2 yr-1
    enddo
    forcingData => climateData
    write(*,*)"runID:  ", runID
    write(*,*)"forcing: hours,days,years", datalines,data_days,data_yrs

  end subroutine read_NACPforcing
#endif

!========================================================================
#ifdef DroughtPaleo
  ! Set up forcing data with paleo precipitation and temperature (monthly)
  subroutine set_PaleoForcing(fdata,fPaleoP,fPaleoT,iDraw, &
    forcingData,datalines,data_days,data_yrs,step_hour)
    implicit none
    character(len=*),intent(in) :: fdata ! Base climate data
    character(len=*),intent(in) :: fPaleoP,fPaleoT ! Paleo inversion data
    integer,intent(in) :: iDraw
    type(climate_data_type),pointer,intent(inout) :: forcingData(:)
    integer,intent(inout) :: datalines,data_days,data_yrs
    real, intent(inout)   :: step_hour

    !------------local var -------------------
    integer, parameter :: N_draws = 1000
    integer, parameter :: N_months = 12
    integer, parameter :: PaleoYears  = 900 ! 901
    integer, parameter :: PaleoMonths = PaleoYears * 12 ! 10812
    integer, parameter :: MonthDays(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    character(len=160)  commts,PaleoPfile,PaleoTfile,fname3
    character(len=10)  mAbv,DrawID
    type(climate_data_type), pointer :: climateData(:)
    real, pointer :: monthlyP(:,:),monthlyT(:,:)
    real, dimension(PaleoMonths,N_draws) :: PaleoP, PaleoT
    real :: fPrcp, dTmp
    integer :: PaleoForcingLines
    integer :: Lines_skip = 3 + 4 ! three lines of comments and 4 lines of data, Sep - Dec
    integer :: istat1, istat2
    integer :: i,j,n
    real :: iYear
    integer :: iLine,iBase,iBY,iY,iM,iD ! Year, Month, Day, Hour

    ! Read in baseline forcing data (1901~1930, 30 years)
    call Read_ForcingData(fdata,forcingData,datalines,data_days,data_yrs,step_hour)
    ! Calculate monthely P and T
    allocate(monthlyP(data_yrs,12),monthlyT(data_yrs,12))
    monthlyP = 0.0
    monthlyT = 0.0
    iBase = 0
    do iY =1,data_yrs
      do iM=1,12
        n =  0
        do iD=1, MonthDays(iM) * int(24.0/step_hour)
          n = n + 1
          iBase = iBase + 1
          if(iBase > datalines) exit
          monthlyP(iY,iM) = monthlyP(iY,iM) + forcingData(iBase)%rain
          monthlyT(iY,iM) = monthlyT(iY,iM) + forcingData(iBase)%Tair
        enddo
        if(monthlyP(iY,iM)<1.0E-9)then ! Assign a very small value for zero rainfall month
          monthlyP(iY,iM) = 1.0E-9     ! just for put the paleo rainfall at the last hour
          forcingData(iBase)%rain = 1.0E-9
        endif
        monthlyP(iY,iM) = monthlyP(iY,iM) * (step_hour * 3600) ! Monthly total
        monthlyT(iY,iM) = monthlyT(iY,iM) /n - 273.16 ! K to C
      enddo
    enddo

    ! Read in Paleo precipitation and temperature data, monthly, 1001~1901
    PaleoPfile=trim(filepath_in)//trim(fPaleoP)
    PaleoTfile=trim(filepath_in)//trim(fPaleoT)
    inquire (file=PaleoPfile, iostat=istat1)
    if (istat1 /= 0) then
      write (*, '("Error: input file ", a, " does not exist")') PaleoPfile
      stop
    endif
    inquire (file=PaleoTfile, iostat=istat2)
    if (istat2 /= 0) then
      write (*, '("Error: input file ", a, " does not exist")') PaleoTfile
      stop
    endif
    open(21,file=PaleoPfile,status='old',ACTION='read',IOSTAT=istat2)
    open(22,file=PaleoTfile,status='old',ACTION='read',IOSTAT=istat2)
    do i=1,Lines_skip
      read(21,*) commts
      read(22,*) commts
    enddo

    do i=1,PaleoMonths
      read(21,*,IOSTAT=istat1)iYear, mAbv,(PaleoP(i,j),j=1,1000)
      read(22,*,IOSTAT=istat2)iYear, mAbv,(PaleoT(i,j),j=1,1000)
    enddo

    ! Replace base data's P and T
    PaleoForcingLines = INT(PaleoYears*365*24/step_hour)
    allocate(climateData(PaleoForcingLines))
    iBase = 0
    iLine = 0
    do iY =1, PaleoYears ! 901
      iBY = MOD(iY-1,data_yrs)+1 ! Corresponding base data year
      do iM=1,12
        ! Calculate ratios of Paleo P and T to the base data's
        fPrcp = PaleoP((iY-1)*12+iM,iDraw) / monthlyP(iBY,iM)
        dTmp  = PaleoT((iY-1)*12+iM,iDraw) - monthlyT(iBY,iM)
        do iD=1, MonthDays(iM) * int(24.0/step_hour)
          iBase = MOD(iLine,datalines) + 1
          iLine = iLine + 1
          climateData(iline)       = forcingData(iBase)
          climateData(iline)%rain  = forcingData(iBase)%rain * fPrcp
          climateData(iline)%Tair  = forcingData(iBase)%Tair + dTmp
          climateData(iline)%Tsoil = forcingData(iBase)%Tsoil+ dTmp
        enddo ! month hours
      enddo   ! Months
      !if(iBY==data_yrs)write(*,*)'iY,iLine:',iY,iLine
    enddo     ! years
    deallocate(monthlyP,monthlyT)
    deallocate(forcingdata)
    ! Update data array for model run
    forcingData => climateData
    datalines = iLine
    data_days = PaleoYears * 365
    data_yrs   = PaleoYears

    if(WriteForcing) then
      ! Write climateData to a csv file, for checking only
      write(DrawID, '(I0)')iDraw
      fname3 = trim(filepath_out)//trim(fPaleoP(1:3))//'_Hourly_'//trim(DrawID)//'.csv'
      open(15,file=trim(fname3))
      write(15,*)"YEAR,DOY,HoD,PAR,Swdown,Tair,Tsoil,RH,RAIN,WIND,PRESSURE,CO2"
      do i=1,PaleoForcingLines
        write(15,'(2(I4,","),12(E15.4,","))') &
        forcingData(i)%year, forcingData(i)%doy,       forcingData(i)%hod,  &
        forcingData(i)%PAR,  forcingData(i)%radiation, forcingData(i)%Tair, &
        forcingData(i)%Tsoil,forcingData(i)%RH,        forcingData(i)%rain, &
        forcingData(i)%windU,forcingData(i)%P_air,     forcingData(i)%CO2
      enddo
      close(15)
    endif

  end subroutine set_PaleoForcing
#endif

!==============================================================================
!======================= Setup output files ===================================
  subroutine setup_output_files()

    ! ----------Local vars ------------
    implicit none
    character(len=150) :: YearlyCohort2, DailyPatch2  ! For DroughtMIP only
    character(len=120) :: filesuffix, fpath
    character(len=6)   :: LonLat
    integer :: istat1, istat2, istat3
#ifdef Hydro_test
    integer :: k_h
    character(len=8) :: farea_hdr(Ysw_max)

    do k_h = 1, Ysw_max
      write(farea_hdr(k_h), '(a,i0)') 'farea', k_h
    enddo
#endif

    call init_annual_diagnostics_buffers()

    ! File path and names
    fpath = trim(filepath_out)
    filesuffix   = trim(runID) ! tag for simulation experiments
#ifdef GlobalRun
    write(LonLat, GridIDFMT) GridID
    filesuffix = trim(filesuffix)//trim(LonLat)
#endif
    file_out(1) = trim(fpath)//trim(filesuffix)//'_Cohort_hourly.csv'       ! hourly
    file_out(2) = trim(fpath)//trim(filesuffix)//'_Ecosystem_hourly.csv'    ! hourly
    file_out(3) = trim(fpath)//trim(filesuffix)//'_Cohort_daily.csv'        ! daily
    file_out(4) = trim(fpath)//trim(filesuffix)//'_Ecosystem_daily.csv'     ! Daily
    file_out(5) = trim(fpath)//trim(filesuffix)//'_Cohort_yearly.csv'       ! Yearly
    file_out(6) = trim(fpath)//trim(filesuffix)//'_Ecosystem_yearly.csv'    ! Yearly

#ifdef DroughtMIP
    ! For DroughtMIP
    YearlyCohort2 = trim(fpath)//trim(filesuffix)//'2_Cohort_yearly.csv'    ! Yearly
    DailyPatch2   = trim(fpath)//trim(filesuffix)//'2_Ecosystem_daily.csv'  ! Daily
#endif

    ! Open files
    if(outputhourly)then
      open(fno1,file=trim(file_out(1)),ACTION='write', IOSTAT=istat1)
      write(fno1,'(5(a8,","),30(a12,","))')'G'//LonLat, &       ! Hourly cohort
      'year','doy','hour','cID','sp','layer', &
      'density','dbh','height','Acrown',      &
      'bl','LAI','GPP', 'NPP', 'Transp',      &
#ifdef Hydro_test
      'Psi_L','Psi_W','W_lf','W_sw'
#else
      'W_supply','W_scale'
#endif

      open(fno2,file=trim(file_out(2)), ACTION='write', IOSTAT=istat1)
      write(fno2,'(5(a8,","),30(a12,","))')    &       ! Hourly tile
      'G'//LonLat,'year','doy','hour','rad',&
      'Tair','Prcp', 'GPP', 'Resp',         &
      'Transp','Evap','Runoff','Soilwater', &
      'wcl', 'psi_soil','k_soil',           &
      'bl','Psi_L','Psi_W','W_lf','W_sw','Transp'
    endif

    if(outputdaily)then
      open(fno3,file=trim(file_out(3)), ACTION='write', IOSTAT=istat2)
      if(istat2 /= 0)then
        write(*,*) 'fno3 open error. Stopped!'
        stop
      endif

      open(fno4,file=trim(file_out(4)),  ACTION='write', IOSTAT=istat2)
      if(istat2 /= 0)then
        write(*,*) 'fno4 open error. Stopped!'
        stop
      endif

      ! Write in file headers
      write(fno3,'(60(a8,","))')'yr'//LonLat,'doy',    &  ! Cohort daily
      'c_No','PFT','layer','Pheno','ndm','ncd',     &
      'density','Acrown','LAI','LeafAge',           &
      'gpp','resp','transp','NPPL','NPPR','NPPW',   &
      'NSC','seedC','leafC','rootC','SW-C','HW-C',  &
      'NSN','seedN','leafN','rootN','SW-N','HW-N',  &
      'W_LF','W_SW','W_HW',                         &
      'GDD','ALT','AWD'

      write(fno4,'(2(a8,","),55(a10,","))')'Yr'//LonLat, 'doy',   &  ! Tile daily
      'Tc','Prcp','Trsp','Evap','Roff','WaterS','thetaS',&
      'WC1_5','WC2_25','WC3_50','WC4_100','WC5_120',     &
      'LAI','GPP','Rauto','Rh',                          &
      'fineL', 'strucL', 'McrbC', 'fastSOC', 'slowSOC',  &
      'fineN', 'strucN', 'McrbN', 'fastSON', 'slowSON',  &
      'mineralN', 'N_uptk','CH4', 'dNorg','dNgas','dNmin' !,'Kappa'
    endif

    ! Open yearly output files
    open(fno5,file=trim(file_out(5)),ACTION='write', IOSTAT=istat3)
    if(istat3 /= 0)then
      write(*,*) 'fno5 open error. Stopped!'
      stop
    endif

    open(fno6,file=trim(file_out(6)), ACTION='write', IOSTAT=istat3)
    if(istat3 /= 0)then
      write(*,*) 'fno6 open error. Stopped!'
      stop
    endif

#ifdef DO_ANIMAL
    open(fno7, file=trim(fpath)//trim(filesuffix)//'_Animal_yearly.csv', &
         ACTION='write', IOSTAT=istat3)
    if(istat3 /= 0)then
      write(*,*) 'fno7 (Animal_yearly) open error. Stopped!'
      stop
    endif
    write(fno7,'(4(a8,","),12(a14,","))') &
      'G'//LonLat,'year','cohortNo','AFT',                        &
      'nindivs','age',                                            &
      'IntakePlant','IntakePrey','C_removed',                     &
      'FecalC','FecalN','CarcassC','CarcassN',                    &
      'mu_starve','deaths'
#endif

#ifdef DroughtMIP
    !For baseline runs
    open(fno4,file=trim(file_out(4)),  ACTION='write', IOSTAT=istat2)
    write(fno4,'(3(a5,","),55(a10,","))')'YEAR', 'Month','DAY',   &  ! Daily tile, 'tile',
    'GPP','NPP','ET','LAI','LFLIT','SW1','SW2','SW3','SW4'

    ! for scenario runs
    open(fno4+10,file=trim(DailyPatch2),  ACTION='write', IOSTAT=istat2)
    write(fno4+10,'(2(a5,","),55(a10,","))')'YEAR', 'Month','DAY',   &  ! Daily tile, 'tile',
    'GPP','NPP','ET','LAI','LFLIT','SW1','SW2','SW3','SW4'

    write(fno5,'(3(a5,","),55(a10,","))')'YEAR', 'SP','ID',   &
    'NLIVE','DBH','HT','TB','AGB','WD','SLA','Acrown'

    open(fno5+10,file=trim(YearlyCohort2),ACTION='write', IOSTAT=istat3)
    write(fno5+10,'(3(a5,","),55(a10,","))')'YEAR', 'SP','ID',   &
    'NLIVE','DBH','HT','TB','AGB','WD','SLA','Acrown'

#elif DBEN_run
    write(fno5,'(4(a5,","),40(a9,","))')'tile',       &    ! Yearly cohort
    'yr','cNo.','cID','PFT','Woody','Layer',          &
    'Density','f_L','dbh','height','Acrown','Aleaf',  &
    'bl','br','bSW','bHW','seed','nsc',               &
    'GPP','NPP','dDBH','dBA','dCA',                   &
    'Gtree','f_sd','f_lf','f_fr','f_wd','mu'

    write(fno6,'(1(a5,","),80(a12,","))')'year',           &  ! Yearly tile
    'CAI', 'LAI', 'GPP', 'Rauto', 'Rh',               &
    'rain','SoilWater','Transp','Evap','Runoff',      &
    'plantC', 'soilC', 'plantN', 'soilN',             &
    'leafC', 'rootC', 'swC', 'hwC', 'SeedC', 'NSC',   &
    'leafN', 'rootN', 'swN', 'hwN', 'SeedN', 'NSN',   &
    'fineL', 'strucL', 'McrbC', 'fastSOC', 'slowSOC', &
    'fineN', 'strucN', 'McrbN', 'fastSON', 'slowSON','CH4'

#elif FACE_run
    write(fno5,'(4(a5,","),40(a7,","))')              &    ! Yearly cohort
    'yr','cNo.','PFT','layer','f_L','N_ha','mu',      &
    'dD','dCA','dbh','ht','Acrown','Aleaf',           &
    'bl','br','bSW','bHW','seed','nsc',               &
    'N_lf','N_fr','N_SW','N_HW','N_sd','NSN','N_up',  &
    'GPP','NPP','NPPl','NPPfr','NPPw','Trsp',         &
    'demandW','Asap','Ktree','treeHU','treeW0'
    write(fno6,'(1(a5,","),80(a12,","))')'year',      &  ! Yearly tile
    'CAI', 'LAI', 'GPP', 'Rauto', 'Rh',               &
    'rain','SoilWater','Transp','Evap','Runoff',      &
    'plantC', 'soilC', 'plantN', 'soilN',             &
    'leafC', 'rootC', 'swC', 'hwC', 'SeedC', 'NSC',   &
    'leafN', 'rootN', 'swN', 'hwN', 'SeedN', 'NSN',   &
    'fineL', 'strucL', 'McrbC', 'fastSOC', 'slowSOC', &
    'fineN', 'strucN', 'McrbN', 'fastSON', 'slowSON', &
    'mineralN','Nm_SL', 'N_up', 'Nm_FR', 'N_loss',    &
    'CO2','CH4'

#else
    write(fno5,'(4(a8,","),380(a8,","))')             &    ! Yearly cohort
    'G'//LonLat,'yr','cNo.','cID', 'PFT','layer',     &
    'N_ha','f_L','dD','dBA','dCA','dbh','ht','Acrown',&
    'Aleaf','bl','br','bSW','bHW','seed','nsc','NSN', &
    'GPP','NPP','Gtree','f_sd','f_lf','f_fr','f_wd',  &
    'mu','Trsp','dmdW','Nup','Nfix','gddON','TcOFF',  &
    'Atrunk','Asap','Ktree','treeHU',                 &
#ifdef Hydro_test
    'treeW0', (farea_hdr(k_h), k_h=1,Ysw_max)
#else
    'treeW0'
#endif

    write(fno6,'(1(a8,","),80(a12,","))')'G'//LonLat,'year',     &  ! Yearly tile
    'CAI','LAI','GPP', 'Rauto', 'Rh', 'burned',                  &
    'Tmp','rain','SoilWater','Transp','Evap','Runoff',           &
    'plantC', 'soilC', 'plantN', 'soilN',                        &
    'NSC', 'SeedC', 'leafC', 'rootC', 'swC', 'hwC',              &
    'NSN', 'SeedN', 'leafN', 'rootN', 'swN', 'hwN',              &
    'fineL', 'strucL', 'McrbC', 'fastSOC', 'slowSOC',            &
    'fineN', 'strucN', 'McrbN', 'fastSON', 'slowSON','mineralN', &
    'WC1_5','WC2_25','WC3_50','WC4_100','WC5_120',               &
    'N_fxed','N_uptk','Nm_SL','Nm_FR','dNorg','dNgas','dNmin',   &
    'TreeCA','GrassCA','BMgrass','PET','Frisk','Pfire','CH4',    &
    'HarvC','HarvN'

#endif

  end subroutine setup_output_files

!================================================
  subroutine zip_output_files()
    implicit none
    character(len=256) :: command
    integer :: N_files, i, iostat
    integer :: idx(6) = [6,5,4,3,2,1]

    if (          outputhourly .and. outputdaily )then
      N_files = 6
    elseif((.not. outputhourly).and. outputdaily )then
      N_files = 4
    elseif((.not. outputhourly).and.(.not.outputdaily))then
      N_files = 2
    elseif (      outputhourly .and. (.not. outputdaily))then
      N_files = 4
      idx(3:4) = [1,2]
    endif
    do i = 1, N_files
      command = 'gzip -f ' // trim(file_out(idx(i)))
      call execute_command_line(command, exitstat=iostat)
      if (iostat /= 0) then
        print *, 'Error zipping: ', trim(file_out(idx(i))), ' (Exit status: ', iostat, ')'
      end if
    enddo
  end subroutine zip_output_files

!=========================================================================
!================= Diagnostics============================================
! Hourly fluxes sum to daily
  subroutine hourly_diagnostics(vegn,forcing,iyears,idoy,ihour,iday)
    implicit none
    type(vegn_tile_type), intent(inout) :: vegn
    type(climate_data_type),intent(in):: forcing
    integer, intent(in) :: iyears,idoy,ihour,iday

    !-------local var ------
    type(cohort_type), pointer :: cc    ! current cohort
    integer :: i

    !! Output horly diagnostics
    If(outputhourly .and. iday > totdays-366*5 ) then !  .and. ihour==12
      !write(fno1,'(4(I8,","))')vegn%n_cohorts
      do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        write(fno1,'(7(I8,","),40(F12.4,","))')vegn%tileID, &
        iyears,idoy,ihour,cc%ccID,cc%species,cc%layer,    &
        cc%nindivs*10000,cc%dbh,cc%height,cc%Acrown,      &
        cc%bl,cc%LAI,cc%gpp,cc%npp,cc%transp,             &
#ifdef Hydro_test
        cc%psi_leaf,cc%psi_stem,cc%W_lf,cc%W_sw
#else
        cc%W_supply,cc%W_scale
#endif

      enddo
      ! Hourly tile
      associate ( cc1 => vegn%cohorts(1))
        write(fno2,'(4(I5,","),60(E12.4,","))') vegn%tileID,   &
        iyears,idoy,ihour,forcing%radiation,forcing%Tair,    &
        forcing%rain,vegn%GPP,vegn%resp,vegn%transp,         &
        vegn%evap,vegn%runoff,vegn%soilwater,                &
        vegn%wcl(2),vegn%psi_soil(2),vegn%K_soil(2),         &
        cc1%bl,cc1%psi_leaf,cc1%psi_stem,cc1%W_lf,         &
        cc1%W_sw,cc1%transp
      end associate
    endif

  end subroutine hourly_diagnostics

!========================================================================================
  subroutine daily_diagnostics(vegn,iyears,idoy,iday,MonthDays)
    implicit none
    type(vegn_tile_type), intent(inout) :: vegn
    integer, intent(in) :: iyears,idoy,iday
    integer, intent(in) :: MonthDays(0:12)
    !-------local var ------
    type(cohort_type), pointer :: cc    ! current cohort
    integer :: i,j
#ifdef DroughtMIP
    integer :: f_eco,iyr_out
    integer :: iMonth, iDate

    ! Output daily cohorts
    if(iyears > 900)then
      !Write to two files
      if (iyears <= 1000) then
        f_eco = fno4
        iyr_out = iyears - 900
      else
        f_eco = fno4 + 10
        iyr_out = iyears - 1000
      endif

      !Convert doy to Month and Date
      !MonthDays(0:12): (/0,31,59,90,120,151,181,212,243,273,304,334,365)
      do i=1,12
        if(idoy <= MonthDays(i))then
          iMonth = i
          iDate  = idoy - MonthDays(i-1)
          exit
        endif
      enddo

      !! Tile daily
      write(f_eco,'(3(I5,","),65(F12.4,","))')iyr_out,iMonth,iDate,          &
      vegn%dailyGPP*1000.,vegn%dailyNPP*1000.,vegn%dailyTrsp+vegn%dailyEvap, &
      vegn%LAI,vegn%dailyLFLIT*1000., (vegn%wcl(i),i=2,5)
    endif

#elif DroughtFMT
    if(outputdaily.and. iday > skipped_days)then      !! Tile daily
      write(fno4,'(2(I5,","),70(E12.6,","))')iyears,idoy,            &
      vegn%tc_pheno, vegn%dailyPrcp,vegn%dailyTrsp,vegn%dailyEvap,   &
      vegn%dailyRoff,vegn%SoilWater,vegn%thetaS,(vegn%wcl(j),j=1,5), &
      vegn%LAI,vegn%dailyGPP,vegn%dailyResp,vegn%dailyRh,vegn%dailyCH4
    endif
#else
    if(outputdaily .and. iday > skipped_days)then
      !write(fno3,'(3(I6,","))')iyears, idoy,vegn%n_cohorts
      !! Cohort daily
      do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        write(fno3,'(8(I5,","),60(E12.6,","))')iyears,idoy,i, &
        cc%species,cc%layer,cc%status,cc%ndm,cc%ncd,     &
        cc%nindivs*10000.,cc%Acrown,cc%LAI,cc%leafage,   &
        cc%dailygpp,cc%dailyresp,cc%dailytrsp,           &
        cc%NPPleaf,cc%NPProot,cc%NPPwood,                &
        cc%NSC,cc%seedC,cc%bl,cc%br,cc%bsw,cc%bHW,       &
        cc%NSN*1000,cc%seedN*1000, cc%leafN*1000,        &
        cc%rootN*1000,cc%swN*1000,cc%hwN*1000,           &
        cc%W_lf,cc%W_sw,cc%W_hw,cc%gdd,cc%ALT,cc%AWD
      enddo
      !! Tile daily
      write(fno4,'(2(I5,","),70(E12.6,","))')iyears,idoy,      &
      vegn%Tc_daily, vegn%dailyPrcp,vegn%dailyTrsp,            &
      vegn%dailyEvap,vegn%dailyRoff,                           &
      vegn%SoilWater,vegn%thetaS,(vegn%wcl(j),j=1,5),          &
      vegn%LAI,vegn%dailyGPP, vegn%dailyResp, vegn%dailyRh,    &
      (vegn%SOC(j),j=1,5), (vegn%SON(j)*1000,j=1,5),           &
      vegn%mineralN*1000,vegn%dailyNup*1000, vegn%dailyCH4,    &
      vegn%dNorg_Daily*1000, vegn%dNgas_Daily*1000, vegn%dNmin_Daily*1000 !,vegn%kp(1)
    endif
#endif
  end subroutine daily_diagnostics

!==================================================================================================
  subroutine annual_diagnostics(vegn, iyears)
    implicit none
    type(vegn_tile_type), intent(inout) :: vegn
    integer, intent(in) :: iyears

    ! --------local var --------
    type(cohort_type), pointer :: cc
    character(len=annual_line_len) :: line
    real :: treeG, fseed, fleaf, froot,fwood
    real :: dDBH, dBA, dCA
    real :: plantC, plantN, soilC, soilN, N_loss_yr
    integer :: i,j,iyr_out
#ifdef DroughtMIP
    integer :: f_cht,yr_Eq,yr_Sc
#endif

    ! Yearly mean temperature
    vegn%YearlyTmp = vegn%YearlyTmp/365.0
    ! Max LAI
    vegn%LAImax = 0.0
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      vegn%LAImax = vegn%LAImax + cc%Aleafmax * cc%nindivs
    enddo
#ifdef ScreenOutput
    write(*,'(2(I6,","),3(F9.3,","))')iyears,vegn%n_cohorts
    write(*,'(3(a4,","),30(a9,","))')'cc','PFT','L',      &
    'n','f_CA','dD','DBH','NSC','Atrunk','Asap','Ktree', &
    'GPP','mu','W_scale','treeHU','treeW0'
#endif
    ! Cohotrs ouput
    iyr_out = iyears-yr_ResetVeg+30
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      associate ( sp => spdata(cc%species))
        treeG = MAX(1.0E-6, cc%seedC + cc%NPPleaf + cc%NPProot + cc%NPPwood)
        fseed = cc%seedC/treeG
        fleaf = cc%NPPleaf/treeG
        froot = cc%NPProot/treeG
        fwood = cc%NPPwood/treeG
        dDBH  = (cc%DBH - cc%DBH_ys) * 1000.0 ! mm
        dBA   = 3.1415926 * (cc%DBH**2 - cc%DBH_ys**2)/4.0
        dCA   = cc%Acrown - DBH2CA(cc%DBH_ys,cc%species)

#ifdef DroughtMIP
        yr_Sc = yr_Baseline
        yr_Eq = yr_Sc - yr_Baseline ! 100
        if(iyears > yr_Eq)then
          if (iyears <= yr_Sc) then
            f_cht = fno5
            iyr_out = iyears - yr_Eq
          else
            f_cht = fno5 + 10
            iyr_out = iyears - yr_Sc
          endif
          write(line,'(3(I8,","),300(E15.4,","))')               &
          iyr_out,cc%species,i, cc%nindivs*10000*(1.0-cc%mu),    &
          cc%dbh*100.,cc%height,TreeTotalC(cc),TreeTotalC(cc)*0.7, &
          2.0*sp%rho_wood,1.0/(2.0*sp%LMA), cc%Acrown
          call append_annual_line(f_cht, line)

        endif

#elif DBEN_run
        if(iyr_out > 0) then
          write(line,'(7(I8,","),300(E15.4,","))')vegn%tileID, &
          iyr_out,i,cc%ccID,cc%species,sp%lifeform,    &
          cc%layer,cc%nindivs*10000,cc%layerfrac,      &
          cc%dbh,cc%height,cc%Acrown,cc%Aleafmax,      &
          cc%bl,cc%br,cc%bsw,cc%bHW,cc%seedC,cc%nsc,   &
          cc%annualGPP,cc%annualNPP,dDBH,dBA,dCA,      &
          treeG,fseed,fleaf,froot,fwood,cc%mu
          call append_annual_line(fno5, line)
        endif
#elif FACE_run
        write(line,'(4(I8,","),300(E15.6,","))')iyears,i,   &
        cc%species,cc%layer,cc%layerfrac,cc%nindivs*10000,&
        cc%mu,dDBH,dCA,cc%dbh,cc%height,cc%Acrown,        &
        cc%Aleafmax,cc%bl,cc%br,cc%bsw,cc%bHW,cc%seedC,   &
        cc%nsc,cc%leafN*1000,cc%rootN*1000,cc%swN*1000,   &
        cc%hwN*1000,cc%seedN*1000, cc%NSN*1000,           &
        cc%NupYr*1000,cc%annualGPP,cc%annualNPP,          &
        cc%NPPleaf,cc%NPProot,cc%NPPwood,cc%annualTrsp,   &
        cc%totDemand,cc%Asap,cc%Ktrunk,cc%treeHU,cc%treeW0
        call append_annual_line(fno5, line)

#else
        write(line,'(6(I8,","),300(E15.6,","))')vegn%tileID, &
        iyears,i,cc%ccID,cc%species,cc%layer,            &
        cc%nindivs*10000,cc%layerfrac,dDBH,dBA,dCA,      &
        cc%dbh,cc%height,cc%Acrown,cc%Aleafmax,cc%bl,    &
        cc%br,cc%bsw,cc%bHW,cc%seedC,cc%nsc,cc%NSN,      &
        cc%annualGPP,cc%annualNPP,treeG,fseed,fleaf,     &
        froot,fwood,cc%mu,cc%annualTrsp,cc%totDemand,    &
        cc%NupYr,cc%NfixedYr,cc%gdd_ON,cc%Tc_OFF,        &
        cc%Atrunk,cc%Asap,cc%Ktrunk,cc%treeHU,           &
#ifdef Hydro_test
        cc%treeW0,(cc%farea(j),j=1,Ysw_max)
#else
        cc%treeW0
#endif
        call append_annual_line(fno5, line)

#endif

#ifdef ScreenOutput
        ! Screen output
        write(*,'(3(I4,","),1(F9.1,","),10(F9.3,","),10(F9.1,","))') &
        i,cc%species,cc%layer, &
        cc%nindivs*10000,cc%layerfrac,dDBH,cc%dbh,cc%nsc, &
        cc%Atrunk,cc%Asap,cc%Ktrunk,cc%annualGPP,cc%mu,   &
        cc%annualTrsp/cc%totDemand,cc%treeHU,cc%treeW0
#endif

      end associate
    enddo

    ! tile pools output

    if(iyr_out > 0) then
      call vegn_sum_tile(vegn)
      plantC = vegn%NSC + vegn%SeedC + vegn%leafC + vegn%rootC +  &
               vegn%SwC + vegn%HwC
      plantN = vegn%NSN + vegn%SeedN + vegn%leafN + vegn%rootN +  &
               vegn%SwN + vegn%HwN
      
      soilC  = sum(vegn%SOC(:))
      soilN  = sum(vegn%SON(:)) + vegn%mineralN

      N_loss_yr = (vegn%dNorg_Yr + vegn%dNgas_Yr + vegn%dNmin_Yr)*1000.
#ifdef FACE_run
      write(line,'(1(I5,","),85(E15.6,","))') iyears, &
      vegn%CAI,vegn%LAImax,vegn%annualGPP,vegn%annualResp,vegn%annualRh,  &
      vegn%annualPrcp, vegn%SoilWater, vegn%annualTrsp, vegn%annualEvap,  &
      vegn%annualRoff, plantC, soilC, plantN*1000, soilN*1000,            &
      vegn%leafC, vegn%rootC, vegn%SwC, vegn%HwC, vegn%SeedC,             &
      vegn%NSC, vegn%leafN*1000,vegn%rootN*1000,vegn%SwN*1000,            &
      vegn%HwN*1000, vegn%SeedN*1000, vegn%NSN*1000,                      &
      (vegn%SOC(j),j=1,5), (vegn%SON(j)*1000,j=1,5),                      &
      vegn%mineralN*1000, vegn%annualN*1000, vegn%NupYr*1000,             &
      vegn%Nm_Fire*1000, N_loss_yr, vegn%CO2_c,vegn%annualCH4
      call append_annual_line(fno6, line)
#elif DroughtMIP
      if (iyears > yr_Eq) then
        write(line,'(2(I5,","),80(E15.6,","))')&
        vegn%tileID,iyears - yr_Sc,vegn%CAI,vegn%LAI,                   &
        vegn%annualGPP,vegn%annualResp,vegn%annualRh,vegn%C_burned,     &
        vegn%annualPrcp,vegn%SoilWater,vegn%annualTrsp,vegn%annualEvap, &
        vegn%annualRoff,plantC,soilC,plantN*1000,soilN*1000,vegn%NSC,   &
        vegn%SeedC,vegn%leafC,vegn%rootC,vegn%SwC,vegn%HwC,             &
        vegn%NSN*1000,vegn%SeedN*1000,vegn%leafN*1000,vegn%rootN*1000,  &
        vegn%SwN*1000,vegn%HwN*1000,(vegn%SOC(j),j=1,5),                &
        (vegn%SON(j)*1000,j=1,5),vegn%mineralN*1000,vegn%annualCH4,     &
        (vegn%wcl(j),j=1,soil_L)
        call append_annual_line(fno6, line)
      endif

#else
      write(line,'(2(I5,","),120(E15.6,","))')  &
      vegn%tileID,iyears,vegn%CAI,vegn%LAI,vegn%annualGPP,            &
      vegn%annualResp,vegn%annualRh,vegn%C_burned,vegn%YearlyTmp,     &
      vegn%annualPrcp,vegn%SoilWater,vegn%annualTrsp,vegn%annualEvap, &
      vegn%annualRoff,plantC,soilC,plantN*1000,soilN*1000,vegn%NSC,   &
      vegn%SeedC,vegn%leafC,vegn%rootC,vegn%SwC,vegn%HwC,             &
      vegn%NSN*1000,vegn%SeedN*1000,vegn%leafN*1000,vegn%rootN*1000,  &
      vegn%SwN*1000,vegn%HwN*1000,(vegn%SOC(j),j=1,5),                &
      (vegn%SON(j)*1000,j=1,5),vegn%mineralN*1000,                    &
      (vegn%wcl(j),j=1,soil_L),vegn%NfixedYr*1000,vegn%NupYr*1000,    &
      vegn%Nm_Soil*1000,vegn%Nm_Fire*1000,                            &
      vegn%dNorg_Yr*1000, vegn%dNgas_Yr*1000, vegn%dNmin_Yr*1000,     &
      vegn%TreeCA,vegn%GrassCA,vegn%GrassBM,vegn%annualPET,           &
      vegn%Frisk,vegn%Pfire,vegn%annualCH4,vegn%HarvestC,vegn%HarvestN
      call append_annual_line(fno6, line)
#endif
    endif
  end subroutine annual_diagnostics

!=============================================================================
#ifdef DO_ANIMAL
  subroutine ani_annual_diagnostics(vegn, iyears)
    ! Write one CSV row per animal cohort to fno7 (_Animal_yearly.csv),
    ! then reset annual accumulators. Follows the pattern of annual_diagnostics.
    implicit none
    type(vegn_tile_type), intent(inout) :: vegn
    integer,              intent(in)    :: iyears

    integer :: i
    character(len=annual_line_len) :: line

    do i = 1, vegn%n_ani_cohorts
      associate (ac => vegn%ani_cohorts(i))
        write(line,'(4(I8,","),11(E15.6,","))')  &
          vegn%tileID, iyears, i, ac%aft,        &
          ac%nindivs, ac%age,                    &
          ac%annualIntakePlant,                  &
          ac%annualIntakePrey,                   &
          ac%annualC_removed,                    &
          ac%annualC_feces,  ac%annualN_feces,   &
          ac%annualC_carcass,ac%annualN_carcass, &
          ac%mu_starve,      ac%deaths
        call append_annual_line(fno7, line)

#ifdef ScreenOutput
        write(*,'(A,I3,A,I2,2(A,ES11.3),2(A,ES11.3))') &
          '  AFT cohort', i, '  AFT=', ac%aft,         &
          '  nindivs=', ac%nindivs, '  age=', ac%age,  &
          '  PlantIntake=', ac%annualIntakePlant,      &
          '  PreyIntake=', ac%annualIntakePrey
#endif

        ! Reset annual accumulators
        ac%annualIntakePlant = 0.0
        ac%annualIntakePrey  = 0.0
        ac%annualC_removed   = 0.0
        ac%annualC_feces     = 0.0
        ac%annualN_feces     = 0.0
        ac%annualC_carcass   = 0.0
        ac%annualN_carcass   = 0.0
      end associate
    end do

  end subroutine ani_annual_diagnostics
#endif /* DO_ANIMAL */

!================================================
end module io_mod
