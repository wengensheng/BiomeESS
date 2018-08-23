!========================================================================
!==== Biome Ecological Strategy Simulator (BiomeESS) ====================
!============   Main program   ==========================================
!=============== 10-21-2017 =============================================
!========================================================================
!
! This work was financially supported by US Forest Service and Princeton
! Environment Institute. The technical details of this model are in:
!
! Weng, E. S., Farrior, C. E., Dybzinski, R., Pacala, S. W., 2017.
! Predicting vegetation type through physiological and environmental 
! interactions with leaf traits: evergreen and deciduous forests in an 
! earth system modeling framework. Global Change Biology, 
! doi: 10.1111/gcb.13542.
!
! Weng, E. S., Malyshev, S., Lichstein, J. W., Farrior, C. E., 
! Dybzinski, R., Zhang, T., Shevliakova, E., Pacala, S. W., 2015. 
! Scaling from individual trees to forests in an Earth system modeling 
! framework using a mathematically tractable model of height-structured 
! competition. Biogeosciences, 12: 2655–2694, doi:10.5194/bg-12-2655-2015.
!
!
! Contact Ensheng Weng (wengensheng@gmail.com) for qeustions.
!                      (02/03/2017)
!
!------------------------------------------------------------------------
!
! This simulator can simulate evolutionarily stable strategy (ESS) of LMA
! and reproduce the forest succession patterns. But, since it
! does not include the models of photosynthesis, leaf stomatal
! conductance, transpiration, soil water dynamics, and energy balance, it 
! cannot simulate the ESS of allocation as reported in Weng et al. 2015 
! Biogeosciences.
!
! Processes included in this simulator are:
!     photosynthesis, transpiration, plant respiration
!     soil respraition,soil water dynamics
!     Phenology
!     Plant growth: Allometry and allocation
!     Reproduction
!     Mortality
!     Population dynamics
!     Soil C-N dynamics
!
!
!----------------------------------------
! Subroutine call structure:

!----- END -----------------------------------------------------------
!

program BiomeESS
   use datatypes
   use esdvm
   use soil_mod
   implicit none
   type(vegn_tile_type),  pointer :: vegn
   type(soil_tile_type),  pointer :: soil
   type(cohort_type),     pointer :: cp,cc

   integer,parameter :: rand_seed = 86456
   integer,parameter :: totalyears = 10
   integer,parameter :: nCohorts = 1
   integer :: datalines ! the total lines in forcing data file
   integer :: yr_data   ! Years of the forcing data
   integer :: days_data ! days of the forcing data
   integer :: steps_per_day ! 24 or 48
   real    :: timestep  ! hour, Time step of forcing data, usually hourly (1.0)
   real    :: tsoil, soil_theta
   real    :: NPPtree,fseed, fleaf, froot, fwood ! for output
   real    :: dDBH ! yearly growth of DBH, mm
   real    :: plantC,plantN, soilC, soilN
   real    :: dSlowSOM  ! for multiple tests only
   character(len=150) :: plantcohorts,plantCNpools,soilCNpools,allpools,faststepfluxes  ! output file names
   logical :: new_annual_cycle = .False.
   logical :: switch = .True.
   integer :: istat1,istat2,istat3
   integer :: year0, year1, iyears
   integer :: fno1,fno2,fno3,fno4,fno5 ! output files
   integer :: totyears, totdays
   integer :: i, j, k, idays, idoy
   integer :: simu_steps,idata
   character(len=50) :: filepath_out,filesuffix
   character(len=50) :: parameterfile(10),chaSOM(10)
   character(len=50) :: namelistfile = 'parameters_WC_biodiversity.nml' ! 'parameters_CN.nml'
   ! 'parameters_Allocation.nml' !'parameters_Konza.nml' !
   !

   ! create output files
   filepath_out='output/'
   filesuffix  = 'test.csv' ! tag for simulation experiments
   plantcohorts = trim(filepath_out)//'Annual_cohorts'//trim(filesuffix)
   plantCNpools = trim(filepath_out)//'Cohorts_daily'//trim(filesuffix)  ! daily
   soilCNpools  = trim(filepath_out)//'Ecosystem_daily'//trim(filesuffix)
   allpools     = trim(filepath_out)//'Ecosystem_yearly'//trim(filesuffix)
   faststepfluxes = trim(filepath_out)//'PhotosynthesisDynamics'//trim(filesuffix) ! hourly

   fno1=91; fno2=101; fno3=102; fno4=103; fno5=104
   open(fno1, file=trim(faststepfluxes),ACTION='write', IOSTAT=istat1)
   open(fno2,file=trim(plantcohorts),   ACTION='write', IOSTAT=istat1)
   open(fno3,file=trim(plantCNpools),   ACTION='write', IOSTAT=istat2)
   open(fno4,file=trim(soilCNpools),    ACTION='write', IOSTAT=istat3)
   open(fno5,file=trim(allpools),       ACTION='write', IOSTAT=istat3)
   ! head
   write(fno1,'(5(a8,","),25(a12,","))')      &
        'year','doy','hour','rad',            &
        'Tair','Prcp', 'GPP', 'Resp',         &
        'Transp','Evap','Runoff','Soilwater', &
        'wcl','FLDCAP','WILTPT'
   write(fno2,'(3(a5,","),25(a9,","))')            &
        'cID','PFT','layer','density', 'f_layer',  &
        'dDBH','dbh','height','Acrown',            &
        'wood','nsc', 'NSN','NPPtr','seed',        &
        'NPPL','NPPR','NPPW','GPP-yr','NPP-yr',    &
        'N_uptk','N_fix','maxLAI'

   write(fno3,'(5(a5,","),25(a8,","))')              &
        'year','doy','hour','cID','PFT',             &
        'layer','density', 'f_layer', 'LAI',         &
        'gpp','resp','transp',                       &
        'NSC','seedC','leafC','rootC','SW-C','HW-C', &
        'NSN','seedN','leafN','rootN','SW-N','HW-N'

   write(fno4,'(2(a5,","),55(a10,","))')  'year','doy',    &
        'Tc','Prcp', 'totWs',  'Trsp', 'Evap','Runoff',    &
        'ws1','ws2','ws3', 'LAI','GPP', 'Rauto', 'Rh',     &
        'NSC','seedC','leafC','rootC','SW-C','HW-C',       &
        'NSN','seedN','leafN','rootN','SW-N','HW-N',       &
        'McrbC', 'fastSOM',   'slowSOM',                   &
        'McrbN', 'fastSoilN', 'slowSoilN',                 &
        'mineralN', 'N_uptk'

   write(fno5,'(1(a5,","),80(a12,","))')  'year',              &
        'CAI','LAI','GPP', 'Rauto',   'Rh',                    &
        'rain','SiolWater','Transp','Evap','Runoff',           &
        'plantC','soilC',    'plantN', 'soilN','totN',         &
        'NSC', 'SeedC', 'leafC', 'rootC', 'SapwoodC', 'WoodC', &
        'NSN', 'SeedN', 'leafN', 'rootN', 'SapwoodN', 'WoodN', &
        'McrbC','fastSOM',   'SlowSOM',                        &
        'McrbN','fastSoilN', 'slowSoilN',                      &
        'mineralN', 'N_fxed','N_uptk','N_yrMin','N_P2S','N_loss', &
        'seedC','seedN','Seedling-C','Seedling-N'


   ! Parameter initialization: Initialize PFT parameters
   call initialize_PFT_data(namelistfile)
   ! Initialize vegetation tile and plant cohorts
   allocate(vegn)
   call initialize_vegn_tile(vegn,nCohorts,namelistfile)
   ! Sort and relayer cohorts
   call relayer_cohorts(vegn)
   call Zero_diagnostics(vegn)

   ! Read in forcing data
   !call read_FACEforcing(forcingData,datalines,days_data,yr_data,timestep)
   call read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
   steps_per_day = int(24.0/timestep)
   dt_fast_yr = 1.0/(365.0 * steps_per_day)
   step_seconds = 24.0*3600.0/steps_per_day ! seconds_per_year * dt_fast_yr
   write(*,*)steps_per_day,dt_fast_yr,step_seconds
   ! total years of model run
   totyears = model_run_years
   totdays  = INT(totyears/yr_data+1)*days_data
   equi_days = totdays - days_data

   ! ----- model run ---------- ! Model run starts here !!
   year0 = forcingData(1)%year
   iyears = 1
   idoy   = 0
   simu_steps = 0
   do idays =1, totdays ! 1*days_data ! days for the model run
        idoy = idoy + 1
        ! get daily mean temperature
        vegn%Tc_daily = 0.0
        tsoil         = 0.0
        do i=1,steps_per_day
             idata = MOD(simu_steps, datalines)+1
             year0 = forcingData(idata)%year  ! Current year
             vegn%Tc_daily = vegn%Tc_daily + forcingData(idata)%Tair
             tsoil         = forcingData(idata)%tsoil
             simu_steps = simu_steps + 1

             !! fast-step calls, hourly or half-hourly
             call vegn_CNW_budget_fast(vegn,forcingData(idata))
             ! diagnostics
             call hourly_diagnostics(vegn,forcingData(idata),iyears,idoy,i,idays,fno1)
        enddo ! hourly or half-hourly
        vegn%Tc_daily = vegn%Tc_daily/steps_per_day
        tsoil         = tsoil/steps_per_day
        soil_theta    = vegn%thetaS
        !write(*,*)idays,equi_days
        call daily_diagnostics(vegn,forcingData(idata),iyears,idoy,idays,fno3,fno4)
        !write(*,*)iyears,idoy
        ! daily calls
        call vegn_phenology(vegn,j)
        !call vegn_starvation(vegn)
        call vegn_growth_EW(vegn)

        !! annual calls
        idata = MOD(simu_steps+1, datalines)+1 !
        year1 = forcingData(idata)%year  ! Check if it is the last day of a year
        new_annual_cycle = ((year0 /= year1).OR. & ! new year
                (idata == steps_per_day .and. simu_steps > datalines)) ! last line
        if(new_annual_cycle)then

            idoy = 0
            !call annual_calls(vegn)
            if(update_annualLAImax) call vegn_annualLAImax_update(vegn)

            ! mortality

            call annual_diagnostics(vegn,iyears,fno2,fno5)
            call vegn_annual_starvation(vegn)
            call vegn_nat_mortality(vegn, real(seconds_per_year))

            ! Reproduction and Re-organize cohorts
            call vegn_reproduction(vegn)
            call kill_lowdensity_cohorts(vegn)
            call relayer_cohorts(vegn)
            call vegn_mergecohorts(vegn)

            ! set annual variables zero
            call Zero_diagnostics(vegn)

            ! update the years of model run
            iyears = iyears + 1
        endif
   enddo

   !deallocate(cc)
   close(91)
   close(101)
   close(102)
   close(103)
   close(104)
   deallocate(vegn%cohorts)
   deallocate(forcingData)

  contains

!========================================================================
! read in forcing data (Users need to write their own data input procedure)
subroutine read_FACEforcing(forcingData,datalines,days_data,yr_data,timestep)
  type(climate_data_type),pointer,intent(inout) :: forcingData(:)
  integer,intent(inout) :: datalines,days_data,yr_data
  real, intent(inout)   :: timestep
  !------------local var -------------------
  type(climate_data_type), pointer :: climateData(:)
  character(len=80)  commts
  integer, parameter :: niterms=9       ! MDK data for Oak Ridge input
  integer, parameter :: ilines=22*366*24 ! the maxmum records of Oak Ridge FACE, 1999~2007
  integer,dimension(ilines) :: year_data
  real,   dimension(ilines) :: doy_data,hour_data
  real input_data(niterms,ilines)
  real inputstep
  integer :: istat1,istat2,istat3
  integer :: doy,idays
  integer :: i,j,k
  integer :: m,n

  climfile=trim(filepath_in)//trim(climfile)

! open forcing data
  open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
  write(*,*)istat2
! skip 2 lines of input met data file
  read(11,'(a160)') commts
! read(11,'(a160)') commts ! MDK data only has one line comments
  m       = 0  ! to record the lines in a file
  idays   = 1  ! the total days in a data file
  yr_data = 0 ! to record years of a dataset
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
      !write(*,*)year_data(m),doy_data(m),hour_data(m)
  enddo ! end of reading the forcing file

  timestep = hour_data(2) - hour_data(1)
  write(*,*)"forcing",datalines,yr_data,timestep,dt_fast_yr
  if (timestep==1.0)then
      write(*,*)"the data freqency is hourly"
  elseif(timestep==0.5)then
      write(*,*)"the data freqency is half hourly"
  else
      write(*,*)"Please check time step!"
      stop
  endif
  close(11)    ! close forcing file
! Put the data into forcing 
  datalines = m - 1
  days_data = idays
  yr_data  = year_data(datalines-1) - year_data(1) + 1

  allocate(climateData(datalines))
  do i=1,datalines
     climateData(i)%year      = year_data(i)          ! Year
     climateData(i)%doy       = doy_data(i)           ! day of the year
     climateData(i)%hod       = hour_data(i)          ! hour of the day
     climateData(i)%PAR       = input_data(1,i)       ! umol/m2/s
     climateData(i)%radiation = input_data(2,i)       ! W/m2
     climateData(i)%Tair      = input_data(3,i) + 273.16  ! air temperature, K
     climateData(i)%Tsoil     = input_data(4,i) + 273.16  ! soil temperature, K
     climateData(i)%RH        = input_data(5,i) * 0.01    ! relative humidity (0.xx)
     climateData(i)%rain      = input_data(6,i)/(timestep * 3600)! ! kgH2O m-2 s-1
     climateData(i)%windU     = input_data(7,i)        ! wind velocity (m s-1)
     climateData(i)%P_air     = input_data(8,i)        ! pa
     climateData(i)%CO2       = input_data(9,i) * 1.0e-6       ! mol/mol
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  write(*,*)"forcing", datalines,days_data,yr_data
end subroutine read_FACEforcing

!=============================================================
! for reading in NACP site synthesis forcing
subroutine read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
  type(climate_data_type),pointer,intent(inout) :: forcingData(:)
  integer,intent(inout) :: datalines,days_data,yr_data
  real, intent(inout)   :: timestep
  !------------local var -------------------
  type(climate_data_type), pointer :: climateData(:)
  character(len=80)  commts
  integer, parameter :: niterms=15       ! NACP site forcing
  integer, parameter :: ilines=22*366*48 ! the maxmum records
  integer,dimension(ilines) :: year_data
  real,   dimension(ilines) :: doy_data,hour_data
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
  yr_data = 0 ! to record years of a dataset
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
      !write(*,*)year_data(m),doy_data(m),hour_data(m)
      ! discard one line
      !read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      !                        (input_data(n,m),n=1,niterms)
  enddo ! end of reading the forcing file

  timestep = hour_data(2) - hour_data(1)
  write(*,*)"forcing",datalines,yr_data,timestep,dt_fast_yr
  if (timestep==1.0)then
      write(*,*)"the data freqency is hourly"
  elseif(timestep==0.5)then
      write(*,*)"the data freqency is half hourly"
  else
      write(*,*)"Please check time step!"
      stop
  endif
  close(11)    ! close forcing file
! Put the data into forcing 
  datalines = m - 1
  days_data = idays
  yr_data  = year_data(datalines-1) - year_data(1) + 1

  allocate(climateData(datalines))
  do i=1,datalines
     climateData(i)%year      = year_data(i)          ! Year
     climateData(i)%doy       = doy_data(i)           ! day of the year
     climateData(i)%hod       = hour_data(i)          ! hour of the day
     climateData(i)%PAR       = input_data(11,i)*2.0  ! umol/m2/s
     climateData(i)%radiation = input_data(11,i)      ! W/m2
     climateData(i)%Tair      = input_data(1,i)       ! air temperature, K
     climateData(i)%Tsoil     = input_data(1,i)       ! soil temperature, K
     climateData(i)%rain      = input_data(7,i)       ! kgH2O m-2 s-1
     climateData(i)%windU     = input_data(5,i)        ! wind velocity (m s-1)
     climateData(i)%P_air     = input_data(9,i)        ! pa
     climateData(i)%RH        = input_data(3,i)/mol_h2o*mol_air* & ! relative humidity (0.xx)
                                climateData(i)%P_air/esat(climateData(i)%Tair-273.16)
     climateData(i)%CO2       = input_data(15,i) * 1.0e-6       ! mol/mol
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  write(*,*)"forcing", datalines,days_data,yr_data
  
end subroutine read_NACPforcing


!=====================================================
end program BiomeESS



