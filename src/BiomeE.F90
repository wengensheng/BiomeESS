!========================================================================
!========= Biome Ecological strategy simulator (BiomeE) =================
!=======================   Main program   ===============================
!========================== 12-30-2019 ==================================
!========================================================================
!
! This work was financially supported by US Forest Service and Princeton
! Environment Institute. The technical details of this model can be found
! in:
!
! Weng, E., Dybzinski, R., Farrior, C. E., and Pacala, S. W.: Competition
! alters predicted forest carbon cycle responses to nitrogen availability
! and elevated CO2: simulations using an explicitly competitive, game-
! theoretic vegetation demographic model, Biogeosciences, 16, 4577–4599,
! https://doi.org/10.5194/bg-16-4577-2019, 2019.
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
! Contact Ensheng Weng (wengensheng@gmail.com) for questions.
!
! (Lase updated 12/30/2017, 07/24/2022, 01/03/2023)
!========================================================================

!------------------------------------------------------------------------
! Processes included in this simulator are:
!     Photosynthesis, transpiration, plant respiration, Plant hydraulics
!     Plant growth: Allometry and allocation; Phenology
!     Demography: Reproduction, Mortality
!     Population dynamics
!     Soil respiration, Soil C-N dynamics
!     Soil water dynamics: soil surface evaporation, infiltration, runoff
!
!----------------------------- END ----------------------------------

!---------------
module BiomeE_mod
 use datatypes
 use io_mod
 use esdvm

 implicit none
 private

 ! ------ public subroutines ---------
 public :: BiomeE_main
 public :: BiomeE_Initialization, BiomeE_run, BiomeE_end

 !-------------Vars for the model -------------
 type(land_grid_type), pointer :: land ! Land grid
 ! Input forcing data
 type(climate_data_type), pointer :: forcingData(:)
 ! output files
 integer :: fno1, fno2, fno3, fno4, fno5, fno6

 contains

!================== BiomeE Driver ===========================================
 subroutine BiomeE_main() ! Weng 03/20/2023, main BiomeE module
  implicit none
  !=============== Namelist file (must be hardwired) ===============
  character(len=80)  :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'

  call read_namelist(fnml)
  call setup_forcingdata()
  call BiomeE_initialization()
  call BiomeE_run()
  call BiomeE_end()
end subroutine BiomeE_main

! --------- Setup forcing data and step lenght ----------------------
subroutine setup_forcingdata()

#ifdef DroughtPaleo
  call set_PaleoForcing(climfile,PaleoPfile,PaleoTfile,iDraw, &
        forcingData,datalines,days_data,yr_data,step_hour)
#else
  call read_FACEforcing(climfile,forcingData,datalines,days_data,yr_data,step_hour)
  !call read_NACPforcing(forcingData,datalines,days_data,yr_data,step_hour)
  !call read_CRUforcing(forcingData,datalines,days_data,yr_data,step_hour)
#endif

  ! ------ Setup steps for model run ------
  steps_per_day = int(24.0/step_hour)
  dt_fast_yr    = step_hour/(365.0 * 24.0)
  step_seconds  = step_hour*3600.0
  write(*,*)'steps/day,dt_fast,s/step',steps_per_day,dt_fast_yr,step_seconds

  stop ! just testing data setup

end subroutine

!----------------------------------------------------------------------------
subroutine BiomeE_initialization()
  ! Weng 08/08/2022, for model initialization
  implicit none
  type(vegn_tile_type), pointer :: vegn => NULL()
  type(vegn_tile_type), pointer :: pveg => NULL()
  integer :: timeArray(3),i

  ! ---------- Time stamp -------------
  call itime(timeArray)     ! Get current time
  write(*,'(a12,3(I2,":"))')'Start time: ', timeArray

  ! Hack for closedN setting
  if(do_closedN_run) then
    K_nitrogen = 0.0 ! rate of a year, 2.5
    rho_SON    = 0.0 ! organic nitrogen release rate
    etaN       = 0.0 ! Nitrogen loss rate with runoff
    fDON       = 0.0 ! DON fraction
    N_input    = 0.0 ! N input, kg N m-2 yr-1
  endif

  ! Setup total days of model run
  totdays   = INT(model_run_years/yr_data+1)*days_data
  equi_days = Max(0, totdays - days_data)

  ! ------ Soil and PFT parameters ------
  call initialize_soilpars()
  call initialize_PFT_data()

  ! ------ Land grid, vegetation tiles, and plant cohorts ------
  allocate(land)
  land%nTiles = 0
  do i =1, N_VegTile
    allocate(vegn)
    call initialize_vegn_tile(vegn)
    ! Sort and relayer cohorts
    call relayer_cohorts(vegn)
    call Zero_diagnostics(vegn)
    vegn%Tc_pheno = forcingData(1)%Tair
    vegn%tileID = i
    land%nTiles = land%nTiles + 1
    if(i==1)then
      land%firstVegn => vegn
      pveg => vegn
    else
      pveg%next => vegn
      vegn%prev => pveg
      pveg      => vegn
    endif
    vegn => NULL()
  enddo
  vegn => land%firstVegn
  pveg => NULL()

  ! --------- Setup output files ---------------
  call set_up_output_files(fno1,fno2,fno3,fno4,fno5,fno6)

  ! ------ Start a new random number series ------
  call RANDOM_SEED()

end subroutine BiomeE_initialization

!----------------------------------------------------------------------------
subroutine BiomeE_run()
  ! Weng 08/08/2022, for model run
  implicit none
  type(vegn_tile_type), pointer :: vegn => NULL()
  type(climate_data_type) :: climateData
  integer :: i, k, idays, idata, jdata, idoy
  integer :: n_steps, n_yr, year0, year1
  integer :: MonthDays(0:12)
  integer :: iCO2_hist, CO2_start_yr, CO2_end_yr, CO2_yrs, skip_yrs
  integer :: tot_yrs,spin_yrs,hist_yrs ! for FACE MDS III
  real    :: r_d
  logical :: new_annual_cycle
  logical :: BaseLineClimate = .True.

#ifdef FACE_run
  ! History CO2 concentration, data from 1700 to 2024 (325 years)
  CO2_start_yr = 1850 ! Minimum 1700
  CO2_end_yr   = 1997
  CO2_yrs      = Max(CO2_end_yr - CO2_start_yr + 1, 1)
  skip_yrs     = Max(CO2_start_yr - 1700, 0)
  iCO2_hist= skip_yrs + 1
  !Total years, CO2-history years, and experiment years
  tot_yrs  = INT(model_run_years/yr_data+1)*yr_data
  hist_yrs = CO2_yrs
  spin_yrs = tot_yrs ! in case there is unrecognizable CO2Tag
  if(CO2Tag == 'Init')then
    spin_yrs = tot_yrs
  elseif(CO2Tag == 'Hist')then
    spin_yrs = tot_yrs - hist_yrs
  elseif(CO2Tag == 'aCO2' .OR. CO2Tag == 'eCO2')then
    spin_yrs = INT(750/yr_data)*yr_data - hist_yrs
  endif
#endif

  !----------------------
  n_yr    = 1
  idoy    = 0
  MonthDays = MonthDOY
  n_steps = StartLine - 1 ! steps skipped acc. the starting line, for UFL only
  do idays =1, totdays - (StartLine - 1)/steps_per_day ! 1*days_data ! days for the model run
    idoy = idoy + 1
    ! Leap year or not
    if(idoy == 1)then
      jdata = MIN(datalines,MOD(n_steps, datalines) + 59*steps_per_day+1)
      if(forcingData(jdata)%doy>28)then
        MonthDays(2:12) = MonthDOY(2:12)+1!leap year
      else
        MonthDays = MonthDOY!Non-leap year
      endif
    endif

    land%Tc_daily = 0.0 ! Zero daily mean temperature
    ! Fast-step update (hourly or half-hourly)
    do i=1,steps_per_day
      n_steps = n_steps + 1
      idata = MOD(n_steps-1, datalines) + 1
      climateData = forcingData(idata)
      ! Set up scenarios for rainfall and CO2 concentration
      climateData%rain = forcingData(idata)%rain !* Sc_prcp
      climateData%CO2  = CO2_c ! ppm
#ifdef FACE_run
      if(n_yr <= spin_yrs)then
        climateData%CO2  = CO2_c ! ppm
      elseif(n_yr > spin_yrs .and. n_yr <= spin_yrs+hist_yrs)then ! Treatment CO2 (Hist, aCO2, eCO2)
        climateData%CO2  = CO2_Hist(iCO2_hist)
      elseif(n_yr > spin_yrs+hist_yrs)then
        if(CO2Tag == 'aCO2')then
          climateData%CO2 = forcingData(idata)%CO2
        elseif(CO2Tag == 'eCO2')then
          climateData%CO2 = forcingData(idata)%eCO2
        endif
      endif
#endif
#ifdef UFL_Monoculture
      if(n_yr > 500)then
        climateData%rain = 0.5 * forcingData(idata)%rain * Sc_prcp
      endif
#endif
      land%Tc_daily = land%Tc_daily + climateData%Tair - 273.16

      vegn => land%firstVegn
      do while(ASSOCIATED(vegn))
        if(i == INT(steps_per_day/2) .and. idoy == 181)vegn%CO2_c = climateData%CO2 ! * 1.0e6
        call vegn_CNW_budget_fast(vegn,climateData)
        call hourly_diagnostics(vegn,climateData, &
                                n_yr,idoy,i,idays,fno1,fno2)
        vegn => vegn%next
      enddo
    enddo ! steps_per_day
    land%Tc_daily = land%Tc_daily/steps_per_day

    ! Daily update
    vegn => land%firstVegn
    do while(ASSOCIATED(vegn))
      vegn%Tc_daily = land%Tc_daily
      call vegn_daily_update(vegn,dt_daily_yr)
      call daily_diagnostics(vegn,n_yr,idoy,idays,MonthDays,fno3,fno4)
      vegn => vegn%next
    enddo

    ! Annual update
    ! Check if the next step is a new year
    year0 = forcingData(idata)%year  ! Current step year
    idata = MOD(n_steps, datalines) + 1 ! Next step idata
    year1 = forcingData(idata)%year  ! Nex step year
    new_annual_cycle = ((year0 /= year1) .OR. (MOD(n_steps,datalines)==0))
    if(new_annual_cycle)then
      idoy = 0
      vegn => land%firstVegn
      do while(ASSOCIATED(vegn))
        ! Update plant hydraulic states, for the last year
        call vegn_hydraulic_states(vegn,real(seconds_per_year))
        ! Fire disturbance
        if(do_fire) call vegn_fire (vegn, real(seconds_per_year))
#ifdef SingleTreeTest
        call vegn_SingleCohort_annual_update(vegn)
        call annual_diagnostics(vegn,n_yr,fno5,fno6)
#else
        call annual_diagnostics(vegn,n_yr,fno5,fno6)
        call vegn_demographics(vegn,real(seconds_per_year))
#endif

        ! Case studies
        if(do_migration) call vegn_migration(vegn) ! for competition
        ! if(update_annualLAImax) call vegn_annualLAImax_update(vegn)

        ! Cohort management
        call kill_lowdensity_cohorts(vegn)
        ! calculate the number of cohorts with indivs>mindensity
         k = 0
         do i = 1, vegn%n_cohorts
            if (vegn%cohorts(i)%nindivs > 0.5*min_nindivs) k=k+1
         enddo
         if(k==0)then
           write(*,*)"zero cohorts, reset!"
           call reset_vegn_initial(vegn)
         endif
         call kill_old_grass(vegn)
         !call vegn_gap_fraction_update(vegn) !for CROWN_GAP_FILLING
         call relayer_cohorts(vegn)
         call vegn_mergecohorts(vegn)
        ! Summarize tile and zero annual reporting variables
        call vegn_sum_tile(vegn)
        call Zero_diagnostics(vegn)

#ifdef DBEN_runs
        !! Reset vegetation to initial conditions, for DBEN
        CALL RANDOM_NUMBER(r_d)
        if((n_yr==yr_ResetVeg).or.(n_yr>yr_ResetVeg .and. r_d<envi_fire_prb)) &
          call reset_vegn_initial(vegn)
#endif
        vegn => vegn%next
      enddo

      ! update the years of model run
      n_yr = n_yr + 1

#ifdef FACE_run
      if(n_yr > spin_yrs .and. n_yr <= spin_yrs+hist_yrs)then
        iCO2_hist = Min(iCO2_hist + 1, CO2_end_yr)
      endif
#endif

#ifdef DroughtMIP
      if(n_yr == yr_Baseline + 1 .and. BaseLineClimate)then
        call read_forcingdata(Scefile)
        n_steps = 0
        BaseLineClimate = .False.
      endif
#endif

#ifdef UFL_Monoculture
      if(n_yr > 505)exit
#endif
    endif
  enddo
end subroutine BiomeE_run

!----------------------------------------------------------------------------
subroutine BiomeE_end
  type(vegn_tile_type), pointer :: vegn => NULL()
  type(vegn_tile_type), pointer :: pveg => NULL()
  integer :: timeArray(3)

  !------------ Close output files and release memory
  close(fno1); close(fno2); close(fno3)
  close(fno4); close(fno5); close(fno6)

  vegn => land%firstVegn
  do while(ASSOCIATED(vegn))
    pveg => vegn%next
    deallocate(vegn%cohorts)
    deallocate(vegn)
    vegn => pveg
  enddo
  deallocate(land)
  deallocate(forcingData)

  ! ---------- Time stamp -------------
  call itime(timeArray)     ! Get current time
  write(*,'(a12,3(I2,":"))')'End time: ', timeArray
end subroutine BiomeE_end

!----------------------------------------------------------------------------
end module BiomeE_mod
