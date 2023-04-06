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

!#define DemographyOFF

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
 ! Main vegn unit
 type(land_grid_type), pointer :: land
 type(vegn_tile_type), pointer :: vegn
 ! Input forcing data
 type(climate_data_type), pointer :: forcingData(:)
 ! output files
 integer :: fno1, fno2, fno3, fno4, fno5, fno6

 contains

!================== BiomeE Driver ===========================================
 subroutine BiomeE_main() ! Weng 03/20/2023, main BiomeE module
  implicit none

  call BiomeE_initialization()
  call BiomeE_run()
  call BiomeE_end()
end subroutine BiomeE_main

!----------------------------------------------------------------------------
subroutine BiomeE_initialization()
  ! Weng 08/08/2022, for model initialization
  implicit none
  type(vegn_tile_type), pointer :: pveg => NULL()
  integer :: timeArray(3),i
  !=============== Namelist file (must be hardwired) ===============
  character(len=80)  :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'

  ! ---------- Time stamp -------------
  call itime(timeArray)     ! Get current time
  write(*,'(a6,3(I2,":"))')'Time: ', timeArray

  ! ---------------------- Read the namelist file -----------------
  call read_namelist(fnml)

  ! --------- Read forcing data ----------------------
  call read_FACEforcing(forcingData,datalines,days_data,yr_data,step_hour)
  !call read_NACPforcing(forcingData,datalines,days_data,yr_data,step_hour)
  !call read_CRUforcing(forcingData,datalines,days_data,yr_data,step_hour)
  steps_per_day = int(24.0/step_hour)
  dt_fast_yr    = step_hour/(365.0 * 24.0)
  step_seconds  = step_hour*3600.0
  write(*,*)'steps/day,dt_fast,s/step',steps_per_day,dt_fast_yr,step_seconds
  ! total years of model run
  totyears  = model_run_years
  totdays   = INT(totyears/yr_data+1)*days_data
  equi_days = Max(0, totdays - days_data)

  ! Set up scenarios for rainfall and CO2 concentration
  forcingData%rain = forcingData%rain * Sc_prcp
  forcingData%CO2  = CO2_c * 1.0e-6
  !stop

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
  write(*,*)'total tiles:', N_VegTile

  ! --------- Setup output files ---------------
  call set_up_output_files(fno1,fno2,fno3,fno4,fno5,fno6)

  ! ------ Start a new random number series ------
  call RANDOM_SEED()

end subroutine BiomeE_initialization

!----------------------------------------------------------------------------
subroutine BiomeE_run()
  ! Weng 08/08/2022, for model run
  implicit none
  integer :: i, idays, idata, idoy
  integer :: n_steps, n_yr, year0, year1
  real    :: r_d
  logical :: new_annual_cycle

  !----------------------
  vegn%Tc_pheno = forcingData(1)%Tair
  year0   = forcingData(1)%year
  n_yr    = 1
  idoy    = 0
  n_steps = 0
  do idays =1, totdays ! 1*days_data ! days for the model run
    idoy = idoy + 1
    land%Tc_daily = 0.0 ! Zero daily mean temperature

    ! Fast-step calls, hourly or half-hourly
    do i=1,steps_per_day
      n_steps = n_steps + 1
      idata = MOD(n_steps-1, datalines) + 1
      year0 = forcingData(idata)%year  ! Current step year
      land%Tc_daily = land%Tc_daily + forcingData(idata)%Tair
      vegn => land%firstVegn
      do while(ASSOCIATED(vegn))
        call vegn_CNW_budget_fast(vegn,forcingData(idata))
        call hourly_diagnostics(vegn,forcingData(idata), &
                                n_yr,idoy,i,idays,fno1,fno2)
        vegn => vegn%next
      enddo
    enddo ! steps_per_day
    land%Tc_daily = land%Tc_daily/steps_per_day

    ! Daily run
    vegn => land%firstVegn
    do while(ASSOCIATED(vegn))
      ! Update vegn age
      call vegn_age(vegn,dt_daily_yr)
      call daily_diagnostics(vegn,n_yr,idoy,idays,fno3,fno4)
      ! Daily calls
      vegn%Tc_daily = land%Tc_daily
      call vegn_phenology(vegn)
      call vegn_growth(vegn)
      !call vegn_daily_starvation(vegn)
      vegn => vegn%next
    enddo

    !! Check if the next step is a new year
    idata = MOD(n_steps, datalines) + 1
    year1 = forcingData(idata)%year  ! Nex step year
    new_annual_cycle = ((year0 /= year1) .OR. & ! new year
                (MOD(n_steps, datalines) == 0)) ! last line

    ! Annual calculations
    if(new_annual_cycle)then
      idoy = 0
      vegn => land%firstVegn
      do while(ASSOCIATED(vegn))
        ! Update plant hydraulic states, for the last year
        call vegn_hydraulic_states(vegn,real(seconds_per_year))
        call annual_diagnostics(vegn,n_yr,fno5,fno6)

        ! Case studies
        ! N is losing after changing the soil pool structure. Hack !!!!!
        ! if(do_closedN_run) call Recover_N_balance(vegn)
        ! if(do_fire) call vegn_fire (vegn, real(seconds_per_year))
        if(do_migration) call vegn_migration(vegn) ! for competition
        ! if(update_annualLAImax) call vegn_annualLAImax_update(vegn)
#ifndef DemographyOFF
        ! For the incoming year
        call vegn_annual_starvation(vegn) ! turn it off for grass run
        call vegn_nat_mortality(vegn, real(seconds_per_year))
        call vegn_reproduction(vegn)
#endif

        ! Cohort management
        call kill_lowdensity_cohorts(vegn)
        call kill_old_grass(vegn)
        !call vegn_gap_fraction_update(vegn) !for CROWN_GAP_FILLING
        call relayer_cohorts(vegn)
        call vegn_mergecohorts(vegn)

        ! zero annual reporting variables
        call Zero_diagnostics(vegn)
        !! Reset vegetation to initial conditions
        CALL RANDOM_NUMBER(r_d)
        if((n_yr==yr_ResetVeg).or.(n_yr>yr_ResetVeg .and. r_d<envi_fire_prb)) &
          call reset_vegn_initial(vegn)
        vegn => vegn%next
      enddo

      ! update the years of model run
      n_yr = n_yr + 1
    endif
  enddo
end subroutine BiomeE_run

!----------------------------------------------------------------------------
subroutine BiomeE_end
  type(vegn_tile_type), pointer :: pveg => NULL()
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
end subroutine BiomeE_end

!----------------------------------------------------------------------------
end module BiomeE_mod
