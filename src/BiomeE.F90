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
  use model_utils
  use io_mod
  use esdvm
  use restart_mod

  implicit none
  private

  ! ------ public subroutines ---------
  public :: BiomeE_main

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
    type(vegn_tile_type), pointer :: vegn => NULL()
    type(vegn_tile_type), pointer :: pveg => NULL()
    integer :: i

    ! Setup total days of model run (forcing data have been read in)
    totdays  = INT(model_run_years/data_yrs+1)*data_days
    totyears = INT(model_run_years/data_yrs+1)*data_yrs
    if(output_days > 0)then
      skipped_days = totdays - output_days
    else
      skipped_days = totdays - data_days
    endif

    ! ------ Land grid, vegetation tiles, and plant cohorts ------
    allocate(land)
    land%nTiles = 0

    if (do_restart_read) then
      ! ---- Restore full model state from restart file ----
      call read_restart(land)
      ! Rebuild canopy layers and zero flux accumulators
      vegn => land%firstVegn
      do while (associated(vegn))
        call vegn_RelayerCohorts(vegn)
        call vegn_sum_tile(vegn)
        call Zero_diagnostics(vegn)
        vegn => vegn%next
      enddo
    else
      ! ---- Original (spin-up) initialization path ----

      ! Setup initial PFTs and cohorts information (not cohorts per se)
#ifdef DO_Climate_VEG
      ! Update init_cohort_* arrays, 09/09/2025
      call Climate_envelope_vars(forcingData,steps_per_day)
      call Set_PFTs_from_Data()
#endif

      if(init_cohort_Indiv(1)<0.0) &
      call Assign_Std_Cohorts (init_cohort_sps,init_cohort_N)

      do i =1, N_VegTile
        allocate(vegn)
        call initialize_vegn_tile(vegn)
        call vegn_RelayerCohorts(vegn)
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

    endif ! do_restart_read

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
    real    :: r_d
    logical :: new_annual_cycle
#ifdef DroughtMIP
    logical :: BaseLineClimate = .True.
#endif

#ifdef HistCO2
    integer :: spin_yrs,hist_yrs,i_hist ! for FACE MDS III

    ! Total model run years and spin-up years
    hist_yrs = Max(CO2_end_yr - CO2_start_yr + 1, 1)
    spin_yrs = totyears - hist_yrs - post_yrs
    i_hist   = Max(CO2_start_yr - 1700, 0) + 1
#endif

    !----------------------
    n_yr    = 1
    idoy    = 0
    MonthDays = MonthDOY
    n_steps = 0
    do idays = 1, totdays ! Total days for model run
      idoy = idoy + 1
      ! Leap year or not (CRU data has 365 days/yr; use year number for calendar check)
      if(idoy == 1)then
        jdata = MOD(n_steps, datalines) + 1
        year0 = forcingData(jdata)%year
        if(MOD(year0,4)==0 .and. (MOD(year0,100)/=0 .or. MOD(year0,400)==0)) then
          MonthDays(2:12) = MonthDOY(2:12)+1 ! leap year
        else
          MonthDays = MonthDOY ! non-leap year
        endif
      endif

      land%Tc_daily = 0.0 ! Zero daily mean temperature
      ! Fast-step update (hourly or half-hourly)
      do i=1,steps_per_day
        n_steps = n_steps + 1
        idata = MOD(n_steps-1, datalines) + 1
        climateData = forcingData(idata)
        ! Set up scenarios for rainfall and CO2 concentration
        climateData%rain = forcingData(idata)%rain * Sc_prcp
        climateData%Tair = forcingData(idata)%Tair + Sc_dT
#ifdef HistCO2
        climateData%CO2 = CO2_Hist(i_hist)
        if(n_yr > spin_yrs .and. CO2Tag == 'eCO2') &
          climateData%CO2 = climateData%CO2 + dCO2
#else
        if(fixedCO2) climateData%CO2 = CO2_c ! ppm
        if(CO2Tag == 'eCO2') climateData%CO2 = forcingData(idata)%eCO2
#endif
        land%Tc_daily = land%Tc_daily + climateData%Tair - 273.16

        vegn => land%firstVegn
        do while(ASSOCIATED(vegn))
          vegn%CO2_c = climateData%CO2 ! * 1.0e6
          call vegn_CNW_budget_fast(vegn,climateData)
          call hourly_diagnostics(vegn,climateData,n_yr,idoy,i,idays)
          vegn => vegn%next
        enddo
      enddo ! steps_per_day
      land%Tc_daily = land%Tc_daily/steps_per_day

      ! Daily update
      vegn => land%firstVegn
      do while(ASSOCIATED(vegn))
        vegn%Tc_daily = land%Tc_daily
        call vegn_daily_update(vegn,dt_daily_yr)
        call daily_diagnostics(vegn,n_yr,idoy,idays,MonthDays)
        vegn => vegn%next
      enddo

      ! Annual update
      ! Southern hemisphere has been shifted up for 182 days with "ShiftSHdata"
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
          if(do_fire) call vegn_fire(vegn,real(seconds_per_year))

#ifdef SingleTreeTest
          call vegn_SingleCohort_annual_update(vegn)
          call annual_diagnostics(vegn,n_yr)
#else
          call annual_diagnostics(vegn,n_yr)
          call vegn_demographics(vegn,real(seconds_per_year))
#endif

          ! Case studies
          if(do_RecoverSP .and. MOD(n_yr, FreqY0)==0) &
          call vegn_species_recovery(vegn) ! for competition
          ! if(update_annualLAImax) call vegn_annualLAImax_update(vegn)

          ! --------- Cohort management ---------
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
          call vegn_RelayerCohorts(vegn)
          call vegn_mergecohorts(vegn)
          call vegn_Remove_empty_cc(vegn)
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

#ifdef HistCO2
        ! CO2 concentration for this year
        write(*,*)'i_hist, CO2Yrs, spin_yrs',i_hist, CO2Yrs, spin_yrs
        write(*,*)'Used CO2 concentration:',climateData%CO2, CO2_Hist(i_hist)
        ! Next Year's i_hist
        if(n_yr > spin_yrs .and. n_yr <= spin_yrs+hist_yrs) i_hist = Min(i_hist + 1, CO2Yrs)
#endif

#ifdef DroughtMIP
        if(n_yr == yr_Baseline + 1 .and. BaseLineClimate)then
          call setup_forcingdata(Scefile)
          n_steps = 0
          BaseLineClimate = .False.
        endif
#endif

      endif
    enddo
  end subroutine BiomeE_run

!----------------------------------------------------------------------------
subroutine BiomeE_end
  type(vegn_tile_type), pointer :: vegn => null()
  type(vegn_tile_type), pointer :: pveg => null()
  logical :: is_open

  call flush_annual_diagnostics_buffers()

  !------------ Write restart checkpoint (before closing output files)
  if (do_restart_write) call write_restart(land)

  !------------ Close output files (only if opened)
  inquire(unit=fno1, opened=is_open); if (is_open) close(fno1)
  inquire(unit=fno2, opened=is_open); if (is_open) close(fno2)
  inquire(unit=fno3, opened=is_open); if (is_open) close(fno3)
  inquire(unit=fno4, opened=is_open); if (is_open) close(fno4)
  inquire(unit=fno5, opened=is_open); if (is_open) close(fno5)
  inquire(unit=fno6, opened=is_open); if (is_open) close(fno6)

  !------------ Release vegetation tiles/cohorts
  if (associated(land)) then
    vegn => land%firstVegn
    do while (associated(vegn))
      pveg => vegn%next

      if (associated(vegn%cohorts)) then
        deallocate(vegn%cohorts)
        nullify(vegn%cohorts)
      endif

      if (associated(vegn%initialCC)) then
        deallocate(vegn%initialCC)
        nullify(vegn%initialCC)
      endif

      deallocate(vegn)
      vegn => pveg
    enddo
    deallocate(land)
    nullify(land)
    nullify(vegn)
  endif

  !------------ Release forcing data
  if (associated(forcingData)) then
    deallocate(forcingData)
    nullify(forcingData)
  endif
end subroutine BiomeE_end

!----------------------------------------------------------------------------
end module BiomeE_mod
