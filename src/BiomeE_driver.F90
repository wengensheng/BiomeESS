#define Do_demography

!---------------
module BiomeE_mod
 use datatypes
 use io_mod
 use esdvm

 implicit none
 private

 ! ------ public subroutines ---------
 public :: BiomeE_Initialization, BiomeE_run, BiomeE_end

 ! Main vegn unit
 type(vegn_tile_type),public,pointer :: vegn
 character(len=50) :: runID, fnamelist, fpath_out
 integer :: nml_unit, iostat,io ! for reading the namelist file
 integer :: fno1, fno2, fno3, fno4, fno5, fno6 ! output files
 integer,parameter :: N_initCC = 1

 contains

!================== BiomeE Driver ===========================================
!----------------------------------------------------------------------------
subroutine BiomeE_Initialization()
  ! Weng 08/08/2022, for model initialization
  implicit none
  integer :: timeArray(3), rand_seed
  real    :: r_rand

  ! ---------------------- Setup run files ----------------------
  runID = 'BCI_hydro' ! 'OR_phiRL' ! 'Konza2' !  'OR_GAPLUE' !  'FACE_hydro' !
  fnamelist = 'parameters_'//trim(runID)//'.nml'

  ! --------- Read forcing data ----------------------
  nml_unit = 901
  open(nml_unit, file=fnamelist, form='formatted', action='read', status='old')
  read (nml_unit, nml=initial_state_nml, iostat=io)
  close (nml_unit)

  call read_FACEforcing(forcingData,datalines,days_data,yr_data,step_hour)
  !call read_NACPforcing(forcingData,datalines,days_data,yr_data,step_hour)
  steps_per_day = int(24.0/step_hour)
  dt_fast_yr   = step_hour/(365.0 * 24.0)
  step_seconds = step_hour*3600.0
  write(*,*)'steps/day,dt_fast,s/step',steps_per_day,dt_fast_yr,step_seconds
  ! total years of model run
  totyears  = model_run_years
  totdays   = INT(totyears/yr_data+1)*days_data
  equi_days = Max(0, totdays - days_data)

  ! Set up rainfall scenario for phiRL test runs
  forcingData%rain = forcingData%rain * Sc_prcp

  ! --------- Model initialization: ---------------
  ! create output files
  fpath_out='output/'
  call set_up_output_files(runID,fpath_out,fno1,fno2,fno3,fno4,fno5,fno6)

  ! Setup random numbers ! call random_seed()
  call itime(timeArray)     ! Get the current time
  rand_seed = timeArray(1)+timeArray(2)+timeArray(3)
  r_rand    = rand(rand_seed)

  ! ------Initialize soil and PFT parameters
  call initialize_soilpars(fnamelist)
  call initialize_PFT_data(fnamelist)

  ! Initialize vegetation tile and plant cohorts
  allocate(vegn)
  call initialize_vegn_tile(vegn,N_initCC,fnamelist)
  ! Sort and relayer cohorts
  call relayer_cohorts(vegn)
  call Zero_diagnostics(vegn)

end subroutine BiomeE_Initialization

!----------------------------------------------------------------------------
subroutine BiomeE_run()
  ! Weng 08/08/2022, for model run
  implicit none
  logical :: new_annual_cycle
  real    :: tsoil, soil_theta
  integer :: year0, year1, idays, simu_steps,idata
  integer :: idoy, iyears, i, j, k

  !----------------------
  year0 = forcingData(1)%year
  iyears = 1
  idoy   = 0
  simu_steps = 0
  vegn%Tc_pheno = forcingData(1)%Tair
  do idays =1, totdays ! 1*days_data ! days for the model run
    idoy = idoy + 1
    ! get daily mean temperature
    vegn%Tc_daily = 0.0
    tsoil         = 0.0
    do i=1,steps_per_day ! fast-step calls, hourly or half-hourly
      idata = MOD(simu_steps, datalines)+1
      year0 = forcingData(idata)%year  ! Current year
      vegn%Tc_daily = vegn%Tc_daily + forcingData(idata)%Tair
      tsoil         = forcingData(idata)%tsoil
      simu_steps = simu_steps + 1
      call vegn_CNW_budget_fast(vegn,forcingData(idata))
      call hourly_diagnostics(vegn,forcingData(idata),iyears,idoy,i,idays,fno1,fno2)
    enddo ! steps_per_day
    vegn%Tc_daily = vegn%Tc_daily/steps_per_day
    tsoil         = tsoil/steps_per_day
    soil_theta    = vegn%thetaS

    call daily_diagnostics(vegn,iyears,idoy,idays,fno3,fno4)
    ! daily calls
    call vegn_phenology(vegn,j)
    !call vegn_daily_starvation(vegn)
    call vegn_growth(vegn)

    !! annual calls
    idata = MOD(simu_steps+1, datalines)+1 !
    year1 = forcingData(idata)%year  ! Check if it is the last day of a year
    new_annual_cycle = ((year0 /= year1).OR. & ! new year
            (idata == steps_per_day .and. simu_steps > datalines)) ! last line
    if(new_annual_cycle)then ! Annual calculations
      idoy = 0
      ! Update plant hydraulic states, for the last year
      call vegn_hydraulic_states(vegn,real(seconds_per_year))
      call annual_diagnostics(vegn,iyears,fno5,fno6)

      ! Case studies
      ! N is losing after changing the soil pool structure. Hack !!!!!
      ! if(do_closedN_run) call Recover_N_balance(vegn)
      if(do_fire) call vegn_fire (vegn, real(seconds_per_year))
      if(do_migration) call vegn_migration(vegn) ! for competition
      if(update_annualLAImax) call vegn_annualLAImax_update(vegn)

#ifdef Do_demography
      ! For the incoming year
      call vegn_annual_starvation(vegn) ! turn it off for grass run
      call vegn_nat_mortality(vegn, real(seconds_per_year))
      ! Reproduction and Reorganize cohorts
      call vegn_reproduction(vegn)
#endif

      ! Cohort management
      call kill_lowdensity_cohorts(vegn)
      !call vegn_gap_fraction_update(vegn) !for CROWN_GAP_FILLING
      call relayer_cohorts(vegn)
      call vegn_mergecohorts(vegn)

      ! zero annual reporting variables
      call Zero_diagnostics(vegn)

      ! update the years of model run
      iyears = iyears + 1
    endif
  enddo
end subroutine BiomeE_run

!----------------------------------------------------------------------------
subroutine BiomeE_end
  !------------ Close output files and release memory
  close(fno1); close(fno2); close(fno3)
  close(fno4); close(fno5); close(fno6)

  deallocate(vegn%cohorts)
  deallocate(forcingData)
end subroutine BiomeE_end

!----------------------------------------------------------------------------

end module BiomeE_mod
