! BiomeE global version 0.1
! Ensheng Weng, 07/18/2025

program BiomeE_global_driver
  !use omp_lib
  use datatypes
  use netcdf_io
  use io_mod, only: read_namelist,setup_forcingdata, &
                    setup_output_files, zip_output_files
  use BiomeE_mod, only: BiomeE_main

  implicit none

  !-----------
  !type(climate_data_type), pointer :: forcing1(:),forcing2(:),forcing3(:)
  character(len=80)  :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'
  real :: start_time, end_time, last_time, elapsed_time
  real, pointer :: GridData(:,:)    ! N_yr*Ntime, N_vars
  integer :: N_yrs, totL
  integer :: m, n
  integer :: timeArray(3)

  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  last_time = start_time
  call itime(timeArray)     ! Get current time
  write(*,'(a12,3(I2,":"))')'Start time: ', timeArray

#ifdef GlobalRun
  ! Read in namelist file and setup output files
  call read_namelist(fnml) ! Namelist file (must be hardwired)
  call read_global_setting(fnml)
  ! ------ Setup steps for model run ------
  N_yrs = yr_end - yr_start + 1
  step_hour = 1.0 ! hourly time step
  steps_per_hour = 1.0/step_hour
  steps_per_day = int(24.0*steps_per_hour)
  datalines = N_yrs * 365 * int(24.0/step_hour)
  days_data = N_yrs * 365
  yr_data   = N_yrs
  call ReadNCfiles(ncfilepath, ncfields, yr_start, yr_end)

  !------------ Forcing data interpolation and model run
  !$omp parallel do private(i, GPP, NPP, Tr) shared(climate, veg_state, fluxes)
  if(start_grid > N_VegGrids)then
    print '(A, I8, I8)', "start_grid is beyond the total grids: ", start_grid, N_VegGrids
    stop
  endif

  do m = start_grid, N_VegGrids
    GridID = GridLonLat(m) ! for file names
    write(*,'(a30,3(I8,","))')'Running at grid (iLon, iLat): ', &
               GridLonLat(m), CRUgrid(m)%iLon, CRUgrid(m)%iLat
    print '(A, I8, A, I8)', 'Grid ', m, ' of ', N_VegGrids
    ! Data interpolated to hourly
    call CRU_Interpolation(CRUgrid(m),steps_per_hour,forcingData)
    call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
    call BiomeE_main()
    call zip_output_files()

    ! ---------- Time stamp -------------
    call cpu_time(end_time)
    elapsed_time = end_time - last_time
    write(*,'(a22,3(f12.1,","))')'Time taken (seconds): ', elapsed_time
    last_time = end_time
  enddo
  !$omp end parallel do

  ! Release netcdf-related allocatable data arrays
  deallocate(tswrfH)
  !deallocate(CRUData)
  deallocate(GridClimateData, CRUgrid)
  deallocate(GridLonLat)

#else
  ! Read in namelist file and setup output files
  call read_namelist(fnml) ! Namelist file (must be hardwired)
  ! Single site run
  call setup_forcingdata()
  call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
  call BiomeE_main()
#endif

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(a24,3(f12.1,","))')'Time taken (seconds):', elapsed_time

end program BiomeE_global_driver
