! ------------ Model driver -------------
! ------------- 08/30/2022 --------------
! BiomeE global version 0.1, Ensheng Weng, 07/18/2025

program BiomeE
  !use omp_lib
  use datatypes
#ifdef GlobalRun
  use netcdf_io
#endif
  use io_mod, only: setup_forcingdata,setup_output_files,zip_output_files
  use BiomeE_mod, only: BiomeE_main

  implicit none
  !-----------
  character(len=80)  :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'
  real :: start_time, end_time, last_time, elapsed_time
  real, pointer :: GridData(:,:)    ! N_yr*Ntime, N_vars
  integer :: N_yrs, totL
  integer :: m, n
  integer :: timeArray(3)

  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  call itime(timeArray)     ! Get current time
  write(*,'(a20,3(I2,":"))')'Model run start at: ', timeArray

  ! Read in model initialization and plant & soil parameters
  call read_init_namelist(fnml)
  call initialize_soilpars(fnml)
  call initialize_PFT_data(fnml)

#ifdef GlobalRun
  ! ---------- Global or regional run with netcdf input files -----
  ! Read in the namelist for global data and model settings
  call read_global_setting(fnml)

  ! ------ Setup steps for model run ------
  N_yrs = yr_end - yr_start + 1
  step_hour = 1.0 ! hourly time step
  steps_per_hour = 1.0/step_hour
  steps_per_day = int(24.0*steps_per_hour)
  datalines = N_yrs * 365 * int(24.0/step_hour)
  days_data = N_yrs * 365
  yr_data   = N_yrs

  ! Read in netCDF global data files
  call ReadNCfiles(ncfilepath, ncfields, yr_start, yr_end)

  !------------ Forcing data interpolation and model run
  !$omp parallel do private(GridID,forcingData,fno1,fno2,fno3,fno4,fno5,fno6) shared(GridLonLat, LandGrid)
  do m = start_grid, N_VegGrids
    GridID = GridLonLat(m) ! for file names
    fno1=GridID+1; fno2=GridID+2; fno3=GridID+3
    fno4=GridID+4; fno5=GridID+5; fno6=GridID+6
    write(*,'(a30,3(I8,","))')'Running at grid (iLon, iLat): ', &
               GridLonLat(m), LandGrid(m)%iLon, LandGrid(m)%iLat
    print '(A, I8, A, I8)', 'Grid ', m, ' of ', N_VegGrids
    call cpu_time(last_time) ! Record time needed for one grid simulation
    ! Data interpolated to hourly
    call CRU_Interpolation(LandGrid(m),steps_per_hour,forcingData)
    call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
    call Set_PFTs_from_Climate(forcingData,datalines,days_data,yr_data)
    !cycle ! for testing set_pfts_from_climate
    call BiomeE_main()
    call zip_output_files()

    ! ---------- Time stamp -------------
    call cpu_time(end_time)
    elapsed_time = end_time - last_time
    write(*,'(a6,f6.1, a12, I6, a14, f7.1, a9)') &
            'Used: ', elapsed_time, ' seconds at ', GridID, &
            '. Total time: ', (end_time-start_time)/60.0, ' minutes.'

    last_time = end_time
  enddo
  !$omp end parallel do

  ! Release netcdf-related allocatable data arrays
  deallocate(tswrfH)
  !deallocate(CRUData)
  deallocate(ClimData)
  deallocate(LandGrid)
  deallocate(GridLonLat)
#else
  ! ---------- Single site run with csv/txt forcing data input ----------
  call setup_forcingdata()
  call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
  call BiomeE_main()

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(a24,3(f8.2,","))')'Total time (minutes):', elapsed_time/60.
#endif

end program BiomeE
