! BiomeE global run
! Ensheng Weng, 07/18/2025

program BiomeE_global_driver
  use BiomeE_mod, only: BiomeE_main
  use io_mod, only: read_namelist,setup_forcingdata, &
                    setup_output_files, zip_output_files
  use netcdf_io
  use datatypes

  implicit none

  !-----------
  !type(climate_data_type), pointer :: forcing1(:),forcing2(:),forcing3(:)
  character(len=80)  :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'
  real :: start_time, end_time, elapsed_time
  real, pointer :: GridData(:,:)    ! N_yr*Ntime, N_vars
  integer :: N_yrs, totL
  integer :: timeArray(3)

  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  call itime(timeArray)     ! Get current time
  write(*,'(a12,3(I2,":"))')'Start time: ', timeArray

  ! Read in namelist file and setup output files
  call read_namelist(fnml) ! Namelist file (must be hardwired)

#ifdef GlobalRun
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
  do iLon = LowerLon, UpperLon
    do iLat = LowerLat, UpperLat
      if(CRUData(1, 3, iLon, iLat)< 9999.0)then
        write(*,*)'model run at grid (iLon, iLat): ', iLon, iLat
        GridID = iLon * 1000 + iLat
        GridData => CRUData(:,:,iLon, iLat)
        call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
        call CRU_Interpolation(GridData,iLon,iLat,steps_per_hour,forcingData)
        call BiomeE_main()
        call zip_output_files()
      endif
    enddo
  enddo

  ! Release netcdf-related allocatable data arrays
  deallocate(tswrfH)
  deallocate(CRUData)

#else
  ! Single site run
  call setup_forcingdata()
  call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
  call BiomeE_main()
#endif

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(a24,3(f8.2,","))')'Time taken (seconds):', elapsed_time

end program BiomeE_global_driver
