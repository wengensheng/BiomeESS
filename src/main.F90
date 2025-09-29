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
  integer :: timeArray(3), m
  logical :: file_exists

  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  call itime(timeArray)     ! Get current time
  write(*,'(a20,3(I2,":"))')'Model run start at: ', timeArray

  ! Read in model initialization and plant & soil parameters
  call model_para_init(fnml)

#ifdef GlobalRun
  ! ---------- Global or regional run with netcdf input files -----
  ! Read in the namelist for global data and model settings
  call read_global_setting(fnml)

  ! Read in netCDF global data files
  call ReadNCfiles(ncfilepath, ncfields, yr_start, yr_end)

  !------------ Forcing data interpolation and model run
  !$omp parallel do private(GridID,forcingData,fno1,fno2,fno3,fno4,fno5,fno6) shared(GridLonLat, LandGrid)
  do m = start_grid, N_VegGrids
    call cpu_time(last_time) ! Record time needed for one grid simulation

    ! ------ Get a grid's forcingData
    GridID = GridLonLat(m) ! for file names
#ifdef Use_InterpolatedData
    call read_interpolatedCRU(GridID,yr_start,yr_end,forcingData,file_exists)
    if(.not. file_exists)then
      print '(A, I8, A)', 'Grid ', GridID, ' is skipped b/c of no input file'
      cycle
    endif
#else
    ! Data interpolated to hourly
    call CRU_Interpolation(LandGrid(m),forcingData)
    if(WriteForcing)then
      deallocate(forcingData)
      cycle ! Skip model runs
    endif
#endif

    ! Output file grid ID
    fno1=GridID+1; fno2=GridID+2; fno3=GridID+3
    fno4=GridID+4; fno5=GridID+5; fno6=GridID+6
    write(*,'(a20,3(I6,","))')'Running at grid: ', GridID
    print '(A, I8, A, I8)', 'Grid ', m, ' of ', N_VegGrids

    ! ------- Run model -----------
    call setup_output_files()
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
  deallocate(GridLonLat)
#ifndef Use_InterpolatedData
  !deallocate(CRUData)
  deallocate(CRUtime)
  deallocate(ClimData)
  deallocate(LandGrid)
#endif

#else
  ! ---------- Single site run with csv/txt forcing data input ----------
  call setup_forcingdata(climfile)
  call setup_output_files()
  call BiomeE_main()
#endif
  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(a24,3(f8.2,","))')'Total time (minutes):', elapsed_time/60.

end program BiomeE
