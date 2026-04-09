! ------------ Model driver -------------
! ------------- 08/30/2022 --------------
! BiomeE global version 0.1, Ensheng Weng, 07/18/2025, almost done, 01/15/2026

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
  character(len=256) :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'
  real :: start_time, end_time, last_time, elapsed_time
  integer :: timeArray(3), m
  ! Wall time counting
  INTEGER(kind=8) :: start_count, end_count, count_rate, count_max
  REAL(kind=8)    :: wall_time
  logical :: file_exists

  ! If a namelist file is provided on the command line, use it; otherwise use default.
  if (command_argument_count() >= 1) then
    call get_command_argument(1, fnml)
    fnml = adjustl(fnml)
  endif
  inquire(file=trim(fnml), exist=file_exists)
  if (.not. file_exists) then
    write(*,*) 'Namelist file does not exist: ', trim(fnml)
    stop
  endif


  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  call itime(timeArray)     ! Get current time
  write(*,'(a20,3(I2,":"))')'Model run start at: ', timeArray

  ! Get initial clock information
  CALL SYSTEM_CLOCK(COUNT_RATE=count_rate, COUNT_MAX=count_max)
  ! Get the starting time
  CALL SYSTEM_CLOCK(COUNT=start_count)

  ! Read in model setting and control namelist (for both global and site runs)
  call read_init_namelist(fnml)

#ifdef GlobalRun
  ! ---------- Global or regional run with netcdf input files -----
  call read_global_setting(fnml) ! Read in global data and model run settings
  call Preset_GlobalPFTs()   ! Hardwired global pft parameters
  call model_para_init(fnml) ! Update hardwired global PFT & soil parameters from namelists

#ifdef Use_InterpolatedData
  call read_GridLonLat(GridListFile,file_exists)
  if(.not. file_exists)then
    write(*,*)trim(GridListFile), 'does not exist. Stopped!'
    stop !
  endif
#else
  ! Read in netCDF global data files
  call ReadNCfiles(ncfilepath, veg_path, ndp_path)
#endif

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(A,3(f8.2,","))')'Data reading CPU time (minutes): ', elapsed_time/60.

  ! Get the ending time
  CALL SYSTEM_CLOCK(COUNT=end_count)

  ! Calculate wall time
  IF (end_count < start_count) THEN
    ! Handle clock wrap-around if it occurs
    wall_time = REAL(end_count + count_max - start_count, kind=8) / count_rate/60.0
  ELSE
    wall_time = REAL(end_count - start_count, kind=8) / count_rate / 60.0
  END IF
  PRINT *, "Data reading wall time:", wall_time, " minutes"

  !------------ Forcing data interpolation and model run ---------------------------
  !$omp parallel do private(GridID,forcingData,fno1,fno2,fno3,fno4,fno5,fno6) shared(GridLonLat, LandGrid)
  do m = grid_No1, grid_No2  ! Grids in GridLonLat
    ! Get this grid's ID
    GridID = GridLonLat(m) ! for file names
    GridVC = GridVegCov(:,m)
    GridFR = GridFarm(:,m)
    call cpu_time(last_time) ! Record time needed for one grid simulation
    print '(A, I6, A, I6)', 'Working at grid: ', GridID, '. Grid No. ', m
    print '(A, I6, A, I6)', 'The ', m - grid_No1 + 1, 'th grid of ', grid_No2 - grid_No1 + 1

    ! --------- Get the forcingData for this grid -----------
#ifdef Use_InterpolatedData
    ! Read interpolated data from disk files. Moved here to avoid reading errors
    call read_interpolatedCRU(int_fpath,int_prefix,GridID,yr_start,yr_end,forcingData,file_exists)
    if(.not. file_exists)then
      print '(A, I8, A)', 'Grid ', GridID, ' is skipped b/c of no input file or shorter than needed.'
      cycle
    endif
#else
    ! Interpolate grid data to hourly
    call CRU_Interpolation(LandGrid(m),forcingData)
    if(WriteForcing)then
      deallocate(forcingData)
      cycle ! Skip model run. Only output interpolated grid forcing data
    endif
#endif

    ! ------ Set up output files for this grid ----------------
    fno1 = GridID + 1000000; fno2 = GridID + 2000000
    fno3 = GridID + 3000000; fno4 = GridID + 4000000
    fno5 = GridID + 5000000; fno6 = GridID + 6000000
    call setup_output_files() ! Setup output files before reading forcing data (?)

    ! ------- Run model -----------
    call BiomeE_main()

#ifdef Zip_outputs
    call zip_output_files()
#endif
    ! ---------- Time stamp -------------
    call cpu_time(end_time)
    elapsed_time = end_time - last_time
    write(*,'(a6,f6.1, a12, I6, a17, f7.1, a9)') &
        'Used: ', elapsed_time, ' seconds at ', GridID, &
        '. Total CPU time: ', (end_time-start_time)/60.0, ' minutes.'

    last_time = end_time
  enddo
  !$omp end parallel do

  ! Release netcdf-related allocatable data arrays and close files
  call CRU_end()

#else

  ! ---------- Single site run with csv/txt forcing data input ----------
  call model_para_init(fnml) ! Read in PFT & soil parameters from namelists
  call setup_forcingdata(climfile)
  call setup_output_files()
  call BiomeE_main()

#endif

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(A,3(f7.2,","))')'Total CPU time (minutes): ', elapsed_time/60.

  ! Get the ending time
  CALL SYSTEM_CLOCK(COUNT=end_count)
  ! Calculate wall time
  IF (end_count < start_count) THEN
    ! Handle clock wrap-around if it occurs
    wall_time = REAL(end_count + count_max - start_count, kind=8) / count_rate/60.0
  ELSE
    wall_time = REAL(end_count - start_count, kind=8) / count_rate / 60.0
  END IF
  PRINT '(A, f7.2, A, f6.2, A)', "Total wall time:", wall_time, " minutes (", wall_time/60., ' hours).'

end program BiomeE
