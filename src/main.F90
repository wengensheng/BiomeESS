! ------------ Model driver -------------
! ------------- 08/30/2022 --------------

program BiomeE
  use BiomeE_mod, only: BiomeE_main
  use io_mod, only: read_namelist,setup_forcingdata,setup_output_files
  use datatypes

  implicit none
  !-----------
  character(len=80)  :: fnml = './para_files/input.nml' ! 'parameters_ORNL_test.nml'
  real :: start_time, end_time, elapsed_time
  integer :: timeArray(3)

  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  call itime(timeArray)     ! Get current time
  write(*,'(a12,3(I2,":"))')'Start time: ', timeArray

  ! Input and output file setting, and model run
  call read_namelist(fnml) ! Namelist file (must be hardwired)
  call setup_forcingdata()
  call setup_output_files(fno1,fno2,fno3,fno4,fno5,fno6)
  call BiomeE_main()

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(a24,3(f8.2,","))')'Time taken (seconds):', elapsed_time
end program BiomeE
