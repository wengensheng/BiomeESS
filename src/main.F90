! ------------ Model driver -------------
! ------------- 08/30/2022 --------------

program BiomeE
  use BiomeE_mod, only: BiomeE_main
  !-----------
  real :: start_time, end_time, elapsed_time
  integer :: timeArray(3)

  ! ---------- Time stamp -------------
  call cpu_time(start_time)
  call itime(timeArray)     ! Get current time
  write(*,'(a12,3(I2,":"))')'Start time: ', timeArray

  ! ---------- Model Run -------------
  call BiomeE_main()

  ! ---------- Time stamp -------------
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  write(*,'(a24,3(f8.2,","))')'Time taken (seconds):', elapsed_time
end program BiomeE
