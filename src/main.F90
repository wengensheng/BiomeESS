! ------------ Model driver -------------
! ------------- 08/30/2022 --------------

program BiomeE
  use BiomeE_mod
  !-----------
  call BiomeE_initialization()
  call BiomeE_run()
  call BiomeE_end()
end program BiomeE
