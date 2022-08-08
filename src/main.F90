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
!             (Lase updated 12/30/2017, 07/24/2022)
!
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

program BiomeE
  use BiomeE_mod

  !-----------
  call BiomeE_initialization()
  call BiomeE_run()
  call BiomeE_end()
end program BiomeE
