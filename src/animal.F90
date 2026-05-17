!======================================================================
! animal.F90 - Animal Functional Type (AFT) dynamics for BiomeE
!
! Follows vegetation.F90 conventions:
!   - AFT parameters live in aftdata(0:N_AFT), populated by initialize_AFT_pars
!   - Animal cohorts live in vegn%ani_cohorts(:), allocated by initialize_ani_tile
!   - Subroutines address aftdata via  associate(sp => aftdata(ac%aft))
!
! SOC/SON pool routing (mirrors plant2soil in vegetation.F90):
!   Feces, undigested prey remains, and carcasses all go to SOC/SON(4) fast SOM 
!
! Calling sequence in main loop:
!   call initialize_AFT_pars()             ! once at model startup
!   call initialize_ani_tile(vegn)         ! once per tile
!   ...
!   call ani_daily_update(vegn, deltat)    ! each day (after vegn_daily_update)
!   ...
!   call ani_annual_diagnostics(vegn)      ! once per year
!
! Weng Lab - 2026
!======================================================================

#ifdef DO_ANIMAL
module animal_mod
  use datatypes
  implicit none
  private

  public :: initialize_AFT_pars
  public :: initialize_ani_tile
  public :: ani_daily_update
  public :: ani_reproduction

contains

!=============================================================================
subroutine initialize_AFT_pars()
  ! Copy module-level AFT parameter arrays -> aftdata(:).
  ! Analogous to initialize_PFT_pars / spdata population in model_utils.F90.
  ! Call once at model startup, after reading ani_parameters_nml.
  implicit none
  integer :: i

  aftdata(:)%diet_class    = aft_diet_class
  aftdata(:)%body_mass     = aft_body_mass
  aftdata(:)%f_C_body      = aft_f_C_body
  aftdata(:)%f_N_body      = aft_f_N_body
  aftdata(:)%I_max_plant   = aft_I_max_plant
  aftdata(:)%K_half_plant  = aft_K_half_plant
  aftdata(:)%browse_height = aft_browse_height
  aftdata(:)%f_plant_diet  = aft_f_plant_diet
  aftdata(:)%litter_pref   = aft_litter_pref
  aftdata(:)%I_max_prey    = aft_I_max_prey
  aftdata(:)%K_half_prey   = aft_K_half_prey
  aftdata(:)%f_prey_diet   = aft_f_prey_diet
  aftdata(:)%digestibility = aft_digestibility
  aftdata(:)%f_N_feces     = aft_f_N_feces
  aftdata(:)%I_maint       = aft_I_maint
  aftdata(:)%mu_starve_max = aft_mu_starve_max
  aftdata(:)%mu_background = aft_mu_background
  aftdata(:)%r_max         = aft_r_max
  ! Palatability: copy each AFT column from the 2-D module array
  ! aft_palatability is (0:MSPECIES, 0:N_AFT) so the first index is species
  do i = 0, N_AFT
    aftdata(i)%palatability = aft_palatability(:, i)
  end do

end subroutine initialize_AFT_pars

!=============================================================================
subroutine initialize_ani_tile(vegn)
  ! Allocate vegn%ani_cohorts and set initial state from init_ani_cohort_* vars.
  ! Analogous to initialize_vegn_tile for plant cohorts.
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn
  integer :: i, n

  n = init_ani_cohort_N
  vegn%n_ani_cohorts = n
  if (n <= 0) return

  allocate(vegn%ani_cohorts(n))
  do i = 1, n
    vegn%ani_cohorts(i)%aft     = init_ani_cohort_aft(i)
    vegn%ani_cohorts(i)%nindivs = init_ani_cohort_nindivs(i)
    vegn%ani_cohorts(i)%age     = 0.0
    vegn%ani_cohorts(i)%body_C  = aftdata(init_ani_cohort_aft(i))%body_mass &
                                * aftdata(init_ani_cohort_aft(i))%f_C_body
  end do

end subroutine initialize_ani_tile

!=============================================================================
subroutine ani_daily_update(vegn, deltat)
  ! Advance all animal cohorts one day. Analogous to vegn_daily_update.
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn
  real,                 intent(in)    :: deltat  ! dt_daily_yr = 1/365

  integer :: i

  do i = 1, vegn%n_ani_cohorts
    if (vegn%ani_cohorts(i)%nindivs <= 0.0) cycle
    call ani_cohort_daily(vegn%ani_cohorts(i), vegn, deltat)
  end do

end subroutine ani_daily_update

!=============================================================================
subroutine ani_cohort_daily(ac, vegn, deltat)
  ! One daily step for a single animal cohort.
  implicit none
  type(ani_cohort_type), target, intent(inout) :: ac
  type(vegn_tile_type),          intent(inout) :: vegn
  real,                          intent(in)    :: deltat

  real :: intake_plant_ind, intake_prey_ind
  real :: intake_plant_tot, intake_prey_tot
  real :: intake_eff, C_removed
  real :: C_feces, N_feces, C_carcass, N_carcass
  real :: mu_starve, mu_total, dN

  associate (sp => aftdata(ac%aft))

  !------------------------------------------------------------------
  ! 1. Plant feeding (herbivore / omnivore)
  !------------------------------------------------------------------
  if (sp%f_plant_diet > 0.0 .and. sp%I_max_plant > 0.0) then
    call ani_plant_intake(ac, vegn, intake_plant_ind, intake_plant_tot, C_removed)
  else
    intake_plant_ind = 0.0 ; intake_plant_tot = 0.0 ; C_removed = 0.0
  end if

  !------------------------------------------------------------------
  ! 2. Prey feeding (carnivore / omnivore)
  !------------------------------------------------------------------
  if (sp%f_prey_diet > 0.0 .and. sp%I_max_prey > 0.0) then
    call ani_prey_intake(ac, vegn, intake_prey_ind, intake_prey_tot)
  else
    intake_prey_ind = 0.0 ; intake_prey_tot = 0.0
  end if

  !------------------------------------------------------------------
  ! 3. Feces -> metabolic litter SOC(4) / SON(4)
  !    Both plant and prey intake are in kg C m-2 day-1.
  !    Undigested fraction = (1 - digestibility) * intake.
  !    f_N_feces is N:C ratio of feces (kg N kg-1 C).
  !------------------------------------------------------------------
  C_feces = (intake_plant_tot + intake_prey_tot) * (1.0 - sp%digestibility)
  N_feces = intake_plant_tot * sp%f_N_feces &
          + intake_prey_tot  * (1.0 - sp%digestibility) &
          * safe_div(sp%f_N_body, sp%f_C_body)
  vegn%SOC(4) = vegn%SOC(4) + C_feces
  vegn%SON(4) = vegn%SON(4) + N_feces

  !------------------------------------------------------------------
  ! 4. Starvation mortality
  !    intake_eff: individual C intake rate, kg C ind-1 day-1
  !    maintenance threshold: I_maint * body_C
  !------------------------------------------------------------------
  select case (sp%diet_class)
    case (AFT_HERBIVORE)
      intake_eff = intake_plant_ind
    case (AFT_CARNIVORE)
      intake_eff = intake_prey_ind
    case default   ! omnivore: sum of plant and prey C
      intake_eff = intake_plant_ind + intake_prey_ind
  end select

  if (sp%I_maint > 0.0) then
    mu_starve = sp%mu_starve_max * max(0.0, 1.0 - intake_eff / &
                (sp%I_maint * ac%body_C))
  else
    mu_starve = 0.0
  end if
  mu_total = sp%mu_background + mu_starve * (1.0 - sp%mu_background)

  dN = min(ac%nindivs * (1.0 - exp(-mu_total)), ac%nindivs)
  ac%nindivs = max(0.0, ac%nindivs - dN)

  !------------------------------------------------------------------
  ! 5. Carcasses -> fast SOM, SOC(4) / SON(4)
  !------------------------------------------------------------------
  C_carcass = dN * ac%body_C
  N_carcass = dN * sp%body_mass * sp%f_N_body
  vegn%SOC(4) = vegn%SOC(4) + C_carcass
  vegn%SON(4) = vegn%SON(4) + N_carcass

  !------------------------------------------------------------------
  ! 6. Store daily fluxes and advance age
  !------------------------------------------------------------------
  ac%intake_plant  = intake_plant_ind
  ac%intake_prey   = intake_prey_ind
  ac%C_removed_veg = C_removed
  ac%C_feces       = C_feces
  ac%N_feces       = N_feces
  ac%C_carcass     = C_carcass
  ac%N_carcass     = N_carcass
  ac%mu_starve     = mu_starve
  ac%deaths        = dN
  ac%age           = ac%age + deltat

  !------------------------------------------------------------------
  ! 7. Accumulate annual diagnostics
  !------------------------------------------------------------------
  ac%annualIntakePlant = ac%annualIntakePlant + intake_plant_ind
  ac%annualIntakePrey  = ac%annualIntakePrey  + intake_prey_ind
  ac%annualC_removed   = ac%annualC_removed   + C_removed
  ac%annualC_feces     = ac%annualC_feces     + C_feces
  ac%annualN_feces     = ac%annualN_feces     + N_feces
  ac%annualC_carcass   = ac%annualC_carcass   + C_carcass
  ac%annualN_carcass   = ac%annualN_carcass   + N_carcass

  end associate

end subroutine ani_cohort_daily

!=============================================================================
subroutine ani_plant_intake(ac, vegn, intake_ind, intake_tot, C_removed_tot)
  ! Type II functional response on plant leaf C + metabolic litter C (SOC pool 1).
  ! All biomass quantities are in kg C m-2 (no DM conversion needed).
  ! Leaf access is height-limited: fraction = browse_height / (height + browse_height).
  ! Litter access is controlled by sp%litter_pref (0 = no litter eating).
  ! C removed from cc%bl and vegn%SOC(1) is distributed by palatability-weighted forage.
  ! cc%brsC and cc%brsN record the daily leaf C and N eaten from each plant cohort.
  implicit none
  type(ani_cohort_type), target, intent(inout) :: ac
  type(vegn_tile_type),          intent(inout) :: vegn
  real,                          intent(out)   :: intake_ind    ! kg C ind-1 day-1
  real,                          intent(out)   :: intake_tot    ! kg C m-2  day-1
  real,                          intent(out)   :: C_removed_tot ! kg C m-2  day-1 (leaves only)

  type(cohort_type),     pointer :: cc
  real :: B_cc(vegn%n_cohorts)   ! available leaf C per plant cohort, kg C m-2
  real :: w_cc(vegn%n_cohorts)   ! palatability-weighted leaf C per cohort, kg C m-2
  real :: B_vegn, W_vegn
  real :: B_SOM1, W_SOM1         ! litter C and weighted litter C, kg C m-2
  real :: B_avail, w_total       ! total available C and weighted C, kg C m-2
  real :: C_removed_cc, C_removed_all, C_removed_litter, frac_remain, K_plant
  real :: N_removed_cc
  integer :: i

  associate (sp => aftdata(ac%aft))

  B_vegn = 0.0 ; w_vegn = 0.0
  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    ! cc%bl is kg C ind-1; scale by density and height-limited browse fraction -> kg C m-2
    B_cc(i) = cc%bl * cc%nindivs * sp%browse_height / (cc%height + sp%browse_height)
    w_cc(i) = sp%palatability(cc%species) * B_cc(i)
    B_vegn  = B_vegn + B_cc(i)
    w_vegn  = w_vegn + w_cc(i)
  end do

  !-- Add metabolic litter (SOC pool 1) as supplemental forage --
  ! vegn%SOC(1) is already kg C m-2
  B_SOM1  = vegn%SOC(1)
  W_SOM1  = sp%litter_pref * B_SOM1
  B_avail = B_vegn + B_SOM1
  w_total = w_vegn + W_SOM1

  !-- Type II functional response: B_avail in kg C m-2, intake in kg C ind-1 day-1 --
  ! I_max_plant is fraction of body_C day-1
  K_plant = B_avail / (sp%K_half_plant + B_avail)
  intake_ind    = sp%f_plant_diet * sp%I_max_plant * ac%body_C * K_plant
  intake_tot    = intake_ind * ac%nindivs        ! kg C m-2 day-1
  C_removed_all = min(intake_tot, B_avail)       ! kg C m-2 day-1

  !-- Distribute C removal to plant cohort leaves --
  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    if (w_total > 0.0 .and. w_cc(i) > 0.0) then
      C_removed_cc = C_removed_all * (w_cc(i) / w_total)
    else
      C_removed_cc = 0.0
    end if
    ! cap at available leaf C for this cohort
    C_removed_cc = min(C_removed_cc, B_cc(i))
    if (cc%nindivs > 0.0 .and. C_removed_cc > 0.0) then
      frac_remain = max(0.0, 1.0 - C_removed_cc / max(cc%bl * cc%nindivs, 1.e-12))
      ! Record leaf C and N eaten from this cohort (kg m-2 day-1) before scaling
      N_removed_cc = safe_div(cc%leafN, cc%bl) * C_removed_cc
      cc%brsC = cc%brsC   + C_removed_cc
      cc%brsN = cc%brsN   + N_removed_cc
      cc%bl    = cc%bl    * frac_remain
      cc%leafN = cc%leafN * frac_remain
    end if
  end do

  !-- Distribute C removal to metabolic litter --
  if (W_SOM1 > 0.0) then
    C_removed_litter = min(C_removed_all * (W_SOM1/w_total), vegn%SOC(1))
    vegn%SOC(1) = max(0.0, vegn%SOC(1) - C_removed_litter)
  else
    C_removed_litter = 0.0
  end if

  C_removed_tot = C_removed_all - C_removed_litter  ! leaf-only C removal returned to caller

  end associate

end subroutine ani_plant_intake

!=============================================================================
subroutine ani_prey_intake(ac, vegn, intake_ind, intake_tot)
  ! Type II functional response on total body C of all OTHER animal cohorts.
  ! All biomass quantities in kg C m-2.
  ! Prey density is reduced proportionally to prey C removed.
  implicit none
  type(ani_cohort_type), target, intent(inout) :: ac
  type(vegn_tile_type),          intent(inout) :: vegn
  real,                          intent(out)   :: intake_ind   ! kg C ind-1 day-1
  real,                          intent(out)   :: intake_tot   ! kg C m-2  day-1

  type(ani_cohort_type), pointer :: ac_j
  real :: B_prey(vegn%n_ani_cohorts)
  real :: B_prey_total, C_eaten_prey, C_removed_cc, frac_remain, K_prey
  integer :: j

  associate (sp => aftdata(ac%aft))

  B_prey_total = 0.0 ; B_prey = 0.0
  do j = 1, vegn%n_ani_cohorts
    ac_j => vegn%ani_cohorts(j)
    if (associated(ac_j, ac) .or. ac_j%nindivs <= 0.0) cycle
    ! Prey body C density: nindivs * body_C, kg C m-2
    B_prey(j)    = ac_j%nindivs * ac_j%body_C
    B_prey_total = B_prey_total + B_prey(j)
  end do

  if (B_prey_total > 0.0) then
    ! I_max_prey is fraction of body_C day-1 -> kg C ind-1 day-1
    K_prey = B_prey_total / (sp%K_half_prey + B_prey_total)
    intake_ind = sp%f_prey_diet * sp%I_max_prey * ac%body_C * K_prey
  else
    intake_ind = 0.0
  end if
  intake_tot   = intake_ind * ac%nindivs       ! kg C m-2 day-1
  C_eaten_prey = min(intake_tot, B_prey_total) ! kg C m-2 day-1 actually consumed

  do j = 1, vegn%n_ani_cohorts
    ac_j => vegn%ani_cohorts(j)
    if (associated(ac_j, ac) .or. ac_j%nindivs <= 0.0 .or. B_prey_total <= 0.0) cycle
    ! C removed from prey cohort j, proportional to its share of prey C
    C_removed_cc = min(C_eaten_prey * safe_div(B_prey(j), B_prey_total), B_prey(j))
    frac_remain  = max(0.0, 1.0 - C_removed_cc / max(B_prey(j), 1.e-12))
    ac_j%nindivs = ac_j%nindivs * frac_remain
  end do

  end associate

end subroutine ani_prey_intake

!=============================================================================
subroutine ani_reproduction(vegn)
  ! Annual reproduction: increase nindivs for each cohort in proportion to
  ! food surplus above maintenance. No new cohorts are created; newborns join
  ! the existing cohort at adult body size (density-only increase).

  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  integer :: i
  real    :: annual_intake_ind, r_repro, C_maint

  do i = 1, vegn%n_ani_cohorts
    associate (ac => vegn%ani_cohorts(i), sp => aftdata(vegn%ani_cohorts(i)%aft))
    ! Migration when individuals are low
    if( ac%nindivs < zero_thld) ac%nindivs = 1.0E-7

    ! Calculate reproduction
    if (sp%r_max > 0.0) then
      ! Per-individual annual intake in the same units as I_maint
      select case (sp%diet_class)
        case (AFT_HERBIVORE)
          annual_intake_ind = ac%annualIntakePlant
        case (AFT_CARNIVORE)
          annual_intake_ind = ac%annualIntakePrey
        case default  ! omnivore
          annual_intake_ind = ac%annualIntakePlant + ac%annualIntakePrey
      end select

      ! Reproduction rate: surplus C intake above maintenance, scaled by body C
      ! I_maint * body_C = daily maintenance C threshold per individual
      C_maint = sp%I_maint * ac%body_C * days_per_year
      r_repro = max(0.0, min(sp%r_max, (annual_intake_ind - C_maint) / ac%body_C))
      ac%nindivs = (1.0 + r_repro) * ac%nindivs
      write(*,*)'annual_intake,cost, r_repro,nindivs',annual_intake_ind,C_maint, r_repro,ac%nindivs

    end if

    end associate
  end do

end subroutine ani_reproduction

!=============================================================================
!  Internal helpers
!=============================================================================

pure real function safe_div(a, b) result(r)
  real, intent(in) :: a, b
  if (abs(b) > 1.e-30) then
    r = a / b
  else
    r = 0.0
  end if
end function safe_div

!=============================================================================
end module animal_mod
#endif /* DO_ANIMAL */
