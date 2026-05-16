!======================================================================
! animal.F90 - Animal Functional Type (AFT) dynamics for BiomeE
!
! Follows vegetation.F90 conventions:
!   - AFT parameters live in aftdata(0:N_AFT), populated by initialize_AFT_pars
!   - Animal cohorts live in vegn%ani_cohorts(:), allocated by initialize_ani_tile
!   - Subroutines address aftdata via  associate(sp => aftdata(ac%aft))
!
! SOC/SON pool routing (mirrors plant2soil in vegetation.F90):
!   SOC/SON(1) metabolic litter  -- feces and undigested prey remains
!   SOC/SON(2) structural litter -- carcasses (slow-decomposing)
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
  aftdata(:)%DM_to_C       = aft_DM_to_C
  aftdata(:)%digestibility = aft_digestibility
  aftdata(:)%f_N_feces     = aft_f_N_feces
  aftdata(:)%I_maint       = aft_I_maint
  aftdata(:)%mu_starve_max = aft_mu_starve_max
  aftdata(:)%mu_background = aft_mu_background
  aftdata(:)%r_max         = aft_r_max
  ! f_C_feces derived: undigested fraction times DM-to-C conversion
  aftdata(:)%f_C_feces = (1.0 - aft_digestibility) * aft_DM_to_C
  ! Palatability: copy each AFT row from the 2-D module array
  do i = 0, N_AFT
    aftdata(i)%palatability = aft_palatability(i, :)
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
    call ani_cohort_daily(i, vegn, deltat)
  end do

end subroutine ani_daily_update

!=============================================================================
subroutine ani_cohort_daily(idx, vegn, deltat)
  ! One daily step for cohort idx in vegn%ani_cohorts.
  implicit none
  integer,               intent(in)    :: idx
  type(vegn_tile_type),  intent(inout) :: vegn
  real,                  intent(in)    :: deltat

  type(ani_cohort_type), pointer :: ac
  real :: intake_plant_ind, intake_prey_ind
  real :: intake_plant_tot, intake_prey_tot
  real :: intake_eff, C_removed
  real :: C_feces, N_feces, C_carcass, N_carcass
  real :: mu_starve, mu_total, dN

  ac => vegn%ani_cohorts(idx)

  associate (sp => aftdata(ac%aft))

  !------------------------------------------------------------------
  ! 1. Plant feeding (herbivore / omnivore)
  !------------------------------------------------------------------
  if (sp%f_plant_diet > 0.0 .and. sp%I_max_plant > 0.0) then
    call ani_plant_intake(idx, vegn, intake_plant_ind, intake_plant_tot, C_removed)
  else
    intake_plant_ind = 0.0 ; intake_plant_tot = 0.0 ; C_removed = 0.0
  end if

  !------------------------------------------------------------------
  ! 2. Prey feeding (carnivore / omnivore)
  !------------------------------------------------------------------
  if (sp%f_prey_diet > 0.0 .and. sp%I_max_prey > 0.0) then
    call ani_prey_intake(idx, vegn, intake_prey_ind, intake_prey_tot)
  else
    intake_prey_ind = 0.0 ; intake_prey_tot = 0.0
  end if

  !------------------------------------------------------------------
  ! 3. Feces -> metabolic litter SOC(4) / SON(4)
  !    Plant-derived: undigested DM fraction
  !    Prey-derived:  undigested prey body C and N
  !------------------------------------------------------------------
  C_feces = intake_plant_tot * sp%f_C_feces &
          + intake_prey_tot  * (1.0 - sp%digestibility)
  N_feces = intake_plant_tot * sp%f_N_feces &
          + intake_prey_tot  * (1.0 - sp%digestibility) &
            * safe_div(sp%f_N_body, sp%f_C_body)
  vegn%SOC(4) = vegn%SOC(4) + C_feces
  vegn%SON(4) = vegn%SON(4) + N_feces

  !------------------------------------------------------------------
  ! 4. Starvation mortality
  !    intake_eff is in the same units as I_maint for this diet class
  !------------------------------------------------------------------
  select case (sp%diet_class)
    case (AFT_HERBIVORE)
      intake_eff = intake_plant_ind
    case (AFT_CARNIVORE)
      intake_eff = intake_prey_ind
    case default   ! omnivore: weighted sum, convert prey C to DM equiv.
      intake_eff = sp%f_plant_diet * intake_plant_ind &
                 + sp%f_prey_diet  * safe_div(intake_prey_ind, sp%DM_to_C)
  end select

  if (sp%I_maint > 0.0) then
    mu_starve = sp%mu_starve_max * max(0.0, 1.0 - intake_eff / sp%I_maint)
  else
    mu_starve = 0.0
  end if
  mu_total = mu_starve + sp%mu_background

  dN = min(ac%nindivs * (1.0 - exp(-mu_total)), ac%nindivs)
  ac%nindivs = max(0.0, ac%nindivs - dN)

  !------------------------------------------------------------------
  ! 5. Carcasses -> fast SOM, SOC(4) / SON(4)
  !------------------------------------------------------------------
  C_carcass = dN * sp%body_mass * sp%f_C_body
  N_carcass = dN * sp%body_mass * sp%f_N_body
  vegn%SOC(2) = vegn%SOC(2) + C_carcass
  vegn%SON(2) = vegn%SON(2) + N_carcass

  !------------------------------------------------------------------
  ! 6. Store daily fluxes and advance age
  !------------------------------------------------------------------
  ac%intake_plant  = intake_plant_tot
  ac%intake_prey   = intake_prey_tot
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
  ac%annualIntakePlant = ac%annualIntakePlant + intake_plant_tot
  ac%annualIntakePrey  = ac%annualIntakePrey  + intake_prey_tot
  ac%annualC_removed   = ac%annualC_removed   + C_removed
  ac%annualC_feces     = ac%annualC_feces     + C_feces
  ac%annualN_feces     = ac%annualN_feces     + N_feces
  ac%annualC_carcass   = ac%annualC_carcass   + C_carcass
  ac%annualN_carcass   = ac%annualN_carcass   + N_carcass

  end associate

end subroutine ani_cohort_daily

!=============================================================================
subroutine ani_plant_intake(idx, vegn, intake_ind, intake_tot, C_removed_tot)
  ! Type II functional response on plant leaf biomass + metabolic litter (SOC pool 1).
  ! Leaf access is height-limited: fraction = browse_height / (height + browse_height).
  ! Litter access is controlled by sp%litter_pref (0 = no litter eating).
  ! C removed from cc%bl and vegn%SOC(1) is distributed by palatability-weighted forage.
  implicit none
  integer,               intent(in)    :: idx
  type(vegn_tile_type),  intent(inout) :: vegn
  real,                  intent(out)   :: intake_ind    ! kg DM ind-1 day-1
  real,                  intent(out)   :: intake_tot    ! kg DM m-2  day-1
  real,                  intent(out)   :: C_removed_tot ! kg C  m-2  day-1 (leaves only)

  type(ani_cohort_type), pointer :: ac
  type(cohort_type),     pointer :: cc
  real :: B_cc(vegn%n_cohorts)   ! available leaf DM per plant cohort, kg DM m-2
  real :: w_cc(vegn%n_cohorts)   ! palatability-weighted forage per cohort, kg DM m-2
  real :: B_avail, w_total       ! total available DM and weighted DM, kg DM m-2
  real :: B_litter, w_litter     ! litter DM and weighted litter, kg DM m-2
  real :: C_removed_cc, C_removed_all, C_removed_litter, frac_remain
  integer :: i

  ac => vegn%ani_cohorts(idx)
  associate (sp => aftdata(ac%aft))

  B_avail = 0.0 ; w_total = 0.0 ; B_cc = 0.0 ; w_cc = 0.0

  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    ! cc%bl is kg C ind-1; convert to kg DM ind-1 via DM_to_C, then scale by density and browse fraction
    B_cc(i) = (cc%bl / sp%DM_to_C) * cc%nindivs &
            * sp%browse_height / (cc%height + sp%browse_height)
    w_cc(i)  = sp%palatability(cc%species) * B_cc(i)
    B_avail  = B_avail + B_cc(i)
    w_total  = w_total + w_cc(i)
  end do

  !-- Add metabolic litter (SOC pool 1) as supplemental forage --
  ! vegn%SOC(1) is kg C m-2; convert to kg DM m-2
  B_litter = vegn%SOC(1) / sp%DM_to_C
  w_litter = sp%litter_pref * B_litter
  B_avail  = B_avail + B_litter
  w_total  = w_total + w_litter

  !-- Type II functional response: B_avail in kg DM m-2, intake in kg DM ind-1 day-1 --
  if (B_avail > 0.0) then
    intake_ind = sp%I_max_plant * B_avail / (sp%K_half_plant + B_avail)
  else
    intake_ind = 0.0
  end if
  intake_tot    = intake_ind * ac%nindivs               ! kg DM m-2 day-1
  C_removed_all = min(intake_tot, B_avail) * sp%DM_to_C ! kg C  m-2 day-1

  !-- Distribute C removal to plant cohort leaves --
  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    if (w_total > 0.0 .and. w_cc(i) > 0.0) then
      C_removed_cc = C_removed_all * (w_cc(i) / w_total)
    else
      C_removed_cc = 0.0
    end if
    ! cap at available leaf C for this cohort (B_cc is DM; multiply back to C)
    C_removed_cc = min(C_removed_cc, B_cc(i) * sp%DM_to_C)
    if (cc%nindivs > 0.0 .and. C_removed_cc > 0.0) then
      frac_remain = max(0.0, 1.0 - C_removed_cc / max(cc%bl * cc%nindivs, 1.e-12))
      cc%bl    = cc%bl    * frac_remain
      cc%leafN = cc%leafN * frac_remain
    end if
  end do

  !-- Distribute C removal to metabolic litter --
  C_removed_litter = 0.0
  if (w_total > 0.0 .and. w_litter > 0.0) then
    C_removed_litter = C_removed_all * (w_litter / w_total)
    C_removed_litter = min(C_removed_litter, vegn%SOC(1))
    vegn%SOC(1) = max(0.0, vegn%SOC(1) - C_removed_litter)
  end if

  C_removed_tot = C_removed_all - C_removed_litter  ! leaf-only C removal returned to caller

  end associate

end subroutine ani_plant_intake

!=============================================================================
subroutine ani_prey_intake(idx, vegn, intake_ind, intake_tot)
  ! Type II functional response on total body C of all OTHER animal cohorts.
  ! Prey density is reduced proportionally to prey C removed.
  implicit none
  integer,               intent(in)    :: idx
  type(vegn_tile_type),  intent(inout) :: vegn
  real,                  intent(out)   :: intake_ind   ! kg C ind-1 day-1
  real,                  intent(out)   :: intake_tot   ! kg C m-2  day-1

  type(ani_cohort_type), pointer :: ac
  real :: B_prey(vegn%n_ani_cohorts)
  real :: B_prey_total, C_removed_prey, C_removed_cc, frac_remain
  integer :: j

  ac => vegn%ani_cohorts(idx)
  associate (sp => aftdata(ac%aft))

  B_prey_total = 0.0 ; B_prey = 0.0

  do j = 1, vegn%n_ani_cohorts
    if (j == idx) cycle
    if (vegn%ani_cohorts(j)%nindivs <= 0.0) cycle
    B_prey(j)    = vegn%ani_cohorts(j)%nindivs &
                 * aftdata(vegn%ani_cohorts(j)%aft)%body_mass &
                 * aftdata(vegn%ani_cohorts(j)%aft)%f_C_body
    B_prey_total = B_prey_total + B_prey(j)
  end do

  if (B_prey_total > 0.0) then
    intake_ind = sp%I_max_prey * B_prey_total / (sp%K_half_prey + B_prey_total)
  else
    intake_ind = 0.0
  end if
  intake_tot     = intake_ind * ac%nindivs
  C_removed_prey = min(intake_tot, B_prey_total)

  do j = 1, vegn%n_ani_cohorts
    if (j == idx) cycle
    if (vegn%ani_cohorts(j)%nindivs <= 0.0) cycle
    if (B_prey_total <= 0.0) cycle
    C_removed_cc = C_removed_prey * safe_div(B_prey(j), B_prey_total)
    C_removed_cc = min(C_removed_cc, B_prey(j))
    if (C_removed_cc > 0.0) then
      frac_remain = max(0.0, 1.0 - C_removed_cc / max(B_prey(j), 1.e-12))
      vegn%ani_cohorts(j)%nindivs = vegn%ani_cohorts(j)%nindivs * frac_remain
    end if
  end do

  end associate

end subroutine ani_prey_intake

!=============================================================================
subroutine ani_reproduction(vegn)
  ! Annual reproduction: increase nindivs for each cohort in proportion to
  ! food surplus above maintenance. No new cohorts are created; newborns join
  ! the existing cohort at adult body size (density-only increase).
  !
  ! Birth rate:
  !   f_repro = max(0, min(1, annual_intake_per_ind / (I_maint*365) - 1))
  !   births  = r_max * f_repro * nindivs
  !
  ! f_repro = 0 when intake <= I_maint (survival only, no surplus for breeding)
  ! f_repro = 1 when intake >= 2 * I_maint (full reproduction at r_max)
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  integer :: i
  real    :: annual_intake_ind, f_repro, births

  do i = 1, vegn%n_ani_cohorts
    associate (ac => vegn%ani_cohorts(i), sp => aftdata(vegn%ani_cohorts(i)%aft))

    if (sp%r_max > 0.0 .and. ac%nindivs > 0.0) then

      ! Per-individual annual intake in the same units as I_maint
      select case (sp%diet_class)
        case (AFT_HERBIVORE)
          annual_intake_ind = safe_div(ac%annualIntakePlant, ac%nindivs)
        case (AFT_CARNIVORE)
          annual_intake_ind = safe_div(ac%annualIntakePrey,  ac%nindivs)
        case default  ! omnivore
          annual_intake_ind = safe_div(sp%f_plant_diet * ac%annualIntakePlant &
                            + sp%f_prey_diet  * safe_div(ac%annualIntakePrey, sp%DM_to_C), &
                              ac%nindivs)
      end select

      ! Surplus above maintenance drives reproduction (0 at maintenance, 1 at 2x maintenance)
      f_repro = max(0.0, min(1.0, annual_intake_ind / max(sp%I_maint * days_per_year, 1.e-30) - 1.0))

      births = sp%r_max * f_repro * ac%nindivs
      ac%nindivs = ac%nindivs + births

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
