!---------------
module esdvm
 use datatypes
 use io_mod
 use soil_mod
 implicit none
 private

 !Core functions
 public :: initialize_vegn_tile
 public :: vegn_phenology
 public :: vegn_CNW_budget_fast,vegn_daily_update, vegn_demographics
 public :: vegn_SingleCohort_annual_update
 public :: relayer_cohorts, vegn_hydraulic_states, vegn_SW2HW_hydro
 public :: kill_lowdensity_cohorts,kill_old_grass,vegn_mergecohorts

 !For specific experiments
 public :: vegn_fire, vegn_species_recovery, vegn_species_switch
 public :: vegn_annualLAImax_update, vegn_gap_fraction_update
 public :: reset_vegn_initial

 integer :: MaxCohortID = 0

 contains

!========================================================================
!==================== BiomeE surbroutines ===============================
!========================================================================
subroutine vegn_CNW_budget_fast(vegn, forcing)
  ! hourly carbon, nitrogen, and water dynamics, Weng 2016-11-25
  ! include Nitrogen uptake and carbon budget
  ! C_growth is calculated here to drive plant growth and reproduciton
  type(vegn_tile_type), intent(inout) :: vegn
  type(climate_data_type),intent(in):: forcing

  !-------local var
  type(cohort_type), pointer :: cc  ! current cohort
  integer :: i
  real :: tair, tsoil ! temperature of soil, degC
  real :: thetaS ! soil wetness, unitless
  real :: NSC_supply,LR_demand,LR_deficit
  real :: LeafGrowthMin, RootGrowthMin,NSCtarget,v
  real :: LR_growth,WS_growth
  real :: R_days,fNSC,fLFR,fStem
  integer :: layer

  ! Climatic variable
  tair   = forcing%Tair  - 273.16   ! degC
  tsoil  = forcing%tsoil - 273.16  ! degC
  thetaS = (vegn%wcl(2)-vegn%WILTPT)/(vegn%FLDCAP-vegn%WILTPT)
  vegn%annualPET = vegn%annualPET + potentialET(forcing) * step_seconds ! Potential ET, kg m-2 step-1

  ! Water supply for leaves
  ! Soil water parameters (psi and conductivity for each layer)
  call SoilWater_psi_K(vegn)

#ifdef Hydro_test
  ! Update plant hydraulic status, fluxes, and water supply for transpiration
  call Plant_water_dynamics_linear(vegn)
  ! Photosynsthesis
  call vegn_photosynthesis(forcing, vegn)

  ! Phloem transport, Mazen Nakad, 10/08/2023
  call vegn_Phloem_transport(forcing,vegn)
#else
  ! Water supply from soil directly
  call SoilWaterSupply(vegn)
  call vegn_photosynthesis(forcing, vegn)
  call SoilWaterTranspUpdate(vegn)
#endif

  ! Soil water dynamics: infiltration and surface evap.
  call SoilWaterDynamics(forcing,vegn)

  ! Plant Respiration
  call vegn_respiration(forcing,vegn)

  ! Soil organic matter decomposition
  call Soil_BGC(vegn, forcing%tsoil, thetaS)

  !! Nitrogen uptake
  call vegn_N_uptake(vegn, forcing%tsoil)

end subroutine vegn_CNW_budget_fast

!==========================================================================
subroutine vegn_daily_update(vegn, deltat)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! dt_daily_yr, 1.0/365.0

  call vegn_age(vegn,deltat) ! Update vegn age
  call vegn_phenology(vegn)
  call vegn_growth(vegn)
  ! Turnover of leaves and fine roots
  call vegn_tissue_turnover(vegn)
  call SoilWater_psi_K(vegn)
  call grass_thinning(vegn) ! Added 12/03/2025
  !call vegn_daily_starvation(vegn)
  call vegn_cohort_update(vegn)
  call Vegn_N_deposition(vegn,deltat) ! Daily N deposition
  call vegn_sum_tile(vegn)  ! Update tile variables
end subroutine vegn_daily_update

!==========================================================================
subroutine vegn_demographics(vegn, deltat)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! seconds of a year
  !-------- local vars ----------
  !real totN0

  !! Total N balance checking
  !call vegn_sum_tile(vegn)
  !totN0 = TotalN(vegn)
  ! For the incoming year
  call vegn_annual_starvation(vegn) ! turn it off for grass run
  call vegn_nat_mortality(vegn, deltat)
  call vegn_reproduction(vegn)
  !call check_N_conservation(vegn,totalN1,'annual')

end subroutine vegn_demographics

!========================================================================
subroutine vegn_SingleCohort_annual_update(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  call vegn_reprod_samesized(vegn)
  call vegn_Wood_turnover(vegn)
  call vegn_selfthinning(vegn)
  call vegn_cohort_update(vegn)
end subroutine vegn_SingleCohort_annual_update

!================ For single cohort test 02/21/2024 =====================
!========================================================================
subroutine vegn_reprod_samesized(vegn)
  ! This subroutine only generates the exact same individuals as parents
  ! Annual time step. (02/21/2024)
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  !-----local--------
  type(cohort_type), pointer :: cc ! parent and child cohort pointers
  integer :: i
  real :: plantC, plantN
  real :: n_new, n_newC, n_newN
  real :: N_demand, N_left, C_demand, C_left

  !--------- Reproduction ---------------
  !looping through cohorts and convert their C_seed to plants
  !by increasing cc%n
  ! Looping through all cohorts
  do i=1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     !Carbon content of a current individual of this cohort
     plantC = cc%bl + cc%br + cc%bsw + cc%bHW + cc%nsc
     plantN = cc%leafN+cc%rootN+cc%swN+cc%hwN+cc%NSN
     n_newC  = cc%seedC * cc%nindivs / plantC
     n_newN  = cc%seedN * cc%nindivs / plantN
     n_new = min(n_newC, n_newN)

     N_demand = n_new * plantN
     C_demand = n_new * plantC
     N_left = cc%seedN * cc%nindivs - N_demand
     C_left = cc%seedC * cc%nindivs - C_demand

     ! Update density and seed pools
     if(n_new > 0.0)then
        cc%nindivs  = cc%nindivs + n_new
        cc%NSN = cc%NSN + N_left/cc%nindivs ! put the left N back to NSN pool
        cc%NSC = cc%NSC + C_left/cc%nindivs
        cc%seedC = 0.0
        cc%seedN = 0.0
     endif
  enddo
end subroutine vegn_reprod_samesized

!========================================================================
subroutine vegn_Wood_turnover(vegn)
  ! Mortality for the single cohort test, 02/21/2024
  type(vegn_tile_type), intent(inout) :: vegn

  !-------local var----------
  type(cohort_type), pointer :: cc    ! current cohort
  real :: loss_coarse, lossN_coarse, loss_fine, lossN_fine
  real :: alpha_WD ! turnover rate, fraction per day
  real :: dCSW, dCHW, dNHW, dNSW  ! Turnover of tissues
  real :: dAleaf, dBL,dBR,dNL,dNR
  integer :: i

  ! update plant carbon and nitrogen for all cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )

     ! Wood turnover
     alpha_WD = mortality_rate(cc)
     dCSW = cc%bsw * alpha_WD
     dNSW = cc%swN * alpha_WD
     dCHW = cc%bHW * alpha_WD
     dNHW = cc%hwN * alpha_WD

     ! Update plant C and N pools
     cc%bsw  = cc%bsw - dCSW
     cc%bHW  = cc%bHW - dCHW
     cc%swN  = cc%swN - dNSW
     cc%hwN  = cc%hwN - dNHW

     !! update plant architecture given increase of bsw
     call BM2Architecture(cc, cc%bsw+cc%bHW)

     ! Update bl_max, br_max, and NSNmax with shifts from understory to the top layer
     call update_max_LFR_NSN(cc)

     ! Update leaves and fine roots
     dBL = max(0.0, cc%bl - cc%bl_max)
     dBR = max(0.0, cc%br - cc%br_max)
     dNL = dBL/sp%CNleaf0
     dNR = dBR/sp%CNroot0

     cc%bl    = cc%bl    - dBL
     cc%br    = cc%br    - dBR
     cc%leafN = cc%leafN - dNL
     cc%rootN = cc%rootN - dNR

     ! Update leaf area and LAI
     cc%Aleaf = BL2Aleaf(cc%bl,cc)

     ! put C and N into soil pools
     dAleaf = BL2Aleaf(dBL,cc)
     loss_coarse  = (1.-l_fract) * cc%nindivs * (dBL - dAleaf * LMAmin + dCSW + dCHW)
     loss_fine    = (1.-l_fract) * cc%nindivs * (dBR + dAleaf * LMAmin)
     lossN_coarse = (1.-retransN)* cc%nindivs * (dNL - dAleaf * sp%LNbase + dNSW + dNHW)
     lossN_fine   = (1.-retransN)* cc%nindivs * (dNR + dAleaf * sp%LNbase)

     vegn%SOC(1) = vegn%SOC(1)   +  &
                   fsc_fine * loss_fine + fsc_wood * loss_coarse
     vegn%SOC(2) = vegn%SOC(2) +  &
                   ((1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse)
     vegn%SON(1)  = vegn%SON(1) +    &
                   fsc_fine * lossN_fine + fsc_wood * lossN_coarse
     vegn%SON(2) = vegn%SON(2) + &
                   (1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse

     !    annual N from plants to soil
     vegn%NorgP2S = vegn%NorgP2S + lossN_fine + lossN_coarse

#ifdef Hydro_test
     ! Assume water in plants goes to first layer of soil
     vegn%wcl(1) = vegn%wcl(1) + cc%nindivs * &
                   (cc%W_leaf * dBL /(cc%bl +dBL) + &
                    cc%W_stem * dCSW/(cc%bsw+dCSW)+ &
                    cc%W_dead * dCHW/(cc%bhw+dCHW))/(thksl(1)*1000.0)
#endif

    END ASSOCIATE
  enddo
end subroutine vegn_Wood_turnover

!========================================================================
subroutine vegn_selfthinning(vegn)
  ! Selfthinning removes the individuals that cannot be hold.
  ! Daily time step. 02/21/2024
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  !-----local--------
  type(cohort_type), pointer :: cc ! parent and child cohort pointers
  real  :: totalCA, CAcrown, dn
  integer :: L = 1

  totalCA  = 0.d0
  do i=1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate (sp => spdata(cc%species) )
       CAcrown  = (1.0 + sp%f_cGap) * cc%Acrown
       totalCA  = totalCA + cc%nindivs * CAcrown
       cc%layer = L
       if(totalCA .gt. L) then
           ! selfthinning
           dn = MIN(cc%nindivs,(totalCA - L) / CAcrown)
           ! Record the actual mortality rate for output
           cc%mu = mortality_rate(cc) + dn/cc%nindivs
           ! Kill some trees in this cohort to let totalCA = L
           cc%nindivs = cc%nindivs - dn
           totalCA = totalCA - dn * CAcrown
           ! Put the dead trees carbon to soil pools, pp%Tpool (1,12,1)
           call plant2soil(vegn,cc,dn)
       endif
     end associate
  enddo
end subroutine vegn_selfthinning
!========================================================================

!=============== Plant physiology =======================================
!=============== Hourly subroutines =====================================
! Weng 2017-10-18:compute stomatal conductance, photosynthesis and respiration
! updates cc%An_op and cc%An_cl, from LM3
subroutine vegn_photosynthesis (forcing, vegn)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  !use datatypes
  implicit none

  type(climate_data_type), intent(in)    :: forcing
  type(vegn_tile_type), intent(inout)    :: vegn

  !----- local var --------------
  type(cohort_type), pointer :: cc
  real :: rad_top      ! downward PAR radiation at the top of the canopy, W/m2
  real :: rad_net      ! net PAR radiation absorbed by the canopy, W/m2
  real :: Tair, TairK  ! air temperature, degC and degK
  real :: cana_q       ! specific humidity in canopy air space, kg/kg
  real :: cana_co2     ! co2 concentration in canopy air space, mol CO2/mol dry air
  real :: p_surf       ! surface pressure, Pa
  real :: water_supply ! water supply per m2 of leaves (mol H2O m-2 leaf s-1)
  real :: fc           ! wet and snow-covered fraction of leaves
  real :: psyn         ! net photosynthesis, mol C/(m2 leaf s)
  real :: resp         ! leaf respiration, mol C/(m2 leaf s)
  real :: wd, transp   ! mol H2O per m2 of leaf per second
  real :: wd1          ! Kg H2O per tree per time step
  real :: kappa        ! light extinction coefficient of crown layers
  real, dimension(10) :: f_light
  integer :: i, layer, nlayers, ncoh

  ncoh = vegn%n_cohorts
  ! Determine number of light layers to use for f_light calculation.
  nlayers = min(int(vegn%CAI + 1.0), 9)

  !! Light supply for photosynthesis
  ! Initialize kp (ensure vegn%kp has appropriate bounds)
  vegn%kp = 0.0
  do i = 1, ncoh
     cc => vegn%cohorts(i)
     layer = max(1, min(cc%layer, 9))
     ! accumulate per-layer projected crown area * extinction
     vegn%kp(layer) = vegn%kp(layer) + cc%extinct * cc%Acrown * cc%nindivs
  end do

  ! Light fraction -- compute up to nlayers
  f_light = 0.0
  f_light(1) = 1.0
  do i = 2, nlayers
      f_light(i) = f_light(i-1) * (exp(-vegn%kp(i-1) * vegn%LAI_L(i-1)) + vegn%f_gap(i-1))
  end do

  ! Photosynthesis loop across cohorts
  do i = 1, ncoh
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
       if (cc%status == LEAF_ON .and. cc%Aleaf > tiny) then
          ! Convert forcing data / units
          layer = max(1, min(cc%layer, nlayers))
          ! ensure forcing%radiation is PAR or convert accordingly
          rad_top = f_light(layer) * f_PAR * forcing%radiation
          rad_net = f_light(layer) * f_PAR * forcing%radiation * 0.9

          p_surf = forcing%P_air          ! Pa
          TairK  = forcing%Tair          ! K
          Tair   = forcing%Tair - 273.16 ! degC (esat expects degC)
          cana_q = (esat(Tair) * forcing%RH * mol_h2o) / (p_surf * mol_air) ! specific humidity
          cana_co2 = forcing%CO2 * 1.0e-6  ! ppm -> mol/mol
          water_supply = cc%W_supply / (cc%Aleaf * step_seconds * mol_h2o) ! water supply: mol H2O m-2 leaf s-1

          fc = 0.0 ! Assume it is zero
          call gs_Leuning(rad_top, rad_net, TairK, cana_q, cc%LAI, &
                          p_surf, water_supply, cc%species, sp%pt, &
                          cana_co2, cc%extinct, fc, cc%layer,      &
                          psyn, resp, wd, transp)

          ! store outputs (document sign conventions)
          cc%An_op = psyn          ! mol C s-1 m-2 leaf
          cc%An_cl = -resp         ! (kept negative here as original — verify convention)
          cc%gpp   = (psyn - resp) * mol_C * cc%Aleaf * step_seconds
          cc%transp = transp * mol_h2o * cc%Aleaf * step_seconds
          if (abs(wd) > tiny) then
             cc%w_scale = transp / wd
          else
             cc%w_scale = -9999.0
          endif

          ! For UFL drought mortality (accumulate potential demand)
          wd1 = wd * mol_h2o * cc%Aleaf * step_seconds
          cc%totDemand = cc%totDemand + wd1
          cc%dailyWdmd = cc%dailyWdmd + wd1
       else
          ! no leaves => zero photosynthesis & transp
          cc%An_op   = 0.0
          cc%An_cl   = 0.0
          cc%gpp     = 0.0
          cc%transp  = 0.0
          cc%w_scale = -9999.0
       endif

       ! NaN checks using ieee_is_nan for portability
       if (ieee_is_nan(cc%gpp)) then
          write(*,*) 'Error: cc%gpp is NaN for cohort ', i, ' species ', cc%species
          stop 1
       endif
       if (ieee_is_nan(cc%transp)) then
          write(*,*) 'Error: transp is NaN, wd, transp, lai = ', wd, transp, cc%LAI
          stop 1
       endif
     end associate
  end do  ! cohorts loop

end subroutine vegn_photosynthesis

!============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, &
                   p_surf, ws, pft, pt, ca, kappa, f_w, layer, &
                   apot, acl, wd, transp)
  ! Note: included Yang Qi's corrections on:
  ! lai_eq, Ag_l, and gs_w (02/11/2026)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  ! ---- arguments
  real,    intent(in)   :: rad_top   ! PAR dn on top of the canopy, w/m2
  real,    intent(in)   :: rad_net   ! PAR net on top of the canopy, w/m2
  real,    intent(in)   :: tl        ! leaf temperature, degK
  real,    intent(in)   :: ea        ! specific humidity in the canopy air, kg/kg
  real,    intent(in)   :: lai       ! leaf area index
  real,    intent(in)   :: p_surf    ! surface pressure, Pa
  real,    intent(in)   :: ws        ! water supply, mol H20/(m2 of leaf s)
  integer, intent(in)   :: pft       ! species
  integer, intent(in)   :: pt        ! physiology type (C3 or C4)
  real,    intent(in)   :: ca        ! CO2 in canopy air, mol/mol dry air
  real,    intent(in)   :: kappa     ! canopy extinction coefficient
  real,    intent(in)   :: f_w       ! fraction of leaf wet/snow-covered
  integer, intent(in)   :: layer     ! canopy layer index
  real,    intent(out)  :: apot      ! net photosynthesis, mol C/(m2 leaf s)
  real,    intent(out)  :: acl       ! leaf respiration, mol C/(m2 leaf s)
  real,    intent(out)  :: wd, transp ! mol H2O/(m2 leaf s)

  ! ---- local vars
  ! photosynthesis
  real :: vm
  real :: kc, ko          ! Michaelis-Menten constants for CO2 and O2
  real :: ci
  real :: capgam          ! CO2 compensation point
  real :: f2, f3
  real :: coef0, coef1
  real :: Resp

  ! conductance related
  real :: gs
  real :: b
  real :: ds              ! humidity deficit, kg/kg
  real :: hl              ! saturated specific humidity at leaf temp, kg/kg
  real :: do1

  ! miscellaneous / numerics
  real :: dum2
  real, parameter :: light_crit = 0.0
  real, parameter :: gs_lim     = 0.25
  real, parameter :: lai_min    = 1.0e-6
  real, parameter :: kappa_min  = 1.0e-6

  ! new average computations
  real :: lai_eq
  real :: light_top
  real :: par_net
  real :: Ag
  real :: An
  real :: Ag_l
  real :: Ag_rb
  real :: anbar
  real :: gsbar
  real :: w_scale

  ! soil water stress
  real :: Ed, an_w, gs_w

  ! protected/clamped values
  real :: kappa_eff, f_w_eff, arg

  b   = 0.01
  do1 = 0.09 ! kg/kg

  ! Clamp wet fraction and kappa (extra safety)
  f_w_eff   = min(1.0, max(0.0, f_w))
  kappa_eff = max(kappa, kappa_min)

  ! Convert Solar influx from W/(m^2) to mol_of_quanta/(m^2 s) PAR
  light_top = rad_top * rad_phot ! for lai_eq
  par_net   = rad_net * rad_phot ! for leaf photosynthesis rate

  ! Humidity deficit, kg/kg
  call qscomp(tl, p_surf, hl)
  ds = max(hl - ea, 0.0)

  ! Guard LAI to avoid division by zero
  if (lai <= lai_min) then
     apot   = 0.0
     acl    = 0.0
     wd     = 0.0
     transp = 0.0
     return
  endif

  associate ( sp => spdata(pft) )
    !  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf
    !  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf
    !  vm=sp%Vmax*exp(3000.0*(1.0/288.2-1.0/tl))

    ! Weng, 2013-01-10
    ko = 0.248    * exp(35948.0/Rgas*(1.0/298.2-1.0/tl)) * p_sea/p_surf
    kc = 0.000404 * exp(59356.0/Rgas*(1.0/298.2-1.0/tl)) * p_sea/p_surf
    vm = sp%Vmax  * exp(24920.0/Rgas*(1.0/298.2-1.0/tl))

    !decrease Vmax due to aging of temperate deciduous leaves
    !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
    !! Turned off by Weng, 2013-02-01, since we can't trace new leaves
    !  if (sp%leafage_tau>0 .and. leafage>sp%leafage_onset) then
    !     vm=vm*exp(-(leafage-sp%leafage_onset)/sp%leafage_tau)
    !  endif

    ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986
    capgam = 0.5 * kc / ko * 0.21 * 0.209 ! Farquhar & Caemmerer 1982

    ! Find respiration for the whole canopy layer
    !  Resp=sp%gamma_resp*vm*lai /((layer-1)*1.0+1.0)  ! Weng, 2013-01-17 add '/ ((layer-1)*1.0+1.0)'

    ! 2014-09-03, for Nitrogen model: resp = D*(A + B*LMA)
    ! (A+B*LMA) = LNA, D=Vmax/LNA = 25E-6/0.0012 = 0.02 for a standard deciduous species
    !! Leaf resp as a function of nitrogen
    !  Resp=sp%gamma_resp*0.04*sp%LNA  & ! basal rate, mol m-2 s-1
    !       * exp(24920/Rgas*(1.0/298.2-1.0/tl))         & ! temperature scaled
    !       * lai                                        & ! whole canopy
    !       /((layer-1)*1.0+1.0)                         !
    !! as a function of LMA
    !  Resp=(sp%gamma_LNbase*sp%LNbase+sp%gamma_LMA*sp%LMA)  & ! basal rate, mol m-2 s-1
    !  Resp=sp%gamma_LNbase*(2.5*sp%LNA-1.5*sp%LNbase)     & ! basal rate, mol m-2 s-1

    ! Respiration for whole canopy (scales with LAI)
    Resp = sp%gamma_LN/seconds_per_year          & ! per seconds,  mol m-2 s-1
         * sp%LNA * lai / mol_c                 & ! whole canopy, mol m-2 s-1
         * exp(24920.0/Rgas*(1.0/298.2-1.0/tl))   ! temperature scaled

    ! Temperature effects (suppression at cold/hot)
    Resp = Resp / ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                 * (1.0+exp(0.4*(tl-45.0-TFREEZE))))

    ! initialize
    Ag_l  = 0.0
    Ag_rb = 0.0
    Ag    = 0.0
    An    = -Resp
    anbar = -Resp / lai
    gsbar = b

    ! find the LAI level at which gross photosynthesis rates are equal
    if (light_top > light_crit) then

       if (pt == PT_C4) then
          coef0 = (1.0 + ds/do1) / sp%m_cond
          ci    = (ca + 1.6*coef0*capgam) / (1.0 + 1.6*coef0)

          if (ci > capgam) then
             f2   = vm
             f3   = 18000.0 * vm * ci
             dum2 = min(f2, f3)

             ! lai_eq = -log(dum2/(sp%alpha_ps*light_top))/kappa
             arg = dum2 / max(sp%alpha_ps*light_top, tiny)
             if (arg > tiny) then
                lai_eq = -log(arg) / kappa_eff
             else
                lai_eq = 0.0
             endif
             lai_eq = min(max(0.0, lai_eq), lai)

             Ag_l  = sp%alpha_ps * par_net * (exp(-lai_eq*kappa_eff) - exp(-lai*kappa_eff)) / kappa_eff
             Ag_rb = dum2 * lai_eq

             Ag = (Ag_l + Ag_rb) / ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                                  * (1.0+exp(0.4*(tl-45.0-TFREEZE))))

             An    = Ag - Resp
             anbar = An / lai

             if (anbar > 0.0) then
                gsbar = anbar / (ci - capgam) / coef0
             endif
          endif

       else  ! C3
          coef0 = (1.0 + ds/do1) / sp%m_cond
          coef1 = kc * (1.0 + 0.209/ko)
          ci    = (ca + 1.6*coef0*capgam) / (1.0 + 1.6*coef0)

          f2   = vm * (ci - capgam) / (ci + coef1)
          f3   = vm / 2.0
          dum2 = min(f2, f3)

          if (ci > capgam) then
             ! lai_eq=-log(dum2*(ci+2*capgam)/(ci-capgam)/(sp%alpha_ps*light_top))/kappa
             arg = dum2 * (ci + 2.0*capgam) / max((ci - capgam) * sp%alpha_ps * light_top, tiny)
             !(sp%alpha_ps*light_top*kappa) --> (sp%alpha_ps*light_top) ! Yang Qi: no kappa)
             if (arg > tiny) then
                lai_eq = -log(arg) / kappa_eff
             else
                lai_eq = 0.0
             endif
             lai_eq = min(max(0.0, lai_eq), lai)

             Ag_l  = sp%alpha_ps * (ci-capgam)/(ci+2.0*capgam) * par_net  &
                   * (exp(-lai_eq*kappa_eff) - exp(-lai*kappa_eff)) / kappa_eff
                   ! Yang Qi: 1.0-exp(-lai*kappa)) --> kappa (kappa_eff)
             Ag_rb = dum2 * lai_eq

             Ag = (Ag_l + Ag_rb) / ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                                  * (1.0+exp(0.4*(tl-45.0-TFREEZE))))

             An    = Ag - Resp
             anbar = An / lai

             if (anbar > 0.0) then
                gsbar = anbar / (ci - capgam) / coef0
             endif
          endif
       endif

    endif ! light_top > light_crit

    ! Wet/snow-covered leaf scaling
    an_w = anbar * (1.0 - sp%ps_wet * f_w_eff)
    gs_w = gsbar * (1.0 - sp%ps_wet * f_w_eff) ! ! Yang Qi's correction (09/24, no 1.56)
    if (gs_w > gs_lim) then
       if (an_w > 0.0) an_w = an_w * gs_lim / gs_w
       gs_w = gs_lim
    endif
  end associate

  ! Water availability diagnostic demand
  Ed = gs_w * ds * mol_air / mol_h2o

  if (Ed > ws) then
     w_scale = ws / Ed
     gs_w = w_scale * gs_w
     if (an_w > 0.0) an_w = an_w * w_scale
     if (an_w < 0.0 .and. gs_w > b) gs_w = b
  endif

  apot   = an_w
  acl    = -Resp / lai
  transp = min(ws, Ed)   ! mol H2O/(m2 leaf s)
  wd     = Ed

  ! Convert stomatal conductance to m/s (kept local; not returned)
  gs = gs_w
  gs = gs * Rgas * Tl / p_surf

  ! Error check
  if (ieee_is_nan(transp)) then
     write(*,*) 'ws, Ed = ', ws, Ed
     stop '"transp" is a NaN'
  endif
end subroutine gs_Leuning

! ============================================================================
! Copied from LM4-ESM4_2021.02: vegn_radiation.F90
subroutine twostream(mu,mu_bar,LAI,albedo_g,phi1,phi2,rl,tl, &
                     transm_dir, scatter_dir, albedo_dir,    &
                     transm_dif, albedo_dif )

  real, intent(in)  :: mu         ! cosine of direct light zenith angle
  real, intent(in)  :: mu_bar     ! average inverse diffuse optical depth per unit leaf area
  real, intent(in)  :: LAI        ! leaf area index
  real, intent(in)  :: albedo_g   ! ground surface albedo
  real, intent(in)  :: phi1, phi2 ! coefficients of expression for G_mu
  real, intent(in)  :: rl         ! reflectivity of leaves
  real, intent(in)  :: tl         ! transmittance of leaves
  ! output
  real, intent(out) :: transm_dir ! canopy transmittance for direct beam -- that
                                  ! is, the part of the beam that passes through
                                  ! the canopy untouched
  real, intent(out) :: scatter_dir! part of direct beam scattered down, at the
                                  ! bottom of the canopy
  real, intent(out) :: albedo_dir ! overall land surface albedo for direct beam
  real, intent(out) :: transm_dif ! canopy transmittance for diffuse incident light
  real, intent(out) :: albedo_dif ! overall land surface albedo for diffuse incident light

  ! ---- local vars
  real :: G_mu        ! relative projected leaf area in direction of direct beam
  real :: K           ! optical depth for direct beam per unit LAI
  real :: g1,g2,g3,g4 ! coefficients in the two-stream equation
  real :: kappa       ! eigenvalue of free solution
  real :: a_up, b_up, c_up ! coefficients of upward diffuse light flux
  real :: a_dn, b_dn, c_dn ! coefficients of downward diffuse light flux
  real :: x1,x2       ! intermediate coefficients
  real :: a11,a12,a21,a22, d1,d2 ! coefficients of linear system
  real :: D           ! determinant of the matrix
  real :: A,B         ! coefficients of diffuse light function
  real :: dif_dn_bot, dif_up_bot, dif_up_top

  real, parameter :: eta = 6.0  ! this value is suitable only for uniform leaf
                                ! angular distribution !!!

  ! calculate coefficients of optical path
  G_mu = phi1 + phi2 * mu
  K = G_mu/mu

  ! given optical parameters, calculate coefficients of basic equation
  g1 = (1-(rl+tl)/2+(rl-tl)/eta)/mu_bar;
  g2 = (  (rl+tl)/2+(rl-tl)/eta)/mu_bar;
  g3 = G_mu*((rl+tl)/2+mu*(rl-tl)/eta/G_mu);
  g4 = G_mu*((rl+tl)/2-mu*(rl-tl)/eta/G_mu);

  ! calculate eigenvalue of free solution (=exponent coefficient of
  ! free solution, notes 12)
  kappa = sqrt(g1**2-g2**2);

  ! calculate forced term coefficients for diffuse light intensity
  c_up = ( K*g3-g1*g3-g2*g4)/(K*K-g1*g1+g2*g2);
  c_dn = (-K*g4-g1*g4-g2*g3)/(K*K-g1*g1+g2*g2);
  ! calculate intermediate coefficients for solution
  x1 = g1+g2+kappa; x2 = g1+g2-kappa;

  !write(*,*)mu,K,g1,g2,g3,g4,c_up,c_dn

  ! calculate coefficients of the matrix
  a11 = x2;
  a12 = x1;
  d1  = -c_dn;
  a21 = exp(kappa*LAI)*(x1-albedo_g*x2);
  a22 = exp(-kappa*LAI)*(x2-albedo_g*x1);
  d2  = exp(-K*LAI)*(albedo_g*c_dn + albedo_g*mu - c_up);
  ! solve the equation system
  D = a11*a22-a12*a21;
  A = (d1*a22-d2*a12)/D;
  B = (a11*d2-a21*d1)/D;

  ! calculate coefficients of the diffuse light solution
  a_up = A*x1; b_up=B*x2;
  a_dn = A*x2; b_dn=B*x1;

  ! calculate downward diffuse light at the bottom of the canopy
  dif_dn_bot = a_dn*exp(kappa*LAI) + b_dn*exp(-kappa*LAI)+c_dn*exp(-K*LAI);
  dif_up_bot = a_up*exp(kappa*LAI) + b_up*exp(-kappa*LAI)+c_up*exp(-K*LAI);

  ! calculate canopy transmittance and scattered part for direct light
  scatter_dir = dif_dn_bot/mu; ! scatter
  transm_dir  = exp(-K*LAI); ! transmittance for unmolested direct light

  ! calculate canopy reflectance for direct light
  dif_up_top = a_up + b_up + c_up;
  albedo_dir = dif_up_top/mu;

  ! calculate upward diffuse light at the top of the canopy
  ! no need to recalculate D, since the matrix is the same
  d1 = 1.0;
  d2 = 0.0;
  A = (d1*a22-d2*a12)/D;
  B = (a11*d2-a21*d1)/D;

  ! calculate coefficients of the diffuse light solution
  a_up = A*x1; b_up=B*x2;
  a_dn = A*x2; b_dn=B*x1;

  transm_dif = a_dn*exp(kappa*LAI) + b_dn*exp(-kappa*LAI);
  albedo_dif = a_up + b_up;

end subroutine twostream

!============================================================================
subroutine vegn_respiration(forcing,vegn)
  type(climate_data_type),intent(in)  :: forcing
  type(vegn_tile_type), intent(inout) :: vegn

  !---------local var ---------
  type(cohort_type), pointer :: cc
  real :: TairK     ! air temperature, degK
  real :: tf,tfs ! thermal inhibition factors for above- and below-ground biomass
  real :: r_leaf, r_stem, r_root
  real :: Acambium  ! cambium area, m2/tree
  real :: fnsc,NSCtarget ! used to regulation respiration rate
  real :: r_Nfix    ! respiration cost for N fixation
  real :: facuC     ! Carbon available for faculative N fixation
  integer :: i

  !-----------------------
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     TairK = forcing%tair
     associate ( sp => spdata(cc%species) )
       ! Maintenance respiration
       !call plant_respiration(cc,forcing%tair) ! get resp per tree per time step
       ! temperature response function
       tf  = exp(9000.0*(1.0/298.16-1.0/tairK))

       !  tfs = thermal_inhibition(tsoil)  ! original
       tfs = tf ! Rm_T_response_function(tsoil) ! Weng 2014-01-14
       ! With nitrogen model, leaf respiration is a function of leaf nitrogen
       NSCtarget = 3. * (cc%bl_max + cc%br_max)
       fnsc = min(max(0.0,cc%nsc/NSCtarget-0.05),1.0)
       Acambium = PI * cc%DBH * cc%height * 1.2

       ! r_leaf = fnsc*sp%gamma_LN  * cc%leafN * tf * dt_fast_yr  ! tree-1 step-1
       ! LeafN  = sp%LNA * cc%Aleaf
       r_leaf   = cc%An_cl * mol_C * cc%Aleaf * step_seconds
       r_stem   = fnsc*sp%gamma_SW  * Acambium * tf * dt_fast_yr ! kgC tree-1 step-1
       r_root   = fnsc*sp%gamma_FR  * cc%rootN * tf * dt_fast_yr ! root respiration ~ root N

       if(sp%R0_Nfix > 1.0e-6)then ! Nitrogen fixer
          ! Baseline nitrogen fixation (Obligate)
          cc%fixedN = sp%R0_Nfix * cc%br * fnsc * tf * dt_fast_yr ! kgN tree-1 step-1
          r_Nfix    = sp%C0_Nfix * cc%fixedN ! KgC tree-1 step-1

          ! Extra C for N fixation (Facultative)
          facuC = sp%S_facuN * cc%extraC / steps_per_day ! Carbon used for faculative N fixation
          cc%fixedN = cc%fixedN + facuC / sp%C0_Nfix !
          r_Nfix    = r_Nfix    + facuC
          cc%nsc    = cc%nsc    - facuC
       else
          cc%fixedN = 0.0
          r_Nfix    = 0.0
       endif

       ! Total Respiration and NPP
       cc%resl = r_leaf + r_stem ! tree-1 step-1
       cc%resr = r_root + r_Nfix ! tree-1 step-1
       cc%resp = cc%resl + cc%resr + cc%resg/steps_per_day   !kgC tree-1 step-1
       cc%resp = min(cc%resp,max(0.0,cc%nsc+cc%gpp-cc%resp)) ! Hack, Weng 09/24/2023
       cc%npp  = cc%gpp  - cc%resp ! kgC tree-1 step-1

       ! Update NSC and NSN
       cc%nsc = cc%nsc + cc%npp
       cc%NSN = cc%NSN + cc%fixedN
     end associate
  enddo ! all cohorts

end subroutine vegn_respiration

!=====================================================
! Weng, 2016-11-28
subroutine vegn_N_uptake(vegn, tsoil)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil, degK

  !-------local var
  type(cohort_type),pointer :: cc
  real    :: totNup    ! kgN m-2
  real    :: avgNup
  real    :: rho_N_up,N_roots   ! actual N uptake rate
  logical :: NSN_not_full
  integer :: i

  !! Nitrogen uptake parameter
  ! It considers competition here. How much N one can absorp depends on
  ! how many roots it has and how many roots other individuals have.
  N_Roots  = 0.0
  vegn%N_uptake = 0.0
  if(vegn%mineralN > 0.0)then
     do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        associate (sp => spdata(cc%species))
        !A scheme for deciduous to get enough N:
        !cc%NSNmax = sp%fNSNmax*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
        if(cc%NSN < cc%NSNmax) &
          N_Roots = N_Roots + cc%br * cc%nindivs

        end associate
     enddo
     ! M-M equation for Nitrogen absoption, McMurtrie et al. 2012, Ecology & Evolution
     ! rate at given root biomass and period of time
     if(N_roots>0.0)then
        ! Add a temperature response equation herefor rho_N_up0 (Zhu Qing 2016)
        ! rho_N_up = 1.-exp(-rho_N_up0 * N_roots/(N_roots0+N_roots) * hours_per_year * dt_fast_yr) ! rate at given root density and time period
        rho_N_up = rho_N_up0 * N_roots/(N_roots0+N_roots) * hours_per_year * dt_fast_yr
        totNup = rho_N_up * vegn%mineralN * exp(9000.0 * (1./298.16 - 1./tsoil)) ! kgN m-2 time step-1
        avgNup = totNup / N_roots ! kgN time step-1 kg roots-1
        ! Nitrogen uptaken by each cohort, N_uptake
        vegn%N_uptake = 0.0
        do i = 1, vegn%n_cohorts
           cc => vegn%cohorts(i)
           cc%N_uptake  = 0.0
           if(cc%NSN < cc%NSNmax)then
               cc%N_uptake = min(cc%br*avgNup, cc%NSNmax- cc%NSN)
               cc%nsn      = cc%nsn   + cc%N_uptake
               cc%NupYr    = cc%NupYr + cc%N_uptake
               ! subtract N from mineral N
               vegn%mineralN = vegn%mineralN - cc%N_uptake * cc%nindivs
               vegn%N_uptake = vegn%N_uptake + cc%N_uptake * cc%nindivs
           endif
        enddo
        cc =>null()
     endif ! N_roots>0
  endif
end subroutine vegn_N_uptake

!============================================================================
!========================== Daily subroutines ===============================
subroutine vegn_growth(vegn)
  ! updates cohort biomass pools, LAI, and height using accumulated
  ! C_growth and bHW_gain
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc    ! current cohort
  real :: LFR_deficit, LF_deficit, FR_deficit
  real :: G_LFR  ! amount of carbon spent on leaf and root growth
  real :: Cgrowth, Nsupply
  real :: dBL, dBR, dBSW, dSeed ! growth of leaf, root, sapwood, and seeds, kgC/individual
  real :: DBH0 ! the DBH before grwoth
  real :: Ndemand, extraN, N_used
  real :: r_N_SD ! Nitrogen Supply/Demand ratio
  integer :: i,j,k

  !Allocate C_gain to tissues
  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    ! Skip non-growing season
    if(cc%status == LEAF_OFF)then
      cc%resg = 0.0 ! Zero growth respiration
      cycle ! Skip the following growth processes
    endif

    ! Get carbon from NSC pool
    call fetch_CN_for_growth(cc,Cgrowth,Nsupply) ! Weng, 2017-10-19
    associate (sp => spdata(cc%species))
      ! Allocate carbon to the plant pools
      ! calculate the carbon spent on growth of leaves and roots
      LF_deficit = max(0., cc%bl_max - cc%bl)
      FR_deficit = max(0., cc%br_max - cc%br)
      LFR_deficit = LF_deficit + FR_deficit
      G_LFR = max(min(LF_deficit + FR_deficit, f_LFR_max * Cgrowth),0.) ! (1.- Wood_fract_min)
      !! and distribute it between roots and leaves
      dBL = min(max(0.,(G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/ &
                       (cc%bl_max + cc%br_max)), G_LFR)
      if((G_LFR-dBL) > FR_deficit) dBL = G_LFR - FR_deficit
      dBR  = G_LFR - dBL
      ! calculate carbon spent on growth of sapwood growth
      if(cc%layer == 1 .AND. cc%age > sp%AgeRepro)then
          dSeed = sp%v_seed * (Cgrowth - G_LFR)
          dBSW  = (1.0-sp%v_seed)* (Cgrowth - G_LFR)
      else
          dSeed= 0.0
          dBSW = Cgrowth - G_LFR
      endif

      ! For grasses, temporary, added addtional allocation to seeds for grasses in all layers
      if(sp%lifeform == 0) then
          dSeed = dSeed + sp%v_seed * G_LFR
          dBR   = (1.0 - sp%v_seed) * dBR
          dBL   = (1.0 - sp%v_seed) * dBL
          G_LFR = dBL + dBR
      endif

      ! ------------Updated 2026-01-21 from a scheme 2019-05-21 -----------
      !! Nitrogen adjustment on allocations between wood and leaves+roots
      ! same ratio reduction for leaf, root, and seed if(Nsupply < Ndemand),
      !! Nitrogen demand by leaves, roots, and seeds (Their C/N ratios are fixed.)
      Ndemand = dBL/sp%CNleaf0 + dBR/sp%CNroot0 + dSeed/sp%CNseed0 + dBSW/sp%CNwood0
      if(Ndemand > 0.0 .and. Nsupply < Ndemand) then
        r_N_SD = MAX(0.0, Nsupply/Ndemand) ! N supply-demand ratio
        cc%extraC   = (1.0-r_N_SD) * (dBL+dBR+dSeed)
        dBSW =  dBSW + cc%extraC
        dBR  =  r_N_SD * dBR
        dBL  =  r_N_SD * dBL
        dSeed=  r_N_SD * dSeed
      else
        cc%extraC   = 0.0
      endif

      ! Update plant pools
      cc%NSC   = cc%NSC   - dBR - dBL -dSeed - dBSW
      cc%resg  = 0.5 * (dBR+dBL+dSeed+dBSW) !  daily
      !update biomass pools
      cc%bl    = cc%bl    + dBL
      cc%br    = cc%br    + dBR
      cc%bsw   = cc%bsw   + dBSW
      cc%seedC = cc%seedC + dSeed
      ! Update leaf age
      cc%leafage = cc%bl/(dBL+cc%bl)*cc%leafage + 1.0/365.0

      !!update nitrogen pools, Nitrogen allocation
      cc%leafN = cc%leafN + dBL   /sp%CNleaf0
      cc%rootN = cc%rootN + dBR   /sp%CNroot0
      cc%seedN = cc%seedN + dSeed /sp%CNseed0
      cc%swN = cc%swN + f_N_add * cc%NSN + &
         (Nsupply - dBL/sp%CNleaf0 - dBR/sp%CNroot0 - dSeed/sp%CNseed0)
      !extraN = max(0.0,cc%swN+cc%hwN - (cc%bsw+cc%bHW)/sp%CNwood0)
      extraN   = max(0.0,cc%swN - cc%bsw/sp%CNwood0)
      cc%swN = cc%swN - extraN
      cc%NSN   = cc%NSN   + extraN - f_N_add*cc%NSN - Nsupply !! update NSN

      ! accumulated C allocated to leaf, root, and wood
      cc%NPPleaf = cc%NPPleaf + dBL
      cc%NPProot = cc%NPProot + dBR
      cc%NPPwood = cc%NPPwood + dBSW

      ! Update plant architecture
      DBH0 = cc%DBH   ! Keep previous DBH
      call BM2Architecture(cc, cc%bsw+cc%bHW)
      ! Update Ktrunk with new sapwood
      k = Max(MIN(cc%Nrings, Ysw_max),1)
      cc%Kx(k)   = NewWoodKx(cc)
      cc%Lring(k)= HT2Lpath(cc%height)
      cc%Ktrunk  = cc%Ktrunk + 0.25*PI*(cc%DBH**2-DBH0**2)*cc%Kx(k)/cc%Lring(k)
    end associate ! F2003
  enddo
  cc => null()
end subroutine vegn_growth ! daily

!========================================================================
subroutine vegn_cohort_update(vegn)
  ! Update derived variables of cohorts at daily time step
  ! after daily growth and tissue turnover, 12/04/2025
  type(vegn_tile_type), intent(inout) :: vegn

  !-------local var----------
  type(cohort_type), pointer :: cc    ! current cohort
  integer :: i, j, L

  ! Calculate root area and each canopy layer's LAI, CAI and gaps
  vegn%LAI_L  = 0.0
  vegn%CAI_L  = 0.0
  vegn%f_gap  = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
      ! update accumulative LAI for each corwn layer
       L = Max(1, Min(cc%layer,CLmax)) ! between 1 ~ CLmax
       vegn%LAI_L(L) = vegn%LAI_L(L) + cc%Aleaf  * cc%nindivs/(1.0-sp%f_cGap)
       vegn%CAI_L(L) = vegn%CAI_L(L) + cc%Acrown * cc%nindivs / (1.0-sp%f_cGap)
       vegn%f_gap(L) = vegn%f_gap(L) + cc%Acrown * cc%nindivs * sp%f_cGap

      ! Root area and vertical distribution
       cc%rootarea = cc%br * sp%SRA
       do j=1,soil_L
         cc%ArootL(j) = cc%rootarea * sp%root_frac(j)
       enddo
     end associate
  enddo

  ! Update cohort derived variables
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! Update bl_max, br_max, and NSNmax
     call update_max_LFR_NSN(cc)
     ! Update plant hydraulic states
     call Update_plant_hydro_vars(cc)
     call Plant_water2psi_exp(cc)
     ! Update leaf variables
     cc%Aleaf = BL2Aleaf(cc%bl,cc)
     cc%LAI   = Aleaf2LAI(vegn%CAI_L(cc%layer),cc)
     cc%Aleafmax = Max(cc%Aleafmax, cc%Aleaf)
  enddo
end subroutine vegn_cohort_update

!===============================================================================
subroutine grass_thinning(vegn)
  ! Kill extra grasses, 12/03/2025
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars -------------
  type(cohort_type), pointer :: cc    ! current cohort
  real :: TopGC, r_rem ! for grass thinning
  integer :: i

  ! Calculate grass CA (TopGC) in the top layer for grass thinning
  TopGC = 0.0
  do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      if(cc%layer > 1) exit
      associate ( sp => spdata(cc%species))
        if(sp%lifeform == 0 .and. cc%layer == 1 .and. cc%status == LEAF_ON) then
          TopGC = TopGC + cc%Acrown*cc%nindivs
        endif
      end associate
  enddo

  ! Reduce cohort density when TopGC is 0.05 greater than MaxGrassCA
  if(TopGC > MaxGrassCA + 0.05)then ! 0.05 is used for delaying correction
    r_rem = 1.0 - MaxGrassCA/TopGC
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      if(cc%layer > 1) exit
      associate ( sp => spdata(cc%species))
      if(sp%lifeform == 0 .and. cc%layer == 1 .and. cc%status == LEAF_ON)then
        call plant2soil(vegn,cc,r_rem*cc%nindivs)
        cc%nindivs = (1.0 - r_rem) * cc%nindivs
      endif
      end associate
    enddo
  endif
end subroutine grass_thinning ! daily

!========= Calculate carbon and nitrogen supply ==========================
subroutine fetch_CN_for_growth(cc,Cgrowth,Nsupply)
  !@sum Fetch C from labile C pool according to the demand of leaves and fine roots,
  !@+   and the push of labile C pool
  !@+   Daily call.
  !@+   added by Weng, 12-06-2016

  implicit none
  type(cohort_type), intent(inout) :: cc
  real, intent(out):: Cgrowth, Nsupply

  !------local var -----------
  logical :: woody
  logical :: dormant,growing
  real :: NSCtarget, bl_max, br_max
  real :: C_push, C_pull, growthC
  real :: N_push, N_pull, growthN
  real :: LFR_rate

  associate ( sp => spdata(cc%species) )
    ! Fetch C from labile C pool if it is in the growing season
    NSCtarget = 3.0 * (cc%bl_max + cc%br_max)      ! kgC/tree
    LFR_rate = sp%LFR_rate ! 1.0 !  1.0/16.0 ! filling rate/day

    C_pull = LFR_rate * (Max(cc%bl_max - cc%bl,0.0) +   &
              Max(cc%br_max - cc%br,0.0))
    N_pull = LFR_rate * (Max(cc%bl_max - cc%bl,0.0)/sp%CNleaf0 +  &
              Max(cc%br_max - cc%br,0.0)/sp%CNroot0)

    C_push = max(0.0,cc%nsc-0.1*NSCtarget)/(days_per_year*sp%tauNSC) ! max(cc%nsc-NSCtarget, 0.0)/(days_per_year*sp%tauNSC)
    N_push = max(0.0,cc%NSN)/(days_per_year*sp%tauNSC)

    Cgrowth = Min(max(0.02*cc%NSC,0.0), C_pull + C_push)
    Nsupply = Min(max(0.02*cc%NSN,0.0), N_pull + N_push)

  end associate

end subroutine fetch_CN_for_growth

!================================================================
subroutine update_max_LFR_NSN(cc)
  !@sum: Daily call for calculating bl_max and br_max
  !@+   added by Weng, 08-19-2022
  implicit none
  type(cohort_type), intent(inout) :: cc

  !----- local vars ---------
  real :: BL_c, BL_u, f_CO2

  ! Update bl_max and br_max daily, Weng 2014-01-23, 2021-06-04, 08/24/2022
  ! The new updates allow a gradual increase of BLmax when a tree enters
  ! the canopy layer and an abrupt increase for grasses.
  associate ( sp => spdata(cc%species) )
    f_CO2 =Min(Max(SQRT(cc%CO2_c/400.0), 0.5),2.0) ! CO2 effects on max LAI, for FACE-MDS
    BL_c = f_CO2 * sp%LAImax * sp%LMA * cc%Acrown * (1.0-sp%f_cGap)
    if (cc%layer > 1) then
      BL_u = BL_c/cc%layer
    else ! cc%layer = 1
      BL_u = BL_c / (1+cc%layer)            ! Woody plants only
    endif
    cc%bl_max = BL_u + min(1., cc%topyear/sp%transT) * (BL_c - BL_u)
!    if(cc%layer > 1) cc%bl_max = sp%LMA * 1.0 * cc%Acrown * (1.0-sp%f_cGap)
    cc%br_max = BLmax2BRmax(cc)
    cc%NSNmax = ccNSNmax(cc)
  end associate
end subroutine update_max_LFR_NSN

!============================================================================
! Updated by Weng, 06-04-2021
subroutine vegn_phenology(vegn) ! daily step
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc
  integer :: i,j !,var_acc
  real    :: gdd_ON, Tc_OFF
  real    :: totC, totN, ccNSC, ccNSN
  logical :: PhenoON, PhenoOFF

  ! -------------- update vegn GDD and tc_pheno ---------------
  vegn%tc_pheno = vegn%tc_pheno * 0.8 + vegn%Tc_daily * 0.2 ! C
  do i = 1, vegn%n_cohorts
    cc=>vegn%cohorts(i)
    associate (sp => spdata(cc%species) )
      if(sp%phenotype == 0) then ! Deciduous speices
        if(cc%status == LEAF_ON)then
           cc%ngd = Min(366, cc%ngd + 1)
           if(cc%ngd > Days_thld) cc%ALT = cc%ALT + MIN(0.,vegn%tc_pheno-sp%tc0_off)
           if(cc%dailyWdmd > 0.0) cc%AWD = 0.9*cc%AWD + 0.1 * (cc%dailyTrsp/cc%dailyWdmd)
        else  ! cc%status == LEAF_OFF
           cc%ndm = cc%ndm + 1
           if(vegn%tc_pheno<T0_chill) cc%ncd = cc%ncd + 1
           ! Keep gdd as zero in early non-growing season when days < Days_thld
           ! var_acc = merge(1, 0, cc%ndm > Days_thld)
           if(cc%ndm>Days_thld)cc%gdd = cc%gdd + max(0.0,vegn%tc_pheno-T0_gdd) ! * var_acc
        endif ! cc%status
      endif   ! sp%phenotype == 0

    end associate
  enddo

  ! --------- Change pheno status ----------------------------
  ! Turn ON the phenology of deciduous species
  do i = 1,vegn%n_cohorts
    cc => vegn%cohorts(i)
    associate (sp => spdata(cc%species) )
      cc%firstday = .false.
      if(sp%phenotype == 0) then ! Deciduous speices
        gdd_ON  = sp%gdd_par1 + sp%gdd_par2 * exp(sp%gdd_par3*cc%ncd) ! GDD threshold for leaf green-up
        PhenoON = ( cc%status/=LEAF_ON )                              &
         !.and.(cc%gdd>sp%gdd_crit    .and. vegn%tc_pheno>sp%tc0_on)  &  ! Thermal conditions
          .and.(cc%gdd > gdd_ON       .and. vegn%tc_pheno>sp%tc0_on)  &  ! Thermal conditions
          .and.(vegn%thetaS>sp%betaON .and. cc%Ndm > Days_thld)       &  ! Water
          .and.(.NOT.(sp%lifeform==0  .and. cc%layer > MaxGrassLyr))     ! If grasses, layer <= MaxGrassLyr

        if(PhenoON)then
          ! Setup phenology status
          cc%status = LEAF_ON ! Turn on a growing season
          cc%firstday = .True.
          cc%gdd_ON = gdd_ON
          cc%gdd = 0.0
          cc%ncd = 0
          cc%ndm = 0
          cc%AWD = 1.0 ! Accumulative water available ratio

          ! Setup water pressure
          cc%psi_s0   = maxval(vegn%psi_soil(:))
          cc%psi_stem = cc%psi_s0
          cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height)
        endif

        ! Reset deciduous grasses at the first day of a growing season
        if(sp%lifeform ==0 .and. (cc%firstday .and. cc%age>0.5))then
            ccNSC = (cc%NSC +cc%bl +  cc%bsw  +cc%bHW  +cc%br   +cc%seedC) * cc%nindivs
            ccNSN = (cc%NSN +cc%leafN+cc%swN+cc%hwN+cc%rootN+cc%seedN) * cc%nindivs
            cc%nindivs = MIN(ccNSC /sp%s0_plant, ccNSN/(sp%s0_plant/sp%CNroot0))
            totC = ccNSC / cc%nindivs
            totN = ccNSN / cc%nindivs
            call setup_seedling(cc,totC,totN)
        endif
      else
        cc%status=LEAF_ON ! Evergree species
      endif
    end associate
  enddo  ! vegn%n_cohorts
  if(any(vegn%cohorts(:)%firstday)) call relayer_cohorts(vegn)

  ! ---------- OFF of a growing season --------------
  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate (sp => spdata(cc%species) )
       if(sp%phenotype == 0) then
         ! Critical temperature trigering offset of phenology
         Tc_OFF = sp%tc0_off - 5. * exp(-0.05*(cc%ngd-N0_GD))
         PhenoOFF = (cc%status == LEAF_ON .and. cc%NGD > Days_thld .and.  & ! Minimum days of a growing season
              ((cc%ALT < cold_thld .and. vegn%tc_pheno < Tc_OFF) .or.  & ! Cold-deciduous
              (vegn%thetaS < sp%betaOFF)  & ! Drought-deciduous
              !(cc%AWD < sp%AWD_crit)      & ! Drought-deciduous
              ))
         if(PhenoOFF )then
            cc%status = LEAF_OFF  ! Turn off a growing season
            cc%Tc_OFF = Tc_OFF
            cc%gdd = 0.0          ! Start to count a new cycle of GDD
            cc%ngd = 0
            cc%ALT = 0.0
            cc%AWD = 1.0 ! Accumulative water available ratio
         endif
         call Seasonal_fall(cc,vegn) ! leaf fall
       endif
     end associate
  enddo ! vegn%n_cohorts
end subroutine vegn_phenology

! ============================================================================
subroutine vegn_tissue_turnover(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  !-------local var
  type(cohort_type), pointer :: cc    ! current cohort
  real :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real :: alpha_L, alpha_S, alpha_R   ! turnover rates of leaves, stems, and roots
  real :: dBL, dBR, dBStem  ! leaf and fine root carbon tendencies
  real :: dNL, dNR, dNStem  ! leaf and fine root nitrogen tendencies
  real :: dAleaf ! leaf area decrease due to dBL
  integer :: i

  ! update plant carbon and nitrogen for all cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
     ! turnover rate, fraction per day
     alpha_L = MIN(.2, Max(2.*cc%leafage/sp%leafLS - 1., 0.)) ! Leaves
     alpha_R = sp%alpha_FR /days_per_year                     ! Roots
     alpha_S = 0.0                                            ! Woody stems
     ! Grass stems
     if(sp%lifeform == 0) alpha_S = alpha_L

     ! Tissue changes
     dBL    = cc%bl    * alpha_L
     dNL    = cc%leafN * alpha_L
     dBStem = cc%bsw   * alpha_S
     dNStem = cc%swN   * alpha_S
     dBR    = cc%br    * alpha_R
     dNR    = cc%rootN * alpha_R

     !    Retranslocation to NSC and NSN
     cc%nsc = cc%nsc + l_fract  * (dBL + dBR + dBStem)
     cc%NSN = cc%NSN + retransN * (dNL + dNR + dNStem)

     ! Update leaf age
     if(cc%bl>0.0001)then
       cc%leafage = (1.0 - dBL/cc%bl)*cc%leafage
     else
       cc%leafage = 0.d0
     endif
     !    update plant pools
     cc%bl    = cc%bl    - dBL
     cc%bsw   = cc%bsw   - dBStem
     cc%br    = cc%br    - dBR

     cc%leafN = cc%leafN - dNL
     cc%swN   = cc%swN   - dNStem
     cc%rootN = cc%rootN - dNR

     ! update NPP for leaves, fine roots, and wood
     cc%NPPleaf = cc%NPPleaf - l_fract * dBL
     cc%NPProot = cc%NPProot - l_fract * dBR
     cc%NPPwood = cc%NPPwood - l_fract * dBStem

     ! Put C and N into soil pools
     dAleaf = BL2Aleaf(dBL,cc)
     loss_coarse  = (1.-l_fract) * cc%nindivs * (dBL - dAleaf * LMAmin    + dBStem)
     loss_fine    = (1.-l_fract) * cc%nindivs * (dBR + dAleaf * LMAmin)
     lossN_coarse = (1.-retransN)* cc%nindivs * (dNL - dAleaf * sp%LNbase + dNStem)
     lossN_fine   = (1.-retransN)* cc%nindivs * (dNR + dAleaf * sp%LNbase)

     vegn%dailyLFLIT = vegn%dailyLFLIT + (1.-l_fract) * cc%nindivs * dBL
     vegn%SOC(1) = vegn%SOC(1)   +  &
                        fsc_fine * loss_fine + fsc_wood * loss_coarse
     vegn%SOC(2) = vegn%SOC(2) +  &
                         ((1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse)
     vegn%SON(1)  = vegn%SON(1) +    &
                          fsc_fine * lossN_fine + fsc_wood * lossN_coarse
     vegn%SON(2) = vegn%SON(2) + &
                          (1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse

     !    annual N from plants to soil
     vegn%NorgP2S = vegn%NorgP2S + lossN_fine + lossN_coarse

    END ASSOCIATE
  enddo

end subroutine vegn_tissue_turnover

!========================================================================
! Starvation due to low NSC or NSN, daily
subroutine vegn_daily_starvation (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars --------
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  integer :: i, k
  type(cohort_type), pointer :: cc
  type(cohort_type), dimension(:),pointer :: ccold, ccnew

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species))
    !Mortality due to starvation
    deathrate = 0.0
    !if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max .OR.(cc%layer >1 .and. sp%lifeform ==0)) then
    if (cc%nsc < 0.01*cc%bl_max )then ! .OR. cc%NSN < 0.01*cc%bl_max/sp%CNleaf0
         deathrate = 1.0
         deadtrees = cc%nindivs * deathrate !individuals / m2
         write(*,*)"Daily starvation is happening"
         ! Carbon and Nitrogen from plants to soil pools
         call plant2soil(vegn,cc,deadtrees)
         ! update cohort individuals
         cc%nindivs = 0.0 ! cc%nindivs*(1.0 - deathrate)
     else
         deathrate = 0.0
     endif
     end associate
  enddo
  ! Remove the cohorts with 0 individuals
  !call kill_lowdensity_cohorts(vegn)
end subroutine vegn_daily_starvation

!========= Leaf and stem fall ==========================
subroutine Seasonal_fall(cc,vegn)
  !@sum leaf and stem fall for deciduous plants, including deciduous trees and grasses
  !@+   DAILY call.
  !@+   added by Weng, 12-03-2017

  implicit none
  type(cohort_type), intent(inout) :: cc
  type(vegn_tile_type), intent(inout) :: vegn
    !------local var -----------
  real    :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real    :: dAleaf, dBL, dBR, dNL, dNR, dBStem, dNStem, dWLeaf, dWStem      ! per day
  real    :: leaf_fall_rate, root_mort_rate      ! per day

  leaf_fall_rate = 0.05
  root_mort_rate = 0.025
  !End a growing season: leaves fall for deciduous
  associate (sp => spdata(cc%species) )
  if(cc%status == LEAF_OFF .AND. cc%bl > 0.0)then
     dBL = min(leaf_fall_rate * cc%bl_max, cc%bl)
     dBR = min( root_mort_rate * cc%br_max, cc%br)  ! Just for test: keep roots
     if(sp%lifeform == 0)then  ! grasses
         dBStem = MIN(1.0,dBL/cc%bl) * cc%bsw
         dNStem = MIN(1.0,dBL/cc%bl) * cc%swN
         dWStem = MIN(1.0,dBL/cc%bl) * cc%W_stem
     else
         dBStem = 0.0 ! trees
         dNStem = 0.0 ! trees
         dWStem = 0.0
     endif
     ! Nitrogen and water out
     dNL = dBL/cc%bl * cc%leafN !dBL/sp%CNleaf0
     dWLeaf = cc%W_leaf*dBL/cc%bl

     if(cc%br>0)then
        dNR = dBR/cc%br * cc%rootN !dBR/sp%CNroot0
     else
        dNR = 0.0
     endif

     dAleaf = BL2Aleaf(dBL,cc)
#ifdef Hydro_test
     ! Put plant water into the first soil layer
     vegn%wcl(1) = vegn%wcl(1) + cc%nindivs*(dWLeaf+dWStem)/(thksl(1)*1000.0)
#endif

     !Retranslocation to NSC and NSN
     cc%nsc = cc%nsc + l_fract  * (dBL + dBR + dBStem)
     cc%NSN = cc%NSN + retransN * (dNL + dNR + dNStem)
     !update plant pools
     cc%bl    = cc%bl  - dBL
     cc%br    = cc%br  - dBR
     cc%bsw   = cc%bsw - dBStem ! for grass
     if(cc%bl<0.0001) cc%leafage = 0.0

     cc%leafN = cc%leafN - dNL
     cc%rootN = cc%rootN - dNR
     cc%swN   = cc%swN - dNStem

     !update NPP for leaves, fine roots, and wood
     cc%NPPleaf = cc%NPPleaf - l_fract * dBL
     cc%NPProot = cc%NPProot - l_fract * dBR
     cc%NPPwood = cc%NPPwood - l_fract * dBStem
     cc%Aleaf   = BL2Aleaf(cc%bl,cc)
     cc%LAI     = vegn%CAI_L(cc%layer)*cc%Aleaf/(cc%Acrown *(1.0-sp%f_cGap))

     !put C and N into soil pools:  Substraction of C and N from leaf and root pools
     loss_coarse  = (1.-l_fract) * cc%nindivs * (dBStem+dBL - dAleaf * LMAmin)
     loss_fine    = (1.-l_fract) * cc%nindivs * (dBR        + dAleaf * LMAmin)
     lossN_coarse = (1.-retransN)* cc%nindivs * (dNStem+dNL - dAleaf * sp%LNbase)
     lossN_fine   = (1.-retransN)* cc%nindivs * (dNR        + dAleaf * sp%LNbase)

     vegn%dailyLFLIT = vegn%dailyLFLIT + (1.-l_fract) * cc%nindivs * dBL
     vegn%SOC(1) = vegn%SOC(1) +  &
                      fsc_fine * loss_fine + fsc_wood * loss_coarse
     vegn%SOC(2) = vegn%SOC(2) +   &
                      (1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse
     vegn%SON(1)  = vegn%SON(1) +    &
                       fsc_fine * lossN_fine + fsc_wood * lossN_coarse
     vegn%SON(2) = vegn%SON(2) +   &
                       (1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse

     !annual N from plants to soil
     vegn%NorgP2S = vegn%NorgP2S + lossN_fine + lossN_coarse
  endif
  end associate

end subroutine Seasonal_fall

!========================================================================
! Weng, 08/23/2022: time counter, update cohort ages and the time of them
!                   staying in the first layer
subroutine vegn_age (vegn,t_yr) ! daily
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: t_yr ! step length (year), 1.0/365.0

  !----- local var --------------
  type(cohort_type),pointer :: cc
  integer :: i

  vegn%age = vegn%age + t_yr ! Tile age
  ! Update cohort ages
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     cc%CO2_c = vegn%CO2_c ! for FACE-MDS
     ! Update cohort age
     cc%age = cc%age + t_yr
     ! Update time in the top layer
     if (cc%layer == 1) cc%topyear = cc%topyear + t_yr
  enddo
end subroutine vegn_age

! =========================================================================
subroutine vegn_N_deposition(vegn, dt_daily)
  ! Weng, 05/15/2023: Nitrogen deposition, daily
  type(vegn_tile_type), intent(inout) :: vegn
  real                , intent(in)    :: dt_daily

  ! Update mineral N pool (mineralN)
  vegn%mineralN = vegn%mineralN + vegn%N_input * dt_daily
end subroutine vegn_N_deposition

!==========================================================================
!========================= Annual subroutines =============================

!=======================================================================
! the reproduction of each canopy cohort, yearly time step
! calculate the new cohorts added in this step and states:
! tree density, DBH, woddy and fine biomass
subroutine vegn_reproduction (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc ! parent and child cohort pointers
  type(cohort_type), dimension(:),pointer :: ccold, ccnew   ! pointer to old cohort array
  integer,dimension(MSPECIES) :: reproPFTs
  real,   dimension(MSPECIES) :: seedC, seedN ! seed pool of productible PFTs
  real :: failed_seeds, N_failedseed !, prob_g, prob_e
  real :: totC,totN  ! the biomass of a new seedling
  integer :: newcohorts, matchflag, nPFTs ! number of new cohorts to be created
  integer :: nCohorts, istat
  integer :: i, j, k, n ! cohort indices

  ! Looping through all reproducible cohorts and Check if reproduction happens
  reproPFTs = -999 ! the code of reproductive PFT
  seedC = 0.0
  seedN = 0.0
  nPFTs = 0
  do k=1, vegn%n_cohorts
     cc => vegn%cohorts(k)
     if(cohort_can_reproduce(cc))then
        matchflag = 0
        do i=1,nPFTs
           if(cc%species == reproPFTs(i))then
               seedC(i) = seedC(i) + cc%seedC  * cc%nindivs
               seedN(i) = seedN(i) + cc%seedN  * cc%nindivs
               ! reset parent's seed C and N
               cc%seedC = 0.0
               cc%seedN = 0.0
               ! Update flag
               matchflag = 1
               exit
           endif
        enddo
        if(matchflag==0)then ! a new PFT, put the seeds to the next place
            nPFTs            = nPFTs + 1 ! update the number of reproducible PFTs
            reproPFTs(nPFTs) = cc%species ! PFT number
            seedC(nPFTs)     = cc%seedC * cc%nindivs ! seed carbon
            seedN(nPFTs)     = cc%seedN * cc%nindivs ! seed nitrogen
            ! reset parent's seed C and N
            cc%seedC = 0.0
            cc%seedN = 0.0
        endif
     endif ! cohort_can_reproduce
  enddo ! k, vegn%n_cohorts

  ! Generate new cohorts
  newcohorts = nPFTs
  if (newcohorts >= 1) then   ! build new cohorts for seedlings
     ccold => vegn%cohorts ! keep old cohort information
     nCohorts = vegn%n_cohorts + newcohorts
     allocate(ccnew(1:nCohorts), STAT = istat)
     ccnew(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
     vegn%cohorts => ccnew

     ! set up new cohorts
     k = vegn%n_cohorts
     do i = 1, newcohorts
        k = k+1 ! increment new cohort index
        ! Copy old information to new cohort, Weng, 2021-06-02
        do n =1, vegn%n_cohorts ! go through old cohorts
          if(reproPFTs(i) == ccold(n)%species)then
            ccnew(k) = ccold(n) ! Use the information from partent cohort
            exit
          endif
        enddo

        ! Update new cohort information
        cc => vegn%cohorts(k)
        cc%species = reproPFTs(i)
        cc%ccID = MaxCohortID + i ! new cohort ID
        ! update child cohort parameters
        associate (sp => spdata(reproPFTs(i)))
          totC = sp%s0_plant
          cc%nindivs = seedC(i)/totC ! density
          if(seedC(i)>0.0)then
            totN = totC * seedN(i)/seedC(i)
          else
            totN = totC/10.0
          endif
          call setup_seedling(cc,totC,totN)
          ! Take water from the second soil layer
          vegn%wcl(2) = Max(0.0, vegn%wcl(2)-cc%nindivs*cc%W_stem/(thksl(2)*1000.0))

          !! seeds fail
          ! cc%nindivs = cc%nindivs * sp%prob_g * sp%prob_e
          !!       put failed seeds to soil carbon pools
          !        failed_seeds = 0.0 ! (1. - sp%prob_g*sp%prob_e) * seedC(i)!!

          !        vegn%litter = vegn%litter + failed_seeds
          !        vegn%SOC(1) = vegn%SOC(1) +        fsc_fine *failed_seeds
          !        vegn%SOC(2) = vegn%SOC(2) + (1.0 - fsc_fine)*failed_seeds

          !!      Nitrogen of seeds to soil SOMs
          !        N_failedseed= 0.0 ! (1.-sp%prob_g*sp%prob_e)   * seedN(i)
          !        vegn%SON(1)  = vegn%SON(1)   +        fsc_fine * N_failedseed
          !        vegn%SON(2) = vegn%SON(2)  + (1.0 - fsc_fine)* N_failedseed

          !       annual N from plants to soil
          !       vegn%NorgP2S = vegn%NorgP2S + N_failedseed
        end associate
     enddo
     deallocate (ccold)
     MaxCohortID = MaxCohortID + newcohorts
     vegn%n_cohorts = k
     ccnew => null()
  endif ! set up new born cohorts

end subroutine vegn_reproduction

!============================================
subroutine setup_seedling(cc,totC,totN)
  type(cohort_type), intent(inout) :: cc
  real,              intent(in)    :: totC
  real,              intent(in)    :: totN
  ! ----- local var -----------
  integer :: layer
  ! ----------------------------
  associate(sp=>spdata(cc%species))
     layer = max(1, cc%layer)
     if(sp%phenotype == 0)then
        cc%status = LEAF_OFF
     else
        cc%status = LEAF_ON
     endif
     ! Leaf age
     cc%leafage = 0.0
     ! Carbon pools
     cc%bl     = 0.0 * totC
     cc%br     = 0.1 * totC
     cc%bsw    = f_iniBSW * totC
     cc%bHW    = 0.0 * totC
     cc%seedC  = 0.0
     cc%nsc    = totC - cc%bsw -cc%br
     ! Nitrogen pools
     cc%leafN  = cc%bl/sp%CNleaf0
     cc%rootN  = cc%br/sp%CNroot0
     cc%swN  = cc%bsw/sp%CNwood0
     cc%hwN  = cc%bHW/sp%CNwood0
     cc%seedN  = 0.0
     cc%NSN    = totN - (cc%leafN+cc%rootN+cc%swN) !cc%br/sp%CNroot0

     ! Structure
     call BM2Architecture(cc,cc%bsw+cc%bHW )
     call update_max_LFR_NSN(cc)

     ! Zero cohort variables
     cc%topyear = 0.0
     cc%age     = 0.0

     ! Hydraulic states
     cc%transp= 0.0
     cc%WTC0  = 0.0
     cc%Kx    = 0.0
     cc%farea = 1.0
     cc%accH  = 0.0
     cc%plcH  = 0.0
     cc%Rring = 0.0
     cc%Lring = 0.0
     cc%Aring = 0.0

     ! First ring
     cc%Nrings  = 1
     cc%WTC0(1) = NewWoodWTC(cc)
     cc%Kx(1)   = NewWoodKx(cc)
     cc%Lring(1)= HT2Lpath(cc%height)
     cc%Aring(1)= PI * 0.25*cc%DBH**2
     cc%Ktrunk  = cc%Aring(1) * cc%Kx(1)/cc%Lring(1)
     cc%treeHU  = 0.0
     cc%treeW0  = cc%WTC0(1) * cc%Aring(1)
     cc%Asap    = cc%Aring(1)
     cc%Atrunk  = cc%Aring(1)

     ! Cohort hydraulic properties
     call Update_plant_hydro_vars(cc)
     cc%W_leaf = cc%Wmax_l
     cc%W_stem = cc%Wmax_s
     cc%W_dead = 0.0
     call Plant_water2psi_exp(cc)

  end associate
end subroutine setup_seedling

! ============================================================================
function cohort_can_reproduce(cc); logical cohort_can_reproduce
  type(cohort_type), intent(in) :: cc

  associate (sp => spdata(cc%species) )! F2003
  cohort_can_reproduce = (cc%layer == 1 .and. &
                          cc%nindivs > 0.0 .and. &
                          cc%age   > sp%AgeRepro.and. &
                          cc%seedC > sp%s0_plant .and. &
                          cc%seedN > sp%s0_plant/sp%CNseed0)
  end associate

end function

! ============================================================================
!------------------------Mortality------------------------------------
subroutine vegn_nat_mortality (vegn, deltat)
  ! TODO: update background mortality rate as a function of wood density (Weng, Jan. 07 2017)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! seconds since last mortality calculations, s

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  type(spec_data_type),   pointer :: sp
  real :: deadtrees ! number of trees that died over the time step
  integer :: i, k

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     cc%mu = mortality_rate(cc)
     !if(cc%mu>=0.9) write(*,*)"Bulk motality happens!!!",cc%species,i
     !deadtrees = cc%nindivs*(1.0-exp(0.0-cc%mu*deltat/seconds_per_year)) ! individuals / m2
     deadtrees = cc%nindivs * cc%mu * deltat/seconds_per_year ! individuals / m2
     ! Carbon and Nitrogen from dead plants to soil pools
     call plant2soil(vegn,cc,deadtrees)
     ! Update plant density
     cc%nindivs = cc%nindivs - deadtrees
  enddo

end subroutine vegn_nat_mortality

!========================================================================
! Starvation due to low NSC and annual NPP
subroutine vegn_annual_starvation (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars --------
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  integer :: i, k
  type(cohort_type), pointer :: cc
  type(cohort_type), dimension(:),pointer :: ccold, ccnew

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species)  )
       ! Mortality due to starvation
       deathrate = 0.0
       !if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max .OR.(cc%layer >1 .and. sp%lifeform ==0)) then
       !if (cc%nsc < 0.01*cc%bl_max .OR. cc%annualNPP < 0.0) then ! .OR. cc%NSN < 0.01*cc%bl_max/sp%CNleaf0
       if (cc%nsc < 0.0001*cc%bl_max) then !  .OR. cc%annualNPP < 0.0, annualNPP < 0 is for grasses only
           deathrate = 1.0
           deadtrees = cc%nindivs * deathrate !individuals / m2
           ! write(*,*)"Yearly starvation: cNo., PFT:",i,cc%species,deathrate
           ! Carbon and Nitrogen from plants to soil pools
           call plant2soil(vegn,cc,deadtrees)
           ! update cohort individuals
           cc%nindivs = 0.0 ! cc%nindivs*(1.0 - deathrate)
       else
           deathrate = 0.0
       endif
     end associate
  enddo

end subroutine vegn_annual_starvation

! ===============================
subroutine plant2soil(vegn,cc,deadtrees)
  type(vegn_tile_type), intent(inout) :: vegn
  type(cohort_type),    intent(inout) :: cc
  real,                 intent(in)    :: deadtrees ! dead trees/m2

  ! ---- local vars --------
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse

     associate (sp => spdata(cc%species))
    ! Carbon and Nitrogen from plants to soil pools
     loss_coarse  = deadtrees * (cc%bHW + cc%bsw   + cc%bl    - cc%Aleaf*LMAmin)
     loss_fine    = deadtrees * (cc%nsc + cc%seedC + cc%br    + cc%Aleaf*LMAmin)
     lossN_coarse = deadtrees * (cc%hwN + cc%swN   + cc%leafN - cc%Aleaf*sp%LNbase)
     lossN_fine   = deadtrees * (cc%rootN+cc%seedN + cc%NSN   + cc%Aleaf*sp%LNbase)
#ifdef Hydro_test
     ! Assume water in plants goes to first layer of soil
     vegn%wcl(1) = vegn%wcl(1) +  deadtrees * (cc%W_leaf+cc%W_stem+cc%W_dead)/(thksl(1)*1000.0)
#endif
     vegn%SOC(1) = vegn%SOC(1) + fsc_fine *loss_fine + fsc_wood *loss_coarse
     vegn%SOC(2) = vegn%SOC(2) + (1.0-fsc_fine)*loss_fine + (1.0-fsc_wood)*loss_coarse

     vegn%SON(1) = vegn%SON(1) + &
                fsc_fine *lossN_fine +    fsc_wood *lossN_coarse
     vegn%SON(2) = vegn%SON(2) + &
                (1.-fsc_fine)*lossN_fine +(1.-fsc_wood)*lossN_coarse

     ! annual N from plants to soil
     vegn%NorgP2S = vegn%NorgP2S + lossN_fine + lossN_coarse

     end associate
end subroutine plant2soil

!============================================================================
real function mortality_rate(cc) result(mu) ! per year
  !@sum calculate cohort mortality/year, Ensheng Weng, 12/07/2021
  ! Mortality rate should be a function of growth rate, age, and environmental
  ! conditions. Here, we only used used a couple of parameters to calculate
  ! mortality as functions of social status, seedling size, and adult size.
  ! Grass is saprately defined.
  type(cohort_type),intent(inout) :: cc

  !-------local var -------------
  integer :: n ! the latest ring
  real :: f_L, f_S, f_D ! Layer, seeding, and size effects on mortality
  real :: mu_bg      ! Background mortality rate
  real :: mu_hydro   ! Hydraulic failure
  real :: mu_drought ! for UFL drought mortality
  real :: mu_add

  !---------------------
  mu_bg    = 0.0  ! Background mortality rate
  mu_hydro = 0.0
  mu_drought = 0.0
  associate ( sp => spdata(cc%species))
    n = MIN(cc%Nrings, Ysw_max)
    f_L = sp%A_un * SQRT(Max(0.0, cc%layer-1.0)) ! Layer effects (0~ infinite)
    f_S = sp%A_sd * exp(sp%B_sd*cc%dbh) + 1.0    ! Understory seedling
    f_D = 1.0 + sp%A_DBH/(1.+exp((sp%D0mu-cc%dbh)/sp%B_DBH)) ! Size effects (big tees, U-shaped)
    ! Background mortality rate
    mu_bg = Min(0.5,sp%mu0_topL * (1.d0+f_L*f_S)*f_D) ! per year

    if(DO_DroughtMu)then
      ! Annual drought mortality, From Lichstein et al. 2024 (J. Ecology)
      ! It can be turned off by setting a large sp%W_mu0
      if(cc%totDemand>0.00001)then
        cc%w_scale = cc%annualTrsp/cc%totDemand
      else
        cc%w_scale = 1.0
      endif
      mu_drought = 1.0/(1.0 + exp(-20.0*(1.0 - cc%w_scale - sp%W_mu0)))
    endif

#ifdef Hydro_test
    ! Trunk hydraulic failure probability
    mu_hydro = exp(sp%s_hu * (1.0 - min(1.,cc%treeHU/cc%treeW0)))
    !mu_hydro = Max(0., 1. - cc%farea(n))
    !mu_hydro = Max(0., 1. - cc%Asap/cc%Acrown/(sp%LAImax*sp%phiCSA))
#endif

  ! Total mortality rate:
  mu_add = mu_hydro + mu_drought - mu_hydro*mu_drought
  mu = mu_bg + (1.0 - mu_bg) * mu_add
  end associate

end function mortality_rate

!============================================================================
!-----------------------Plant Hydraulics------------------------------
subroutine vegn_hydraulic_states(vegn, deltat)
  ! Update plant hydraulic states, yearly time step
  ! Author: Ensheng Weng, 2021-03-15, updated 2023-10-8
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! seconds since last mortality calculations, s

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  type(spec_data_type),   pointer :: sp
  !real :: Trsp_sap ! Asap normalized tranp rate, m
  real :: trsp_ring ! Water flow amount of a ring, ton
  real :: funcA,Lmax
  integer :: i, j, k

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species))
       ! Set up the first year seedling
       if(cc%Nrings == 1)then
          cc%WTC0(1) = NewWoodWTC(cc)
          cc%Kx(1)   = NewWoodKx(cc)
          cc%farea(1) = 1.0
          cc%accH(1)  = 0.0
          cc%plcH(1)  = 0.0
          cc%Rring(1) = cc%DBH/2.0
          cc%Lring(1) = HT2Lpath(cc%height)
          cc%Aring(1) = PI * cc%Rring(1)**2
       endif

       ! Set up a new ring and move previous years' states inward
       if(cc%DBH > cc%DBH_ys)then
         cc%Nrings = cc%Nrings + 1 ! A new ring
         if(cc%Nrings > Ysw_max)then
            do j=2, Ysw_max
               cc%WTC0(j-1)  = cc%WTC0(j)
               cc%Kx(j-1)    = cc%Kx(j)
               cc%farea(j-1) = cc%farea(j)
               cc%accH(j-1)  = cc%accH(j)
               cc%plcH(j-1)  = cc%plcH(j)
               cc%Rring(j-1) = cc%Rring(j)
               cc%Lring(j-1) = cc%Lring(j)
               cc%Aring(j-1) = cc%Aring(j)
            enddo
         endif
         k = MIN(cc%Nrings, Ysw_max)
         ! WTC0 and Kx represent scientific hypotheses. They can be constant,
         ! or functions of environmental conditions, growht rates, etc.
         cc%WTC0(k) = NewWoodWTC(cc)
         cc%Kx(k)   = NewWoodKx(cc)
         ! Other cohort variables of the new ring
         cc%farea(k) = 1.0
         cc%accH(k)  = 0.0
         cc%plcH(k)  = 0.0
         cc%Rring(k) = cc%DBH/2.0
         cc%Lring(k) = HT2Lpath(cc%height)
         if(k>1)then
            cc%Aring(k) = PI * Max(0.0,cc%Rring(k)**2 - (cc%DBH_ys/2.)**2)
         else
            cc%Aring(k) = PI * cc%Rring(k)**2 ! Only for the first year
         endif
         call calculate_Asap_Ktrunk (cc) ! Trunk conductance and sapwood area
       endif ! Setup the new ring

       ! Update hydrauolic usage in this year
       Lmax = HT2Lpath(cc%height) ! maxval(cc%Lring(:)) !
       cc%treeHU = 0.0
       cc%treeW0 = 0.0
       do k=1, MIN(cc%Nrings, Ysw_max)
         ! Lifetime water transported for functional xylem conduits
         funcA = cc%farea(k) * cc%Aring(k)
         if(funcA > 1.0e-8)then
           trsp_ring = 1.e-3 * cc%annualTrsp * cc%Kring(k)/cc%Ktrunk ! ton, per ring
           cc%accH(k) = cc%accH(k) + trsp_ring/funcA  ! m, for functional conduits only
           !cc%accH(k) = cc%accH(k) + 1.e-3 * cc%annualTrsp * cc%Kx(k)/Lmax/cc%Ktrunk
         endif

         ! Update each ring's functional fraction
         cc%farea(k) = 1. - exp(-r_DF * (1. - MIN(1.0,(cc%accH(k)+cc%plcH(k))/cc%WTC0(k))))
         !cc%farea(k) = 1. - 1./(1. + exp(r_DF * (1. - cc%accH(k)/cc%WTC0(k))))

         ! Update tree hydraulic usage and WTC0
         funcA = cc%farea(k) * cc%Aring(k)
         cc%treeHU   = cc%treeHU + funcA * cc%accH(k)
         cc%treeW0   = cc%treeW0 + funcA * cc%WTC0(k)
       enddo
     end associate

     ! Update Asap and Ktrunk
     call calculate_Asap_Ktrunk (cc)

     ! This year's mortality rate, calculation only, for output
     cc%mu = mortality_rate(cc)
  enddo
  ! Sapwood conversion and tile variables
  call vegn_SW2HW_hydro(vegn)
  call vegn_sum_tile(vegn)
end subroutine vegn_hydraulic_states

!========================================================================
! Weng 2022-03-29 ! Updated 01/13/2023
subroutine Plant_water_dynamics_linear(vegn)     ! forcing,
  !type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real :: Q_air          ! transpiration amount in this step
  real :: psi_ht         ! Gravitational water pressure, MPa
  real :: psi_leaf,psi_stem,psi_sl,plc
  real :: k_rs(soil_L)   ! soil-root water conductance by soil layer
  real :: k_stem         ! The conductance of the tree at current conditions
  real :: sumK, sumPK, dpsi, psi_soil
  real :: W_psi_s0 ! Stem water content at psi_s0
  real :: Q_tot  ! Total soil water uptake by a plant
  real :: W_maxL ! Maximum available water in a soil layer for an individual
  real :: f_soil ! Fraction of water left in soil for next step
  integer :: i,j,k

  ! Plant water potentials, water fluxes and content
  do j = 1, vegn%n_cohorts
     cc => vegn%cohorts(j)
     associate ( sp => spdata(cc%species) )
       ! Refine plant water status based on water content
       cc%psi_s0 = maxval(vegn%psi_soil(:))
       call Plant_water2psi_exp(cc)
       cc%H_leaf = sp%CR_Leaf * cc%W_leaf  ! Leaf Capacitance
       cc%H_stem = sp%CR_Wood * cc%W_stem  ! Stem capacitance

       ! Calculate ! Water flux from stems to leaves
       Q_air  = cc%transp   ! /step_seconds
       psi_sl = (cc%psi_stem + cc%psi_leaf)/2
       plc = plc_function(psi_sl,sp%psi50_WD,sp%Kexp_WD)
       k_stem = cc%Ktrunk * plc
       psi_ht = HT2MPa(cc%height) ! MPa

       !Approximately estimate psi_leaf and Q_leaf
       ! find out "psi_leaf" that makes water flux from stems to leaves equal
       ! to Q_air (Transpiration)
       psi_leaf = (k_stem*step_seconds*(cc%psi_stem-psi_ht)+  &
            cc%H_leaf*cc%psi_leaf - Q_air) /(k_stem*step_seconds+cc%H_leaf)
       cc%Q_leaf = (cc%psi_stem - psi_leaf - psi_ht) * k_stem * step_seconds

       ! Adjust Q_leaf to make it reasonable
       cc%Q_leaf = MIN(Max(cc%Wmax_L - cc%W_leaf + Q_air, 0.0), &
                       Max(cc%W_stem - cc%Wmin_s, 0.0), cc%Q_leaf)

       !Update water content and psi
       cc%W_leaf = cc%W_leaf + cc%Q_leaf - Q_air
       cc%W_stem = cc%W_stem - cc%Q_leaf
       cc%psi_leaf = psi_leaf

       ! Xylem damage when plc is low
       if(plc <= plc_crit)then
         k = MIN(cc%Nrings, Ysw_max)
         cc%plcH(1:k) = cc%plcH(1:k) + sp%f_plc*cc%WTC0(1:k) * step_seconds/seconds_per_day
       endif

       ! Water fluxes from soil to stem base
       W_psi_s0 = cc%Wmax_s * exp(cc%psi_s0 * sp%CR_Wood)
       if(cc%W_stem < W_psi_s0) then
         sumK  = 0.0
         sumPK = 0.0
         do i=1, soil_L
           ! RAI(i) = cc%ArootL(i)/cc%Acrown
           ! k_rs(i)= vegn%K_soil(i)*SQRT(RAI(i))/(3.14159*thksl(i))
           k_rs(i)= vegn%K_soil(i)*cc%ArootL(i) ! per tree
           !if(vegn%psi_soil(i)<cc%psi_stem) k_rs(i)=0.0
           sumK = sumK + k_rs(i)
           sumPK= sumPK+ k_rs(i) * vegn%psi_soil(i)
         enddo
         psi_stem = (sumPK*step_seconds - cc%Q_leaf + cc%psi_stem*cc%H_stem)/ &
                    (sumK *step_seconds + cc%H_stem)
         ! Calculate water uptake by layer in theory
         do i=1, soil_L
           W_maxL =  vegn%freewater(i)*cc%ArootL(i)/max(1.e-6,vegn%ArootL(i))
           dpsi   = max(0.0, vegn%psi_soil(i) - psi_stem)
           cc%Q_soil(i) = Min(k_rs(i)*dpsi*step_seconds, W_maxL)
         enddo

         ! Check if the plant needs this ammount of water
         Q_tot = sum(cc%Q_soil)
         if(Q_tot > 0.0) &
           cc%Q_soil = cc%Q_soil * Min(1.0,(W_psi_s0-cc%W_stem)/Q_tot)

         ! Update stem and soil water content, and stem psi
         do i=1, soil_L
           cc%W_stem = cc%W_stem + cc%Q_soil(i)
           vegn%freewater(i) = vegn%freewater(i) - cc%Q_soil(i)*cc%nindivs
           vegn%wcl(i) = vegn%WILTPT + vegn%freewater(i)/(thksl(i)*1000.0)
         enddo
         cc%psi_stem = psi_stem ! log(cc%W_stem/cc%Wmax_s)/sp%CR_Wood
       endif ! Water uptake from soil

       ! Next step water supply for regulating stomata conductance
       Q_tot = sum(cc%Q_soil)
       f_soil = 1.0 - Q_tot/vegn%soilwater ! only when Q_tot is close to soil water
       cc%W_supply = PlantWaterSupply(cc,step_seconds) ! + f_soil*Q_tot
     end associate
  enddo
  ! Update soil free water to make sure it is the sum of all layers'
  vegn%soilwater = sum(vegn%freewater(:))

end subroutine Plant_water_dynamics_linear

!================================================================
subroutine vegn_SW2HW_hydro(vegn)
  ! Coverstion of sapwood to heartwood, yearly time step
  ! Author: Ensheng Weng, 2023-10-08
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  real :: Atrunk, D_hw, SW1, dSW, r_sw
  integer :: i

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! Woody plants only
     associate (sp => spdata(cc%species) )
     if(sp%lifeform>0 .and. cc%bsw > 0.01)then
       ! Heartwood diameter and new SW
       Atrunk = PI * cc%DBH * cc%DBH * 0.25
       D_hw   = 2.0 * SQRT(Max(0.0, (Atrunk - cc%Asap)/PI) )  ! Diameter, 2*r
       SW1    = sp%alphaBM * (cc%dbh**sp%thetaBM - D_hw**sp%thetaBM)

       ! Convertion of sapwood to heartwood (C, N, and water)
       if( cc%bsw > SW1)then
          dSW       = cc%bsw - SW1
          r_sw      = dSW / cc%bsw
          cc%bsw    = cc%bsw    - dSW
          cc%bhw    = cc%bhw    + dSW
          cc%swN    = cc%swN    - cc%swN    * r_sw
          cc%hwN    = cc%hwN    + cc%swN    * r_sw
          cc%W_stem = cc%W_stem - cc%W_stem * r_sw
          cc%W_dead = cc%W_dead + cc%W_stem * r_sw
       endif
     endif
     end associate

     ! Update other variables
     call Update_plant_hydro_vars(cc)
     call Plant_water2psi_exp(cc)
  enddo
end subroutine vegn_SW2HW_hydro

!--------------------------------------------------------------
subroutine vegn_SW2HW_fixedHv(vegn)
  ! Coverstion of sapwood to heartwood by fixed Hv, yearly time step
  ! from LM3-PPA, Weng, 08-19-2022
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: D_hw   ! diameter of heartwood at breast height, m
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: dSW    ! Sapwood to Heartwood, kgC/individual
  real :: dNS    ! Nitrogen from SW to HW
  integer :: i

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
       if(sp%lifeform>0)then ! woody plants
          CSAsw  = cc%bl_max/sp%LMA * sp%phiCSA * cc%height
          CSAtot = 0.25 * PI * cc%DBH**2
          CSAwd  = max(0.0, CSAtot - CSAsw)
          D_hw   = 2*sqrt(CSAwd/PI)
          BSWmax = sp%alphaBM * (cc%DBH**sp%thetaBM - D_hw**sp%thetaBM)
          dSW    = max(cc%bsw - BSWmax, 0.0)
          dNS    = dSW/cc%bsw *cc%swN
          ! update C and N of sapwood and wood
          cc%bHW   = cc%bHW   + dSW
          cc%bsw   = cc%bsw   - dSW
          cc%swN = cc%swN - dNS
          cc%hwN = cc%hwN + dNS
          cc%W_stem = cc%W_stem - cc%W_stem * dSW/cc%bSW
          cc%W_dead = cc%W_dead + cc%W_stem * dSW/cc%bSW

          !Update Atrunk and Asap
          D_hw = bm2dbh(cc%bHW,cc%species)
          cc%Atrunk = PI*(cc%DBH/2)**2
          cc%Asap = cc%Atrunk - PI*(D_hw/2)**2
       endif
     end associate

     ! Update other variables
     call Update_plant_hydro_vars(cc)
     call Plant_water2psi_exp(cc)
  enddo
end subroutine vegn_SW2HW_fixedHv

!===================================================================
real function PlantWaterSupply(cc,step_seconds) result(pws)
  !@sum: Calculate water supply for transpiration from plant water status
  ! Weng, 07/13/2022
  type(cohort_type),intent(in) :: cc
  real,             intent(in) :: step_seconds

  !----- Local vars ---------------
  real :: f0_sup
  real :: W_stem,W_leaf
  real :: k_stem,psi_stem,psi_leaf
  real :: psi0, dpsi
  real :: S_stem,S_leaf,wflux
  real :: step_base = 300. ! Seconds
  integer :: n_iterations,i
  logical :: do_simple_W_supply = .False. ! .True. !
  !---------------------

  ! Prameters
  associate ( sp => spdata(cc%species) )
    f0_sup = sp%f_supply * step_seconds/3600.0
    ! Stem water supply
    if(do_simple_W_supply)then ! An arbitrary assignment
      S_stem = f0_sup * Max((cc%W_stem - cc%Wmin_s),0.0)
    else ! Calculated as a function of woody properties
      n_iterations = int(step_seconds/step_base)
      psi_leaf = log(cc%W_leaf/cc%Wmax_l)/sp%CR_Leaf
      W_stem = cc%W_stem
      S_stem = 0.0
      do i =1, n_iterations
        psi_stem = log(W_stem/cc%Wmax_s)/sp%CR_Wood
        dpsi     = max(0.0, psi_stem - psi_leaf) ! sp%psi50_WD
        if(dpsi <= 0.02 .or. W_stem <= cc%Wmin_s)exit
        k_stem = cc%Ktrunk * plc_function(psi_stem,sp%psi50_WD,sp%Kexp_WD)
        wflux  = MIN(dpsi * k_stem * step_base, 0.5*(W_stem - cc%Wmin_S))
        S_stem = S_stem + wflux
        W_stem = W_stem - wflux
        !if(S_stem >= f0_sup * (cc%W_stem - cc%Wmin_s))exit
      enddo
    endif
  end associate
  !Leaf water supply
  S_leaf = max(f0_sup * (cc%W_leaf - cc%Wmin_L),0.0)
  ! Total plant water supply
  pws = S_leaf + S_stem

  !------------------------
  if(isnan(pws))then
     write(*,*)'cc%species',cc%species
     write(*,*)'cc%bl,cc%bsw',cc%bl,cc%bsw
     write(*,*)'cc%Wmin_L,cc%Wmax_L',cc%Wmin_L,cc%Wmax_L
     write(*,*)'cc%Wmin_s,cc%Wmax_s',cc%Wmin_s,cc%Wmax_s
     write(*,*)'cc%W_leaf,cc%W_stem',cc%W_leaf,cc%W_stem
     write(*,*)'cc%H_leaf,cc%H_stem',cc%H_leaf,cc%H_stem
     write(*,*)'cc%psi_leaf,cc%psi_stem',cc%psi_leaf,cc%psi_stem
     stop 'pws (PlantWaterSupply) is an NaN!'
  endif
end function PlantWaterSupply

!===================================================================
subroutine calculate_Asap_Ktrunk (cc)
  !@sum: Sapwood cross-sectional area and
  ! Total trunk conductance, Weng, 08/13/2022
  type(cohort_type),intent(inout) :: cc

  !----- Local vars ---------------
  integer :: i

  !---------------------
  !The contribution of each ring to the Ktrunk is calculated as where and
  ! how much water it lifts:
  !Kring=farea * Aringcc%Kx(i)/Lring *(Lring/Lmax)
  ! "Lring/Lmax" is the fraction of water lifts to tree height (or psi_wood)
  cc%Ktrunk = 0.0
  cc%Asap   = 0.0
  do i=1, MIN(cc%Nrings, Ysw_max)
    cc%Asap     = cc%Asap   + cc%Aring(i) * cc%farea(i)
    cc%Kring(i) = cc%farea(i)*cc%Aring(i) * cc%Kx(i)/HT2Lpath(cc%height) ! maxval(cc%Lring(:))
    cc%Ktrunk   = cc%Ktrunk + cc%Kring(i)
  enddo
  cc%Atrunk = 0.25 * PI * cc%DBH**2 ! Just for updating merged cohorts
end subroutine calculate_Asap_Ktrunk

!=================================================
! Weng: update soil root area layers, hydraulic variables, 03/29/2022
subroutine Update_plant_hydro_vars(cc)
  type(cohort_type), intent(inout) :: cc
  !----------local var ----------
  integer :: j

  associate (sp => spdata(cc%species) )
    ! Plant hydraulics-related variables
    cc%Wmax_l = cc%bl *(1.0/sp%rho_leaf - 1.0/rho_cellwall)*rho_H2O ! max leaf water, kg H2O
    cc%Wmax_s = cc%bsw*(1.0/sp%rho_wood - 1.0/rho_cellwall)*rho_H2O ! max stem water, kg H2O
    cc%Wmin_l = cc%Wmax_l * exp(sp%psi0_LF*sp%CR_Leaf)
    cc%Wmin_s = cc%Wmax_s * exp(sp%psi0_WD*sp%CR_Wood)
  end associate
end subroutine Update_plant_hydro_vars

!==============================================================================!
!  SUBROUTINE: PSI2RWC
!> \breif Convert water potential of leaf and wood to relative water content
!> \details Here we assume a constant hydraulic capacitance for both leaf and
!> wood. From the definition of hydraulic capacitance we have \n
!>       hydro_cap = delta_water_content / delta_psi \n
!> Since psi = 0. when water_content is at saturation, we have \n
!>       hydro_cap = (1. - rwc) * water_content_at_saturation / (0. - psi) \n
!> Reorganize the equation above, we can get \n
!>       rwc = 1. + psi * hydro_cap / water_content_at_saturation
!==============================================================================!
subroutine Plant_psi2water(cc)
  implicit none
  type(cohort_type), intent(inout) :: cc

  associate ( sp => spdata(cc%species) )
    !cc%W_leaf = cc%Wmax_L - (cc%Wmax_L - cc%Wmin_L) * cc%psi_leaf/sp%psi0_LF
    !cc%W_stem = cc%Wmax_S - (cc%Wmax_S - cc%Wmin_S) * cc%psi_stem/sp%psi0_WD
    cc%W_leaf = cc%Wmax_l * exp(cc%psi_leaf*sp%CR_Leaf)
    cc%W_stem = cc%Wmax_s * exp(cc%psi_stem*sp%CR_Wood)
  end associate
end subroutine Plant_psi2water

!===========================================================================!
!  SUBROUTINE: RWC2PSI
!> \brief Convert plant water content to water potential
!===========================================================================!
subroutine Plant_water2psi_exp(cc)
  implicit none
  type(cohort_type), intent(inout) :: cc
  !----local vars
  real :: W_status

  associate ( sp => spdata(cc%species))
    if(cc%Wmax_s > 1.0E-4)then
      W_status = MIN(max(1.0E-4,cc%W_stem),cc%Wmax_s)
      cc%psi_stem = log(W_status/cc%Wmax_s)/sp%CR_Wood
    else
      cc%psi_stem = cc%psi_s0
    endif

    if(cc%Wmax_l > 1.0E-4)then
      W_status = MIN(max(1.0E-4,cc%W_leaf),cc%Wmax_l)
      cc%psi_leaf = log(W_status/cc%Wmax_l)/sp%CR_Leaf
    else
      cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height)
    endif
  end associate
end subroutine Plant_water2psi_exp

!====================================
subroutine plant_psi_s0(vegn,cc,psi_s0)
  !@sum calculate stem base water potential when water flux is zero
  !@+   Weng, 03/20/2022
  implicit none
  type(vegn_tile_type), intent(in) :: vegn
  type(cohort_type),    intent(in) :: cc
  real,                 intent(out):: psi_s0
  !------local var -----------
  real :: RAI(soil_L)    ! root area index
  real :: k_rs(soil_L)   ! root-soil layer conductance
  real :: sumK,sumPK
  integer :: i

  ! Calculate plant water potential
  sumK  = 0.0
  sumPK = 0.0
  do i=1, soil_L
    ! RAI(i) = cc%ArootL(i)/cc%Acrown
    ! k_rs(i)= vegn%K_soil(i)*SQRT(RAI(i))/(3.14159*thksl(i))
    k_rs(i)= vegn%K_soil(i)*cc%ArootL(i) ! per tree
    sumK = sumK + k_rs(i)
    sumPK= sumPK+ k_rs(i) * vegn%psi_soil(i)
  enddo
  psi_s0 = sumPK / sumK

end subroutine plant_psi_s0

!==========================================
real function HT2MPa(ht) result(P) ! water gravational pressure
  real, intent(in):: ht ! Plant height, m
  !---------------------
  P = 1000.0 * 9.8 * ht * 1.0e-6 ! MPa
end function HT2MPa

!==========================================
real function HT2Lpath(ht) result(L) ! water path length from root to leaves
  real, intent(in):: ht ! Plant height, m
  !---------------------
  L = ht * 1.5 + 0.2 ! m
end function HT2Lpath

!==========================================
real function plc_function(psi, psi50, expK) result(plc)
  real, intent(in):: psi
  real, intent(in):: psi50
  real, intent(in):: expK
  !---------------------
  plc = 1.0 / (1.0 + (psi / psi50) ** expK)
end function plc_function

!============================================================================
real function NewWoodWTC(cc) result(WTC)
   !@sum calculate wood water transport capacity (WTC0), Ensheng Weng, 02/23/2022
    type(cohort_type),intent(in) :: cc

    !---------------------
    associate ( sp => spdata(cc%species))
    if (do_VariedWTC0)then
      WTC = sp%WTC0 * (1. + m0_WTC * cc%DBH**sp%thetaHT)
    else
      WTC = sp%WTC0 * (1. + m0_WTC * 0.5 ) ! 0.25**0.5 = 0.5
    endif
    end associate
end function NewWoodWTC

!============================================================================
real function NewWoodKx(cc) result(Kx)
    !@sum calculate wood water conductivity (Kx), Ensheng Weng, 02/23/2022
    type(cohort_type),intent(in) :: cc
    !---------------------
    associate ( sp => spdata(cc%species))
    ! Kx scaling with tree height to keep psi_leaf constant
    if (do_VariedKx)then
      Kx = sp%Kx0 * (1. + m0_kx * cc%DBH**sp%thetaHT)  ! (mm/s)/(MPa/m)
    else
      Kx = sp%Kx0 !* (1. + m0_kx * 0.5 ) ! 0.25**0.5 = 0.5
    endif
    end associate
end function NewWoodKx

!==============================================================
!============= Vegetation initializations =====================
subroutine initialize_vegn_tile(vegn)
   type(vegn_tile_type),intent(inout) :: vegn

   !--------local vars -------
   type(cohort_type),dimension(:), pointer :: cc => null()
   type(cohort_type),pointer :: cp
   real    :: btotal
   integer :: nCohorts = 1 ! Randomly generate n Cohorts if not defined
   integer :: i, istat

   ! Setup initial soil and vegetation conditions
   vegn%age = 0.0 ! Set the tile age as zero.
   call initialize_soil(vegn)    ! Initial Soil pools and environmental conditions
   call initialize_cohorts(vegn) ! Initialize plant cohorts
   call relayer_cohorts(vegn)    ! Sorting cohorts

   ! Setup an ID for each cohort
   do i=1,size(vegn%cohorts) ! vegn%n_cohorts
      cp => vegn%cohorts(i)
      cp%ccID = i
   enddo
   MaxCohortID = cp%ccID

   ! tile summary
   call vegn_sum_tile(vegn)
   vegn%initialN0 = totalN(vegn)
   vegn%totN      = vegn%initialN0

   ! Make a copy of the initial cohorts
   allocate(cc(1:init_cohort_N), STAT = istat)
   cc = vegn%cohorts
   vegn%initialCC   => cc
   vegn%n_initialCC = init_cohort_N
   cc => null()

end subroutine initialize_vegn_tile

!============================================================================
!Weng, 12/20/2022, Reset to Initial Vegetation States
subroutine reset_vegn_initial(vegn)
   type(vegn_tile_type),intent(inout) :: vegn

   !--------local vars -------
   type(cohort_type),dimension(:), pointer :: cc1
   type(cohort_type), pointer :: cp
   integer :: i

   ! Keep old cohorts
   cc1 => vegn%cohorts
   !Reset to initial plant cohorts
   call initialize_cohorts(vegn)
   ! Relayering and summary
   call relayer_cohorts(vegn)
   call vegn_sum_tile(vegn)
   ! ID each cohort
   do i=1, vegn%n_cohorts
      cp => vegn%cohorts(i)
      cp%ccID = MaxCohortID + i
   enddo
   MaxCohortID = cp%ccID
   !write(*,*)"Vegetaion resetted to initial conditions!"

   ! Release old cohorts
   deallocate(cc1)
end subroutine reset_vegn_initial

!================================================================
subroutine initialize_cohorts(vegn)
   type(vegn_tile_type),intent(inout) :: vegn

   !--------local vars -------
   type(cohort_type),dimension(:), pointer :: cc
   type(cohort_type),pointer :: cp
   real    :: btotal
   integer :: nCohorts = 1 ! Randomly generate n Cohorts if not defined
   integer :: i, istat

   ! Check initial species and density
   write(*,*)'Initial PFTs: ', init_cohort_sps(1:init_cohort_N)
   write(*,*)'Initial Dens: ', init_cohort_Indiv(1:init_cohort_N)

   ! Initialize plant cohorts
   allocate(cc(1:init_cohort_N), STAT = istat)
   vegn%cohorts => cc
   vegn%n_cohorts = init_cohort_N
   cc => null()
   do i=1,init_cohort_N
      cp => vegn%cohorts(i)
      cp%status  = LEAF_OFF ! ON=1, OFF=0 ! ON
      cp%layer   = 1
      cp%age     = 0
      cp%topyear = 0.0
      cp%species = init_cohort_sps(i)
      cp%nindivs = init_cohort_Indiv(i) ! trees/m2
      cp%nsc     = init_cohort_nsc(i)
      cp%bsw     = init_cohort_bsw(i)
      cp%bHW     = init_cohort_bHW(i)
      btotal     = cp%bsw + cp%bHW  ! kgC /tree
      call initialize_cohort_from_biomass(cp,btotal,maxval(vegn%psi_soil(:)))
   enddo
end subroutine initialize_cohorts

!============================================================
subroutine initialize_soil(vegn)
   type(vegn_tile_type),intent(inout) :: vegn

   ! Initial Soil pools and environmental conditions
   vegn%SOC(4)    = init_fast_SOC ! kgC m-2
   vegn%SOC(5)    = init_slow_SOC ! slow soil carbon pool, (kg C/m2)
   vegn%SON(4)    = vegn%SOC(4)/CN0SOM(4)  ! fast soil nitrogen pool, (kg N/m2)
   vegn%SON(5)    = vegn%SOC(5)/CN0SOM(5)  ! slow soil nitrogen pool, (kg N/m2)
   vegn%N_input   = N_input  ! kgN m-2 yr-1, N input to soil
   vegn%mineralN  = init_mineralN  ! Mineral nitrogen pool, (kg N/m2)
   vegn%previousN = vegn%mineralN
   !Soil water
   vegn%soiltype = soiltype
   vegn%FLDCAP   = soilpars(soiltype)%vwc_fc
   vegn%WILTPT   = soilpars(soiltype)%vwc_wilt
   vegn%wcl      = soilpars(soiltype)%vwc_sat
   vegn%soilZ    = sum(thksl(:))
   vegn%soilWP0  = vegn%WILTPT * vegn%soilZ * 1000.0 ! Minimum soil water
   vegn%thetaS   = 1.0
   call SoilWater_psi_K(vegn)
end subroutine initialize_soil

! ============================================================================
! Initialize a cohort by initial biomass and soil water conditions
subroutine initialize_cohort_from_biomass(cc,btot,psi_s0)
  type(cohort_type), intent(inout) :: cc
  real,intent(in)    :: btot ! total biomass per individual, kg C
  real,intent(in)    :: psi_s0 ! Initial stem water potential

  !---- local vars ------------
  integer :: j

  associate(sp=>spdata(cc%species))
    call BM2Architecture(cc,btot)
    call update_max_LFR_NSN(cc)
    if(sp%leafLS>1.0)then
      cc%leafage = 0.5 ! 1.0
      cc%bl      = cc%bl_max
    else
      cc%leafage = 0.0
      cc%bl      = 0.0
    endif
    cc%br     = cc%br_max
    cc%nsc    = 2.0 * (cc%bl_max + cc%br_max)
    cc%seedC  = 0.0

    ! N pools
    cc%NSN   = 5.0*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
    cc%leafN = cc%bl/sp%CNleaf0
    cc%rootN = cc%br/sp%CNroot0
    cc%swN   = cc%bsw/sp%CNwood0
    cc%hwN   = cc%bHW/sp%CNwood0
    cc%seedN = 0.0

    ! ----------Plant hydraulics----------
    cc%Nrings  = 1
    cc%WTC0(1) = NewWoodWTC(cc)
    cc%Kx(1)   = NewWoodKx(cc)
    cc%Rring(1)= cc%DBH/2.0
    cc%Lring(1)= HT2Lpath(cc%height) !Max(cc%height,0.01)
    cc%Aring(1)= PI * cc%Rring(1)**2
    cc%Ktrunk  = PI * cc%Aring(1) * cc%Kx(1)/cc%Lring(1)
    cc%treeHU = 0.0
    cc%treeW0 = cc%WTC0(1) * cc%Aring(1)
    cc%farea(1)= 1.0
    cc%accH(1) = 0.0
    cc%plcH(1) = 0.0
    ! Initial psi
    cc%psi_stem = psi_s0
    cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height) ! MPa
    call Update_plant_hydro_vars(cc)
    call plant_psi2water(cc)

  end associate
end subroutine initialize_cohort_from_biomass

!=======================================================================
!==================== Cohort management ================================
!=======================================================================
subroutine relayer_cohorts (vegn)
  ! Arrange crowns into canopy layers according to PPA
  type(vegn_tile_type), intent(inout) :: vegn ! input cohorts

  ! ---- local constants
  real, parameter :: tolerance = 1e-4
  real, parameter :: layer_vegn_cover = 1.0
  ! ---- local vars
  integer :: idx(vegn%n_cohorts) ! indices of cohorts in decreasing height order
  integer :: i ! new cohort index
  integer :: k ! old cohort index
  integer :: L ! layer index (top-down)
  integer :: N0,N1 ! initial and final number of cohorts
  real    :: frac ! fraction of the layer covered so far by the canopies
  type(cohort_type), pointer :: cc(:),new(:)
  real    :: nindivs

  ! rank cohorts in descending order by height
  N0 = vegn%n_cohorts
  cc => vegn%cohorts
  call rank_descending(cc(1:N0)%height,idx)

  ! calculate max possible number of new cohorts : it is equal to the number of
  ! old cohorts, plus the number of layers -- since the number of full layers is
  ! equal to the maximum number of cohorts split by the layer boundaries
  N1 = vegn%n_cohorts + int(sum(cc(1:N0)%nindivs*cc(1:N0)%Acrown))
  allocate(new(N1))

  ! copy cohort information to the new cohorts, splitting the old cohorts that
  ! stride the layer boundaries
  i = 1
  k = 1
  L = 1
  frac = 0.0
  nindivs = cc(idx(k))%nindivs
  do
     new(i)         = cc(idx(k))
     new(i)%nindivs = min(nindivs,(layer_vegn_cover-frac)/cc(idx(k))%Acrown)
     new(i)%layer   = L
     frac = frac+new(i)%nindivs*new(i)%Acrown
     nindivs = nindivs - new(i)%nindivs

     if (abs(nindivs*cc(idx(k))%Acrown)<tolerance) then
       new(i)%nindivs = new(i)%nindivs + nindivs ! allocate the remainder of individuals to the last cohort
       if (k==N0) exit ! end of loop
       k = k+1
       nindivs = cc(idx(k))%nindivs  ! go to the next input cohort
     endif

     if (abs(layer_vegn_cover - frac)<tolerance) then
       L = L+1
       frac = 0.0              ! start new layer
     endif
     i = i+1
  enddo

  ! replace the array of cohorts
  vegn%cohorts => new
  vegn%n_cohorts = i

  ! update layer fraction for each cohort
  do i=1, vegn%n_cohorts
     vegn%cohorts(i)%layerfrac = vegn%cohorts(i)%nindivs * vegn%cohorts(i)%Acrown
  enddo

  ! Release old cohort array
  deallocate(cc)
  cc  => null()
  new => null()
end subroutine relayer_cohorts

!============================================================================
! Merge similar cohorts in a tile
! Only merge cohorts within current cohort array, Weng, 06/21/2025
subroutine vegn_mergecohorts(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  integer :: i, j

  merged(:) = .FALSE.
  do i = 1, vegn%n_cohorts
     if(merged(i)) cycle ! skip cohorts that were already merged
     ! Go through rest cohorts for merging them into the current one
     do j = i+1, vegn%n_cohorts
        if (Mergeable_cohorts(vegn%cohorts(j),vegn%cohorts(i))) then
           call merge_cohorts(vegn%cohorts(j),vegn%cohorts(i))
           merged(j) = .TRUE.
        endif
     enddo
  enddo

end subroutine vegn_mergecohorts

! ============================================================================
! kill low density cohorts, a new function seperated from vegn_mergecohorts
! Weng, 2014-07-22
! If MergeLowDenCohorts is True, this routine only remove zero density
! cohorts and the low density cohorts that have no other cohorts with the
! same PFT in a layer. Weng, 06/25/2025

subroutine kill_lowdensity_cohorts(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cp, cc(:) ! array to hold new cohorts
  integer :: i,j,k

 ! calculate the number of cohorts with indivs>min_nindivs
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs > min_nindivs) k=k+1
  enddo

  ! exclude cohorts that have low individuals
  if (k>0 .and. k<vegn%n_cohorts)then
     allocate(cc(k))
     j=0
     do i = 1,vegn%n_cohorts
        cp =>vegn%cohorts(i)
        if (cp%nindivs > min_nindivs) then
           j=j+1
           cc(j) = cp
        else
           ! Carbon and Nitrogen from plants to soil pools
           call plant2soil(vegn,cp,cp%nindivs)
        endif
     enddo
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc
     vegn%n_cohorts = k
  endif
end subroutine kill_lowdensity_cohorts

! ============================================================================
! kill old grass cohorts
! Weng, 01/22/2023
subroutine kill_old_grass(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cp, cc(:) ! array to hold new cohorts
  !logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  logical :: OldGrass
  integer :: i,j,k

 ! calculate the number of cohorts that are not old grass
  k = 0
  do i = 1, vegn%n_cohorts
    cp =>vegn%cohorts(i)
    associate(sp=>spdata(cp%species))
      OldGrass = (sp%lifeform ==0 .and. cp%age > MaxGrassAge)
      if (.not. OldGrass) k=k+1
    end associate
  enddo
  if (k==0)then
     write(*,*)'in kill_old_grass: All cohorts are old grass, No action!'
     !stop
  endif

  ! exclude cohorts that are old grass
  if (k>0 .and. k<vegn%n_cohorts)then
     allocate(cc(k))
     j=0
     do i = 1,vegn%n_cohorts
        cp =>vegn%cohorts(i)
        associate(sp=>spdata(cp%species))
        OldGrass = (sp%lifeform ==0 .and. cp%age > MaxGrassAge)
        if (.not. OldGrass) then
           j=j+1
           cc(j) = cp
        else
           ! Carbon and Nitrogen from plants to soil pools
           call plant2soil(vegn,cp,cp%nindivs)
        endif
        end associate
     enddo
     vegn%n_cohorts = j
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc
  endif
end subroutine kill_old_grass

! ============================================================================
! Just merge cohorts c1 into c2, regardless of their sizes and densities
subroutine merge_cohorts(c1,c2) ! Put c1 into c2
  type(cohort_type), intent(inout) :: c1
  type(cohort_type), intent(inout) :: c2

  real :: x1, x2 ! normalized relative weights
  real :: Ntot, btot, c2H, c2DBH
  integer :: i

     ! Set up merging weights (x1 and x2)
     Ntot = c1%nindivs + c2%nindivs
     x1 = c1%nindivs/Ntot
     x2 = c2%nindivs/Ntot
     c2DBH = c2%dbh
     c2H   = c2%height

     ! Update c2 (merging c1 into c2)
     c2%ccID = c1%ccID
     c2%nindivs = Ntot
     c2%age = x1*c1%age + x2*c2%age
     c2%topyear = x1*c1%topyear + x2*c2%topyear

     !  Carbon
     c2%bl  = x1*c1%bl  + x2*c2%bl
     c2%br  = x1*c1%br  + x2*c2%br
     c2%bsw = x1*c1%bsw + x2*c2%bsw
     c2%bHW = x1*c1%bHW + x2*c2%bHW
     c2%seedC = x1*c1%seedC + x2*c2%seedC
     c2%nsc = x1*c1%nsc + x2*c2%nsc

     !  Nitrogen
     c2%leafN = x1*c1%leafN + x2*c2%leafN
     c2%rootN = x1*c1%rootN + x2*c2%rootN
     c2%swN   = x1*c1%swN   + x2*c2%swN
     c2%hwN   = x1*c1%hwN   + x2*c2%hwN
     c2%seedN = x1*c1%seedN + x2*c2%seedN
     c2%NSN   = x1*c1%NSN   + x2*c2%NSN

     ! Water content
     c2%W_leaf  = x1*c1%W_leaf  + x2*c2%W_leaf
     c2%W_stem  = x1*c1%W_stem  + x2*c2%W_stem
     c2%W_dead  = x1*c1%W_dead  + x2*c2%W_dead

     ! Allometry recalculation
     btot = c2%bsw + c2%bHW
     call BM2Architecture(c2,btot)
     call update_max_LFR_NSN(c2)
     call Update_plant_hydro_vars(c2)
     call Plant_water2psi_exp(c2)

     ! Reset tree rings' hydraulics
     c2%Kx    = x1*c1%Kx    + x2*c2%Kx
     c2%WTC0  = x1*c1%WTC0  + x2*c2%WTC0
     c2%accH  = x1*c1%accH  + x2*c2%accH
     c2%farea = x1*c1%farea + x2*c2%farea

     ! Recalculate ring variables
     c2%Rring   = c2%Rring * c2%dbh/c2DBH
     c2%Lring   = c2%Lring * c2%height/c2H
     do i=2,MIN(c2%Nrings, Ysw_max)
       if (c2%Rring(i)>c2%Rring(i-1)) then
         c2%Aring(i) = PI*(c2%Rring(i)**2 - c2%Rring(i-1)**2)
       endif
     enddo

     ! Update tree trunk sapwood area and conductance
     call calculate_Asap_Ktrunk (c2)

     ! Zero c1's desnity so that it can be removed by kill_lowdensity_cohorts
     c1%nindivs = 0.0

end subroutine merge_cohorts

! ============================================================================
function Mergeable_cohorts(c1,c2); logical Mergeable_cohorts
   type(cohort_type), intent(in) :: c1,c2
   logical :: sameSpecies,sameLayer,sameSize,sameSizeTree,sameSizeGrass
   logical :: LowDensity, Not_ZeroDensity

   LowDensity = .False. ! Default

   ! Mergeable criteria:
   Not_ZeroDensity = c1%nindivs > 0.0 .and. c2%nindivs > 0.0
   sameSpecies  = c1%species == c2%species
   sameLayer    = (c1%layer == c2%layer) .or. &
                  ((spdata(c1%species)%lifeform ==0) .and. &
                   (spdata(c2%species)%lifeform ==0) .and. &
                   (c1%layer>1 .and.c2%layer>1))
   sameSizeTree = (spdata(c1%species)%lifeform > 0).and.  &
                  (spdata(c2%species)%lifeform > 0).and.  &
                 ((abs(c1%DBH - c2%DBH)/c2%DBH < diff_S0 ) .or.  &
                  (abs(c1%DBH - c2%DBH)        < 0.001))  ! it'll be always true for grasses
   sameSizeGrass= (spdata(c1%species)%lifeform ==0) .and. &
                  (spdata(c2%species)%lifeform ==0) .and. &
                  (c1%DBH == c2%DBH)  ! it'll be always true for grasses
   sameSize = sameSizeTree .OR. sameSizeGrass

   if (MergeLowDenCohorts) &
       LowDensity = c1%nindivs < min_nindivs .OR. c2%nindivs < min_nindivs

   Mergeable_cohorts = sameSpecies .and. sameLayer &
                      .and. Not_ZeroDensity       &
                      .and. (sameSize .or. LowDensity)
end function

! ============================================================================
subroutine check_N_conservation(vegn,totN0,tag)
  type(vegn_tile_type),intent(in) :: vegn
  real,                intent(in) :: totN0
  character(len = *),  intent(in) :: tag
  !-------local var --------
  real :: totN1
  ! Total N balance checking
  totN1 = TotalN(vegn)
  if(abs(totN0 - totN1) > 1.0E-6)then ! Precision: 1.19209290E-07
    write(*,*)"Imbalance of nitrogen in: ", tag
    write(*,*)'N0, N1, N0-N1', totN0, totN1, totN0 - totN1
    !stop
  endif
  write(*,*)tag, ': N0, N1, N0-N1', totN0, totN1, totN0 - totN1
end subroutine check_N_conservation

! ============================================================================
function TotalN(vegn)
  real :: TotalN ! returned value
  type(vegn_tile_type), intent(in) :: vegn

  TotalN = vegn%NSN + vegn%SeedN + vegn%leafN + vegn%rootN + &
           vegn%SwN + vegn%HwN + vegn%mineralN + sum(vegn%SON(:))
end function

!======================= Specific experiments ================================

!----------------------- fire disturbance (Konza) ---------------------------
subroutine vegn_fire (vegn, deltat)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! seconds since last mortality calculations, s

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  type(spec_data_type), pointer :: sp
  ! -- fire effects variables --
  real :: P_ET    ! ratio of P to ET, for fire probability calculation
  real :: Frisk   ! Environmental fire risk (a function of P_ET)
  real :: P_Ign, r_Ign  ! Probability of ignition, and a random number for fire occurence
  real :: Ign_G0, Ign_W0 ! Grass and wood ignition probability
  real :: f_grass, f_wood  ! grasses and canopy tree spread probabilities
  real :: flmb_G, flmb_W ! Flamability of grasses and woody plants
  real :: d_tree ! Tree's sensitivity to ground surface fire: max 1.0, min 0.0
  real :: s_fireG ! Grass fire severity
  real :: mu_fire ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  real :: Cfire, Cfast, Cslow, Nfire, Nfast, Nslow ! C and N fates at fires
  integer :: i, k

  !  Parameters (defined in datatypes.F90 and read in from the namelist file):
  !  Frisk: Environmental fire occurrence probability, a function of environmental
  !  conditions that can result in fire if fuel is available
  !  (i.e., (match-dropping probability). It should be function of environmental conditions
  !  Vegetation flammability parameters, Ign_G0, Ign_W0:
  !  Ignition probability for grasses and woody plants once environmental conditions meet Frisk
  !  For grasses: Ign_G0 = 1.0; For woody plants: Ign_W0 = 0.025
  !  mu0_FireW, mu0_FireG: mortality rates of trees and grasses due to fire
  !  r_BK0: shape parameter, for bark resistance, exponential equation,
  !                                  120 --> 0.006 m of bark 0.5 survival

  ! Environmental risk
  if(Do_FixedFrisk) then
    Frisk = EnvF0 ! Fixed environment risk
  else
    P_ET = vegn%annualPrcp / vegn%annualPET
    Frisk = 1.0/(1.0 + exp(A_MI*(P_ET - MI0Fire)))
  endif

  ! Ignition probabilities of grasses and woody PFTs
  Ign_G0 = 0.0
  Ign_W0 = 0.0
  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    associate ( sp => spdata(cc%species))
      if(sp%lifeform==0) then  ! grasses
        Ign_G0 = max(Ign_G0, sp%IgniteP) ! Use the max IgniteP
      else                     ! trees
        Ign_W0 = max(Ign_W0, sp%IgniteP)
      endif
    end associate
  enddo

  ! Burning probability
  f_grass = min(1.0, vegn%GrassCA)
  f_wood  = min(1.0, vegn%TreeCA)
  flmb_G  = Ign_G0 * f_grass
  flmb_W  = Ign_W0 * f_wood
  P_Ign   = 1.0 - (1.- flmb_G * Frisk)*(1. - flmb_W * Frisk)

  ! fire effects on vegetation and soil
  CALL RANDOM_NUMBER(r_Ign) ! r_Ign = rand(0)
  Cfire = 0.0; Cfast = 0.0; Cslow = 0.0
  Nfire = 0.0; Nfast = 0.0; Nslow = 0.0
  if(r_Ign < P_Ign)then ! Fire_ON
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      associate ( sp => spdata(cc%species))
      if(sp%lifeform==0) then  ! grasses
         mu_fire = mu0_FireG
      else                     ! trees
         if(r_Ign < flmb_W * Frisk) then   ! tree canopy fire
            mu_fire = 0.99 * min(1.0, 1.25 * f_wood)
         else                                     ! grass fire
            ! 05/01/2024 (Kelvin), s_fireG should be a function of grass biomass.
            ! At high grass biomass, the fire has higher serverity and kills
            ! more shrubs. So, extreme droughts can trigger expansion of shrubs.
            ! s_fireG is grass fire severity, defined as a function of grass biomass
            s_fireG = max(0.0,min(1.0,(vegn%GrassBM/FSBM0)**2)) ! grass fire severity
            cc%D_bark = f_bk * cc%dbh    ! bark thickness,
            d_tree = exp(r_BK0*cc%D_bark)! Tree's fire sensitivity to grass fire
            !d_tree = 1. - cc%D_bark/(cc%D_bark+D_BK0) !Alternative formulation
            mu_fire = mu0_FireW * d_tree * f_grass * s_fireG
         endif
      endif

      ! Burned vegetation and soils
      deadtrees = cc%nindivs * MIN(1.0,mu_fire * deltat/seconds_per_year) ! individuals / m2

      ! Carbon and Nitrogen release by burning
      Cfire = Cfire + (0.2*cc%NSC + 0.7*cc%bl    + 0.2*(cc%bsw+cc%bHW) + 0.0*cc%br    + 0.0*cc%seedC) * deadtrees
      Cfast = Cfast + (0.8*cc%NSC + 0.3*cc%bl    + 0.0*(cc%bsw+cc%bHW) + 1.0*cc%br    + 1.0*cc%seedC) * deadtrees
      Cslow = Cfast + (0.0*cc%NSC + 0.0*cc%bl    + 0.8*(cc%bsw+cc%bHW) + 0.0*cc%br    + 0.0*cc%seedC) * deadtrees
      Nfire = Nfire + (0.2*cc%NSN + 0.7*cc%leafN + 0.2*(cc%swN+cc%hwN) + 0.0*cc%rootN + 0.0*cc%seedN) * deadtrees
      Nfast = Nfast + (0.8*cc%NSN + 0.3*cc%leafN + 0.0*(cc%swN+cc%hwN) + 1.0*cc%rootN + 1.0*cc%seedN) * deadtrees
      Nslow = Nslow + (0.0*cc%NSN + 0.0*cc%leafN + 0.8*(cc%swN+cc%hwN) + 0.0*cc%rootN + 0.0*cc%seedN) * deadtrees

      ! Update plant density
      cc%nindivs = cc%nindivs - deadtrees
      end associate
    enddo
    ! C and N fluxes due to fire
    vegn%C_burned = Cfire + 0.7*vegn%SOC(1)+0.2*vegn%SOC(2) ! Burned litter: 70% of fine litter and 20% of coarse litter are burned
    vegn%mineralN = Nfire + 0.7*vegn%SON(1)+0.2*vegn%SON(2) + vegn%mineralN
    vegn%Nm_Fire  = Nfire + 0.7*vegn%SON(1)+0.2*vegn%SON(2) + vegn%Nm_Fire

    ! Update Litter C and N pools
    vegn%SOC(1) = (1.0-0.7)*vegn%SOC(1) + Cfast
    vegn%SOC(2) = (1.0-0.2)*vegn%SOC(2) + Cslow
    vegn%SON(1) = (1.0-0.7)*vegn%SON(1) + Nfast
    vegn%SON(2) = (1.0-0.2)*vegn%SON(2) + Nslow

    ! Annual N from plants to soil
    vegn%NorgP2S = vegn%NorgP2S + Cfast + Cslow

#ifdef ScreenOutput
    write(*,*)"fire, TreeCA, GrassCA", &
        r_Ign < P_Ign, vegn%TreeCA, vegn%GrassCA
#endif
  endif

  ! Record Frisk and Pfire for output
  vegn%Frisk = Frisk
  vegn%Pfire = P_Ign

end subroutine vegn_fire

!=======================================================================
! switch the species of the first cohort to another species
! bugs !!!!!!
 subroutine vegn_species_switch(vegn,N_SP,iyears,FREQ)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in):: N_SP  ! total species in model run settings
  integer, intent(in):: iyears
  integer, intent(in):: FREQ  ! frequency of species switching

  ! ---- local vars --------
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse
  integer :: i, k
  type(cohort_type), pointer :: cc

     cc => vegn%cohorts(1)
     associate (sp => spdata(cc%species)) ! F2003
     if(cc%bl > 0.0) then ! remove all leaves to keep mass balance
        loss_coarse  = cc%nindivs * (cc%bl - cc%Aleaf*LMAmin)
        loss_fine    = cc%nindivs *  cc%Aleaf*LMAmin
        lossN_coarse = cc%nindivs * (cc%leafN - cc%Aleaf*sp%LNbase)
        lossN_fine   = cc%nindivs *  cc%Aleaf*sp%LNbase
        ! Carbon to soil pools
        vegn%SOC(1) = vegn%SOC(1) + fsc_fine *loss_fine + &
                fsc_wood *loss_coarse
        vegn%SOC(2) = vegn%SOC(2) + (1.0-fsc_fine)*loss_fine + &
                (1.0-fsc_wood)*loss_coarse
        ! Nitrogen to soil pools
        vegn%SON(1) = vegn%SON(1) + fsc_fine  *lossN_fine +   &
                                        fsc_wood *lossN_coarse
        vegn%SON(2) = vegn%SON(2) +(1.-fsc_fine) *lossN_fine +   &
                                      (1.-fsc_wood)*lossN_coarse
        ! annual N from plants to soil
        vegn%NorgP2S = vegn%NorgP2S + lossN_fine + lossN_coarse
        ! remove leaves
        cc%bl = 0.0
     endif
     end associate
     ! Change species
     cc%species = mod(iyears/FREQ,N_SP)+2

 end subroutine vegn_species_switch

!=======================================================================
! Put missing PFTs back from the initial cohorts (initialCC)
! 10/17/2018, Weng
subroutine vegn_species_recovery (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc ! parent and child cohort pointers
  type(cohort_type), dimension(:),pointer :: ccold, ccnew   ! pointer to old cohort array
  integer,dimension(16) :: initialPFTs, currPFTs, missingPFTs
  real :: addedC, addedN ! seed pool of productible PFTs
  real :: totC, totN
  integer :: newcohorts, matchflag
  integer :: nPFT0 ! number of initial PFTs
  integer :: nPFTs ! number of new cohorts to be created
  integer :: nCohorts, istat
  integer :: i, j, k, n ! cohort indices

  ! Looping through all current cohorts and get the exsiting PFTs
  currPFTs = -999 ! the code of existing PFT
  nPFTs = 0
  do k=1, vegn%n_cohorts
     cc => vegn%cohorts(k)
     matchflag = 0
     do i=1,nPFTs
        if(cc%species == currPFTs(i))then
            matchflag = 1
            exit
        endif
     enddo
     if(matchflag==0)then ! when it is a new PFT, put it to the next place
        nPFTs            = nPFTs + 1 ! update the number of existing PFTs
        currPFTs(nPFTs) = cc%species ! PFT number
     endif
  enddo ! k, vegn%n_cohorts
  !write(*,'(A12,16(I4,","))')'current cc',currPFTs(1:nPFTs)

  ! Looping through all initial cohorts and get the initial PFTs
  initialPFTs = -999 ! the code of initial PFT
  nPFT0 = 0
  do k=1, vegn%n_initialCC
     cc => vegn%initialCC(k)
     !initialPFTs(k) = cc%species
     !nPFT0 = nPFT0 + 1
     matchflag = 0
     do i=1,nPFT0
        if(cc%species == initialPFTs(i))then
            matchflag = 1
            exit
        endif
     enddo
     if(matchflag==0)then ! when it is a new PFT, put it to the next place
        nPFT0      = nPFT0 + 1 ! update the number of existing PFTs
        initialPFTs(nPFT0) = cc%species ! PFT number
     endif
  enddo ! k, vegn%n_initialCC
  !write(*,'(A12,16(I4,","))')'initial cc',initialPFTs(1:nPFT0)

  ! Get missing PFT(s)
  missingPFTs = -999
  n = 0
  if(nPFTs < nPFT0)then
      do i=1, nPFT0
         matchflag = 0
         do j=1, nPFTs ! go through all current PFTs
            if(initialPFTs(i)==currPFTs(j))then
               matchflag = 1
               exit
            endif
         enddo
         if(matchflag ==0)then
             n = n + 1
             missingPFTs(n) = initialPFTs(i)
         endif
      enddo
      !write(*,*)'missing PFTs',missingPFTs(1)
  endif

  ! Generate new cohorts
  addedC = 0.0
  addedN = 0.0
  newcohorts = n
  if (newcohorts >= 1) then   ! build new cohorts for seedlings
     ccold => vegn%cohorts ! keep old cohort information
     nCohorts = vegn%n_cohorts + newcohorts
     allocate(ccnew(1:nCohorts), STAT = istat)
     ccnew(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
     vegn%cohorts => ccnew

     deallocate (ccold)

     ! set up new cohorts
     k = vegn%n_cohorts
     do i = 1, newcohorts
        k = k+1 ! increment new cohort index
        cc => vegn%cohorts(k)
        ! Give the new cohort an ID
        cc%ccID = MaxCohortID + i
        ! find out the equivalent initial cohort with species n
        do n=1, vegn%n_initialCC
            if(vegn%initialCC(n)%species == missingPFTs(i))then
               ! update new cohort parameters
               cc = vegn%initialCC(n) ! Reset the missing PFT to the initial
               !cc%nindivs = 0.01 ! a small number of individuals
               exit
            endif
        enddo
        totC = cc%bl + cc%br + cc%bsw + cc%bHW + cc%nsc
        totN = cc%leafN + cc%rootN + cc%swN + cc%hwN + cc%NSN
        addedC = addedC + cc%nindivs * totC
        addedN = addedN + cc%nindivs * totN
        call setup_seedling(cc,totC,totN)
     enddo

     MaxCohortID = MaxCohortID + newcohorts
     vegn%n_cohorts = k
     ccnew => null()

     ! Make carbon and nitrogen balance
     vegn%SOC(2) = vegn%SOC(2) - min(0.99*vegn%SOC(2),addedC)
     vegn%SON(2) = vegn%SON(2) - min(0.99*vegn%SON(2),addedN)
  endif ! set up newly moved-in cohorts

end subroutine vegn_species_recovery

! =============================================================================
! Added by Weng 2015-02-29
subroutine vegn_annualLAImax_update(vegn)
  ! used for updating LAImax according to mineral N in soil
  ! Potential problems:
  !   1. All species LAImax are updated
  !   2. For evergreen, LAImax can be less than current LAI.
  !  Weng, 2017-08-02
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc
  real   :: LAImin, LAIfixedN, LAImineralN
  real   :: LAI_Nitrogen
  real   :: fixedN, rootN
  logical:: fixedN_based
  integer :: i
  ! Calculating LAI max based on mineral N or mineralN + fixed N
  fixedN_based =  .False. ! .True. !
  LAImin       = 0.5

  !fixedN = 0.0
  !do i = 1,vegn%n_cohorts
  !      cc => vegn%cohorts(i)
  !      fixedN = fixedN + cc%NfixedYr * cc%Acrown * cc%nindivs
  !enddo

  ! Mineral+fixed N-based LAImax
  ! LAI_fixedN = sp%R0_Nfix * sp%LMA * sp%CNleaf0 * sp%leafLS / sp%LMA
  ! cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
  vegn%previousN = 0.8 * vegn%previousN + 0.2 * vegn%Nm_Soil
  do i=0,MSPECIES
      associate (sp => spdata(i) )

      LAIfixedN  = 0.5 * sp%R0_Nfix * sp%CNleaf0 * sp%leafLS
      LAImineralN = 0.5*vegn%previousN*sp%CNleaf0*sp%leafLS/sp%LMA
      !LAImineralN = vegn%previousN/(sp%LMA/(sp%CNleaf0*sp%leafLS)+sp%phiRL*sp%alpha_FR/sp%SRA /sp%CNroot0)
      LAI_nitrogen = LAIfixedN + LAImineralN

      spdata(i)%LAImax = MAX(LAImin, MIN(LAI_nitrogen,sp%LAI_light))
      spdata(i)%LAImax_u = MIN(sp%LAImax,1.2)
      end associate
  enddo

  !  ! update the PFTs in the first layer based on fixed N
  !  if(fixedN_based)then ! based on "cc%NfixedYr + vegn%previousN"
  !!    Reset sp%LAImax
  !     do i = 1,vegn%n_cohorts
  !        cc => vegn%cohorts(i)
  !        associate (sp => spdata(cc%species) )
  !        sp%LAImax    = 0.0  ! max(sp%LAImax,ccLAImax)
  !        sp%n_cc      = 0
  !        end associate
  !     enddo
  !!   Sum ccLAImax in the first layer
  !     do i = 1,vegn%n_cohorts
  !        cc => vegn%cohorts(i)
  !        associate ( sp => spdata(cc%species) )
  !        if(sp%LAImax < LAImin)then
  !           LAI_nitrogen = 0.5*(vegn%previousN+cc%NfixedYr)*sp%CNleaf0*sp%leafLS/sp%LMA
  !           if(sp%R0_Nfix > 0.0)
  !           sp%LAImax    = MAX(LAImin, MIN(LAI_nitrogen,sp%LAI_light))
  !        endif
  !        end associate
  !     enddo
  !  endif

end subroutine vegn_annualLAImax_update

!==============================================================
!Weng, 06-13-2021, Crown gap with biodiversity
subroutine vegn_gap_fraction_update(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real :: totCA, CA(0:MSPECIES), f_keep
  integer :: i

  ! Calculate totCA and each PFT's CA
  totCA     = 0.0
  CA(:)     = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     totCA = totCA + cc%Acrown * cc%nindivs
     CA(cc%species) = CA(cc%species) + cc%Acrown * cc%nindivs
  enddo
  !Calculate f_keep
  f_keep = 0.0
  do i=0, MSPECIES-1
    if(CA(i)>0.0)f_keep = f_keep + (CA(i)/totCA)**3
  enddo
  ! update f_cGap
  do i=0, MSPECIES-1
     if(CA(i)>0.0) spdata(i)%f_cGap = f_cGap(i) * f_keep
  enddo
  write(*,*)'PFT1_gap,PFT2_gap,f_keep',spdata(1)%f_cGap,spdata(2)%f_cGap,f_keep

end subroutine vegn_gap_fraction_update

!========================================================================
!=====================Plant hydraulics testing codes, not used =============

! Weng 2022-02-16 ! Compute water flux via tree trunk (i.e., soil-trunk-leaves)
subroutine Plant_water_dynamics_equi(vegn) ! forcing,
  !type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real :: psi_ht            ! Gravitational water pressure, MPa
  real :: psi_leaf,psi_stem
  real :: psi_soil(soil_L),K_soil(soil_L)
  real :: thetaS(soil_L) ! soil moisture index (0~1)
  real :: dpsiSR(soil_L) ! pressure difference between soil and root, MPa
  real :: k_stem          ! actual conductance of the whole tree
  real :: dW_L,dW_S  ! water demand of leaves and stems  ! mm/s
  integer :: i,j, layer

  ! Equilibrium plant water potentials
  do j = 1, vegn%n_cohorts
     cc => vegn%cohorts(j)
     associate ( sp => spdata(cc%species) )

     ! Calculate plant tissue water potentials ! Hack
     dW_L = Max(0.0, cc%transp + cc%Wmax_L - cc%W_leaf)/step_seconds ! per step -> per second
     dW_S = Max(0.0, cc%Wmax_S - cc%W_stem)/step_seconds ! per step -> per second
     call plant_water_potential_equi(vegn,cc,dW_L,dW_S,psi_leaf,psi_stem)
     cc%psi_stem = psi_stem
     cc%psi_leaf = psi_leaf

     ! Stem-Leaf water flux
     k_stem = cc%Ktrunk * plc_function(cc%psi_stem,sp%psi50_WD,sp%Kexp_WD)
     psi_ht = HT2MPa(cc%height) ! MPa
     cc%Q_leaf = (cc%psi_stem - cc%psi_leaf - psi_ht) * k_stem * step_seconds
     cc%W_leaf = cc%W_leaf - cc%transp + cc%Q_leaf
     cc%W_stem = cc%W_stem - cc%Q_leaf

     ! Water fluxes from soil layers to stem base (and between soil layers via roots)
     do i=1, soil_L ! Calculate water uptake potential layer by layer
        dpsiSR(i) = Max(0.0, vegn%psi_soil(i) - cc%psi_stem)
        cc%Q_soil(i) = cc%ArootL(i) * vegn%K_soil(i) * dpsiSR(i) * step_seconds
        cc%Q_soil(i) = Max(0.0, Min(cc%Q_soil(i), 0.05 * vegn%freewater(i)*cc%Acrown))
        ! Water uptaken by roots, hourly
        vegn%wcl(i) = Max(vegn%wcl(i) - (cc%Q_soil(i)*cc%nindivs/(thksl(i)*1000.0)),vegn%WILTPT)
        cc%W_stem = cc%W_stem + cc%Q_soil(i)
     enddo

     ! Water supply for the next step photosynthesis
     cc%W_supply = 0.05 * vegn%soilwater * cc%Acrown
     end associate
   enddo

end subroutine Plant_water_dynamics_equi
!====================================
subroutine plant_water_potential_equi(vegn,cc,Q_leaf,Q_stem,psi_leaf,psi_stem)
  !@sum leaf and stem equilibrium water potential
  !@+   Weng, 03/20/2022
  implicit none
  type(vegn_tile_type), intent(in) :: vegn
  type(cohort_type),    intent(in) :: cc
  real,                 intent(in) :: Q_leaf, Q_stem    ! mm/s/tree
  real,                 intent(out):: psi_leaf, psi_stem

  !------local var -----------
  real :: psi_ht          ! Gravitational water pressure, MPa
  real :: RAI(soil_L)    ! root area index
  real :: k_rs(soil_L)   ! root-soil layer conductance
  real :: sumK,sumPK
  real :: k_stem         ! with water conditions
  integer :: i

  !! Water supply from each soil layer
  associate ( sp => spdata(cc%species) )
    sumK  = 0.0
    sumPK = 0.0
    do i=1, soil_L
      ! RAI(i) = cc%ArootL(i)/cc%Acrown
      ! k_rs(i)= vegn%K_soil(i)*SQRT(RAI(i))/(3.14159*thksl(i))
      k_rs(i)= vegn%K_soil(i)*cc%ArootL(i) ! per tree
      sumK = sumK + k_rs(i)
      sumPK= sumPK+ k_rs(i) * vegn%psi_soil(i)
    enddo
    psi_ht = HT2MPa(cc%height) ! MPa
    k_stem = cc%Ktrunk * plc_function(cc%psi_stem,sp%psi50_WD,sp%Kexp_WD)
    psi_stem = Max(sp%psi0_WD, (sumPK - Q_stem)/sumK)
    psi_leaf = Max(sp%psi0_LF, psi_stem - Q_leaf/k_stem - psi_ht) ! Ktrunk: (mm/s)/(MPa/m)
  end associate
end subroutine plant_water_potential_equi

!====================================
subroutine plant_water_dynamics_Xiangtao(vegn)
  !@sum leaf and stem water potential and fluxes
  !@+   Weng, 03/25/2022
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  !------local var -----------
  type(cohort_type),pointer :: cc
  real :: transp          ! mm/tree/s
  real :: psi_ht          ! Gravitational water pressure, MPa
  real :: k_rs(soil_L)   ! root to each soil layer conductance
  real :: layer_water_supply(soil_L)
  real :: psi_s0   ! Maximum soil layer water potential
  real :: psi_leaf
  real :: psi_stem
  real :: k_stem ! conductance of tree trunk (with modifications of current states)
  real :: ap, bp, exp_term
  real :: weighted_gw_rate  = 0.d0 ! soil water flux rate from this plant's view
  real :: weighted_gw_cond  = 0.d0 ! soil cond from this plant's view
  real :: wflux_wl  ! Water flux from wood to leaves, kg H2O/tree/step
  real :: wflux_gw  ! Water flux from soil to stems, kg H2O/tree/step
  real :: dpsi
  integer :: i,j

  !  From Xiangtao Xu -------------------------
  !  https://github.com/xiangtaoxu/ED2/blob/master/ED/src/dynamics/plant_hydro.f90
  !  !---------------------------------------------------------------------------!
  !  ! 1.2.3. "Normal case", with positive c_leaf and positive stem_cond.  Check
  !  !        reference X16 for derivation of the equations.
  !  !---------------------------------------------------------------------------!
  !  ap = - stem_cond / c_leaf_d                                           ! [1/s]
  !  bp = ((wood_psi_d - hite_d) * stem_cond - transp_d) / c_leaf_d        ! [m/s]
  !
  !  !----- Project the final leaf psi. -----------------------------------------!
  !  exp_term      = exp(max(ap * dt_d,lnexp_min8))
  !  proj_leaf_psi = max( leaf_psi_lwr_d                                         &
  !                     , ((ap * leaf_psi_d + bp) * exp_term - bp) / ap )
  !  !---------------------------------------------------------------------------!
  !  !----- Calculate the average sapflow rate within the time step [kgH2O/s]. --!
  !  wflux_wl_d = (proj_leaf_psi - leaf_psi_d) * c_leaf_d / dt_d + transp_d
  !    Calculate soil-root water conductance kg H2O/m/s based on reference [K03].
  !---------------------------------------------------------------------------------!
  !gw_cond = soil_cond_d(k) * sqrt(RAI) / (pi18 * dslz8(k))  & ! kg H2O / m3 / s
  !        / nplant_d                                        ! ! conducting area  m2
  !!---------------------------------------------------------------------------------!
  !!      Disable hydraulic redistribution.  Assume roots will shut down if they are
  !! going to lose water to soil.
  !!---------------------------------------------------------------------------------!
  ! if (soil_psi_d(k) <= wood_psi_d) then
  !   gw_cond = 0.d0
  ! end if
  !!---------------------------------------------------------------------------------!
  !!---------------------------------------------------------------------------------!
  !!    Calculate weighted conductance, weighted psi, and water_supply_layer_frac.
  !!---------------------------------------------------------------------------------!
  ! weighted_gw_cond      = weighted_gw_cond + gw_cond                  ! kgH2O/m/s
  ! weighted_soil_psi     = weighted_soil_psi + gw_cond * soil_psi_d(k) ! kgH2O/s
  ! layer_water_supply(k) = gw_cond * (soil_psi_d(k) - wood_psi_d)      ! kgH2O/s
  !!---------------------------------------------------------------------------------!
  !!---------------------------------------------------------------------------------!
  !!     Calculate the average soil water uptake. Check reference X16 for derivation
  !! of the equations.
  !!---------------------------------------------------------------------------------!
  ! ap = - weighted_gw_cond  / c_stem_d  ! ! 1/s
  ! bp = (weighted_soil_psi - wflux_wl_d) / c_stem_d ! m/s
  !!---------------------------------------------------------------------------------!
  !!----- Project the final wood psi, but ensure it will be bounded. ----------------!
  ! exp_term        = exp(max(ap * dt_d,lnexp_min8))
  ! proj_wood_psi   = max( wood_psi_lwr_d                                             &
  !                     , ((ap * wood_psi_d + bp) * exp_term - bp) / ap )
  !!---------------------------------------------------------------------------------!
  !!----- Calculate the average root extraction within the time step [kgH2O/s]. -----!
  ! wflux_gw_d     = (proj_wood_psi - wood_psi_d) * c_stem_d  / dt_d + wflux_wl_d
  !!---------------------------------------------------------------------------------!
  !!------------------------------------------------------------------------------------!
  !!    Now estimate the water uptake from each layer based on layer_water_supply.
  !!------------------------------------------------------------------------------------!
  ! if (sum(layer_water_supply) == 0.d0) then
  !   wflux_gw_layer_d = 0.d0
  ! else
  !   wflux_gw_layer_d = layer_water_supply / sum(layer_water_supply) * wflux_gw_d
  ! end if
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  psi_s0 = maxval(vegn%psi_soil(:))
  do j = 1, vegn%n_cohorts
    cc => vegn%cohorts(j)
    associate ( sp => spdata(cc%species) )

    !-------Cohort specific variables-----------
    transp = cc%transp ! /step_seconds
    psi_ht = HT2MPa(cc%height) ! MPa
    k_stem = cc%Ktrunk * plc_function(cc%psi_stem,sp%psi50_WD,sp%Kexp_WD)

    call Plant_water2psi_exp(cc) ! Refine plant water potential based on water content
    cc%H_leaf = sp%CR_Leaf * cc%W_leaf  ! Leaf Capacitance
    cc%H_stem = sp%CR_Wood * cc%W_stem  ! Stem capacitance

    !! Soil Water psi and K from plant's perspective
    weighted_gw_cond = 0.0     ! kgH2O/m/s
    weighted_gw_rate = 0.0     ! kgH2O/s
    do i=1, soil_L
       k_rs(i)= vegn%K_soil(i)*cc%ArootL(i) ! per tree
       weighted_gw_cond = weighted_gw_cond + k_rs(i)                  ! kgH2O/m/s
       weighted_gw_rate = weighted_gw_rate + k_rs(i) * vegn%psi_soil(i) ! kgH2O/s
       dpsi = max(0.0, vegn%psi_soil(i)-cc%psi_stem)
       layer_water_supply(i) = k_rs(i) * dpsi * step_seconds  ! kgH2O/step
    enddo

    !------------------ Leaf psi and water flux from stems to leaves --------!
    !  ap = - stem_cond / c_leaf_d                                           ! [1/s]
    !  bp = ((wood_psi_d - hite_d) * stem_cond - transp_d) / c_leaf_d        ! [m/s]
    !  exp_term      = exp(max(ap * dt_d,lnexp_min8))
    !  proj_leaf_psi = max( leaf_psi_lwr_d                                         &
    !                     , ((ap * leaf_psi_d + bp) * exp_term - bp) / ap )
    !  wflux_wl_d = (proj_leaf_psi - leaf_psi_d) * c_leaf_d / dt_d + transp_d
    if(cc%Aleaf > 0.0)then
      ap = -k_stem/cc%H_leaf
      bp = ((cc%psi_stem - psi_ht)*k_stem - transp/step_seconds)/cc%H_leaf
      exp_term = exp(ap * step_seconds)
      psi_leaf = max(sp%psi0_LF,((ap*cc%psi_leaf+bp)*exp_term - bp)/ap)
      wflux_wl = (psi_leaf - cc%psi_leaf)*cc%H_leaf + transp
    else
      cc%psi_leaf = cc%psi_stem - psi_ht
      psi_leaf    = cc%psi_leaf
      wflux_wl    = 0.0
      cc%W_leaf   = 0.0
    endif

    !------------------ Stem psi and water flux from soil to leaves------------------
    ! ap = - weighted_gw_cond  / c_stem_d  ! ! 1/s
    ! bp = (weighted_gw_rate - wflux_wl_d) / c_stem_d ! m/s
    ! exp_term        = exp(max(ap * dt_d,lnexp_min8))
    ! proj_wood_psi   = max( wood_psi_lwr_d                                             &
    !                     , ((ap * wood_psi_d + bp) * exp_term - bp) / ap )
    ! wflux_gw_d     = (proj_wood_psi - wood_psi_d) * c_stem_d  / dt_d + wflux_wl_d
    ap = -weighted_gw_cond / cc%H_stem
    bp = (weighted_gw_rate - wflux_wl/step_seconds) / cc%H_stem
    exp_term = exp(ap * step_seconds)
    psi_stem = max( sp%psi0_WD,((ap * cc%psi_stem + bp) * exp_term - bp) / ap)
    wflux_gw = (psi_stem - cc%psi_stem) * cc%H_stem + wflux_wl ! kg H2O/tree/step

    !------------------ Update plant water and hydraulic status ------------------
    cc%Q_leaf = wflux_wl !* step_seconds
    cc%Q_stem = wflux_gw !* step_seconds
    !cc%W_leaf = cc%W_leaf - cc%transp + cc%Q_leaf
    !cc%W_stem = cc%W_stem - cc%Q_leaf + cc%Q_stem
    cc%W_leaf = cc%W_leaf - cc%transp + wflux_wl
    cc%W_stem = cc%W_stem + (psi_stem - cc%psi_stem) * cc%H_stem
    cc%psi_leaf = psi_leaf
    cc%psi_stem = psi_stem

    !------------------ Plant water update from soil layers ------------------
    if (sum(layer_water_supply) == 0.d0) then
       cc%Q_soil = 0.d0
    else
       cc%Q_soil = layer_water_supply/sum(layer_water_supply)*Max(cc%Q_stem,0.0)
    end if

    ! Update soil water in each layer
    do i=1, soil_L
       vegn%wcl(i) = Max(vegn%wcl(i) - (cc%Q_soil(i)*cc%nindivs/(thksl(i)*1000.0)),vegn%WILTPT)
    enddo

    !----------------Next step maximum transpiration---------------
    k_stem = cc%Ktrunk * plc_function(cc%psi_stem,sp%psi50_WD,sp%Kexp_WD)
    cc%W_supply = 0.05 * sum(vegn%freewater(:))*cc%Acrown

    end associate
  enddo

end subroutine plant_water_dynamics_Xiangtao

!===========================================================================
subroutine Plant_water2psi_linear(cc)
  implicit none
  type(cohort_type), intent(inout) :: cc
  !----local vars
  real :: dW_L,dW_S

  associate ( sp => spdata(cc%species) )
    !dW_S =  Max(0.0,MIN(1.0,(cc%Wmax_S - cc%W_stem)/(cc%Wmax_S - cc%Wmin_S)))
    dW_S =  (cc%Wmax_S - cc%W_stem)/(cc%Wmax_S - cc%Wmin_S)
    cc%psi_stem = sp%psi0_WD * dW_S
    if(cc%bl > 0.001)then
      dW_L =  Max(0.0,MIN(1.0,(cc%Wmax_L - cc%W_leaf)/(cc%Wmax_L - cc%Wmin_L)))
      cc%psi_leaf = sp%psi0_LF * dW_L
    else
      cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height)
    endif
  end associate
end subroutine Plant_water2psi_linear

! ====================================

end module esdvm
