#define Hydro_test
!---------------
module esdvm
 use datatypes
 use io_mod
 use soil_mod
 implicit none
 private

 !Core functions
 public :: initialize_vegn_tile, vegn_age,vegn_CNW_budget_fast,vegn_growth
 public :: vegn_phenology,vegn_daily_starvation,vegn_annual_starvation
 public :: vegn_reproduction, vegn_nat_mortality, vegn_hydraulic_states
 public :: relayer_cohorts, vegn_mergecohorts, kill_lowdensity_cohorts
 !For specific experiments
 public :: vegn_fire, vegn_migration, vegn_species_switch, Recover_N_balance
 public :: vegn_annualLAImax_update, vegn_gap_fraction_update

 contains

!========================================================================
!==================== BiomeE surbroutines ===============================
!========================================================================

!=============== Hourly subroutines =====================================
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
  tair   = forcing%Tair -273.16   ! degC
  tsoil  = forcing%tsoil -273.16  ! degC
  thetaS = (vegn%wcl(2)-vegn%WILTPT)/(vegn%FLDCAP-vegn%WILTPT)

  ! Water supply for leaves
  ! Soil water parameters (psi and conductivity for each layer)
  call SoilWater_psi_K(vegn)

#ifdef Hydro_test
  ! Dynamic tree trunk conductivity and water potential
  call Plant_water_dynamics_linear(vegn)  ! A linearized calculation of psi and k
  !call plant_water_dynamics_Xiangtao(vegn) ! Xiangtao's model
  ! Photosynsthesis
  call vegn_photosynthesis(forcing, vegn)
#else
  ! Water supply for photosynthesis from soil layers
  call SoilWaterSupply(vegn)
  call vegn_photosynthesis(forcing, vegn)
  call SoilWaterTranspUpdate(vegn)
#endif

  ! Update soil water: infiltration and surface evap.
  call SoilWaterDynamics(forcing,vegn)

  ! Respiration
  call vegn_respiration(forcing,vegn)

  ! update soil carbon
   call Soil_BGC(vegn, forcing%tsoil, thetaS)

  !! Nitrogen uptake
   call vegn_N_uptake(vegn, forcing%tsoil)

end subroutine vegn_CNW_budget_fast

!========================================================================
! Weng, 08/23/2022: time counter, update cohort ages and the time of them
!                   staying in the first layer
subroutine vegn_age (vegn,t_yr)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: t_yr ! step length (year)

  !----- local var --------------
  type(cohort_type),pointer :: cc
  integer :: i

  ! Tile age
  vegn%age = vegn%age + t_yr
  ! Update cohort ages
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! Update cohort age
     cc%age = cc%age + t_yr
     ! Update time in the top layer
     if (cc%layer == 1) cc%topyear = cc%topyear + t_yr
  enddo
end subroutine vegn_age

!=============== Plant physiology =======================================
! Weng 2017-10-18:compute stomatal conductance, photosynthesis and respiration
! updates cc%An_op and cc%An_cl, from LM3
subroutine vegn_photosynthesis (forcing, vegn)
  type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real :: rad_top  ! downward radiation at the top of the canopy, W/m2
  real :: rad_net  ! net radiation absorbed by the canopy, W/m2
  real :: Tair, TairK     ! air temperature, degC and degK
  real :: cana_q   ! specific humidity in canopy air space, kg/kg
  real :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real :: p_surf   ! surface pressure, Pa
  real :: water_supply ! water supply per m2 of leaves
  real :: fw, fs ! wet and snow-covered fraction of leaves
  real :: psyn   ! net photosynthesis, mol C/(m2 of leaves s)
  real :: resp   ! leaf respiration, mol C/(m2 of leaves s)
  real :: tempLAI,w_scale2, transp ! mol H20 per m2 of leaf per second
  real :: kappa  ! light extinction coefficient of corwn layers
  real :: f_light(10)=0.0      ! light fraction of each layer
  integer :: i, layer

  !! Light supply for photosynthesis
  vegn%kp = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     layer = Max (1, Min(cc%layer,9))
     ! Calculate kappa according to sun zenith angle ! kappa = cc%extinct/max(cosz,0.01) !
     vegn%kp(layer) = vegn%kp(layer)  &  ! -0.75
                    + cc%extinct * cc%Acrown * cc%nindivs
  enddo

  ! Light fraction
  f_light = 0.0
  f_light(1) = 1.0
  do i =2, layer !MIN(int(vegn%CAI+1.0),9)
      f_light(i) = f_light(i-1)  &
           * (exp(0.-vegn%kp(i-1)*vegn%LAIlayer(i-1)) + vegn%f_gap(i-1))
  enddo

  ! Photosynthesis
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
       if(cc%status == LEAF_ON .and. cc%lai > 0.01) then
         ! Convert forcing data
          layer = Max (1, Min(cc%layer,9))
          rad_top = f_light(layer) * forcing%radiation ! downward radiation at the top of the canopy, W/m2
          rad_net = f_light(layer) * forcing%radiation * 0.9 ! net radiation absorbed by the canopy, W/m2
          p_surf  = forcing%P_air  ! Pa
          TairK   = forcing%Tair ! K
          Tair    = forcing%Tair - 273.16 ! degC
          cana_q  = (esat(Tair)*forcing%RH*mol_h2o)/(p_surf*mol_air)  ! air specific humidity, kg/kg
          cana_co2= forcing%CO2 ! co2 concentration in canopy air space, mol CO2/mol dry air
         ! recalculate the water supply to mol H20 per m2 of leaf per second
          water_supply = cc%W_supply/(cc%Aleaf*step_seconds*mol_h2o) ! mol m-2 leafarea s-1

         fw = 0.0; fs = 0.0
         call gs_Leuning(rad_top, rad_net, TairK, cana_q, cc%lai, &
                         p_surf, water_supply, cc%species, sp%pt, &
                         cana_co2, cc%extinct, fs+fw, cc%layer,   &
                         psyn, resp,w_scale2,transp ) ! output
         ! put into cohort data structure for future use in growth
         cc%An_op  = psyn  ! molC s-1 m-2 of leaves
         cc%An_cl  = -resp  ! molC s-1 m-2 of leaves
         cc%w_scale  = w_scale2
         cc%transp = transp * mol_h2o * cc%Aleaf * step_seconds  ! Transpiration (kgH2O/(tree step), Weng, 2017-10-16
         cc%gpp  = (psyn-resp) * mol_C * cc%Aleaf * step_seconds ! kgC step-1 tree-1

         !if(isnan(cc%gpp))cc%gpp=0.0
         if(isnan(cc%gpp))stop '"gpp" is a NaN'
         if(isnan(cc%transp))then
            write(*,*)'w_scale2,transp,lai',w_scale2,transp,cc%lai
            stop '"transp" is a NaN'
         endif
       else
         ! no leaves means no photosynthesis and no stomatal conductance either
         cc%An_op  = 0.0
         cc%An_cl  = 0.0
         cc%gpp    = 0.0
         cc%transp = 0.0
         cc%w_scale  = -9999
       endif
     end associate
  enddo ! vegn, go through all cohorts
end subroutine vegn_photosynthesis

!============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, &
                   p_surf, ws, pft, pt, ca, kappa, f_w, layer, &
                   apot, acl,w_scale2, transp)
  real,    intent(in)   :: rad_top ! PAR dn on top of the canopy, w/m2
  real,    intent(in)   :: rad_net ! PAR net on top of the canopy, w/m2
  real,    intent(in)   :: tl   ! leaf temperature, degK
  real,    intent(in)   :: ea   ! specific humidity in the canopy air, kg/kg
  real,    intent(in)   :: lai  ! leaf area index
  !real,    intent(in)   :: leaf_age ! age of leaf since budburst (deciduos), days
  real,    intent(in)   :: p_surf ! surface pressure, Pa
  real,    intent(in)   :: ws   ! water supply, mol H20/(m2 of leaf s)
  integer, intent(in)   :: pft  ! species
  integer, intent(in)   :: pt   ! physiology type (C3 or C4)
  real,    intent(in)   :: ca   ! concentartion of CO2 in the canopy air space, mol CO2/mol dry air
  real,    intent(in)   :: kappa! canopy extinction coefficient (move inside f(pft))
  real,    intent(in)   :: f_w ! fraction of leaf that's wet or snow-covered
  integer, intent(in)   :: layer  ! the layer of this canopy
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  !real,    intent(out)   :: gs   ! stomatal conductance, m/s
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
  real,    intent(out)   :: w_scale2,transp  ! transpiration, mol H20/(m2 of leaf s)

  ! ---- local vars
  ! photosynthesis
  real :: vm
  real :: kc,ko  ! Michaelis-Menten constants for CO2 and O2, respectively
  real :: ci
  real :: capgam ! CO2 compensation point
  real :: f2,f3
  real :: coef0,coef1
  real :: Resp

  ! conductance related
  real :: gs
  real :: b
  real :: ds  ! humidity deficit, kg/kg
  real :: hl  ! saturated specific humidity at the leaf temperature, kg/kg
  real :: do1

  ! misceleneous
  real :: dum2
  real, parameter :: light_crit = 0
  real, parameter :: gs_lim = 0.25
  real, parameter :: Rgas = 8.314 ! J mol-1 K-1, universal gas constant
  ! new average computations
  real :: lai_eq;
  real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta
  real :: light_top
  real :: par_net
  real :: Ag
  real :: An
  real :: Ag_l
  real :: Ag_rb
  real :: anbar
  real :: gsbar
  real :: w_scale
  real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa
  ! soil water stress
  real :: Ed, an_w, gs_w

  b=0.01
  do1=0.09 ! kg/kg
  if (pft < 2) do1=0.15

  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot
  par_net   = rad_net*rad_phot

  ! Humidity deficit, kg/kg
  call qscomp(tl, p_surf, hl)
  ds = max(hl - ea,0.0)

  associate ( sp => spdata(pft) )
    !  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf
    !  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf
    !  vm=sp%Vmax*exp(3000.0*(1.0/288.2-1.0/tl))

    ! Weng, 2013-01-10
    ko=0.248    * exp(35948/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
    kc=0.000404 * exp(59356/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
    vm=sp%Vmax*exp(24920/Rgas*(1.0/298.2-1.0/tl)) ! / ((layer-1)*1.0+1.0) ! Ea = 33920

    !decrease Vmax due to aging of temperate deciduous leaves
    !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
    !! Turned off by Weng, 2013-02-01, since we can't trace new leaves
    !  if (sp%leaf_age_tau>0 .and. leaf_age>sp%leaf_age_onset) then
    !     vm=vm*exp(-(leaf_age-sp%leaf_age_onset)/sp%leaf_age_tau)
    !  endif

    ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986
    capgam=0.5*kc/ko*0.21*0.209 ! Farquhar & Caemmerer 1982

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
    Resp= sp%gamma_LN/seconds_per_year          & ! per seconds,  mol m-2 s-1
         * sp%LNA * lai / mol_c                 & ! whole canopy, mol m-2 s-1
         * exp(24920/Rgas*(1.0/298.2-1.0/tl))     ! temperature scaled

    ! Temperature effects
    Resp = Resp / ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                 * (1.0+exp(0.4*(tl-45.0-TFREEZE))))

    ! ignore the difference in [CO2] near the leaf and in the canopy air, rb=0.
    Ag_l  = 0.
    Ag_rb = 0.
    Ag    = 0.
    anbar = -Resp/lai
    gsbar = b
    ! find the LAI level at which gross photosynthesis rates are equal
    ! only if PAR is positive
    if(light_top > light_crit)then
       if(pt==PT_C4) then ! C4 species
          coef0=(1+ds/do1)/sp%m_cond;
          ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
          if (ci>capgam) then
             f2=vm
             f3=18000.0*vm*ci ! 18000 or 1800?
             dum2=min(f2,f3)

             ! find LAI level at which rubisco limited rate is equal to light limited rate
             lai_eq = -log(dum2/(kappa*sp%alpha_ps*light_top))/kappa
             lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

             ! gross photosynthesis for light-limited part of the canopy
             Ag_l   = sp%alpha_ps * par_net     &
                    * (exp(-lai_eq*kappa)-exp(-lai*kappa)) &
                    / (1-exp(-lai*kappa))

             ! gross photosynthesis for rubisco-limited part of the canopy
             Ag_rb  = dum2*lai_eq
             Ag=(Ag_l+Ag_rb)/ &
               ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
               *(1.0+exp(0.4*(tl-45.0-TFREEZE))))
             An=Ag-Resp
             anbar=An/lai

             if(anbar>0.0) then
                 gsbar=anbar/(ci-capgam)/coef0;
             endif
          endif ! ci>capgam
       else ! C3 species
          coef0=(1+ds/do1)/sp%m_cond;
          coef1=kc*(1.0+0.209/ko);
          ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
          f2=vm*(ci-capgam)/(ci+coef1);
          f3=vm/2.;
          dum2=min(f2,f3);
          if (ci>capgam) then
             ! find LAI level at which rubisco limited rate is equal to light limited rate
             lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
                         (sp%alpha_ps*light_top*kappa))/kappa
             lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

             ! gross photosynthesis for light-limited part of the canopy
             Ag_l   = sp%alpha_ps              &
                  * (ci-capgam)/(ci+2.*capgam) * par_net   &
                  * (exp(-lai_eq*kappa)-exp(-lai*kappa))  &
                  / (1.0-exp(-lai*kappa))
             ! gross photosynthesis for rubisco-limited part of the canopy
             Ag_rb  = dum2*lai_eq
             Ag = (Ag_l+Ag_rb) /((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                * (1.0+exp(0.4*(tl-45.0-TFREEZE))))
             An = Ag - Resp
             anbar = An/lai
             if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0
             endif
          endif ! ci>capgam
       endif
    endif ! light is available for photosynthesis

    an_w=anbar
    if (an_w > 0.) then
       an_w=an_w*(1-sp%ps_wet*f_w)
    endif
    gs_w = 1.56 * gsbar *(1-sp%ps_wet*f_w) !Weng: 1.56 for H2O?
    if (gs_w > gs_lim) then
        if(an_w > 0.) an_w = an_w*gs_lim/gs_w
        gs_w = gs_lim
    endif
  end associate
  ! find water availability diagnostic demand
  Ed = gs_w * ds*mol_air/mol_h2o ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]
  ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
  if (Ed>ws) then
     w_scale=ws/Ed
     gs_w=w_scale*gs_w
     if(an_w > 0.0) an_w = an_w*w_scale
     if(an_w < 0.0.and.gs_w >b) gs_w=b
  endif
  gs   = gs_w
  apot = an_w
  acl  = -Resp/lai
  transp = min(ws,Ed) ! mol H20/(m2 of leaf s)

  ! Convert units of stomatal conductance to m/s from mol/(m2 s) by
  ! multiplying it by a volume of a mole of gas
  gs = gs * Rugas * Tl / p_surf

  ! for reporting
  w_scale2=min(1.0,ws/Ed)

  ! Error check
  if(isnan(transp))then
    write(*,*)'ws,ed',ws,ed
    stop '"transp" is a NaN'
  endif
end subroutine gs_Leuning

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
  real :: r_Nfix    ! respiration due to N fixation
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
       NSCtarget = 0.5 * (cc%bl_max + cc%br_max)
       fnsc = min(max(0.0,cc%nsc/NSCtarget),1.0)
       Acambium = PI * cc%DBH * cc%height * 1.2

       ! r_leaf = fnsc*sp%gamma_LN  * cc%leafN * tf * dt_fast_yr  ! tree-1 step-1
       ! LeafN  = sp%LNA * cc%Aleaf
       r_leaf   = cc%An_cl * mol_C * cc%Aleaf * step_seconds
       r_stem   = fnsc*sp%gamma_SW  * Acambium * tf * dt_fast_yr ! kgC tree-1 step-1
       r_root   = fnsc*sp%gamma_FR  * cc%rootN * tf * dt_fast_yr ! root respiration ~ root N

       ! Facultive Nitrogen fixation
       !if(cc%NSN < cc%NSNmax .and. cc%NSC > 0.5 * NSCtarget)then
       !   cc%fixedN = sp%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
       !else
       !   cc%fixedN = 0.0 ! sp%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
       !endif

       ! Obligate Nitrogen Fixation
       cc%fixedN = fnsc*sp%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
       r_Nfix    = sp%NfixCost0 * cc%fixedN ! + 0.25*sp%NfixCost0 * cc%N_uptake    ! tree-1 step-1

       ! Total Respiration and NPP
       cc%resl = r_leaf + r_stem ! tree-1 step-1
       cc%resr = r_root + r_Nfix ! tree-1 step-1
       cc%resp = cc%resl + cc%resr + cc%resg/steps_per_day   !kgC tree-1 step-1
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
               cc%N_uptake  = min(cc%br*avgNup, cc%NSNmax- cc%NSN)
               cc%nsn       = cc%nsn + cc%N_uptake
               cc%annualNup = cc%annualNup + cc%N_uptake
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
  real :: dBL, dBR, dBSW, dSeed ! growth of leaf, root, sapwood, and seeds, kgC/individual
  real :: dDBH, dHT, dCA ! tendencies of DBH, height, and Acrown
  real :: BM, DBH1 ! the DBH before grwoth
  real :: N_demand, Nsupplyratio, extraN
  real :: r_N_SD
  integer :: i,j,k

  ! Turnover of leaves and fine roots
  call vegn_tissue_turnover(vegn)
  call SoilWater_psi_K(vegn)

  !Allocate C_gain to tissues
  do i = 1, vegn%n_cohorts
    cc => vegn%cohorts(i)
    associate (sp => spdata(cc%species))
    if (cc%status == LEAF_ON) then
       ! Get carbon from NSC pool
       call fetch_CN_for_growth(cc) ! Weng, 2017-10-19

       ! Allocate carbon to the plant pools
       ! calculate the carbon spent on growth of leaves and roots
       LF_deficit = max(0., cc%bl_max - cc%bl)
       FR_deficit = max(0., cc%br_max - cc%br)
       LFR_deficit = LF_deficit + FR_deficit
       G_LFR = max(min(LF_deficit + FR_deficit, f_LFR_max * cc%C_growth),0.) ! (1.- Wood_fract_min)
       !! and distribute it between roots and leaves
       dBL = min(max(0.,(G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/ &
                        (cc%bl_max + cc%br_max)), G_LFR)
       !! flexible allocation scheme
       !dBL = min(LF_deficit, 0.6*G_LFR)

       if((G_LFR-dBL) > FR_deficit) dBL = G_LFR - FR_deficit
       dBR  = G_LFR - dBL

       ! calculate carbon spent on growth of sapwood growth
       if(cc%layer == 1 .AND. cc%age > sp%AgeRepro)then
           dSeed = sp%v_seed * (cc%C_growth - G_LFR)
           dBSW  = (1.0-sp%v_seed)* (cc%C_growth - G_LFR)
       else
           dSeed= 0.0
           dBSW = cc%C_growth - G_LFR
       endif

       ! For grasses, temporary
       if(sp%lifeform ==0 )then
           dSeed = dSeed + 0.15*G_LFR
           G_LFR = 0.85 * G_LFR
           dBR   = 0.85 * dBR
           dBL   = 0.85 * dBL
       endif
       !! Nitrogen adjustment on allocations between wood and leaves+roots
       !! Nitrogen demand by leaves, roots, and seeds (Their C/N ratios are fixed.)
       N_demand = dBL/sp%CNleaf0 + dBR/sp%CNroot0 + dSeed/sp%CNseed0 + dBSW/sp%CNsw0
       !! Nitrogen available for all tisues, including wood
       if(cc%N_growth < N_demand)then
         ! a new method, Weng, 2019-05-21
         ! same ratio reduction for leaf, root, and seed if(cc%N_growth < N_demand)
         Nsupplyratio = MAX(0.0, MIN(1.0, cc%N_growth/N_demand))
         !r_N_SD = (cc%N_growth-cc%C_growth/sp%CNsw0)/(N_demand-cc%C_growth/sp%CNsw0) ! fixed wood CN
         r_N_SD = cc%N_growth/N_demand ! = Nsupplyratio
         if(sp%lifeform > 0 )then ! for trees
            if(r_N_SD<=1.0 .and. r_N_SD>0.0)then
             dBSW =  dBSW + (1.0-r_N_SD) * (dBL+dBR+dSeed)
             dBR  =  r_N_SD * dBR
             dBL  =  r_N_SD * dBL
             dSeed=  r_N_SD * dSeed
            elseif(r_N_SD <= 0.0)then
             dBSW = cc%N_growth/sp%CNsw0
             dBR  =  0.0
             dBL  =  0.0
             dSeed=  0.0
            endif
         else ! for grasses
            dBSW =  dBSW + (1.0 - Nsupplyratio) * (dBL+dBR+dSeed) ! Nsupplyratio * dBSW !
            dBR  =  Nsupplyratio * dBR
            dBL  =  Nsupplyratio * dBL
            dSeed=  Nsupplyratio * dSeed
         endif
       endif

       !update biomass pools
       cc%bl    = cc%bl    + dBL
       cc%br    = cc%br    + dBR
       cc%bsw   = cc%bsw   + dBSW
       cc%seedC = cc%seedC + dSeed
       cc%NSC   = cc%NSC   - dBR - dBL -dSeed - dBSW
       cc%resg  = 0.5 * (dBR+dBL+dSeed+dBSW) !  daily

       !!update nitrogen pools, Nitrogen allocation
       cc%leafN = cc%leafN + dBL   /sp%CNleaf0
       cc%rootN = cc%rootN + dBR   /sp%CNroot0
       cc%seedN = cc%seedN + dSeed /sp%CNseed0
       cc%sapwN = cc%sapwN + f_N_add * cc%NSN + &
          (cc%N_growth - dBL/sp%CNleaf0 - dBR/sp%CNroot0 - dSeed/sp%CNseed0)
       !extraN = max(0.0,cc%sapwN+cc%woodN - (cc%bsw+cc%bHW)/sp%CNsw0)
       extraN   = max(0.0,cc%sapwN - cc%bsw/sp%CNsw0)
       cc%sapwN = cc%sapwN - extraN
       cc%NSN   = cc%NSN   + extraN - f_N_add*cc%NSN - cc%N_growth !! update NSN
       cc%N_growth = 0.0

       ! accumulated C allocated to leaf, root, and wood
       cc%NPPleaf = cc%NPPleaf + dBL
       cc%NPProot = cc%NPProot + dBR
       cc%NPPwood = cc%NPPwood + dBSW

       ! Keep previous DBH
       DBH1 = cc%DBH

       !! update plant architecture given increase of bsw
       BM        = cc%bsw + cc%bHW
       call BM2Architecture(cc,BM)
       !dDBH = dBSW / (sp%thetaBM * sp%alphaBM * cc%DBH**(sp%thetaBM-1.0))
       !dHT  = sp%thetaHT * sp%alphaHT * cc%DBH**(sp%thetaHT-1) * dDBH
       !dCA  = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * dDBH
       !cc%DBH    = BM2DBH(    BM,cc%species) !cc%DBH    = cc%DBH    + dDBH
       !cc%height = DBH2HT(cc%DBH,cc%species) !cc%height = cc%height + dHT
       !cc%Acrown = DBH2CA(cc%DBH,cc%species) !cc%Acrown = cc%Acrown + dCA

       ! Update bl_max, br_max, and NSNmax with shifts from understory to the top layer
       call update_max_LFR_NSN(cc)

#ifdef Hydro_test
       ! Update plant hydraulic states
       call Update_hydraulic_vars(cc)
       if(cc%firstday)then
         cc%psi_s0   = maxval(vegn%psi_soil(:))
         cc%psi_stem = cc%psi_s0
         cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height)
       endif

       ! Update Ktrunk with new sapwood
       k = Max(MIN(cc%Nrings, Ysw_max),1)
       cc%Kx(k)   = NewWoodKx(cc)
       cc%Lring(k)= HT2Lpath(cc%height)
       cc%Ktrunk  = cc%Ktrunk+ &
             0.25*PI*(cc%DBH**2-DBH1**2)*cc%Kx(k)/cc%Lring(k)

       ! Plant water status update
       call Update_hydraulic_vars(cc)
       call Plant_water2psi_exp(cc)
#else
       !Convert C and N from sapwood to heartwood
       call Sap2HeartWood_fixedHv(cc)
#endif

     elseif(cc%status == LEAF_OFF .and. cc%C_growth > 0.)then
       cc%nsc = cc%nsc + cc%C_growth
       cc%resg = 0.0
     endif ! "cc%status == LEAF_ON"

     ! reset carbon acculmulation terms
     cc%C_growth = 0

  end associate ! F2003
  enddo
  cc => null()

  ! Update tile variables
  call vegn_sum_tile(vegn)

end subroutine vegn_growth ! daily

!========= Calculate carbon and nitrogen supply ==========================
subroutine fetch_CN_for_growth(cc)
  !@sum Fetch C from labile C pool according to the demand of leaves and fine roots,
  !@+   and the push of labile C pool
  !@+   Daily call.
  !@+   added by Weng, 12-06-2016

  implicit none
  type(cohort_type), intent(inout) :: cc

  !------local var -----------
  logical :: woody
  logical :: dormant,growing
  real :: NSCtarget, bl_max, br_max
  real :: C_push, C_pull, growthC
  real :: N_push, N_pull, growthN
  real :: LFR_rate
  real :: resp_growth

  associate ( sp => spdata(cc%species) )
    LFR_rate = sp%LFR_rate ! 1.0 !  1.0/16.0 ! filling rate/day
    NSCtarget = 3.0 * (cc%bl_max + cc%br_max)      ! kgC/tree

    ! Fetch C from labile C pool if it is in the growing season
    if (cc%status == LEAF_ON) then ! growing season
        C_pull = LFR_rate * (Max(cc%bl_max - cc%bl,0.0) +   &
                  Max(cc%br_max - cc%br,0.0))
        N_pull = LFR_rate * (Max(cc%bl_max - cc%bl,0.0)/sp%CNleaf0 +  &
                  Max(cc%br_max - cc%br,0.0)/sp%CNroot0)

        C_push = cc%nsc/(days_per_year*sp%tauNSC) ! max(cc%nsc-NSCtarget, 0.0)/(days_per_year*sp%tauNSC)
        N_push = cc%NSN/(days_per_year*sp%tauNSC)

        cc%C_growth = Min(max(0.02*cc%NSC,0.0), C_pull + C_push)
        cc%N_growth = Min(max(0.02*cc%NSN,0.0), N_pull + N_push)
    else ! non-growing season
        cc%C_growth = 0.0
        cc%N_growth = 0.0
        cc%resg     = 0.0
    endif
  end associate

end subroutine fetch_CN_for_growth

!================================================================
subroutine update_max_LFR_NSN(cc)
  !@sum: Daily call for calculating bl_max and br_max
  !@+   added by Weng, 08-19-2022
  implicit none
  type(cohort_type), intent(inout) :: cc

  !----- local vars ---------
  real :: BL_c, BL_u

  ! Update bl_max and br_max daily, Weng 2014-01-23, 2021-06-04, 08/24/2022
  ! The new updates allow a gradual increase of BLmax when a tree enters
  ! the canopy layer and a abrupt increase for grasses.
  associate ( sp => spdata(cc%species) )
    BL_c = sp%LMA * sp%LAImax * cc%Acrown * (1.0-sp%f_cGap)
    if (cc%layer > 1) then
      BL_u = BL_c/cc%layer
    else ! cc%layer = 1
      BL_u = BL_c / (1+cc%layer)            ! Woody plants only
    endif
    cc%bl_max = BL_u + min(1., cc%topyear/sp%transT) * (BL_c - BL_u)

    ! Root max
    cc%br_max = BLmax2BRmax(cc)
    ! NSN max
    cc%NSNmax = ccNSNmax(cc)
  end associate
end subroutine update_max_LFR_NSN

!============================================================================
! Updated by Weng, 06-04-2021
subroutine vegn_phenology(vegn) ! daily step
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc
  integer :: i,j
  integer :: Days_thld = 60 ! minimum days of the growing or non-growing season
  real    :: cold_thld = -20.  ! threshold of accumulative low temperature
  real    :: GDD_adp, Tc_off_crit
  real    :: Tc_adj ! Tc_critical adjust according to growing season lenght
  real    :: ccNSC, ccNSN
  real    :: totC, totN
  logical :: PhenoON, PhenoOFF

  ! -------------- update vegn GDD and tc_pheno ---------------
  ! environmental factors for each cohort
  vegn%tc_pheno = vegn%tc_pheno * 0.8 + vegn%Tc_daily * 0.2
  do i = 1, vegn%n_cohorts
    cc=>vegn%cohorts(i)
    associate (sp => spdata(cc%species) )
      if(cc%status == LEAF_ON)then
         cc%ncd = 0
         cc%ndm = 0
         cc%gdd = 0.0
         cc%ngd = Min(366, cc%ngd + 1)
         if(cc%ngd > Days_thld)cc%ALT = cc%ALT + MIN(0.,vegn%tc_pheno-sp%tc_crit_off)
      else  ! cc%status == LEAF_OFF
         cc%ngd = 0
         cc%ALT = 0.0
         cc%ndm = cc%ndm + 1
         if(vegn%tc_pheno<sp%tc_crit_off)then
            cc%ncd = cc%ncd + 1
         endif
         ! Keep gdd as zero in early non-growing season when days < 60
         if(cc%ndm>Days_thld)cc%gdd = cc%gdd + max(0.0,vegn%tc_pheno-T0_gdd)
      endif ! cc%status
    end associate
  enddo

  ! --------- Change pheno status ----------------------------
  ! ON and OFF of phenology: change the indicator of growing season for deciduous
  do i = 1,vegn%n_cohorts
    cc => vegn%cohorts(i)
    associate (sp => spdata(cc%species) )
      !for evergreen
      if(sp%phenotype==1 .and. cc%status /= LEAF_ON) cc%status=LEAF_ON
      !for deciduous and grasses
      ! GDD_adp = sp%gdd_crit*exp(gdd_par3*cc%ncd) + gdd_par1 ! for adaptive phenology
      PhenoON = ((sp%phenotype==0 .and. cc%status/=LEAF_ON)     &
         ! Temperature conditions
         .and.(cc%gdd>sp%gdd_crit .and. vegn%tc_pheno>sp%tc_crit_on)  &
         !!!  Woody plants            Grasses in the top layer   !!!
         !.and.(sp%lifeform==1 .OR.(sp%lifeform==0 .and. cc%layer==1))  &
         )

      cc%firstday = .false.
      if(PhenoON)then
          cc%status = LEAF_ON ! Turn on a growing season
          cc%firstday = .True.
      endif

      ! Reset grass density at the first day of a growing season
      if(cc%firstday .and. sp%lifeform ==0 .and. cc%age>1.)then
          !        reset grass density and size for perenials
          ccNSC   = (cc%NSC +cc%bl +  cc%bsw  +cc%bHW  +cc%br   +cc%seedC) * cc%nindivs
          ccNSN   = (cc%NSN +cc%leafN+cc%sapwN+cc%woodN+cc%rootN+cc%seedN) * cc%nindivs
          ! reset
          cc%nindivs = MIN(ccNSC /sp%s0_plant, ccNSN/(sp%s0_plant/sp%CNroot0))
          totC = ccNSC / cc%nindivs
          totN = ccNSN / cc%nindivs
          call setup_seedling(cc,totC,totN)
      endif
    end associate
  enddo

  if(PhenoON) call relayer_cohorts(vegn)

  ! OFF of a growing season
  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate (sp => spdata(cc%species) )
     !Tc_adj = - 5. * exp(-0.05*max(-15.,real(cc%ngd-N0_GD)))
     Tc_off_crit = sp%tc_crit_off - 5. * exp(-0.05*(cc%ngd-N0_GD))
     PhenoOFF = (cc%status==LEAF_ON .and. sp%phenotype == 0 .and. &
                cc%ALT < cold_thld .and. vegn%tc_pheno < Tc_off_crit)
     end associate

     if(PhenoOFF )then
        cc%status = LEAF_OFF  ! Turn off a growing season
        cc%gdd    = 0.0        ! Start to counting a new cycle of GDD
     endif
     ! leaf fall
     call Seasonal_fall(cc,vegn)
  enddo
end subroutine vegn_phenology

! ============================================================================
subroutine vegn_tissue_turnover(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  !-------local var
  type(cohort_type), pointer :: cc    ! current cohort
  real :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real :: alpha_L   ! turnover rate of leaves
  real :: alpha_S   ! turnover rate of stems
  real :: dBL, dBR, dBStem  ! leaf and fine root carbon tendencies
  real :: dNL, dNR, dNStem  ! leaf and fine root nitrogen tendencies
  real :: dAleaf ! leaf area decrease due to dBL
  integer :: i

  ! update plant carbon and nitrogen for all cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
     !    Turnover of leaves and roots regardless of the STATUS of leaf
     !    longevity. Deciduous: 0; Evergreen 0.035/LMa
     !    root turnover
     if(cc%status==LEAF_OFF)then
        alpha_L = sp%alpha_L ! 60.0 ! yr-1, for decuduous leaf fall
     else
        alpha_L = sp%alpha_L
     endif
     ! Stem turnover
     if(sp%lifeform == 0)then
        alpha_S = alpha_L
     else
        alpha_S = 0.0
     endif
     dBL    = cc%bl    *    alpha_L  /days_per_year
     dNL    = cc%leafN *    alpha_L  /days_per_year

     dBStem = cc%bsw   *    alpha_S  /days_per_year
     dNStem = cc%sapwN *    alpha_S  /days_per_year

     dBR    = cc%br    * sp%alpha_FR /days_per_year
     dNR    = cc%rootN * sp%alpha_FR /days_per_year

     dAleaf = BL2Aleaf(dBL,cc)

     !    Retranslocation to NSC and NSN
     cc%nsc = cc%nsc + l_fract  * (dBL + dBR + dBStem)
     cc%NSN = cc%NSN + retransN * (dNL + dNR + dNStem)
     !    update plant pools
     cc%bl    = cc%bl    - dBL
     cc%bsw   = cc%bsw   - dBStem
     cc%br    = cc%br    - dBR

     cc%leafN = cc%leafN - dNL
     cc%sapwN = cc%sapwN - dNStem
     cc%rootN = cc%rootN - dNR

     !    update leaf area and LAI
     cc%Aleaf= BL2Aleaf(cc%bl,cc)
     cc%lai     = cc%Aleaf/(cc%Acrown *(1.0-sp%f_cGap))

     !    update NPP for leaves, fine roots, and wood
     cc%NPPleaf = cc%NPPleaf - l_fract * dBL
     cc%NPProot = cc%NPProot - l_fract * dBR
     cc%NPPwood = cc%NPPwood - l_fract * dBStem

     !    put C and N into soil pools
     loss_coarse  = (1.-l_fract) * cc%nindivs * (dBL - dAleaf * LMAmin    + dBStem)
     loss_fine    = (1.-l_fract) * cc%nindivs * (dBR + dAleaf * LMAmin)
     lossN_coarse = (1.-retransN)* cc%nindivs * (dNL - dAleaf * sp%LNbase + dNStem)
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
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse

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
  real    :: dAleaf, dBL, dBR, dNL, dNR, dBStem, dNStem      ! per day
  real    :: leaf_fall_rate, root_mort_rate      ! per day

  leaf_fall_rate = 0.05
  root_mort_rate = 0.025
  !End a growing season: leaves fall for deciduous
  associate (sp => spdata(cc%species) )
  if(cc%status == LEAF_OFF .AND. cc%bl > 0.0)then
     dBL = min(leaf_fall_rate * cc%bl_max, cc%bl)
     dBR = min( root_mort_rate * cc%br_max, cc%br)  ! Just for test: keep roots
     dBStem = 0.0 ! trees
     dNStem = 0.0 ! trees
     if(sp%lifeform == 0)then  ! grasses
         dBStem = MIN(1.0,dBL/cc%bl) * cc%bsw
         dNStem = MIN(1.0,dBL/cc%bl) * cc%sapwN
     endif
     ! Nitrogen out
     if(cc%bl>0)then
        dNL = dBL/cc%bl * cc%leafN !dBL/sp%CNleaf0
     else
        dNL = 0.0
     endif
     if(cc%br>0)then
        dNR = dBR/cc%br * cc%rootN !dBR/sp%CNroot0
     else
        dNR = 0.0
     endif

     dAleaf = BL2Aleaf(dBL,cc)
     ! Put plant water into the first soil layer
     if(cc%bl >0.0)  &
     vegn%wcl(1) = vegn%wcl(1) + cc%nindivs * cc%W_leaf*dBL/cc%bl/(thksl(1)*1000.0)
     if(cc%bSW>0.0) &
     vegn%wcl(1) = vegn%wcl(1) + cc%nindivs * cc%W_stem*dBStem/cc%bSW/(thksl(1)*1000.0)

     !Retranslocation to NSC and NSN
     cc%nsc = cc%nsc + l_fract  * (dBL + dBR + dBStem)
     cc%NSN = cc%NSN + retransN * (dNL + dNR + dNStem)
     !update plant pools
     cc%bl    = cc%bl  - dBL
     cc%br    = cc%br  - dBR
     cc%bsw   = cc%bsw - dBStem ! for grass

     cc%leafN = cc%leafN - dNL
     cc%rootN = cc%rootN - dNR
     cc%sapwN = cc%sapwN - dNStem

     !update NPP for leaves, fine roots, and wood
     cc%NPPleaf = cc%NPPleaf - l_fract * dBL
     cc%NPProot = cc%NPProot - l_fract * dBR
     cc%NPPwood = cc%NPPwood - l_fract * dBStem
     cc%Aleaf= BL2Aleaf(cc%bl,cc)
     cc%lai     = cc%Aleaf/(cc%Acrown *(1.0-sp%f_cGap))

     !put C and N into soil pools:  Substraction of C and N from leaf and root pools
     loss_coarse  = (1.-l_fract) * cc%nindivs * (dBStem+dBL - dAleaf * LMAmin)
     loss_fine    = (1.-l_fract) * cc%nindivs * (dBR        + dAleaf * LMAmin)
     lossN_coarse = (1.-retransN)* cc%nindivs * (dNStem+dNL - dAleaf * sp%LNbase)
     lossN_fine   = (1.-retransN)* cc%nindivs * (dNR        + dAleaf * sp%LNbase)

     vegn%SOC(1) = vegn%SOC(1) +  &
                      fsc_fine * loss_fine + fsc_wood * loss_coarse
     vegn%SOC(2) = vegn%SOC(2) +   &
                      (1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse
     vegn%SON(1)  = vegn%SON(1) +    &
                       fsc_fine * lossN_fine + fsc_wood * lossN_coarse
     vegn%SON(2) = vegn%SON(2) +   &
                       (1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse

     !annual N from plants to soil
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
  endif
  end associate

end subroutine Seasonal_fall

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
  integer,dimension(16) :: reproPFTs
  real,   dimension(16) :: seedC, seedN ! seed pool of productible PFTs
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
          !       vegn%N_P2S_yr = vegn%N_P2S_yr + N_failedseed
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
     cc%sapwN  = cc%bsw/sp%CNsw0
     cc%woodN  = cc%bHW/sp%CNwood0
     cc%seedN  = 0.0
     cc%NSN    = totN - (cc%leafN+cc%rootN+cc%sapwN) !cc%br/sp%CNroot0

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
     cc%Rring = 0.0
     cc%Lring = 0.0
     cc%Aring = 0.0

     ! First ring
     cc%Nrings  = 1
     cc%WTC0(1) = NewWoodWTC(cc)
     cc%Kx(1)   = NewWoodKx(cc)
     cc%Lring(1)= HT2Lpath(cc%height)
     cc%Aring(1)= PI * 0.25*cc%DBH**2
     cc%Ktrunk  = PI * 0.25*cc%DBH**2 * cc%Kx(1)/cc%Lring(1)
     cc%treeHU  = 0.0
     cc%treeW0  = cc%WTC0(1) * cc%Aring(1)
     cc%Asap    = cc%Aring(1)
     cc%Atrunk  = cc%Aring(1)

     ! Cohort hydraulic properties
     call Update_hydraulic_vars(cc)
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
       if (cc%nsc < 0.01*cc%bl_max .OR. cc%annualNPP < 0.0) then ! annualNPP < 0 is for grasses only
           deathrate = 1.0
           deadtrees = cc%nindivs * deathrate !individuals / m2
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
     lossN_coarse = deadtrees * (cc%woodN+cc%sapwN + cc%leafN - cc%Aleaf*sp%LNbase)
     lossN_fine   = deadtrees * (cc%rootN+cc%seedN + cc%NSN   + cc%Aleaf*sp%LNbase)

     ! Assume water in plants goes to first layer of soil
     vegn%wcl(1) = vegn%wcl(1) +  deadtrees * (cc%W_leaf+cc%W_stem+cc%W_dead)/(thksl(1)*1000.0)

     vegn%SOC(1) = vegn%SOC(1) + fsc_fine *loss_fine + fsc_wood *loss_coarse
     vegn%SOC(2) = vegn%SOC(2) + (1.0-fsc_fine)*loss_fine + (1.0-fsc_wood)*loss_coarse

     vegn%SON(1) = vegn%SON(1) + &
                fsc_fine *lossN_fine +    fsc_wood *lossN_coarse
     vegn%SON(2) = vegn%SON(2) + &
                (1.-fsc_fine)*lossN_fine +(1.-fsc_wood)*lossN_coarse

     ! annual N from plants to soil
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse

     end associate
end subroutine plant2soil

!============================================================================
real function mortality_rate(cc) result(mu) ! per year
  !@sum calculate cohort mortality/year, Ensheng Weng, 12/07/2021
  ! Mortality rate should be a function of growth rate, age, and environmental
  ! conditions. Here, we only used used a couple of parameters to calculate
  ! mortality as functions of social status, seedling size, and adult size.
  ! Grass is saprately defined.
  type(cohort_type),intent(in) :: cc

  !-------local var -------------
  integer :: n ! the latest ring
  real :: f_L, f_S, f_D ! Layer, seeding, and size effects on mortality
  real :: expD, m_S ! Mortality multiplier for size effects
  real :: mu_bg, mu_hydro  ! Background mortality rate and hydraulic failure

  !---------------------
  associate ( sp => spdata(cc%species))
    n = MIN(cc%Nrings, Ysw_max)
    f_L  = sp%A_un * SQRT(Max(0.0, cc%layer - 1.0)) ! Layer effects (0~ infinite)
    f_S  = sp%A_sd * exp(sp%B_sd*cc%dbh) + 1.0      ! Understory seedling

    ! Size effect on the mortality of adult trees
    if(do_U_shaped_mortality)then
       m_S = 5.0
    else
       m_S = 0.0
    endif
    expD = exp(sp%A_D * (cc%dbh - sp%D0mu))
    f_D  = 1. + m_S * expD / (1. + expD) ! Size effects (big tees)

    mu_bg = Min(0.5,sp%r0mort_c * (1.d0+f_L*f_S)*f_D) ! per year

#ifdef Hydro_test
    ! Trunk hydraulic failure probability
    mu_hydro = exp(sp%s_hu * (1.0 - min(1.,cc%treeHU/cc%treeW0)))
    !mu_hydro = Max(0., 1. - cc%farea(n))
    !mu_hydro = Max(0., 1. - cc%Asap/cc%Acrown/(sp%LAImax*sp%phiCSA))

#endif
  end associate

  ! Return mortality rate:
  mu = mu_bg + mu_hydro - mu_bg * mu_hydro ! Add hydraulic failure

end function mortality_rate

!============================================================================
!-----------------------Plant Hydraulics------------------------------
subroutine vegn_hydraulic_states(vegn, deltat)
  ! Update plant hydraulic states and coverstion of sapwood to heartwood
  ! yearly time step
  ! Author: Ensheng Weng, 2021-03-15, updated 2021-12-8
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
          cc%Rring(1) = cc%DBH/2.0
          cc%Lring(1) = HT2Lpath(cc%height)
          cc%Aring(1) = PI * cc%Rring(1)**2
       endif
       ! Set up the space for the new ring : moving previous years' states inward
       if(cc%Nrings >= Ysw_max)then
          do j=2, Ysw_max
             cc%WTC0(j-1)  = cc%WTC0(j)
             cc%Kx(j-1)    = cc%Kx(j)
             cc%farea(j-1) = cc%farea(j)
             cc%accH(j-1)  = cc%accH(j)
             cc%Rring(j-1) = cc%Rring(j)
             cc%Lring(j-1) = cc%Lring(j)
             cc%Aring(j-1) = cc%Aring(j)
          enddo
       endif

       ! Growth: the new ring formed in this year
       cc%Nrings = cc%Nrings + 1 ! A new ring
       k = MIN(cc%Nrings, Ysw_max)
       ! WTC0 and Kx represent scientific hypotheses. They can be constant,
       ! or functions of environmental conditions, growht rates, etc.
       cc%WTC0(k) = NewWoodWTC(cc)
       cc%Kx(k)   = NewWoodKx(cc)
       ! Other cohort variables of the new ring
       cc%farea(k) = 1.0
       cc%accH(k)  = 0.0
       cc%Rring(k) = cc%DBH/2.0
       cc%Lring(k) = HT2Lpath(cc%height)
       if(k>1)then
          cc%Aring(k) = PI * Max(0.0,cc%Rring(k)**2 - (cc%DBH_ys/2.)**2)
       else
          cc%Aring(k) = PI * cc%Rring(k)**2 ! Only for the first year
       endif
       call calculate_Asap_Ktrunk (cc) ! Tree trunk conductance and sapwood area
       ! Update tree total water transport
       ! cc%treeHU   = cc%treeHU + cc%Acrown * cc%annualTrsp * 1.e-3 ! ton/tree
       ! cc%treeW0   = cc%treeW0 + cc%Aring(k) * cc%WTC0(k)

       ! Update trunk hydraulic status
       !Trsp_sap = 1.e-3 * cc%annualTrsp/cc%Asap ! m, usage of functional conduits
       Lmax = HT2Lpath(cc%height) ! maxval(cc%Lring(:)) !
       cc%treeHU = 0.0
       cc%treeW0 = 0.0
       do k=1, MIN(cc%Nrings, Ysw_max)
         !cc%accH(k) = cc%accH(k) + Trsp_sap ! Assume the same

         ! Lifetime water transported for functional xylem conduits
         funcA = cc%farea(k) * cc%Aring(k)
         if(funcA > 1.0e-8)then
           trsp_ring = 1.e-3 * cc%annualTrsp * cc%Kring(k)/cc%Ktrunk ! ton, per ring
           cc%accH(k) = cc%accH(k) + trsp_ring/funcA  ! m, for functional conduits only
           !cc%accH(k) = cc%accH(k) + 1.e-3 * cc%annualTrsp * cc%Kx(k)/Lmax/cc%Ktrunk
         endif

         ! Update each ring's functional fraction
         cc%farea(k) = 1. - exp(-r_DF * (1. - MIN(1.0,cc%accH(k)/cc%WTC0(k))))
         !cc%farea(k) = 1. - 1./(1. + exp(r_DF * (1. - cc%accH(k)/cc%WTC0(k))))

         ! Update tree hydraulic usage and WTC0
         funcA = cc%farea(k) * cc%Aring(k)
         cc%treeHU   = cc%treeHU + funcA * cc%accH(k)
         cc%treeW0   = cc%treeW0 + funcA * cc%WTC0(k)

       enddo
       call calculate_Asap_Ktrunk (cc) ! Update Asap and Ktrunk
     end associate

     ! ----- Conversion of sapwood to heartwood -----
     call Sap2HeartWood_Hydro(cc)

     ! Update other variables
     call Update_hydraulic_vars(cc)
     call Plant_water2psi_exp(cc)
  enddo
end subroutine vegn_hydraulic_states

!========================================================================
! Weng 2022-03-29 ! Updated 07/29/2022
subroutine Plant_water_dynamics_linear(vegn)     ! forcing,
  !type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real :: Q_air          ! transpiration amount in this step
  real :: psi_ht         ! Gravitational water pressure, MPa
  real :: psi_leaf,psi_stem,psi_sl
  real :: k_rs(soil_L)   ! soil-root water conductance by soil layer
  real :: k_stem         ! The conductance of the tree at current conditions
  real :: sumK, sumPK, dpsi, psi_soil, W_psi_soil
  integer :: i,j

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
       k_stem = cc%Ktrunk * plc_function(psi_sl,sp%psi50_WD,sp%Kexp_WD)
       psi_ht = HT2MPa(cc%height) ! MPa

       !Approximately estimate psi_leaf and Q_leaf
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

       ! Water fluxes from soil to stem base (and between soil layers via roots)
       W_psi_soil = cc%Wmax_s * exp(cc%psi_s0 * sp%CR_Wood)
       if(cc%W_stem < W_psi_soil) then
         sumK  = 0.0
         sumPK = 0.0
         do i=1, soil_L
           ! RAI(i) = cc%rootareaL(i)/cc%Acrown
           ! k_rs(i)= vegn%K_soil(i)*SQRT(RAI(i))/(3.14159*thksl(i))
           k_rs(i)= vegn%K_soil(i)*cc%rootareaL(i) ! per tree
           !if(vegn%psi_soil(i)<cc%psi_stem) k_rs(i)=0.0
           sumK = sumK + k_rs(i)
           sumPK= sumPK+ k_rs(i) * vegn%psi_soil(i)
         enddo
         psi_stem = (sumPK*step_seconds - cc%Q_leaf + cc%psi_stem*cc%H_stem)/ &
                    (sumK *step_seconds + cc%H_stem)
         ! Calculate water uptake by layer in theory
         do i=1, soil_L
            dpsi = max(0.0, vegn%psi_soil(i) - psi_stem)
            cc%Q_soil(i) = k_rs(i) * dpsi * step_seconds
            cc%Q_soil(i) = Min(cc%Q_soil(i), 0.05*vegn%freewater(i)*cc%Acrown)
         enddo

         ! Check if the plant needs this ammount of water
         if(sum(cc%Q_soil)>0.0) &
           cc%Q_soil = cc%Q_soil*Min(1.0,(W_psi_soil-cc%W_stem)/sum(cc%Q_soil))

         ! Update stem and soil water content, and stem psi
         do i=1, soil_L
           cc%W_stem  = cc%W_stem   + cc%Q_soil(i)
           vegn%wcl(i)= vegn%wcl(i) - cc%Q_soil(i)*cc%nindivs/(thksl(i)*1000.0)
         enddo
         cc%psi_stem = psi_stem ! log(cc%W_stem/cc%Wmax_s)/sp%CR_Wood
       endif

       ! Water supply for regulating stomata conductance
       cc%W_supply = PlantWaterSupply(cc,step_seconds)

       !------------------------
       if(isnan(cc%W_supply))then
          write(*,*)'cc%species',cc%species
          write(*,*)'cc%bl,cc%bsw',cc%bl,cc%bsw
          write(*,*)'cc%Wmin_L,cc%Wmax_L',cc%Wmin_L,cc%Wmax_L
          write(*,*)'cc%Wmin_s,cc%Wmax_s',cc%Wmin_s,cc%Wmax_s
          write(*,*)'cc%W_leaf,cc%W_stem',cc%W_leaf,cc%W_stem
          write(*,*)'cc%H_leaf,cc%H_stem',cc%H_leaf,cc%H_stem
          write(*,*)'cc%psi_leaf,cc%psi_stem',cc%psi_leaf,cc%psi_stem
          stop '"transp" is a NaN'
       endif
     end associate
  enddo
end subroutine Plant_water_dynamics_linear

!===================================================================
real function PlantWaterSupply(cc,step_seconds) result(pws)
  !@sum: Calculate water supply for transpiration from plant water status
  ! Weng, 07/13/2022
  type(cohort_type),intent(in) :: cc
  real,             intent(in) :: step_seconds

  !----- Local vars ---------------
  real :: expK,psi50, CR0, f0_sup
  real :: W_stem,W_leaf
  real :: k_stem,psi_stem
  real :: psi0, dpsi
  real :: S_stem,S_leaf,wflux
  real :: step_base = 600. ! Seconds
  integer :: n_iterations,i
  logical :: do_simple_W_supply = .True. ! .False. !
  !---------------------

  ! Prameters
  associate ( sp => spdata(cc%species) )
    psi50  = sp%psi50_WD
    expk   = sp%Kexp_WD
    CR0    = sp%CR_Wood
    f0_sup = sp%f_supply * step_seconds/3600.0
  end associate

  !Leaf water supply
  S_leaf = f0_sup * (cc%W_leaf - cc%Wmin_L)
  ! Stem water supply
  if(do_simple_W_supply)then ! An arbitrary assignment
    S_stem = f0_sup * Max((cc%W_stem - cc%Wmin_s),0.0) ! + (cc%psi_stem - cc%psi_leaf)* k_stem * step_seconds
  else ! Calculated as a function of woody properties
    n_iterations = int(step_seconds/step_base)
    W_stem = cc%W_stem
    psi0   = psi50 ! the hypothetical terminal branch water pressure
    S_stem = 0.0
    do i =1, n_iterations
      W_stem = MIN(max(0.5*cc%Wmin_S, W_stem), cc%Wmax_s)
      psi_stem = log(W_stem/cc%Wmax_s)/CR0
      dpsi     = psi_stem - psi0
      if(dpsi <= 0.02)exit
      k_stem = cc%Ktrunk * plc_function(psi_stem,psi50,expK)
      wflux  = MIN(dpsi * k_stem * step_base, 0.5 * (W_stem - cc%Wmin_S))
      W_stem = W_stem - wflux
      S_stem = S_stem + wflux
      if(S_stem >= f0_sup * (cc%W_stem - cc%Wmin_s))exit
    enddo
  endif

  ! Total plant water supply
  pws = S_stem + S_leaf

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
  cc%Atrunk = 0.0
  do i=1, MIN(cc%Nrings, Ysw_max)
    cc%Asap     = cc%Asap   + cc%Aring(i) * cc%farea(i)
    cc%Atrunk   = cc%Atrunk + cc%Aring(i)
    cc%Kring(i) = cc%farea(i)*cc%Aring(i) * cc%Kx(i)/maxval(cc%Lring(:))
    cc%Ktrunk   = cc%Ktrunk + cc%Kring(i)
  enddo
end subroutine calculate_Asap_Ktrunk

!=================================================
! Weng: Update C and N of sapwood and heartwood, 08/22/2022
subroutine Sap2HeartWood_Hydro(cc)
  type(cohort_type), intent(inout) :: cc
  !----------local var ----------
  real :: D_hw, woodC, newSW, woodN, dSW

  ! Woody plants only
  associate (sp => spdata(cc%species) )
  if(sp%lifeform>0)then
    ! Calculate heartwood diameter
    D_hw  = SQRT(Max(0.0, PI*(cc%DBH/2)**2 - cc%Asap)/PI)
    woodC  = cc%bsw   + cc%bHW
    woodN  = cc%sapwN + cc%woodN

    ! Convert sapwood to heart wood
    newSW = sp%alphaBM * (cc%dbh**sp%thetaBM - D_hw**sp%thetaBM)
    dSW   = cc%bSW - newSW
    cc%W_stem = cc%W_stem - cc%W_stem * dSW/cc%bSW
    cc%W_dead = cc%W_dead + cc%W_stem * dSW/cc%bSW
    cc%bHW    = cc%bHW + dSW
    cc%bsw    = woodC - cc%bHW
    cc%sapwN  = woodN * cc%bsw/woodC
    cc%woodN  = woodN * cc%bHW/woodC
  endif
  end associate
end subroutine Sap2HeartWood_Hydro

!================================================================
subroutine Sap2HeartWood_fixedHv(cc)
  !@sum: Daily call for converting sapwood to heartwood with fixed Hv
  !@+   added by Weng, 08-19-2022
  implicit none
  type(cohort_type), intent(inout) :: cc
  !--- Local vars ------
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: dSW    ! Sapwood to Heartwood, kgC/individual
  real :: dNS    ! Nitrogen from SW to HW

  !-------------------------
  associate ( sp => spdata(cc%species) )
    if(sp%lifeform>0)then ! woody plants
       CSAsw  = cc%bl_max/sp%LMA * sp%phiCSA * cc%height ! with Plant hydraulics, Weng, 2016-11-30
       CSAtot = 0.25 * PI * cc%DBH**2
       CSAwd  = max(0.0, CSAtot - CSAsw)
       DBHwd  = 2*sqrt(CSAwd/PI)
       BSWmax = sp%alphaBM * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
       dSW   = max(cc%bsw - BSWmax, 0.0)
       dNS    = dSW/cc%bsw *cc%sapwN
       ! update C and N of sapwood and wood
       cc%bHW   = cc%bHW   + dSW
       cc%bsw   = cc%bsw   - dSW
       cc%sapwN = cc%sapwN - dNS
       cc%woodN = cc%woodN + dNS
       cc%W_stem = cc%W_stem - cc%W_stem * dSW/cc%bSW
       cc%W_dead = cc%W_dead + cc%W_stem * dSW/cc%bSW
    endif
  end associate
end subroutine Sap2HeartWood_fixedHv

!=================================================
! Weng: update soil root area layers, hydraulic variables, 03/29/2022
subroutine Update_hydraulic_vars(cc)
  type(cohort_type), intent(inout) :: cc
  !----------local var ----------
  integer :: j

  associate (sp => spdata(cc%species) )
    cc%Aleaf  = BL2Aleaf(cc%bl,cc)
    cc%lai       = cc%Aleaf/(cc%Acrown *(1.0-sp%f_cGap))
    cc%rootarea  = cc%br * sp%SRA
    do j=1,soil_L
       cc%rootareaL(j) = cc%rootarea * sp%root_frac(j)
    enddo
    ! Plant hydraulics-related variables
    cc%V_leaf = cc%bl/sp%LMA * sp%leafTK ! area * thicknees
    cc%V_stem = cc%bsw/sp%rho_wood ! m^3  ! (cc%bsw+cc%bhw)
    cc%Wmax_l = (cc%V_leaf - cc%bl /rho_cellwall)*rho_H2O ! max leaf water, kg H2O
    cc%Wmax_s = (cc%V_stem - cc%bsw/rho_cellwall)*rho_H2O ! max stem water, kg H2O
    cc%Wmin_l = cc%Wmax_l * exp(sp%psi0_LF*sp%CR_Leaf)
    cc%Wmin_s = cc%Wmax_s * exp(sp%psi0_WD*sp%CR_Wood)
    !cc%H_leaf = (cc%Wmin_l - cc%Wmax_l)/sp%psi0_LF  ! Linearized Leaf Capacitance
    !cc%H_stem = (cc%Wmin_s - cc%Wmax_s)/sp%psi0_WD  ! Linearized Stem capacitance
  end associate
end subroutine Update_hydraulic_vars

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
      W_status = MIN(max(1.0E-6,cc%W_stem),cc%Wmax_s)
      cc%psi_stem = log(W_status/cc%Wmax_s)/sp%CR_Wood
    else
      cc%psi_stem = cc%psi_s0
    endif

    if(cc%Wmax_l > 1.0E-4)then
      W_status = MIN(max(1.0E-6,cc%W_leaf),cc%Wmax_l)
      cc%psi_leaf = log(W_status/cc%Wmax_l)/sp%CR_Leaf
    else
      cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height)
    endif
  end associate
end subroutine Plant_water2psi_exp

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
    ! RAI(i) = cc%rootareaL(i)/cc%Acrown
    ! k_rs(i)= vegn%K_soil(i)*SQRT(RAI(i))/(3.14159*thksl(i))
    k_rs(i)= vegn%K_soil(i)*cc%rootareaL(i) ! per tree
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
real function HT2Lpath(ht) result(L) ! water path lenght from root to leaves
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
       k_rs(i)= vegn%K_soil(i)*cc%rootareaL(i) ! per tree
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

!==============================================================
!============= Vegetation initializations =====================
subroutine initialize_vegn_tile(vegn,nml_file)
   type(vegn_tile_type),intent(inout),pointer :: vegn
   character(len=50),intent(in) :: nml_file

   !--------local vars -------
   type(cohort_type),dimension(:), pointer :: cc, initialCC
   type(cohort_type),pointer :: cp
   integer,parameter :: rand_seed = 86456
   real    :: r
   real    :: btotal
   integer :: nCohorts = 1 ! Randomly generate n Cohorts if not defined
   integer :: i, istat
   integer :: io           ! i/o status for the namelist
   integer :: ierr         ! error code, returned by i/o routines
   integer :: nml_unit

   vegn%age = 0.0 ! Set the tile age as zero.
   !  Read parameters from the parameter file (namelist)
   if(read_from_parameter_file)then
      ! --- Generate cohorts according to "initial_state_nml" ---
      nml_unit = 999
      open(nml_unit, file=nml_file, form='formatted', action='read', status='old')
      read (nml_unit, nml=initial_state_nml, iostat=io, end=20)
      20  close (nml_unit)
      write(*,nml=initial_state_nml)

      ! Initial Soil pools and environmental conditions
      vegn%SOC(4)  = init_fast_soil_C ! kgC m-2
      vegn%SOC(5)  = init_slow_soil_C ! slow soil carbon pool, (kg C/m2)
      vegn%SON(4)  = vegn%SOC(4)/CN0SOM(4)  ! fast soil nitrogen pool, (kg N/m2)
      vegn%SON(5)  = vegn%SOC(5)/CN0SOM(5)  ! slow soil nitrogen pool, (kg N/m2)
      vegn%N_input      = N_input  ! kgN m-2 yr-1, N input to soil
      vegn%mineralN     = init_Nmineral  ! Mineral nitrogen pool, (kg N/m2)
      vegn%previousN    = vegn%mineralN
      !Soil water
      vegn%soiltype = soiltype
      vegn%FLDCAP = FLDCAP
      vegn%WILTPT = WILTPT
      vegn%wcl    = FLDCAP
      vegn%thetaS = 1.0
      call SoilWater_psi_K(vegn)

      ! Initialize plant cohorts
      allocate(cc(1:init_n_cohorts), STAT = istat)
      vegn%cohorts => cc
      vegn%n_cohorts = init_n_cohorts
      cc => null()

      do i=1,init_n_cohorts
         cp => vegn%cohorts(i)
         cp%status  = LEAF_OFF ! ON=1, OFF=0 ! ON
         cp%layer   = 1
         cp%age     = 0
         cp%topyear = 0.0
         cp%species = init_cohort_species(i)
         cp%ccID    = i
         cp%nindivs = init_cohort_nindivs(i) ! trees/m2
         cp%nsc     = init_cohort_nsc(i)
         cp%bsw     = init_cohort_bsw(i)
         cp%bHW     = init_cohort_bHW(i)
         btotal     = cp%bsw + cp%bHW  ! kgC /tree
         call initialize_cohort_from_biomass(cp,btotal,maxval(vegn%psi_soil(:)))
      enddo
      MaxCohortID = cp%ccID
      ! Sorting these cohorts
      call relayer_cohorts(vegn)
   else
     ! ------- Generate one cohort randomly --------
     ! Initialize plant cohorts
      allocate(cc(1:nCohorts), STAT = istat)
      vegn%cohorts => cc
      vegn%n_cohorts = nCohorts
      cc => null()
      r = rand(rand_seed)
      do i=1,nCohorts
         cp => vegn%cohorts(i)
         cp%status  = LEAF_OFF ! ON=1, OFF=0 ! ON
         cp%layer   = 1
         cp%age     = 0
         cp%topyear = 0.0
         cp%species = 3     ! INT(rand()*5)+1
         cp%ccID    = i
         cp%nindivs = 0.001 ! rand()/10. ! trees/m2
         cp%nsc     = 0.005
         cp%bsw     = 0.2   ! kgC /tree
         cp%bHW     = 0.0
         btotal     = cp%bsw + cp%bHW
         call initialize_cohort_from_biomass(cp,btotal,maxval(vegn%psi_soil(:)))
      enddo
      ! Sorting these cohorts
      call relayer_cohorts(vegn)
      ! ID each cohort
      do i=1,nCohorts
         cp => vegn%cohorts(i)
         cp%ccID = MaxCohortID + i
      enddo
      MaxCohortID = cp%ccID
      ! Initial Soil pools and environmental conditions
      vegn%SOC(4) = 0.2 ! kgC m-2
      vegn%SOC(5) = 7.0 ! slow soil carbon pool, (kg C/m2)
      vegn%SON(4) = vegn%SOC(4)/CN0SOM(4)  ! fast soil nitrogen pool, (kg N/m2)
      vegn%SON(5) = vegn%SOC(5)/CN0SOM(5)  ! slow soil nitrogen pool, (kg N/m2)
      vegn%N_input     = N_input  ! kgN m-2 yr-1, N input to soil
      vegn%mineralN    = 0.005  ! Mineral nitrogen pool, (kg N/m2)
      vegn%previousN   = vegn%mineralN
   endif  ! initialization: random or pre-described

   ! tile summary
   call vegn_sum_tile(vegn)
   vegn%initialN0 = vegn%NSN + vegn%SeedN + vegn%leafN +      &
                    vegn%rootN + vegn%SapwoodN + vegn%woodN + &
                    sum(vegn%SON(:)) + vegn%mineralN
   vegn%totN =  vegn%initialN0

   !Setup reserved initial cohorts
   allocate(cc(1:init_n_cohorts), STAT = istat)
   cc = vegn%cohorts
   vegn%initialCC   => cc
   vegn%n_initialCC = init_n_cohorts
   cc => null()

end subroutine initialize_vegn_tile

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
    cc%bl     = 0.0
    cc%br     = cc%bl_max
    cc%nsc    = 2.0 * (cc%bl_max + cc%br_max)
    cc%seedC  = 0.0

    ! N pools
    cc%NSN    = 5.0*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
    cc%leafN  = cc%bl/sp%CNleaf0
    cc%rootN  = cc%br/sp%CNroot0
    cc%sapwN  = cc%bsw/sp%CNsw0
    cc%woodN  = cc%bHW/sp%CNwood0
    cc%seedN  = 0.0

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
    ! Initial psi
    cc%psi_stem = psi_s0
    cc%psi_leaf = cc%psi_stem - HT2MPa(cc%height) ! MPa
    call Update_hydraulic_vars(cc)
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
  cc=>vegn%cohorts
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
  deallocate(vegn%cohorts)
  vegn%cohorts => new
  vegn%n_cohorts = i
  ! update layer fraction for each cohort
  do i=1, vegn%n_cohorts
     vegn%cohorts(i)%layerfrac = vegn%cohorts(i)%nindivs * vegn%cohorts(i)%Acrown
  enddo

end subroutine relayer_cohorts

! ======================================================================
subroutine rank_descending(x,idx)
  ! ranks array x in descending order: on return, idx() contains indices
  ! of elements of array x in descending order of x values. These codes
  ! are from Sergey Malyshev (LM3PPA, Weng et al. 2015 Biogeosciences)

   real,    intent(in)  :: x(:)
   integer, intent(out) :: idx(:)
   integer :: i,n
   integer, allocatable :: t(:)

   n = size(x)
   do i = 1,n
      idx(i) = i
   enddo

   allocate(t((n+1)/2))
   call mergerank(x,idx,n,t)
   deallocate(t)
end subroutine

! =====================================================================
! based on:
! http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
subroutine merge(x,a,na,b,nb,c,nc)
   integer, intent(in) :: na,nb,nc ! Normal usage: NA+NB = NC
   real, intent(in)       :: x(*)
   integer, intent(in)    :: a(na)    ! B overlays C(NA+1:NC)
   integer, intent(in)    :: b(nb)
   integer, intent(inout) :: c(nc)

   integer :: i,j,k

   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) >= x(b(j))) then
         c(k) = a(i) ; i = i+1
      else
         c(k) = b(j) ; j = j+1
      endif
      k = k + 1
   enddo
   do while (i <= na)
      c(k) = a(i) ; i = i + 1 ; k = k + 1
   enddo
end subroutine merge

recursive subroutine mergerank(x,a,n,t)
  integer, intent(in) :: n
  real,    intent(in) :: x(*)
  integer, dimension(n), intent(inout) :: a
  integer, dimension((n+1)/2), intent (out) :: t

  integer :: na,nb
  integer :: v

  if (n < 2) return
  if (n == 2) then
     if ( x(a(1)) < x(a(2)) ) then
        v = a(1) ; a(1) = a(2) ; a(2) = v
     endif
     return
  endif
  na=(n+1)/2
  nb=n-na

  call mergerank(x,a,na,t)
  call mergerank(x,a(na+1),nb,t)

  if (x(a(na)) < x(a(na+1))) then
     t(1:na)=a(1:na)
     call merge(x,t,na,a(na+1),nb,a,n)
  endif
end subroutine mergerank

!============================================================================
! Merge similar cohorts in a tile
subroutine vegn_mergecohorts(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 1.0E-6
  integer :: i,j,k,m

  allocate(cc(vegn%n_cohorts))
  merged(:) = .FALSE.
  k = 0
  do i = 1, vegn%n_cohorts
     if(merged(i)) cycle ! skip cohorts that were already merged
     k = k+1
     cc(k) = vegn%cohorts(i)
     ! try merging the rest of the cohorts into current one
     do j = i+1, vegn%n_cohorts
        if (merged(j)) cycle ! skip cohorts that are already merged
        if (cohorts_can_be_merged(vegn%cohorts(j),cc(k))) then
           call merge_cohorts(vegn%cohorts(j),cc(k))
           merged(j) = .TRUE.

        endif
     enddo
  enddo
  ! at this point, k is the number of new cohorts
  vegn%n_cohorts = k
  deallocate(vegn%cohorts)
  vegn%cohorts=>cc

end subroutine vegn_mergecohorts

! ============================================================================
! kill low density cohorts, a new function seperated from vegn_mergecohorts
! Weng, 2014-07-22
subroutine kill_lowdensity_cohorts(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cp, cc(:) ! array to hold new cohorts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  real, parameter :: mindensity = 0.25E-4
  integer :: i,j,k

 ! calculate the number of cohorts with indivs>mindensity
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs > mindensity) k=k+1
  enddo
  if (k==0)then
     write(*,*)'in kill_lowdensity_cohorts: All cohorts < mindensity, No action!'
     !stop
  endif

  ! exclude cohorts that have low individuals
  if (k>0 .and. k<vegn%n_cohorts)then
     allocate(cc(k))
     j=0
     do i = 1,vegn%n_cohorts
        cp =>vegn%cohorts(i)
        associate(sp=>spdata(cp%species))
        if (cp%nindivs > mindensity) then
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
end subroutine kill_lowdensity_cohorts

! ============================================================================
subroutine merge_cohorts(c1,c2) ! Put c1 into c2
  type(cohort_type), intent(in) :: c1
  type(cohort_type), intent(inout) :: c2

  real :: x1, x2 ! normalized relative weights
  real :: btot, c2H, c2DBH
  integer :: i

  if(c1%nindivs > 0.0 .or. c2%nindivs > 0.0)then
     x1 = c1%nindivs/(c1%nindivs+c2%nindivs)
     x2 = c2%nindivs/(c1%nindivs+c2%nindivs)
     c2%ccID = c1%ccID
     c2DBH = c2%dbh
     c2H   = c2%height

  ! update number of individuals in merged cohort
     c2%nindivs = c1%nindivs + c2%nindivs
     c2%age = x1*c1%age + x2*c2%age
     c2%topyear = x1*c1%topyear + x2*c2%topyear
     c2%C_growth = x1*c1%C_growth + x2*c2%C_growth

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
     c2%sapwN = x1*c1%sapwN + x2*c2%sapwN
     c2%woodN = x1*c1%woodN + x2*c2%woodN
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
     call Update_hydraulic_vars(c2)
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

  endif
end subroutine merge_cohorts

! ============================================================================
function cohorts_can_be_merged(c1,c2); logical cohorts_can_be_merged
   type(cohort_type), intent(in) :: c1,c2

   real, parameter :: mindensity = 1.0E-4
   logical :: sameSpecies,sameLayer,sameSize,sameSizeTree,sameSizeGrass

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
                 ((c1%DBH == c2%DBH).and.c1%age> 2. .and. c2%age>2.)  ! it'll be always true for grasses
   sameSize = sameSizeTree .OR. sameSizeGrass
   cohorts_can_be_merged = sameSpecies .and. sameLayer .and. sameSize
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
  real :: fire_prob, r_fire
  real :: f_bk, r_bk ! coefficient of bark thickness, bark resistance to fire
  real :: s_fire ! Tree's sensitivity to ground surface fire: max 1.0, min 0.0
  real :: flmb_G,flmb_W
  real :: f_grass, f_wood  ! grasses and canopy tree spread probabilities
  real :: mu_fire ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  integer :: i, k

  !  ! Parameters (defined in datatypes.F90 and read in from the namelist file):
  !  envi_fire_prb: Environmental fire occurrence probability, a function of environmental
  !  conditions that can result in fire if fuel is available
  !  (i.e., (match-dropping probability). It should be function of environmental conditions
  ! Vegetation flammability parameters, Ignition_G0, Ignition_W0: Ignition probability
  ! for grasses and woody plants once environmental conditions meet envi_fire_prb
  !  For grasses: Ignition_G0 = 1.0; For woody plants: Ignition_W0 = 0.025
  !  m0_w_fire, m0_g_fire: mortality rates of trees and grasses due to fire
  !  r_BK0: shape parameter ! -480.0  ! for bark resistance, exponential equation,
  !                                  120 --> 0.006 m of bark 0.5 survival
  !  For an old scheme
  !  f_HT0: shape parameter fire resistance (due to growth of bark) as a function of height
  !  h0_escape: tree height that escapes direct burning of grass fires
  !  D_BK0: Bark thickness at half survival rate.
  ! Hoffmann et 2012. shrubs: Y=1.105*X^1.083; trees: Y=0.31*X^1.276 for (Y:mm, X:cm)
  f_bk  = 0.1105
  vegn%treecover = 0.0
  vegn%grasscover = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species))
     if(cc%layer == 1)then
       if(sp%lifeform==0) &
          vegn%grasscover = vegn%grasscover + cc%Acrown*cc%nindivs
       if(sp%lifeform==1) &
          vegn%treecover = vegn%treecover + cc%Acrown*cc%nindivs
     endif
     end associate
  enddo
  ! Fire spread parameters for grasses and trees
  f_grass  = min(1.0, vegn%grasscover)
  f_wood   = min(1.0, vegn%treecover)
  ! Flammability
  flmb_G = Ignition_G0 * f_grass
  flmb_W = Ignition_W0 * f_wood
  fire_prob = 1.-(1.-flmb_G*envi_fire_prb)*(1.-flmb_W*envi_fire_prb)
  ! fire effects on vegetation and soil
  r_fire    = rand(0) !CALL RANDOM_NUMBER(r_fire)
  vegn%C_combusted = 0.0
  if(r_fire < fire_prob)then ! Fire_ON
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      associate ( sp => spdata(cc%species))
      if(sp%lifeform==0) then  ! grasses
         mu_fire = m0_g_fire
      else                     ! trees
         if(r_fire < flmb_W*envi_fire_prb)then
            mu_fire = 0.99 * f_wood   ! tree canopy fire
         else ! grass fire
            cc%D_bark = f_bk * cc%dbh    ! bark thickness,
            s_fire = exp(r_BK0*cc%D_bark)
            !s_fire = 1. - cc%D_bark/(cc%D_bark+D_BK0) !Alternative formulation
            mu_fire = m0_w_fire * s_fire * (1.0-f_wood)
         endif

         !!!!----- Old scheme -------- !!!!!!
         ! s_fire = 1.0 - cc%height/(cc%height+f_HT0)
         !if(cc%height < h0_escape) then ! Short trees
         !   mu_fire = m0_w_fire * s_fire * f_grass
         !   ! Didn't consider the probability of protection of big trees over small trees
         !   ! by excluding grasses and it leads to heterogeinity.
         !else  ! Tall trees
         !   mu_fire = m0_w_fire * s_fire * f_wood
         !endif
      endif

      ! Burned vegetation and soils
      deadtrees = cc%nindivs * MIN(1.0,mu_fire * deltat/seconds_per_year) ! individuals / m2
      ! Carbon and Nitrogen release by burning
      vegn%C_combusted = vegn%C_combusted + &
            (cc%NSC+cc%bl+cc%bsw+cc%bHW+cc%br+cc%seedC)*deadtrees
      vegn%mineralN = vegn%mineralN + &
            (cc%NSN+cc%leafN+cc%rootN+cc%sapwN+cc%woodN+cc%seedN)*deadtrees

      ! Update plant density
      cc%nindivs = cc%nindivs - deadtrees
      end associate
    enddo
    ! Burned litter: 70% of fine litter and 20% of coarse litter are burned
    vegn%C_combusted = vegn%C_combusted + 0.7*vegn%SOC(1)+0.2*vegn%SOC(2)
    vegn%mineralN    = vegn%mineralN    + 0.7*vegn%SON(1)+0.2*vegn%SON(2)
    vegn%SOC(1) = (1.0-0.7)*vegn%SOC(1); vegn%SOC(2) = (1.0-0.2)*vegn%SOC(2)
    vegn%SON(1) = (1.0-0.7)*vegn%SON(1); vegn%SON(2) = (1.0-0.2)*vegn%SON(2)
    write(*,*)"fire, treecover, grasscover", &
        r_fire < fire_prob, vegn%treecover, vegn%grasscover
  endif

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
        vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
        ! remove leaves
        cc%bl = 0.0
     endif
     end associate
     ! Change species
     cc%species = mod(iyears/FREQ,N_SP)+2

 end subroutine vegn_species_switch

!=======================================================================
! Put missing PFTs back with information of initial cohorts (initialCC)
! 10/17/2018, Weng
subroutine vegn_migration (vegn)
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
  write(*,'(A12,16(I4,","))')'current cc',currPFTs(1:nPFTs)

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
  write(*,'(A12,16(I4,","))')'initial cc',initialPFTs(1:nPFT0)

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
      write(*,*)'missing PFTs',missingPFTs(1)
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
               cc = vegn%initialCC(n)
               cc%nindivs = 0.01 ! a small number of individuals
               exit
            endif
        enddo
        totC = cc%bl + cc%br + cc%bsw + cc%bHW + cc%nsc
        totN = cc%leafN + cc%rootN + cc%sapwN + cc%woodN + cc%NSN
        addedC = addedC + cc%nindivs * totC
        addedN = addedN + cc%nindivs * totN
        call setup_seedling(cc,totC,totN)
     enddo

     MaxCohortID = MaxCohortID + newcohorts
     vegn%n_cohorts = k
     ccnew => null()

     ! Make carbon and nitrogen balance
     vegn%SOC(2) = vegn%SOC(2) - min(0.05*vegn%SOC(2),addedC)
     vegn%SON(2) = vegn%SON(2) - min(0.05*vegn%SON(2),addedN)
  endif ! set up newly moved-in cohorts

end subroutine vegn_migration

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
  !      fixedN = fixedN + cc%annualfixedN * cc%Acrown * cc%nindivs
  !enddo

  ! Mineral+fixed N-based LAImax
  ! LAI_fixedN = sp%Nfixrate0 * sp%LMA * sp%CNleaf0 * sp%leafLS / sp%LMA
  ! cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
  vegn%previousN = 0.8 * vegn%previousN + 0.2 * vegn%annualN
  do i=0,MSPECIES
      associate (sp => spdata(i) )

      LAIfixedN  = 0.5 * sp%Nfixrate0 * sp%CNleaf0 * sp%leafLS
      LAImineralN = 0.5*vegn%previousN*sp%CNleaf0*sp%leafLS/sp%LMA
      !LAImineralN = vegn%previousN/(sp%LMA/(sp%CNleaf0*sp%leafLS)+sp%phiRL*sp%alpha_FR/sp%SRA /sp%CNroot0)
      LAI_nitrogen = LAIfixedN + LAImineralN

      spdata(i)%LAImax = MAX(LAImin, MIN(LAI_nitrogen,sp%LAI_light))
      spdata(i)%LAImax_u = MIN(sp%LAImax,1.2)
      end associate
  enddo

  !  ! update the PFTs in the first layer based on fixed N
  !  if(fixedN_based)then ! based on "cc%annualfixedN + vegn%previousN"
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
  !           LAI_nitrogen = 0.5*(vegn%previousN+cc%annualfixedN)*sp%CNleaf0*sp%leafLS/sp%LMA
  !           if(sp%Nfixrate0 > 0.0)
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

!==========================================================================
 ! Hack !!!!!
 subroutine Recover_N_balance(vegn)
   type(vegn_tile_type), intent(inout) :: vegn
      if(abs(vegn%totN-vegn%initialN0)*1000>0.001)then
         vegn%SON(5) = vegn%SON(5) - vegn%totN + vegn%initialN0
         vegn%totN =  vegn%initialN0
      endif

 end subroutine

!========================================================================
!================================Obsolete================================
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
        cc%Q_soil(i) = cc%rootareaL(i) * vegn%K_soil(i) * dpsiSR(i) * step_seconds
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
      ! RAI(i) = cc%rootareaL(i)/cc%Acrown
      ! k_rs(i)= vegn%K_soil(i)*SQRT(RAI(i))/(3.14159*thksl(i))
      k_rs(i)= vegn%K_soil(i)*cc%rootareaL(i) ! per tree
      sumK = sumK + k_rs(i)
      sumPK= sumPK+ k_rs(i) * vegn%psi_soil(i)
    enddo
    psi_ht = HT2MPa(cc%height) ! MPa
    k_stem = cc%Ktrunk * plc_function(cc%psi_stem,sp%psi50_WD,sp%Kexp_WD)
    psi_stem = Max(sp%psi0_WD, (sumPK - Q_stem)/sumK)
    psi_leaf = Max(sp%psi0_LF, psi_stem - Q_leaf/k_stem - psi_ht) ! Ktrunk: (mm/s)/(MPa/m)
  end associate
end subroutine plant_water_potential_equi

! ====================================

end module esdvm
