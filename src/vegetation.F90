module esdvm
 use datatypes
 use soil_mod
 implicit none
 private

! ------ public subroutines ---------
public :: initialize_cohort_from_biomass, initialize_vegn_tile
public :: vegn_phenology,vegn_CNW_budget_fast, vegn_growth_EW,update_layer_LAI
public :: vegn_reproduction, vegn_annualLAImax_update, annual_calls
public :: vegn_starvation, vegn_nat_mortality, vegn_species_switch
public :: relayer_cohorts, vegn_mergecohorts, kill_lowdensity_cohorts

 contains
!=============== ESS subroutines ========================================
!========================================================================
 subroutine vegn_CNW_budget_fast(vegn, forcing)
! hourly carbon, nitrogen, and water dynamics, Weng 2016-11-25
! include Nitrogen uptake and carbon budget
! carbon_gain is calculated here to drive plant growth and reproduciton
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
  tsoil  = forcing%tsoil -273.16 ! degC
  thetaS = (vegn%wcl(2)-WILTPT)/(FLDCAP-WILTPT)

  ! Photosynsthesis
  call vegn_photosynthesis(forcing, vegn)
  !write(*,*)"photosynthesis OK"
  ! Respiration and allocation for growth
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
     ! increment tha cohort age
     cc%age = cc%age + dt_fast_yr

     ! Maintenance respiration
     call plant_respiration(cc,forcing%tair) ! get resp per tree per time step
     cc%resp = cc%resp + cc%resg/24.0 ! put growth respiration into total resp.
     cc%npp  = cc%gpp - cc%resp ! kgC tree-1 hour-1
     ! detach photosynthesis model from plant growth
     !cc%nsc  = cc%nsc + 2.4 * cc%crownarea * dt_fast_yr - cc%resp
     cc%nsc = cc%nsc + cc%npp
     cc%NSN = cc%NSN + cc%fixedN

     END ASSOCIATE
  enddo ! all cohorts

  ! Update soil water 
   !call SoilWaterDynamicsBucket(forcing,vegn)
   call SoilWaterDynamicsLayer(forcing,vegn)
  ! update soil carbon
   call SOMdecomposition(vegn, forcing%tsoil, thetaS)

!! Nitrogen uptake
   call vegn_N_uptake(vegn, forcing%tsoil)

end subroutine vegn_CNW_budget_fast

! ============= Plant physiology ========================================
!========================================================================
! Weng 2017-10-18
! compute stomatal conductance, photosynthesis and respiration
! updates cc%An_op and cc%An_cl, from LM3
subroutine vegn_photosynthesis (forcing, vegn)
  type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

!----- local var --------------
  type(cohort_type),pointer :: cc
  real :: freewater, W_up0 ! potential water uptake by roots, kgH2O m-2 fasttimestep-1
  real :: fWup       ! fraction to the actual soil water
  real :: ws, thetaS ! soil moisture index (0~1)
  real :: dpsiSR     ! pressure difference between soil water and root water, MPa
  real :: rad_top  ! downward radiation at the top of the canopy, W/m2
  real :: rad_net  ! net radiation absorbed by the canopy, W/m2
  real :: Tair, TairK     ! air temperature, degC and degK
  real  :: Dair     ! VPD, pa
  real  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real  :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real  :: p_surf   ! surface pressure, Pa
  real  :: water_supply ! water supply per m2 of leaves
  real  :: fw, fs ! wet and snow-covered fraction of leaves
  real  :: psyn   ! net photosynthesis, mol C/(m2 of leaves s)
  real  :: resp   ! leaf respiration, mol C/(m2 of leaves s)
  real  :: tempLAI,w_scale2, transp ! mol H20 per m2 of leaf per second
  real  :: kappa  ! light extinction coefficient of corwn layers
  real :: f_light(10)=0.0      ! light fraction of each layer
  real :: LAIlayer(10),totCA,f_gap ! additional GPP for lower layer cohorts due to gaps
  integer :: i, layer

  !! Water supply for photosynthesis, Layers
  call water_supply_layer(forcing, vegn)

!! Light supply for photosynthesis
! update accumulative LAI for each corwn layer
  vegn%LAI      = 0.0
  LAIlayer = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     layer = Max (1, Min(cc%layer,9))
     LAIlayer(layer) = LAIlayer(layer) + cc%leafarea * cc%nindivs
     !/(1.0-sp%internal_gap_frac)
     vegn%LAI = vegn%LAI + cc%leafarea  * cc%nindivs
  enddo
  f_gap = 0.1 ! 0.1
  ! Calculate kappa according to sun zenith angle !!
  kappa = 0.75
  ! Light fraction
  f_light = 0.0
  f_light(1) = 1.0
  do i =2, layer !MIN(int(vegn%CAI+1.0),9)
      f_light(i) = f_light(i-1) * &
        (exp(-kappa*LAIlayer(i-1)/(1.-f_gap)) + f_gap)
  enddo

  ! Photosynthesis
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
     if(cc%status == LEAF_ON .and. cc%lai > 0) then
        ! Photosynthesis can be calculated by a photosynthesis model

        ! Convert forcing data
         layer = Max (1, Min(cc%layer,9))
         rad_top = f_light(layer) * forcing%radiation ! downward radiation at the top of the canopy, W/m2
         rad_net = f_light(layer) * forcing%radiation * 0.9 ! net radiation absorbed by the canopy, W/m2
         p_surf  = forcing%P_air  ! Pa
         TairK   = forcing%Tair ! K
         Tair   = forcing%Tair - 273.16 ! degC
         cana_q  = esat(Tair)*forcing%RH/p_surf * mol_h2o/mol_air  ! air specific humidity, kg/kg
         cana_co2= forcing%CO2 ! co2 concentration in canopy air space, mol CO2/mol dry air
        ! recalculate the water supply to mol H20 per m2 of leaf per second
         water_supply = cc%W_supply/(cc%leafarea*step_seconds*mol_h2o) ! mol m-2 leafarea s-1
      
        !call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
        fw = 0.0
        fs = 0.0
        call gs_Leuning(rad_top, rad_net, TairK, cana_q, cc%lai, &
                    p_surf, water_supply, cc%species, sp%pt, &
                    cana_co2, cc%extinct, fs+fw, cc%layer, &
             ! output:
                    psyn, resp,w_scale2,transp )

        ! store the calculated photosynthesis, photorespiration, and transpiration for future use
        ! in growth
        cc%An_op  = psyn  ! molC s-1 m-2 of leaves
        cc%An_cl  = resp  ! molC s-1 m-2 of leaves
        cc%w_scale  = w_scale2
        cc%transp = transp * mol_h2o * cc%leafarea * step_seconds ! Transpiration (kgH2O/(tree hour), Weng, 2017-10-16
        cc%gpp  = (psyn-resp) * mol_C * cc%leafarea * step_seconds ! kgC hour-1 tree-1
     else
        ! no leaves means no photosynthesis and no stomatal conductance either
        cc%An_op  = 0.0;  cc%An_cl  = 0.0
        cc%gpp    = 0.0;  cc%transp = 0.0
        cc%w_scale  = -9999
     endif

     end associate
  enddo ! vegn, go through all cohorts
end subroutine vegn_photosynthesis

! ============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, &
                   p_surf, ws, pft, pt, ca, kappa, leaf_wet, layer, &
                   apot, acl,w_scale2, transp)
  real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
  real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
  real,    intent(in)    :: tl   ! leaf temperature, degK
  real,    intent(in)    :: ea   ! specific humidity in the canopy air, kg/kg
  real,    intent(in)    :: lai  ! leaf area index
  !real,    intent(in)    :: leaf_age ! age of leaf since budburst (deciduos), days
  real,    intent(in)    :: p_surf ! surface pressure, Pa
  real,    intent(in)    :: ws   ! water supply, mol H20/(m2 of leaf s)
  integer, intent(in)    :: pft  ! species
  integer, intent(in)    :: pt   ! physiology type (C3 or C4)
  real,    intent(in)    :: ca   ! concentartion of CO2 in the canopy air space, mol CO2/mol dry air
  real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
  real,    intent(in)    :: leaf_wet ! fraction of leaf that's wet or snow-covered
  integer, intent(in)    :: layer  ! the layer of this canopy
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  !real,    intent(out)   :: gs   ! stomatal conductance, m/s
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
  real,    intent(out)   :: w_scale2,transp  ! transpiration, mol H20/(m2 of leaf s)

  ! ---- local vars     
  ! photosynthesis
  real :: vm;
  real :: kc,ko; ! Michaelis-Menten constants for CO2 and O2, respectively
  real :: ci;
  real :: capgam; ! CO2 compensation point
  real :: f2,f3;
  real :: coef0,coef1;

  real :: Resp;

  ! conductance related
  real :: gs;
  real :: b;
  real :: ds;  ! humidity deficit, kg/kg
  real :: hl;  ! saturated specific humidity at the leaf temperature, kg/kg
  real :: do1;
  
  ! misceleneous
  real :: dum2;
  real, parameter :: light_crit = 0;
  real, parameter :: gs_lim = 0.25;
  real, parameter :: Rgas = 8.314 ! J mol-1 K-1, universal gas constant
  ! new average computations
  real :: lai_eq;
  real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta 
  real :: light_top;
  real :: par_net;
  real :: Ag;
  real :: An;
  real :: Ag_l;
  real :: Ag_rb;
  real :: anbar;
  real :: gsbar;
  real :: w_scale;
  real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa
  ! soil water stress
  real :: Ed,an_w,gs_w;

  b=0.01;
  do1=0.09 ; ! kg/kg
  if (pft < 2) do1=0.15;


  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot;
  par_net   = rad_net*rad_phot;
  
  ! calculate humidity deficit, kg/kg
  call qscomp(tl, p_surf, hl)
  ds = max(hl - ea,0.0)

  ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986

!  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
!  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
!  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));

! corrected by Weng, 2013-01-17
! Weng, 2013-01-10
  ko=0.248    * exp(35948/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
  kc=0.000404 * exp(59356/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
  vm=spdata(pft)%Vmax*exp(24920/Rgas*(1.0/298.2-1.0/tl)) ! / ((layer-1)*1.0+1.0) ! Ea = 33920

  !decrease Vmax due to aging of temperate deciduous leaves 
  !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
!! Turned off by Weng, 2013-02-01, since we can't trace new leaves
!  if (spdata(pft)%leaf_age_tau>0 .and. leaf_age>spdata(pft)%leaf_age_onset) then
!     vm=vm*exp(-(leaf_age-spdata(pft)%leaf_age_onset)/spdata(pft)%leaf_age_tau)
!  endif

  capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

  ! Find respiration for the whole canopy layer
  
!  Resp=spdata(pft)%gamma_resp*vm*lai /((layer-1)*1.0+1.0)  ! Weng, 2013-01-17 add '/ ((layer-1)*1.0+1.0)'

! 2014-09-03, for Nitrogen model: resp = D*(A + B*LMA)
! (A+B*LMA) = LNA, D=Vmax/LNA = 25E-6/0.0012 = 0.02 for a standard deciduous species
!! Leaf resp as a function of nitrogen
!  Resp=spdata(pft)%gamma_resp*0.04*spdata(pft)%LNA  & ! basal rate, mol m-2 s-1
!       * exp(24920/Rgas*(1.0/298.2-1.0/tl))         & ! temperature scaled
!       * lai                                        & ! whole canopy
!       /((layer-1)*1.0+1.0)                         !
!! as a function of LMA
!  Resp=(spdata(pft)%gamma_LNbase*spdata(pft)%LNbase+spdata(pft)%gamma_LMA*spdata(pft)%LMA)  & ! basal rate, mol m-2 s-1
!  Resp=spdata(pft)%gamma_LNbase*(2.5*spdata(pft)%LNA-1.5*spdata(pft)%LNbase)     & ! basal rate, mol m-2 s-1
  Resp= spdata(pft)%gamma_LN/seconds_per_year & ! per seconds, ! basal rate, mol m-2 s-1
       * spdata(pft)%LNA * lai / mol_c    &  ! whole canopy, ! basal rate, mol m-2 s-1
       * exp(24920/Rgas*(1.0/298.2-1.0/tl))        ! temperature scaled
               !                                  &
!       /((layer-1)*1.0+1.0)
! Temperature effects
   Resp=Resp/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
  
  
  ! ignore the difference in concentrations of CO2 near
  !  the leaf and in the canopy air, rb=0.
  Ag_l=0.;
  Ag_rb=0.;
  Ag=0.;
  anbar=-Resp/lai;
  gsbar=b;
  ! find the LAI level at which gross photosynthesis rates are equal
  ! only if PAR is positive
  if ( light_top > light_crit ) then
     if (pt==PT_C4) then ! C4 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        if (ci>capgam) then
           f2=vm;
           f3=18000.0*vm*ci; ! 18000 or 1800?
           dum2=min(f2,f3)
           
           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq = -log(dum2/(kappa*spdata(pft)%alpha_phot*light_top))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           Ag_l   = spdata(pft)%alpha_phot * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           Ag=(Ag_l+Ag_rb)/ &
             ((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
           An=Ag-Resp;
           anbar=An/lai;
     
           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif
        endif ! ci>capgam
     else ! C3 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        coef1=kc*(1.0+0.209/ko);
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        f2=vm*(ci-capgam)/(ci+coef1);
        f3=vm/2.;
        dum2=min(f2,f3);
        if (ci>capgam) then
           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
                       (spdata(pft)%alpha_phot*light_top*kappa))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           Ag_l   = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           Ag=(Ag_l+Ag_rb) /((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
           An=Ag-Resp;
           anbar=An/lai;

           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif
        endif ! ci>capgam
     endif
  endif ! light is available for photosynthesis
  !write(898,'(1(I4,","),10(E10.4,","))') &
  !     layer, light_top, par_net, kappa, lai, lai_eq, ci, capgam, Ag_l, Ag_rb, Ag
  
  an_w=anbar;
  if (an_w > 0.) then
     an_w=an_w*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);
  endif
  gs_w = 1.56 * gsbar *(1-spdata(pft)%wet_leaf_dreg*leaf_wet); !Weng: 1.56 for H2O?
  if (gs_w > gs_lim) then
      if(an_w > 0.) an_w = an_w*gs_lim/gs_w;
      gs_w = gs_lim;
  endif
  ! find water availability diagnostic demand
  Ed = gs_w * ds*mol_air/mol_h2o ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]
  ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
  if (Ed>ws) then
     w_scale=ws/Ed;
     gs_w=w_scale*gs_w;
     if(an_w > 0.0) an_w = an_w*w_scale;
     if(an_w < 0.0.and.gs_w >b) gs_w=b;
  endif
  gs=gs_w
  apot=an_w
  acl=-Resp/lai
  transp = min(ws,Ed) ! mol H20/(m2 of leaf s)
! just for reporting
  w_scale2=min(1.0,ws/Ed)
   ! finally, convert units of stomatal conductance to m/s from mol/(m2 s) by
   ! multiplying it by a volume of a mole of gas
   gs = gs * Rugas * Tl / p_surf
   !write(899, '(25(E12.4,","))') rad_net,par_net,apot*3600*12,acl*3600*12,Ed
end subroutine gs_Leuning

!============================================================================
subroutine plant_respiration(cc, tairK)
  type(cohort_type), intent(inout) :: cc
  real, intent(in) :: tairK ! degK
  !---------local var ---------
  real :: tf,tfs ! thermal inhibition factors for above- and below-ground biomass
  real :: r_leaf, r_stem, r_root
  real :: Acambium  ! cambium area, m2/tree
  ! real :: LeafN     ! leaf nitrogen, kgN/Tree
  real :: fnsc,NSCtarget ! used to regulation respiration rate
  real :: r_Nfix    ! respiration due to N fixation
  
  integer :: sp ! shorthand for cohort species
  sp = cc%species
  ! temperature response function
  tf  = exp(9000.0*(1.0/298.16-1.0/tairK))
!  tfs = thermal_inhibition(tsoil)  ! original
  tfs = tf ! Rm_T_response_function(tsoil) ! Weng 2014-01-14
! With nitrogen model, leaf respiration is a function of leaf nitrogen
  NSCtarget = 3.0 * (cc%bl_max + cc%br_max)
  fnsc = 1.0 ! min(max(0.0,cc%nsc/NSCtarget),1.0)
  Acambium = PI * cc%DBH * cc%height * 1.2

  ! Facultive Nitrogen fixation
  !if(cc%NSN < cc%NSNmax .and. cc%NSC > 0.5 * NSCtarget)then
  !   cc%fixedN = spdata(sp)%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
  !else
  !   cc%fixedN = 0.0 ! spdata(sp)%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 step-1
  !endif

  ! Obligate Nitrogen Fixation
  cc%fixedN = fnsc*spdata(sp)%NfixRate0 * cc%br * tf * dt_fast_yr ! kgN tree-1 hour-1
  r_Nfix    = spdata(sp)%NfixCost0 * cc%fixedN                  ! tree-1 hour-1
  ! LeafN    = spdata(sp)%LNA * cc%leafarea
  r_stem   = fnsc*spdata(sp)%gamma_SW  * Acambium * tf * dt_fast_yr ! kgC tree-1 hour-1
  r_root   = fnsc*spdata(sp)%gamma_FR  * cc%rootN * tf * dt_fast_yr ! root respiration ~ root N
  r_leaf   = fnsc*spdata(sp)%gamma_LN  * cc%leafN * tf * dt_fast_yr  ! tree-1 hour-1

  cc%resp = (r_leaf + r_stem + r_root + r_Nfix) !kgC tree-1 hour-1
  cc%resl = r_leaf !tree-1 hour-1
  cc%resr = r_root + r_Nfix ! tree-1 hour-1

end subroutine plant_respiration

!========= Plant growth ==========================
subroutine fetch_Carbon_for_growth(cc)
!@sum Fetch C from labile C pool according to the demand of leaves and fine roots,
!@+   and the push of labile C pool
!@+   DAILY call.
!@+   added by Weng, 12-06-2016
    implicit none
    type(cohort_type), intent(inout) :: cc
    !------local var -----------
    logical :: woody
    logical :: dormant,growing

    real*8 :: NSCtarget, growthC
    real*8 :: C_push, C_demand
    real*8 :: LFR_rate ! make these two variables to PFT-specific parameters
    real*8 :: bl_max, br_max
    real*8 :: resp_growth

    ! make these two variables to PFT-specific parameters
    LFR_rate = 1.0/16.0 ! filling rate/day

    associate ( sp => spdata(cc%species) )

    NSCtarget = 2.0 * (cc%bl_max + cc%br_max)      ! kgC/tree
    ! Fetch C from labile C pool if it is in the growing seasaon
    if (cc%status == LEAF_ON) then ! growing season
        C_demand = (Max(cc%bl_max - cc%bl,0.0) +   &
                    Max(cc%br_max - cc%br,0.0))* LFR_rate
        C_push = max(cc%nsc-0.5*NSCtarget, 0.0)/(365.0*sp%tauNSC)
        growthC= Max(0.0,MIN(0.025*cc%nsc, C_demand+C_push))

        cc%resg        = 0.25 * growthC
        cc%carbon_gain = 0.75 * growthC    ! kgC/tree

        ! Update NSC pool
        cc%nsc = cc%nsc - cc%carbon_gain ! growthC !

    else ! non-growing season
        cc%carbon_gain = 0.0
        cc%resg    = 0.0
    endif
    end associate

 end subroutine fetch_Carbon_for_growth

! ============================================================================

subroutine vegn_growth_EW(vegn)
! updates cohort biomass pools, LAI, and height using accumulated 
! carbon_gain and bHW_gain
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(cohort_type), pointer :: cc    ! current cohort
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: dB_LRS, G_LFR  ! amount of carbon spent on leaf and root growth
  real :: dSeed ! allocation to seeds, Weng, 2016-11-26
  real :: dBL, dBR ! tendencies of leaf and root biomass, kgC/individual
  real :: dBSW ! tendency of sapwood biomass, kgC/individual
  real :: dBHW ! tendency of wood biomass, kgC/individual
  real :: dDBH ! tendency of breast height diameter, m
  real :: dCA ! tendency of crown area, m2/individual
  real :: dHeight ! tendency of vegetation height
  real :: dNS    ! Nitrogen from SW to HW
  real :: sw2nsc = 0.0 ! conversion of sapwood to non-structural carbon
  real :: b,BL_u,BL_c
  real :: alphaBL, alphaBR
  real :: DBHtp
  real :: N_supply, N_demand,fNr,Nsupplyratio,extrasapwN
  integer :: i,j

  ! Turnover of leaves and fine roots
  call vegn_tissue_turnover(vegn)

  !Allocate C_gain to tissues
  DBHtp = 0.8
  fNr   = 0.25
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
 !    call biomass_allocation(cc)
     associate (sp => spdata(cc%species)) ! F2003
     if (cc%status == LEAF_ON) then
        ! Get carbon from NSC pool
        call fetch_Carbon_for_growth(cc) ! Weng, 2017-10-19

        ! Allocate carbon to the plant pools
        ! calculate the carbon spent on growth of leaves and roots
        G_LFR = max(0.0, min(cc%bl_max+cc%br_max-cc%bl-cc%br,  &
                            (1.-Wood_fract_min)*cc%carbon_gain))
        ! and distribute it between roots and leaves
        dBL  = min(G_LFR, max(0.0, &
          (G_LFR*cc%bl_max + cc%bl_max*cc%br - cc%br_max*cc%bl)/(cc%bl_max + cc%br_max) &
          ))
        dBR  = G_LFR - dBL
        ! calculate carbon spent on growth of sapwood growth
        if(cc%layer == 1 .AND. cc%age > sp%maturalage)then
            dSeed=      sp%v_seed * (cc%carbon_gain - G_LFR)
            dBSW = (1.0-sp%v_seed)* (cc%carbon_gain - G_LFR)
        else
            dSeed= 0.0
            dBSW = cc%carbon_gain - G_LFR
        endif
!       Specially for grasses, temporary
        if(sp%lifeform ==0 )then
            dSeed = dSeed + 0.15*G_LFR
            G_LFR = 0.85 * G_LFR
            dBR   = 0.85 * dBR
            dBL   = 0.85 * dBL
        endif
!!       Nitrogen effects on allocations between wood and leaves+roots
!
!!       Nitrogen demand by leaves, roots, and seeds (Their C/N ratios are fixed.)
        N_demand = dBL/sp%CNleaf0 + dBR/sp%CNroot0 + dSeed/sp%CNseed0
!!       Nitrogen available for all tisues, including wood
        N_supply= MAX(0.0, fNr*cc%NSN)
!!       same ratio reduction for leaf, root, and seed if(N_supply < N_demand)
        IF(N_demand > N_supply )then
            dB_LRS = dBL+dBR+dSeed
            Nsupplyratio = N_supply / N_demand
            dBR  =  Nsupplyratio * dBR
            dBL  =  Nsupplyratio * dBL
            dSeed=  Nsupplyratio * dSeed
            dBSW =  dBSW + (1.0 - Nsupplyratio) * dB_LRS ! dB_LRS: = dBL+dBR+dSeed, original guess)
            N_demand = N_supply ! dBL/sp%CNleaf0 + dBR/sp%CNroot0 + dSeed/sp%CNseed0
        ENDIF
!       update biomass pools
        cc%bl     = cc%bl    + dBL
        cc%br     = cc%br    + dBR
        cc%bsw    = cc%bsw   + dBSW
        cc%seedC  = cc%seedC + dSeed
!!      update nitrogen pools, Nitrogen allocation
        cc%NSN   = cc%NSN   - N_supply
        cc%leafN = cc%leafN + dBL   /sp%CNleaf0
        cc%rootN = cc%rootN + dBR   /sp%CNroot0
        cc%seedN = cc%seedN + dSeed /sp%CNseed0
        cc%sapwN = cc%sapwN + MAX((N_supply - N_demand),0.0)
!       Return excessiive Nitrogen in SW back to NSN
        if(cc%sapwN > cc%bsw/sp%CNsw0)then
           extrasapwN = cc%sapwN - cc%bsw/sp%CNsw0
           cc%NSN     = cc%NSN   + extrasapwN ! MAX(0.0, cc%sapwN - cc%bsw/sp%CNsw0)
           cc%sapwN   = cc%sapwN - extrasapwN ! MAX(0.0, cc%sapwN - cc%bsw/sp%CNsw0)
        endif

!       accumulated C allocated to leaf, root, and wood
        cc%NPPleaf = cc%NPPleaf + dBL
        cc%NPProot = cc%NPProot + dBR
        cc%NPPwood = cc%NPPwood + dBSW

!       update breast height diameter given increase of bsw
        dDBH   = dBSW / (sp%thetaBM * sp%alphaBM * cc%DBH**(sp%thetaBM-1))
        dHeight= sp%thetaHT * sp%alphaHT * cc%DBH**(sp%thetaHT-1) * dDBH
        dCA    = sp%thetaCA * sp%alphaCA * cc%DBH**(sp%thetaCA-1) * dDBH
!       update plant architecture
        cc%DBH       = cc%DBH       + dDBH
        cc%height    = cc%height    + dHeight
        cc%crownarea = cc%crownarea + dCA
        cc%leafarea  = leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
        cc%lai       = cc%leafarea/cc%crownarea !(cc%crownarea *(1.0-sp%internal_gap_frac))
        call rootarea_and_verticalprofile(cc)
!       convert sapwood to heartwood ! Nitrogen from sapwood to heart wood
        if(sp%lifeform>0)then
           CSAsw  = cc%bl_max/sp%LMA * sp%phiCSA * cc%height ! with Plant hydraulics, Weng, 2016-11-30
           CSAtot = 0.25 * PI * cc%DBH**2
           CSAwd  = max(0.0, CSAtot - CSAsw)
           DBHwd  = 2*sqrt(CSAwd/PI)
           BSWmax = sp%alphaBM * (cc%DBH**sp%thetaBM - DBHwd**sp%thetaBM)
           dBHW   = max(cc%bsw - BSWmax, 0.0)
           dNS   = dBHW/cc%bsw *cc%sapwN
           ! update C and N of sapwood and wood
           cc%bHW   = cc%bHW   + dBHW
           cc%bsw   = cc%bsw   - dBHW
           cc%sapwN = cc%sapwN - dNS
           cc%woodN = cc%woodN + dNS
        endif

!       update bl_max and br_max daily
        BL_c = sp%LMA * sp%LAImax * cc%crownarea * &
               (1.0-sp%internal_gap_frac)
        BL_u = sp%LMA*cc%crownarea*(1.0-sp%internal_gap_frac)* &
                    sp%underLAImax
        if (cc%layer == 1) cc%topyear = cc%topyear + 1.0 /365.0
        if (cc%layer > 1 .and. cc%firstlayer == 0) then ! changed back, Weng 2014-01-23
            cc%bl_max = BL_u
!           Keep understory tree's root low and constant
            cc%br_max = 0.8*cc%bl_max/(sp%LMA*sp%SRA) ! sp%phiRL
        else
            cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
            cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
        endif
        ! Grasses have the saem bl_max regardless of their layer position
        if(sp%lifeform == 0) then
           cc%bl_max = BL_c
           cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
        endif ! for grasses
     elseif(cc%status == LEAF_OFF .and. cc%carbon_gain>0.)then
        cc%nsc = cc%nsc + cc%carbon_gain
     endif ! "cc%status == LEAF_ON"
     ! reset carbon acculmulation terms
     cc%carbon_gain = 0
  end associate ! F2003
  enddo
  cc => null()

end subroutine vegn_growth_EW ! daily

!=================================================
! Weng, 2017-10-26
subroutine update_layer_LAI(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

!----- local var --------------
  type(cohort_type),pointer :: cc
  integer :: i, layer

! update accumulative LAI for each corwn layer
  vegn%LAI      = 0.0
  vegn%LAIlayer = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     layer = Max (1, Min(cc%layer,9)) ! between 1~9
     vegn%LAIlayer(layer) = vegn%LAIlayer(layer) + cc%leafarea * cc%nindivs !/(1.0-sp%internal_gap_frac)
  enddo
 end subroutine update_layer_LAI

!=================================================
! Weng: partioning root area into layers, 10-24-2017
subroutine rootarea_and_verticalprofile(cc)
  type(cohort_type), intent(inout) :: cc
  !----------local var ----------
  integer :: j

  associate (sp => spdata(cc%species) )
  cc%rootarea  = cc%br * sp%SRA
  do j=1,max_lev
     cc%rootareaL(j) = cc%rootarea * sp%root_frac(j)
  enddo
  end associate
 end subroutine rootarea_and_verticalprofile
!============================================================================
subroutine vegn_phenology(vegn,doy) ! daily step
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: doy

  ! ---- local vars
  type(cohort_type), pointer :: cc
  integer :: i,j
  real    :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real    :: stem_fall, grassdensity   ! for grasses only
  real    :: dAleaf, dBL, dBR, dNL, dNR, dNS      ! per day
  real    :: leaf_fall_rate, root_mort_rate      ! per day
  real    :: BL_u,BL_c
  real    :: retransN  ! retranslocation coefficient of Nitrogen
  logical :: cc_firstday = .false.
  logical :: growingseason
  logical :: TURN_ON_life, TURN_OFF_life

  retransN = 0.5
  leaf_fall_rate = 0.075; root_mort_rate = 0.0

  vegn%litter = 0   ! daily litter

  ! update vegn GDD and tc_pheno
  vegn%gdd      = vegn%gdd + max(0.0, vegn%tc_daily - 278.15)
  vegn%tc_pheno = vegn%tc_pheno * 0.85 + vegn%Tc_daily * 0.15

! ON and OFF of phenology: change the indicator of growing season for deciduous
  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! update GDD for each cohort
     cc%gdd = cc%gdd + max(0.0, vegn%tc_daily - 278.15) ! GDD5
     
     associate (sp => spdata(cc%species) )
!    for evergreen
     if(sp%phenotype==1 .and. cc%status==LEAF_OFF) cc%status=LEAF_ON
!    for deciduous and grasses
     TURN_ON_life = (sp%phenotype == 0              .and. &
                    cc%status    == LEAF_OFF       .and. &
                    cc%gdd        > sp%gdd_crit    .and. &
                    vegn%tc_pheno > sp%tc_crit_on) .and. &
             (sp%lifeform /=0 .OR.(sp%lifeform ==0 .and.cc%layer==1))

     TURN_OFF_life = (sp%phenotype  == 0 .and.       &
                    cc%status     == LEAF_ON .and. &
                    vegn%tc_pheno < sp%tc_crit)

     cc_firstday = .false.
     if(TURN_ON_life)then
         cc%status = LEAF_ON ! Turn on a growing season
         cc_firstday = .True.

     elseif(TURN_OFF_life )then
        cc%status = LEAF_OFF  ! Turn off a growing season
        cc%gdd   = 0.0        ! Start to counting a new cycle of GDD
        vegn%gdd = 0.0
     endif

!    setup grass density at the first day of a growing season
     if(cc_firstday .and. sp%lifeform ==0 .and. &
        cc%age>3 .and.(cc%nsc>2*sp%seedlingsize .or. cc%nsc<1*sp%seedlingsize) ) then
!        reset grass density and size for perenials
         grassdensity = cc%nindivs ! last year's density
         cc%nindivs = grassdensity*cc%nsc /(2*sp%seedlingsize)
         cc%br  = grassdensity/cc%nindivs * cc%br
         cc%nsc = grassdensity/cc%nindivs * cc%nsc
         cc%bsw = 0.2 *sp%seedlingsize  ! for setting up a initial size
         cc%nsc = cc%nsc - cc%bsw

        call rootarea_and_verticalprofile(cc)

!!        Nitrogen pools
         cc%rootN = grassdensity/cc%nindivs * cc%rootN
         cc%NSN    = grassdensity/cc%nindivs * cc%NSN
         cc%woodN   = cc%bsw/sp%CNsw0
         cc%NSN    = cc%NSN - cc%woodN

         call init_cohort_allometry(cc)
     endif

!    End a growing season: leaves fall for deciduous
     if(cc%status == LEAF_OFF .AND. cc%bl > 0.0)then
        dBL = cc%bl ! min(leaf_fall_rate * cc%bl_max, cc%bl)
        dBR = 0.0   ! min( root_mort_rate * cc%br_max, cc%br)  ! Just for test: keep roots
        stem_fall = 0.0 ! trees
        if(sp%lifeform==0)then  ! grasses
            stem_fall = MIN(1.0,dBL/cc%bl) * cc%bsw
            dNS = MIN(1.0,dBL/cc%bl) * cc%woodN
        else
            stem_fall = 0.0
        endif
        ! Nitrogen out
        dNL = dBL/sp%CNleaf0
        dNR = dBR/sp%CNroot0

        dAleaf = leaf_area_from_biomass(dBL,cc%species,cc%layer,cc%firstlayer)


!       Retranslocation to NSC and NSN
        cc%nsc = cc%nsc + l_fract  * (dBL + dBR + stem_fall)
        cc%NSN = cc%NSN + retransN * (dNL + dNR + dNS)
!       update plant pools
        cc%bl    = cc%bl  - dBL
        cc%br    = cc%br  - dBR
        cc%bsw   = cc%bsw - stem_fall
        !cc%bsw   = cc%bsw - stem_fall ! for grasses
        cc%leafN = cc%leafN - dNL ! cc%leafN * (1. - dBL    / cc%bl)
        cc%rootN = cc%rootN - dNR ! cc%rootN * (1.- dBR / cc%br)
        cc%woodN  = cc%woodN  - dNS
!       update NPP for leaves, fine roots, and wood

        cc%NPPleaf = cc%NPPleaf - l_fract * dBL
        cc%NPProot = cc%NPProot - l_fract * dBR
        cc%NPPwood = cc%NPPwood - l_fract * stem_fall
        ! cc%NPPwood = cc%NPPwood - l_fract * stem_fall  ! for grasses
        cc%leafarea= leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
        cc%lai     = cc%leafarea/(cc%crownarea *(1.0-sp%internal_gap_frac))

!       put C and N into soil pools:  Substraction of C and N from leaf and root pools
        loss_coarse  = cc%nindivs * (stem_fall+dBL - dAleaf * LMAmin)
        loss_fine    = cc%nindivs * (dBR + dAleaf * LMAmin) ! + stem_fall ! for grasses
        lossN_coarse = cc%nindivs * (dNS + dNL - dAleaf * sp%LNbase)
        lossN_fine   = cc%nindivs * (dNR + dAleaf * sp%LNbase)  !  + stem_fall/sp%CNwood
        vegn%metabolicL = vegn%metabolicL + (1.-l_fract) *  &
                         (fsc_fine * loss_fine + fsc_wood * loss_coarse)
        vegn%structuralL = vegn%structuralL + (1.-l_fract) *     &
                         ((1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse)

!       Nitrogen to soil SOMs
        vegn%metabolicN  = vegn%metabolicN + (1.-retransN) *    &
                          (fsc_fine * lossN_fine + fsc_wood * lossN_coarse)
        vegn%structuralN = vegn%structuralN + (1.-retransN) *    &
                          ((1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse)

!       annual N from plants to soil
        vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
     endif
     end associate
     cc => null()
  enddo
end subroutine vegn_phenology

!============================================================================
!------------------------Mortality------------------------------------
subroutine vegn_nat_mortality (vegn, deltat)
! TODO: update background mortality rate as a function of wood density (Weng, Jan. 07 2017)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! time since last mortality calculations, s

  ! ---- local vars
  type(cohort_type), pointer :: cc => null()
  type(spec_data_type),   pointer :: sp
  real :: loss_fine, loss_coarse
  real :: lossN_fine,lossN_coarse
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  real :: DBHtp, tmp
  integer :: i, k

  real, parameter :: min_nindivs = 1e-5 ! 2e-15 ! 1/m. If nindivs is less than this number, 
  ! then the entire cohort is killed; 2e-15 is approximately 1 individual per Earth 
  ! surface area
  !write(*,*)'total cohorts:', vegn%n_cohorts
  !write(*,'(a81)')'i,PFT,layer,density,layerfrac,dDBH,dbh,height,Acrown,wood,nsc,NPPL,NPPW,aGPP,aNPP'
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species))

     ! mortality rate can be a function of growth rate, age, and environmental
     ! conditions. Here, we only used two constants for canopy layer and under-
     ! story layer (mortrate_d_c and mortrate_d_u)
     ! for trees
     if(cc%layer > 1) then
            tmp = (1 + A_mort*exp(B_mort*(DBH_mort-cc%dbh)) &
                       /(1.0 + exp(B_mort*(DBH_mort-cc%dbh))))
            deathrate = spdata(cc%species)%mortrate_d_u * tmp
     else
            deathrate = spdata(cc%species)%mortrate_d_c !sp%mortrate_d_c

     endif
     deadtrees = cc%nindivs * (1.0-exp(-deathrate*deltat/seconds_per_year)) ! individuals / m2
     ! add dead C from leaf and root pools to fast soil carbon
     loss_coarse  = deadtrees * (cc%bHW + cc%bsw   + cc%bl - cc%leafarea*LMAmin)
     loss_fine    = deadtrees * (cc%nsc + cc%seedC + cc%br + cc%leafarea*LMAmin)
     lossN_coarse = deadtrees * (cc%woodN + cc%sapwN + cc%leafN - sp%LNbase*cc%leafarea)
     lossN_fine   = deadtrees * (cc%rootN + cc%seedN + cc%NSN   + sp%LNbase*cc%leafarea)

     vegn%metabolicL  = vegn%metabolicL +    fsc_fine *loss_fine +    fsc_wood *loss_coarse
     vegn%structuralL = vegn%structuralL + (1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse
     ! Nitrogen to soil pools

     vegn%metabolicN = vegn%metabolicN +  fsc_fine * lossN_fine +   &
                                          fsc_wood * lossN_coarse
     vegn%structuralN = vegn%structuralN +(1.-fsc_fine) * lossN_fine +   &
                                          (1.-fsc_wood) * lossN_coarse

!    annual N from plants to soil
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
     ! Update plant density
     cc%nindivs = cc%nindivs-deadtrees
     end associate
  enddo
end subroutine vegn_nat_mortality

!========================================================================
! Starvation due to low NSC
subroutine vegn_starvation (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars --------
  real :: deathrate ! mortality rate, 1/year
  real :: deadtrees ! number of trees that died over the time step
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse
  integer :: i, k
  type(cohort_type), pointer :: cc
  type(cohort_type), dimension(:),pointer :: ccold, ccnew

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species)  )

!   Mortality due to starvation
    deathrate = 0.0
!   if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max .OR.(cc%layer >1 .and. sp%lifeform ==0)) then
    if (cc%bsw<0 .or. cc%nsc < 0.00001*cc%bl_max) then
         deathrate = 1.0
         deadtrees = cc%nindivs * deathrate !individuals / m2

         ! Carbon and Nitrogen from plants to soil pools
         loss_coarse  = deadtrees * (cc%bHW + cc%bsw   + cc%bl - cc%leafarea*LMAmin)
         loss_fine    = deadtrees * (cc%nsc + cc%seedC + cc%br + cc%leafarea*LMAmin)
         lossN_coarse = deadtrees * (cc%woodN + cc%sapwN + cc%leafN - sp%LNbase*cc%leafarea)
         lossN_fine   = deadtrees * (cc%rootN + cc%seedN + cc%NSN   + sp%LNbase*cc%leafarea)

         vegn%metabolicL  = vegn%metabolicL + fsc_fine *loss_fine + fsc_wood *loss_coarse
         vegn%structuralL = vegn%structuralL + (1.0-fsc_fine)*loss_fine + (1.0-fsc_wood)*loss_coarse

         vegn%metabolicN = vegn%metabolicN + &
                           fsc_fine *lossN_fine +    fsc_wood *lossN_coarse
         vegn%structuralN = vegn%structuralN + &
                       (1.-fsc_fine)*lossN_fine +(1.-fsc_wood)*lossN_coarse

         ! annual N from plants to soil
         vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse

!        update cohort individuals
         cc%nindivs = cc%nindivs*(1.0-deathrate)
     else
         deathrate = 0.0
     endif

     end associate
  enddo
end subroutine vegn_starvation

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
  integer :: newcohorts, matchflag, nPFTs ! number of new cohorts to be created
  integer :: nCohorts, istat
  integer :: i, j, k ! cohort indices

! Looping through all reproductable cohorts and Check if reproduction happens
  reproPFTs = -999 ! the code of reproductive PFT
  vegn%totseedC = 0.0
  vegn%totseedN = 0.0
  vegn%totNewCC = 0.0
  vegn%totNewCN = 0.0
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
               vegn%totSeedC = vegn%totSeedC + cc%seedC  * cc%nindivs
               vegn%totSeedN = vegn%totSeedN + cc%seedN  * cc%nindivs
               cc%seedC = 0.0
               cc%seedN = 0.0

               matchflag = 1
               exit
           endif
        enddo
        if(matchflag==0)then ! when it is a new PFT, put it to the next place
            nPFTs            = nPFTs + 1 ! update the number os reproducible PFTs
            reproPFTs(nPFTs) = cc%species ! PFT number
            seedC(nPFTs)     = cc%seedC * cc%nindivs ! seed carbon
            seedN(nPFTs)     = cc%seedN * cc%nindivs ! seed nitrogen
            vegn%totSeedC = vegn%totSeedC + cc%seedC  * cc%nindivs
            vegn%totSeedN = vegn%totSeedN + cc%seedN  * cc%nindivs
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

     deallocate (ccold)

     ! set up new cohorts
     k = vegn%n_cohorts
     do i = 1, newcohorts
        k = k+1 ! increment new cohort index
        cc => vegn%cohorts(k)
        ! Give the new cohort an ID
        cc%ccID = MaxCohortID + i
        ! update child cohort parameters
        associate (sp => spdata(reproPFTs(i))) ! F2003
        ! density
        cc%nindivs = seedC(i)/sp%seedlingsize

        cc%species = reproPFTs(i)
        cc%status  = LEAF_OFF
        cc%firstlayer = 0
        cc%topyear = 0.0
        cc%age     = 0.0

        ! Carbon pools
        cc%bl      = 0.0 * sp%seedlingsize
        cc%br      = 0.1 * sp%seedlingsize
        cc%bsw     = 0.4 * sp%seedlingsize
        cc%bHW     = 0.0 * sp%seedlingsize
        cc%seedC   = 0.0
        cc%nsc     = sp%seedlingsize - cc%bsw !
        call rootarea_and_verticalprofile(cc)

!!      Nitrogen pools
        cc%leafN  = cc%bl/sp%CNleaf0
        cc%rootN  = cc%br/sp%CNroot0
        cc%sapwN  = cc%bsw/sp%CNsw0
        cc%woodN  = cc%bHW/sp%CNwood0
        cc%seedN  = 0.0
        cc%NSN    = sp%seedlingsize * seedN(i) / seedC(i) -  &
                    (cc%leafN + cc%rootN + cc%sapwN + cc%woodN)

        vegn%totNewCC = vegn%totNewCC + cc%nindivs*(cc%bl + cc%br + cc%bsw + cc%bHW + cc%nsc)
        vegn%totNewCN = vegn%totNewCN + cc%nindivs*(cc%leafN + cc%rootN + cc%sapwN + cc%woodN + cc%NSN)

        call init_cohort_allometry(cc)
!!        !! seeds fail
        !cc%nindivs = cc%nindivs * sp%prob_g * sp%prob_e
!!       put failed seeds to soil carbon pools
!        failed_seeds = 0.0 ! (1. - sp%prob_g*sp%prob_e) * seedC(i)!!

!        vegn%litter = vegn%litter + failed_seeds
!        vegn%metabolicL = vegn%metabolicL +        fsc_fine *failed_seeds
!        vegn%structuralL = vegn%structuralL + (1.0 - fsc_fine)*failed_seeds

!!      Nitrogen of seeds to soil SOMs
!        N_failedseed= 0.0 ! (1.-sp%prob_g*sp%prob_e)   * seedN(i)
!        vegn%metabolicN  = vegn%metabolicN   +        fsc_fine * N_failedseed
!        vegn%structuralN = vegn%structuralN  + (1.0 - fsc_fine)* N_failedseed

!       annual N from plants to soil
 !   vegn%N_P2S_yr = vegn%N_P2S_yr + N_failedseed

        end associate   ! F2003
     enddo
     MaxCohortID = MaxCohortID + newcohorts
     vegn%n_cohorts = k
     ccnew => null()
     call zero_diagnostics(vegn)
  endif ! set up new born cohorts

end subroutine vegn_reproduction

! ============================================================================
function cohort_can_reproduce(cc); logical cohort_can_reproduce
  type(cohort_type), intent(in) :: cc

  associate (sp => spdata(cc%species) )! F2003
  cohort_can_reproduce = (cc%layer == 1 .and. &
                          cc%age   > sp%maturalage.and. &
                          cc%seedC > sp%seedlingsize .and. &
                          cc%seedN > sp%seedlingsize/sp%CNseed0)
  end associate

end function


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
        loss_coarse  = cc%nindivs * (cc%bl - cc%leafarea*LMAmin)
        loss_fine    = cc%nindivs *  cc%leafarea*LMAmin
        lossN_coarse = cc%nindivs * (cc%leafN - sp%LNbase*cc%leafarea)
        lossN_fine   = cc%nindivs *  sp%LNbase*cc%leafarea
        ! Carbon to soil pools
        vegn%metabolicL  = vegn%metabolicL  + fsc_fine *loss_fine + &
                                              fsc_wood *loss_coarse
        vegn%structuralL = vegn%structuralL + (1.0-fsc_fine)*loss_fine + &
                                              (1.0-fsc_wood)*loss_coarse
        ! Nitrogen to soil pools
        vegn%metabolicN = vegn%metabolicN + fsc_fine  *lossN_fine +   &
                                        fsc_wood *lossN_coarse
        vegn%structuralN = vegn%structuralN +(1.-fsc_fine) *lossN_fine +   &
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

! ============================================================================
! Arrange crowns into canopy layers according to their height and crown areas.
subroutine relayer_cohorts (vegn)
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

!  rand_sorting = .TRUE. ! .False.
  
  ! rank cohorts in descending order by height. For now, assume that they are 
  ! in order
  N0 = vegn%n_cohorts; cc=>vegn%cohorts
  call rank_descending(cc(1:N0)%height,idx)
  
  ! calculate max possible number of new cohorts : it is equal to the number of
  ! old cohorts, plus the number of layers -- since the number of full layers is 
  ! equal to the maximum number of times an input cohort can be split by a layer 
  ! boundary.
  N1 = vegn%n_cohorts + int(sum(cc(1:N0)%nindivs*cc(1:N0)%crownarea))
  allocate(new(N1))

  ! copy cohort information to the new cohorts, splitting the old cohorts that 
  ! stride the layer boundaries
  i = 1 ; k = 1 ; L = 1 ; frac = 0.0 ; nindivs = cc(idx(k))%nindivs
  do 
     new(i)         = cc(idx(k))
     new(i)%nindivs = min(nindivs,(layer_vegn_cover-frac)/cc(idx(k))%crownarea)
     new(i)%layer   = L
     if (L==1) new(i)%firstlayer = 1
!    if (L>1)  new(i)%firstlayer = 0  ! switch off "push-down effects"
     frac = frac+new(i)%nindivs*new(i)%crownarea
     nindivs = nindivs - new(i)%nindivs
     
     if (abs(nindivs*cc(idx(k))%crownarea)<tolerance) then
       new(i)%nindivs = new(i)%nindivs + nindivs ! allocate the remainder of individuals to the last cohort
       if (k==N0) exit ! end of loop
       k = k+1 ; nindivs = cc(idx(k))%nindivs  ! go to the next input cohort
     endif
     
     if (abs(layer_vegn_cover - frac)<tolerance) then
       L = L+1 ; frac = 0.0              ! start new layer
     endif
!     write(*,*)i, new(i)%layer
     i = i+1
  enddo
  
  ! replace the array of cohorts
  deallocate(vegn%cohorts)
  vegn%cohorts => new ; vegn%n_cohorts = i
  ! update layer fraction for each cohort
  do i=1, vegn%n_cohorts
     vegn%cohorts(i)%layerfrac = vegn%cohorts(i)%nindivs * vegn%cohorts(i)%crownarea
  enddo

end subroutine relayer_cohorts

! ============================================================================

 subroutine vegn_tissue_turnover(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  !-------local var
  type(cohort_type), pointer :: cc    ! current cohort
  real :: loss_coarse, loss_fine, lossN_coarse, lossN_fine
  real :: dBL, dBR  ! leaf and fine root carbon tendencies
  real :: dNL, dNR  ! leaf and fine root nitrogen tendencies
  real :: dAleaf ! leaf area decrease due to dBL
  real :: retransN = 0.5  ! retranslocation coefficient of Nitrogen
  integer :: i

  ! update plant carbon and nitrogen for all cohorts
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species) )
!    Turnover of leaves and roots regardless of the STATUS of leaf
!    longevity. Deciduous: 0; Evergreen 0.035/LMa
!    root turnover
     dBL = cc%bl    * sp%alpha_L  /days_per_year
     dBR = cc%br    * sp%alpha_FR /days_per_year
     dNL = cc%leafN * sp%alpha_L  /days_per_year
     dNR = cc%rootN * sp%alpha_FR /days_per_year
     dAleaf = leaf_area_from_biomass(dBL,cc%species,cc%layer,cc%firstlayer)

!    Retranslocation to NSC and NSN
     cc%nsc = cc%nsc + l_fract  * (dBL + dBR)
     cc%NSN = cc%NSN + retransN * (dNL + dNR)
!    update plant pools
     cc%bl    = cc%bl    - dBL
     cc%br    = cc%br    - dBR
     cc%leafN = cc%leafN - dNL
     cc%rootN = cc%rootN - dNR

!    update leaf area and LAI
     cc%leafarea= leaf_area_from_biomass(cc%bl,cc%species,cc%layer,cc%firstlayer)
     cc%lai     = cc%leafarea/(cc%crownarea *(1.0-sp%internal_gap_frac))

!    update NPP for leaves, fine roots, and wood
     cc%NPPleaf = cc%NPPleaf - l_fract * dBL
     cc%NPProot = cc%NPProot - l_fract * dBR

!    put C and N into soil pools
     loss_coarse  = cc%nindivs * (dBL - dAleaf * LMAmin)
     loss_fine    = cc%nindivs * (dBR + dAleaf * LMAmin)
     lossN_coarse = cc%nindivs * (dNL - dAleaf * sp%LNbase)
     lossN_fine   = cc%nindivs * (dNR + dAleaf * sp%LNbase)

     vegn%metabolicL = vegn%metabolicL + (1.-l_fract) *  &
                        (fsc_fine * loss_fine + fsc_wood * loss_coarse)
     vegn%structuralL = vegn%structuralL + (1.-l_fract) *     &
                         ((1.-fsc_fine)*loss_fine + (1.-fsc_wood)*loss_coarse)

!    Nitrogen to soil SOMs
     vegn%metabolicN  = vegn%metabolicN + (1.-retransN) *    &
                          (fsc_fine * lossN_fine + fsc_wood * lossN_coarse)
     vegn%structuralN = vegn%structuralN + (1.-retransN) *    &
                          ((1.-fsc_fine) * lossN_fine + (1.-fsc_wood) * lossN_coarse)

!    annual N from plants to soil
     vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse

    END ASSOCIATE
  enddo

end subroutine vegn_tissue_turnover
!=====================================================
! Weng, 2016-11-28
subroutine vegn_N_uptake(vegn, tsoil)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil, degK

  !-------local var
  type(cohort_type),pointer :: cc
  real    :: rho_N_up0 = 0.02 ! hourly N uptake rate, fraction of the total mineral N
  real    :: N_roots0  = 0.1  ! root biomass at half max N-uptake rate,kg C m-2
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
!!       A scheme for deciduous to get enough N:
        cc%NSNmax = 0.2 * cc%crownarea  ! 5*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)) !
        NSN_not_full = (cc%NSN < cc%NSNmax) !
        if(NSN_not_full) N_Roots = N_Roots + cc%br * cc%nindivs

        end associate
     enddo
     ! M-M equation for Nitrogen absoption, McMurtrie et al. 2012, Ecology & Evolution
     ! rate at given root biomass and period of time
     if(N_roots>0.0)then
        ! Add a temperature response equation herefor rho_N_up0 (Zhu Qing 2016)
        rho_N_up = 1.-exp(-rho_N_up0 * N_roots/(N_roots0+N_roots) * hours_per_year * dt_fast_yr) ! rate at given root density and time period
        totNup = rho_N_up * vegn%mineralN  &
                * exp(9000.0 * (1./298.16 - 1./tsoil)) ! kgN m-2 time step-1
        avgNup = totNup / N_roots ! kgN time step-1 kg roots-1
        ! Nitrogen uptaken by each cohort, N_uptake
        vegn%N_uptake = 0.0
        do i = 1, vegn%n_cohorts
           cc => vegn%cohorts(i)
           associate ( sp => spdata(cc%species) )
           NSN_not_full = (cc%NSN < cc%NSNmax)
           if(NSN_not_full)then
               cc%N_uptake = cc%br  * avgNup  !
               cc%nsn      = cc%nsn + cc%N_uptake
               cc%annualNup  = cc%annualNup + cc%N_uptake/cc%crownarea

               ! subtract N from mineral N
               vegn%mineralN = vegn%mineralN - cc%N_uptake * cc%nindivs
               vegn%N_uptake = vegn%N_uptake + cc%N_uptake * cc%nindivs
               vegn%accu_Nup  = vegn%accu_Nup + cc%N_uptake * cc%nindivs

           endif
           end associate
        enddo
        cc =>null()
     endif ! N_roots>0

  endif
  
end subroutine vegn_N_uptake
! ============================================================================
! Nitrogen mineralization and immoblization with microbial C & N pools
! it's a new decomposition model with coupled C & N pools and variable 
! carbon use efficiency
subroutine SOMdecomposition(vegn, tsoil, thetaS)
  type(vegn_tile_type), intent(inout) :: vegn
  real                , intent(in)    :: tsoil ! soil temperature, deg K 
  real                , intent(in)    :: thetaS

  real :: CUE0=0.4  ! default microbial CUE
  real :: phoMicrobial = 2.5 ! turnover rate of microbes (yr-1)
  real :: CUEfast,CUEslow
  real :: CNm = 10.0  ! Microbial C/N ratio
  real :: NforM, fNM=0.0  ! mineral N available for microbes
  real :: micr_C_loss, fast_L_loss, slow_L_loss
  real :: runoff ! kg m-2 /step
  real :: N_loss
  real :: DON_fast,DON_slow,DON_loss ! Dissolved organic N loss, kg N m-2 step-1
  real :: fDON=0.0   ! 0.02     ! fractio of DON production in decomposition
  real :: fast_N_free 
  real :: slow_N_free 
  real :: CNfast, CNslow
  real :: A  ! decomp rate reduction due to moisture and temperature
  
!  runoff = vegn%Wrunoff * 365*24*3600 *dt_fast_yr !kgH2O m-2 s-1 ->kg m-2/time step
  runoff = vegn%runoff  !* dt_fast_yr !kgH2O m-2 yr-1 ->kgH2O m-2/time step, weng 2017-10-15
! CN ratios of soil C pools
  CNfast = vegn%metabolicL/vegn%metabolicN
  CNslow = vegn%structuralL/vegn%structuralN

!! C decomposition
!  A=A_function(tsoil,thetaS)
!  micr_C_loss = vegn%microbialC *A*phoMicrobial* dt_fast_yr
!  fast_L_loss = vegn%metabolicL*A*K1           * dt_fast_yr
!  slow_L_loss = vegn%structuralL*A*K2          * dt_fast_yr

! C decomposition
  A=A_function(tsoil,thetaS)
  micr_C_loss = vegn%microbialC * (1.0 - exp(-A*phoMicrobial* dt_fast_yr))
  fast_L_loss = vegn%metabolicL * (1.0 - exp(-A*K1          * dt_fast_yr))
  slow_L_loss = vegn%structuralL* (1.0 - exp(-A*K2          * dt_fast_yr))

! Carbon use efficiencies of microbes
  NforM = fNM * vegn%mineralN
  CUEfast = MIN(CUE0,CNm*(fast_L_loss/CNfast + NforM)/fast_L_loss)
  CUEslow = MIN(CUE0,CNm*(slow_L_loss/CNslow + NforM)/slow_L_loss)

! update C and N pools
! Carbon pools
  vegn%microbialC  = vegn%microbialC - micr_C_loss &
                    + fast_L_loss * CUEfast &
                    + slow_L_loss * CUEslow
  vegn%metabolicL = vegn%metabolicL - fast_L_loss
  vegn%structuralL = vegn%structuralL - slow_L_loss

! Find papers about soil DON losses
! DON loss, revised by Weng. 2016-03-03  ??
  fDON        = 0.25 ! 0.25 ! * dt_fast_yr ! 0.05 !* dt_fast_yr
  runoff      = 0.2 ! 0.2 ! mm day-1
  ! Assume it is proportional to decomposition rates
  ! Find some papers!!
  DON_fast    = fDON * fast_L_loss/CNfast * (etaN*runoff)
  DON_slow    = fDON * slow_L_loss/CNslow * (etaN*runoff)
  DON_loss    = DON_fast + DON_slow

! Update Nitrogen pools
  vegn%microbialN= vegn%microbialC/CNm
  vegn%metabolicN  = vegn%metabolicN  - fast_L_loss/CNfast - DON_fast
  vegn%structuralN = vegn%structuralN - slow_L_loss/CNslow - DON_slow

! Mixing of microbes to litters
  vegn%metabolicL   = vegn%metabolicL + MLmixRatio*fast_L_loss * CUEfast
  vegn%metabolicN   = vegn%metabolicN + MLmixRatio*fast_L_loss * CUEfast/CNm

  vegn%structuralL = vegn%structuralL + MLmixRatio*slow_L_loss * CUEslow
  vegn%structuralN = vegn%structuralN + MLmixRatio*slow_L_loss * CUEslow/CNm

  vegn%microbialC  = vegn%microbialC  - MLmixRatio*(fast_L_loss*CUEfast+slow_L_loss*CUEslow)
  vegn%microbialN  = vegn%microbialC/CNm
    
! update mineral N pool (mineralN)
  fast_N_free = MAX(0.0, fast_L_loss*(1./CNfast - CUEfast/CNm))
  slow_N_free = MAX(0.0, slow_L_loss*(1./CNslow - CUEslow/CNm))

! N_loss = MAX(0.,vegn%mineralN)        * A * K_nitrogen * dt_fast_yr
  N_loss = MAX(0.,vegn%mineralN) * (1.0-exp(-etaN*runoff - A*K_nitrogen * dt_fast_yr))
  vegn%Nloss_yr = vegn%Nloss_yr + N_loss + DON_loss

  vegn%mineralN = vegn%mineralN - N_loss       &
                  + vegn%N_input * dt_fast_yr  &
                  + fast_N_free + slow_N_free  &
                  + micr_C_loss/CNm
  vegn%annualN   = vegn%annualN - N_loss       &
                  + vegn%N_input * dt_fast_yr  &
                  + fast_N_free + slow_N_free  &
                  + micr_C_loss/CNm

! Check if soil C/N is above CN0
  fast_N_free = MAX(0., vegn%metabolicN  - vegn%metabolicL/CN0metabolicL)
  slow_N_free = MAX(0., vegn%structuralN - vegn%structuralL/CN0structuralL)
  vegn%metabolicN  = vegn%metabolicN  - fast_N_free
  vegn%structuralN = vegn%structuralN - slow_N_free
  vegn%mineralN    = vegn%mineralN + fast_N_free + slow_N_free
  vegn%annualN     = vegn%annualN  + fast_N_free + slow_N_free

! Heterotrophic respiration: decomposition of litters and SOM, kgC m-2 hour-1
  vegn%rh =  (micr_C_loss + fast_L_loss*(1.-CUEfast)+ slow_L_loss*(1.-CUEslow))

end subroutine SOMdecomposition

! ============================================================================
! The combined reduction in decomposition rate as a funciton of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
function A_function(tsoil, thetaS) result(A)
  real :: A                 ! return value, resulting reduction in decomposition rate
  real, intent(in) :: tsoil ! effective temperature for soil carbon decomposition
  real, intent(in) :: thetaS

  real :: soil_temp ! temperature of the soil, deg C
  real :: Td        ! rate multiplier due to temp
  real :: Wd        ! rate reduction due to mositure

  ! coefficeints and terms used in temperaturex term
  real :: Topt,Tmax,t1,t2,tshl,tshr

  soil_temp = tsoil-273.16

  ! EFFECT OF TEMPERATURE , ! from Bolker's century code
  Tmax=45.0 
  if (soil_temp > Tmax) soil_temp = Tmax 
  Topt=35.0 
  tshr=0.2 
  tshl=2.63 
  t1=(Tmax-soil_temp)/(Tmax-Topt) 
  t2=exp((tshr/tshl)*(1.-t1**tshl)) 
  Td=t1**tshr*t2 

  if (soil_temp > -10) Td=Td+0.05 
  if (Td > 1.) Td=1. 

  ! EFFECT OF MOISTURE
  ! Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272
  ! This differs from the Century Wd
  ! was modified by slm/ens based on the figures from the above paper 
  !     (not the reported function)

  if(thetaS <= 0.3) then
     Wd = 0.2 
  else if(thetaS <= 0.6) then
     Wd = 0.2+0.8*(thetaS-0.3)/0.3
  else 
     Wd = 1.0 ! exp(2.3*(0.6-thetaS)); ! Weng, 2016-11-26
  endif

  A = (Td*Wd)  ! the combined (multiplicative) effect of temp and water
               ! on decomposition rates
end function A_function

!=======================================================================


! =================== Cohort management ================================
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
  integer :: i,j,k

  allocate(cc(vegn%n_cohorts))
  merged(:)=.FALSE. ; k = 0
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
  real :: loss_fine,loss_coarse
  real :: lossN_fine,lossN_coarse
  integer :: i,j,k


 ! calculate the number of cohorts with indivs>mindensity
  k = 0
  do i = 1, vegn%n_cohorts
     if (vegn%cohorts(i)%nindivs > mindensity) k=k+1
  enddo
  if (k==0) write(*,*)'kill_lowdensity_cohorts: ','All cohorts have died'

  ! exclude cohorts that have low individuals
  if (k < vegn%n_cohorts) then
     allocate(cc(k))
     k=0
     do i = 1,vegn%n_cohorts
        cp =>vegn%cohorts(i)
        associate(sp=>spdata(cp%species))
        if (cp%nindivs > mindensity) then
           k=k+1
           cc(k) = cp
        else
           ! Carbon and Nitrogen from plants to soil pools
           loss_coarse  = cp%nindivs * (cp%bHW + cp%bsw   + cp%bl - cp%leafarea*LMAmin)
           loss_fine    = cp%nindivs * (cp%nsc + cp%seedC + cp%br + cp%leafarea*LMAmin)
           lossN_coarse = cp%nindivs * (cp%woodN + cp%sapwN + cp%leafN - sp%LNbase*cp%leafarea)
           lossN_fine   = cp%nindivs * (cp%rootN + cp%seedN + cp%NSN   + sp%LNbase*cp%leafarea)

           vegn%metabolicL  = vegn%metabolicL + fsc_fine *loss_fine + fsc_wood *loss_coarse
           vegn%structuralL = vegn%structuralL + (1.0-fsc_fine)*loss_fine + (1.0-fsc_wood)*loss_coarse

           vegn%metabolicN = vegn%metabolicN + &
                           fsc_fine *lossN_fine +    fsc_wood *lossN_coarse
           vegn%structuralN = vegn%structuralN + &
                       (1.-fsc_fine)*lossN_fine +(1.-fsc_wood)*lossN_coarse
           ! annual N from plants to soil
           vegn%N_P2S_yr = vegn%N_P2S_yr + lossN_fine + lossN_coarse
        endif
        end associate
     enddo
     vegn%n_cohorts = k
     deallocate (vegn%cohorts)
     vegn%cohorts=>cc
  endif
end subroutine kill_lowdensity_cohorts


! ============================================================================
subroutine merge_cohorts(c1,c2)
  type(cohort_type), intent(in) :: c1
  type(cohort_type), intent(inout) :: c2
  
  real :: x1, x2 ! normalized relative weights

  if(c1%nindivs > 0.0 .or. c2%nindivs > 0.0)then
     x1 = c1%nindivs/(c1%nindivs+c2%nindivs)
     x2 = c2%nindivs/(c1%nindivs+c2%nindivs)
  !else
  !   x1 = 0.5
  !   x2 = 0.5
  !endif
  ! update number of individuals in merged cohort
     c2%nindivs = c1%nindivs + c2%nindivs
  !  Carbon
     c2%bl  = x1*c1%bl  + x2*c2%bl
     c2%br  = x1*c1%br  + x2*c2%br
     c2%bsw = x1*c1%bsw + x2*c2%bsw
     c2%bHW = x1*c1%bHW + x2*c2%bHW
     c2%seedC = x1*c1%seedC + x2*c2%seedC
     c2%nsc = x1*c1%nsc + x2*c2%nsc
  !   Allometry
     c2%dbh = x1*c1%dbh + x2*c2%dbh
     c2%height = x1*c1%height + x2*c2%height
     c2%crownarea = x1*c1%crownarea + x2*c2%crownarea
     c2%age = x1*c1%age + x2*c2%age
     c2%carbon_gain = x1*c1%carbon_gain + x2*c2%carbon_gain
     c2%topyear = x1*c1%topyear + x2*c2%topyear

  !  Nitrogen
     c2%leafN = x1*c1%leafN + x2*c2%leafN
     c2%rootN = x1*c1%rootN + x2*c2%rootN
     c2%sapwN = x1*c1%sapwN + x2*c2%sapwN
     c2%woodN = x1*c1%woodN + x2*c2%woodN
     c2%seedN = x1*c1%seedN + x2*c2%seedN
     c2%NSN   = x1*c1%NSN   + x2*c2%NSN

  !  calculate the resulting dry heat capacity
     c2%leafarea = leaf_area_from_biomass(c2%bl, c2%species, c2%layer, c2%firstlayer)
  endif
end subroutine merge_cohorts

! ============================================================================
function cohorts_can_be_merged(c1,c2); logical cohorts_can_be_merged
   type(cohort_type), intent(in) :: c1,c2

   real, parameter :: mindensity = 1.0E-4
   logical :: sameSpecies, sameLayer, sameSize, sameSizeTree, sameSizeGrass, lowDensity

   sameSpecies  = c1%species == c2%species
   sameLayer    = (c1%layer == c2%layer) ! .and. (c1%firstlayer == c2%firstlayer)
   sameSizeTree = (spdata(c1%species)%lifeform > 0).and.  &
                  (spdata(c2%species)%lifeform > 0).and.  &
                 ((abs(c1%DBH - c2%DBH)/c2%DBH < 0.2 ) .or.  &
                  (abs(c1%DBH - c2%DBH)        < 0.001))  ! it'll be always true for grasses
   sameSizeGrass= (spdata(c1%species)%lifeform ==0) .and. &
                  (spdata(c2%species)%lifeform ==0) .and. &
                 (((c1%DBH == c2%DBH) .And.(c1%nsc == c2%nsc)) .OR. &
                    c1%topyear==c2%topyear)  ! it'll be always true for grasses
   sameSize = sameSizeTree .OR. sameSizeGrass
   lowDensity  = .FALSE. ! c1%nindivs < mindensity 
                         ! Weng, 2014-01-27, turned off
   cohorts_can_be_merged = sameSpecies .and. sameLayer .and. sameSize
end function

! ============================================================================
! calculate tree height, DBH, height, and crown area by initial biomass
! The allometry equations are from Ray Dybzinski et al. 2011 and Forrior et al. in review
!         HT = alphaHT * DBH ** (gamma-1)   ! DBH --> Height
!         CA = alphaCA * DBH ** gamma       ! DBH --> Crown Area
!         BM = alphaBM * DBH ** (gamma + 1) ! DBH --> tree biomass
subroutine initialize_cohort_from_biomass(cc,btot)
  type(cohort_type), intent(inout) :: cc
  real,intent(in)    :: btot ! total biomass per individual, kg C

  associate(sp=>spdata(cc%species))
     cc%DBH        = (btot / sp%alphaBM) ** ( 1.0/sp%thetaBM )
     cc%height     = sp%alphaHT * cc%dbh ** sp%thetaHT
     cc%crownarea  = sp%alphaCA * cc%dbh ** sp%thetaCA

     cc%bl_max = sp%LMA   * sp%LAImax        * cc%crownarea
     cc%br_max = sp%phiRL * sp%LAImax/sp%SRA * cc%crownarea
     cc%NSNmax = 0.2 * cc%crownarea ! 5.0*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
     cc%nsc    = 2.0 * (cc%bl_max + cc%br_max)
     call rootarea_and_verticalprofile(cc)
!    N pools
     cc%NSN    = 5.0*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
     cc%leafN  = cc%bl/sp%CNleaf0
     cc%rootN  = cc%br/sp%CNroot0
     cc%sapwN  = cc%bsw/sp%CNsw0
     cc%woodN  = cc%bHW/sp%CNwood0
  end associate
end subroutine initialize_cohort_from_biomass

! ============================================================================
subroutine annual_calls(vegn)
   type(vegn_tile_type), intent(inout) :: vegn

    ! ---------- annual call -------------
    ! update the LAImax of each PFT according to available N for next year
    ! call vegn_annualLAImax_update(vegn)

    ! Reproduction and mortality
    call vegn_reproduction(vegn)
    call vegn_nat_mortality(vegn, real(seconds_per_year))
    call vegn_starvation(vegn)  ! called daily

    ! Re-organize cohorts
    call relayer_cohorts(vegn)
    call kill_lowdensity_cohorts(vegn)
    call vegn_mergecohorts(vegn)

    ! set annual variables zero
    call Zero_diagnostics(vegn)
 end subroutine annual_calls

!============================================
subroutine init_cohort_allometry(cc)
  type(cohort_type), intent(inout) :: cc
  ! ----- local var -----------
  real    :: btot ! total biomass per individual, kg C

  btot = max(0.0001,cc%bHW+cc%bsw)
  associate(sp=>spdata(cc%species))
     cc%DBH        = (btot / sp%alphaBM) ** ( 1.0/sp%thetaBM )
!    cc%treeBM     = sp%alphaBM * cc%dbh ** sp%thetaBM
     cc%height     = sp%alphaHT * cc%dbh ** sp%thetaHT
     cc%crownarea  = sp%alphaCA * cc%dbh ** sp%thetaCA

     ! calculations of bl_max and br_max are here only for the sake of the
     ! diagnostics, because otherwise those fields are inherited from the 
     ! parent cohort and produce spike in the output, even though these spurious
     ! values are not used by the model
     cc%bl_max = sp%LMA   * sp%LAImax        * cc%crownarea
     cc%br_max = sp%phiRL * sp%LAImax/sp%SRA * cc%crownarea
     cc%NSNmax = 0.2 * cc%crownarea ! 5.0*(cc%bl_max/sp%CNleaf0 + cc%br_max/sp%CNroot0)
  end associate
end subroutine

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
  real   :: fixedN
  logical:: fixedN_based
  integer :: i
  ! Calculating LAI max based on mineral N or mineralN + fixed N
  fixedN_based =  .False. ! .True. !
  LAImin       = 0.5

  !fixedN = 0.0
  !do i = 1,vegn%n_cohorts
  !      cc => vegn%cohorts(i)
  !      fixedN = fixedN + cc%annualfixedN * cc%crownarea * cc%nindivs
  !enddo

 ! Mineral+fixed N-based LAImax
 ! LAI_fixedN = sp%Nfixrate0 * sp%LMA * sp%CNleaf0 * sp%leafLS / sp%LMA
 ! cc%br_max = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
  vegn%previousN = 0.8 * vegn%previousN + 0.2 * vegn%annualN
  do i=0,MSPECIES
      associate (sp => spdata(i) )

      LAIfixedN  = 0.5 * sp%Nfixrate0 * sp%CNleaf0 * sp%leafLS
      LAImineralN = 0.5*vegn%previousN*sp%CNleaf0*sp%leafLS/sp%LMA
      LAI_nitrogen = LAIfixedN + LAImineralN

      spdata(i)%LAImax = MAX(LAImin, MIN(LAI_nitrogen,sp%LAI_light))
      spdata(i)%underLAImax = MIN(sp%LAImax,1.2)
      end associate
  enddo

!  ! update the PFTs in the first layer based on fixed N
!  if(fixedN_based)then ! based on "cc%annualfixedN + vegn%previousN"
!!    Reset sp%LAImax
!     do i = 1,vegn%n_cohorts
!        cc => vegn%cohorts(i)
!        associate (sp => spdata(cc%species) )
!        sp%LAImax    = 0.0  ! max(sp%LAImax,ccLAImax)
!        sp%layerfrac = 0.0
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

! ============================================================================
function leaf_area_from_biomass(bl,species,layer,firstlayer) result (area)
  real :: area ! returned value
  real,    intent(in) :: bl      ! biomass of leaves, kg C/individual
  integer, intent(in) :: species ! species
  integer, intent(in) :: layer, firstlayer

! modified by Weng (2014-01-09), 07-18-2017
  area = bl/spdata(species)%LMA
  !if(layer > 1.AND. firstlayer == 0)then
  !   area = bl/(0.5*spdata(species)%LMA) ! half thickness for leaves in understory
  !else
  !   area = bl/spdata(species)%LMA
  !endif
end function
! ============================================================================

!============= Vegetation initializations =====================
subroutine initialize_vegn_tile(vegn,nCohorts,namelistfile)
   type(vegn_tile_type),intent(inout),pointer :: vegn
   integer,intent(in) :: nCohorts
   character(len=50),intent(in) :: namelistfile
!--------local vars -------

   type(cohort_type),dimension(:), pointer :: cc
   type(cohort_type),pointer :: cp
   integer,parameter :: rand_seed = 86456
   real    :: r
   real    :: btotal
   integer :: i, istat
   integer :: io           ! i/o status for the namelist
   integer :: ierr         ! error code, returned by i/o routines
   integer :: nml_unit

!  Read parameters from the parameter file (namelist)
   if(read_from_parameter_file)then
      ! --- Generate cohorts according to "initial_state_nml" ---
      nml_unit = 999
      open(nml_unit, file=namelistfile, form='formatted', action='read', status='old')
      read (nml_unit, nml=initial_state_nml, iostat=io, end=20)
20    close (nml_unit)
      write(*,nml=initial_state_nml)
      ! Initialize plant cohorts
      allocate(cc(1:init_n_cohorts), STAT = istat)
      vegn%cohorts => cc
      vegn%n_cohorts = init_n_cohorts
      cc => null()

      do i=1,init_n_cohorts
         cp => vegn%cohorts(i)
         cp%status  = LEAF_OFF ! ON=1, OFF=0 ! ON
         cp%layer   = 1
         cp%species = init_cohort_species(i)
         cp%ccID =  i
         cp%nsc     = init_cohort_nsc(i)
         cp%nindivs = init_cohort_nindivs(i) ! trees/m2
         cp%bsw     = init_cohort_bsw(i)
         cp%bHW   = init_cohort_bHW(i)
         btotal     = cp%bsw + cp%bHW  ! kgC /tree
         call initialize_cohort_from_biomass(cp,btotal)
      enddo
      MaxCohortID = cp%ccID
      ! Sorting these cohorts
      call relayer_cohorts(vegn)
      ! Initial Soil pools and environmental conditions
      vegn%metabolicL   = init_fast_soil_C ! kgC m-2
      vegn%structuralL  = init_slow_soil_C ! slow soil carbon pool, (kg C/m2)
      vegn%metabolicN   = vegn%metabolicL/CN0metabolicL  ! fast soil nitrogen pool, (kg N/m2)
      vegn%structuralN  = vegn%structuralL/CN0structuralL  ! slow soil nitrogen pool, (kg N/m2)
      vegn%N_input      = N_input  ! kgN m-2 yr-1, N input to soil
      vegn%mineralN     = init_Nmineral  ! Mineral nitrogen pool, (kg N/m2)
      vegn%previousN    = vegn%mineralN
      !Soil water
      vegn%FLDCAP = FLDCAP
      vegn%WILTPT = WILTPT
      vegn%wcl = FLDCAP
      ! Update soil water
      vegn%SoilWater = 0.0
      do i=1, max_lev
         vegn%SoilWater = vegn%SoilWater + vegn%wcl(i)*thksl(i)*1000.0
      enddo
      vegn%thetaS = 1.0

   else
     ! ------- Generate cohorts randomly --------
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
         cp%species = INT(rand()*5)+1
         cp%nindivs = rand()/10. ! trees/m2
         btotal     = rand()*100.0  ! kgC /tree
         call initialize_cohort_from_biomass(cp,btotal)
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
      vegn%metabolicL  = 0.2 ! kgC m-2
      vegn%structuralL = 7.0 ! slow soil carbon pool, (kg C/m2)
      vegn%metabolicN  = vegn%metabolicL/CN0metabolicL  ! fast soil nitrogen pool, (kg N/m2)
      vegn%structuralN = vegn%structuralL/CN0structuralL  ! slow soil nitrogen pool, (kg N/m2)
      vegn%N_input     = N_input  ! kgN m-2 yr-1, N input to soil
      vegn%mineralN    = 0.005  ! Mineral nitrogen pool, (kg N/m2)
      vegn%previousN   = vegn%mineralN

   endif  ! initialization: random or pre-described
end subroutine initialize_vegn_tile

! ====================================

end module esdvm



