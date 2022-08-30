! The subroutines are from LM3PPA, the version used in Weng et al. 2016.
! This simulator can simulate evolutionarily stable strategy (ESS) of LMA
! and reproduce the forest succession patterns shown in Weng et al.,
! 2016 Global Change Biology along the graidient of temperature.

module soil_mod
 use datatypes
 implicit none
 private

! ------ public subroutines ---------
public :: Soil_BGC, SoilWater_psi_K, SoilWaterDynamics, &
          SoilWaterTranspUpdate, SoilWaterSupply
public :: soil_data_beta

!---------------------------------
! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

! soil layer depth
real     :: dz(soil_L) = thksl   ! thicknesses of layers
real     :: zfull(soil_L)
real     :: zhalf(soil_L+1)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! =========================================================================
! =============== Soil Biogeochemical Cycle ===============================
! =========================================================================
! Nitrogen mineralization and immoblization with microbial C & N pools
! it's a new decomposition model with coupled C & N pools and variable
! carbon use efficiency
subroutine Soil_BGC (vegn, tsoil, thetaS)
  type(vegn_tile_type), intent(inout) :: vegn
  real                , intent(in)    :: tsoil ! soil temperature, deg K
  real                , intent(in)    :: thetaS

  !---- local var -------------
  real :: d_C(5), d_N(5), newM(5)
  real :: CUEf0, CUEs0
  real :: extraN, N_m    ! Mineralized nitrogen
  real :: N_loss ! Mineral Nitrogen loss, kg N m-2 step-1
  real :: runoff ! kg m-2 /step
  real :: DON_fast,DON_slow,DON_loss ! Dissolved organic N loss, kg N m-2 step-1
  real :: A  ! decomp rate reduction due to moisture and temperature
  integer :: i

  ! Default microbial CUE for fast and slow SOM
  CUEf0 = CUEmax0       ! 0.4
  CUEs0 = CUEmax0 * 0.5 ! 0.2

  ! Environmental scalar
  A=A_function(tsoil,thetaS)
  runoff = vegn%runoff  !mm/step, weng 2017-10-15

  ! Put litters into soil to start decomposition processes
  do i=1, 2
    d_C(i) = vegn%SOC(i) * K0SOM(i) * dt_fast_yr
    d_N(i) = vegn%SON(i) * K0SOM(i) * dt_fast_yr
    vegn%SOC(i)  = vegn%SOC(i) - d_C(i)
    vegn%SON(i)  = vegn%SON(i) - d_N(i)
    vegn%SOC(3+i) = vegn%SOC(3+i) + d_C(i)
    vegn%SON(3+i) = vegn%SON(3+i) + d_N(i)   
  enddo

  ! Turnover in SOM4 and SOM5
  do i=3, 5
     !d_C(i) = vegn%SOC(i)*(1. - exp(-A*K0SOM(i)*dt_fast_yr))
     d_C(i) = vegn%SOC(i) * A * K0SOM(i) * dt_fast_yr
     d_N(i) = vegn%SON(i) * A * K0SOM(i) * dt_fast_yr
  enddo

  ! New microbes grown from SOM decomposition
  newM(4) = Min(CUEf0*d_C(4), d_N(4)/CN0SOM(3))
  newM(5) = Min(CUEs0*d_C(5), d_N(5)/CN0SOM(3))
  newM(3) = (newM(4)+newM(5)) * (1.-f_M2SOM)

  ! Update C and N pools
  vegn%SOC(3) = vegn%SOC(3) - d_C(3) + newM(3)
  vegn%SOC(4) = vegn%SOC(4) - d_C(4) + newM(4) * f_M2SOM
  vegn%SOC(5) = vegn%SOC(5) - d_C(5) + newM(5) * f_M2SOM

  vegn%SON(3) = vegn%SON(3) - d_N(3) + newM(3) / CN0SOM(3)
  vegn%SON(4) = vegn%SON(4) - d_N(4) + newM(4) * f_M2SOM/CN0SOM(3)
  vegn%SON(5) = vegn%SON(5) - d_N(5) + newM(5) * f_M2SOM/CN0SOM(3)

  ! Mineralized nitrogen and Heterotrophic respiration, kg m-2 step-1
  vegn%rh = d_C(3) + d_C(4) + d_C(5) - newM(4) - newM(5) !
  N_m = d_N(3)+d_N(4)+d_N(5) - (newM(4)+newM(5))/CN0SOM(3)

  ! Organic and mineral nitrogen losses: Assume it is proportional to decomposition rates
  ! Find papers about these processes!! Experimental!
  DON_fast = fDON * d_N(4) * (etaN*runoff) + vegn%SON(4) * rho_SON * A * dt_fast_yr
  DON_slow = fDON * d_N(5) * (etaN*runoff) + vegn%SON(5) * rho_SON * A * dt_fast_yr
  DON_loss = DON_fast + DON_slow
  vegn%SON(4) = vegn%SON(4) - DON_fast
  vegn%SON(5) = vegn%SON(5) - DON_slow

  ! Mineral nitrogen loss
  !N_loss = MAX(0.,vegn%mineralN) * A * K_nitrogen * dt_fast_yr
  !N_loss = MAX(0.,vegn%mineralN) * (1. - exp(0.0 - etaN*runoff - A*K_nitrogen*dt_fast_yr))
  N_loss = vegn%mineralN * MIN(0.25, (A * K_nitrogen * dt_fast_yr + etaN*runoff))
  vegn%Nloss_yr = vegn%Nloss_yr + N_loss + DON_loss

  ! Update mineral N pool (mineralN)
  vegn%mineralN = vegn%mineralN + vegn%N_input * dt_fast_yr    &
                + N_m - N_loss
  vegn%annualN  = vegn%annualN  + vegn%N_input * dt_fast_yr    &
                + N_m - N_loss

 ! Check if soil C/N is lower than CN0
  do i=4, 5
     extraN = vegn%SON(i) - vegn%SOC(i)/CN0SOM(i)
     if (extraN > 0.0)then
        vegn%SON(i)   = vegn%SON(i)   - extraN
        vegn%mineralN = vegn%mineralN + extraN
        vegn%annualN  = vegn%annualN  + extraN
     endif
  enddo

end subroutine Soil_BGC

! =========================================================================
! =============== soil water subroutines ==================================

! ============================================================================
! Weng, 2017-10-27
subroutine SoilWaterDynamics(forcing,vegn)    !outputs
  ! All of inputs, the unit of water is 'mm',
  ! Soil moisture (soil water content) is a ratio
  type(vegn_tile_type), intent(inout) :: vegn
  type(climate_data_type),intent(in):: forcing

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real    :: WaterBudgetL(soil_L)
  real    :: rainwater,W_deficit(soil_L),W_add(soil_L)
  real    :: kappa  ! light extinction coefficient of corwn layers
  real    :: Esoil      ! soil surface evaporation, kg m-2 s-1
  real    :: Hsoil      ! sensible heat from soil
  real    :: Rsoilabs   ! W/m2
  real    :: Hgrownd    ! Ground heat flux, W/m2
  real    :: TairK,Tair      ! temperature, K and C, respectively
  real    :: RH         ! relative humidity, ratio to the saturated (0~1)
  real    :: Dair       ! VPD, pa
  real    :: rhocp !
  real    :: H2OLv
  real    :: slope
  real    :: psyc
  real    :: Cmolar ! mole density of air (mol/m3)
  real    :: rsoil  ! s m-1
  real    :: raero
  real    :: rLAI
  integer :: i,j,k

  WaterBudgetL = 0.0
  !! Soil surface evaporation
  !    calculate kappa  ! light extinction coefficient of corwn layers
  kappa = 0.75
  ! thermodynamic parameters for air
  Rsoilabs = forcing%radiation * exp(-kappa*vegn%LAI)
  Hgrownd = 0.0
  TairK = forcing%Tair
  Tair  = forcing%Tair - 273.16
  rhocp = cpair * forcing%P_air * mol_air / (Rugas*TairK)
  H2OLv =H2oLv0 - 2.365e3*Tair
  RH = forcing%RH  ! Check forcing's unit of humidity
  Dair  = esat(Tair)*(1.0 - RH)
  slope = (esat(Tair+0.1)-esat(Tair))/0.1
  psyc=forcing%P_air*cpair*mol_air/(H2OLv*mol_h2o)
  Cmolar=forcing%P_air/(Rugas*TairK) ! mole density of air (mol/m3)
  rsoil = exp(8.206-4.255*vegn%fldcap) ! s m-1, Liu Yanlan et al. 2017, PNAS
  !Rsoil=3.0E+10 * (FILDCP-vegn%wcl(1))**16 ! Kondo et al. 1990
  !rsoil=7500 * exp(-50.0*vegn%wcl(1))  ! s m-1
  raero=50./(forcing%windU + 0.2)
  rLAI=exp(vegn%LAI)

  !latent heat flux into air from soil
  !Eleaf(ileaf)=1.0*  &
  !     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/  &  !2* Weng 0215
  !     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))

  Esoil=(slope*Rsoilabs + rhocp*Dair/raero)/ &
        (slope + psyc*(1.0+rsoil/raero)) *   &
        max(vegn%wcl(1),0.0)/vegn%FLDCAP ! (vegn%wcl(1)-ws0)/(vegn%FLDCAP-ws0)
  !sensible heat flux into air from soil
  !Hsoil = Rsoilabs - Esoil - Hgrownd

  !Calculate Esoil, kg m-2 step-1
  vegn%evap = min(Esoil/H2OLv * step_seconds, &
                  0.2*vegn%wcl(1) * thksl(1) *1000.) ! kg m-2 step-1
  !vegn%wcl(1) = vegn%wcl(1) - vegn%evap/(thksl(1) *1000.)
  WaterBudgetL(1) = WaterBudgetL(1) - vegn%evap

  !! soil water refill by precipitation
  rainwater =  forcing%rain * step_seconds
  if(rainwater > 0.0)then
     do i=1, soil_L
        W_deficit(i) = (vegn%FLDCAP - vegn%wcl(i)) * thksl(i)*1000.0
        W_add(i) = min(rainwater, W_deficit(i))
        rainwater = rainwater - W_add(i)
        !vegn%wcl(i) = vegn%wcl(i) + W_add(i)/(thksl(i)*1000.0)
        WaterBudgetL(i) = WaterBudgetL(i) + W_add(i)

        if(rainwater<=0.0)exit
     enddo
  endif
  vegn%runoff = rainwater ! mm step-1

  ! Total soil water
  do i=1,soil_L
     vegn%wcl(i) = vegn%wcl(i) +  WaterBudgetL(i)/(thksl(i)*1000.0)
     vegn%freewater(i) = max(0.0,((vegn%wcl(i)-vegn%WILTPT)*thksl(i)*1000.0)) ! kg/m2, or mm
  enddo
  vegn%soilwater = sum(vegn%freewater(:))

end subroutine SoilWaterDynamics

!======================================================================
subroutine SoilWater_psi_K(vegn)
  !@sum leaf and stem equilibrium water potential
  !@+   Weng, 03/20/2022
  implicit none
  type(vegn_tile_type), intent(inout) :: vegn

  !------local var -----------
  real :: psi_sat_ref,k_sat_ref,v_sat,chb
  integer :: i

  !! Soil Water psi and K
  psi_sat_ref = soilpars(vegn%soiltype)%psi_sat_ref
  k_sat_ref   = soilpars(vegn%soiltype)%k_sat_ref
  v_sat       = soilpars(vegn%soiltype)%vwc_sat
  chb         = soilpars(vegn%soiltype)%chb
  do i=1, soil_L
     vegn%psi_soil(i) = calc_soil_psi(psi_sat_ref,chb,v_sat,vegn%wcl(i))
     vegn%K_soil(i)   = calc_soil_K  (k_sat_ref,  chb,v_sat,vegn%wcl(i))
  enddo
end subroutine SoilWater_psi_K

!==========================================
! Weng, 03/21/2022
real function calc_soil_psi(psi_sat_ref,chb,V_sat,Vwc) result(psi)
    !Calculate soil water potential (psi, MPa)
    real,intent(in):: psi_sat_ref
    real,intent(in):: chb          ! Parameter for scaling
    real,intent(in):: V_sat        ! Field capacity
    real,intent(in):: Vwc          ! Soil water content
    !---------------------
    psi = psi_sat_ref*1.0e-6 * ((V_sat/Vwc)**chb)    ! MPa
end function calc_soil_psi

!==========================================
! Weng, 03/21/2022
real function calc_soil_K(k_sat_ref,chb,V_sat,Vwc) result(K)
    !Calculate soil water conductivity (K, kg H2O/(m2 MPa s))
    real,intent(in):: k_sat_ref
    real,intent(in):: chb          ! Parameter for scaling
    real,intent(in):: V_sat        ! Field capacity
    real,intent(in):: Vwc          ! Soil water content
    !---------------------
    k = 18./1000. * k_sat_ref * (Vwc/V_sat)**(2*chb+3)! kg H2O/(m2 MPa s)

end function calc_soil_K

!========================================================================
! Weng 2017-10-18 ! compute available water for photosynthesis
!========================================================================
! Calculate soil water uptake by plants and supply for transpiration
subroutine SoilWaterSupply(vegn) ! forcing,
  !type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real :: freewater(soil_L)
  real :: thetaS(soil_L) ! soil moisture index (0~1)
  real :: LayerTot(soil_L) ! potential water uptake, kg H2O s-1 m-2
  real :: dpsiSR(soil_L) ! pressure difference between soil and root, MPa
  real :: fWup(soil_L)      ! fraction to the actual soil water
  integer :: i,j

  ! Calculating soil water availability for transpiration in next step
  do i=1, soil_L ! Calculate water uptake potential layer by layer
     freewater(i) = max(0.0,((vegn%wcl(i)-vegn%WILTPT) * thksl(i) * 1000.0)) ! kg/m2, or mm
     thetaS(i)    = max(0.0, (vegn%wcl(i)-vegn%WILTPT)/(vegn%FLDCAP-vegn%WILTPT))
     ! The difference of water potential between roots and soil
     dpsiSR(i) = 1.5 * thetaS(i)**2 ! *1.0e6  MPa

     ! Water uptake capacity
     LayerTot(i) = 0.0 ! Potential water uptake per layer by all cohorts
     do j = 1, vegn%n_cohorts
        cc => vegn%cohorts(j)
        ! Potential water uptake per soil layer by all cohorts
        cc%WupL(i) = cc%rootareaL(i) * vegn%K_soil(i) * dpsiSR(i) * step_seconds
        LayerTot(i) = LayerTot(i) + cc%WupL(i) * cc%nindivs
     enddo

     ! Adjust cc%WupL(i) according to soil available water
     do j = 1, vegn%n_cohorts
        cc => vegn%cohorts(j)
        if(LayerTot(i)>0.0) &
            fWup(i) = Min(0.2 * freewater(i) / LayerTot(i),1.0)! ratio of available soil water
        cc%WupL(i) = fWup(i) * cc%WupL(i) ! kg tree-1 step-1
     enddo ! cohort for each layer
  enddo    ! all layers

  ! total water suplly for next step's transpiration
  do j = 1, vegn%n_cohorts
     cc => vegn%cohorts(j)
     cc%W_supply = sum(cc%WupL(:))
  enddo

end subroutine SoilWaterSupply

! ============================================================================
! 07/07/2022, Weng
! Update soil water by deducting transpiration
subroutine SoilWaterTranspUpdate(vegn)
  type(vegn_tile_type),intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  real    :: WaterBudgetL(soil_L)
  real    :: W_supply, fsupply ! fraction of transpiration from a soil layer
  integer :: i,j

  ! Soil water conditions
  ! call SoilWaterSupply(vegn)

  ! Water uptaken by roots in each soil layer
  WaterBudgetL = 0.0
  do j = 1, vegn%n_cohorts
      cc => vegn%cohorts(j)
      W_supply = sum(cc%WupL)
      if(W_supply > 0.0)then
         do i=1,soil_L
            fsupply = cc%WupL(i)/W_supply
            WaterBudgetL(i) = WaterBudgetL(i) - fsupply * cc%transp * cc%nindivs
         enddo
      endif
  enddo ! all cohorts
  ! Deduct soil water
  do i=1,soil_L
     vegn%wcl(i) = vegn%wcl(i) + WaterBudgetL(i)/(thksl(i)*1000.0)
  enddo
end subroutine SoilWaterTranspUpdate


! ==============Used in LM3, but not here =============================
! =============== just for reference =========================================
! ============================================================================
! compute uptake-related properties
subroutine soil_data_beta(soil, vegn, soil_beta, soil_water_supply, &
                            soil_uptake_T )
  type(soil_tile_type), intent(in)    :: soil
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(out) :: soil_beta(:) ! relative water availability, used only in VEGN_PHOT_SIMPLE treatment
  real, intent(out) :: soil_water_supply(:) ! max rate of water supply to roots, kg/(indiv s)
  real, intent(out) :: soil_uptake_T(:) ! an estimate of temperature of the water
             ! taken up by transpiration. In case of 'linear' uptake it is an exact
             ! value; in case of 'darcy*' treatments the actual uptake profile
             ! is calculated only in step 2, so the value returned is an estimate

  ! ---- local vars -----------
  integer :: k, l
  real    :: dz(soil_L)    ! thicknesses of layers
  real, dimension(num_l) :: &
       uptake_frac_max, & ! normalized root distribution
       vegn_uptake_term, &
       vlc, vsc, & ! volumetric fractions of water and ice in the layer
       root_length, & ! vertical distribution of volumetric root length, m/m3
       VRL, & ! volumetric root length
       u, du ! uptake and its derivative (the latter is not used)
  real :: z  !  soil depth
  !real :: psi_wilt ! added by Weng, 2017-10-29
  logical :: uptake_oneway = .TRUE. ! added by Weng
  logical :: uptake_from_sat = .true.
  type (cohort_type), pointer :: cc

  dz = thksl ! Weng
  do l = 1, num_l
    vlc(l) = max(0., soil%prog(l)%wl / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%prog(l)%ws / (dens_h2o*dz(l)))
  enddo

  ! calculate volumetric root length for the entire tile
  VRL(:) = 0.0
  do k = 1, vegn%n_cohorts
     cc=>vegn%cohorts(k)
     call cohort_root_properties(cc, dz(1:num_l), cc%root_length(1:num_l), &
                                 cc%K_r, cc%r_r)
     VRL(:) = VRL(:)+cc%root_length(1:num_l)*cc%nindivs
  enddo

  ! calculate characteristic half-distance between roots, m
  where (VRL(:) > 0)
     vegn%root_distance(1:num_l) = 1.0/sqrt(PI*VRL(:))
  elsewhere
     vegn%root_distance(1:num_l) = 1.0 ! the value doesn't matter since uptake is 0 anyway
  end where

  do k = 1, vegn%n_cohorts
     cc=>vegn%cohorts(k)
     call cohort_uptake_profile (cc, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

     do l = 1, num_l
        cc%uptake_frac(l) = uptake_frac_max(l) &
             * max(0.0, min(1.0,(vlc(l)-soil%w_wilt(l))/&
                  (0.75*(soil%w_fc(l)-soil%w_wilt(l)))))
     enddo
     soil_beta(k) = sum(cc%uptake_frac(:))
     if (soil_beta(k) /= 0) then
          cc%uptake_frac(:) = cc%uptake_frac(:) / soil_beta(k)
     else
          cc%uptake_frac(:) = uptake_frac_max(:)
     endif

     ! calculate total water supply
     call darcy2d_uptake_lin ( soil, psi_wilt, vegn%root_distance, cc%root_length, &
             cc%K_r, cc%r_r, uptake_oneway, uptake_from_sat, u, du)
     soil_water_supply(k) = max(0.0,sum(u))
     soil_uptake_T(k) = soil%uptake_T

  enddo
end subroutine soil_data_beta
! ============================================================================
subroutine darcy2d_uptake_lin ( soil, psi_x0, R, VRL, K_r, r_r,uptake_oneway, &
    uptake_from_sat, u, du )
  type(soil_tile_type), intent(in) :: soil
  real, intent(in) :: &
       psi_x0,    & ! water potential inside roots (in xylem) at zero depth, m
       R(:),      & ! characteristic half-distance between roots, m
       VRL(:),    & ! Volumetric Root Length (root length per unit volume), m/m3
       K_r,       & ! permeability of the root membrane per unit area, kg/(m3 s)
       r_r          ! radius of fine roots, m
  logical, intent(in) :: &
       uptake_oneway, & ! if true, then the roots can only take up water, but
                   ! never loose it to the soil
       uptake_from_sat   ! if false, uptake from saturated soil is prohibited
  real, intent(out) :: &
       u(:), &      ! layer-by-layer distribution of uptake, kg/(m2 s)
       du(:)        ! derivative of u w.r.t. root water potential, kg/(m3 s)
  ! ---- local vars
  integer :: k
  real :: psi_x     ! water potential inside roots (psi_x0+z), m
  real :: psi_soil  ! water potential of soil, m
  real :: psi_sat   ! saturation soil water potential, m
  real :: K_sat     ! hydraulic conductivity of saturated soil, kg/(m2 s)

  real :: psi_root  ! water potential at the root/soil interface, m
  real :: psi_root0 ! initial guess of psi_root, m


  ! calculate some hydraulic properties common for all soil layers
  psi_sat = soil%pars%psi_sat_ref/soil%pars%alpha
  K_sat   = soil%pars%k_sat_ref*soil%pars%alpha**2

  u = 0; du = 0
  do k = 1, num_l
     psi_x    = psi_x0 + zfull(k)
     psi_soil = soil%psi(k)
     psi_root0= soil%psi(k) ! change it later to prev. time step value
     if ( soil%prog(k)%ws > 0 ) &
          cycle ! skip layers with ice
     if ( uptake_oneway.and.psi_x > soil%psi(k) ) &
          cycle ! skip layers where roots would loose water
     if ( .not.(uptake_from_sat).and.psi_soil >= psi_sat ) &
          cycle ! skip layers where the soil is saturated

     ! calculates soil term of uptake expression
     call darcy2d_flow_lin (psi_x, psi_soil, psi_root0, K_sat, psi_sat, soil%pars%chb, &
          K_r, r_r, R(k), u(k), du(k), psi_root)

     ! scale by volumetric root length and thickness of layer to get total
     ! uptake from the current soil layer
     u(k)  = VRL(k)*dz(k)*u(k)
     du(k) = VRL(k)*dz(k)*du(k)
  enddo

end subroutine darcy2d_uptake_lin

! ============================================================================
! given soil and root parameters, calculate the flux of water toward root
! per unit root length, and its derivative w.r.t. xylem water potential
! this version calculates fluxes linearized around psi_root0
subroutine darcy2d_flow_lin (psi_x, psi_soil, psi_root0, K_sat, psi_sat, b, K_r, &
     r_r, R, u, du, psi_root)
  real, intent(in) :: &
       psi_x,    & ! xylem water potential, m
       psi_soil, & ! soil water potential, m
       psi_root0,& ! value of psi_root we linearize around, m
       K_sat,    & ! saturated soil hydraulic conductivity, kg/(m2 s)
       psi_sat,  & ! saturates soil water potential, m
       b,        & ! power of soil moisture characteristic function
       K_r,      & ! root membrane permeability per unit area, kg/(m3 s)
       r_r,      & ! radius of root, m
       R           ! characteristic radial half-distance between roots, m
  real, intent(out) :: &
       u,        & ! uptake, kg/(m s)
       du,       & ! derivative of uptake w.r.t psi_x, kg/(m2 s)
       psi_root    ! water potential at the root/soil interface, m

  ! ---- local vars
  real :: u_soil0 ! flux through soil for psi_root = psi_root0
  real :: du_soil ! its derivative w.r.t. psi_root
  real :: C_r !
  real :: n
  real :: K_root  ! root membrane permeability per unit length, kg/(m2 s)

  C_r=2*PI/(log(R/r_r))
  n = -(1+3/b)
  K_root = 2*PI*r_r*K_r

  ! calculate flux through soil for psi_root = psi_root0
  u_soil0 = C_r*K_sat*&
       (psi_sat/n* &
          (  (min(psi_soil ,psi_sat)/psi_sat)**n   &
            -(min(psi_root0,psi_sat)/psi_sat)**n ) &
          + max(0.0, psi_soil  - psi_sat)          &
          - max(0.0, psi_root0 - psi_sat)          )
  ! and its derivative w.r.t. psi_root at psi_root0
  du_soil=-C_r*K_sat*(min(psi_root0,psi_sat)/psi_sat)**(n-1)

  ! flux through soil+membrane
  u  = K_root/(-du_soil+K_root)*(u_soil0+du_soil*(psi_x-psi_root0))
  ! and its derivative w.r.t. psi_x
  du = K_root/(-du_soil+K_root)*du_soil
  ! water potential at the root-soil interface
  psi_root = psi_x + u/K_root
end subroutine

! ============================================================================
! returns properties of the fine roots
subroutine cohort_root_properties(cohort, dz, vrl, K_r, r_r)
  type(cohort_type), intent(in)  :: cohort
  real, intent(in)  :: dz(:)
  real, intent(out) :: &
       vrl(:), & ! volumetric fine root length, m/m3
       K_r,    & ! root membrane permeability per unit area, kg/(m3 s)
       r_r       ! radius of fine roots, m

  integer :: sp, l
  real :: factor, z
  real :: vbr ! volumetric biomass of fine roots, kg C/m3

  sp = cohort%species

  factor = 1.0/(1.0-exp(-sum(dz)/cohort%root_zeta))
  z = 0
  do l = 1, size(dz)
     ! calculate the volumetric fine root biomass density [kgC/m3] for current layer
     ! NOTE: sum(brv*dz) must be equal to cohort%br, which is achieved by normalizing
     ! factor
     vbr = cohort%br * &
          (exp(-z/cohort%root_zeta) - exp(-(z+dz(l))/cohort%root_zeta))*factor/dz(l)
     ! calculate the volumetric fine root length
     vrl(l) = vbr*spdata(sp)%srl

     z = z + dz(l)
  enddo

  K_r = spdata(sp)%root_perm
  r_r = spdata(sp)%root_r

end subroutine

! ============================================================================
! calculates vertical distribution of active roots: given layer thicknesses,
! returns fraction of active roots per level
subroutine cohort_uptake_profile(cohort, dz, uptake_frac_max, vegn_uptake_term)
  type(cohort_type), intent(in)  :: cohort
  real, intent(in)  :: dz(:)
  real, intent(out) :: uptake_frac_max(:)
  real, intent(out) :: vegn_uptake_term(:)

  real, parameter :: res_scaler = mol_air/mol_h2o  ! scaling factor for water supply
  ! NOTE: there is an inconsistency there between the
  ! units of stomatal conductance [mol/(m2 s)], and the units of humidity deficit [kg/kg],
  ! in the calculations of water demand. Since the uptake options other than LINEAR can't
  ! use res_scaler, in this code the units of humidity deficit are converted to mol/mol,
  ! and the additional factor is introduced in res_scaler to ensure that the LINEAR uptake
  ! gives the same results.

  integer :: l
  real    :: z, sum_rf

     !linear scaling, LM3V
     z = 0
     do l = 1, size(dz)
        uptake_frac_max(l) = (exp(-z/cohort%root_zeta)    &
                - exp(-(z+dz(l))/cohort%root_zeta))
        uptake_frac_max(l) = &
                max( uptake_frac_max(l), 0.0)
        z = z + dz(l)
     enddo

  sum_rf = sum(uptake_frac_max)
  if(sum_rf>0) &
       uptake_frac_max(:) = uptake_frac_max(:)/sum_rf

  if (cohort%br <= 0) then
     vegn_uptake_term(:) = 0.0
  else
     vegn_uptake_term(:) = uptake_frac_max(:) * &
          res_scaler * spdata(cohort%species)%root_r * cohort%br
  endif

end subroutine
! ================================================

end module soil_mod
