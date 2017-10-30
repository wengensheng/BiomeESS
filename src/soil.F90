! The subroutines are from LM3PPA, the version used in Weng et al. 2016.
! This simulator can simulate evolutionarily stable strategy (ESS) of LMA
! and reproduce the forest succession patterns shown in Weng et al.,
! 2016 Global Change Biology along the graidient of temperature. 
!
! This simulator includes:
! photosynthesis, leaf stomatal conductance, transpiration,
! plant respiration, plant growth, reproduction, mortality,
! PPA, soil water dynamics, and soil organic matter decomposition.


module soil_mod
 use datatypes
 implicit none
 private

! ------ public subroutines ---------
public :: SoilWaterDynamicsBucket, water_supply_Bucket
public :: SoilWaterDynamicsLayer, water_supply_layer
public :: soil_data_beta

!---------------------------------

 contains
!=============== soil water subroutines ==================================
!=========================================================================

!========================================================================
! Weng 2017-10-18 ! compute available water for photosynthesis
subroutine water_supply_layer(forcing, vegn)
  type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

!----- local var --------------
  type(cohort_type),pointer :: cc
  real :: fWup(max_lev)      ! fraction to the actual soil water
  real :: freewater(max_lev)
  real :: W_up0(max_lev) ! potential water uptake, mol s-1 m-2
  real :: thetaS(max_lev) ! soil moisture index (0~1)
  real :: dpsiSR(max_lev) ! pressure difference between soil water and root water, MPa
  real  :: step_seconds ! seconds per step
  integer :: i,j, layer

  step_seconds = seconds_per_year * dt_fast_yr
!! Water pressure difference
  do i=1, max_lev ! Calculate water uptake potential layer by layer
     freewater(i) = max(0.0,((vegn%wcl(i)-WILTPT)*thksl(i)*1000.0))
     thetaS(i)    = max(0.0, (vegn%wcl(i)-WILTPT)/(FLDCAP-WILTPT))
     dpsiSR(i) = 1.5 *1.0e6 * thetaS(i)**2 ! Pa

     W_up0(i) = 0.0 ! Potential water uptake per layer
     do j = 1, vegn%n_cohorts
        cc => vegn%cohorts(j)
        associate ( sp => spdata(cc%species) )

        !cc%WupL(i) = cc%rootareaL(i)*sp%Kw_root*dpsiSR(i) * (step_seconds*mol_h2o) ! kg tree-1 step-1
        cc%WupL(i) = 0.00016*freewater(i)*cc%crownarea ! for test
        W_up0(i) = W_up0(i) + cc%WupL(i) * cc%nindivs ! water uptake per layer by all
        end associate
     enddo
     if(W_up0(i)>0.0)then
         fWup(i) = Min(0.2*freewater(i)/W_up0(i),1.0)! ratio of available soil water
     else
         fWup(i) = 0.0
     endif
     do j = 1, vegn%n_cohorts
        cc => vegn%cohorts(j)
        cc%WupL(i) = fWup(i) * cc%WupL(i) ! kg tree-1 step-1
     enddo
  enddo

! actual W_suplly
  do j = 1, vegn%n_cohorts
     cc => vegn%cohorts(j)
     cc%W_supply = 0.0
     do i=1,max_lev
        cc%W_supply = cc%W_supply + cc%WupL(i) ! kg tree-1 step-1
     enddo
  enddo

 end subroutine water_supply_layer

!========================================================================
! Weng 2017-10-18
! Calculate water supply by assuming soil is a bucket (no layers)
subroutine water_supply_bucket(forcing, vegn)
  type(climate_data_type),intent(in):: forcing
  type(vegn_tile_type), intent(inout) :: vegn

!----- local var --------------
  type(cohort_type),pointer :: cc
  real :: freewater, W_up0 ! potential water uptake by roots, kgH2O m-2 fasttimestep-1
  real :: fWup       ! fraction to the actual soil water
  real :: ws, thetaS ! soil moisture index (0~1)
  real :: dpsiSR     ! pressure difference between soil water and root water, MPa
  real :: step_seconds ! seconds per step
  integer :: i, layer

  step_seconds = seconds_per_year * dt_fast_yr

  ! Update soil water
   vegn%SoilWater = 0.0
   do i=1, max_lev
      vegn%SoilWater = vegn%SoilWater + vegn%wcl(i)*thksl(i)*1000.0
   enddo
!! Water supply for photosynthesis, a bucket
  ws = vegn%SoilWater/(rzone*1000.0)
  vegn%wcl = ws ! vol/vol
  thetaS = max(0.0, (ws - WILTPT)/(FLDCAP - WILTPT))
  vegn%thetaS = thetaS
  dpsiSR = 7.5 * SQRT(thetaS) *1.0e6

  W_up0 = 0.0
  freewater = max(0.0,(vegn%SoilWater-WILTPT*rzone*1000.0) )
  do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        associate ( sp => spdata(cc%species) )
        cc%W_supply = cc%rootarea * sp%Kw_root * dpsiSR ! mol s-1 tree-1
        W_up0 = W_up0 + cc%W_supply * cc%nindivs  * seconds_per_year * dt_fast_yr ! water uptake per unit ground area, kgH2O
        end associate

  enddo
  if(W_up0>0.0)then
     fWup = Min((0.5 * freewater)/W_up0,1.0)! rate at given root density and time period
  else
     fWup = 0.0
  endif
  ! Actual water uptake
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     cc%W_supply = fWup * cc%W_supply
     !write(*,*)cc%W_supply,cc%rootarea,fWup,freewater,W_up0,cc%nindivs
  enddo

end subroutine water_supply_bucket


! ============================================================================
! Weng, 2017-09-15, a bucket
subroutine SoilWaterDynamicsBucket(forcing,vegn)    !outputs
!     All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
  type(vegn_tile_type), intent(inout) :: vegn
  type(climate_data_type),intent(in):: forcing

!----- local var --------------
  type(cohort_type),pointer :: cc
  real    :: Edepth = 0.05     ! the depth of soil surface evaporation, m
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

! soil water refill by precipitation
   vegn%SoilWater = vegn%SoilWater + forcing%rain * seconds_per_year * dt_fast_yr
   vegn%runoff = max(0.0, vegn%SoilWater - FLDCAP*rzone*1000.0) ! mm hour-1
   vegn%SoilWater = vegn%SoilWater - vegn%runoff

! Water uptaken by roots, hourly
  vegn%transp = 0.0
  do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      vegn%transp = vegn%transp + cc%transp * cc%nindivs
  enddo
  vegn%SoilWater = vegn%SoilWater - vegn%transp

! Soil surface evaporation
!    thermodynamic parameters for air
      Rsoilabs = forcing%radiation * exp(-0.5*vegn%LAI)
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
      !Rsoil=3.0E+10 * (FILDCP-vegn%wcl(1))**16 ! Kondo et al. 1990
      rsoil=7500 * exp(-50.0*vegn%wcl(1))  ! s m-1
      raero=50./(forcing%windU + 0.2)
      rLAI=exp(vegn%LAI)
!     latent heat flux into air from soil
!           Eleaf(ileaf)=1.0*
!     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
!     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      Esoil=(slope*(Rsoilabs-Hgrownd)+rhocp*Dair/(raero+rLAI))/ &
            (slope+psyc*(rsoil/(raero+rLAI)+1.0))
!     sensible heat flux into air from soil
!      Hsoil = Rsoilabs - Esoil - Hgrownd

  !Calculate Esoil, kg m-2 s-1
  Esoil = 0.0
  vegn%evap = 0.0 ! min(Esoil/H2OLv * 3600.0, &
                  ! vegn%SoilWater * Edepth/rzone) ! kg m-2 hour-1
  ! Update soil water
  vegn%SoilWater = vegn%SoilWater - vegn%evap

end subroutine SoilWaterDynamicsBucket

! ============================================================================
! Weng, 2017-10-27
subroutine SoilWaterDynamicsLayer(forcing,vegn)    !outputs
!     All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
  type(vegn_tile_type), intent(inout) :: vegn
  type(climate_data_type),intent(in):: forcing

!----- local var --------------
  type(cohort_type),pointer :: cc
  real    :: rainwater,W_deficit(max_lev),W_add(max_lev)
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
  real    :: ftransp ! transp/W_supply
  integer :: i,j,k

! soil water refill by precipitation
   rainwater =  forcing%rain * seconds_per_year * dt_fast_yr
   do i=1, max_lev
      W_deficit(i) = (FLDCAP - vegn%wcl(i))*thksl(i)*1000.0
      W_add(i) = min(rainwater, W_deficit(i))
      vegn%wcl(i) = vegn%wcl(i) + W_add(i)/(thksl(i)*1000.0)
      rainwater = rainwater - W_add(i)
   enddo
   vegn%runoff = max(0.0, rainwater) ! mm hour-1

! Water uptaken by roots, hourly
  vegn%transp = 0.0
  do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      if(cc%W_supply>0.0)then
         ftransp = cc%transp/cc%W_supply
      else
         ftransp = 0.0
      endif
      do j=1,max_lev
         vegn%wcl(i) = vegn%wcl(i) - &
             (ftransp * cc%WupL(i) * cc%nindivs/(thksl(i)*1000))
      enddo
  enddo

! Soil surface evaporation
!    thermodynamic parameters for air
      Rsoilabs = forcing%radiation * exp(-0.65*vegn%LAI)
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
      !Rsoil=3.0E+10 * (FILDCP-vegn%wcl(1))**16 ! Kondo et al. 1990
      rsoil=7500 * exp(-50.0*vegn%wcl(1))  ! s m-1
      raero=50./(forcing%windU + 0.2)
      rLAI=exp(vegn%LAI)
!     latent heat flux into air from soil
!           Eleaf(ileaf)=1.0*
!     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
!     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      Esoil=(slope*(Rsoilabs-Hgrownd)+rhocp*Dair/(raero+rLAI))/ &
            (slope+psyc*(rsoil/(raero+rLAI)+1.0))
!     sensible heat flux into air from soil
!      Hsoil = Rsoilabs - Esoil - Hgrownd

  !Calculate Esoil, kg m-2 s-1
  ! Esoil = 0.0
  vegn%evap = min(Esoil/H2OLv * 3600.0, &
                  0.2*vegn%wcl(1) * thksl(1) *1000.) ! kg m-2 step-1
  vegn%wcl(1) = vegn%wcl(1) - vegn%evap/(thksl(1) *1000.)

  ! Update soil water
  vegn%SoilWater = 0.0
   do i=1, max_lev
      vegn%SoilWater = vegn%SoilWater + vegn%wcl(i)*thksl(i)*1000.0
   enddo

end subroutine SoilWaterDynamicsLayer

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

  ! ---- local vars
  integer :: k, l
  real    :: dz(max_lev)    ! thicknesses of layers
  real, dimension(num_l) :: &
       uptake_frac_max, & ! normalized root distribution
       vegn_uptake_term, &
       vlc, vsc, & ! volumetric fractions of water and ice in the layer
       root_length, & ! vertical distribution of volumetric root length, m/m3
       VRL, & ! volumetric root length
       u, du ! uptake and its derivative (the latter is not used)
  real :: z  !  soil depth
  real :: psi_wilt ! added by Weng, 2017-10-29
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

end subroutine

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
! ================================================

end module soil_mod



