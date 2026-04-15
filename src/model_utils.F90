module model_utils
  use datatypes

  implicit none
  private

  public :: model_para_init
  public :: read_init_namelist, read_vegn_namelist, read_soil_namelist
  public :: read_global_setting, Preset_GlobalPFTs
  public :: initialize_soilpars, initialize_PFT_pars, init_derived_species_data
  public :: Set_PFTs_from_map, Climate_envelope_vars
  public :: Set_PFTs_from_Data, Set_PFTs_from_Climate, Set_PFTs_from_LandCover
  public :: Assign_Std_Cohorts
  public :: vegn_sum_tile, Zero_diagnostics
  public :: qscomp, BM2Architecture
  public :: DBH2HT, DBH2CA, DBH2BM, BM2DBH
  public :: ccNSNmax, CA2BLmax, BLmax2BRmax, BL2Aleaf, Aleaf2LAI
  public :: A_function, calc_solarzen, PotentialET, c3c4
  public :: esat
  public :: rank_descending, merge, mergerank

contains

  !========================== Summarize tile variables =====================
  ! Weng, 2021-06-02
  subroutine vegn_sum_tile(vegn)
    implicit none
    type(vegn_tile_type), intent(inout) :: vegn

    type(cohort_type),pointer :: cc
    real :: BMG
    integer :: i

    vegn%NSC     = 0.0
    vegn%SeedC   = 0.0
    vegn%leafC   = 0.0
    vegn%rootC   = 0.0
    vegn%SwC     = 0.0
    vegn%HwC     = 0.0

    vegn%NSN     = 0.0
    vegn%SeedN   = 0.0
    vegn%leafN   = 0.0
    vegn%rootN   = 0.0
    vegn%SwN     = 0.0
    vegn%HwN     = 0.0

    vegn%W_stem = 0.0
    vegn%W_dead = 0.0
    vegn%W_leaf = 0.0

    vegn%LAI    = 0.0
    vegn%CAI    = 0.0
    vegn%ArootL = 0.0
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      associate ( sp => spdata(cc%species))
        vegn%NSC     = vegn%NSC     + cc%NSC    * cc%nindivs
        vegn%SeedC   = vegn%SeedC   + cc%seedC  * cc%nindivs
        vegn%leafC   = vegn%leafC   + cc%bl     * cc%nindivs
        vegn%rootC   = vegn%rootC   + cc%br     * cc%nindivs
        vegn%SwC     = vegn%SwC     + cc%bsw    * cc%nindivs
        vegn%HwC     = vegn%HwC     + cc%bHW    * cc%nindivs
        vegn%CAI     = vegn%CAI     + cc%Acrown * cc%nindivs
        vegn%LAI     = vegn%LAI     + cc%Aleaf  * cc%nindivs
        vegn%ArootL  = vegn%ArootL  + cc%ArootL * cc%nindivs

        vegn%NSN     = vegn%NSN     + cc%NSN    * cc%nindivs
        vegn%SeedN   = vegn%SeedN   + cc%seedN  * cc%nindivs
        vegn%leafN   = vegn%leafN   + cc%leafN  * cc%nindivs
        vegn%rootN   = vegn%rootN   + cc%rootN  * cc%nindivs
        vegn%SwN     = vegn%SwN     + cc%swN    * cc%nindivs
        vegn%HwN     = vegn%HwN     + cc%hwN    * cc%nindivs

        vegn%W_stem  = vegn%W_stem  + cc%W_stem * cc%nindivs
        vegn%W_dead  = vegn%W_dead  + cc%W_dead * cc%nindivs
        vegn%W_leaf  = vegn%W_leaf  + cc%W_leaf * cc%nindivs
      end associate
    enddo

    if(all(vegn%cohorts(:)%status == LEAF_ON))then
      vegn%TreeCA  = 0.0
      vegn%GrassCA = 0.0
      BMG          = 0.0
      do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        associate ( sp => spdata(cc%species))
          if(sp%lifeform==0) BMG = BMG + (cc%bl+cc%br+cc%bsw)*cc%nindivs
          if(cc%layer == 1)then
            if(sp%lifeform==0) then
              vegn%GrassCA = vegn%GrassCA + cc%Acrown*cc%nindivs
            else
              vegn%TreeCA  = vegn%TreeCA  + cc%Acrown*cc%nindivs
            endif
          endif
        end associate
      enddo
      vegn%GrassBM = max(vegn%GrassBM, BMG)
    endif

  end subroutine vegn_sum_tile

  !==================================================================================================
  ! Weng, 2016-11-28
  subroutine Zero_diagnostics(vegn)
    type(vegn_tile_type), intent(inout) :: vegn
    type(cohort_type),pointer :: cc
    integer :: i

    vegn%NfixDaily = 0.0
    vegn%dailyPrcp = 0.0
    vegn%dailyTrsp = 0.0
    vegn%dailyEvap = 0.0
    vegn%dailyRoff = 0.0
    vegn%dailyNup  = 0.0
    vegn%dailyGPP  = 0.0
    vegn%dailyNPP  = 0.0
    vegn%dailyResp = 0.0
    vegn%dailyRh   = 0.0
    vegn%dailyCH4  = 0.0
    vegn%dNorg_daily = 0.0
    vegn%dNgas_daily = 0.0
    vegn%dNmin_daily = 0.0

    vegn%NfixedYr   = 0.0
    vegn%annualPrcp = 0.0
    vegn%annualTrsp = 0.0
    vegn%annualEvap = 0.0
    vegn%annualRoff = 0.0
    vegn%annualGPP  = 0.0
    vegn%annualNPP  = 0.0
    vegn%annualResp = 0.0
    vegn%annualRh   = 0.0
    vegn%annualCH4  = 0.0
    vegn%NorgP2S    = 0.0
    vegn%Nm_Soil    = 0.0
    vegn%Nm_Fire    = 0.0
    vegn%C_burned   = 0.0
    vegn%NupYr      = 0.0
    vegn%dNorg_Yr   = 0.0
    vegn%dNgas_Yr   = 0.0
    vegn%dNmin_Yr   = 0.0
    vegn%GrassBM    = 0.0
    vegn%annualPET  = 0.0
    vegn%YearlyTmp  = 0.0

    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      cc%gpp       = 0.0
      cc%npp       = 0.0
      cc%resp      = 0.0
      cc%resl      = 0.0
      cc%resr      = 0.0
      cc%resg      = 0.0
      cc%transp    = 0.0
      cc%dailyWdmd = 0.0
      cc%dailyTrsp = 0.0
      cc%dailyGPP  = 0.0
      cc%dailyNPP  = 0.0
      cc%dailyResp = 0.0
      cc%dailyNup  = 0.0
      cc%NfixDaily = 0.0
      cc%annualTrsp= 0.0
      cc%annualGPP = 0.0
      cc%annualNPP = 0.0
      cc%annualResp= 0.0
      cc%NupYr     = 0.0
      cc%NfixedYr  = 0.0
      cc%totDemand = 0.0
      cc%NPPleaf   = 0.0
      cc%NPProot   = 0.0
      cc%NPPwood   = 0.0
      cc%DBH_ys    = cc%DBH
      cc%Aleafmax  = 0.0
    enddo
  end subroutine Zero_diagnostics

  !==========================================================================
  ! Initialize model processes/settings and model parameters
  !----------------------------------------------------------------
  subroutine model_para_init(fnml)
    character(len=*),intent(in) :: fnml
    call initialize_soilpars(fnml)
    call initialize_PFT_pars(fnml)

    ! Hack for closedN setting
    if(do_closedN_run) then
      K_DeNitr  = 0.0 ! rate of a year, 2.5
      rho_SON   = 0.0 ! organic nitrogen release rate
      fdsvN     = 0.0 ! Maximum nitrogen loss rate with runoff
      N_input   = 0.0 ! N input, kg N m-2 yr-1
    endif

  end subroutine model_para_init

  !============================ Subroutines =================================
  ! Read in parameters in the namelist file, Weng, 09/07/2025

  !----------------------------------------------------------------
  subroutine read_init_namelist(fnml)
    character(len=*),intent(in) :: fnml
    !--------local vars -----------
    integer :: rc, fu

    ! Check whether file exists
    inquire (file=fnml, iostat=rc)
    if (rc /= 0) then
      write (*, '("Error: input file ", a, " does not exist")') fnml
      stop
    end if

    ! Open and read Namelist file.
    open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
    read (nml=initial_state_nml, iostat=rc, unit=fu)
    if (rc == 0) then
      write(*,*)'Namelist initial_state_nml read successfully.'
    else
      write(*,*)'Namelist initial_state_nml error', rc
      stop
    endif
    !write(*,nml=initial_state_nml)
    close (fu)

  end subroutine read_init_namelist

  !----------------------------------------------------------------
  subroutine read_vegn_namelist(fnml)
    character(len=*),intent(in) :: fnml
    !--------local vars -----------
    integer :: rc, fu

    open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
    read (nml=vegn_parameters_nml, iostat=rc, unit=fu)
    if (rc == 0) then
      write(*,*)'Namelist vegn_parameters_nm read successfully.'
    else
      write(*,*)'Namelist vegn_parameters_nml error', rc
      stop
    endif
    !write(*,nml=vegn_parameters_nml)
    close (fu)

    ! Ajusting parameters if they are not in the resonable ranges
    CH4_wfps0 = max(0.0, min(CH4_wfps0, 0.8))
    CH4_wfps1 = max(CH4_wfps1, CH4_wfps0 + 0.2)

  end subroutine read_vegn_namelist

  !----------------------------------------------------------------
  subroutine read_soil_namelist(fnml)
    character(len=*),intent(in) :: fnml
    !--------local vars -----------
    integer :: rc, fu

    ! Open and read Namelist file.
    open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
    read (nml=soil_data_nml, iostat=rc, unit=fu)
    if (rc == 0) then
      write(*,*)'Namelist soil_data_nml read successfully.'
    else
      write(*,*)'Namelist soil_data_nml error', rc
      stop
    endif
    !write(*,nml=soil_data_nml)
    close (fu)

  end subroutine read_soil_namelist

  !----------------------------------------------------------------
  subroutine read_global_setting(fnml)
    character(len=*),intent(in) :: fnml
    !--------local vars -----------
    integer :: rc, fu
    ! Open the namelist file and read global_setting_nml
    open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
    read (nml=global_setting_nml, iostat=rc, unit=fu)
    if (rc == 0) then
      write(*,*)'Namelist global_setting_nml read successfully.'
    else
      write(*,*)'Namelist global_setting_nml error', rc
      stop
    endif
    !write(*,nml=initial_state_nml)
    close (fu)
    ! In case the values are not in reasonable ranges
    LowerLon = max(LowerLon, 1)
    UpperLon = min(UpperLon, Nlon)
    LowerLat = max(LowerLat, 1)
    UpperLat = min(UpperLat, Nlat)
  end subroutine read_global_setting

  !=============================================================================
  ! For global testing of tree-grass-desert shrub and evergreen-deciduous forests
  ! Eight PFTs: 0: C4G, 1: C3G, 2: TrE, 3: TrD, 4: TmE, 5: TmD, 6: Nfx, 7: DeS
  ! Weng, 09/06/2025
  subroutine Preset_GlobalPFTs()
    implicit none
    !------------------------ Pre-setup global PFT parameters --------------------
    !---------------------0: C4G, 1: C3G, 2: TrE, 3: TrD, 4: TmE, 5: TmD, 6: Nfx, 7: DeS
    pt(0:N_EST)        = [1,      0,      0,      0,      0,      0,      0,      0     ] ! 0 for C3, 1 for C4
    phenotype(0:N_EST) = [0,      0,      1,      0,      1,      0,      0,      0     ] ! 0: Deciduous, 1: evergreen
    lifeform(0:N_EST)  = [0,      0,      1,      1,      1,      1,      1,      1     ] ! life form of PFTs: 0 for grasses, 1 for trees
    s0_plant(0:N_EST)  = [0.005,  0.005,  0.02,   0.02,   0.02,   0.02,   0.02,   0.01  ] ! kgC, initial seedling size
    LAImax(0:N_EST)    = [2.5,    2.5,    4.5,    4.5,    3.5,    3.5,    3.0,    2.0   ] ! maximum LAI for a tree
    LMA(0:N_EST)       = [0.025,  0.02,   0.05,   0.025,  0.14,   0.02,   0.02,   0.03  ] ! leaf mass per unit area, kg C/m2
    LNbase(0:N_EST)    = [1.0E-3, 1.0E-3, 0.9E-3, 1.2E-3, 0.9E-3, 1.2E-3, 0.9E-3, 1.0E-3] !functional nitrogen per unit leaf area, kg N/m2, 1.1E-3 for Acer, 1.5E-3 for Populus
    alphaHT(0:N_EST)   = [30.,    30.,    36.,    36.,    36.,    36.,    36.,    18.   ] ! Height parameter (HT = alphatHT * D^0.5)
    alphaCA(0:N_EST)   = [120.,   120.,   150.,   150.,   150.,   150.,   150.,   120.  ] ! Crown area parameter (CA = alphaCA * D^1.5)
    phiRL(0:N_EST)     = [1.5,    1.5,    1.5,    1.5,    1.2,    1.2,    1.8,    1.5   ] ! ratio of fine root area to leaf area
    tauNSC(0:N_EST)    = [2.0,    2.0,    1.5,    1.5,    1.5,    1.5,    1.5,    2.0   ] ! NSC residence time,years
    m_cond(0:N_EST)    = [7.0,    9.0,    9.0,    9.0,    9.0,    9.0,    9.0,    9.0   ] !
    rho_wood(0:N_EST)  = [90.,    90.,    320.,   320.,   330.,   350.,   280.,   450.  ] ! kgC m-3
    mu0_topL(0:N_EST)  = [.02,    .02,    .025,   .025,   .02,    .015,   .06,    .01   ] ! 0.01 ! yearly ! 0.012 for Acer, 0.0274 for Populus
    D0mu(0:N_EST)      = [0.0,    0.0,    0.8,    0.8,    1.2,    1.2,    0.5,    0.25  ] ! m, Mortality curve parameter
    A_sd(0:N_EST)      = [0.0,    0.0,    8.0,    8.0,    8.0,    8.0,    8.0,    2.0   ] ! Max multiplier for seedling mortality
    B_sd(0:N_EST)      = [-60.,   -60.,   -25.,   -25.,   -25.,   -25.,   -25.,   -40.  ] ! Mortality sensitivity for seedlings
    A_DBH(0:N_EST)     = [4.0,    4.0,    4.0,    4.0,    4.0,    4.0,    4.0,    4.0   ] ! Max multiplier for DBH-based mortality
    B_DBH(0:N_EST)     = [.125,   .125,   .125,   .125,   .125,   .125,   .125,   .125  ] ! Size-based Mortality sensitivity, m
    W_mu0(0:N_EST)     = [2.5,    2.5,    0.65,   0.85,   0.85,   0.85,   0.75,   1.5   ] ! Jeremy's half-mortality transp deficit, high:0.5, low: 0.75, No effects: 2.5
    IgniteP(0:N_EST)   = [1.0,    1.0,    .01,    .02,    .02,    .02,    .02,    .02   ] ! Intrinsic flammability
    mu0fire(0:N_EST)   = [0.2,    0.2,    .80,    .99,    .99,    .80,    .90,    .60   ] ! Intrinsic flammability
    gamma_SW(0:N_EST)  = [0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02  ] ! Wood Acambium respiration rate (kgC/m2/yr
    Tc0_OFF(0:N_EST)   = [12.,    8.,     15.,    15.,    -50.,   12.,    15.0,   12.   ] ! 283.16 ! OFF ! C for convenience
    Tc0_ON(0:N_EST)    = [10.,    5.,     10.,    10.,    -50.,   8.,     10.,     8.   ] ! 280.16 ! ON  ! C for convenience
    betaON(0:N_EST)    = [0.3,    0.3,    0.0,    0.6,    0.0,    0.3,    0.6,    0.2   ] ! Critical soil moisture for phenology ON
    betaOFF(0:N_EST)   = [0.2,    0.2,    0.0,    0.4,    0.0,    0.2,    0.4,    0.1   ] ! Critical soil moisture for phenology OFF
    gdd_par1(0:N_EST)  = [50.,    20.,    0.0,    0.0,    0.0,    50.,    50.,    50.   ] ! 50.d0   ! These three parameters are used to calculate gdd_crit
    gdd_par2(0:N_EST)  = [800.,   600.,   0.0,    600.,   0.0,    600.,   600.,   600.  ] ! 650.d0  !800.d0  ! 638.d0
    gdd_par3(0:N_EST)  = [-0.02,  -0.02,  -0.02,  -0.02,  -0.02,  -0.02,  -0.02,  -0.02 ] ! -0.01d0
    R0_Nfix(0:N_EST)   = [0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.03,   0.0   ] ! Nitrogen fixation rate, 0.03 kgN kgRootC-1 yr-1
    S_facuN(0:N_EST)   = [0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.2,    0.0   ] ! Intensity of faculative N fixation
    ! Not used in current model setting (Global ESS PFTs)
    gdd_crit(0:N_EST)  = [300.,   300.,   300.,   300.,   300.,   300.,   300.,   300.  ] ! 280.0 !
    s_hu(0:N_EST)      = [-25.0,  -25.0,  -25.0,  -25.0,  -25.0,  -25.0,  -25.0,  -25.0 ] ! hydraulic mortality sensitivity
    AWD_crit(0:N_EST)  = [0.3,    0.3,    0.7,    0.7,    0.7,    0.7,    0.7,    0.2   ] ! Critical plant water availability factor (0~1)

    write(*,*)"GlobalPFT parameters implemented for Biome ESS"
  end subroutine Preset_GlobalPFTs

  !========================Parameter initialization =========================
  subroutine initialize_soilpars(fnml)
    character(len=*),intent(in) :: fnml

    ! Read in parameters in soil_data_nml
    call read_soil_namelist(fnml)
    ! initialize soil parameters
    soilpars%GMD         = GMD ! geometric mean partice diameter, mm
    soilpars%GSD         = GSD ! geometric standard deviation of particle size
    soilpars%vwc_sat     = vwc_sat
    soilpars%k_sat_ref   = k_sat_ref ! hydraulic conductivity of saturated soil, kg/(m2 s)
    soilpars%psi_sat_ref = psi_sat_ref ! saturation soil water potential, m
    soilpars%chb         = chb       ! Soil texture parameter
    soilpars%alpha       = alphaSoil       ! *** REPLACE LATER BY alpha(layer)
    soilpars%heat_capacity_dry = heat_capacity_dry

    ! ---- derived constant soil parameters
    ! w_fc (field capacity) set to w at which hydraulic conductivity equals
    ! a nominal drainage rate "rate_fc". w_wilt set to w at which psi is psi_wilt
    soilpars%vwc_fc = soilpars%vwc_sat !&
    !*(rate_fc/(soilpars%k_sat_ref*soilpars%alpha**2))**(1/(3+2*soilpars%chb))
    soilpars%vwc_wilt = soilpars%vwc_sat *0.12 !&
    !*(soilpars%psi_sat_ref/(psi_wilt*soilpars%alpha))**(1/soilpars%chb)
    soilpars%vlc_min = soilpars%vwc_sat*K_rel_min**(1/(3+2*soilpars%chb))

    ! Original LM3PPA codes
    !  soil%w_wilt(:) = soil%pars%vwc_sat &
    !       *(soil%pars%psi_sat_ref/(psi_wilt*soil%pars%alpha))**(1/soil%pars%chb)
    !  soil%w_fc  (:) = soil%pars%vwc_sat &
    !       *(rate_fc/(soil%pars%k_sat_ref*soil%pars%alpha**2))**(1/(3+2*soil%pars%chb))
    !  soil%pars%vwc_wilt = soil%w_wilt(1)
    !  soil%pars%vwc_fc   = soil%w_fc  (1)
    !  soil%pars%vlc_min = soil%pars%vwc_sat*K_rel_min**(1/(3+2*soil%pars%chb))

  end subroutine initialize_soilpars

  ! ================================================
  subroutine initialize_PFT_pars(fnml)
    character(len=*),intent(in) :: fnml

    ! ---- local vars
    integer :: i

    ! Update parameters in vegn_parameters_nml
    call read_vegn_namelist(fnml)

    ! Update derived vegn & soil BGC parameters
    etaN = max(etaN,1.0e-6)

    ! initialize vegetation data structure
    spdata%pt       = pt
    spdata%phenotype= phenotype
    spdata%f_cGap   = f_cGap
    spdata%Vmax     = Vmax
    spdata%Vannual  = Vannual
    spdata%m_cond   = m_cond
    spdata%alpha_ps = alpha_ps
    spdata%ps_wet   = ps_wet
    spdata%LFR_rate = LFR_rate
    spdata%gamma_L  = gamma_L
    spdata%gamma_LN = gamma_LN
    spdata%gamma_SW = gamma_SW
    spdata%gamma_FR = gamma_FR
    spdata%rho_wood = rho_wood
    spdata%f_taper  = f_taper
    spdata%rho_FR   = rho_FR
    spdata%root_r   = root_r
    spdata%root_zeta= root_zeta
    spdata%root_perm= root_perm
    spdata%Kw_root  = Kw_root
    !  spdata%rho_N_up0 = rho_N_up0
    !  spdata%N_roots0  = N_roots0

    ! Plant traits
    spdata%LMA      = LMA      ! leaf mass per unit area, kg C/m2
    spdata%LNbase   = LNbase   ! Basal leaf nitrogen per unit area, kg N/m2
    spdata%CN0leafST= CN0leafST    ! Supportive tissues
    spdata%lifeform = lifeform
    spdata%leaf_size= leaf_size
    spdata%alphaHT  = alphaHT
    spdata%thetaHT  = thetaHT
    spdata%alphaCA  = alphaCA
    spdata%thetaCA  = thetaCA
    spdata%alphaBM  = alphaBM
    spdata%thetaBM  = thetaBM

    spdata%AgeRepro = AgeRepro
    spdata%v_seed   = v_seed
    spdata%s0_plant = s0_plant
    spdata%prob_g   = prob_g
    spdata%prob_e   = prob_e
    spdata%mu0_topL = mu0_topL
    spdata%D0mu     = D0mu
    spdata%A_un     = A_un
    spdata%A_sd     = A_sd
    spdata%B_sd     = B_sd
    spdata%A_DBH    = A_DBH
    spdata%B_DBH    = B_DBH
    spdata%s_hu     = s_hu
    spdata%W_mu0    = W_mu0
    spdata%kx0      = kx0
    spdata%WTC0     = WTC0
    spdata%CR_Leaf  = CR_Leaf
    spdata%CR_Wood  = CR_Wood
    spdata%psi0_LF  = psi0_LF
    spdata%psi0_WD  = psi0_WD
    spdata%psi50_WD = psi50_WD
    spdata%Kexp_WD  = Kexp_WD
    spdata%f_plc    = f_plc

    LAI_light = LAImax
    spdata%LAImax   = LAImax
    spdata%LAImax_u = LAImax/3.0
    spdata%LAI_light= LAI_light
    spdata%tauNSC   = tauNSC
    spdata%fNSNmax  = fNSNmax
    spdata%transT   = transT
    spdata%phiRL    = phiRL
    spdata%phiCSA   = phiCSA
    spdata%alpha_FR = alpha_FR  ! root urnover rate

    !! Nitrogen Weng 2012-10-24
    !spdata%CNleaf0 = CNleaf0
    spdata%CNwood0  = CNwood0
    spdata%CNroot0  = CNroot0
    spdata%CNseed0  = CNseed0
    spdata%f_N_add  = f_N_add
    spdata%R0_Nfix  = R0_Nfix
    spdata%C0_Nfix  = C0_Nfix ! scalar, not an array, same cost for all PFTs
    spdata%S_facuN  = S_facuN

    ! Phenology
    spdata%Tc0_OFF  = Tc0_OFF ! C
    spdata%Tc0_ON   = Tc0_ON  ! C
    spdata%gdd_crit = gdd_crit
    spdata%gdd_par1 = gdd_par1
    spdata%gdd_par2 = gdd_par2
    spdata%gdd_par3 = gdd_par3
    spdata%betaON   = betaON
    spdata%betaOFF  = betaOFF
    spdata%AWD_crit = AWD_crit

    ! Plant flammability
    spdata%IgniteP = IgniteP
    spdata%mu0fire = mu0fire

    !write(*,*)'  kx0,    WTC0,    CR_Wood,    psi50_WD,    psi0_WD,    Kexp_WD,    f_supply,    mu0_topL'
    do i = 0, MSPECIES
      call init_derived_species_data(spdata(i))
    enddo
  end subroutine initialize_PFT_pars

  !------------------------------------------
  subroutine init_derived_species_data(sp)
    type(spec_data_type), intent(inout) :: sp
    ! ---- local vars ------
    integer :: j
    real :: rdepth(0:soil_L)
    real :: residual
    real :: R_WD

    ! specific root area
    sp%SRA = 2.0/(sp%root_r*sp%rho_FR) ! m2/kgC
    ! root vertical profile
    rdepth=0.0
    do j=1,soil_L
      rdepth(j) = rdepth(j-1)+thksl(j)
      sp%root_frac(j) = exp(-rdepth(j-1)/sp%root_zeta)- &
      exp(-rdepth(j)  /sp%root_zeta)
    enddo
    residual = exp(-rdepth(soil_L)/sp%root_zeta)
    do j=1,soil_L
      sp%root_frac(j) = sp%root_frac(j) + residual*thksl(j)/rdepth(soil_L)
    enddo

    ! calculate alphaBM parameter of allometry. note that rho_wood was re-introduced for this calculation
    sp%alphaBM    = sp%rho_wood * sp%f_taper * PI/4. * sp%alphaHT ! 5200

    !  Vmax as a function of LNbase
    sp%Vmax = 0.02 * sp%LNbase ! 0.03125 * sp%LNbase ! Vmax/LNbase= 25E-6/0.8E-3 = 0.03125 !
    !  CN0 of leaves
    sp%LNA     = sp%LNbase +  sp%LMA/sp%CN0leafST
    sp%CNleaf0 = sp%LMA/sp%LNA
    !  Leaf life span as a function of LMA
    sp%leafLS = c_LLS * sp%LMA
    if(sp%leafLS>1.0)then
      sp%phenotype = 1 ! Everygreen
    else
      sp%phenotype = 0 ! Deciduous
    endif
    !  Leaf turnover rate, (leaf longevity as a function of LMA)
    sp%alpha_L = 1.0/sp%leafLS * sp%phenotype
    !sp%leafTK = 4.0e-4 * SQRT(sp%LMA/0.02) ! Niinemets 2001, Ecology
    sp%leafTK = TK0_leaf * SQRT(sp%LMA) ! Leaf thickness, m
    sp%rho_leaf= sp%LMA/sp%leafTK
    sp%CR_Leaf = CR0_LF  * (0.02/sp%LMA)
    ! Wood hydraulic traits as functions of wood density, 06/30/2022, Weng
    R_WD        = sp%rho_wood/WDref0
    sp%kx0      = kx0_WD  * R_WD**(-1)  ! (mm/s)/(Mpa/m)
    sp%WTC0     = WTC0_WD * R_WD**1.5
    sp%CR_Wood  = CR0_WD  * R_WD**(-1.67)  ! Compress ratio per MPa, Santiago et al. 2018
    sp%psi50_WD = p50_WD  * R_WD**1.73 - 1.0 !- 1.09 - 3.57 * (sp%rho_wood/500.) ** 1.73
    sp%psi0_WD  = p50_WD  * R_WD**1.73 - 2.0
    sp%Kexp_WD  = ths0_WD * R_WD**(-0.5) ! rho_wood cannot be 200 for ths0_WD * R_WD**(-1)
    sp%f_supply = fsup0_WD/(R_WD+1.0)
    sp%f_plc    = fplc0_WD * exp(A_plc0_WD*R_WD)

    ! Mortality rate as a function of wood density
    if(Do_mu0_F_WDen)then
      sp%mu0_topL = 0.048 - 0.024 * R_WD
    endif
    !write(*,'(40(F10.4,","))')sp%kx0,sp%WTC0,sp%CR_Wood,sp%psi50_WD,sp%psi0_WD,sp%Kexp_WD,sp%f_supply,sp%mu0_topL

    ! -------- Check parameter boundaries
    if(sp%lifeform==0)then
      sp%transT = dt_daily_yr ! Grass BLmax transition in canopy layer is short
    else
      sp%transT = max(sp%transT, dt_daily_yr) !Minimum for woody plants
    endif

  end subroutine init_derived_species_data

  !=============================================================================
  subroutine Set_PFTs_from_map(fPFT)
    real, intent(in) :: fPFT(:)

    !--------- local vars ------------
    integer :: GridPFTs(N_PFTs)
    integer :: i
    real :: f_min = 0.01 ! coverage fraction threshold

    ! Sorting PFT numbers according to fPFT
    call rank_descending(fPFT,GridPFTs)
    !PFTID = [character(len=3) :: 'C4G','C3G','TEB','TDB','EGN','CDB','CDN','CAS','AAS']
    GridPFTs = GridPFTs - 1 ! PFT No. starts from 0.

    ! Find out PFTs in this grid
    init_cohort_N = min(M_initialCH,Max(1, COUNT(fPFT > f_min)))
    do i=1, init_cohort_N
      init_cohort_sps(i)   = GridPFTs(i)
      init_cohort_Indiv(i) = 0.2  ! initial individual density, individual/m2
      init_cohort_bsw(i)   = 0.01 ! initial biomass of sapwood, kg C/individual
      init_cohort_nsc(i)   = 0.01 ! initial non-structural biomass, kg C/individual
    enddo

    ! Initial soil Carbon and Nitrogen for a vegn tile, Weng 2012-10-24
    init_fast_SOC  = 0.5  ! initial fast soil C, kg C/m2
    init_slow_SOC  = 20.0 ! initial slow soil C, kg C/m2
    init_mineralN  = 0.02 ! Mineral nitrogen pool, (kg N/m2)
    N_input        = 2.0  ! 0.0008 ! annual N input to soil N pool, kgN m-2 yr-1
  end subroutine Set_PFTs_from_map

  !=============================================================================
  ! Calculate vars of climate envelopes
  subroutine Climate_envelope_vars (forcingData,steps_per_day)
    implicit none
    type(climate_data_type), intent(in) :: forcingData(:)
    integer, intent(in) :: steps_per_day

    !--------- local vars ------------
    integer, parameter :: mw = 15
    real, parameter :: pet_floor = 1.0e-6
    real, allocatable :: dailyET(:),dailyTc(:),dailyPr(:),meanTc(:),TminYr(:)
    real :: tmpL(mw*2+1)
    real :: totPrcp, totPET
    integer ::  N_Yrs, N_days,i,j,k,m,n,L,w

    ! ------ Calculate days and years of the data ------
    N_days = size(forcingData)/steps_per_day
    N_yrs = N_days/365
    if (N_yrs < 1) then
      error stop "Set_PFTs_from_Climate: need at least one year's data."
    endif

    ! Allocate variables
    allocate(dailyET(N_days),dailyTc(N_days),dailyPr(N_days),meanTc(N_days))
    allocate(TminYr(N_yrs))

    ! Calculate daily temperature and yearly mean Precipitation
    dailyTc = 0.0
    dailyPr = 0.0
    dailyET = 0.0
    totPrcp = 0.0
    totPET  = 0.0
    k = 0
    do i = 1, N_days
      do j = 1, steps_per_day
        k = k + 1
        dailyTc(i) = dailyTc(i) + forcingData(k)%Tair - 273.15
        dailyPr(i) = dailyPr(i) + forcingData(k)%rain * step_seconds
        dailyET(i) = dailyET(i) + potentialET(forcingData(k)) * step_seconds
      enddo
      dailyTc(i) = dailyTc(i)/steps_per_day
      totPrcp = totPrcp + dailyPr(i)
      totPET  = totPET  + dailyET(i)
    enddo
    Mst_IDX  = totPrcp / max(totPET, pet_floor) ! Moisture Index
    meanPrcp = totPrcp / N_yrs
    meanPET  = totPET / N_yrs

    ! Calculate running mean temperature (meanTc) with a window of L
    L = mw * 2 + 1
    do i = 1, N_days
      m = i - mw
      n = i + mw
      if (m >= 1 .and. n <= N_days)then
        tmpL(1:L) = dailyTc(m:n)
      else if (m < 1) then
        w = 1 - m                     ! number from the end
        tmpL(1:w) = dailyTc(N_days+m:N_days)
        tmpL(w+1:L) = dailyTc(1:n)
      else ! n > N_days: take tail then wrap to start
        w = L - (n - N_days)      ! number we can take before hitting the end
        tmpL(1:w) = dailyTc(m:N_days)
        tmpL(w+1:L) = dailyTc(1:n-N_days)
      endif
      meanTc(i) = sum(tmpL)/real(size(tmpL))
    enddo

    ! Find out yearly minimum Tm and calculate yearly mean
    do i = 1, N_yrs
      TminYr(i) = 999.0
      do j = 1, 365
        k = (i-1)*365 + j
        TminYr(i) = min(meanTc(k),TminYr(i))
      enddo
    enddo
    meanTmin = sum(TminYr)/real(size(TminYr))

    !! Write out
    !write(99,'(6(a8,","))')'i','j','dailyET','dailyTC','meanTc','dailyPr'
    !do i = 1, N_yrs
    !   do j = 1, 365
    !     k = (i-1)*365 + j
    !     write(99,'(2(I8,","),4(f12.4,","))')i,j,dailyET(k),dailyTC(k),meanTc(k),dailyPr(k)
    !   enddo
    !enddo

    ! Screen output
    write(*,'(4(a6,f8.2,";"))') 'Prcp: ', meanPrcp, 'PET: ', meanPET, &
                                'P/ET: ', Mst_IDX,  'Tmin: ', meanTmin

    ! Release allocatable variables
    deallocate(dailyET,dailyTc,dailyPr,meanTc,TminYr)
  end subroutine Climate_envelope_vars

  !=============================================================================
  ! An interface for setting PFTs from climate or land cover data, Weng, 04/08/2026
  subroutine Set_PFTs_from_Data()
    implicit none

    !--------- local vars ------------
    !PFTs: 0:C4G, 1:C3G, 2:TrE, 3:TrD, 4:TmE, 5:TmD, 6:Nfx, 7:DeS
    integer,allocatable :: PFTID(:)
    integer :: N_PFTID

#ifdef WIEMIP_setting
    call Set_PFTs_from_LandCover (N_PFTID,PFTID)
#else
    call Set_PFTs_from_Climate (N_PFTID,PFTID)
#endif

    ! Assign standard initial cohorts by updating init_cohort_*
    call Assign_Std_Cohorts (PFTID,N_PFTID)

    deallocate(PFTID)
  end subroutine Set_PFTs_from_Data

  !=============================================================================
  ! for testing tree-grass-desert shrub and evergreen-deciduous forests only
  ! Weng, 09/06/2025
  subroutine Set_PFTs_from_Climate (N_PFTID,PFTID)
    implicit none
    integer, intent(out) :: N_PFTID
    integer, allocatable, intent(out) :: PFTID(:)

    !--------- local vars ------------

    ! Assign allocatable PFTID
    N_PFTID = 4
    allocate(PFTID(N_PFTID))

    ! Assign PFT groups according to climate data at each grid
    !if(meanPrcp > Pr_thld)then
    if(Mst_IDX > MI0DeSB) then
      if(meanTmin > TcrTREE)then ! Tropical vs. Temperate trees
        PFTID = [1,2,3,6]
      else
        PFTID = [1,4,5,6]
      endif
    else  ! Mst_IDX <= MI0DeSB
      ! Replace evergreen with desert shrubs
      if(meanTmin > TcrTREE)then ! Tropical vs. Temperate trees
        PFTID = [0,3,6,7]
      else
        PFTID = [1,5,6,7]
      endif
    endif
    ! Replace C3 with C4 grasses in dry and warm climates
    if(Mst_IDX < MI0C3C4 .and. meanTmin > TcrC3C4) PFTID(1) = 0

    ! Screen output
    write(*,'(a12, 4(I6,","))')'Grid PFTs: ', PFTID
  end subroutine Set_PFTs_from_Climate

  !=============================================================================
  ! for WIEMIP project, Weng, 04/08/2026
  subroutine Set_PFTs_from_LandCover(N_PFTID,PFTID)
    ! Two options: 1 ESACCI PFT map (2011), 2 Hurtt Land cover
    implicit none
    integer, intent(out) :: N_PFTID
    integer, allocatable, intent(out) :: PFTID(:)

    !--------- local vars ------------
    integer :: idx
    logical :: is_crop = .False.
    logical :: use_Hurtt_data = .True.

    ! Assign PFT groups according to climate and land cover data at each grid
    idx = maxloc(GridVC, dim=1)
    write(*,*)'Max vegetation id', idx

    ! Igor & Paul: "GridFR" is from Hurtt's cropland file (states4.nc).
    if(use_Hurtt_data)then
      is_crop = GridFR(LC_year0) > 0.5
    else
      is_crop = idx >=9 ! C3 or C4 grasses
    endif
    ! Assign PFTs
    if (is_crop) then ! Use land cover data for crops
      N_PFTID = 1
      allocate(PFTID(N_PFTID))
      if(meanTmin > TcrC3C4)then ! C3 vs C4 grasses
        PFTID = [0]
      else
        PFTID = [1]
      endif
    else  ! Assign PFTs according to climate data
      N_PFTID = 2
      allocate(PFTID(N_PFTID))
      if(Mst_IDX > MI0DeSB) then
        if(meanTmin > TcrTREE)then ! Tropical vs. Temperate trees
          PFTID = [2,3]
        else
          PFTID = [4,5]
        endif
      else  ! Mst_IDX <= MI0DeSB
        ! Replace evergreen with desert shrubs
        if(meanTmin > TcrTREE)then ! Tropical vs. Temperate trees
          PFTID = [0,7]
        else
          PFTID = [1,7]
        endif
      endif
      ! Replace C3 with C4 grasses in dry and warm climates
      if(Mst_IDX < MI0C3C4 .and. meanTmin > TcrC3C4) PFTID(1) = 0
    endif
    write(*,'(a12, 4(I6,","))')'Grid PFTs: ', PFTID
  end subroutine Set_PFTs_from_LandCover

  !=============================================================================
  ! Assign a standard cohort for each PFT
  subroutine Assign_Std_Cohorts (PFTID,totPFT)
    implicit none
    integer, intent(in) :: PFTID(:), totPFT

    ! Assign initial cohorts
    ! 0:C4, 1:C3, 2:TrE, 3:TrD, 4:TmE, 5:TmD, 6:N-fixer, 7:Desert shrub
    init_cohort_N = totPFT
    init_cohort_sps(1:totPFT)   = PFTID(1:totPFT)
    init_cohort_Indiv(1:totPFT) = std_den(PFTID(1:totPFT))
    init_cohort_bsw(1:totPFT)   = std_bsw(PFTID(1:totPFT))
    init_cohort_nsc(1:totPFT)   = std_nsc(PFTID(1:totPFT))

  end subroutine Assign_Std_Cohorts

  ! ============================================================
  subroutine qscomp(T, p, qsat)
    real, intent(in) :: T    ! temperature, degK
    real, intent(in) :: p    ! pressure, Pa
    real, intent(out):: qsat ! saturated specific humidity, kg/kg
    !--------local var
    real :: esat ! sat. water vapor pressure
    real :: Temp ! degC

    ! calculate saturated specific humidity
    Temp = T - 273.16 ! degC
    esat=MIN(610.78*exp(17.27*Temp/(Temp+237.3)), p) ! Pa
    qsat = 0.622*esat /(p - 0.378*esat )
  end subroutine qscomp

  !===========================
  FUNCTION esat(T) ! pressure, Pa
    IMPLICIT NONE
    REAL :: esat
    REAL, INTENT(IN) :: T ! degC
    esat=610.78*exp(17.27*T/(T+237.3))
  END FUNCTION esat

  ! ============================================================================
  subroutine BM2Architecture(cc,BM)
    type(cohort_type), intent(inout) :: cc
    real,intent(in)    :: BM ! wood biomass, kg C
    !------------------------------------
    associate(sp=>spdata(cc%species))
      cc%DBH    = BM2DBH(   BM,cc%species)
      cc%height = DBH2HT(cc%DBH,cc%species)
      cc%Acrown = DBH2CA(cc%DBH,cc%species)
      cc%bl_max = CA2BLmax(cc) !sp%LMA  * sp%LAImax * cc%Acrown * (1.0-sp%f_cGap)/max(1,cc%layer)
      cc%br_max = BLmax2BRmax(cc) !sp%phiRL* cc%bl_max/(sp%LMA*sp%SRA)
      cc%NSNmax = ccNSNmax(cc) ! sp%fNSNmax*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
    end associate
  end subroutine BM2Architecture

  !============================================================================
  ! Allometry equations, Weng (2014-01-09),07-18-2017, 06-04-2021,2022-08-20
  !-------------------------------------------
  function DBH2HT(DBH,SP) result (HT)
    real :: HT ! returned value
    real,   intent(in) :: DBH
    integer,intent(in) :: SP
    HT = spdata(SP)%alphaHT * DBH ** spdata(SP)%thetaHT
  end function
  !-------------------------------------------
  function DBH2CA(DBH,SP) result (CA)
    real :: CA ! returned value
    real,   intent(in) :: DBH
    integer,intent(in) :: SP
    CA = spdata(SP)%alphaCA * DBH ** spdata(SP)%thetaCA
  end function
  !-------------------------------------------
  function DBH2BM(DBH,SP) result (BM)
    real :: BM ! returned value
    real,   intent(in) :: DBH
    integer,intent(in) :: SP
    BM = spdata(SP)%alphaBM * DBH ** spdata(SP)%thetaBM
  end function
  !-------------------------------------------
  function BM2DBH(BM,SP) result (DBH)
    real :: DBH ! returned value
    real,   intent(in) :: BM
    integer,intent(in) :: SP
    DBH = (BM/spdata(SP)%alphaBM) ** ( 1.0/spdata(SP)%thetaBM )
  end function

  !-------------------------------------------
  function ccNSNmax(cc) result (NSNmax)
    real :: NSNmax ! returned value
    type(cohort_type), intent(in) :: cc    ! cohort to update

    associate(sp=>spdata(cc%species))
      NSNmax = sp%fNSNmax*(cc%bl_max/(sp%CNleaf0*sp%leafLS)+cc%br_max/sp%CNroot0)
    end associate
  end function

  !-------------------------------------------
  function CA2BLmax(cc) result (BLmax)
    real :: BLmax ! returned value
    type(cohort_type), intent(in) :: cc    ! cohort to update

    associate(sp=>spdata(cc%species))
      BLmax = sp%LMA  * sp%LAImax * cc%Acrown * (1.0-sp%f_cGap)/max(1,cc%layer)
#ifdef UFL_test
      if(cc%layer > 1) BLmax = sp%LMA  * 1.0 * cc%Acrown * (1.0-sp%f_cGap)
#endif
    end associate
  end function

  !-------------------------------------------
  function BLmax2BRmax(cc) result (BRmax)
    real :: BRmax ! returned value
    type(cohort_type), intent(in) :: cc    ! cohort to update

    associate(sp=>spdata(cc%species))
      BRmax = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
#ifdef UFL_test
      if(cc%layer > 1) BRmax = 0.5*cc%bl_max/(sp%LMA*sp%SRA)
#endif
    end associate
  end function

  !-------------------------------------------
  function BL2Aleaf(bl,cc) result (area)
    real :: area ! returned value
    real, intent(in) :: bl      ! biomass of leaves, kg C/individual
    type(cohort_type), intent(in) :: cc    ! cohort to update

    area = bl/spdata(cc%species)%LMA
  end function

  !-------------------------------------------
  function Aleaf2LAI(CAI_L,cc) result (LAI)
    ! Squeeze leaves within canopy when CAI > 1.0. Weng, 12/04/2025
    ! Since we do realyering once per year, it is possible CAI>1 with growth
    real :: LAI ! returned value
    real, intent(in) :: CAI_L           ! biomass of leaves, kg C/individual
    type(cohort_type), intent(in) :: cc ! cohort to update
    LAI = Max(CAI_L,1.0)*cc%Aleaf/(cc%Acrown*(1.0-spdata(cc%species)%f_cGap))
  end function

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
    WD = 1.0/(1.0 + exp(12.0*(0.5-thetaS))) ! Weng, 2025-02-11

    A = (Td*Wd)  ! the combined (multiplicative) effect of temp and water
    ! on decomposition rates
  end function A_function

  !============================================================================
  ! Weng, 05/24/2018
  subroutine calc_solarzen(td,latdegrees,cosz,solarelev,solarzen)
    !* Calculate solar zenith angle **in radians**
    !* From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
    implicit none
    real,intent(in) :: td            ! day(to minute fraction)
    real,intent(in) :: latdegrees    ! latitude in degrees
    real :: hour,latrad
    real :: delta    ! declination angle
    real :: pi, rad
    real,intent(out) :: cosz         ! cosz=cos(zen angle)=sin(elev angle)
    real,intent(out) :: solarelev    ! solar elevation angle (rad)
    real,intent(out) :: solarzen     ! solar zenith angle (rad)
    pi  = 3.1415926
    rad = pi / 180.0 ! Conversion from degrees to radians.
    hour = (td-floor(td))*24.0
    latrad = latdegrees*rad
    delta  = asin(-sin(rad*23.450)*cos(2.0*pi*(td+10.0)/365.0))
    cosz = sin(latrad)*sin(delta) + &
    cos(latrad)*cos(delta)*cos(rad* 15.0*(hour-12.0))
    cosz = max (cosz, 0.0)  ! Sun's angular is 0.01

    ! compute the solar elevation and zenth angles below
    solarelev = asin(cosz)/pi*180.0  !since asin(cos(zen))=pi/2-zen=elev
    solarzen = 90.0 - solarelev ! pi/2.d0 - solarelev
  end subroutine calc_solarzen

  !=================================================================
  ! Weng, 2025-09-07
  ! ------- Potential ET -----------
  ! from ChapGPT for short grass (FAO reference crop, 0.12 m), raero = 100
  ! Aero conductance ! https://www.fao.org/4/x0490e/x0490e06.htm
  ! ChatGPT: Short grass (FAO reference crop, 0.12 m): ≈100 s m⁻¹.
  ! Tall crops (e.g., maize, ~2 m): 20-50 s m-1 because taller, rougher canopies enhance turbulence.
  ! Forests (>20 m): 5–30 s m-1 due to strong turbulence above canopy.
  ! Smooth bare soil / desert: >150–300 s m-1
  function PotentialET(forcing)
    implicit none
    real :: PotentialET  ! returned value, kg H2O/m2/second
    type(climate_data_type),intent(in):: forcing

    !----- local var --------------
    ! Meteorological variables
    real :: Rnet          ! net radiation, and soil surface radiation, W/m2
    real :: RH, Uwind     ! relative humidity (0~1); wind speed (m/s)
    real :: TairK, Tair   ! temperature (K and C)
    real :: P_air, Dair   ! Air pressure and VPD, (pa)
    ! Assumed vegetation states
    real :: Karman = 0.41 ! von Kármán constant (~0.41)
    real :: dV     = 0.1  ! a small number
    real :: Zmh    = 2.0   ! Measurement height (m)
    real :: Zvg    = 0.12  ! FAO, short-crop height (m), for PET
    real :: Z0m    = 0.014 ! roughness length for momentum [m]
    real :: Z0h    = 0.02  ! roughness length for heat and vapour [m]
    real :: dz             ! zero plane displacement height (zero wind height) [m]
    ! For calculation
    real :: rhocp, H2OLv, slope, psyc
    real :: rAero         ! resistance, s m-1
    real :: Esoil         ! Soil surface evaporation, W/m2

    ! Forcing variables
    Rnet  = forcing%radiation * 0.9 ! assuming 10% reflectance
    Uwind = forcing%windU + dV      ! in case U is zero.
    TairK = forcing%Tair
    Tair  = forcing%Tair - 273.16
    P_air = forcing%P_air
    RH    = Max(0.01, min(0.99, forcing%RH))  ! Check forcing's unit of humidity

    ! --------- Potential evapotranspiration (PET) --------
    rhocp = cpair * P_air * mol_air / (Rgas*TairK)
    H2OLv = H2oLv0 - 2.365e3 * Tair
    Dair  = esat(Tair) * (1.0 - RH)
    slope = (esat(Tair + dV) - esat(Tair)) / dV
    psyc  = P_air * cpair * mol_air / (H2OLv*mol_h2o)
    dz    = Zvg*2/3 ! Zero wind speed height
    raero = (log((Zmh - dz)/Z0m)*log((Zmh - dz)/Z0h))/(Karman*Karman*Uwind)
    Esoil = (slope*Rnet + rhocp*Dair/raero)/(slope + psyc)
    !(slope + psyc * (1.0+rsoil/raero)) ! AET, Liqing Peng et al. 2019 GCB

    PotentialET = Esoil/H2OLv ! kg H2O m-2 second-1
  end function PotentialET

  ! ============================================================================
  function c3c4(temp, precip) result (pt)
    integer :: pt
    real,              intent(in) :: temp   ! temperatire, degK
    real,              intent(in) :: precip ! precipitation, ???

    real :: pc4
    ! Rule based on analysis of ED global output; equations from JPC, 2/02
    pc4=exp(-0.0421*(273.16+25.56-temp)-(0.000048*(273.16+25.5-temp)*precip))

    if(pc4>0.5) then
      pt=PT_C4
    else
      pt=PT_C3
    endif
  end function

  !================= utilities ===========================================
  !=======================================================================
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

  !=======================================================================
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

!================================================
end module model_utils
