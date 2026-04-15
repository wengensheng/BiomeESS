module datatypes
  ! define data types, constants, and some basic functions
  implicit none

  !===============constants===============
  integer, parameter :: days_per_year    = 365
  integer, parameter :: hours_per_year   = 365 * 24  ! 8760
  integer, parameter :: MonthDOY(0:12)   =(/0,31,59,90,120,151,181,212,243,273,304,334,366/)
  real,    parameter :: seconds_per_year = 365. * 24. * 3600.
  real,    parameter :: seconds_per_day  = 24. * 3600. ! 86400

  ! Physical constants
  real, parameter :: PI      = 3.1415926
  real, parameter :: rho_H2O = 1000.0 ! water density (kg m-3)
  real, parameter :: TFREEZE = 273.16
  real, parameter :: Rgas    = 8.314472 ! universal gas constant, J K-1 mol-1
  real, parameter :: mol_C   = 12.0e-3 ! molar mass of carbon, kg
  real, parameter :: mol_air = 28.96440e-3 ! molar mass of air, kg
  real, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
  real, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg
  real, parameter :: cpair   = 1010.     ! air heat capapcity (J/kg/K)
  real, parameter :: H2OLv0  = 2.501e6   ! latent heat H2O (J/kg)
  real, parameter :: p_sea   = 101325.  ! atmospheric pressure  (Pa)
  real, parameter :: solarC  = 1361.0   ! Solar constant, W/m2
  real, parameter :: f_PAR   = 0.5  ! Fraction of PAR in total solar radiation
  real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta

  ! Vegetation and soil types
  integer, parameter :: n_dim_soil_types = 9
  integer, parameter :: MSPECIES   = 30 ! 15
  integer, parameter :: N_EST      = 7 ! 8 ESS PFTs since started from 0.
  integer, parameter :: Nsoiltypes = 7
  integer, parameter :: LEAF_ON    = 1
  integer, parameter :: LEAF_OFF   = 0
  integer, parameter :: NBANDS     = 2, & ! number of spectral bands for short-wave radiation calculations
                        BAND_VIS   = 1, & ! visible radiation (wavelenght range?)
                        BAND_NIR   = 2    ! near infra-red radiation (wavelenght range?)

  real,    parameter :: min_nindivs= 0.1E-4 ! 2e-15 ! 1/m. 2e-15 is approximately 1 individual per Earth
  real,    parameter :: zero_thld = 1.0e-12

  ! For vegn_photosynthesis
  integer, parameter :: CLmax   = 5 ! Maximum crown layers
  real, parameter :: light_crit = 1.0e-12
  real, parameter :: gs_lim     = 0.25
  real, parameter :: lai_min    = 1.0e-6
  real, parameter :: lai_max    = 8.0        ! cc%LAI ranges from 0 to ~8
  real, parameter :: ext_min    = 0.10       ! cc%extinct is usually 0.2~0.9
  real, parameter :: ext_max    = 0.90

  ! Plant hydraulics-mortality
  integer, parameter :: Ysw_max      = 210 ! Maximum function years of xylems
  real, parameter    :: WDref0       = 300.0   ! Reference wood density, kgC m-3
  real, parameter    :: rho_cellwall = 750.0 ! kgC m-3, Kellogg & Wangaard 1969 1.5 g/cc
  real, parameter    :: LMAmin       = 0.02    ! minimum LMA (kgC/m2), boundary condition

  ! Soil SOM reference C/N ratios
  integer, parameter :: N_SOM = 5
  !  fineL, structuralL, microbial, fast, slow
  real, parameter :: CN0SOM(5) = (/50., 150., 10., 15., 40./) ! target CN ratios of SOM

  ! Soil water layers
  integer, parameter :: soil_L = 5 ! Soil layers, for soil water dynamics
  integer, parameter :: num_l  = soil_L ! Soil layers,
  real, parameter :: psi_wilt  = -150.0  ! matric head at wilting
  real, parameter :: K_rel_min = 1.e-12
  real, parameter :: rate_fc   = 0.1/86400 ! 0.1 mm/d threshold drainage rate at Field Capacity
  real, parameter :: ws0 = 0.02 ! hygroscopic point
  !real, parameter :: Edepth = 0.05 !m, the depth of soil for surface evaporation
  integer, parameter :: & ! soil types
                     Sand        = 1,  LoamySand   = 2, &
                     SandyLoam   = 3,  SiltLoam    = 4, &
                     FrittedClay = 5,  Loam        = 6, &
                     Clay        = 7

  ! LM3 vegetation and soil types
  integer, parameter :: PHEN_DECIDIOUS = 0, PHEN_EVERGREEN = 1 ! phenology type
  integer, parameter :: PT_C3 = 0, PT_C4 = 1  ! photosynthesis types
  integer, parameter :: NSPECIES = 5,   & ! number of species
                        SP_C4GRASS = 0, & ! c4 grass
                        SP_C3GRASS = 1, & ! c3 grass
                        SP_TEMPDEC = 2, & ! temperate deciduous
                        SP_TROPICAL= 3, & ! non-grass tropical
                        SP_EVERGR  = 4    ! non-grass evergreen
  integer, parameter :: N_LU_TYPES = 4, & ! number of land use types
                        LU_PAST    = 1, & ! pasture
                        LU_CROP    = 2, & ! crops
                        LU_NTRL    = 3, & ! natural vegetation
                        LU_SCND    = 4    ! secondary vegetation
  integer, parameter :: CO2Yrs = 325      ! 325 = the years from 1701 to 2025
  real, parameter :: CO2_Hist(CO2Yrs) = & ! CO2 concentration 1700~2024, ppm
  ! 1,     2,     3,     4,     5,     6,     7,     8,     9,     10
  (/276.59,276.62,276.65,276.67,276.70,276.72,276.75,276.78,276.80,276.83,  &
  276.86,276.89,276.92,276.95,276.98,277.02,277.05,277.09,277.13,277.17,  &
  277.21,277.25,277.30,277.35,277.40,277.45,277.50,277.54,277.58,277.62,  &
  277.64,277.66,277.66,277.65,277.63,277.59,277.54,277.49,277.42,277.35,  &
  277.28,277.22,277.15,277.09,277.03,276.98,276.93,276.89,276.85,276.81,  &
  276.78,276.74,276.71,276.68,276.65,276.63,276.61,276.59,276.58,276.58,  &
  276.58,276.60,276.62,276.66,276.70,276.76,276.83,276.91,276.99,277.09,  &
  277.19,277.29,277.40,277.51,277.63,277.73,277.86,277.98,278.11,278.24,  &
  278.39,278.54,278.72,278.90,279.10,279.30,279.52,279.74,279.97,280.20,  &
  280.44,280.67,280.90,281.13,281.35,281.57,281.78,281.98,282.18,282.37,  &
  282.54,282.71,282.86,283.01,283.14,283.27,283.38,283.49,283.57,283.65,  &
  283.71,283.76,283.79,283.81,283.81,283.80,283.77,283.73,283.69,283.64,  &
  283.60,283.55,283.51,283.48,283.46,283.46,283.47,283.50,283.53,283.58,  &
  283.64,283.70,283.76,283.84,283.92,284.01,284.11,284.23,284.37,284.53,  &
  284.71,284.92,285.14,285.36,285.58,285.79,285.98,286.14,286.28,286.39,  &
  286.46,286.51,286.52,286.52,286.49,286.46,286.44,286.41,286.40,286.40,  &
  286.42,286.45,286.51,286.59,286.71,286.85,287.02,287.21,287.42,287.64,  &
  287.87,288.11,288.34,288.56,288.77,288.97,289.16,289.35,289.56,289.77,  &
  290.00,290.25,290.51,290.80,291.13,291.48,291.86,292.25,292.65,293.05,  &
  293.46,293.86,294.25,294.62,294.98,295.31,295.60,295.86,296.11,296.34,  &
  296.57,296.80,297.04,297.29,297.57,297.87,298.17,298.49,298.82,299.15,  &
  299.49,299.84,300.20,300.58,300.97,301.36,301.76,302.16,302.55,302.92,  &
  303.28,303.63,303.97,304.32,304.66,305.01,305.37,305.73,306.10,306.47,  &
  306.84,307.21,307.58,307.96,308.32,308.69,309.04,309.38,309.70,310.01,  &
  310.27,310.50,310.68,310.82,310.93,311.03,311.12,311.24,311.36,311.52,  &
  311.71,311.92,312.16,312.44,312.75,313.10,313.48,313.93,314.69,315.47,  &
  316.38,317.09,317.88,318.40,319.00,319.40,320.71,321.50,322.35,323.90,  &
  324.94,325.55,326.67,328.87,329.36,330.25,331.13,332.73,334.58,335.81,  &
  337.80,339.31,340.54,342.11,343.64,345.15,346.43,348.25,350.45,352.08,  &
  353.27,354.56,355.50,356.27,357.78,359.71,361.39,362.70,365.44,367.26,  &
  368.42,370.02,372.15,374.62,376.33,378.45,380.51,382.53,384.42,386.05,  &
  388.23,390.03,392.13,394.88,396.94,399.24,402.52,404.71,406.94,409.63,  &
  414.21,416.41,418.53,421.08,424.61/)

  !===============data types ==============================
  !-----------PFT data type----------------
  type spec_data_type
    integer :: lifeform  ! 0 for grasses, 1 for trees
    integer :: phenotype ! phenology type: 0 for deciduous, 1 for evergreen
    integer :: pt        ! photosynthetic physiology of species
    ! leaf traits
    real :: LMA          ! leaf mass per unit area, kg C/m2
    real :: leafLS       ! leaf life span
    real :: alpha_L      ! leaf turn over rate
    real :: LNA          ! leaf Nitrogen per unit area, kg N/m2
    real :: LNbase       ! basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
    real :: CN0leafST    ! leaf structural tissues, 175
    real :: leaf_size    ! characteristic leaf size
    real :: leafTK       ! leaf thickness, m
    real :: rho_leaf     ! leaf mass density (kgC/m3)
    real :: alpha_ps     ! photosynthesis efficiency
    real :: m_cond       ! Stomatal Sensitivity (g1): how "cautious" a plant is. low m indicates a water-conservative strategy
    real :: Vmax         ! max rubisco rate, mol m-2 s-1
    real :: Vannual      ! annual productivity per unit area at full fun (kgC m-2 yr-1)
    real :: gamma_L      ! leaf respiration coeficient (per yr)
    real :: gamma_LN     ! leaf respiration coeficient per unit N
    real :: ps_wet       ! wet leaf photosynthesis down-regulation
    real :: LFR_rate     ! Leaf filling rate per day
    ! radiation parameters for 2 bands, VIS and NIR
    real :: leaf_refl (NBANDS) = (/ 0.10, 0.50/) ! reflectance of leaf
    real :: leaf_tran (NBANDS) = (/ 0.05, 0.25/) ! transmittance of leaf
    real :: leaf_emis          = 1.00            ! emissivity of leaf
    ! root traits
    real :: rho_FR       ! material density of fine roots (kgC m-3)
    real :: root_r       ! radius of the fine roots, m
    real :: root_zeta    ! e-folding parameter (zeta) of root vertical distribution (m)
    real :: root_frac(soil_L)    ! root fraction
    real :: SRA          ! speific fine root area, m2/kg C
    real :: SRL          ! specific root lenght
    real :: gamma_FR     ! Fine root respiration rate, kgC kgC-1 yr-1
    real :: alpha_FR     ! Turnover rate of Fine roots, fraction yr-1
    real :: Kw_root      ! fine root water conductivity mol m m-2 s−1 Pa−1 !
    real :: root_perm    ! assume it is a fraction of root area active for water flow
    !  real :: rho_N_up0   ! maximum N uptake rate
    !  real :: N_roots0    ! root biomass at half of max. N-uptake rate
    real :: R0_Nfix    ! Reference N fixation rate (kgN kgC-1 root)
    real :: C0_Nfix    ! Carbon cost of N fixation (kgC kgN-1)
    real :: S_facuN    ! Intensity of faculative N fixation, 0~1 usage of extraC
    ! wood traits
    real :: rho_wood     ! woody density, kg C m-3 wood
    real :: gamma_SW     ! sapwood respiration rate, kgC m-2 Acambium yr-1
    real :: f_taper

    ! Fire related
    real :: IgniteP ! Probability of ignition when climatic conditions are met (i.e., flammability).
    real :: mu0fire ! PFT-specific fire sensitivity

    ! Plant hydraulics
    real :: kx0  ! xylem conductivity, (mm/s)/(Mpa/m)
    real :: WTC0 ! xylem water transfer capacity, m/lifetime
    real :: CR_Leaf ! leaf compression ratio per MPa
    real :: CR_Wood ! Wood compression ratio per MPa
    real :: psi0_LF ! minimum leaf water potential
    real :: psi0_WD ! minimum stem water potential
    real :: psi50_WD !wood potential at which 50% conductivity lost, MPa
    real :: Kexp_WD  ! exponent of the PLC curve
    real :: f_supply ! fraction of stem water available for leaves per hour
    real :: f_plc    ! fraction of WTC loss due to low water potential (per day)

    ! Allometry
    real :: alphaHT, thetaHT ! height = alphaHT * DBH ** thetaHT
    real :: alphaCA, thetaCA ! crown area = alphaCA * DBH ** thetaCA
    real :: alphaBM, thetaBM ! biomass = alphaBM * DBH ** thetaBM
    real :: phiRL            ! ratio of fine root to leaf area
    real :: phiCSA           ! ratio of sapwood CSA to target leaf area
    real :: tauNSC           ! residence time of C in NSC (to define storage capacity)
    real :: fNSNmax          ! multiplier for NSNmax
    real :: f_N_add
    real :: transT           ! Structural transitional time for canopy layer trees
    ! Default C/N ratios
    real :: CNleaf0
    real :: CNroot0
    real :: CNwood0
    real :: CNseed0
    ! phenology
    real :: Tc0_OFF    ! C, for turning OFF a growth season
    real :: Tc0_ON     ! C, for turning ON a growth season
    real :: gdd_crit       ! C, critical value of GDD5 for turning ON growth season
    real :: gdd_par1
    real :: gdd_par2
    real :: gdd_par3
    real :: betaON       ! Critical soil moisture for phenology ON
    real :: betaOFF      ! Critical soil moisture for phenology OFF
    real :: AWD_crit     ! Critical plant water deficit, transp/W_demand
    !  vital rates
    real :: AgeRepro     ! the age that can reproduce
    real :: v_seed       ! fracton of G_SF to G_F
    real :: s0_plant     ! size of the seedlings, kgC/indiv
    real :: prob_g       ! seed germination probability
    real :: prob_e       ! seed establishment probability
    real :: mu0_topL     ! yearly mortality rate of canopy layer trees
    real :: D0mu         ! Reference diameter for size-dependent mortality
    real :: A_un         ! Parameter for understory mortality affected by layers
    real :: A_sd         ! Max multiplier for seedling mortality
    real :: B_sd         ! Mortality sensitivity for seedlings
    real :: A_DBH        ! Max mulitplier for DBH-based mortality
    real :: B_DBH        ! Sensitivity to dbh
    real :: s_hu         ! hydraulic mortality sensitivity
    real :: W_mu0        ! Half-mortality transp deficit ratio, 0.5, 0.75, 2.5
    ! Population level variables
    real :: f_cGap    ! fraction of internal gaps in the canopy
    real :: LAImax    ! max. LAI
    real :: LAImax_u  ! max. LAI understorey
    real :: LAI_light ! light controlled maximum LAI
    integer :: n_cc   ! for calculating LAImax via cc%LAImax derived from cc%NSN
  end type

  !----------cohort-----------------
  type :: cohort_type
    ! for climate-vegetation type
    integer :: phenotype    ! phenology type: 0 for deciduous, 1 for evergreen
    logical :: firstday     ! First day of a growing season
    integer :: pt           ! photosynthetic physiology of species

    ! ---- Cohort properties
    integer :: ccID   = 0   ! cohort ID
    integer :: species= 0   ! PFT type
    integer :: layer  = 1   ! the layer of this cohort (numbered from top, top layer=1)
    real :: nindivs   = 1.0 ! plant density, individuals/m2
    real :: layerfrac = 0.0 ! fraction of layer area occupied by this cohort
    real :: topyear   = 0.0 ! the years that a plant in top layer

    ! Phenology-related
    integer :: status = 0    ! growth status of plant: 1 for ON, 0 for OFF
    integer :: Ngd    = 0    ! growing days
    integer :: Ndm    = 0    ! dormant days
    integer :: Ncd    = 0    ! number of cold days in non-growing season
    real    :: GDD    = 0.0  ! for phenology
    real    :: ALT    = 0.0  ! growing season accumulative cold temperature
    real    :: AWD    = 1.0  ! Accumulative water deficit (Demand - Transp)/Demand
    real    :: Tc_OFF = -99. ! Critical Tc for PhenoOFF
    real    :: GDD_ON = -99. ! Threshold of GDD for PhenoON

    ! Individual properties
    real :: mu        = 0.02! Cohort mortality rate
    real :: age       = 0.0 ! age of cohort, years
    real :: leafage   = 0.0 ! leaf age (year)
    real :: dbh       = 0.0 ! diameter at breast height, m
    real :: height    = 0.0 ! vegetation height, m
    real :: Acrown    = 1.0 ! crown area, m2/individual
    real :: Aleaf     = 0.0 ! total area of leaves, m2/individual
    real :: lai       = 0.0 ! crown leaf area index, m2/m2
    real :: D_bark    = 0.0 ! thickness of bark
    real :: bl_max    = 0.0 ! Max. leaf biomass, kg C/individual
    real :: br_max    = 0.0 ! Max. fine root biomass, kg C/individual
    real :: CSAsw     = 0.0
    real :: DBH_ys        ! DBH at the begining of a year (growing season)

    ! Photosynthesis
    real :: An_op     = 0.0 ! mol C/(m2 of leaf per year)
    real :: An_cl     = 0.0 ! mol C/(m2 of leaf per year)
    real :: w_scale   = -9999.0
    real :: extinct   = 0.75     ! light extinction coefficient in the canopy for photosynthesis

    ! Carbon fluxes
    real :: gpp  = 0.0 ! gross primary productivity kg C/step
    real :: npp  = 0.0 ! net primary productivity kg C/step
    real :: resp = 0.0 ! plant respiration
    real :: resl = 0.0 ! leaf respiration
    real :: resr = 0.0 ! root respiration
    real :: resn = 0.0 ! N fixation respiration
    real :: resg = 0.0 ! growth respiration
    real :: NPPleaf = 0.0 ! C allocated to leaf, root, and wood
    real :: NPProot = 0.0 !
    real :: NPPwood = 0.0 !

    ! Carbon pools
    real :: bl     = 0.0 ! biomass of leaves, kg C/individual
    real :: br     = 0.0 ! biomass of fine roots, kg C/individual
    real :: bsw    = 0.0 ! biomass of sapwood, kg C/individual
    real :: bHW    = 0.0 ! biomass of heartwood, kg C/individual
    real :: seedC  = 0.0 ! biomass put aside for future progeny, kg C/individual
    real :: nsc    = 0.0 ! non-structural carbon, kg C/individual

    ! ---- Nitrogen model related parameters
    real :: NSNmax = 0.
    real :: NSN    = 0. ! non-structural N pool
    real :: leafN  = 0.
    real :: swN    = 0.
    real :: hwN    = 0. ! N of heart wood
    real :: rootN  = 0. ! N of fine roots
    real :: seedN  = 0. !
    real :: fixedN = 0. ! fixed N at each stem per tree
    real :: extraC = 0. ! The carbon amount that is re-allocated to woody tissues due to N deficit
    real :: N_uptake = 0.

    ! ---- water uptake-related variables
    real :: root_length(soil_L) ! m
    real :: rootarea ! total fine root area per tree
    real :: ArootL(soil_L) = 0.0 ! Root area per layer
    real :: WupL(soil_L) = 0.0 ! normalized vertical distribution of uptake
    real :: Q_soil(soil_L) = 0.0 ! Soil to roots water flux (kg H2O/tree/step)
    real :: W_supply  ! potential water uptake rate per unit time per tree
    real :: totDemand  = 0.0 ! Total water demand in a growing season
    real :: transp   ! transpiration rate per tree per time step
    real :: uptake_frac(soil_L) ! for LM3 soil water uptake, Weng, 2017-10-28
    real :: K_r,r_r
    real :: root_zeta

    ! Hydraulics
    integer :: Nrings = 1
    real :: psi_s0   ! Equilibrium stem base water potential (soil-stem flux=0)
    real :: psi_leaf ! MPa, leaf water potential
    real :: psi_stem ! MPa, stem water potential
    real :: H_leaf ! Leaf capacitance, kgH2O MPa-1 (per tree)
    real :: H_stem ! Stem capacitance, kgH2O MPa-1 (per tree)
    real :: W_leaf ! Leaf water content, kgH2O (per tree)
    real :: W_stem ! Stem water content, kgH2O (per tree)
    real :: W_dead ! water storage in heartwood, just for balance counting.
    real :: Wmax_l ! Leaf max water content, kgH2O (per tree)
    real :: Wmax_s ! Stem max water content, kgH2O (per tree)
    real :: Wmin_l ! Leaf min water content, kgH2O (per tree)
    real :: Wmin_s ! Stem min water content, kgH2O (per tree)
    real :: Q_stem ! water flux from soil to stems (kg/tree/step)
    real :: Q_leaf ! water flux from stems to leaves (kg/tree/step)
    real :: Ktrunk ! trunk water conductance, m/(s MPa)
    real :: Asap   ! Functional cross sectional area
    real :: Atrunk ! Sum of all rings
    real :: treeHU ! total water transported by the functional sapwood, m^3
    real :: treeW0 ! total WTC0 of the sapwood, m^3
    real :: Kx(Ysw_max)    = 0.0 ! Initial conductivity of the woody generated in each year
    real :: WTC0(Ysw_max)  = 0.0 ! lifetime water transfer capacity
    real :: accH(Ysw_max)  = 0.0 ! m, total water transport for functional conduits
    real :: plcH(Ysw_max)  = 0.0 ! m, WTC cost at the low xylem water potential
    real :: farea(Ysw_max) = 0.0 ! fraction of functional area, 1.0/(exp(r_DF*(1.0-accH[j]/W0[j]))+1.0)
    real :: Rring(Ysw_max) = 0.0 ! Radius to the outer edge
    real :: Lring(Ysw_max) = 0.0 ! Length of xylem conduits
    real :: Aring(Ysw_max) = 0.0 ! Area of each ring
    real :: Kring(Ysw_max) = 0.0 ! Conductance of each ring

    ! Diagnostics
    real :: Aleafmax = 0.0  ! Yearly maximum leaf area
    real :: dailyTrsp ! Daily transpiration
    real :: dailyWdmd ! Plant water demand
    real :: dailyGPP   ! kgC/tree day-1
    real :: dailyNPP
    real :: dailyResp
    real :: dailyNup
    real :: annualTrsp
    real :: annualGPP ! C flux/tree
    real :: annualNPP
    real :: annualResp
    real :: CO2_c ! ppm
    real :: NfixDaily= 0.0
    real :: NfixedYr = 0.0 ! annual N fixation per unit crown area
    real :: NupYr= 0.0

  end type cohort_type

  !---------------------------
  type :: vegn_tile_type
    integer :: tileID
    integer :: landuse = LU_NTRL
    integer :: n_cohorts = 0
    integer :: n_initialCC = 0
    integer :: n_years   = 0
    integer :: n_canopycc = 0
    type(cohort_type), pointer :: cohorts(:)=>NULL()
    type(cohort_type), pointer :: initialCC(:)=>NULL()
    type(vegn_tile_type), pointer :: prev => null() ! Pointer to the older vegn tile
    type(vegn_tile_type), pointer :: next => null() ! Pointer to the younger vegn tile
    real :: area              ! m2
    real :: age = 0.0         ! tile age
    real :: LAI               ! leaf area index
    real :: LAImax            ! growing season max
    real :: CAI               ! crown area index
    real :: LAI_L(CLmax) = 0.0 ! LAI of each crown layer, max. 9
    real :: f_gap(CLmax) = 0.0 ! gap fraction of each crown layer
    real :: CAI_L(CLmax) = 1.0 ! crown overlap of each crown layer (squeeze factor)
    real :: kp(CLmax)    = 0.0 ! light extinction coefficient for each layer
    ! uptake-related variables
    real :: root_distance(soil_L) ! characteristic half-distance between fine roots, m
    real :: ArootL(soil_L) = 0.0 ! Root are per layer
    ! averaged quantities for PPA phenology
    real :: tc_daily = 0.0
    real :: tc_pheno = 0.0 ! smoothed canopy air temperature for phenology

    ! litter and soil carbon pools
    real :: litter = 0.0 ! litter flux
    real :: SOC(5) = 0. ! metabolicL, structuralL, microbial, fastSOM, slowSOM
    real :: SON(5) = 0.

    !!  Nitrogen pools, Weng 2014-08-08
    real :: mineralN= 0.  ! Mineral nitrogen pool, (kg N/m2)
    real :: totN    = 0.
    real :: N_uptake= 0.0  ! kg N m-2 hour-1
    real :: fixedN  = 0.0  ! kg N/step
    real :: Nm_Soil = 0.0  ! annual available N in a year
    real :: Nm_Fire = 0.0  ! Mineralized N due to burning
    real :: NorgP2S = 0.0  ! annual N from plants to soil
    real :: previousN      ! an weighted annual available N
    real :: initialN0

    ! Soil water
    integer :: soiltype = 3
    real :: FLDCAP  ! soil field capacity
    real :: WILTPT  ! soil wilting point (0.xx)
    real :: evap    ! kg m-2 per unit fast time step (mm/hour)
    real :: transp  ! kg m-2 hour-1
    real :: runoff  ! Water runoff of the veg tile, unit?
    real :: thetaS  ! moisture index (ws - wiltpt)/(fldcap - wiltpt)
    real :: wcl(soil_L)       ! volumetric soil water content for each layer
    real :: freewater(soil_L) ! Available water in each layer
    real :: psi_soil(soil_L)  ! MPa
    real :: K_soil(soil_L)    ! Kg H2O/(m2 s MPa)
    real :: soilZ             ! Soil depth
    real :: soilWater         ! kg m-2 in root zone
    real :: soilWP0           ! kg m-2, minimum soil water (at WILTPT)

    ! Vegetation water content
    real :: W_leaf  ! Leaves
    real :: W_stem  !
    real :: W_dead  ! Heartwood ?

    ! water uptake-related variables
    real :: RAI                ! root area index
    real :: RAIL(soil_L) = 0.0 ! Root length per layer, m of root/m
    real :: W_uptake           ! water uptake rate per unit time per m2

    !  Carbon fluxes
    real :: gpp  = 0  ! gross primary production, kgC m-2 yr-1
    real :: npp  = 0  ! net primary productivity
    real :: resp = 0 ! auto-respiration of plants
    real :: rh   = 0  ! soil carbon lost to the atmosphere

    ! Methane fluxes (soil biogeochemistry)
    real :: ch4_prod = 0.0  ! CH4 production (as C), kgC m-2 step-1
    real :: ch4_oxid = 0.0  ! CH4 oxidation (as C),  kgC m-2 step-1
    real :: ch4_emit = 0.0  ! CH4 emission (as C),   kgC m-2 step-1

    !  fire disturbance
    real :: C_burned  = 0.0 ! Carbon released to atmosphere via fire
    real :: TreeCA    = 0.0 ! tree CAI in the top layer, for fire spread
    real :: GrassCA   = 0.0 ! grass CAI, for the initial fire severity
    real :: GrassBM   = 0.0 ! Grass biomass at the end of growing season, for fire severity
    real :: annualPET = 0.0 ! Potential ET, yearly, for fire risk calculation
    real :: YearlyTmp = 0.0 ! Yearly mean air temperature, Celcius degree
    real :: Frisk     = 0.0 ! Probability of climatic fire risk
    real :: Pfire     = 0.0 ! Probability of burning

    ! Daily diagnostics
    real :: dailyGPP
    real :: dailyNPP
    real :: dailyResp
    real :: dailyRh
    real :: dailyCH4 = 0.0
    real :: dailyNup
    real :: NfixDaily
    real :: dailyLFLIT = 0.0   !kgC day-1, leaf litter flux

    ! Annual diagnostics
    real :: CO2_c ! ppm, annual atmospheric CO2 concentration
    real :: dailyPrcp = 0.0, annualPrcp = 0.0 ! mm m-2 yr-1
    real :: dailyTrsp = 0.0, dailyEvap  = 0.0, dailyRoff = 0.0 ! mm m-2 yr-1
    real :: dNorg_daily = 0.0, dNgas_daily = 0.0, dNmin_daily = 0.0
    real :: dNorg_Yr    = 0.0, dNgas_Yr    = 0.0, dNmin_Yr    = 0.0
    real :: annualTrsp= 0.0, annualEvap = 0.0, annualRoff= 0.0 ! mm m-2 yr-1
    real :: annualGPP = 0.0 ! kgC m-2 ground yr-1
    real :: annualNPP = 0.0
    real :: annualResp= 0.0
    real :: annualRh  = 0.0
    real :: annualCH4 = 0.0
    real :: NupYr     = 0.0 ! accumulated N uptake kgN m-2 yr-1
    real :: NfixedYr  = 0.0 ! fixe N in a tile
    ! for annual reporting at tile level
    real :: NSC, SeedC, leafC, rootC, SwC, HwC
    real :: NSN, SeedN, leafN, rootN, SwN, HwN
    real :: totSeedC,totSeedN
    ! for cohort plant types (climate-vegetation relationship, Biome, LM3)
    real :: t_ann  = 0.0 ! annual mean T, degK
    real :: t_cold = 0.0 ! average temperature of the coldest month, degK
    real :: p_ann  = 0.0 ! annual mean precip
    real :: ncm    = 0.0 ! number of cold months
  end type vegn_tile_type

  type :: land_grid_type
    type(vegn_tile_type), pointer :: firstVegn => NULL() ! first vegn tile
    integer :: nTiles = 0 ! The total number of tiles
    real    :: area       ! are of this grid, m2
    real    :: Tc_daily   ! Daily mean temperature
    real    :: Tc_pheno   ! Running mean for phenology
  end type land_grid_type

  !----------------------------------------
  type :: soil_pars_type
    real :: GMD ! geometric mean partice diameter, mm
    real :: GSD ! geometric standard deviation of particle size
    real :: vwc_wilt
    real :: vwc_fc
    real :: vwc_sat
    real :: vlc_min
    real :: k_sat_ref ! hydraulic conductivity of saturated soil, kg/(m2 s)
    real :: psi_sat_ref ! saturation soil water potential, m
    real :: chb         ! Soil texture parameter
    real :: alpha       ! vertical changes of soil property, 1: no change
    real :: heat_capacity_dry ! J m-3 k-1
    real :: tfreeze
  end type soil_pars_type

  type :: soil_prog_type
    real wl
    real ws
    real T
  end type soil_prog_type

  ! Soil tile data types
  type :: soil_tile_type
    integer :: tag ! kind of the soil
    type(soil_pars_type) :: pars
    type(soil_prog_type), pointer :: prog(:)
    real,                 pointer :: w_fc(:)
    real,                 pointer :: w_wilt(:)
    !real :: Eg_part_ref
    !real :: z0_scalar
    ! data that were local to soil.f90
    real, pointer :: heat_capacity_dry(:)
    real, pointer :: e(:),f(:)
    ! added to avoid recalculation of soil hydraulics in case of Darcy uptake
    real          :: uptake_T ! update temperature from previous time step
    real, pointer :: psi(:) ! soil water potential
  end type soil_tile_type

  !Forcing data type
  type :: climate_data_type
    integer :: year       ! Year
    integer :: doy        ! day of the year
    real :: hod           ! hour of the day
    real :: PAR           ! umol m-2 s-1
    real :: radiation     ! W/m2
    real :: Tair          ! air temperature,  K
    real :: Tsoil         ! soil temperature, K
    real :: RH            ! relative humidity, 0~1
    real :: rain          ! kgH2O m-2 s-1
    real :: windU         ! wind velocity (m s-1)
    real :: P_air         ! pa
    real :: CO2           ! mol/mol
    real :: eCO2          ! mol/mol
    real :: soilwater     ! soil moisture, vol/vol
    real :: N_input       ! N deposition (kgN/m2/yr)
  end type climate_data_type

  !-------------Vars for the model -------------
  type(land_grid_type), pointer :: land ! Land grid
  ! Input forcing data
  type(climate_data_type), pointer :: forcingData(:)
  ! output files
  integer :: fno1=211, fno2=212, fno3=213,fno4=214, fno5=215, fno6=216

  ! -------------------------------------------
  ! Soil water parameters
  real :: WaterLeakRate = 0.0 ! Soil water leak rate, fraction per day
  real :: thksl(soil_L)=(/0.05,0.25,0.5,1.0,1.2/) ! m, thickness of soil layers

  ! Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
  real :: GMD(n_dim_soil_types) = & ! geometric mean partice diameter, mm
  (/ 0.7, 0.4, 0.3, 0.1, 0.1, 0.07, 0.007, 0.3, 0.3 /)
  real :: GSD(n_dim_soil_types) = & ! geometric standard deviation of particle size
  (/5.0, 5.3, 7.4, 6.1, 6.1, 14.0, 15.0, 7.4, 7.4 /)
  real :: vwc_sat(n_dim_soil_types)= &
  (/ 0.380, 0.445, 0.448, 0.412, 0.414, 0.446, 0.424, 0.445, 0.445   /)
  !real :: vlc_min(n_dim_soil_types)
  real :: k_sat_ref(n_dim_soil_types)= & ! mol/(s MPa m) , hydraulic conductivity of saturated soil,
  (/ 70.8, 75.1, 28.2, 12.1, 11.1, 12.7, 1.69, 28.2, 28.2 /)
  !(/ 130.8, 75.1, 53.2, 12.1, 11.1, 12.7, 1.69, 53.2, 53.2 /)
  real :: psi_sat_ref(n_dim_soil_types) = & ! Pa
  (/ -600., -790., -910., -1580., -1680., -1880., -5980., -790., -790./)
  real :: chb(n_dim_soil_types) = &         ! Soil texture parameter
  (/   3.5,   6.4,  11.0,   4.8,   6.3,   8.4,   6.3,   6.4,   6.4   /)
  real :: alphaSoil(n_dim_soil_types) = 1.0       ! *** REPLACE LATER BY alpha(layer)
  real :: heat_capacity_dry(n_dim_soil_types) = &
  (/ 1.2e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.4e6,   1.0   /)

  ! -------------------------------------------
  ! Vegetation and soil BGC parameters
  ! Cohort management
  real :: diff_S0 = 0.2 ! percentage of the difference between cohorts for merging

  ! Growth parameters:
  integer :: MaxGrassLyr = 2  ! Maximal layers that grasses can survive
  real :: MaxGrassAge = 3.0   ! Maximum grass  age (years)
  real :: MaxGrassCA = 1.2    ! Maximum grass CA in the top layer
  real :: fsc_fine  = 1.0     ! fraction of fast turnover carbon in fine biomass
  real :: fsc_wood  = 0.0     ! fraction of fast turnover carbon in wood biomass
  real :: GR_factor = 0.33    ! growth respiration factor
  real :: l_fract   = 0.0     ! 0.25  ! 0.5 ! fraction of the carbon retained after leaf drop
  real :: retransN  = 0.0     ! retranslocation coefficient of Nitrogen
  real :: f_iniBSW  = 0.1     ! 0.01, 0.2
  real :: f_N_add   = 0.02    ! re-fill of N for sapwood
  real :: f_LFR_max = 0.85    ! max allocation to leaves and fine roots each step
  real :: c_LLS     = 28.5714 ! yr/(kg C m-2), 1/LMAs, ! Leaf life span: leafLS = c_LLS * LMA, (LMAs = 0.035., leafLS = 1.0)
  real :: rho_N_up0 = 0.1     ! 0.05 ! hourly N uptake rate, fraction of the total mineral N
  real :: N_roots0  = 0.4     ! root biomass at half max N-uptake rate,kg C m-2
  real :: C0_Nfix   = 12.0    ! gC/gN, carbon cost of N fixation, FUN model, Fisher et al. 2010, GBC; Kim

  ! Plant hydraulics
  real :: psi0_osm = 0.5     ! MPa, leaf osmotic pressure
  real :: TK0_leaf = 0.003   ! leaf thickness at reference LMA (1.0 LMA), m
  real :: kx0_WD   = 5.0     ! kx0 at reference WD
  real :: WTC0_WD  = 1200.0  ! Water transporat capacity at reference WD
  real :: CR0_LF   = 0.50    ! Compression ratio of leaf tissues at reference LMA
  real :: CR0_WD   = 0.35    ! Compression ratio of wood tissues at reference WD
  real :: p50_WD   = -1.565  ! stem psi50 at reference WD
  real :: r_DF     = 100.0   ! sensitivity of defunction due to water transport usage
  real :: m0_WTC   = 8.0     !  DBH-WTC0 Radial variations, 12000/300 = 40,
  real :: m0_kx    = 8.0     ! DBH-Kx0 Radial variations
  real :: ths0_WD  = 3.0     ! exponential of the PLC function for (psi/psi50)
  real :: fsup0_WD = 1.0     ! Fraction of plant water for transpiration per hour at zero WD
  real :: fplc0_WD = 1.0     ! Fraction of WTC loss at low PLC at zero WD
  real :: A_plc0_WD= -3.0    ! Parameter in f_plc = fplc0_WD * exp(A*R_WD)
  real :: plc_crit = 0.5     ! Critical value of plc for making a damage to xylems

  ! Phenology parameters
  ! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)
  integer :: N0_GD     = 90    ! base growing days, 90 days, with a -5 substraction of Tc_crit
  integer :: Days_thld = 60    ! minimum days of the growing or non-growing season
  real    :: cold_thld = -20.  ! threshold of accumulative low temperature (sum(dT*day)
  real    :: T0_gdd    = 5.0   ! Celsius degree
  real    :: T0_chill  = 10.0  ! Celsius degree

  ! Fire regimes
  real :: EnvF0        = 0.0   ! Fixed environmental fire risk, 11/25/2025
  real :: MI0Fire      = 0.5   ! Frisk = 1.0/(1.0 + exp(A_MI*(P/PET - MI0Fire)))
  real :: s0_max       = 1.0   ! maximum fire severity for canopy fires (0.0~1.0)
  real :: A_MI         = 20.0  ! shape parameter of Fire risk vs. P/PET curve
  real :: FSBM0        = 0.2   ! kgC m-2, grass fire severity parameter, as a function of grass BM
  real :: f_bk         = 0.1105! coefficient of bark thickness, Hoffmann et al. 2012.
  ! shrubs: Y=1.105*X^1.083; trees: Y=0.31*X^1.276 for (Y:mm, X:cm)
  real :: r_BK0        = -240.0! bark resistance, exponential equation, 120 --> 0.006 m of bark

  ! Soil organic matter decomposition
  real :: K0SOM(5)     = [0.8, 0.25, 3.0, 1.5, 0.05] ! (/0.8, 0.25, 2.5, 1.0, 0.2/) ! turnover rate of SOM pools (yr-1)
  real :: CUEmax0      = 0.1 ! Maximum carbon use efficiency of microbes
  real :: K_DeNitr     = 8.0     ! mineral Nitrogen turnover rate
  real :: fDON         = 0.02    ! fraction of DON production in decomposition
  real :: rho_SON      = 0.05    ! SON release rate per year
  real :: f_M2SOM      = 0.8     ! the ratio of C and N returned to litters from microbes
  real :: etaN         = 0.025   ! Coefficient of N loss through runoff (etaN*runoff is a fraction of organic or mineral N)
  real :: fdsvN        = 0.30    ! Max fraction of soluble N taken out by runoff (01/17/2026, Weng)

  ! --- Methane (CH4) parameters for soil BGC (MVP) ---
  real :: CH4_alpha = 0.20   ! fraction of Rh routed to CH4 under fully anaerobic conditions (0-1)
  real :: CH4_beta_ox = 0.50 ! fraction of produced CH4 oxidized in oxic fraction (0-1)
  real :: CH4_wfps0 = 0.60   ! wetness threshold for onset of anoxia (WFPS-like, 0-1)
  real :: CH4_wfps1 = 0.95   ! wetness for fully anaerobic (WFPS-like, 0-1)

  ! -------- PFT-specific parameters ----------
  ! Define parameter arrays with the same value. Preset of global PFTs is in Preset_GlobalPFTs
  integer :: pt(0:MSPECIES)       = 0 ! 0 for C3, 1 for C4
  integer :: phenotype(0:MSPECIES)= 0 ! 0 for Deciduous, 1 for evergreen
  integer :: lifeform(0:MSPECIES) = 1 ! life form of PFTs: 0 for grasses, 1 for trees

  ! Allometry and whole tree parameters
  real :: alphaHT(0:MSPECIES) = 36.0
  real :: alphaCA(0:MSPECIES) = 150.0
  real :: alphaBM(0:MSPECIES) = 5200.0
  real :: thetaHT(0:MSPECIES) = 0.5 !
  real :: thetaCA(0:MSPECIES) = 1.5
  real :: thetaBM(0:MSPECIES) = 2.5
  real :: phiRL(0:MSPECIES)   = 3.5 ! ratio of fine root area to leaf area
  real :: phiCSA(0:MSPECIES)  = 0.25E-4 ! ratio of sapwood area to leaf area
  real :: tauNSC(0:MSPECIES)  = 6 ! 3 ! NSC residence time,years
  real :: fNSNmax(0:MSPECIES) = 5 ! 5 ! multiplier for NSNmax as sum of potential bl and br
  real :: transT(0:MSPECIES)  = 3 ! Years
  real :: f_cGap(0:MSPECIES)  = 0.1  ! The gaps between trees
  real :: LFR_rate(0:MSPECIES)= 1.0/21.0

  ! Leaf parameters
  real :: leaf_size(0:MSPECIES)= 0.04 !
  real :: LAImax(0:MSPECIES)   = 3.5    ! maximum LAI for a tree
  real :: LAI_light(0:MSPECIES)= 4.0    ! maximum LAI limited by light
  real :: LMA(0:MSPECIES)      = 0.035  ! leaf mass per unit area, kg C/m2
  real :: leafLS(0:MSPECIES)   = 1.0
  real :: LNbase(0:MSPECIES)   = 1.3E-3 !functional nitrogen per unit leaf area, kg N/m2, 1.1E-3 for Acer, 1.5E-3 for Populus
  real :: CN0leafST(0:MSPECIES)= 40.0 ! 80.0 ! CN ratio of leaf supporting tissues

  ! photosynthesis parameters
  real :: Vmax(0:MSPECIES)     = 35.0E-6 ! mol m-2 s-1
  real :: m_cond(0:MSPECIES)   = 9.0 ! 7.0 !
  real :: alpha_ps(0:MSPECIES) =  0.06 !
  real :: Vannual(0:MSPECIES)  = 1.2 ! kgC m-2 yr-1
  real :: ps_wet(0:MSPECIES)   = 0.3 ! wet leaf photosynthesis down-regulation: 0.3 means
  ! photosynthesis of completely wet leaf will be 30% reduction

  ! Wood parameters
  real :: rho_wood(0:MSPECIES) = 300.0 ! kgC m-3
  real :: f_taper(0:MSPECIES)  = 0.75 ! taper factor, from a cylinder to a tree
  real :: IgniteP(0:MSPECIES)  = 0.02 ! Ignition probability at fire-friendly climates
  real :: mu0fire(0:MSPECIES)  = 0.99 ! Mortality rate at fire (sensitivity to fire)

  ! root parameters
  real :: alpha_FR(0:MSPECIES)  = 1.2 ! Fine root turnover rate yr-1
  real :: rho_FR(0:MSPECIES)    = 200 ! woody density, kgC m-3
  real :: root_r(0:MSPECIES)    = 2.9E-4 !(/1.1e-4, 1.1e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4/)
  real :: root_zeta(0:MSPECIES) = 0.6 ! 0.29 !
  real :: root_perm(0:MSPECIES) = 0.5 ! kg H2O m-2 hour-1, defined by Weng
  real :: Kw_root(0:MSPECIES)   = 6.3E-8 * 1.e3 ! (kg m-2 s−1 MPa−1) ! Ref: 6.3±3.1×10−8 (m s−1 MPa−1)
  ! * (1000000.0/18.0)*1.e-6 ! mol /(s m2 Pa)
  !Ref added by Weng, 2021-11-15
  ! Sutka et al. 2011 Natural Variation of Root Hydraulics in Arabidopsis Grown
  ! in Normal and Salt-Stressed Conditions. Plant Physiol. 155(3): 1264–1276.
  ! doi: 10.1104/pp.110.163113
  ! Miyamoto et al. 2001. Hydraulic conductivity of rice roots. J. Exp. Bot., 52: 1835–1846,
  ! doi: 10.1093/jexbot/52.362.1835

  !(/1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5/)
  ! fine root membrane permeability per unit membrane area, kg/(m3 s).
  ! Root membrane permeability is "high" for the value from Siqueira et al., 2008,
  ! Water Resource Research Vol. 44, W01432, converted to mass units

  !real :: rho_N_up0(0:MSPECIES) = 0.5 ! fraction of mineral N per hour
  !real :: N_roots0(0:MSPECIES) = 0.3 ! kgC m-2

  ! Respiration rates
  real :: gamma_L(0:MSPECIES)  = 0.02 !
  real :: gamma_LN(0:MSPECIES) = 70.5 ! 25.0  ! kgC kgN-1 yr-1
  real :: gamma_SW(0:MSPECIES) = 0.02 ! 0.08 ! kgC m-2 Acambium yr-1
  real :: gamma_FR(0:MSPECIES) = 0.6 ! 12 !kgC kgN-1 yr-1 ! 0.6: kgC kgN-1 yr-1

  ! Phenology parameters
  real :: Tc0_OFF(0:MSPECIES)  = 15. ! Critical temperature (C) for phenology OFF
  real :: Tc0_ON(0:MSPECIES)   = 10. ! Critical temperature (C) for phenology ON
  real :: gdd_crit(0:MSPECIES) = 300. ! Replaced by an equation: sp%gdd_par1 + sp%gdd_par2 * exp(sp%gdd_par3*cc%ncd)
  real :: AWD_crit(0:MSPECIES) = 0.7  ! Critical plant water availability factor (0~1)
  real :: betaON(0:MSPECIES)   = 0.2  ! Critical soil moisture for phenology ON
  real :: betaOFF(0:MSPECIES)  = 0.1  ! Critical soil moisture for phenology OFF
  real :: gdd_par1(0:MSPECIES) = 30.0   !50.d0   ! -68.d0
  real :: gdd_par2(0:MSPECIES) = 800. ! 650.d0  !800.d0  ! 638.d0
  real :: gdd_par3(0:MSPECIES) = -0.02 ! -0.01d0
  ! Reproduction prarameters
  real :: AgeRepro(0:MSPECIES) = 5.0  ! year
  real :: v_seed(0:MSPECIES)   = 0.1  ! fraction of allocation to wood+seeds
  real :: s0_plant(0:MSPECIES) = 0.05 ! kgC, initial seedling size
  real :: prob_g(0:MSPECIES)   = 1.0
  real :: prob_e(0:MSPECIES)   = 1.0

  ! Mortality parameters
  real :: mu0_topL(0:MSPECIES) = 0.012 ! 0.01 ! yearly ! 0.012 for Acer, 0.0274 for Populus
  real :: D0mu(0:MSPECIES)     = 1.2     ! m, Mortality curve parameter
  real :: A_un(0:MSPECIES)     = 3.0     ! Multiplier for understory mortality
  real :: A_sd(0:MSPECIES)     = 9.0     ! Max multiplier for seedling mortality
  real :: B_sd(0:MSPECIES)     = -20.    ! Mortality sensitivity for seedlings
  real :: A_DBH(0:MSPECIES)    = 4.0     ! Max multiplier for DBH-based mortality
  real :: B_DBH(0:MSPECIES)    = 0.125   ! 0.25   ! Size-based Mortality sensitivity, m
  real :: s_hu(0:MSPECIES)     = -25.0   ! hydraulic mortality sensitivity
  real :: W_mu0(0:MSPECIES)    = 1.0     ! Jeremy's half-mortality transp deficit, high:0.5, low: 0.75, No effects: 2.5

  ! Plant hydraulics parameters
  real :: kx0(0:MSPECIES)      = 5.0 ! (mm/s)/(MPa/m) !132000.0 ! 6000.0   ! (m/yr-1)/(MPa/m)
  real :: WTC0(0:MSPECIES)     = 1200.0  ! 2000, m /lifetime
  real :: CR_Leaf(0:MSPECIES)  = 0.5 ! leaf compression ratio per MPa
  real :: CR_Wood(0:MSPECIES)  = 0.2 ! Wood compression ratio per MPa
  real :: psi0_LF(0:MSPECIES)  = -3.0  ! MPa
  real :: psi0_WD(0:MSPECIES)  = -3.0  ! MPa
  real :: psi50_WD(0:MSPECIES) = -1.5  ! MPa !wood potential at which 50% conductivity lost, MPa
  real :: Kexp_WD(0:MSPECIES)  = 3.0
  real :: f_supply(0:MSPECIES) = 0.5
  real :: f_plc(0:MSPECIES)    = 0.05  ! fraction of WTC loss due to low water potential (per day)

  ! C/N ratios for plant pools
  real :: CNleaf0(0:MSPECIES)  = 25.  ! C/N ratios for leaves
  real :: CNwood0(0:MSPECIES)  = 350. ! C/N ratios for woody biomass
  real :: CNroot0(0:MSPECIES)  = 40.0 ! C/N ratios for leaves ! Gordon & Jackson 2000
  real :: CNseed0(0:MSPECIES)  = 20.0 ! C/N ratios for seeds
  real :: R0_Nfix(0:MSPECIES)  = 0.0  ! Reference N fixation rate (0.03 kgN kg rootC-1 yr-1)
  real :: S_facuN(0:MSPECIES)  = 0.0  ! Faculative N fixation intensity, 0~1. 0: Non; 1.0 full extraC usage

  ! Standard cohorts for the ESS PFTs, Weng, 09/12/2025
  !--------------------------------------0:C4G, 1:C3G, 2:TrE, 3:TrD, 4:TmE, 5:TmD, 6:Nfx, 7:DeS
  real, parameter :: std_den(0:N_EST) = [8.0,   8.0,   .05,   .05,   .05,   .05,   .05,   .05]  ! initial individual density, individual/m2
  real, parameter :: std_bsw(0:N_EST) = [.005,  .005,  .02,   .02,   .02,   .02,   .02,   .02]  ! initial biomass of sapwood, kg C/individual
  real, parameter :: std_nsc(0:N_EST) = [.005,  .005,  .02,   .02,   .02,   .02,   .02,   .02]  ! initial non-structural biomass, kg C/individual

  !----- Initial conditions and model control -------------
  integer :: I
  integer, parameter :: M_initialCH = MSPECIES ! Weng, 2014-10-01
  integer :: init_cohort_N = 1 ! M_initialCH
  integer :: init_cohort_sps(M_initialCH) = (/ (I, I = 0, M_initialCH-1) /)
  real :: init_cohort_Indiv(M_initialCH)  = -9.  ! initial individual density, individual/m2
  real :: init_cohort_bl(M_initialCH)     = 0.0  ! initial biomass of leaves, kg C/individual
  real :: init_cohort_br(M_initialCH)     = 0.0  ! initial biomass of fine roots, kg C/individual
  real :: init_cohort_bsw(M_initialCH)    = .01  ! initial biomass of sapwood, kg C/individual
  real :: init_cohort_bHW(M_initialCH)    = 0.0  ! initial biomass of heartwood, kg C/tree
  real :: init_cohort_seedC(M_initialCH)  = 0.0  ! initial biomass of seeds, kg C/individual
  real :: init_cohort_nsc(M_initialCH)    = .01  ! initial non-structural biomass, kg C/individual

  ! Initial soil type, carbon and nitrogen at a vegn tile, Weng 2012-10-24
  integer :: soiltype    = SandyLoam  ! lookup table for soil hydrologic parameters
  real :: init_fast_SOC  = 0.5  ! initial fast soil C, kg C/m2
  real :: init_slow_SOC  = 2.0  ! initial slow soil C, kg C/m2
  real :: init_mineralN  = 0.005  ! Mineral nitrogen pool, (kg N/m2)
  real :: N_input        = 0.002 ! annual N input to soil N pool, kgN m-2 yr-1

  ! Climate-vegetation initialization, 09/20/2025
  real :: Pr_thld = 300.0  ! Desert shrub vs trees, not used!
  real :: MI0DeSB = 0.25   ! Desert shrub vs trees, P/PET
  real :: MI0C3C4 = 0.50   ! Moisture threshold for C4 vs. C3 grasses
  real :: TcrTREE = 12.0   ! Tropical trees vs Temperate/boreal trees
  real :: TcrC3C4 = 0.0    ! Temperature threshold for C4 vs C3 grasses
  real :: meanPrcp, meanPET, meanTmin, Mst_IDX ! Climate statistics for PFT envelopes

  ! Input files
  character(len=80)  :: filepath_in = './input/'
  character(len=80)  :: filepath_out = './output/'
  character(len=80)  :: runID = 'test'
  character(len=160) :: climfile = 'ORNL_forcing.txt'
  character(len=160) :: Scefile = 'ORNL_forcing.txt'
  character(len=160) :: PaleoPfile = 'RMA_P.csv' ! for DroughtPaleo
  character(len=160) :: PaleoTfile = 'RMA_T.csv'
  integer  :: iDraw = 1 ! Sample number
  integer  :: N_VegTile = 1 ! Initial vegn tiles
  integer  :: datalines ! the total lines in forcing data file
  integer  :: StartLine = 1 ! the first step model run start with, for UFL only
  integer  :: yr_data   ! Years of the forcing data
  integer  :: days_data ! days of the forcing data
  real     :: siteLAT = 36.01 !site latitude, ORNL

  ! For global/regional forcing data, soil conditions, and initial conditions
  ! CRU NetCDF file dimensions
  integer, parameter :: N_PFTs  = 8  ! For ESS Biome
  integer, parameter :: N_Crop  = 6
  integer, parameter :: N_Vegs  = 10 ! pft_2011 vegetation types
  integer, parameter :: FM_Yrs  = 1176 ! Total years of Hurtt's crop data
  integer, parameter :: NDIMS = 3, Nlon = 720, Nlat = 360, Ntime = 1460 ! NC file dimensions
  integer, parameter :: Hours_NCstep = hours_per_year/Ntime ! 6
  real,    parameter :: Lon0 = -180.0, Lat0 = -90.0 ! The first grid's coordinates
  real,    parameter :: Wlat = 360.0/Nlon, Wlon = 180.0/Nlat

  ! Data interpolation (hourly)
  integer, parameter :: TargetSteps_per_hour = 1
  integer, parameter :: SHshift = 182*24*TargetSteps_per_hour ! Southern Hemisphere data shift (to July 2nd)
  logical, parameter :: ShiftSHdata = .True.

  type :: grid_initial_type
    integer :: iLon ! grid number along Longitude (from -180 to 180)
    integer :: iLat ! grid number along Latitude (from -90 to 90)
    real    :: SOM(5) ! Soil organic matter (kgC m-2)
    real    :: SON(5) ! Soil organic nitrogen
    real    :: mineralN ! Soil mineral N
    real    :: soiltexture(3)
    real    :: WLTPT, FLDCP ! soil wilting point and field capacity (0.xx)
    real    :: N_input      ! kgN m-2 yr-1
    real, pointer :: Farm(:)
    real, pointer :: VegCover(:)        ! for pft2011_0.5x0.5.nc only
    real, pointer :: climate(:,:)       ! Ntimes, N_vars
  end type grid_initial_type

  ! Model mechanisms setting
  logical  :: Do_DroughtMu        = .True.  ! Drought-induced mortality, Lichstein (2024)
  logical  :: MergeLowDenCohorts  = .True.  ! Merge low density cohorts to their closest ones
  logical  :: Do_RecoverSP        = .False. ! Species recovery from initial conditions
  logical  :: Do_Fire             = .False. ! Allow fire disturbances
  logical  :: Do_FixedFrisk       = .False. ! Fixed fire risk (regardless of climatic condtions)
  logical  :: Do_FixedFireS       = .True.  ! Fixed fire severity for woody plants in canopy fires
  logical  :: Do_ClosedN_run      = .True.  ! Nitrogen input and output are zero
  logical  :: Do_VariedKx         = .True.  ! trunk new xylem has the same kx or not
  logical  :: Do_VariedWTC0       = .True.  ! WTC0 changes with trunk size
  logical  :: Do_mu0_F_WDen       = .False. ! mu0 as a function of wood density
  logical  :: Do_CH4              = .False. ! Methane emission modeling

  ! For global/regional run, Weng, 2025-07-22
  character (len = 256) :: ncfilepath   = '/media/eweng/HD2/weng/Data/CRU/TRENDY2023/1HX1/'
  character (len = 256) :: veg_path     = '/Users/eweng/Documents/Data/Vegetation/'
  character (len = 256) :: Ndp_path     = '/Users/eweng/Documents/BiomeESS/forcingData/N_deposition/'
  character (len = 20)  :: ncversion    = 'crujra.v2.4.5d.'
  character (len = 20)  :: veg_file     = 'pft2011_0.5x0.5.nc' ! Vegetation coverage

  character (len = 256) :: int_fpath    = '/media/eweng/HD2/weng/Data/CRU/TRENDY2023/1HX1/interpolated/'
  character (len = 80)  :: int_prefix   = 'crujra.v2.4.5d.'
  character (len = 50)  :: GridListFile = 'GlobalVegGridList.csv' ! in int_fpath
  character (len = 5)   :: ncfields(7)  = [character(len=5):: 'tmp','pre','dswrf','spfh','pres','ugrd','vgrd']
  character (len = 6)   :: GridIDFMT    = '(I6.6)' ! For the file name string (GridID)
  character(len=30)     :: Vegstr       = 'TOTAL_VEG'
  character(len=9)      :: VegID(N_Vegs)= [character(len=9) :: &
                           'SHRUBS_BD','SHRUBS_BE','SHRUBS_ND','SHRUBS_NE', &
                           'TREES_BD ','TREES_BE ','TREES_ND ','TREES_NE ', &
                           'GRASS_MAN','GRASS_NAT'] ! pft2011_0.5x0.5.nc
  character(len=5)   :: CropID(N_Crop) = [character(len=5) :: &
                           'c3ann', 'c3per', 'c3nfx',  'c4ann', 'c4per', 'pastr']
  character(len=6)   :: NdpID(4) = [character(len=6) :: &
                           'wetnoy','wetnhx','drynoy','drynhx']
  integer :: LowerLon   = 1
  integer :: UpperLon   = 720 ! Grid number from -179.75 (latitude)
  integer :: LowerLat   = 61
  integer :: UpperLat   = 320 ! Grid number from -89.75 (longitude)
  integer :: yr_start   = 2010
  integer :: yr_end     = 2019
  integer :: grid_No1   = 1     ! the first grid in the grid list file
  integer :: grid_No2   = 56395 ! the last grid in vegetated land
  integer :: N_VegGrids = 1 ! Minimum
  integer :: StepLatLon = 1 ! Skip grids. 1: all; 2: one per 2x2 grids
  integer :: GridID     = 999999 ! 216264                ! = iLon*1000 + iLat
  integer :: HemiSP     = 1 ! 1: North hemisphere; 0: South hemisphere
  integer :: LC_year0   = 1176 ! initialization with the last year's cropland cover (WIEMIP setting)
  real    :: GridVC(N_Vegs) = 0.0 ! For WIE-MIP, vegetation cover percentage
  real    :: GridFR(FM_Yrs) = 0.0 ! Farm land ratio, 1176 years

  type(grid_initial_type), pointer :: LandGrid(:) => null()
  integer,allocatable     :: GridLonLat(:)     ! LonLat
  real,allocatable,target :: GridVegCov(:,:)   ! Grid vegetation composition
  real,allocatable,target :: GridFarm(:,:)     ! (FM_Yrs,N_VegGrids)
  real,allocatable,target :: ClimData(:,:,:)   ! N_yr*Ntime, N_vars, N_VegGrids
  real,allocatable        :: CRUtime(:)        ! Days since 1901-01-01 in CRU data
  !real,    allocatable,target :: CRUData(:,:,:,:)  ! N_yr*Ntime, N_vars, Nlon, Nlat

  ! Output interpolated grid climate data files
  logical :: WriteForcing = .False. ! .True. ! Write interpolated forcing data
  integer :: Grids_UN1   = 99       ! Vegetated grids list file
  integer :: Grids_UN2   = 98       ! Grid cropland cover file

  ! Model run control
  character(len=256) :: file_out(6) ! Output file names
  integer  :: model_run_years = 100
  integer  :: totyears, totdays
  integer  :: steps_per_day  = 24 ! 24 or 48
  integer  :: yr_ResetVeg    = 0 ! reseting vegetation to the initial, clearcut
  integer  :: yr_Baseline    = 1000 ! for DroughtMIP baseline model run years
  integer  :: FreqY0         = 5 ! Interval (yrs) of recovering initial species
  integer  :: equi_days      = 0 ! 100 * 365
  integer  :: steps_per_hour = 1
  real     :: step_hour      = 1.0  ! hour, Time step of forcing data, usually hourly (1.0)
  real     :: step_seconds   = 3600.0
  real     :: dt_fast_yr     = 1.0 / (365.0 * 24.0) ! Hourly
  real     :: dt_daily_yr    = 1.0/365.0 ! Daily

  ! Model output
  logical  :: outputhourly = .True.
  logical  :: outputdaily  = .True.

  ! Scenarios
  character(len=4)  :: CO2Tag = 'aCO2' ! only takes 'aCO2' or 'eCO2', for FACE-MDS-3
  logical  :: Sc_CO2  = .True. ! Use CO2_C if true
  real     :: Sc_prcp = 1.0 ! Scenario of rainfall changes
  real     :: Sc_dT   = 0.0 ! Scenario of temperature changes
  real     :: CO2_c   = 375.0 ! 412 ! PPM, CO2 concentration at 2020

  !-------------Plant and soil parameter types -----------------------
  type(spec_data_type), save :: spdata(0:MSPECIES)         ! PFT-specific parameters
  type(soil_pars_type), save :: soilpars(n_dim_soil_types) ! Soil hydraulics parameters

  ! ------------- Global run model setting name list ------------
  namelist /global_setting_nml/ ncfilepath, Ndp_path, veg_path, &
  int_fpath, ncversion, int_prefix, veg_file, GridListFile,     &
  grid_No1, grid_No2, yr_start, yr_end, LowerLon, UpperLon,     &
  LowerLat, UpperLat, StepLatLon, LC_year0, WriteForcing

  ! ------------- Model initialization name list ------------
  namelist /initial_state_nml/ N_VegTile,                      &
  ! initial vegetation and soil states
  init_cohort_N, init_cohort_sps, init_cohort_Indiv,           &
  init_cohort_bl, init_cohort_br, init_cohort_bsw,             &
  init_cohort_bHW, init_cohort_seedC, init_cohort_nsc,         &
  init_fast_SOC, init_slow_SOC, init_mineralN, N_input,        &
  ! Climate envelopes for initializing PFTs
  Pr_thld, MI0DeSB, MI0C3C4, TcrTREE, TcrC3C4,                 &
  ! Model run controls
  filepath_in,filepath_out,climfile,model_run_years,runID,     &
  outputhourly,outputdaily,Sc_prcp,Sc_dT,Sc_CO2,CO2_c, &
  ! Model components
  MergeLowDenCohorts, Do_DroughtMu, Do_RecoverSP, FreqY0,      &
  Do_ClosedN_run, Do_VariedKx, Do_variedWTC0, Do_mu0_F_WDen,   &
  Do_Fire, Do_FixedFrisk, Do_FixedFireS, Do_CH4,               &
  ! Specific test
  siteLAT,Scefile,StartLine,yr_ResetVeg,yr_Baseline,CO2Tag,    &
  PaleoPfile, PaleoTfile, iDraw

  ! ---------- Soil hydraulic and heat parameter name list ---------
  namelist /soil_data_nml/ soiltype, WaterLeakRate, thksl, GMD,   &
  GSD, vwc_sat, k_sat_ref, psi_sat_ref, chb,alphaSoil,       &
  heat_capacity_dry

  ! --------- Vegetation parameter name list ---------
  namelist /vegn_parameters_nml/  diff_S0,                      &
  pt, phenotype, lifeform, f_taper,f_cGap,                      &
  alphaHT,alphaCA,alphaBM,thetaHT,thetaCA,thetaBM,              &
  ! Leaf
  LMA,LAImax,Vmax,m_cond,Vannual,ps_wet,c_LLS,leaf_size,        &
  ! Wood and root
  rho_wood,rho_FR,root_r,root_zeta,root_perm, Kw_root,          &
  rho_N_up0, N_roots0,                                          &
  ! Growth & respiration
  f_iniBSW,f_LFR_max,GR_factor,LFR_rate,tauNSC,phiRL,phiCSA,    &
  R0_Nfix, C0_Nfix, S_facuN, f_N_add, fNSNmax, retransN,        &
  transT, l_fract, gamma_L, gamma_LN, gamma_SW, gamma_FR,       &
  MaxGrassLyr, MaxGrassAge, MaxGrassCA,                         &
  ! Phenology
  Tc0_OFF, Tc0_ON, T0_chill, betaON, betaOFF, AWD_crit,         &
  Days_thld, cold_thld, N0_GD,                                  &
  gdd_crit, T0_gdd, gdd_par1, gdd_par2, gdd_par3,               &
  ! Reproduction and Mortality
  AgeRepro,v_seed,s0_plant,prob_g,prob_e,                       &
  mu0_topL,D0mu,A_un,A_sd,B_sd,A_DBH,B_DBH,s_hu,W_mu0,          &
  ! Tisue C/N ratios
  LNbase, CN0leafST, CNleaf0, CNwood0, CNroot0, CNseed0,        &
  ! Plant hydraulics
  WTC0_WD,kx0_WD,psi0_WD,p50_WD,ths0_WD,fsup0_WD,CR0_LF,CR0_WD, &
  TK0_leaf,kx0, WTC0, psi0_LF,psi0_osm,r_DF,m0_WTC,m0_kx,       &
  fplc0_WD,A_plc0_WD,f_plc,plc_crit,                            &
  ! Soil BGC
  K0SOM,CUEmax0,fsc_fine,fsc_wood,f_M2SOM,                      &
  K_DeNitr, rho_SON, fDON, etaN, fdsvN,                         &
  CH4_alpha, CH4_beta_ox, CH4_wfps0, CH4_wfps1,                 &
  ! Fire model parameters, updated 11/25/2025
  EnvF0,MI0Fire,FSBM0,A_MI,f_bk,r_BK0,IgniteP,mu0fire,s0_max

  !---------------------------------
end module datatypes
