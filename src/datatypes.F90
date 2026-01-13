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
 integer, parameter :: NBANDS   = 2, & ! number of spectral bands for short-wave radiation calculations
                       BAND_VIS = 1, & ! visible radiation (wavelenght range?)
                       BAND_NIR = 2    ! near infra-red radiation (wavelenght range?)

 real,    parameter :: min_nindivs= 0.1E-4 ! 2e-15 ! 1/m. 2e-15 is approximately 1 individual per Earth

 ! Plant hydraulics-mortality
 integer, parameter :: Ysw_max      = 210 ! Maximum function years of xylems
 integer, parameter :: CLmax        = 5 ! Maximum crown layers
 real, parameter    :: WDref0       = 300.0   ! Reference wood density, kgC m-3
 real, parameter    :: rho_cellwall = 750.0 ! kgC m-3, Kellogg & Wangaard 1969 1.5 g/cc
 real, parameter    :: LMAmin       = 0.02    ! minimum LMA (kgC/m2), boundary condition

 ! Soil SOM reference C/N ratios
 integer, parameter :: N_SOM = 5
 !  fineL, structuralL, microbial, fast, slow
 real, parameter :: CN0SOM(5) = (/50., 150., 10., 15., 40./) ! target CN ratios of SOM
 real, parameter :: CUEmax0 = 0.4 ! CN0fastSOC  = 15.0 ! 25.0 ! 15.0

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
 integer, parameter :: & ! phenology type
          PHEN_DECIDIOUS = 0, PHEN_EVERGREEN = 1
 integer, parameter :: & ! physiology types
                   PT_C3 = 0,  PT_C4 = 1
 integer, parameter :: NSPECIES = 5, & ! number of species
          SP_C4GRASS   = 0, & ! c4 grass
          SP_C3GRASS   = 1, & ! c3 grass
          SP_TEMPDEC   = 2, & ! temperate deciduous
          SP_TROPICAL  = 3, & ! non-grass tropical
          SP_EVERGR    = 4    ! non-grass evergreen
 integer, parameter :: N_LU_TYPES = 4, & ! number of land use types
          LU_PAST    = 1, & ! pasture
          LU_CROP    = 2, & ! crops
          LU_NTRL    = 3, & ! natural vegetation
          LU_SCND    = 4    ! secondary vegetation
 integer, parameter :: CO2Yrs = 325
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
  real :: m_cond       ! factor of stomatal conductance
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
  ! wood traits
  real :: rho_wood     ! woody density, kg C m-3 wood
  real :: gamma_SW     ! sapwood respiration rate, kgC m-2 Acambium yr-1
  real :: f_taper

  ! Fire related
  real :: IgniteP ! Probability of ignition when climatic conditions are met (i.e., flammability).

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
  real :: CNsw0
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
  real :: NSN    = 0.    ! non-structural N pool
  real :: leafN  = 0.
  real :: swN    = 0.
  real :: hwN    = 0. ! N of heart wood
  real :: rootN  = 0. ! N of fine roots
  real :: seedN  = 0. !
  real :: N_uptake = 0.
  real :: fixedN ! fixed N at each stem per tree

  ! ---- water uptake-related variables
  real :: root_length(soil_L) ! m
  real :: rootarea ! total fine root area per tree
  real :: rootdepth  ! maximum depth of fine roots
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
  real :: N_input = 0.   ! annual N input (kgN m-2 yr-1)
  real :: N_uptake= 0.0  ! kg N m-2 hour-1
  real :: fixedN  = 0.0  ! kg N/step
  real :: Nm_Soil = 0.0  ! annual available N in a year
  real :: Nm_Fire = 0.0  ! Mineralized N due to burning
  real :: N_OutYr = 0.0  ! annual N loss
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
  real :: soilWater         ! kg m-2 in root zone

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
  real :: dailyNup
  real :: NfixDaily
  real :: dailyLFLIT = 0.0   !kgC day-1, leaf litter flux

  ! Annual diagnostics
  real :: CO2_c ! ppm, annual atmospheric CO2 concentration
  real :: dailyPrcp = 0.0, annualPrcp = 0.0 ! mm m-2 yr-1
  real :: dailyTrsp = 0.0, dailyEvap  = 0.0, dailyRoff = 0.0 ! mm m-2 yr-1
  real :: annualTrsp= 0.0, annualEvap = 0.0, annualRoff= 0.0 ! mm m-2 yr-1
  real :: annualGPP = 0.0 ! kgC m-2 ground yr-1
  real :: annualNPP = 0.0
  real :: annualResp = 0.0
  real :: annualRh   = 0.0
  real :: NupYr  = 0.0   ! accumulated N uptake kgN m-2 yr-1
  real :: NfixedYr = 0.  ! fixe N in a tile
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
real :: fsc_fine  = 1.0     ! fraction of fast turnover carbon in fine biomass
real :: fsc_wood  = 0.0     ! fraction of fast turnover carbon in wood biomass
real :: GR_factor = 0.33    ! growth respiration factor
real :: l_fract   = 0.0     ! 0.25  ! 0.5 ! fraction of the carbon retained after leaf drop
real :: retransN  = 0.0     ! retranslocation coefficient of Nitrogen
real :: f_iniBSW  = 0.1     !0.01, 0.2
real :: f_N_add   = 0.02    ! re-fill of N for sapwood
real :: f_LFR_max = 0.85    ! max allocation to leaves and fine roots each step
real :: c_LLS     = 28.5714 ! yr/(kg C m-2), 1/LMAs, ! Leaf life span: leafLS = c_LLS * LMA, (LMAs = 0.035., leafLS = 1.0)
real :: rho_N_up0 = 0.1     ! 0.05 ! hourly N uptake rate, fraction of the total mineral N
real :: N_roots0  = 0.4     ! root biomass at half max N-uptake rate,kg C m-2
real :: MaxGrassCA = 1.2 ! Maximum grass CA in the top layer
real :: MaxGrassAge = 3.0   ! Maximum grass  age (years)
integer :: MaxGrassLyr = 2    ! Maximal layers that grasses can survive

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
real :: A_MI         = 20.0  ! shape parameter of Fire risk vs. P/PET curve
real :: FSBM0        = 0.2   ! kgC m-2, grass fire severity parameter, as a function of grass BM
real :: mu0_FireG    = 0.2   ! mortality rates of grasses due to fire
real :: mu0_FireW    = 0.99  ! mortality rates of trees due to fire
real :: f_bk         = 0.1105! coefficient of bark thickness, Hoffmann et al. 2012.
                             ! shrubs: Y=1.105*X^1.083; trees: Y=0.31*X^1.276 for (Y:mm, X:cm)
real :: r_BK0        = -240.0! bark resistance, exponential equation, 120 --> 0.006 m of bark

! Soil organic matter decomposition
real :: K0SOM(5)     = (/0.8, 0.25, 2.5, 1.0, 0.2/) ! turnover rate of SOM pools (yr-1)
real :: K_DeNitr     = 8.0     ! mineral Nitrogen turnover rate
real :: fDON         = 0.02    ! fraction of DON production in decomposition
real :: rho_SON      = 0.05    ! SON release rate per year
real :: f_M2SOM      = 0.8     ! the ratio of C and N returned to litters from microbes
real :: etaN         = 0.025   ! Coefficient of N loss through runoff (etaN*runoff is a fraction of organic or mineral N)

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

! root parameters
real :: alpha_FR(0:MSPECIES) = 1.2 ! Fine root turnover rate yr-1
real :: rho_FR(0:MSPECIES)   = 200 ! woody density, kgC m-3
real :: root_r(0:MSPECIES)   = 2.9E-4 !(/1.1e-4, 1.1e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4/)
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
real :: gamma_L(0:MSPECIES) = 0.02 !
real :: gamma_LN(0:MSPECIES)= 70.5 ! 25.0  ! kgC kgN-1 yr-1
real :: gamma_SW(0:MSPECIES)= 0.02 ! 0.08 ! kgC m-2 Acambium yr-1
real :: gamma_FR(0:MSPECIES)= 0.6 ! 12 !kgC kgN-1 yr-1 ! 0.6: kgC kgN-1 yr-1

! Phenology parameters
real :: Tc0_OFF(0:MSPECIES) = 15. ! Critical temperature (C) for phenology OFF
real :: Tc0_ON(0:MSPECIES)  = 10. ! Critical temperature (C) for phenology ON
real :: gdd_crit(0:MSPECIES)= 300. ! Replaced by an equation: sp%gdd_par1 + sp%gdd_par2 * exp(sp%gdd_par3*cc%ncd)
real :: AWD_crit(0:MSPECIES)= 0.7  ! Critical plant water availability factor (0~1)
real :: betaON(0:MSPECIES)  = 0.2  ! Critical soil moisture for phenology ON
real :: betaOFF(0:MSPECIES) = 0.1  ! Critical soil moisture for phenology OFF
real :: gdd_par1(0:MSPECIES)= 30.0   !50.d0   ! -68.d0
real :: gdd_par2(0:MSPECIES)= 800. ! 650.d0  !800.d0  ! 638.d0
real :: gdd_par3(0:MSPECIES)= -0.02 ! -0.01d0
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
real :: CNleaf0(0:MSPECIES)  = 25. ! C/N ratios for leaves
real :: CNsw0(0:MSPECIES)    = 350.0 ! C/N ratios for woody biomass
real :: CNwood0(0:MSPECIES)  = 350.0 ! C/N ratios for woody biomass
real :: CNroot0(0:MSPECIES)  = 40.0 ! C/N ratios for leaves ! Gordon & Jackson 2000
real :: CNseed0(0:MSPECIES)  = 20.0 ! C/N ratios for seeds
real :: R0_Nfix(0:MSPECIES)  = 0.0  ! Reference N fixation rate (0.03 kgN kg rootC-1 yr-1)
real :: C0_Nfix(0:MSPECIES)  = 12.0 ! Carbon cost of N fixation, FUN model, Fisher et al. 2010, GBC; Kim

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
integer, parameter :: N_PFTs  = 9 ! Global PFTs from Ent Vegetation Map
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
   real    :: fPFT(N_PFTs) = 0.0 ! Fraction of each PFT's coverage
   real    :: SOM(5) ! Soil organic matter (kgC m-2)
   real    :: SON(5) ! Soil organic nitrogen
   real    :: mineralN ! Soil mineral N
   real    :: soiltexture(3)
   real    :: WLTPT, FLDCP ! soil wilting point and field capacity (0.xx)
   real, pointer :: climate(:,:)       ! Ntimes, N_vars
end type grid_initial_type

! Model mechanisms setting
logical  :: Do_DroughtMu        = .True. ! Drought-induced mortality, Lichstein (2024)
logical  :: MergeLowDenCohorts  = .True.
logical  :: update_annualLAImax = .False.
logical  :: do_RecoverSP        = .False.
logical  :: do_fire             = .False.
logical  :: do_closedN_run      = .True.
logical  :: do_VariedKx         = .True. ! trunk new xylem has the same kx or not
logical  :: do_VariedWTC0       = .True.
logical  :: do_WD_mort_function = .False.

! For global/regional run, Weng, 2025-07-22
character (len = 256) :: veg_path     = '/Users/eweng/Documents/Data/Vegetation/'
character (len = 20)  :: veg_file     = 'pft2011_0.5x0.5.nc' ! Vegetation coverage
character (len = 256) :: ncfilepath   = '/media/eweng/HD2/weng/Data/CRU/TRENDY2023/1HX1/'
character (len = 20)  :: ncversion    = 'crujra.v2.4.5d.'
character (len = 256) :: int_fpath    = '/media/eweng/HD2/weng/Data/CRU/TRENDY2023/1HX1/interpolated/'
character (len = 80)  :: int_prefix   = 'crujra.v2.4.5d.'
character (len = 50)  :: GridListFile = 'GlobalVegGridList.csv' ! in int_fpath
character (len = 5)   :: ncfields(7)  = [character(len=5):: 'tmp','pre','dswrf','spfh','pres','ugrd','vgrd']
character (len = 6)   :: GridIDFMT    = '(I6.6)' ! For the file name string (GridID)
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

type(grid_initial_type), pointer :: LandGrid(:) => null()
integer, pointer :: GridLonLat(:)    => null() ! LonLat
real,    pointer :: CRUData(:,:,:,:) => null() ! N_yr*Ntime, N_vars, Nlon, Nlat
real,    pointer :: ClimData(:,:,:)  => null() ! N_yr*Ntime, N_vars, N_VegGrids
real,    pointer :: CRUtime(:)       => null() ! Days since 1901-01-01 in CRU data

! Output interpolated grid climate data files
logical :: WriteForcing = .False. ! .True. ! Write interpolated forcing data
integer :: Grids_Unit   = 99    ! Vegetation grids list file

! Model run control
character(len=256) :: file_out(6) ! Output file names
integer  :: model_run_years = 100
integer  :: totyears, totdays
integer  :: steps_per_day  = 24 ! 24 or 48
integer  :: yr_ResetVeg    = 0 ! reseting vegetation to the initial, clearcut
integer  :: yr_Baseline    = 1000 ! for DroughtMIP baseline model run years
integer  :: F_Recovery     = 5 ! Interval (yrs) of recovering initial species
integer  :: equi_days      = 0 ! 100 * 365
integer  :: steps_per_hour = 1
real     :: step_hour      = 1.0  ! hour, Time step of forcing data, usually hourly (1.0)
real     :: step_seconds   = 1.0 * 3600.0
real     :: dt_fast_yr     = 1.0 / (365.0 * 24.0) ! Hourly
real     :: dt_daily_yr    = 1.0/365.0 ! Daily

! Model output
logical  :: outputhourly = .True.
logical  :: outputdaily  = .True.

! Scenarios
character(len=4)  :: CO2Tag = 'aCO2' ! only takes 'aCO2' or 'eCO2', for FACE-MDS-3
real     :: Sc_prcp = 1.0 ! Scenario of rainfall changes
real     :: Sc_dT   = 0.0 ! Scenario of temperature changes
real     :: CO2_c   = 375.0 ! 412 ! PPM, CO2 concentration at 2020

!-------------Plant and soil parameter types -----------------------
type(spec_data_type), save :: spdata(0:MSPECIES)         ! PFT-specific parameters
type(soil_pars_type), save :: soilpars(n_dim_soil_types) ! Soil hydraulics parameters

! ------------- Global run model setting name list ------------
namelist /global_setting_nml/ ncfilepath, ncversion,            &
    veg_path, veg_file, int_fpath, int_prefix, GridListFile,    &
    grid_No1, grid_No2, yr_start, yr_end, LowerLon, UpperLon,   &
    LowerLat, UpperLat, StepLatLon, WriteForcing

! ------------- Model initialization name list ------------
namelist /initial_state_nml/ &
    ! initial vegetation states
    init_cohort_N, init_cohort_sps, init_cohort_Indiv,          &
    init_cohort_bl, init_cohort_br, init_cohort_bsw,            &
    init_cohort_bHW, init_cohort_seedC, init_cohort_nsc,        &
    init_fast_SOC, init_slow_SOC, init_mineralN, N_input,       &
    Pr_thld, MI0DeSB, MI0C3C4, TcrTREE, TcrC3C4,                &
    ! Model run controls
    filepath_in,filepath_out,runID,climfile,Scefile,StartLine,  &
    PaleoPfile, PaleoTfile, iDraw,                              &
    N_VegTile,siteLAT,model_run_years,yr_ResetVeg,yr_Baseline,  &
    outputhourly,outputdaily,Sc_prcp,Sc_dT,CO2_c,CO2Tag,        &
    update_annualLAImax, MergeLowDenCohorts,F_Recovery,         &
    Do_DroughtMu, do_RecoverSP, do_closedN_run, do_fire,        &
    do_VariedKx, do_variedWTC0, do_WD_mort_function

! ---------- Soil hydraulic and heat parameter name list ---------
namelist /soil_data_nml/ soiltype, WaterLeakRate, thksl, GMD,   &
     GSD, vwc_sat, k_sat_ref, psi_sat_ref, chb,alphaSoil,       &
     heat_capacity_dry

! --------- Vegetation parameter name list ---------
namelist /vegn_parameters_nml/  diff_S0,                        &
  pt, phenotype, lifeform, f_taper,f_cGap,                      &
  alphaHT,alphaCA,alphaBM,thetaHT,thetaCA,thetaBM,              &
  ! Leaf
  LMA,LAImax,Vmax,m_cond,Vannual,ps_wet,c_LLS,leaf_size,        &
  ! Wood and root
  rho_wood,rho_FR,root_r,root_zeta,root_perm, Kw_root,          &
  rho_N_up0, N_roots0,                                          &
  ! Growth & respiration
  f_iniBSW,f_LFR_max,GR_factor,LFR_rate,tauNSC,phiRL,phiCSA,    &
  R0_Nfix, C0_Nfix, f_N_add, fNSNmax, transT, l_fract,retransN, &
  gamma_L, gamma_LN, gamma_SW, gamma_FR,                        &
  MaxGrassLyr, MaxGrassAge, MaxGrassCA,                         &
  ! Phenology
  Tc0_OFF, Tc0_ON, T0_chill, betaON, betaOFF, AWD_crit,         &
  Days_thld, cold_thld, N0_GD,                                  &
  gdd_crit, T0_gdd, gdd_par1, gdd_par2, gdd_par3,               &
  ! Reproduction and Mortality
  AgeRepro,v_seed,s0_plant,prob_g,prob_e,                       &
  mu0_topL,D0mu,A_un,A_sd,B_sd,A_DBH,B_DBH,s_hu,W_mu0,          &
  ! Tisue C/N ratios
  LNbase,CN0leafST,CNleaf0,CNsw0,CNwood0,CNroot0,CNseed0,       &
  ! Plant hydraulics
  WTC0_WD,kx0_WD,psi0_WD,p50_WD,ths0_WD,fsup0_WD,CR0_LF,CR0_WD, &
  TK0_leaf,kx0, WTC0, psi0_LF,psi0_osm,r_DF,m0_WTC,m0_kx,       &
  fplc0_WD,A_plc0_WD,f_plc,plc_crit,                            &
  ! Soil
  K0SOM,fsc_fine,fsc_wood,K_DeNitr,rho_SON,f_M2SOM,fDON,etaN,   &
  ! Fire model parameters, updated 11/25/2025
  EnvF0,MI0Fire,FSBM0,A_MI,mu0_FireW,mu0_FireG,f_bk,r_BK0,IgniteP

!---------------------------------
 contains

!==========================================================================
! Initialize model processes/settings and model parameters
!----------------------------------------------------------------
subroutine model_para_init(fnml)
  character(len=*),intent(in) :: fnml
  call initialize_soilpars(fnml)
  call initialize_PFT_pars(fnml)
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
  UpperLat = max(UpperLat, Nlat)
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
   gamma_SW(0:N_EST)  = [0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02,   0.02  ] ! Wood Acambium respiration rate (kgC/m2/yr
   Tc0_OFF(0:N_EST)   = [12.,    8.,     15.,    15.,    -50.,   12.,    15.0,   12.   ] ! 283.16 ! OFF ! C for convenience
   Tc0_ON(0:N_EST)    = [10.,    5.,     10.,    10.,    -50.,   8.,     10.,     8.   ] ! 280.16 ! ON  ! C for convenience
   betaON(0:N_EST)    = [0.3,    0.3,    0.0,    0.6,    0.0,    0.3,    0.6,    0.2   ] ! Critical soil moisture for phenology ON
   betaOFF(0:N_EST)   = [0.2,    0.2,    0.0,    0.4,    0.0,    0.2,    0.4,    0.1   ] ! Critical soil moisture for phenology OFF
   gdd_par1(0:N_EST)  = [50.,    20.,    0.0,    0.0,    0.0,    50.,    50.,    50.   ] ! 50.d0   ! These three parameters are used to calculate gdd_crit
   gdd_par2(0:N_EST)  = [800.,   600.,   0.0,    600.,   0.0,    600.,   600.,   600.  ] ! 650.d0  !800.d0  ! 638.d0
   gdd_par3(0:N_EST)  = [-0.02,  -0.02,  -0.02,  -0.02,  -0.02,  -0.02,  -0.02,  -0.02 ] ! -0.01d0
   R0_Nfix(0:N_EST)   = [0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.03,   0.0   ] ! Nitrogen fixation rate, 0.03 kgN kgRootC-1 yr-1
   C0_Nfix(0:N_EST)   = [12.0,   12.0,   12.0,   12.0,   12.0,   12.0,   12.0,   12.0  ] ! N fixation carbon cost: 12 gC/gN

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
  spdata%CNsw0    = CNsw0
  spdata%CNwood0  = CNwood0
  spdata%CNroot0  = CNroot0
  spdata%CNseed0  = CNseed0
  spdata%R0_Nfix  = R0_Nfix
  spdata%C0_Nfix  = C0_Nfix

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
   if(do_WD_mort_function)then
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
 subroutine Set_PFTs_from_map(LandGrid)
   type(grid_initial_type), intent(in) :: LandGrid

   !--------- local vars ------------
   integer :: GridPFTs(N_PFTs)
   integer :: i
   real :: f_min = 0.01 ! coverage fraction threshold

   ! Sorting PFT numbers according to fPFT
   call rank_descending(LandGrid%fPFT,GridPFTs)
   !PFTID = [character(len=3) :: 'C4G','C3G','TEB','TDB','EGN','CDB','CDN','CAS','AAS']
   GridPFTs = GridPFTs - 1 ! PFT No. starts from 0.

   ! Find out PFTs in this grid
   init_cohort_N = min(M_initialCH,Max(1, COUNT(LandGrid%fPFT > f_min)))
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
 ! for testing tree-grass-desert shrub and evergreen-deciduous forests only
 ! Weng, 09/06/2025
 subroutine Vegn_PFTs_from_Climate (forcingData,steps_per_day)
   implicit none
   type(climate_data_type), intent(in) :: forcingData(:)
   integer, intent(in) :: steps_per_day

   !--------- local vars ------------
   integer, parameter :: N_PFTID = 4, mw = 15
   integer :: PFTID(N_PFTID) ! 0:C4, 1:C3, 2:TrE, 3:TrD, 4:TmE, 5:TmD, 6:N-fixer, 7:Desert shrub
   real, allocatable :: dailyET(:),dailyTc(:),dailyPr(:),meanTc(:),TminYr(:)
   real :: tmpL(mw*2+1)
   real :: totPrcp, totPET, meanTmin, meanPrcp, Mst_IDX

   integer :: i,j,k,m,n,L,w
   integer :: N_Yrs, N_days

   ! ------ Calculate days and years of the data ------
   N_days = size(forcingData)/steps_per_day
   N_yrs = N_days/365
   if (N_yrs < 1) then
      error stop "Vegn_PFTs_from_Climate: need at least one year's data."
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
   Mst_IDX  = totPrcp / totPET ! Moisture Index
   meanPrcp = totPrcp / N_yrs

   ! Calculate running mean temperature (meanTc) with a windlow of L
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
   write(*,'(a15, 2(f8.2,","))')'Prcp, PET: ', totPrcp/N_yrs,totPET/N_yrs
   write(*,'(2(a6,f8.2,";"), a12, 4(I6,","))')   &
      'P/ET: ', Mst_IDX, 'Tmin: ', meanTmin, 'Grid PFTs: ', PFTID

   ! Assign standard initial cohorts by updating init_cohort_*
   call Assign_Std_Cohorts (PFTID,N_PFTID)

   ! Release allocatable variables
   deallocate(dailyET,dailyTc,dailyPr,meanTc,TminYr)
 end subroutine Vegn_PFTs_from_Climate

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
     real,intent(in) :: td             ! day(to minute fraction)
     real,intent(in) :: latdegrees     ! latitude in degrees
     real :: hour,latrad
     real :: delta    ! declination angle
     real :: pi, rad
     real,intent(out) :: cosz        ! cosz=cos(zen angle)=sin(elev angle)
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
end module datatypes
