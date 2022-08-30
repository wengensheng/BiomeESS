module datatypes
 ! define data types, constants, and some basic functions
 implicit none

 ! ---- public types -------
 public :: spec_data_type, cohort_type, vegn_tile_type

 ! ------ public subroutines and functions ---------
 public :: initialize_PFT_data, initialize_soilpars
 public :: BM2Architecture,DBH2HT,DBH2CA,DBH2BM,BM2DBH, &
           CA2BLmax,BLmax2BRmax,ccNSNmax,BL2Aleaf,      &
           qscomp,calc_solarzen,esat, A_function

 ! ------ public namelists ---------
 public :: vegn_parameters_nml, soil_data_nml, initial_state_nml

 ! ---- public variables ---------
 public :: forcingData, spdata, soilpars

!===============constants===============
 logical, public, parameter :: read_from_parameter_file = .TRUE.
 integer, public, parameter :: days_per_year  = 365
 integer, public, parameter :: hours_per_year = 365 * 24  ! 8760
 real,    public, parameter :: seconds_per_year = 365. * 24. * 3600.
 real,    public, parameter :: seconds_per_day = 24. * 3600.

 ! Physical constants
 real, public, parameter :: TFREEZE = 273.16
 real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1
 real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg
 real, public, parameter :: mol_air = 28.96440e-3 ! molar mass of air, kg
 real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
 real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg
 real, public, parameter :: cpair = 1010.
 real, public, parameter :: H2OLv0=2.501e6   !latent heat H2O (J/kg)
 real, public, parameter :: p_sea = 101325.  ! atmospheric pressure  (Pa)
 real, public, parameter :: DENS_H2O = 1000. ! kg m-3
 real, public, parameter :: PI = 3.1415926

 ! Vegetation and soil types
 integer, public, parameter :: MSPECIES = 15
 integer, public, parameter :: Nsoiltypes = 7
 integer, public, parameter :: n_dim_soil_types = 9
 integer, public, parameter :: LEAF_ON  = 1
 integer, public, parameter :: LEAF_OFF = 0
 real,    public, parameter :: min_nindivs = 1e-5 ! 2e-15 ! 1/m. 2e-15 is approximately 1 individual per Earth

 ! Plant hydraulics-mortality
 real, public, parameter    :: rho_H2O = 1000.0 ! water density (kg m-3)
 real, public, parameter    :: WDref0 = 300.0   ! Reference wood density, kgC m-3
 real, public, parameter    :: rho_cellwall = 750.0 ! kgC m-3, Kellogg & Wangaard 1969 1.5 g/cc
 integer, public, parameter :: Ysw_max = 210 ! Maximum function years of xylems

 ! Soil SOM reference C/N ratios
 integer, public, parameter :: N_SOM = 5
 !  fineL, structuralL, microbial, fast, slow
 real, public, parameter :: CN0SOM(5) = (/50., 150., 10., 15., 40./) ! target CN ratios of SOM
 real, public, parameter :: CUEmax0 = 0.4 ! CN0fastSOC  = 15.0 ! 25.0 ! 15.0

 ! Soil water hydrualics
 integer,public,parameter :: soil_L = 5 ! Soil layers, for soil water dynamics
 integer,public,parameter :: num_l   = soil_L ! Soil layers,
 real, public, parameter :: thksl(soil_L)=(/0.05,0.45,1.5,1.5,1.5/) ! m, thickness of soil layers
 real, public, parameter :: rzone = sum(thksl) !m
 real, public, parameter :: psi_wilt  = -150.0  ! matric head at wilting
 real, public, parameter :: K_rel_min = 1.e-12
 real, public, parameter :: rate_fc   = 0.1/86400 ! 0.1 mm/d drainage rate at FC
 real, public, parameter :: ws0 = 0.02 ! hygroscopic point
 real, public, parameter :: Edepth = 0.05 !m, the depth of soil for surface evaporation
 integer, public, parameter :: & ! soil types
              Sand        = 1,  LoamySand   = 2, &
              SandyLoam   = 3,  SiltLoam    = 4, &
              FrittedClay = 5,  Loam        = 6, &
              Clay        = 7

 ! LM3 vegetation and soil types
 integer, public, parameter :: & ! phenology type
          PHEN_DECIDIOUS = 0, PHEN_EVERGREEN = 1
 integer, public, parameter :: & ! physiology types
                   PT_C3 = 0,  PT_C4 = 1
 integer, public, parameter :: NSPECIES = 5, & ! number of species
          SP_C4GRASS   = 0, & ! c4 grass
          SP_C3GRASS   = 1, & ! c3 grass
          SP_TEMPDEC   = 2, & ! temperate deciduous
          SP_TROPICAL  = 3, & ! non-grass tropical
          SP_EVERGR    = 4    ! non-grass evergreen
 integer, public, parameter :: N_LU_TYPES = 4, & ! number of land use types
          LU_PAST    = 1, & ! pasture
          LU_CROP    = 2, & ! crops
          LU_NTRL    = 3, & ! natural vegetation
          LU_SCND    = 4    ! secondary vegetation

!===============data types ==============================
!-----------PFT data type----------------
type spec_data_type
  integer :: lifeform     ! 0 for grasses, 1 for trees
  integer :: phenotype    ! phenology type: 0 for deciduous, 1 for evergreen
  integer :: pt           ! photosynthetic physiology of species
  ! leaf traits
  real :: LMA          ! leaf mass per unit area, kg C/m2
  real :: leafLS       ! leaf life span
  real :: alpha_L      ! leaf turn over rate
  real :: LNA          ! leaf Nitrogen per unit area, kg N/m2
  real :: LNbase       ! basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
  real :: CN0leafST    ! leaf structural tissues, 175
  real :: leaf_size    ! characteristic leaf size
  real :: leafTK      ! leaf thickness, m
  real :: rho_leaf     ! leaf mass density (kgC/m3)
  real :: alpha_ps   ! photosynthesis efficiency
  real :: m_cond       ! factor of stomatal conductance
  real :: Vmax         ! max rubisco rate, mol m-2 s-1
  real :: Vannual      ! annual productivity per unit area at full fun (kgC m-2 yr-1)
  real :: gamma_L      ! leaf respiration coeficient (per yr)
  real :: gamma_LN     ! leaf respiration coeficient per unit N
  real :: ps_wet ! wet leaf photosynthesis down-regulation
  real :: LFR_rate      ! Leaf filling rate per day
  ! root traits
  real :: rho_FR       ! material density of fine roots (kgC m-3)
  real :: root_r       ! radius of the fine roots, m
  real :: root_zeta    ! e-folding parameter of root vertical distribution (m)
  real :: root_frac(soil_L)    ! root fraction
  real :: SRA          ! speific fine root area, m2/kg C
  real :: SRL          ! specific root lenght
  real :: gamma_FR     ! Fine root respiration rate, kgC kgC-1 yr-1
  real :: alpha_FR     ! Turnover rate of Fine roots, fraction yr-1
  real :: Kw_root      ! fine root water conductivity mol m m-2 s−1 Pa−1 !
  real :: root_perm
  !  real :: rho_N_up0   ! maximum N uptake rate
  !  real :: N_roots0    ! root biomass at half of max. N-uptake rate
  real :: NfixRate0    ! Reference N fixation rate (kgN kgC-1 root)
  real :: NfixCost0    ! Carbon cost of N fixation (kgC kgN-1)
  ! wood traits
  real :: rho_wood     ! woody density, kg C m-3 wood
  real :: gamma_SW     ! sapwood respiration rate, kgC m-2 Acambium yr-1
  real :: f_taper

  ! Plant hydraulics
  real :: kx0  ! xylem conductivity, (mm/s)/(Mpa/m)
  real :: WTC0 ! xylem water transfer capacity, m/lifetime
  real :: CR_Leaf ! leaf compression ratio per MPa
  real :: CR_Wood ! Wood compression ratio per MPa
  real :: w0L_max ! leaf maximum water/carbon ratio
  real :: w0S_max ! stem maximum water/carbon ratio
  real :: w0L_min ! leaf minimum water/carbon ratio
  real :: w0S_min ! stem minimum water/carbon ratio
  real :: psi0_LF ! minimum leaf water potential
  real :: psi0_WD ! minimum stem wood potential
  real :: psi50_WD !wood potential at which 50% conductivity lost, MPa
  real :: Kexp_WD  ! exponent of the PLC curve
  real :: f_supply ! fraction of stem water available for leaves per hour

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
  real :: tc_crit_off     ! K, for turning OFF a growth season
  real :: tc_crit_on      ! K, for turning ON a growth season
  real :: gdd_crit        ! K, critical value of GDD5 for turning ON growth season
  !  vital rates
  real :: AgeRepro       ! the age that can reproduce
  real :: v_seed           ! fracton of G_SF to G_F
  real :: s0_plant     ! size of the seedlings, kgC/indiv
  real :: prob_g,prob_e    ! germination and establishment probabilities
  real :: r0mort_c     ! yearly mortality rate in canopy
  real :: D0mu         ! Reference diameter for size-dependent mortality
  real :: A_un         ! Parameter for understory mortality affected by layers
  real :: A_sd         ! Max multiplier for seedling mortality
  real :: B_sd         ! Mortality sensitivity for seedlings
  real :: A_D          ! Sensitivity to dbh
  real :: s_hu         ! hydraulic mortality sensitivity
  ! Population level variables
  real :: LAImax    ! max. LAI
  real :: LAImax_u  ! max. LAI understorey
  real :: LAI_light ! light controlled maximum LAI
  integer :: n_cc   ! for calculating LAImax via cc%LAImax derived from cc%NSN
  real :: f_cGap    ! fraction of internal gaps in the canopy
  ! "internal" gaps are the gaps that are created within the canopy by the
  ! branch fall processes.
end type

!----------cohort-----------------
type :: cohort_type
  ! for climate-vegetation type
  integer :: phenotype    ! phenology type: 0 for deciduous, 1 for evergreen
  logical :: firstday     ! First day of a growing season
  integer :: pt           ! photosynthetic physiology of species

  ! ---- biological prognostic variables
  integer :: ccID   = 0   ! cohort ID
  integer :: species= 0   ! vegetation species
  real :: gdd       = 0.0   ! for phenology
  real :: ALT       = 0.0  ! growing season accumulative cold temperature
  integer :: Ngd    = 0   ! growing days
  integer :: Ndm    = 0   ! dormant days
  integer :: Ncd    = 0   ! number of cold days in non-growing season
  integer :: status = 0   ! growth status of plant: 1 for ON, 0 for OFF
  integer :: layer  = 1   ! the layer of this cohort (numbered from top, top layer=1)
  real :: layerfrac = 0.0 ! fraction of layer area occupied by this cohort
  real :: leaf_age  = 0.0 ! leaf age (year)

  ! for populatin structure
  real :: nindivs= 1.0 ! density of vegetation, individuals/m2
  real :: mu     = 0.02 ! Cohort mortality rate
  real :: age    = 0.0 ! age of cohort, years
  real :: dbh    = 0.0 ! diameter at breast height, m
  real :: height = 0.0 ! vegetation height, m
  real :: Acrown = 1.0 ! crown area, m2/individual
  real :: Aleaf  = 0.0 ! total area of leaves, m2/individual
  real :: lai    = 0.0 ! crown leaf area index, m2/m2
  real :: D_bark = 0.0 ! thickness of bark
  ! carbon pools
  real :: bl     = 0.0 ! biomass of leaves, kg C/individual
  real :: br     = 0.0 ! biomass of fine roots, kg C/individual
  real :: bsw    = 0.0 ! biomass of sapwood, kg C/individual
  real :: bHW    = 0.0 ! biomass of heartwood, kg C/individual
  real :: seedC  = 0.0 ! biomass put aside for future progeny, kg C/individual
  real :: nsc    = 0.0 ! non-structural carbon, kg C/individual

  ! ----- carbon fluxes
  real :: gpp  = 0.0 ! gross primary productivity kg C/step
  real :: npp  = 0.0 ! net primary productivity kg C/step
  real :: resp = 0.0 ! plant respiration
  real :: resl = 0.0 ! leaf respiration
  real :: resr = 0.0 ! root respiration
  real :: resg = 0.0 ! growth respiration
  real :: NPPleaf,NPProot,NPPwood ! to record C allocated to leaf, root, and wood

  ! for hydraulics-mortality
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
  real :: V_stem ! Volumn of stems (including trunk)
  real :: V_leaf ! Volumn of leaves
  real :: Q_stem ! water flux from soil to stems (kg/tree/step)
  real :: Q_leaf ! water flux from stems to leaves (kg/tree/step)

  real :: Ktrunk ! trunk water conductance, m/(s MPa)
  real :: Asap ! Functional cross sectional area
  real :: Atrunk ! Sum of all rings
  real :: treeHU ! total water transported by the functional sapwood, m^3
  real :: treeW0 ! total WTC0 of the sapwood, m^3
  real :: Kx(Ysw_max) = 0.0 ! Initial conductivity of the woody generated in each year
  real :: WTC0(Ysw_max) = 0.0 ! lifetime water transfer capacity
  real :: accH(Ysw_max) = 0.0 ! m, total water transport for functional conduits
  real :: farea(Ysw_max) = 0.0 ! fraction of functional area, 1.0/(exp(r_DF*(1.0-accH[j]/W0[j]))+1.0)
  real :: Rring(Ysw_max) = 0.0 ! Radius to the outer edge
  real :: Lring(Ysw_max) = 0.0 ! Length of xylem conduits
  real :: Aring(Ysw_max) = 0.0 ! Area of each ring
  real :: Kring(Ysw_max) = 0.0 ! Conductance of each ring

  ! for diagnostics
  real :: dailyTrsp
  real :: dailyGPP   ! kgC/tree day-1
  real :: dailyNPP
  real :: dailyResp
  real :: dailyNup
  real :: annualTrsp
  real :: annualGPP ! C flux/tree
  real :: annualNPP
  real :: annualResp

  ! ---- Nitrogen model related parameters
  real :: NSNmax = 0.
  real :: NSN = 0.    ! non-structural N pool
  real :: leafN = 0.
  real :: sapwN= 0.
  real :: woodN = 0. ! N of heart wood
  real :: rootN = 0. ! N of fine roots
  real :: seedN = 0. !
  real :: N_uptake = 0.
  real :: annualNup  = 0.0
  real :: fixedN ! fixed N at each stem per tree
  real :: dailyfixedN
  real :: annualfixedN = 0.0 ! annual N fixation per unit crown area

  real :: bl_max  = 0.0 ! Max. leaf biomass, kg C/individual
  real :: br_max  = 0.0 ! Max. fine root biomass, kg C/individual
  real :: CSAsw   = 0.0
  real :: topyear = 0.0 ! the years that a plant in top layer
  real :: DBH_ys        ! DBH at the begining of a year (growing season)

  ! ---- water uptake-related variables
  real :: root_length(soil_L) ! m
  real :: rootarea ! total fine root area per tree
  real :: rootdepth  ! maximum depth of fine roots
  real :: rootareaL(soil_L) = 0.0 ! Root length per layer, m of root/m
  real :: WupL(soil_L) = 0.0 ! normalized vertical distribution of uptake
  real :: Q_soil(soil_L) = 0.0 ! Soil to roots water flux (kg H2O/tree/step)
  real :: W_supply  ! potential water uptake rate per unit time per tree
  real :: transp   ! transpiration rate per tree per hour
  real :: uptake_frac(soil_L) ! for LM3 soil water uptake, Weng, 2017-10-28
  real :: K_r,r_r
  real :: root_zeta
  ! for photosynthesis
  real :: An_op = 0.0 ! mol C/(m2 of leaf per year)
  real :: An_cl = 0.0 ! mol C/(m2 of leaf per year)
  real :: w_scale =-9999
  real :: C_growth = 0.0 ! carbon gain since last growth, kg C/individual
  real :: N_growth = 0.0 ! Nitrogen used for plant tissue growth
  real :: extinct = 0.75     ! light extinction coefficient in the canopy for photosynthesis

end type cohort_type

!---------------------------
type :: vegn_tile_type
  integer :: tag ! kind of the tile
  integer :: landuse = LU_NTRL
  integer :: n_cohorts = 0
  integer :: n_initialCC = 0
  integer :: n_years   = 0
  integer :: n_canopycc = 0
  type(cohort_type), pointer :: cohorts(:)=>NULL()
  type(cohort_type), pointer :: initialCC(:)=>NULL()
  real :: area      ! m2
  real :: age = 0.0 ! tile age
  real :: LAI  ! leaf area index
  real :: CAI  ! crown area index
  real :: LAIlayer(9) = 0.0 ! LAI of each crown layer, max. 9
  real :: f_gap(9)    = 0.0 ! gap fraction of each crown layer
  real :: kp(9)       = 0.0 ! light extinction coefficient fro each layer
  ! uptake-related variables
  real :: root_distance(soil_L) ! characteristic half-distance between fine roots, m
  ! averaged quantities for PPA phenology
  real :: tc_daily = 0.0
  real :: tc_pheno = 0.0 ! smoothed canopy air temperature for phenology

  ! litter and soil carbon pools
  real :: litter = 0.0 ! litter flux
  real :: SOC(5) = 0. ! metabolicL, structuralL, microbial, fastSOM, slowSOM
  real :: SON(5) = 0.

  !!  Nitrogen pools, Weng 2014-08-08
  real :: mineralN = 0.  ! Mineral nitrogen pool, (kg N/m2)
  real :: totN = 0.
  real :: N_input        ! annual N input (kgN m-2 yr-1)
  real :: N_uptake= 0.0  ! kg N m-2 hour-1
  real :: fixedN  = 0.0  ! kg N/step
  real :: annualN = 0.0  ! annual available N in a year
  real :: Nloss_yr= 0.0  ! annual N loss
  real :: N_P2S_yr= 0.0  ! annual N from plants to soil
  real :: previousN      ! an weighted annual available N
  real :: initialN0

  ! Soil water
  integer :: soiltype    ! lookup table for soil hydrologic parameters
  real :: FLDCAP  ! soil field capacity
  real :: WILTPT  ! soil wilting point (0.xx)
  real :: evap           ! kg m-2 per unit fast time step (mm/hour)
  real :: transp         ! kg m-2 hour-1
  real :: runoff        ! Water runoff of the veg tile, unit?
  real :: thetaS     ! moisture index (ws - wiltpt)/(fldcap - wiltpt)
  real :: wcl(soil_L)   ! volumetric soil water content for each layer
  real :: freewater(soil_L) ! Available water in each layer
  real :: psi_soil(soil_L) ! MPa
  real :: K_soil(soil_L)   ! Kg H2O/(m2 s MPa)
  real :: soilWater      ! kg m-2 in root zone

  ! Vegetation water
  real :: W_leaf
  real :: W_stem
  real :: W_dead

  ! water uptake-related variables
  real :: RAI ! root area index
  real :: RAIL(soil_L) = 0.0 ! Root length per layer, m of root/m
  real :: W_uptake  ! water uptake rate per unit time per m2

  !  Carbon fluxes
  real :: gpp =0 ! gross primary production, kgC m-2 yr-1
  real :: npp =0 ! net primary productivity
  real :: resp = 0 ! auto-respiration of plants
  real :: rh  =0 ! soil carbon lost to the atmosphere

  !  fire disturbance
  real :: C_combusted = 0.0 ! Carbon released to atmosphere via fire
  real :: treecover   = 0.0 ! tree CAI in the top layer, for fire spread
  real :: grasscover  = 0.0 ! grass CAI, for the initial fire severity
  ! daily diagnostics
  real :: dailyGPP
  real :: dailyNPP
  real :: dailyResp
  real :: dailyRh
  real :: dailyNup
  real :: dailyfixedN
  ! for annual diagnostics
  real :: dailyPrcp=0.0, annualPrcp = 0.0 ! mm m-2 yr-1
  real :: dailyTrsp=0.0,dailyEvap=0.0, dailyRoff=0.0 ! mm m-2 yr-1
  real :: annualTrsp=0.0,annualEvap=0.0, annualRoff=0.0 ! mm m-2 yr-1
  real :: annualGPP = 0.0 ! kgC m-2 ground yr-1
  real :: annualNPP = 0.0
  real :: annualResp = 0.0
  real :: annualRh   = 0.0
  real :: annualNup  = 0.0   ! accumulated N uptake kgN m-2 yr-1
  real :: annualfixedN = 0.  ! fixe N in a tile
  ! for annual reporting at tile level
  real :: NSC, SeedC, leafC, rootC, SapwoodC, WoodC
  real :: NSN, SeedN, leafN, rootN, SapwoodN, WoodN
  real :: totSeedC,totSeedN
  ! for cohort plant types (climate-vegetation relationship, Biome, LM3)
  real :: t_ann  = 0.0 ! annual mean T, degK
  real :: t_cold = 0.0 ! average temperature of the coldest month, degK
  real :: p_ann  = 0.0 ! annual mean precip
  real :: ncm    = 0.0 ! number of cold months
end type vegn_tile_type

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
  real :: heat_capacity_dry
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
   integer :: year          ! Year
   integer :: doy           ! day of the year
   real :: hod           ! hour of the day
   real :: PAR           ! umol m-2 s-1
   real :: radiation     ! W/m2
   real :: Tair          ! air temperature,  K
   real :: Tsoil         ! soil temperature, K
   real :: RH            ! relative humidity
   real :: rain          ! kgH2O m-2 s-1
   real :: windU         ! wind velocity (m s-1)
   real :: P_air         ! pa
   real :: CO2           ! ppm
   real :: soilwater     ! soil moisture, vol/vol
end type climate_data_type
!---------------------------

! Input forcing data
type(climate_data_type),pointer, save :: forcingData(:)
! PFT-specific parameters
type(spec_data_type), save :: spdata(0:MSPECIES) ! define PFTs
! Soil hydraulic parameters
type(soil_pars_type), save :: soilpars(n_dim_soil_types) ! soil parameters

!---------------------------------------
integer :: MaxCohortID = 0

! Vegetation and soil BGC parameters
! Cohort management
real :: diff_S0 = 0.2 ! percentage of the difference between cohorts for merging

! Growth parameters:
real :: LMAmin     = 0.02    ! minimum LMA, boundary condition
real :: fsc_fine   = 1.0     ! fraction of fast turnover carbon in fine biomass
real :: fsc_wood   = 0.0     ! fraction of fast turnover carbon in wood biomass
real :: GR_factor  = 0.33 ! growth respiration factor
real :: l_fract    = 0.0 ! 0.25  ! 0.5 ! fraction of the carbon retained after leaf drop
real :: retransN   = 0.0   ! retranslocation coefficient of Nitrogen
real :: f_iniBSW   = 0.1  !0.01, 0.2
real :: f_N_add    = 0.02 ! re-fill of N for sapwood
real :: f_LFR_max  = 0.85 ! max allocation to leaves and fine roots each step
real :: c_LLS      = 28.57143 ! yr/(kg C m-2), 1/LMAs, ! for leaf life span and LMA leafLS = c_LLS * LMA, when LMAs = 0.035.
real :: rho_N_up0  = 0.1 ! 0.05 ! hourly N uptake rate, fraction of the total mineral N
real :: N_roots0   = 0.4  ! root biomass at half max N-uptake rate,kg C m-2

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
real :: f0_WD    = 0.8     ! Fraction of stem water for transpiration per hour at zero WD

! Mortality as a function of wood density
real :: A_mort   = 0.2    ! mu = A_mort *exp(B_mort*WD/WDref)
real :: B_mort   = -2.1

! Phenology parameters
! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)
real :: T0_gdd   = 273.15 + 5.0 ! 5.d0
real :: T0_chill = 273.15 + 10.0
integer :: N0_GD = 90 ! base growing days, 90 days, with a -5 substraction of Tc_crit
real :: gdd_par1 = 30.0   !50.d0   ! -68.d0
real :: gdd_par2 = 800. ! 650.d0  !800.d0  ! 638.d0
real :: gdd_par3 = -0.02 ! -0.01d0

! Fire regimes
real :: envi_fire_prb = 0.5 ! fire probability due to environment
real :: Ignition_G0 = 1.0   ! Ignition probability for grasses once meets envi_fire_prb
real :: Ignition_W0 = 0.025 ! Ignition probability for woody plants once meets envi_fire_prb
real :: m0_g_fire = 0.2     ! mortality rates of grasses due to fire
real :: m0_w_fire = 0.99    ! mortality rates of trees adue to fire
real :: r_BK0 = -240.0 ! -480.0  ! for bark resistance, exponential equation, 120 --> 0.006 m of bark
! An old scheme
real :: f_HT0 = 10.0 ! shape parameter fire resistence (due to growth of bark) as a function of height
real :: h0_escape = 5.0 ! tree height that escapes direct burning of grass fires
real :: D_BK0   = 5.9/1000.0 ! half survival bark thickness, m

! Soil water properties
real :: soiltype = SandyLoam  ! 1 Sand; 2
real :: FLDCAP = 0.4  ! vol/vol
real :: WILTPT = 0.05 ! vol/vol

! Soil organic matter decomposition
real :: K0SOM(5)  = (/0.8, 0.25, 2.5, 2.0, 0.05/) ! turnover rate of SOM pools (yr-1)
real :: K_nitrogen = 8.0     ! mineral Nitrogen turnover rate
real :: fDON       = 0.02    ! fraction of DON production in decomposition
real :: rho_SON    = 0.05    ! SON release rate
real :: f_M2SOM    = 0.8     ! the ratio of C and N returned to litters from microbes
real :: etaN       = 0.025   ! N loss through runoff (organic and mineral)

! Climate and vegetation types (for LM3)
real :: phen_ev1 = 0.5, phen_ev2 = 0.9 ! thresholds for evergreen/decidious
      ! differentiation (see phenology_type in cohort.F90)
real :: tg_c3_thresh = 1.5 ! threshold biomass between tree and grass for C3 plants
real :: tg_c4_thresh = 2.0 ! threshold biomass between tree and grass for C4 plants

! -------- PFT-specific parameters ----------
! c4grass  c3grass  temp-decid  tropical  evergreen  BE  BD  BN  NE  ND  G  D  T  A
integer :: pt(0:MSPECIES) = 0
!(/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0/) ! 0 for C3, 1 for C4
integer :: phenotype(0:MSPECIES)= 0
! (/0,  0,  0,  0,  1,  1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0 /) ! 0 for Deciduous, 1 for evergreen
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
real :: tauNSC(0:MSPECIES)  = 3 ! 3 ! NSC residence time,years
real :: fNSNmax(0:MSPECIES) = 5 ! 5 ! multiplier for NSNmax as sum of potential bl and br
real :: transT(0:MSPECIES)  = 3 ! Years
real :: f_cGap(0:MSPECIES)  = 0.1  ! The gaps between trees
real :: LFR_rate(0:MSPECIES)= 1.0

! Leaf parameters
real :: leaf_size(0:MSPECIES)= 0.04 !
real :: LAImax(0:MSPECIES)   = 3.5 ! maximum LAI for a tree
real :: LAI_light(0:MSPECIES)= 4.0 ! maximum LAI limited by light
real :: LMA(0:MSPECIES)      = 0.035  !  leaf mass per unit area, kg C/m2
!(/0.04,    0.04,    0.035,   0.035,   0.140,  0.032, 0.032,  0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036  /)
real :: leafLS(0:MSPECIES)   = 1.0
real :: LNbase(0:MSPECIES)   = 0.8E-3 !functional nitrogen per unit leaf area, kg N/m2
real :: CN0leafST(0:MSPECIES)= 80.0 ! CN ratio of leaf supporting tissues

! photosynthesis parameters
real :: Vmax(0:MSPECIES)= 35.0E-6 ! mol m-2 s-1
real :: m_cond(0:MSPECIES)= 7.0 !
real :: alpha_ps(0:MSPECIES)=  0.06 !
real :: Vannual(0:MSPECIES) = 1.2 ! kgC m-2 yr-1
real :: ps_wet(0:MSPECIES) = 0.3 ! wet leaf photosynthesis down-regulation: 0.3 means
        ! photosynthesis of completely wet leaf will be 30% reduction

! Wood parameters
real :: rho_wood(0:MSPECIES) = 300.0 ! kgC m-3
real :: f_taper(0:MSPECIES)  = 0.75 ! taper factor, from a cylinder to a tree

! root parameters
real :: alpha_FR(0:MSPECIES) = 1.2 ! Fine root turnover rate yr-1
!(/0.8, 0.8,0.8, 0.8, 0.8,0.8,0.8,0.8,1.0,1.0,0.6, 1.0, 0.55, 0.9, 0.55, 0.55/)
real :: rho_FR(0:MSPECIES) = 200 ! woody density, kgC m-3
real :: root_r(0:MSPECIES) = 2.9E-4
!(/1.1e-4, 1.1e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4/)
real :: root_zeta(0:MSPECIES) = 0.29 !
real :: Kw_root(0:MSPECIES)= 6.3E-8 * 1.e3 ! (kg m-2 s−1 MPa−1) ! Ref: 6.3±3.1×10−8 (m s−1 MPa−1)
! * (1000000.0/18.0)*1.e-6 ! mol /(s m2 Pa)
!Ref added by Weng, 2021-11-15
! Sutka et al. 2011 Natural Variation of Root Hydraulics in Arabidopsis Grown in Normal and Salt-Stressed Conditions.
! Plant Physiol. 155(3): 1264–1276. doi: 10.1104/pp.110.163113
! Miyamoto et al. 2001. Hydraulic conductivity of rice roots. J. Exp. Bot., 52: 1835–1846,
! doi: 10.1093/jexbot/52.362.1835

!(/1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5/)
   ! fine root membrane permeability per unit membrane area, kg/(m3 s).
   ! Root membrane permeability is "high" for the value from Siqueira et al., 2008,
   ! Water Resource Research Vol. 44, W01432, converted to mass units

!real :: rho_N_up0(0:MSPECIES) = 0.5 ! fraction of mineral N per hour
!real :: N_roots0(0:MSPECIES) = 0.3 ! kgC m-2

! Respiration rates
real :: gamma_L(0:MSPECIES)= 0.02 !
real :: gamma_LN(0:MSPECIES)= 70.5 ! 25.0  ! kgC kgN-1 yr-1
real :: gamma_SW(0:MSPECIES)= 0.08 ! 5.0e-4 ! kgC m-2 Acambium yr-1
real :: gamma_FR(0:MSPECIES)= 12.0 ! 15 !kgC kgN-1 yr-1 ! 0.6: kgC kgN-1 yr-1

! Phenology parameters
real :: tc_crit_off(0:MSPECIES)= 273.15 + 15. ! 283.16 ! OFF
real :: tc_crit_on(0:MSPECIES) = 273.15 + 10. ! 280.16 ! ON
real :: gdd_crit(0:MSPECIES)= 300. ! 280.0 !

! Reproduction prarameters
real :: AgeRepro(0:MSPECIES) = 5.0  ! year
real :: v_seed(0:MSPECIES)   = 0.1  ! fraction of allocation to wood+seeds
real :: s0_plant(0:MSPECIES) = 0.05 ! kgC, initial seedling size
real :: prob_g(0:MSPECIES)   = 1.0
real :: prob_e(0:MSPECIES)   = 1.0

! Mortality parameters
real :: r0mort_c(0:MSPECIES) = 0.01 ! yearly
real :: D0mu(0:MSPECIES)     = 2.0     ! m, Mortality curve parameter
real :: A_un(0:MSPECIES)     = 1.0     ! Multiplier for understory mortality
real :: A_sd(0:MSPECIES)     = 9.0     ! Max multiplier for seedling mortality
real :: B_sd(0:MSPECIES)     = -20.    ! Mortality sensitivity for seedlings
real :: A_D(0:MSPECIES)      = 4.0   ! Sensitivity to dbh
real :: s_hu(0:MSPECIES)     = -25.0 ! hydraulic mortality sensitivity

! Plant hydraulics parameters
real :: kx0(0:MSPECIES)      = 5.0 ! (mm/s)/(MPa/m) !132000.0 ! 6000.0   ! (m/yr-1)/(MPa/m)
real :: WTC0(0:MSPECIES)     = 1200.0  ! 2000, m /lifetime
real :: CR_Leaf(0:MSPECIES)  = 0.5 ! leaf compression ratio per MPa
real :: CR_Wood(0:MSPECIES)  = 0.2 ! Wood compression ratio per MPa
real :: w0L_max(0:MSPECIES)  = 18.0  ! leaf maximum water/carbon ratio ()
real :: w0S_max(0:MSPECIES)  = 2.0   ! stem maximum water/carbon ratio
real :: psi0_LF(0:MSPECIES)  = -3.0  ! MPa
real :: psi0_WD(0:MSPECIES)  = -3.0  ! MPa
real :: psi50_WD(0:MSPECIES) = -1.5  ! MPa !wood potential at which 50% conductivity lost, MPa
real :: Kexp_WD(0:MSPECIES)  = 3.0
real :: f_supply(0:MSPECIES) = 0.5

! C/N ratios for plant pools
real :: CNleaf0(0:MSPECIES)   = 25. ! C/N ratios for leaves
real :: CNsw0(0:MSPECIES)     = 350.0 ! C/N ratios for woody biomass
real :: CNwood0(0:MSPECIES)   = 350.0 ! C/N ratios for woody biomass
real :: CNroot0(0:MSPECIES)   = 40.0 ! C/N ratios for leaves ! Gordon & Jackson 2000
real :: CNseed0(0:MSPECIES)   = 20.0 ! C/N ratios for seeds
real :: NfixRate0(0:MSPECIES) = 0.0  ! Reference N fixation rate (0.03 kgN kg rootC-1 yr-1)
real :: NfixCost0(0:MSPECIES) = 12.0 ! FUN model, Fisher et al. 2010, GBC; Kim

namelist /vegn_parameters_nml/  diff_S0,                              &
  pt, phenotype, lifeform,                                            &
  alphaHT,alphaCA,alphaBM,thetaHT,thetaCA,thetaBM,f_taper,f_cGap,     &
  ! Leaf
  LAImax,LAI_light,LMA,Vmax,m_cond,Vannual,ps_wet,c_LLS,leaf_size,    &
  ! Wood and root
  rho_wood,rho_FR,root_r,root_zeta,Kw_root,rho_N_up0, N_roots0,       &
  ! Growth & respiration
  f_iniBSW,f_LFR_max,GR_factor,LFR_rate,tauNSC,phiRL,phiCSA,          &
  NfixRate0, NfixCost0,f_N_add,fNSNmax,transT, l_fract,               &
  retransN,gamma_L, gamma_LN, gamma_SW, gamma_FR,                     &
  ! Phenology
  gdd_crit,tc_crit_off,tc_crit_on,T0_gdd,T0_chill,                    &
  gdd_par1,gdd_par2,gdd_par3,                                         &
  ! Reproduction and Mortality
  AgeRepro,v_seed,s0_plant,prob_g,prob_e,                             &
  r0mort_c,D0mu,A_un,A_sd,B_sd,A_mort,B_mort,A_D,s_hu,                &
  ! Tisue C/N ratios
  LNbase,CN0leafST,CNleaf0,CNsw0,CNwood0,CNroot0,CNseed0,             &
  ! Plant hydraulics
  WTC0_WD,kx0_WD,psi0_WD,p50_WD,ths0_WD,f0_WD,CR0_LF,CR0_WD,          &
  TK0_leaf,kx0, WTC0, psi0_LF,psi0_osm,r_DF,m0_WTC,m0_kx,             &
  ! Soil
  FLDCAP,WILTPT,LMAmin,fsc_fine,fsc_wood,                             &
  K0SOM,K_nitrogen,rho_SON,f_M2SOM,fDON,etaN,                         &
  ! Fire model parameters, Weng, 01/13/2021
  envi_fire_prb,Ignition_G0, Ignition_W0,m0_w_fire, m0_g_fire, r_BK0, &
  f_HT0 , h0_escape, D_BK0,                    & ! for an old scheme
  phen_ev1, phen_ev2, tg_c3_thresh, tg_c4_thresh ! LM3 PFT transitions

! -------------------------------------------
! Soil parameters
! Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
  real :: GMD(n_dim_soil_types) = & ! geometric mean partice diameter, mm
  (/ 0.7, 0.4, 0.3, 0.1, 0.1, 0.07, 0.007, 0.3, 0.3 /)
  real :: GSD(n_dim_soil_types) = & ! geometric standard deviation of particle size
  (/5.0, 5.3, 7.4, 6.1, 6.1, 14.0, 15.0, 7.4, 7.4 /)
  real :: vwc_sat(n_dim_soil_types)= &
   (/ 0.380, 0.445, 0.448, 0.412, 0.414, 0.446, 0.424, 0.445, 0.445   /)
  !real :: vlc_min(n_dim_soil_types)
  real :: k_sat_ref(n_dim_soil_types)= & ! mol/(s MPa m) , hydraulic conductivity of saturated soil,
  (/ 130.8, 75.1, 53.2, 12.1, 11.1, 12.7, 1.69, 53.2, 53.2 /)
  real :: psi_sat_ref(n_dim_soil_types) = & ! Pa
  (/ -600., -790., -910., -1580., -1680., -1880., -5980., -790., -790./)
  real :: chb(n_dim_soil_types) = &         ! Soil texture parameter
  (/   3.5,   6.4,  11.0,   4.8,   6.3,   8.4,   6.3,   6.4,   6.4   /)
  real :: alphaSoil(n_dim_soil_types) = 1.0       ! *** REPLACE LATER BY alpha(layer)
  real :: heat_capacity_dry(n_dim_soil_types) = &
  (/ 1.2e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.4e6,   1.0   /)

! Soil hydraulic and heat parameter name list
namelist /soil_data_nml/ soiltype,             &
     GMD, GSD, vwc_sat,k_sat_ref, psi_sat_ref, &
     chb, alphaSoil,heat_capacity_dry

!----- Initial conditions -------------
integer, parameter :: MAX_INIT_COHORTS = 10 ! Weng, 2014-10-01
integer :: init_n_cohorts                        = MAX_INIT_COHORTS
integer :: init_cohort_species(MAX_INIT_COHORTS) = 2
real :: init_cohort_nindivs(MAX_INIT_COHORTS) = 1.0  ! initial individual density, individual/m2
real :: init_cohort_bl(MAX_INIT_COHORTS)      = 0.0  ! initial biomass of leaves, kg C/individual
real :: init_cohort_br(MAX_INIT_COHORTS)      = 0.0  ! initial biomass of fine roots, kg C/individual
real :: init_cohort_bsw(MAX_INIT_COHORTS)     = 0.05 ! initial biomass of sapwood, kg C/individual
real :: init_cohort_bHW(MAX_INIT_COHORTS)     = 0.0  ! initial biomass of heartwood, kg C/tree
real :: init_cohort_seedC(MAX_INIT_COHORTS)   = 0.0  ! initial biomass of seeds, kg C/individual
real :: init_cohort_nsc(MAX_INIT_COHORTS)     = 0.05 ! initial non-structural biomass, kg C/

! Initial soil Carbon and Nitrogen for a vegn tile, Weng 2012-10-24
real :: init_fast_soil_C  = 0.0  ! initial fast soil C, kg C/m2
real :: init_slow_soil_C  = 0.0  ! initial slow soil C, kg C/m2
real :: init_Nmineral     = 0.015  ! Mineral nitrogen pool, (kg N/m2)
real :: N_input           = 0.0008 ! annual N input to soil N pool, kgN m-2 yr-1

! Input files
character(len=80) :: filepath_in = '/Users/eweng/Documents/BiomeESS/forcingData/'
character(len=160) :: climfile = 'US-Ha1forcing.txt'
integer  :: datalines ! the total lines in forcing data file
integer  :: yr_data   ! Years of the forcing data
integer  :: days_data ! days of the forcing data

! Model run control
integer  :: model_run_years = 100
integer  :: totyears, totdays, steps_per_day ! 24 or 48
integer  :: equi_days    = 0 ! 100 * 365
real     :: step_hour    = 1.0  ! hour, Time step of forcing data, usually hourly (1.0)
real     :: dt_fast_yr   = 1.0 / (365.0 * 24.0) ! Hourly
real     :: dt_daily_yr  = 1.0/365.0 ! Daily
real     :: step_seconds = 1.0 * 3600.0

!Special test controls
logical  :: outputhourly = .True.
logical  :: outputdaily  = .True.
logical  :: do_U_shaped_mortality = .False.
logical  :: update_annualLAImax = .False.
logical  :: do_migration = .False.
logical  :: do_fire = .False.
logical  :: do_closedN_run = .True. !.False.
logical  :: do_VariedKx   = .True. ! trunk new xylem has the same kx or not
logical  :: do_VariedWTC0 = .True.
logical  :: do_WD_mort_function = .False.

! Scenarios
real     :: Sc_prcp = 1.0 ! Scenario of rainfall changes

! Model initialization name list
namelist /initial_state_nml/ &
    ! initial vegetation states
    init_n_cohorts, init_cohort_species, init_cohort_nindivs,   &
    init_cohort_bl, init_cohort_br, init_cohort_bsw,            &
    init_cohort_bHW, init_cohort_seedC, init_cohort_nsc,        &
    init_fast_soil_C, init_slow_soil_C, init_Nmineral, N_input, &
    ! Model run controls
    filepath_in,climfile, model_run_years, outputhourly, outputdaily,  &
    do_U_shaped_mortality,update_annualLAImax, do_fire, do_migration,  &
    do_closedN_run, do_VariedKx, do_variedWTC0,do_WD_mort_function,    &
    Sc_prcp

!---------------------------------
 contains
!============================ Subroutines =================================

!========================Parameter initialization =========================
!==========================================================================
subroutine initialize_soilpars(namelistfile)
   character(len=50),intent(in) :: namelistfile

  ! ---- local vars
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  integer :: nml_unit

!  Read parameters from the parameter file (namelist)
  if(read_from_parameter_file)then
     nml_unit = 999
     open(nml_unit, file=namelistfile, form='formatted', action='read', status='old')
     read (nml_unit, nml=soil_data_nml, iostat=io, end=10)
10   close (nml_unit)
     write (*, nml=soil_data_nml)
  endif

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
  soilpars%vwc_wilt = soilpars%vwc_sat &
          *(soilpars%psi_sat_ref/(psi_wilt*soilpars%alpha))**(1/soilpars%chb)
  soilpars%vwc_fc = soilpars%vwc_sat &
              *(rate_fc/(soilpars%k_sat_ref*soilpars%alpha**2))**(1/(3+2*soilpars%chb))
  soilpars%vlc_min = soilpars%vwc_sat*K_rel_min**(1/(3+2*soilpars%chb))

end subroutine initialize_soilpars

! ================================================
subroutine initialize_PFT_data(namelistfile)
  ! Initialize PFT parameters
  character(len=50),intent(in) :: namelistfile
  ! ---- local vars
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  integer :: nml_unit

  !  Read parameters from the parameter file (namelist)
  if(read_from_parameter_file)then
      nml_unit = 999
      open(nml_unit, file=namelistfile, form='formatted', action='read', status='old')
      read (nml_unit, nml=vegn_parameters_nml, iostat=io, end=10)
      10 close (nml_unit)
   endif
   write(*,nml=vegn_parameters_nml)

  ! initialize vegetation data structure
  spdata%pt         = pt
  spdata%phenotype  = phenotype
  spdata%Vmax       = Vmax
  spdata%Vannual    = Vannual
  spdata%m_cond     = m_cond
  spdata%alpha_ps = alpha_ps
  spdata%ps_wet = ps_wet
  spdata%LFR_rate = LFR_rate
  spdata%gamma_L  = gamma_L
  spdata%gamma_LN = gamma_LN
  spdata%gamma_SW = gamma_SW
  spdata%gamma_FR = gamma_FR

  spdata%rho_FR    = rho_FR
  spdata%root_r    = root_r
  spdata%root_zeta = root_zeta
  spdata%Kw_root   = Kw_root
  !  spdata%rho_N_up0 = rho_N_up0
  !  spdata%N_roots0  = N_roots0

  spdata%leaf_size   = leaf_size
  spdata%tc_crit_off = tc_crit_off
  spdata%tc_crit_on  = tc_crit_on
  spdata%gdd_crit    = gdd_crit

  ! Plant traits
  spdata%LMA        = LMA      ! leaf mass per unit area, kg C/m2
  spdata%LNbase     = LNbase   ! Basal leaf nitrogen per unit area, kg N/m2
  spdata%CN0leafST  = CN0leafST    ! Supportive tissues
  spdata%lifeform   = lifeform
  spdata%alphaHT    = alphaHT
  spdata%thetaHT    = thetaHT
  spdata%alphaCA    = alphaCA
  spdata%thetaCA    = thetaCA
  spdata%alphaBM    = alphaBM
  spdata%thetaBM    = thetaBM

  spdata%AgeRepro = AgeRepro
  spdata%v_seed   = v_seed
  spdata%s0_plant = s0_plant
  spdata%prob_g   = prob_g
  spdata%prob_e   = prob_e
  spdata%r0mort_c = r0mort_c
  spdata%D0mu     = D0mu
  spdata%A_un     = A_un
  spdata%A_sd     = A_sd
  spdata%B_sd     = B_sd
  spdata%A_D      = A_D
  spdata%s_hu     = s_hu
  spdata%rho_wood = rho_wood
  spdata%f_taper  = f_taper
  spdata%kx0      = kx0
  spdata%WTC0     = WTC0
  spdata%CR_Leaf  = CR_Leaf
  spdata%CR_Wood  = CR_Wood
  spdata%w0L_max  = w0L_max
  spdata%w0S_max  = w0S_max
  spdata%psi0_LF  = psi0_LF
  spdata%psi0_WD  = psi0_WD
  spdata%psi50_WD = psi50_WD
  spdata%Kexp_WD  = Kexp_WD

  spdata%LAImax       = LAImax
  spdata%LAImax_u     = 1.2 ! LAImax
  spdata%LAI_light    = LAI_light
  spdata%tauNSC       = tauNSC
  spdata%fNSNmax      = fNSNmax
  spdata%transT       = transT
  spdata%phiRL        = phiRL
  spdata%phiCSA       = phiCSA
  ! root urnover rate
  spdata%alpha_FR = alpha_FR

  !! Nitrogen Weng 2012-10-24
  ! spdata%CNleaf0 = CNleaf0
  spdata%CNsw0   = CNsw0
  spdata%CNwood0 = CNwood0
  spdata%CNroot0 = CNroot0
  spdata%CNseed0 = CNseed0
  spdata%NfixRate0 = NfixRate0
  spdata%NfixCost0 = NfixCost0
  spdata%f_cGap = f_cGap

  write(*,*)'  kx0,    WTC0,    CR_Wood,    psi50_WD,    psi0_WD,    Kexp_WD,    f_supply,    r0mort_c'
  do i = 0, MSPECIES
     call init_derived_species_data(spdata(i))
  enddo
  end subroutine initialize_pft_data
!------------------------------------------
 subroutine init_derived_species_data(sp)
   type(spec_data_type), intent(inout) :: sp
   ! ---- local vars ------
   integer :: i,j
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
   sp%rho_leaf = sp%LMA/sp%leafTK
   sp%CR_Leaf  = CR0_LF  * (0.02/sp%LMA)
   sp%w0L_max = rho_H2O*(1/sp%rho_leaf - 1/rho_cellwall) ! 18.0  ! leaf max. water/carbon ratio
   sp%w0S_max = rho_H2O*(1/sp%rho_wood - 1/rho_cellwall) ! 2.0   ! stem max. water/carbon ratio

   ! Wood hydraulic traits as functions of wood density, 06/30/2022, Weng
   R_WD        = sp%rho_wood/WDref0
   sp%kx0      = kx0_WD  * R_WD**(-1)  ! (mm/s)/(Mpa/m)
   sp%WTC0     = WTC0_WD * R_WD**1.5
   sp%CR_Wood  = CR0_WD  * R_WD**(-1.67)  ! Compress ratio per MPa, Santiago et al. 2018
   sp%psi50_WD = p50_WD  * R_WD**1.73 - 1.0 !- 1.09 - 3.57 * (sp%rho_wood/500.) ** 1.73
   sp%psi0_WD  = p50_WD  * R_WD**1.73 - 2.0
   sp%Kexp_WD  = ths0_WD * R_WD**(-0.5) ! rho_wood cannot be 200 for ths0_WD * R_WD**(-1)
   sp%f_supply = f0_WD   /(R_WD+1.0)

   ! Mortality rate as a function of wood density
   if(do_WD_mort_function)then
      !sp%r0mort_c = A_mort * exp(B_mort*R_WD)
      sp%r0mort_c = 0.048 - 0.024 * R_WD
   endif
   write(*,'(40(F10.4,","))')sp%kx0,sp%WTC0,sp%CR_Wood,sp%psi50_WD,sp%psi0_WD,sp%Kexp_WD,sp%f_supply,sp%r0mort_c

   ! -------- Check parameter boundaries
   if(sp%lifeform==0)then
     sp%transT = dt_daily_yr ! Grass BLmax transition in canopy layer is short
   else
     sp%transT = max(sp%transT, dt_daily_yr) !Minimum for woody plants
   endif

 end subroutine init_derived_species_data

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
   end associate
 end function

!-------------------------------------------
function BLmax2BRmax(cc) result (BRmax)
  real :: BRmax ! returned value
  type(cohort_type), intent(in) :: cc    ! cohort to update

  associate(sp=>spdata(cc%species))
    BRmax = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
  end associate
end function

 !-------------------------------------------
 function BL2Aleaf(bl,cc) result (area)
 real :: area ! returned value
 real, intent(in) :: bl      ! biomass of leaves, kg C/individual
 type(cohort_type), intent(in) :: cc    ! cohort to update

 area = bl/spdata(cc%species)%LMA
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
       cosz = max (cosz, 0.01)  ! Sun's angular is 0.01

       ! compute the solar elevation and zenth angles below
       solarelev = asin(cosz)/pi*180.0  !since asin(cos(zen))=pi/2-zen=elev
       solarzen = 90.0 - solarelev ! pi/2.d0 - solarelev
  end subroutine calc_solarzen

  ! ===== Get PFTs based on climate conditions =================
  ! Author: Elena and Sergey, LM3

  ! =============================================================================
  subroutine vegn_biogeography(vegn)
    type(vegn_tile_type), intent(inout) :: vegn

    ! ---- local vars
    integer :: i

    do i = 1, vegn%n_cohorts
       call update_species(vegn%cohorts(i), vegn%t_ann, vegn%t_cold, &
            vegn%p_ann*seconds_per_year, vegn%ncm, vegn%landuse)
    enddo
  end subroutine

  ! ============================================================================
  ! given a cohort, climatology, and land use type, determines and updates
  ! physiology type, phenology type, and species of the cohort
  subroutine update_species(c, t_ann, t_cold, p_ann, cm, landuse)
    type(cohort_type), intent(inout) :: c    ! cohort to update
    real,              intent(in) :: t_ann   ! annual-mean temperature, degK
    real,              intent(in) :: t_cold  ! average temperature of the coldest month, degK
    real,              intent(in) :: p_ann   ! annual-mean precipitation, mm/yr
    real,              intent(in) :: cm      ! number of cold months
    integer,           intent(in) :: landuse ! land use type
    ! ---- local var -----------
    integer :: spp

    c%pt    = c3c4(c,t_ann,p_ann)
    c%phenotype = phenology_type(cm)

    if(landuse == LU_CROP) c%phenotype = 0  ! crops can't be evergreen

    if(c%pt==PT_C4) then
       spp=SP_C4GRASS  ! c4 grass
    elseif(c%phenotype==1) then
       spp=SP_EVERGR   ! evergreen non-grass
    elseif(btotal(c) < tg_c3_thresh) then
       spp=SP_C3GRASS  ! c3 grass
    elseif ( t_cold > 278.16 ) then  ! ens,slm Jun 21 2003 to prohibit tropical forest in coastal cells
       spp=SP_TROPICAL ! tropical deciduous non-grass
    else
       spp=SP_TEMPDEC  ! temperate deciduous non-grass
    endif

    ! reset leaf age to zero if species are chnaged
    if (spp/=c%species) c%leaf_age = 0.0

    c%species = spp
  end subroutine

  ! ============================================================================
  function btotal(c)
    real :: btotal ! returned value
    type(cohort_type), intent(in) :: c

    btotal = c%NSC + c%seedC + c%bl + c%br + c%bsw + c%bHW
  end function

  ! ============================================================================
  function c3c4(c, temp, precip) result (pt)
    integer :: pt
    type(cohort_type), intent(in) :: c
    real,              intent(in) :: temp   ! temperatire, degK
    real,              intent(in) :: precip ! precipitation, ???

    real :: pc4

    ! Rule based on analysis of ED global output; equations from JPC, 2/02
    if(btotal(c) < tg_c4_thresh) then
      pc4=exp(-0.0421*(273.16+25.56-temp)-(0.000048*(273.16+25.5-temp)*precip))
    else
      pc4=0.0
    endif

    if(pc4>0.5) then
      pt=PT_C4
    else
      pt=PT_C3
    endif

  end function

  ! ============================================================================
  ! given current conditions, returns type of phenology.
  function phenology_type(cm) ! result(phenology_type)
    integer :: phenology_type ! returned value
    real, intent(in) :: cm ! number of cold months

    real :: pe  ! prob evergreen

    ! GCH, Rule based on analysis of ED global output; equations from JPC, 2/02
    ! GCH, Parameters updated 2/9/02 from JPC
    pe = 1.0/(1.0+((1.0/0.00144)*exp(-0.7491*cm)))

    if(pe>phen_ev1 .and. pe<phen_ev2) then
       phenology_type = PHEN_EVERGREEN ! its evergreen
    else
       phenology_type = PHEN_DECIDIOUS ! its deciduous
    endif
  end function

!================================================
end module datatypes
