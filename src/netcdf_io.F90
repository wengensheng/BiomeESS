! 07/18/2025
! Ensheng Weng, read netcdf files for BiomeE global run

!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is a very simple example which writes a 2D array of
!     sample data. To handle this in netCDF we create two shared
!     dimensions, "x" and "y", and a netCDF variable, called "data".

!     This example demonstrates the netCDF Fortran 90 API. This is part
!     of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_xy_wr.f90,v 1.7 2006/12/09 18:44:58 russ Exp $

!==============================================================
!==============================================================
module netcdf_io
  use netcdf
  use datatypes

  ! CRU NetCDF file dimensions
  integer, parameter :: NDIMS = 3
  integer, parameter :: Nlon = 720, Nlat = 360, Ntime = 1460

  public ReadNCfiles
  public unzip_gzip_file

contains

!===================================================
subroutine read_global_setting(fnml)
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
  read (nml=global_setting_nml, iostat=rc, unit=fu)
  if (rc /= 0) then
    write(*,*)'Namelist global_setting_nml error', rc
    stop
  endif
  !write(*,nml=initial_state_nml)
  close (fu)

  ! Open and read Namelist file.
  open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
  read (nml=soil_data_nml, iostat=rc, unit=fu)
  if (rc /= 0) then
    write(*,*)'Namelist soil_data_nml error', rc
    stop
  endif
  !write(*,nml=soil_data_nml)
  close (fu)

  open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
  read (nml=vegn_parameters_nml, iostat=rc, unit=fu)
  if (rc /= 0) then
    write(*,*)'Namelist vegn_parameters_nml error', rc
    stop
  endif
  !write(*,nml=vegn_parameters_nml)
  close (fu)

  open (action='read', file=fnml, status='old', iostat=rc, newunit=fu)
  read (nml=initial_state_nml, iostat=rc, unit=fu)
  if (rc /= 0) then
    write(*,*)'Namelist initial_state_nml error', rc
    stop
  endif
  !write(*,nml=initial_state_nml)
  close (fu)

end subroutine read_global_setting

!===================================================
subroutine read_PFT_parameters(fnml)
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
  read (nml=vegn_parameters_nml, iostat=rc, unit=fu)
  if (rc /= 0) then
    write(*,*)'Namelist vegn_parameters_nml error', rc
    stop
  endif
  !write(*,nml=vegn_parameters_nml)
  close (fu)

end subroutine read_PFT_parameters

!===================================================
subroutine read_initial_state(fnml)
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
  if (rc /= 0) then
    write(*,*)'Namelist initial_state_nml error', rc
    stop
  endif
  !write(*,nml=initial_state_nml)
  close (fu)

end subroutine read_initial_state

!===================================================
  subroutine ReadNCfiles (fpath,fields,yr_start, yr_end)
    implicit none
    character(len=*), intent(in) :: fpath
    character(len=*), intent(in) :: fields(:)
    integer, intent(in) :: yr_start, yr_end

    !-------- local vars -----------------
    character(len=256) :: command
    character (len =256) :: fname,fnc,fgz
    character (len =20) :: field_idx
    character (len =3) :: PFTID(9)
    character (len =4) :: yr_str

    integer, pointer :: GridMask(:,:)
    integer :: N_yrs,totL,N_vars
    integer :: iostat,i,j,k,m,iLon,iLat
    real :: dataarray(Nlon,Nlat,Ntime),timearray(Ntime)
    real :: PFTdata(144,90,9),VegFraction(144,90)

    ! Read in PFT map
    PFTID = [character(len=3) :: 'TEB','EGN','CDB','TDB','CDN','CAS','AAS','C3G','C4G']
    fnc = trim(fpath)//'BiomeE-PFTs.nc'
    write(*,*)trim(fnc)
    do i=1, 9
      call nc_read_2D(fnc,PFTID(i),144,90,PFTdata(:,:,i))
      write(*,*)"read PFT: ", PFTID(i)
    enddo
    do i =1, 144
      do j=1, 90
        VegFraction(i,j) = max(0.0,min(1.0, sum(PFTdata(i,j,:))))
      enddo
    enddo

    !fields = [character(len=5) :: 'tmp', 'pre', 'tswrf', 'spfh']
    N_vars = size(fields)
    N_yrs = yr_end - yr_start + 1
    totL = N_yrs*Ntime

    ! Open a file for identifying valid grids
    fgz = trim(fpath)//'tswrf_v12_2010.nc.gz'
    fnc = trim(fpath)//'tswrf_v12_2010.nc'
    write(*,*)trim(fgz)

    call unzip_gzip_file(trim(fgz))
    call nc_read_3D(fnc,'tswrf',Nlon,Nlat,Ntime,dataarray)
    ! Remove unziped file
    command = 'rm ' // trim(fnc)
    call execute_command_line(command)

    ! Calculate number of vegetated grids
    allocate(GridMask(LowerLon:UpperLon, LowerLat:UpperLat))
    GridMask = 0
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        i = iLon/5 + 1
        j = iLat/4 + 1
        if(dataarray(ilon,ilat,1) < 9999.0 .and. VegFraction(i,j)>0.1)then
          m = m + 1
          GridMask(ilon,ilat) = 1
        endif
      enddo
    enddo
    N_VegGrids = m
    write(*,*)"Valid grids: ",N_VegGrids

    ! Allocate data arrays
    allocate(tswrfH(totL))
    !allocate(CRUData(totL, 4, LowerLon:UpperLon, LowerLat:UpperLat))
    allocate(GridClimateData(totL, N_vars, N_VegGrids))
    allocate(GridLonLat(N_VegGrids))
    allocate(CRUgrid(N_VegGrids))

    ! Record LonLat in GridLonLat
    m = 0
    do iLon = LowerLon, UpperLon
      do iLat = LowerLat, UpperLat
        if(GridMask(ilon,ilat) > 0)then
          m = m + 1
          GridLonLat(m) = iLon * 1000 + iLat
          CRUgrid(m)%iLon = iLon
          CRUgrid(m)%iLat = iLat
          CRUgrid(m)%climate => GridClimateData(:,:,m)
          ! Assigne PFT coverage for each grid
          i = (iLon-1)/5 + 1
          j = (iLat-1)/4 + 1
          do k = 1, 9
            CRUgrid(m)%fPFT(k) = Min(1.0, Max(0.0, PFTdata(i,j,k)))
          enddo
        endif
      enddo
    enddo

    ! Read in all data
    do j= 1, 4
      field_idx = fields(j) ! 'tmp'
      do i =1, N_yrs
        write(yr_str, '(I4)') yr_start + i - 1
        if(field_idx == 'tswrf')then
          fnc = trim(field_idx)//'_v12_'//trim(yr_str)//'.nc'
        else
          fnc = 'crujra.v2.4.5d.'//trim(field_idx)//'.'//trim(yr_str)//'.365d.noc.nc'
        endif
        fgz = trim(fnc)//'.gz'

        ! Unzip a file
        fgz = trim(fpath)//trim(fgz)
        call unzip_gzip_file(trim(fgz))

        fnc = trim(fpath)//trim(fnc)

        write(*,*)'read nc file:', fnc
        call nc_read_3D(fnc,trim(field_idx),Nlon,Nlat,Ntime,dataarray)
        m = 0
        do iLon = LowerLon, UpperLon
          do iLat = LowerLat, UpperLat
            !CRUData((i-1)*Ntime+1:i*Ntime,j,iLon,iLat) = dataarray(iLon,iLat,:)
            if(GridMask(iLon,iLat) > 0)then
              m = m + 1
              GridClimateData((i-1)*Ntime+1:i*Ntime,j,m)  = dataarray(iLon,iLat,:)
            endif
          enddo
        enddo

        if(field_idx == 'tswrf')then
          call nc_read_1D(fnc,'time',Ntime,timearray)
          tswrfH((i-1)*Ntime+1: i*Ntime) = timearray
          !write(*,*)"Time Data"
          !write(*,'(360(E15.4,","))')  timearray(1:360)
        endif

        ! Remove unziped file
        command = 'rm ' // trim(fnc)
        call execute_command_line(command)
      enddo
    enddo
    deallocate(GridMask)
  end subroutine ReadNCfiles

!=============================================================================
subroutine Assign_global_PFT_parameters()

  ! -------- PFT-specific parameters ----------
  !                     'TEB','EGN','CDB','TDB','CDN','CAS','AAS','C3G','C4G'
   pt(1:N_PFTs)       = (/0,   0,    0,    0,    0,    0,    0,    0,    1/) ! 0 for C3, 1 for C4
   phenotype(1:N_PFTs)= (/1,   1,    0,    0,    0,    1,    1,    0,    0/) ! 0 for Deciduous, 1 for evergreen
   lifeform(1:N_PFTs) = (/1,   1,    1,    1,    1,    1,    1,    0,    0/) ! life form of PFTs: 0 for grasses, 1 for trees

  ! Allometry and whole tree parameters
   alphaHT(1:N_PFTs) = [ 35., 30., 30., 35., 30., 20., 20., 10., 10. ]
   alphaCA(1:N_PFTs) = [ 120.,120.,120.,120.,120.,150.,150., 60., 60. ]
   alphaBM(1:N_PFTs) = 5200.0
   thetaHT(1:N_PFTs) = 0.5 !
   thetaCA(1:N_PFTs) = 1.5
   thetaBM(1:N_PFTs) = 2.5
   phiRL(1:N_PFTs)   = 3.5 ! ratio of fine root area to leaf area
   phiCSA(1:N_PFTs)  = 0.25E-4 ! ratio of sapwood area to leaf area
   tauNSC(1:N_PFTs)  = [ 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 3.0, 3.0 ] ! NSC residence time,years
   fNSNmax(1:N_PFTs) = 5 ! 5 ! multiplier for NSNmax as sum of potential bl and br
   transT(1:N_PFTs)  = 3 ! Years
   f_cGap(1:N_PFTs)  = 0.1  ! The gaps between trees
   LFR_rate(1:N_PFTs)= 1.0/21.0

  ! Leaf parameters
  !                     'TEB','EGN','CDB','TDB','CDN','CAS','AAS','C3G','C4G'
   LAImax(1:N_PFTs)   = [ 4.8, 4.8, 4.5, 4.5, 4.0, 3.0, 3.0, 2.5, 2.5 ]    ! maximum LAI for a tree
   LAI_light(1:N_PFTs)= 5.0    ! maximum LAI limited by light
   LMA(1:N_PFTs)      = [ 7.0e-2, 14.e-2, 2.5e-2, 3.0e-2,3.0e-2, 2.5e-2, 3.0e-2, 2.5e-2, 2.5e-2 ]  ! leaf mass per unit area, kg C/m2
   leafLS(1:N_PFTs)   = 1.0
   LNbase(1:N_PFTs)   = 1.3E-3 !functional nitrogen per unit leaf area, kg N/m2, 1.1E-3 for Acer, 1.5E-3 for Populus
   CN0leafST(1:N_PFTs)= 40.0 ! 80.0 ! CN ratio of leaf supporting tissues
   leaf_size(1:N_PFTs)= 0.04 !

  ! photosynthesis parameters
   Vmax(1:N_PFTs)= 35.0E-6 ! mol m-2 s-1
   m_cond(1:N_PFTs)= 9.0 ! 7.0 !
   alpha_ps(1:N_PFTs)=  0.06 !
   Vannual(1:N_PFTs) = 1.2 ! kgC m-2 yr-1
   ps_wet(1:N_PFTs) = 0.3 ! wet leaf photosynthesis down-regulation: 0.3 means
          ! photosynthesis of completely wet leaf will be 30% reduction

  ! Wood parameters
   rho_wood(1:N_PFTs) = [360.,300.,350.,250.,300.,400.,400.,90.,90.] ! kgC m-3
   f_taper(1:N_PFTs)  = 0.75 ! taper factor, from a cylinder to a tree

  ! root parameters
   alpha_FR(1:N_PFTs) = 1.2 ! Fine root turnover rate yr-1
  !(/0.8, 0.8,0.8, 0.8, 0.8,0.8,0.8,0.8,1.0,1.0,0.6, 1.0, 0.55, 0.9, 0.55, 0.55/)
   rho_FR(1:N_PFTs) = 200 ! woody density, kgC m-3
   root_r(1:N_PFTs) = 2.9E-4
  !(/1.1e-4, 1.1e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4/)
   root_zeta(1:N_PFTs) = 0.6 ! 0.29 !
   root_perm(1:N_PFTs) = 0.5 ! kg H2O m-2 hour-1, defined by Weng
   Kw_root(1:N_PFTs)   = 6.3E-8 * 1.e3 ! (kg m-2 s−1 MPa−1) ! Ref: 6.3±3.1×10−8 (m s−1 MPa−1)

  ! Respiration rates
   gamma_L(1:N_PFTs)= 0.02 !
   gamma_LN(1:N_PFTs)= 70.5 ! 25.0  ! kgC kgN-1 yr-1
   gamma_SW(1:N_PFTs)= 0.02 ! 0.08 ! kgC m-2 Acambium yr-1
   gamma_FR(1:N_PFTs)= 0.6 ! 12 !kgC kgN-1 yr-1 ! 0.6: kgC kgN-1 yr-1

  ! Phenology parameters

   tc_crit_off(1:N_PFTs)= [ 15.0, -80.0, 15.0, 15.0, 15.0, 15.0, 15.0, 5.0, 5.0 ] ! 283.16 ! OFF ! C for convenience
   tc_crit_on(1:N_PFTs) = 10. ! 280.16 ! ON  ! C for convenience
   gdd_crit(1:N_PFTs)= 300. ! 280.0 !
   AWD_crit(1:N_PFTs)= 0.7  ! Critical plant water availability factor (0~1)
   betaON(1:N_PFTs)  = [ 0.0, 0.0, 0.0, 0.4, 0.2, 0.4, 0.4, 0.4, 0.4 ]  ! Critical soil moisture for phenology ON
   betaOFF(1:N_PFTs) = [ 0.0, 0.0, 0.0, 0.2, 0.0, 0.1, 0.1, 0.2, 0.2 ]  ! Critical soil moisture for phenology OFF
   gdd_par1(1:N_PFTs) = [ 20., 0.0, 50., 20., 50., 30., 30., 30., 30. ]   !50.d0   ! -68.d0
   gdd_par2(1:N_PFTs) = [ 200., 0.0, 800., 200., 800., 800., 800., 800., 800. ] ! 650.d0  !800.d0  ! 638.d0
   gdd_par3(1:N_PFTs) = -0.02 ! -0.01d0
  ! Reproduction prarameters
   AgeRepro(1:N_PFTs) = 5.0  ! year
   v_seed(1:N_PFTs)   = 0.1  ! fraction of allocation to wood+seeds
   s0_plant(1:N_PFTs) = [ 0.1,   0.1,  0.1,  0.1,  0.1, 0.05, 0.05, 0.005, 0.005 ] ! kgC, initial seedling size
   prob_g(1:N_PFTs)   = 1.0
   prob_e(1:N_PFTs)   = 1.0

  ! Mortality parameters
   r0mort_c(1:N_PFTs) = [ .025, .025, .025, .025, .025, .02, .02, .02, .02 ] ! 0.01 ! yearly ! 0.012 for Acer, 0.0274 for Populus
   D0mu(1:N_PFTs)     = [ 1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.0, 0.0 ]     ! m, Mortality curve parameter
   A_un(1:N_PFTs)     = 3.0     ! Multiplier for understory mortality
   A_sd(1:N_PFTs)     = [ 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 19.d0, 19.d0, 0.d0, 0.d0 ]     ! Max multiplier for seedling mortality
   B_sd(1:N_PFTs)     = [ -25., -25., -25., -25., -25., -25., -25., -60., -60. ]    ! Mortality sensitivity for seedlings
   A_DBH(1:N_PFTs)    = 4.0     ! Max multiplier for DBH-based mortality
   B_DBH(1:N_PFTs)    = 0.125   ! 0.25   ! Size-based Mortality sensitivity, m
   s_hu(1:N_PFTs)     = -25.0   ! hydraulic mortality sensitivity
   W_mu0(1:N_PFTs)    = 1.0     ! Jeremy's half-mortality transp deficit, high:0.5, low: 0.75, No effects: 2.5

  ! Plant hydraulics parameters
   kx0(1:N_PFTs)      = 5.0 ! (mm/s)/(MPa/m) !132000.0 ! 6000.0   ! (m/yr-1)/(MPa/m)
   WTC0(1:N_PFTs)     = 1200.0  ! 2000, m /lifetime
   CR_Leaf(1:N_PFTs)  = 0.5 ! leaf compression ratio per MPa
   CR_Wood(1:N_PFTs)  = 0.2 ! Wood compression ratio per MPa
   psi0_LF(1:N_PFTs)  = -3.0  ! MPa
   psi0_WD(1:N_PFTs)  = -3.0  ! MPa
   psi50_WD(1:N_PFTs) = -1.5  ! MPa !wood potential at which 50% conductivity lost, MPa
   Kexp_WD(1:N_PFTs)  = 3.0
   f_supply(1:N_PFTs) = 0.5
   f_plc(1:N_PFTs)    = 0.05  ! fraction of WTC loss due to low water potential (per day)

  ! C/N ratios for plant pools
   CNleaf0(1:N_PFTs)   = 25. ! C/N ratios for leaves
   CNsw0(1:N_PFTs)     = 350.0 ! C/N ratios for woody biomass
   CNwood0(1:N_PFTs)   = 350.0 ! C/N ratios for woody biomass
   CNroot0(1:N_PFTs)   = 40.0 ! C/N ratios for leaves ! Gordon & Jackson 2000
   CNseed0(1:N_PFTs)   = 20.0 ! C/N ratios for seeds
   NfixRate0(1:N_PFTs) = 0.0  ! Reference N fixation rate (0.03 kgN kg rootC-1 yr-1)
   NfixCost0(1:N_PFTs) = 12.0 ! FUN model, Fisher et al. 2010, GBC; Kim

end subroutine Assign_global_PFT_parameters

!=============================================================================
subroutine Setup_Grid_Initial_States(CRUgrid)
  type(grid_initial_type), intent(in) :: CRUgrid

  !--------- local vars ------------
  integer :: GridPFTs(N_PFTs)
  integer :: i, totPFT
  real :: f_min = 0.01 ! coverage fraction threshold

  ! Sorting PFT numbers according to fPFT
  call rank_descending(CRUgrid%fPFT,GridPFTs)

  ! Find out PFTs in this grid
  totPFT = Max(1, COUNT(CRUgrid%fPFT > f_min))

  do i=1, totPFT ! init_n_cohorts
    init_cohort_species(i) = GridPFTs(i)
    init_cohort_nindivs(i) = 0.2  ! initial individual density, individual/m2
    init_cohort_bl(i)      = 0.0  ! initial biomass of leaves, kg C/individual
    init_cohort_br(i)      = 0.0  ! initial biomass of fine roots, kg C/individual
    init_cohort_bsw(i)     = 0.3  ! initial biomass of sapwood, kg C/individual
    init_cohort_bHW(i)     = 0.0  ! initial biomass of heartwood, kg C/tree
    init_cohort_seedC(i)   = 0.0  ! initial biomass of seeds, kg C/individual
    init_cohort_nsc(i)     = 0.3  ! initial non-structural biomass, kg C/
  enddo

 ! Initial soil Carbon and Nitrogen for a vegn tile, Weng 2012-10-24
  init_fast_soil_C  = 0.5  ! initial fast soil C, kg C/m2
  init_slow_soil_C  = 20.0  ! initial slow soil C, kg C/m2
  init_Nmineral     = 0.15  ! Mineral nitrogen pool, (kg N/m2)
  N_input           = 0.0 !0.0008 ! annual N input to soil N pool, kgN m-2 yr-1

end subroutine Setup_Grid_Initial_States

!=============================================================================
subroutine CRU_Interpolation(CRUgrid,steps_per_hour,forcingData)
  implicit none
  type(grid_initial_type), pointer, intent(in) :: CRUgrid
  real, intent(in) :: steps_per_hour
  type(climate_data_type), pointer :: forcingData(:) ! output

  !---------- local variables ------------------
  real, pointer :: GridData(:,:)
  integer :: iLon, iLat ! Column and Lines (started from -179.75 and -89.75)
  type(climate_data_type), pointer :: climateData(:) ! will be linked to forcingData
  character(len=120) :: fname ! For testing output
  character(len=6) :: LonLat  ! Used in output file name

  real, allocatable :: fdSW(:)
  real, allocatable :: timecols(:,:)
  real, allocatable :: hourly_data(:,:)
  real    :: ShiftData(12,10) ! To make the first hour is hour 0
  real    :: Lati, Longi
  real    :: steps_in_6H ! Temporary variable, steps interpolated
  real    :: td,cosz,solarelev,solarzen,r_light
  integer :: year0, year1 ! Start and end year
  integer :: yr,doy,iday,ihour,iyr
  integer :: ndays,nyear,totalL
  integer :: m,i,j,k,steps_per_day
  integer :: totyr, Nlines
  logical :: WriteSample = .False.

  ! Assigne CRUgrid data to local variables
  GridData => CRUgrid%climate
  iLon = CRUgrid%iLon
  iLat = CRUgrid%iLat

  ! Latitude and Longitude of this grid
  Longi = (iLon-1)*0.5 - 179.75
  Lati = (iLat-1) *0.5 - 89.75

  ! Data lines
  Nlines = SIZE(tswrfH)
  steps_per_day = int(steps_per_hour * 24)
  totalL = Nlines * int(6 * steps_per_hour) ! for hourly data
  totyr = int(Nlines/(365*4))
  year0 = int(tswrfH(1)/(365*24)) + 1850
  year1 = year0 + totyr - 1

  ! Allocate allocatable variables
  allocate(fdSW(Nlines))
  allocate(timecols(Nlines,3)) ! year, doy, hour
  allocate(hourly_data(totalL,10)) ! Interpolated data (hourly)
  timecols(:,:) = 0
  do m=1, Nlines
    timecols(m,1) = int(tswrfH(m)/(365*24)) + 1850.0       ! Year
    timecols(m,2) = MODULO(int(tswrfH(m)/24),365) + 1.0    ! Day of the year
    timecols(m,3) = MODULO(tswrfH(m), 24.0) + int(Longi/15.0) ! Local time
    timecols(m,3) = MODULO(timecols(m,3), 24.0) ! 0~23
  enddo

  ! Fraction of solar radiation
  do m=1, Nlines
    td = timecols(m,2) + timecols(m,3)/24.0
    call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
    if(cosz>0.001)then
      fdSW(m) = max(min(1.0,GridData(m,3)/(cosz*1380.0)),0.0)
    else
      fdSW(m) = 0.0
    endif
  enddo

  ! ------------- Data interpolation ------------------
  m = 0
  do i = 1, Nlines -1
    do j=0, 5
      m=m+1
      steps_in_6H = 6.0*steps_per_hour
      hourly_data(m,4) = GridData(i,1) + (GridData(i+1,1)-GridData(i,1))*j/steps_in_6H ! Tmp, K
      hourly_data(m,5) = GridData(i,2)/(6.0 * 3600.0)                  ! Precipitation, mm/second
      hourly_data(m,6) = GridData(i,3) + (GridData(i+1,3)-GridData(i,3))*j/steps_in_6H ! Shortwave radiation,need to calculate zenith angle
      hourly_data(m,7) = GridData(i,4) + (GridData(i+1,4)-GridData(i,4))*j/steps_in_6H ! spfh
      hourly_data(m,8) = 101325.0 ! Pa                                       ! pres
      hourly_data(m,9) = 1.2        ! Wind speed m/s
      hourly_data(m,10)= fdSW(i) + (fdSW(i+1)-fdSW(i))*j/steps_in_6H          ! fraction of Shortwave radiation
    enddo
  enddo
  ! Last 6 hours' data
  do j=0, 5
    hourly_data(totalL-j,4) = GridData(Nlines,1) + (GridData(1,1)-GridData(Nlines,1))*j/steps_in_6H ! Tmp, K
    hourly_data(totalL-j,5) = GridData(Nlines,2)/(6.0 * 3600.0)                    ! Precipitation, mm/second
    hourly_data(totalL-j,6) = GridData(Nlines,3) + (GridData(1,3)-GridData(Nlines,3))*j/steps_in_6H ! Shortwave radiation,need to calculate zenith angle
    hourly_data(totalL-j,7) = GridData(Nlines,4) + (GridData(1,4)-GridData(Nlines,4))*j/steps_in_6H ! spfh
    hourly_data(totalL-j,8) = 101325.0 ! Pa
    hourly_data(totalL-j,9)   = 1.2 ! m/s, wind spped
    hourly_data(totalL-j,10)  = fdSW(Nlines) + (fdSW(1)-fdSW(Nlines))*j/steps_in_6H
  enddo

  ! Shift data to fit a whole day
  if(timecols(1,3)>0.9999 .and. timecols(1,3)<11.9999) then ! Move down
    m = max(1,int(timecols(1,3)))
    ShiftData(1:m,:) = hourly_data(totalL-m+1:totalL,:)
    hourly_data(m+1:totalL,:) = hourly_data(1:totalL-m,:)
    hourly_data(1:m,:) = ShiftData(1:m,:)
  elseif (timecols(1,3)>=11.9999) then ! Move up
    m = 24 - max(12,int(timecols(1,3)))
    ShiftData(1:m,:) = hourly_data(1:m,:)
    hourly_data(1:totalL-m,:) = hourly_data(m+1:totalL,:)
    hourly_data(totalL-m+1:totalL,:) = ShiftData(1:m,:)
  endif

  !Time columns
  m=0
  do iyr = year0, year1
    do iday=1, 365
      do ihour = 1, steps_per_day ! 24
        m =  m + 1
        hourly_data(m,1) = iyr
        hourly_data(m,2) = iday
        hourly_data(m,3) = 24.0 * (ihour - 1)/steps_per_day
      enddo
    enddo
  enddo

  ! -------------- Put the data into forcing -------------------
  allocate(climateData(totalL))
  do i=1,totalL
     td = hourly_data(i,2) + hourly_data(i,3)/24.0
     call calc_solarzen(td,Lati,cosz,solarelev,solarzen)
     if(cosz>0.001)then
       climateData(i)%radiation = 1380.0 * cosZ * hourly_data(i,10)      ! W/m2
       climateData(i)%PAR       = climateData(i)%radiation * 2.0 ! umol/m2/s
     else
       climateData(i)%radiation = 0.0 ! W/m2
       climateData(i)%PAR       = 0.0 ! umol/m2/s
     endif
     !'tmp','pre','tswrf','spfh','pres','windU'
     climateData(i)%year  = int(hourly_data(i,1))     ! Year
     climateData(i)%doy   = int(hourly_data(i,2))     ! day of the year
     climateData(i)%hod   = hourly_data(i,3)          ! hour of the day
     climateData(i)%Tair  = hourly_data(i,4)          ! air temperature, K
     climateData(i)%Tsoil = (climateData(i)%Tair - 273.16) * 0.8 + 273.16 ! soil temperature, C
     climateData(i)%rain  = hourly_data(i,5)          ! kgH2O m-2 second-1
     climateData(i)%P_air = hourly_data(i,8)          ! pa
     climateData(i)%windU = hourly_data(i,9)          ! wind velocity (m s-1)
     climateData(i)%RH    = hourly_data(i,7) / mol_h2o * mol_air * &
                            climateData(i)%P_air/  &
                            esat(climateData(i)%Tair-273.16) ! relative humidity (0.xx)
     climateData(i)%CO2   = CO2_Hist(climateData(i)%year-1700)       ! ppm
     climateData(i)%eCO2  = climateData(i)%CO2 + 200.       ! ppm
     climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
  enddo
  forcingData => climateData
  datalines = totalL
  days_data = totyr * 365
  yr_data   = totyr
  steps_per_day = int(24.0/step_hour)
  dt_fast_yr    = step_hour/(365.0 * 24.0)
  step_seconds  = step_hour*3600.0

  ! Write out a sample file
  if(WriteSample)then
    write(LonLat, '(I6)') GridID
    fname = 'CRU_'//trim(LonLat)//'_forcing.csv'
    open(15,file=trim(fname))
    write(15,*)"YEAR,DOY,HOUR,PAR,Swdown,Tair,Tsoil,RH,RAIN,WIND,PRESSURE,aCO2,eCO2"
    do i=1,totalL
        write(15,'(2(I4,","),1(f8.2,","),30(E15.6,","))') &
          climateData(i)%year, climateData(i)%doy,      climateData(i)%hod,  &
          climateData(i)%PAR,  climateData(i)%radiation,climateData(i)%Tair, &
          climateData(i)%Tsoil,climateData(i)%RH,       climateData(i)%rain, &
          climateData(i)%windU,climateData(i)%P_air, &
          climateData(i)%CO2,  climateData(i)%eCO2
    enddo
    close(15)
    WriteSample = .False.
  endif

  !Release memory
  deallocate(hourly_data)
  deallocate(fdSW)
  deallocate(timecols)

end subroutine CRU_Interpolation

!==============================================================
  subroutine unzip_gzip_file(filename_gz)
    character (len = *), intent(in) :: filename_gz
    !----------local vars
    character(len=256) :: command
    integer :: iostat

    ! Construct the gunzip command. The -k option keeps the original .gz file.
    command = 'gunzip -k ' // trim(filename_gz)
    call execute_command_line(command, exitstat=iostat)
    if (iostat == 0) then
      print *, 'Successfully unzipped ', trim(filename_gz)
    else
      print *, 'Error unzipping ', trim(filename_gz), ' (Exit status: ', iostat, ')'
    end if
  end subroutine unzip_gzip_file

!==============================================
  subroutine nc_read_3D(FILE_NAME,field_idx,NX,NY,Ntime,DA)
    ! This is the name of the data file we will create.
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: NX, NY, Ntime
    real, intent(inout) :: DA(:,:,:)

    !----- Local vars ----------------
    integer :: ncid, varid  ! IDs were created with netCDF files

    ! Open the file with NF90_NOWRITE as read-only access
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
    call check( nf90_get_var(ncid, varid, DA) )
    call check( nf90_close(ncid) )
    print *, 'Read file: ncid=',ncid, 'varid=',varid

  end subroutine nc_read_3D

  !==============================================
  subroutine nc_read_2D(FILE_NAME,field_idx,NX,NY,DA)
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: NX, NY
    real, intent(inout) :: DA(:,:)

    !----- Local vars ----------------
    integer :: ncid, varid ! IDs were created with netCDF files

    ! Open the file with NF90_NOWRITE as read-only access
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) ) ! Get the varid of the data variable
    call check( nf90_get_var(ncid, varid, DA) )  ! Read the data.
    call check( nf90_close(ncid) ) ! Close the file, freeing all resources.
    print *, 'Read file: ncid=',ncid, 'varid=',varid
  end subroutine nc_read_2D
  !==============================================
  subroutine nc_read_1D(FILE_NAME,field_idx,Ntime,DA)
    ! This is the name of the data file we will read
    character (len = *), intent(in) :: FILE_NAME,field_idx
    integer, intent(in) :: Ntime
    real, intent(inout) :: DA(:)

    !----- Local vars ----------------
    integer :: ncid, varid ! IDs were created with netCDF files

    ! Open the file with NF90_NOWRITE as read-only access
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
    call check( nf90_get_var(ncid, varid, DA) )
    call check( nf90_close(ncid) )
    print *, 'Read file: ncid=',ncid, 'varid=',varid
  end subroutine nc_read_1D

!===================================================
  subroutine nc_write(FILE_NAME,NDIMS,NX,NY)
    ! This is the name of the data file we will create.
    character (len = *), intent(in) :: FILE_NAME
    integer, intent(in) :: NDIMS, NX, NY ! , Ntime

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    integer :: ncid, varid
    integer :: x_dimid, y_dimid

    ! This is the data array we will write. It will just be filled with
    ! a progression of integers for this example.
    integer,allocatable :: data_out(:,:) !(NY, NX)
    integer,allocatable :: dimids(:)

    ! Loop indexes, and error handling.
    integer :: x, y

    ! Create some pretend data. If this wasn't an example program, we
    ! would have some real data to write, for example, model output.
    allocate(dimids(NDIMS))
    allocate(data_out(NY,NX))
    do x = 1, NX
       do y = 1, NY
          data_out(y, x) = (x - 1) * NY + (y - 1)
       end do
    end do

    ! Always check the return code of every netCDF function call. In
    ! this example program, wrapping netCDF calls with "call check()"
    ! makes sure that any return which is not equal to nf90_noerr (0)
    ! will print a netCDF error message and exit.

    ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
    ! overwrite this file, if it already exists.
    call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )

    ! Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, "x", NX, x_dimid) )
    call check( nf90_def_dim(ncid, "y", NY, y_dimid) )

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids =  (/ y_dimid, x_dimid /)

    ! Define the variable. The type of the variable in this case is
    ! NF90_INT (4-byte integer).
    call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )

    ! End define mode. This tells netCDF we are done defining metadata.
    call check( nf90_enddef(ncid) )

    ! Write the pretend data to the file. Although netCDF supports
    ! reading and writing subsets of data, in this case we write all the
    ! data in one operation.
    call check( nf90_put_var(ncid, varid, data_out) )

    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    call check( nf90_close(ncid) )
    write(*,*)data_out
    print *, "*** SUCCESS writing ncfile simple_xy.nc! "

    deallocate(data_out)
    deallocate(dimids)

  end subroutine nc_write

!===================================================
  subroutine nc_read(FILE_NAME,NDIMS,NX,NY)
    ! This is the name of the data file we will create.
    character (len = *), intent(in) :: FILE_NAME
    integer, intent(in) :: NDIMS, NX, NY

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    integer :: ncid, varid
    ! We are reading 2D data, a 6 x 12 grid.
    integer,allocatable :: data_in(:,:) !(NY, NX)

    ! Loop indexes, and error handling.
    integer :: x, y

    allocate(data_in(NY,NX))
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
    ! the file.
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
    print *, 'ncid=',ncid
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, "data", varid) )
    print *, 'varid=',varid
    ! Read the data.
    call check( nf90_get_var(ncid, varid, data_in) )

    ! Check the data.
    do x = 1, NX
       do y = 1, NY
          print *, "data_in(", y, ", ", x, ") = ", data_in(y, x)
       end do
    end do

    ! Close the file, freeing all resources.
    call check( nf90_close(ncid) )

    print *,"*** SUCCESS reading ncfile ", FILE_NAME, "! "

    deallocate(data_in)
  end subroutine nc_read

!=================================================================
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end module
!=================================================================
