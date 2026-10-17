!============================================================================
! restart_mod.F90 -- checkpoint/restart I/O for BiomeE-Global
!
! Usage:
!   call write_restart(land)   -- at end of run when do_restart_write=.True.
!   call read_restart(land)    -- at start of run when do_restart_read=.True.
!
! Binary format: Fortran unformatted sequential (default record-length units).
! The file starts with a header (magic number, format version, array sizes) that is
! checked on read. The cohort records are written/read by write_cohort/read_cohort,
! which are also used for the initial cohorts (initialCC); the READ and WRITE record
! order is identical. Animal cohorts keep their state variables (aft, nindivs, age,
! body_C); their fluxes and annual accumulators restart from zero.
!
! For global runs where GridID > 0 the file name is constructed by inserting
! the GridID before the extension, e.g.:
!   restart_file = './output/BiomeE_restart.bin'
!   -> actual path: './output/BiomeE_restart_205249.bin'
!============================================================================

module restart_mod
  use datatypes
  use model_utils, only: ccNSCmax, ccNSNmax
  implicit none
  private

  public :: write_restart, read_restart

  integer, parameter :: RESTART_MAGIC   = 20260617 ! identifies a BiomeE restart file
  integer, parameter :: RESTART_VERSION = 3        ! 2: header, ccID, facuC, HarvYrs, initialCC, animals; 3: thetaS, soilwater

contains

!----------------------------------------------------------------------------
! Helper: build the actual file path (appends GridID when > 0)
!----------------------------------------------------------------------------
  function restart_path() result(path)
    character(len=512) :: path
    character(len=20)  :: gid_str
    integer :: dot_pos, n

    ! GridID default is 999999 (site run). Only append it for real global grids.
    if (GridID /= 999999) then
      ! Find the last '.' in restart_file to insert GridID before extension
      n = len_trim(restart_file)
      dot_pos = 0
      do while (n > 0)
        if (restart_file(n:n) == '.') then
          dot_pos = n
          exit
        endif
        n = n - 1
      enddo
      write(gid_str, '("_", I6.6)') GridID
      if (dot_pos > 0) then
        path = trim(restart_file(1:dot_pos-1)) // trim(gid_str) // &
               trim(restart_file(dot_pos:len_trim(restart_file)))
      else
        path = trim(restart_file) // trim(gid_str) // '.bin'
      endif
    else
      path = trim(restart_file)
    endif
  end function restart_path

!============================================================================
  subroutine write_cohort(iunit, cc)
!============================================================================
    integer, intent(in) :: iunit
    type(cohort_type), intent(in) :: cc


  !-- integer fields --
  write(iunit) cc%species
  write(iunit) cc%ccID
  write(iunit) cc%layer
  write(iunit) cc%status
  write(iunit) cc%Ngd
  write(iunit) cc%Ndm
  write(iunit) cc%Ncd
  write(iunit) cc%Nrings
  write(iunit) cc%phenotype
  write(iunit) cc%pt

  !-- structure real scalars --
  write(iunit) cc%nindivs
  write(iunit) cc%layerfrac
  write(iunit) cc%topyear
  write(iunit) cc%age
  write(iunit) cc%leafage
  write(iunit) cc%mu
  write(iunit) cc%dbh
  write(iunit) cc%height
  write(iunit) cc%Acrown
  write(iunit) cc%Aleaf
  write(iunit) cc%lai
  write(iunit) cc%D_bark
  write(iunit) cc%bl_max
  write(iunit) cc%br_max
  write(iunit) cc%CSAsw
  write(iunit) cc%DBH_ys
  write(iunit) cc%extinct
  write(iunit) cc%w_scale

  !-- phenology --
  write(iunit) cc%GDD
  write(iunit) cc%ALT
  write(iunit) cc%AWD
  write(iunit) cc%Tc_OFF
  write(iunit) cc%GDD_ON

  !-- C pools --
  write(iunit) cc%bl
  write(iunit) cc%br
  write(iunit) cc%bsw
  write(iunit) cc%bHW
  write(iunit) cc%seedC
  write(iunit) cc%nsc

  !-- N pools --
  write(iunit) cc%NSN
  write(iunit) cc%leafN
  write(iunit) cc%swN
  write(iunit) cc%hwN
  write(iunit) cc%rootN
  write(iunit) cc%seedN
  write(iunit) cc%facuC

  !-- root / water supply --
  write(iunit) cc%K_r
  write(iunit) cc%r_r
  write(iunit) cc%root_zeta
  write(iunit) cc%W_supply

  !-- root arrays --
  write(iunit) cc%root_length(1:soil_L)
  write(iunit) cc%ArootL(1:soil_L)
  write(iunit) cc%WupL(1:soil_L)
  write(iunit) cc%uptake_frac(1:soil_L)

  !-- hydraulics scalars --
  write(iunit) cc%psi_s0
  write(iunit) cc%psi_leaf
  write(iunit) cc%psi_stem
  write(iunit) cc%H_leaf
  write(iunit) cc%H_stem
  write(iunit) cc%W_lf
  write(iunit) cc%W_sw
  write(iunit) cc%W_hw
  write(iunit) cc%Wmax_l
  write(iunit) cc%Wmax_s
  write(iunit) cc%Wmin_l
  write(iunit) cc%Wmin_s
  write(iunit) cc%Q_stem
  write(iunit) cc%Q_leaf
  write(iunit) cc%Ktrunk
  write(iunit) cc%Asap
  write(iunit) cc%Atrunk
  write(iunit) cc%treeHU
  write(iunit) cc%treeW0
  write(iunit) cc%An_op
  write(iunit) cc%An_cl

  !-- xylem ring arrays --
  write(iunit) cc%Kx(1:Ysw_max)
  write(iunit) cc%WTC0(1:Ysw_max)
  write(iunit) cc%accH(1:Ysw_max)
  write(iunit) cc%plcH(1:Ysw_max)
  write(iunit) cc%farea(1:Ysw_max)
  write(iunit) cc%Rring(1:Ysw_max)
  write(iunit) cc%Lring(1:Ysw_max)
  write(iunit) cc%Aring(1:Ysw_max)
  write(iunit) cc%Kring(1:Ysw_max)

  end subroutine write_cohort

!============================================================================
  subroutine read_cohort(iunit, cc)
!============================================================================
    integer, intent(in) :: iunit
    type(cohort_type), intent(inout) :: cc


  !-- integer fields --
  read(iunit, err=900, end=900) cc%species
  read(iunit, err=900, end=900) cc%ccID
  read(iunit, err=900, end=900) cc%layer
  read(iunit, err=900, end=900) cc%status
  read(iunit, err=900, end=900) cc%Ngd
  read(iunit, err=900, end=900) cc%Ndm
  read(iunit, err=900, end=900) cc%Ncd
  read(iunit, err=900, end=900) cc%Nrings
  read(iunit, err=900, end=900) cc%phenotype
  read(iunit, err=900, end=900) cc%pt

  !-- structure real scalars --
  read(iunit, err=900, end=900) cc%nindivs
  read(iunit, err=900, end=900) cc%layerfrac
  read(iunit, err=900, end=900) cc%topyear
  read(iunit, err=900, end=900) cc%age
  read(iunit, err=900, end=900) cc%leafage
  read(iunit, err=900, end=900) cc%mu
  read(iunit, err=900, end=900) cc%dbh
  read(iunit, err=900, end=900) cc%height
  read(iunit, err=900, end=900) cc%Acrown
  read(iunit, err=900, end=900) cc%Aleaf
  read(iunit, err=900, end=900) cc%lai
  read(iunit, err=900, end=900) cc%D_bark
  read(iunit, err=900, end=900) cc%bl_max
  read(iunit, err=900, end=900) cc%br_max
  read(iunit, err=900, end=900) cc%CSAsw
  read(iunit, err=900, end=900) cc%DBH_ys
  read(iunit, err=900, end=900) cc%extinct
  read(iunit, err=900, end=900) cc%w_scale

  !-- phenology --
  read(iunit, err=900, end=900) cc%GDD
  read(iunit, err=900, end=900) cc%ALT
  read(iunit, err=900, end=900) cc%AWD
  read(iunit, err=900, end=900) cc%Tc_OFF
  read(iunit, err=900, end=900) cc%GDD_ON

  !-- C pools --
  read(iunit, err=900, end=900) cc%bl
  read(iunit, err=900, end=900) cc%br
  read(iunit, err=900, end=900) cc%bsw
  read(iunit, err=900, end=900) cc%bHW
  read(iunit, err=900, end=900) cc%seedC
  read(iunit, err=900, end=900) cc%nsc

  !-- N pools --
  ! NSCmax and NSNmax are not saved; recomputed from bl_max and br_max
  cc%NSCmax = ccNSCmax(cc)
  cc%NSNmax = ccNSNmax(cc)
  read(iunit, err=900, end=900) cc%NSN
  read(iunit, err=900, end=900) cc%leafN
  read(iunit, err=900, end=900) cc%swN
  read(iunit, err=900, end=900) cc%hwN
  read(iunit, err=900, end=900) cc%rootN
  read(iunit, err=900, end=900) cc%seedN
  read(iunit, err=900, end=900) cc%facuC

  !-- root / water supply --
  read(iunit, err=900, end=900) cc%K_r
  read(iunit, err=900, end=900) cc%r_r
  read(iunit, err=900, end=900) cc%root_zeta
  read(iunit, err=900, end=900) cc%W_supply

  !-- root arrays --
  read(iunit, err=900, end=900) cc%root_length(1:soil_L)
  read(iunit, err=900, end=900) cc%ArootL(1:soil_L)
  read(iunit, err=900, end=900) cc%WupL(1:soil_L)
  read(iunit, err=900, end=900) cc%uptake_frac(1:soil_L)

  !-- hydraulics scalars --
  read(iunit, err=900, end=900) cc%psi_s0
  read(iunit, err=900, end=900) cc%psi_leaf
  read(iunit, err=900, end=900) cc%psi_stem
  read(iunit, err=900, end=900) cc%H_leaf
  read(iunit, err=900, end=900) cc%H_stem
  read(iunit, err=900, end=900) cc%W_lf
  read(iunit, err=900, end=900) cc%W_sw
  read(iunit, err=900, end=900) cc%W_hw
  read(iunit, err=900, end=900) cc%Wmax_l
  read(iunit, err=900, end=900) cc%Wmax_s
  read(iunit, err=900, end=900) cc%Wmin_l
  read(iunit, err=900, end=900) cc%Wmin_s
  read(iunit, err=900, end=900) cc%Q_stem
  read(iunit, err=900, end=900) cc%Q_leaf
  read(iunit, err=900, end=900) cc%Ktrunk
  read(iunit, err=900, end=900) cc%Asap
  read(iunit, err=900, end=900) cc%Atrunk
  read(iunit, err=900, end=900) cc%treeHU
  read(iunit, err=900, end=900) cc%treeW0
  read(iunit, err=900, end=900) cc%An_op
  read(iunit, err=900, end=900) cc%An_cl

  !-- xylem ring arrays --
  read(iunit, err=900, end=900) cc%Kx(1:Ysw_max)
  read(iunit, err=900, end=900) cc%WTC0(1:Ysw_max)
  read(iunit, err=900, end=900) cc%accH(1:Ysw_max)
  read(iunit, err=900, end=900) cc%plcH(1:Ysw_max)
  read(iunit, err=900, end=900) cc%farea(1:Ysw_max)
  read(iunit, err=900, end=900) cc%Rring(1:Ysw_max)
  read(iunit, err=900, end=900) cc%Lring(1:Ysw_max)
  read(iunit, err=900, end=900) cc%Aring(1:Ysw_max)
  read(iunit, err=900, end=900) cc%Kring(1:Ysw_max)

    return

900 continue
    write(*,'(A,A)') 'read_restart: error or end of file while reading cohort data from ', trim(restart_path())
    error stop 'read_restart failed: the restart file is incomplete or from another version'
  end subroutine read_cohort

!============================================================================
  subroutine write_restart(land)
!============================================================================
    type(land_grid_type), pointer, intent(in) :: land

    type(vegn_tile_type), pointer :: vegn
    integer :: iunit, i, ios, total_cohorts, ntiles
    character(len=512) :: fpath

    fpath = restart_path()
    open(newunit=iunit, file=trim(fpath), form='unformatted', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A,A,I0)') 'write_restart: cannot open ', trim(fpath), ', iostat=', ios
      error stop 'write_restart failed'
    endif

    !--- header: checked when the file is read ---
    write(iunit) RESTART_MAGIC, RESTART_VERSION
    write(iunit) soil_L, CLmax, Ysw_max

    ntiles = land%nTiles
    write(iunit) ntiles

    total_cohorts = 0
    vegn => land%firstVegn
    do while (associated(vegn))

      !--- tile integer scalars ---
      write(iunit) vegn%tileID
      write(iunit) vegn%landuse
      write(iunit) vegn%n_cohorts
      write(iunit) vegn%n_years
      write(iunit) vegn%soiltype
      write(iunit) vegn%HarvYrs

      !--- tile real scalars ---
      write(iunit) vegn%age
      write(iunit) vegn%tc_daily
      write(iunit) vegn%tc_pheno
      write(iunit) vegn%FLDCAP
      write(iunit) vegn%WILTPT
      write(iunit) vegn%W0topSL
      write(iunit) vegn%thetaS     ! soil wetness of the last step, used by Soil_BGC and phenology
      write(iunit) vegn%soilwater

      !--- tile arrays ---
      write(iunit) vegn%SOC(1:5)
      write(iunit) vegn%SON(1:5)
      write(iunit) vegn%wcl(1:soil_L)
      write(iunit) vegn%freewater(1:soil_L)
      write(iunit) vegn%psi_soil(1:soil_L)
      write(iunit) vegn%K_soil(1:soil_L)
      write(iunit) vegn%LAI_L(1:CLmax)
      write(iunit) vegn%f_gap(1:CLmax)
      write(iunit) vegn%CAI_L(1:CLmax)
      write(iunit) vegn%kp(1:CLmax)

      !--- tile N state ---
      write(iunit) vegn%mineralN
      write(iunit) vegn%previousN
      write(iunit) vegn%initialN0

      !--- tile vegetation water ---
      write(iunit) vegn%W_lf
      write(iunit) vegn%W_sw
      write(iunit) vegn%W_hw

      !--- tile fire/climate ---
      write(iunit) vegn%annualPET
      write(iunit) vegn%YearlyTmp
      write(iunit) vegn%Frisk
      write(iunit) vegn%Pfire
      write(iunit) vegn%TreeCA
      write(iunit) vegn%GrassCA
      write(iunit) vegn%GrassBM
      write(iunit) vegn%t_ann
      write(iunit) vegn%t_cold
      write(iunit) vegn%p_ann
      write(iunit) vegn%ncm

      !--- cohorts ---
      do i = 1, vegn%n_cohorts
        call write_cohort(iunit, vegn%cohorts(i))
        total_cohorts = total_cohorts + 1
      enddo

      !--- initial cohorts (used for species recovery) ---
      write(iunit) vegn%n_initialCC
      do i = 1, vegn%n_initialCC
        call write_cohort(iunit, vegn%initialCC(i))
      enddo

      !--- animal cohorts (state variables only) ---
      write(iunit) vegn%n_ani_cohorts
      do i = 1, vegn%n_ani_cohorts
        write(iunit) vegn%ani_cohorts(i)%aft
        write(iunit) vegn%ani_cohorts(i)%nindivs
        write(iunit) vegn%ani_cohorts(i)%age
        write(iunit) vegn%ani_cohorts(i)%body_C
      enddo

      vegn => vegn%next
    enddo ! tiles

    close(iunit)
    write(*,'(A,I0,A,I0,A,A)') 'Restart written: ', total_cohorts, &
      ' cohorts in ', ntiles, ' tiles -> ', trim(fpath)

  end subroutine write_restart

!============================================================================
  subroutine read_restart(land)
!============================================================================
    type(land_grid_type), pointer, intent(inout) :: land

    type(vegn_tile_type), pointer :: vegn, pveg
    integer :: iunit, i, ios, ntiles, itile, total_cohorts
    integer :: magic, version, n_soil, n_cl, n_ysw
    character(len=512) :: fpath

    fpath = restart_path()
    open(newunit=iunit, file=trim(fpath), form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A,A,I0)') 'read_restart: cannot open ', trim(fpath), ', iostat=', ios
      error stop 'read_restart failed: restart file not found'
    endif

    !--- header ---
    read(iunit, iostat=ios) magic, version
    if (ios /= 0 .or. magic /= RESTART_MAGIC) then
      write(*,'(A,A)') 'read_restart: not a BiomeE restart file of this format: ', trim(fpath)
      error stop 'read_restart failed: unknown file format (old files must be regenerated)'
    endif
    if (version /= RESTART_VERSION) then
      write(*,'(A,I0,A,I0)') 'read_restart: file version ', version, ', code version ', RESTART_VERSION
      error stop 'read_restart failed: version mismatch'
    endif
    read(iunit, err=900, end=900) n_soil, n_cl, n_ysw
    if (n_soil /= soil_L .or. n_cl /= CLmax .or. n_ysw /= Ysw_max) then
      write(*,'(A,3I5,A,3I5)') 'read_restart: array sizes (soil_L,CLmax,Ysw_max) in file:', &
            n_soil, n_cl, n_ysw, ' vs code:', soil_L, CLmax, Ysw_max
      error stop 'read_restart failed: array sizes differ from this build'
    endif

    read(iunit, err=900, end=900) ntiles
    land%nTiles  = ntiles
    nullify(land%firstVegn)
    pveg => NULL()
    total_cohorts = 0

    do itile = 1, ntiles

      allocate(vegn)
      ! Initialise all pointer members to null so deallocation is safe
      nullify(vegn%cohorts)
      nullify(vegn%initialCC)
      nullify(vegn%ani_cohorts)
      nullify(vegn%prev)
      nullify(vegn%next)
      vegn%n_initialCC = 0
      vegn%n_ani_cohorts = 0

      !--- tile integer scalars ---
      read(iunit, err=900, end=900) vegn%tileID
      read(iunit, err=900, end=900) vegn%landuse
      read(iunit, err=900, end=900) vegn%n_cohorts
      read(iunit, err=900, end=900) vegn%n_years
      read(iunit, err=900, end=900) vegn%soiltype
      read(iunit, err=900, end=900) vegn%HarvYrs

      !--- tile real scalars ---
      read(iunit, err=900, end=900) vegn%age
      read(iunit, err=900, end=900) vegn%tc_daily
      read(iunit, err=900, end=900) vegn%tc_pheno
      read(iunit, err=900, end=900) vegn%FLDCAP
      read(iunit, err=900, end=900) vegn%WILTPT
      read(iunit, err=900, end=900) vegn%W0topSL
      read(iunit, err=900, end=900) vegn%thetaS
      read(iunit, err=900, end=900) vegn%soilwater

      !--- tile arrays ---
      read(iunit, err=900, end=900) vegn%SOC(1:5)
      read(iunit, err=900, end=900) vegn%SON(1:5)
      read(iunit, err=900, end=900) vegn%wcl(1:soil_L)
      read(iunit, err=900, end=900) vegn%freewater(1:soil_L)
      read(iunit, err=900, end=900) vegn%psi_soil(1:soil_L)
      read(iunit, err=900, end=900) vegn%K_soil(1:soil_L)
      read(iunit, err=900, end=900) vegn%LAI_L(1:CLmax)
      read(iunit, err=900, end=900) vegn%f_gap(1:CLmax)
      read(iunit, err=900, end=900) vegn%CAI_L(1:CLmax)
      read(iunit, err=900, end=900) vegn%kp(1:CLmax)

      !--- tile N state ---
      read(iunit, err=900, end=900) vegn%mineralN
      read(iunit, err=900, end=900) vegn%previousN
      read(iunit, err=900, end=900) vegn%initialN0

      !--- tile vegetation water ---
      read(iunit, err=900, end=900) vegn%W_lf
      read(iunit, err=900, end=900) vegn%W_sw
      read(iunit, err=900, end=900) vegn%W_hw

      !--- tile fire/climate ---
      read(iunit, err=900, end=900) vegn%annualPET
      read(iunit, err=900, end=900) vegn%YearlyTmp
      read(iunit, err=900, end=900) vegn%Frisk
      read(iunit, err=900, end=900) vegn%Pfire
      read(iunit, err=900, end=900) vegn%TreeCA
      read(iunit, err=900, end=900) vegn%GrassCA
      read(iunit, err=900, end=900) vegn%GrassBM
      read(iunit, err=900, end=900) vegn%t_ann
      read(iunit, err=900, end=900) vegn%t_cold
      read(iunit, err=900, end=900) vegn%p_ann
      read(iunit, err=900, end=900) vegn%ncm

      !--- cohorts ---
      allocate(vegn%cohorts(vegn%n_cohorts))
      do i = 1, vegn%n_cohorts
        call read_cohort(iunit, vegn%cohorts(i))
        total_cohorts = total_cohorts + 1
      enddo

      !--- initial cohorts ---
      read(iunit, err=900, end=900) vegn%n_initialCC
      if (vegn%n_initialCC > 0) then
        allocate(vegn%initialCC(vegn%n_initialCC))
        do i = 1, vegn%n_initialCC
          call read_cohort(iunit, vegn%initialCC(i))
        enddo
      endif

      !--- animal cohorts ---
      read(iunit, err=900, end=900) vegn%n_ani_cohorts
      if (vegn%n_ani_cohorts > 0) then
        allocate(vegn%ani_cohorts(vegn%n_ani_cohorts))
        do i = 1, vegn%n_ani_cohorts
          read(iunit, err=900, end=900) vegn%ani_cohorts(i)%aft
          read(iunit, err=900, end=900) vegn%ani_cohorts(i)%nindivs
          read(iunit, err=900, end=900) vegn%ani_cohorts(i)%age
          read(iunit, err=900, end=900) vegn%ani_cohorts(i)%body_C
        enddo
      endif

      !--- link tile into the list ---
      if (itile == 1) then
        land%firstVegn => vegn
      else
        pveg%next => vegn
        vegn%prev => pveg
      endif
      pveg => vegn
      vegn => NULL()

    enddo ! tiles

    close(iunit)
    write(*,'(A,I0,A,I0,A,A)') 'Restart read: ', total_cohorts, &
      ' cohorts in ', ntiles, ' tiles <- ', trim(fpath)
    return

900 continue
    write(*,'(A,A)') 'read_restart: error or end of file while reading ', trim(fpath)
    error stop 'read_restart failed: the restart file is incomplete'
  end subroutine read_restart

!----------------------------------------------------------------------------
end module restart_mod
