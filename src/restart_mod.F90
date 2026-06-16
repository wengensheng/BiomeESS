!============================================================================
! restart_mod.F90 -- checkpoint/restart I/O for BiomeE-Global
!
! Usage:
!   call write_restart(land)   -- at end of run when do_restart_write=.True.
!   call read_restart(land)    -- at start of run when do_restart_read=.True.
!
! Binary format: Fortran unformatted sequential (default record-length units).
! The READ and WRITE record order is byte-for-byte identical.
!
! For global runs where GridID > 0 the file name is constructed by inserting
! the GridID before the extension, e.g.:
!   restart_file = './output/BiomeE_restart.bin'
!   -> actual path: './output/BiomeE_restart_205249.bin'
!============================================================================

module restart_mod
  use datatypes
  implicit none
  private

  public :: write_restart, read_restart

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
  subroutine write_restart(land)
!============================================================================
    type(land_grid_type), pointer, intent(in) :: land

    type(vegn_tile_type), pointer :: vegn
    integer :: iunit, i, total_cohorts, ntiles
    character(len=512) :: fpath

    fpath  = restart_path()
    iunit  = 99 ! scratch unit; not used elsewhere during this call

    open(unit=iunit, file=trim(fpath), form='unformatted', &
         status='replace', action='write')

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

      !--- tile real scalars ---
      write(iunit) vegn%age
      write(iunit) vegn%tc_daily
      write(iunit) vegn%tc_pheno
      write(iunit) vegn%FLDCAP
      write(iunit) vegn%WILTPT
      write(iunit) vegn%W0topSL

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
        associate(cc => vegn%cohorts(i))

          !-- integer fields --
          write(iunit) cc%species
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
          write(iunit) cc%NSNmax
          write(iunit) cc%NSN
          write(iunit) cc%leafN
          write(iunit) cc%swN
          write(iunit) cc%hwN
          write(iunit) cc%rootN
          write(iunit) cc%seedN

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

        end associate
        total_cohorts = total_cohorts + 1
      enddo ! cohorts

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
    integer :: iunit, i, ntiles, itile, total_cohorts
    character(len=512) :: fpath
    integer :: tileID_tmp, landuse_tmp, n_cohorts_tmp, n_years_tmp, soiltype_tmp

    fpath = restart_path()
    iunit = 99

    open(unit=iunit, file=trim(fpath), form='unformatted', &
         status='old', action='read')

    read(iunit) ntiles
    land%nTiles  = ntiles
    nullify(land%firstVegn)
    pveg => NULL()
    total_cohorts = 0

    do itile = 1, ntiles

      allocate(vegn)
      ! Initialise all pointer members to null so deallocation is safe
      nullify(vegn%cohorts)
      nullify(vegn%initialCC)
      nullify(vegn%prev)
      nullify(vegn%next)
      vegn%n_initialCC = 0

      !--- tile integer scalars ---
      read(iunit) vegn%tileID
      read(iunit) vegn%landuse
      read(iunit) vegn%n_cohorts
      read(iunit) vegn%n_years
      read(iunit) vegn%soiltype

      !--- tile real scalars ---
      read(iunit) vegn%age
      read(iunit) vegn%tc_daily
      read(iunit) vegn%tc_pheno
      read(iunit) vegn%FLDCAP
      read(iunit) vegn%WILTPT
      read(iunit) vegn%W0topSL

      !--- tile arrays ---
      read(iunit) vegn%SOC(1:5)
      read(iunit) vegn%SON(1:5)
      read(iunit) vegn%wcl(1:soil_L)
      read(iunit) vegn%freewater(1:soil_L)
      read(iunit) vegn%psi_soil(1:soil_L)
      read(iunit) vegn%K_soil(1:soil_L)
      read(iunit) vegn%LAI_L(1:CLmax)
      read(iunit) vegn%f_gap(1:CLmax)
      read(iunit) vegn%CAI_L(1:CLmax)
      read(iunit) vegn%kp(1:CLmax)

      !--- tile N state ---
      read(iunit) vegn%mineralN
      read(iunit) vegn%previousN
      read(iunit) vegn%initialN0

      !--- tile vegetation water ---
      read(iunit) vegn%W_lf
      read(iunit) vegn%W_sw
      read(iunit) vegn%W_hw

      !--- tile fire/climate ---
      read(iunit) vegn%annualPET
      read(iunit) vegn%YearlyTmp
      read(iunit) vegn%Frisk
      read(iunit) vegn%Pfire
      read(iunit) vegn%TreeCA
      read(iunit) vegn%GrassCA
      read(iunit) vegn%GrassBM
      read(iunit) vegn%t_ann
      read(iunit) vegn%t_cold
      read(iunit) vegn%p_ann
      read(iunit) vegn%ncm

      !--- cohorts ---
      allocate(vegn%cohorts(vegn%n_cohorts))

      do i = 1, vegn%n_cohorts
        associate(cc => vegn%cohorts(i))

          !-- integer fields --
          read(iunit) cc%species
          read(iunit) cc%layer
          read(iunit) cc%status
          read(iunit) cc%Ngd
          read(iunit) cc%Ndm
          read(iunit) cc%Ncd
          read(iunit) cc%Nrings
          read(iunit) cc%phenotype
          read(iunit) cc%pt

          !-- structure real scalars --
          read(iunit) cc%nindivs
          read(iunit) cc%layerfrac
          read(iunit) cc%topyear
          read(iunit) cc%age
          read(iunit) cc%leafage
          read(iunit) cc%mu
          read(iunit) cc%dbh
          read(iunit) cc%height
          read(iunit) cc%Acrown
          read(iunit) cc%Aleaf
          read(iunit) cc%lai
          read(iunit) cc%D_bark
          read(iunit) cc%bl_max
          read(iunit) cc%br_max
          read(iunit) cc%CSAsw
          read(iunit) cc%DBH_ys
          read(iunit) cc%extinct
          read(iunit) cc%w_scale

          !-- phenology --
          read(iunit) cc%GDD
          read(iunit) cc%ALT
          read(iunit) cc%AWD
          read(iunit) cc%Tc_OFF
          read(iunit) cc%GDD_ON

          !-- C pools --
          read(iunit) cc%bl
          read(iunit) cc%br
          read(iunit) cc%bsw
          read(iunit) cc%bHW
          read(iunit) cc%seedC
          read(iunit) cc%nsc

          !-- N pools --
          read(iunit) cc%NSNmax
          read(iunit) cc%NSN
          read(iunit) cc%leafN
          read(iunit) cc%swN
          read(iunit) cc%hwN
          read(iunit) cc%rootN
          read(iunit) cc%seedN

          !-- root / water supply --
          read(iunit) cc%K_r
          read(iunit) cc%r_r
          read(iunit) cc%root_zeta
          read(iunit) cc%W_supply

          !-- root arrays --
          read(iunit) cc%root_length(1:soil_L)
          read(iunit) cc%ArootL(1:soil_L)
          read(iunit) cc%WupL(1:soil_L)
          read(iunit) cc%uptake_frac(1:soil_L)

          !-- hydraulics scalars --
          read(iunit) cc%psi_s0
          read(iunit) cc%psi_leaf
          read(iunit) cc%psi_stem
          read(iunit) cc%H_leaf
          read(iunit) cc%H_stem
          read(iunit) cc%W_lf
          read(iunit) cc%W_sw
          read(iunit) cc%W_hw
          read(iunit) cc%Wmax_l
          read(iunit) cc%Wmax_s
          read(iunit) cc%Wmin_l
          read(iunit) cc%Wmin_s
          read(iunit) cc%Q_stem
          read(iunit) cc%Q_leaf
          read(iunit) cc%Ktrunk
          read(iunit) cc%Asap
          read(iunit) cc%Atrunk
          read(iunit) cc%treeHU
          read(iunit) cc%treeW0
          read(iunit) cc%An_op
          read(iunit) cc%An_cl

          !-- xylem ring arrays --
          read(iunit) cc%Kx(1:Ysw_max)
          read(iunit) cc%WTC0(1:Ysw_max)
          read(iunit) cc%accH(1:Ysw_max)
          read(iunit) cc%plcH(1:Ysw_max)
          read(iunit) cc%farea(1:Ysw_max)
          read(iunit) cc%Rring(1:Ysw_max)
          read(iunit) cc%Lring(1:Ysw_max)
          read(iunit) cc%Aring(1:Ysw_max)
          read(iunit) cc%Kring(1:Ysw_max)

        end associate
        total_cohorts = total_cohorts + 1
      enddo ! cohorts

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

  end subroutine read_restart

!----------------------------------------------------------------------------
end module restart_mod
