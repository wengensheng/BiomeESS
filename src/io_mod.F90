!#define USE_NETCDF

module io_mod
! Data input and output
 use datatypes
#ifdef USE_NETCDF
 use netcdf
#endif

 implicit none

! ------ public subroutines ---------
public :: set_up_output_files,read_FACEforcing,read_NACPforcing
public :: Zero_diagnostics, hourly_diagnostics, daily_diagnostics, &
          annual_diagnostics
public :: vegn_sum_tile
!---------------------------------

 contains
!====================== Subroutines ======================================

!=================================================
! Weng, 2021-06-02
subroutine vegn_sum_tile(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  !----- local var --------------
  type(cohort_type),pointer :: cc
  integer :: i, layer

  vegn%NSC     = 0.0
  vegn%SeedC   = 0.0
  vegn%leafC   = 0.0
  vegn%rootC   = 0.0
  vegn%SapwoodC= 0.0
  vegn%WoodC   = 0.0

  vegn%NSN     = 0.0
  vegn%SeedN   = 0.0
  vegn%leafN   = 0.0
  vegn%rootN   = 0.0
  vegn%SapwoodN= 0.0
  vegn%WoodN   = 0.0

  vegn%W_stem = 0.0
  vegn%W_dead = 0.0
  vegn%W_leaf = 0.0

  vegn%LAI     = 0.0
  vegn%CAI     = 0.0

  vegn%LAIlayer = 0.0
  vegn%f_gap    = 0.0
  vegn%treecover = 0.0
  vegn%grasscover = 0.0
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     associate ( sp => spdata(cc%species))

     ! update accumulative LAI for each corwn layer
     layer = Max (1, Min(cc%layer,9)) ! between 1~9
     vegn%LAIlayer(layer) = vegn%LAIlayer(layer) + &
                            cc%Aleaf * cc%nindivs/(1.0-sp%f_cGap)
     vegn%f_gap(layer)    = vegn%f_gap(layer)    +  &
                            cc%Acrown * cc%nindivs * sp%f_cGap

    ! For reporting
    ! Vegn C pools:
     vegn%NSC     = vegn%NSC     + cc%NSC    * cc%nindivs
     vegn%SeedC   = vegn%SeedC   + cc%seedC  * cc%nindivs
     vegn%leafC   = vegn%leafC   + cc%bl     * cc%nindivs
     vegn%rootC   = vegn%rootC   + cc%br     * cc%nindivs
     vegn%SapwoodC= vegn%SapwoodC+ cc%bsw    * cc%nindivs
     vegn%woodC   = vegn%woodC   + cc%bHW    * cc%nindivs
     vegn%CAI     = vegn%CAI     + cc%Acrown * cc%nindivs
     vegn%LAI     = vegn%LAI     + cc%Aleaf  * cc%nindivs
    ! Vegn N pools
     vegn%NSN     = vegn%NSN   + cc%NSN      * cc%nindivs
     vegn%SeedN   = vegn%SeedN + cc%seedN    * cc%nindivs
     vegn%leafN   = vegn%leafN + cc%leafN    * cc%nindivs
     vegn%rootN   = vegn%rootN + cc%rootN    * cc%nindivs
     vegn%SapwoodN= vegn%SapwoodN + cc%sapwN * cc%nindivs
     vegn%woodN   = vegn%woodN    + cc%woodN * cc%nindivs
     ! Vegn water pools
     vegn%W_stem = vegn%W_stem   + cc%W_stem * cc%nindivs
     vegn%W_dead = vegn%W_dead   + cc%W_dead * cc%nindivs
     vegn%W_leaf = vegn%W_leaf   + cc%W_leaf * cc%nindivs

     ! Update tree and grass cover
     if(sp%lifeform==0) then
         if(cc%layer == 1)vegn%grasscover = vegn%grasscover + cc%Acrown*cc%nindivs
     elseif(sp%lifeform==1 .and. cc%height > 4.0)then ! for trees in the top layer
         vegn%treecover = vegn%treecover + cc%Acrown*cc%nindivs
     endif

     end associate
  enddo

end subroutine vegn_sum_tile

!================= Diagnostics============================================
! Weng, 2016-11-28
subroutine Zero_diagnostics(vegn)
  ! for annual update
  type(vegn_tile_type), intent(inout) :: vegn
  !-------local var
  type(cohort_type),pointer :: cc
  integer :: i
  !daily
  vegn%dailyfixedN = 0.
  vegn%dailyPrcp = 0.0
  vegn%dailyTrsp = 0.0
  vegn%dailyEvap = 0.0
  vegn%dailyRoff = 0.0
  vegn%dailyNup  = 0.0
  vegn%dailyGPP = 0.0
  vegn%dailyNPP = 0.0
  vegn%dailyResp = 0.0
  vegn%dailyRh   = 0.0

  !annual
  vegn%annualfixedN = 0.
  vegn%annualPrcp = 0.0
  vegn%annualTrsp = 0.0
  vegn%annualEvap = 0.0
  vegn%annualRoff = 0.0
  vegn%annualGPP = 0.0
  vegn%annualNPP = 0.0
  vegn%annualResp = 0.0
  vegn%annualRh   = 0.0
  vegn%N_P2S_yr  = 0.
  vegn%annualN   = 0.
  vegn%Nloss_yr  = 0.
  vegn%annualNup  = 0.0

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     cc%C_growth = 0.0
     cc%N_growth = 0.0
     cc%gpp      = 0.0
     cc%npp      = 0.0
     cc%resp     = 0.0
     cc%resl     = 0.0
     cc%resr     = 0.0
     cc%resg     = 0.0
     cc%transp   = 0.0
     !daily
     cc%dailyTrsp = 0.0
     cc%dailyGPP = 0.0
     cc%dailyNPP = 0.0
     cc%dailyResp= 0.0
     cc%dailyNup   = 0.0
     cc%dailyfixedN = 0.0
     ! annual
     cc%annualTrsp = 0.0
     cc%annualGPP = 0.0
     cc%annualNPP = 0.0
     cc%annualResp= 0.0
     cc%annualNup   = 0.0
     cc%annualfixedN = 0.0
     cc%NPPleaf   = 0.0
     cc%NPProot   = 0.0
     cc%NPPwood   = 0.0
     cc%DBH_ys    = cc%DBH
  enddo
end subroutine Zero_diagnostics

!=========================================================================
! Hourly fluxes sum to daily
 subroutine hourly_diagnostics(vegn,forcing,iyears,idoy,ihour,iday,fno1,fno2)
  type(vegn_tile_type), intent(inout) :: vegn
  type(climate_data_type),intent(in):: forcing
  integer, intent(in) :: iyears,idoy,ihour,iday
  integer, intent(in) :: fno1, fno2

  !-------local var ------
  type(cohort_type), pointer :: cc    ! current cohort
  integer :: i

  ! Tile summary
  vegn%GPP    = 0.; vegn%fixedN = 0.
  vegn%NPP    = 0.; vegn%Resp   = 0.
  vegn%transp = 0.
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! cohort daily
     cc%dailyTrsp = cc%dailyTrsp + cc%transp ! kg day-1
     cc%dailyGPP  = cc%dailygpp  + cc%gpp ! kg day-1
     cc%dailyNPP  = cc%dailyNpp  + cc%Npp ! kg day-1
     cc%dailyResp = cc%dailyResp + cc%Resp ! kg day-1
     cc%dailyfixedN  = cc%dailyfixedN  + cc%fixedN ! kg day-1

     ! Tile hourly
     vegn%GPP    = vegn%GPP    + cc%gpp    * cc%nindivs
     vegn%NPP    = vegn%NPP    + cc%Npp    * cc%nindivs
     vegn%Resp   = vegn%Resp   + cc%Resp   * cc%nindivs
     vegn%transp = vegn%transp + cc%transp * cc%nindivs
     vegn%fixedN = vegn%fixedN + cc%fixedN * cc%nindivs
  enddo

  !! Output horly diagnostics
  If(outputhourly .and. iday > totdays-366*5 ) then !  .and. ihour==12
    write(fno1,'(4(I8,","))')vegn%n_cohorts
    do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        write(fno1,'(6(I8,","),40(F12.4,","))')           &
          iyears,idoy,ihour,cc%ccID,cc%species,cc%layer,  &
          cc%nindivs*10000,cc%dbh,cc%height,cc%Acrown,    &
          cc%bl,cc%LAI,cc%gpp,cc%npp,cc%transp,           &
          cc%psi_leaf,cc%psi_stem,cc%W_leaf,cc%W_stem
    enddo
    ! Hourly tile
    associate ( cc1 => vegn%cohorts(1))
    write(fno2,'(3(I5,","),30(E12.4,","),25(F12.4,","))')  &
      iyears,idoy,ihour,forcing%radiation,forcing%Tair,    &
      forcing%rain,vegn%GPP,vegn%resp,vegn%transp,         &
      vegn%evap,vegn%runoff,vegn%soilwater,                &
      vegn%wcl(2),vegn%psi_soil(2),vegn%K_soil(2),         &
      cc1%bl,cc1%psi_leaf,cc1%psi_stem,cc1%W_leaf,         &
      cc1%W_stem,cc1%transp
    end associate
  endif
  ! Daily summary:
  vegn%dailyNup  = vegn%dailyNup  + vegn%N_uptake
  vegn%dailyGPP  = vegn%dailyGPP  + vegn%gpp
  vegn%dailyNPP  = vegn%dailyNPP  + vegn%npp
  vegn%dailyResp = vegn%dailyResp + vegn%resp
  vegn%dailyRh   = vegn%dailyRh   + vegn%rh
  vegn%dailyTrsp = vegn%dailyTrsp + vegn%transp
  vegn%dailyEvap = vegn%dailyEvap + vegn%evap
  vegn%dailyRoff = vegn%dailyRoff + vegn%runoff
  vegn%dailyPrcp = vegn%dailyPrcp + forcing%rain * step_seconds
  vegn%dailyfixedN  = vegn%dailyfixedN  + vegn%fixedN

end subroutine hourly_diagnostics

!============================================
subroutine daily_diagnostics(vegn,iyears,idoy,iday,fno3,fno4)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: iyears,idoy,iday,fno3,fno4

  !-------local var ------
  type(cohort_type), pointer :: cc    ! current cohort
  integer :: i,j

  ! Output and zero daily variables
  !!! daily !! cohorts output
  if(outputdaily.and. iday>equi_days) &
  write(fno3,'(3(I6,","))')iyears, idoy,vegn%n_cohorts
  do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      if(outputdaily.and. iday>equi_days) &
      write(fno3,'(9(I5,","),1(F12.4,","),50(F12.4,","))')  &
            iyears,idoy,i, cc%ccID,cc%species,              &
            cc%layer,cc%status, cc%ndm, cc%ncd,             &
            cc%nindivs*10000, cc%layerfrac, cc%LAI,         &
            cc%dailygpp,cc%dailyresp,cc%dailytrsp,          &
            cc%NPPleaf,cc%NPProot,cc%NPPwood,               &
            cc%W_leaf,cc%W_stem,cc%W_dead,                  &
            cc%NSC,cc%seedC,cc%bl,cc%br,cc%bsw,cc%bHW,      &
            cc%NSN*1000, cc%seedN*1000, cc%leafN*1000,      &
            cc%rootN*1000,cc%sapwN*1000,cc%woodN*1000,      &
            cc%gdd,cc%ALT

      ! annual sum
      cc%annualGPP  = cc%annualGPP  + cc%dailyGPP
      cc%annualNPP  = cc%annualNPP  + cc%dailyNPP
      cc%annualResp = cc%annualResp + cc%dailyResp
      cc%annualTrsp = cc%annualTrsp + cc%dailyTrsp
      cc%annualfixedN = cc%annualfixedN + cc%dailyfixedN
      ! Zero Daily variables
      cc%dailyTrsp = 0.0
      cc%dailyGPP = 0.0
      cc%dailyNPP = 0.0
      cc%dailyResp = 0.0
      cc%dailyfixedN = 0.0
  enddo
  !! Tile level, daily
  if(outputdaily.and. iday>equi_days) then
     write(fno4,'(2(I5,","),65(F12.4,","))') iyears, idoy,       &
        vegn%tc_pheno, vegn%dailyPrcp, vegn%soilwater,           &
        vegn%dailyTrsp, vegn%dailyEvap,vegn%dailyRoff,           &
        vegn%wcl(1)*thksl(1)*1000.,vegn%wcl(2)*thksl(2)*1000.,   &
        vegn%wcl(3)*thksl(3)*1000.,                              &
        vegn%LAI,vegn%dailyGPP, vegn%dailyResp, vegn%dailyRh,    &
        vegn%W_leaf,vegn%W_stem,vegn%W_dead,                     &
        vegn%NSC, vegn%SeedC, vegn%leafC, vegn%rootC,            &
        vegn%SapwoodC, vegn%woodC,                               &
        vegn%NSN*1000, vegn%SeedN*1000, vegn%leafN*1000,         &
        vegn%rootN*1000, vegn%SapwoodN *1000, vegn%WoodN *1000,  &
        (vegn%SOC(j),j=1,5), (vegn%SON(j)*1000,j=1,5),           &
        vegn%mineralN*1000,vegn%dailyNup*1000,vegn%kp(1)
  endif
  !annual tile
  ! Annual summary:
  vegn%annualNup  = vegn%annualNup  + vegn%dailyNup
  vegn%annualGPP  = vegn%annualGPP  + vegn%dailygpp
  vegn%annualNPP  = vegn%annualNPP  + vegn%dailynpp
  vegn%annualResp = vegn%annualResp + vegn%dailyresp
  vegn%annualRh   = vegn%annualRh   + vegn%dailyrh
  vegn%annualPrcp = vegn%annualPrcp + vegn%dailyPrcp
  vegn%annualTrsp = vegn%annualTrsp + vegn%dailytrsp
  vegn%annualEvap = vegn%annualEvap + vegn%dailyevap
  vegn%annualRoff = vegn%annualRoff + vegn%dailyRoff
  vegn%annualfixedN  = vegn%annualfixedN  + vegn%dailyfixedN

  ! zero:
  vegn%dailyNup  = 0.0
  vegn%dailyGPP  = 0.0
  vegn%dailyNPP  = 0.0
  vegn%dailyResp = 0.0
  vegn%dailyRh   = 0.0
  vegn%dailyPrcp = 0.0
  vegn%dailyTrsp = 0.0
  vegn%dailyEvap = 0.0
  vegn%dailyRoff = 0.0
  vegn%dailyfixedN = 0.0

end subroutine daily_diagnostics

!======================================================
 subroutine annual_diagnostics(vegn, iyears,f1,f2)
   type(vegn_tile_type), intent(inout) :: vegn
   integer, intent(in) :: f1,f2, iyears

    ! --------local var --------
    type(cohort_type), pointer :: cc
    real treeG, fseed, fleaf, froot,fwood,dDBH
    real :: plantC, plantN, soilC, soilN
    integer :: i,j

    write(*, '(2(I6,","),3(F9.3,","))')iyears, vegn%n_cohorts,vegn%age,vegn%cohorts(1)%age
    write(*,'(1(a6,","),2(a4,","),25(a9,","))')       &
            'cID','PFT','L', 'n','f_CA','dDBH',       &
            'DBH','Height','Acrown','NSC','GPP','mu', &
            'Atrunk','Asap','Ktree','treeHU','treeW0'

    ! Cohotrs ouput
    write(f1,'(2(I6,","),1(F9.2,","))')iyears, vegn%n_cohorts
    do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)
        treeG = cc%seedC + cc%NPPleaf + cc%NPProot + cc%NPPwood
        fseed = cc%seedC/treeG
        fleaf = cc%NPPleaf/treeG
        froot = cc%NPProot/treeG
        fwood = cc%NPPwood/treeG
        dDBH  = (cc%DBH - cc%DBH_ys) * 1000.
        write(f1,'(4(I8,","),300(E15.4,","))')           &
            iyears,cc%ccID,cc%species,cc%layer,          &
            cc%nindivs*10000, cc%layerfrac,dDBH,         &
            cc%dbh,cc%height,cc%Acrown,                  &
            cc%bsw+cc%bHW,cc%nsc,cc%NSN,                 &
            treeG,fseed,fleaf,froot,fwood,               &
            cc%annualGPP,cc%annualNPP,cc%annualTrsp,     &
            cc%annualNup,cc%annualfixedN,cc%mu,          &
            cc%Atrunk,cc%Asap,cc%Ktrunk,                 &
            cc%treeHU,cc%treeW0,(cc%farea(j),j=1,Ysw_max)

        ! Screen output
        write(*,'(1(I6,","),2(I4,","),30(F9.3,","))') &
          cc%ccID,cc%species,cc%layer,                &
          cc%nindivs*10000, cc%layerfrac,             &
          dDBH,cc%dbh,cc%height,cc%Acrown,            &
          cc%nsc,cc%annualGPP,cc%mu,cc%Atrunk,        &
          cc%Asap,cc%Ktrunk,cc%treeHU,cc%treeW0
    enddo

    ! tile pools output
    !do i = 1, vegn%n_cohorts
    !    cc => vegn%cohorts(i)
    !    vegn%annualfixedN  = vegn%annualfixedN  + cc%annualfixedN * cc%nindivs
    !enddo
    plantC = vegn%NSC + vegn%SeedC + vegn%leafC + vegn%rootC +   &
            vegn%SapwoodC + vegn%woodC
    plantN = vegn%NSN + vegn%SeedN + vegn%leafN +                &
            vegn%rootN + vegn%SapwoodN + vegn%woodN
    soilC  = sum(vegn%SOC(:))
    soilN  = sum(vegn%SON(:)) + vegn%mineralN
    vegn%totN = plantN + soilN
    write(f2,'(1(I5,","),30(F12.4,","),6(F12.4,","),30(F12.4,","))') &
        iyears,       &
        vegn%CAI,vegn%LAI, vegn%treecover, vegn%grasscover, &
        vegn%annualGPP, vegn%annualResp, vegn%annualRh, vegn%C_combusted, &
        vegn%annualPrcp, vegn%SoilWater,vegn%annualTrsp, vegn%annualEvap, vegn%annualRoff, &
        plantC,soilC,plantN *1000, soilN * 1000, (plantN+soilN)*1000,&
        vegn%NSC, vegn%SeedC, vegn%leafC, vegn%rootC,  &
        vegn%SapwoodC, vegn%woodC,                     &
        vegn%NSN*1000, vegn%SeedN*1000, vegn%leafN*1000, vegn%rootN*1000, &
        vegn%SapwoodN *1000,  vegn%WoodN *1000,  &
        (vegn%SOC(j),j=1,5), (vegn%SON(j)*1000,j=1,5), &
        vegn%mineralN*1000,   vegn%annualfixedN*1000, vegn%annualNup*1000, &
        vegn%annualN*1000,vegn%N_P2S_yr*1000, vegn%Nloss_yr*1000

 end subroutine annual_diagnostics

!========================================================================
! read in forcing data (Users need to write their own data input procedure)
subroutine read_FACEforcing(forcingData,datalines,days_data,yr_data,timestep)
   type(climate_data_type),pointer,intent(inout) :: forcingData(:)
   integer,intent(inout) :: datalines,days_data,yr_data
   real, intent(inout)   :: timestep
   !------------local var -------------------
   type(climate_data_type), pointer :: climateData(:)
   character(len=80)  commts
   integer, parameter :: niterms=9       ! MDK data for Oak Ridge input

   integer, allocatable :: doy_data(:),year_data(:)
   real,    allocatable :: hour_data(:),input_data(:,:)
   real    :: hr, clim(niterms)
   integer :: yr,dy
   integer :: istat1,istat2,istat3
   integer :: ndays,nyear,totlines
   integer :: m,n,i

   ! Open forcing data
   climfile=trim(filepath_in)//trim(climfile)
   open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
   write(*,*)'istat2',istat2

   ! Skip 1 line of input met data file
   read(11,'(a160)') commts ! MDK data only has one line comments
   ! Count total lines
   totlines = 0  ! to record the lines in a file
   do
     read(11,*,IOSTAT=istat3)yr,dy !,hr,(clim(n),n=1,niterms)
     if(istat3 < 0)exit
     totlines = totlines + 1
   enddo ! end of reading the forcing file
   write(*,*)'total lines:',totlines

   ! Allocate arrays for reading in data
   allocate(doy_data(totlines),year_data(totlines),hour_data(totlines))
   allocate(input_data(niterms,totlines))

   ! Read forcing files
   rewind 11
   read(11,'(a160)') commts
   ndays = 0 ! the total days in this data file
   nyear = 0 ! the total years of this data file
   dy    = -1  ! Initial value
   yr    = -1
   do m = 1, totlines
     read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
                         (input_data(n,m),n=1,niterms)
     ! Count days
     if(m > 1) then
       dy = doy_data(m-1)
       yr = year_data(m-1)
     endif
     if(dy /= doy_data(m)) ndays = ndays + 1
     if(yr /= year_data(m))nyear = nyear + 1
   enddo
   ! Check fast time step
   timestep = hour_data(2) - hour_data(1)
   if (timestep==1.0)then
       write(*,*)"the data freqency is hourly"
   elseif(timestep==0.5)then
       write(*,*)"the data freqency is half hourly"
   else
       write(*,*)'hour data:',hour_data(1),hour_data(2),hour_data(3)
       write(*,*)"Please check time step!"
       stop
   endif

   ! Put the data into forcing
   allocate(climateData(totlines))
   do i=1,totlines
      climateData(i)%year      = year_data(i)          ! Year
      climateData(i)%doy       = doy_data(i)           ! day of the year
      climateData(i)%hod       = hour_data(i)          ! hour of the day
      climateData(i)%PAR       = input_data(1,i)       ! umol/m2/s
      climateData(i)%radiation = input_data(2,i)       ! W/m2
      climateData(i)%Tair      = input_data(3,i) + 273.16  ! air temperature, K
      climateData(i)%Tsoil     = input_data(4,i) + 273.16  ! soil temperature, K
      climateData(i)%RH        = input_data(5,i) * 0.01    ! relative humidity (0.xx)
      climateData(i)%rain      = input_data(6,i)/(timestep * 3600)! ! kgH2O m-2 s-1
      climateData(i)%windU     = input_data(7,i)        ! wind velocity (m s-1)
      climateData(i)%P_air     = input_data(8,i)        ! pa
      climateData(i)%CO2       = input_data(9,i) * 1.0e-6       ! mol/mol
      climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
   enddo
   forcingData => climateData
   datalines = totlines
   days_data = ndays
   yr_data   = nyear
   write(*,*)"forcing", datalines,days_data,yr_data

   !Close opened file and release memory
   close(11)    ! close forcing file
   deallocate(doy_data,year_data,hour_data,input_data)
end subroutine read_FACEforcing

 !=============================================================
 ! for reading in NACP site synthesis forcing
subroutine read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
   type(climate_data_type),pointer,intent(inout) :: forcingData(:)
   integer,intent(inout) :: datalines,days_data,yr_data
   real, intent(inout)   :: timestep
   !------------local var -------------------
   type(climate_data_type), pointer :: climateData(:)
   character(len=80)  commts
   integer, parameter :: niterms=15       ! NACP site forcing
   integer, parameter :: ilines=22*366*48 ! the maxmum records
   integer,dimension(ilines) :: year_data
   real,   dimension(ilines) :: doy_data,hour_data
   real input_data(niterms,ilines)
   real inputstep
   integer :: istat1,istat2,istat3
   integer :: doy,idays
   integer :: i,j,k
   integer :: m,n

   climfile=trim(filepath_in)//trim(climfile)
   write(*,*)'inputfile: ',climfile
   ! open forcing data
   open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
   write(*,*)istat2
   ! skip 2 lines of input met data file
   read(11,'(a160)') commts
   read(11,'(a160)') commts
   m       = 0  ! to record the lines in a file
   idays   = 1  ! the total days in a data file
   yr_data = 0 ! to record years of a dataset
   do    ! read forcing files
       m=m+1
       read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
                               (input_data(n,m),n=1,niterms)
       if(istat3<0)exit
       if(m == 1) then
           doy = doy_data(m)
       else
           doy = doy_data(m-1)
       endif
       if(doy /= doy_data(m)) idays = idays + 1
       !write(*,*)year_data(m),doy_data(m),hour_data(m)
       ! discard one line
       !read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
       !                        (input_data(n,m),n=1,niterms)
   enddo ! end of reading the forcing file

   timestep = hour_data(2) - hour_data(1)
   write(*,*)"forcing",datalines,yr_data,timestep,dt_fast_yr
   if (timestep==1.0)then
       write(*,*)"the data freqency is hourly"
   elseif(timestep==0.5)then
       write(*,*)"the data freqency is half hourly"
   else
       write(*,*)"Please check time step!"
       stop
   endif
   close(11)    ! close forcing file
   ! Put the data into forcing
   datalines = m - 1
   days_data = idays
   yr_data  = year_data(datalines-1) - year_data(1) + 1

   allocate(climateData(datalines))
   do i=1,datalines
      climateData(i)%year      = year_data(i)          ! Year
      climateData(i)%doy       = doy_data(i)           ! day of the year
      climateData(i)%hod       = hour_data(i)          ! hour of the day
      climateData(i)%PAR       = input_data(11,i)*2.0  ! umol/m2/s
      climateData(i)%radiation = input_data(11,i)      ! W/m2
      climateData(i)%Tair      = input_data(1,i)       ! air temperature, K
      climateData(i)%Tsoil     = input_data(1,i)       ! soil temperature, K
      climateData(i)%rain      = input_data(7,i)       ! kgH2O m-2 s-1
      climateData(i)%windU     = input_data(5,i)        ! wind velocity (m s-1)
      climateData(i)%P_air     = input_data(9,i)        ! pa
      climateData(i)%RH        = input_data(3,i)/mol_h2o*mol_air* & ! relative humidity (0.xx)
                                 climateData(i)%P_air/esat(climateData(i)%Tair-273.16)
      climateData(i)%CO2       = input_data(15,i) * 1.0e-6       ! mol/mol
      climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
   enddo
   forcingData => climateData
   write(*,*)"forcing", datalines,days_data,yr_data

end subroutine read_NACPforcing

!=========== Write output file header ====================
subroutine set_up_output_files(runID,fpath,fno1,fno2,fno3,fno4,fno5,fno6)
   character(len=50),intent(in):: runID,fpath
   integer,intent(inout):: fno1,fno2,fno3,fno4,fno5,fno6

   ! ----------Local vars ------------
   character(len=150) :: YearlyCohort,DailyCohort,HourlyCohort ! Output file names
   character(len=150) :: YearlyPatch, DailyPatch, HourlyPatch  ! output file names
   character(len=50)  :: filesuffix
   integer :: istat1, istat2, istat3

    ! File path and names
    filesuffix   = '_'//trim(runID)//'.csv' ! tag for simulation experiments
    HourlyCohort = trim(fpath)//'Cohort_hourly'//trim(filesuffix)   ! hourly
    HourlyPatch  = trim(fpath)//'Patch_hourly'//trim(filesuffix)    ! hourly
    DailyCohort  = trim(fpath)//'Cohort_daily'//trim(filesuffix)    ! daily
    DailyPatch   = trim(fpath)//'Patch_daily'//trim(filesuffix)     ! Daily
    YearlyCohort = trim(fpath)//'Cohort_yearly'//trim(filesuffix)   ! Yearly
    YearlyPatch  = trim(fpath)//'Patch_yearly'//trim(filesuffix)    ! Yearly

    ! Open files
    fno1=91; fno2=92; fno3=103; fno4=104; fno5=105; fno6=106
    open(fno1,file=trim(HourlyCohort),ACTION='write', IOSTAT=istat1)
    open(fno2,file=trim(HourlyPatch), ACTION='write', IOSTAT=istat1)
    open(fno3,file=trim(DailyCohort), ACTION='write', IOSTAT=istat2)
    open(fno4,file=trim(DailyPatch),  ACTION='write', IOSTAT=istat2)
    open(fno5,file=trim(YearlyCohort),ACTION='write', IOSTAT=istat3)
    open(fno6,file=trim(YearlyPatch), ACTION='write', IOSTAT=istat3)

    ! Add file header
    write(fno1,'(5(a8,","),30(a12,","))')        &       ! Hourly cohort
         'year','doy','hour','cID','sp','layer', &
         'density','dbh','height','Acrown',      &
         'bl','LAI','GPP', 'NPP', 'Transp',      &
         'Psi_L','Psi_W','W_leaf','W_stem'

    write(fno2,'(5(a8,","),30(a12,","))')      &       ! Hourly tile
         'year','doy','hour','rad',            &
         'Tair','Prcp', 'GPP', 'Resp',         &
         'Transp','Evap','Runoff','Soilwater', &
         'wcl', 'psi_soil','k_soil',           &
         'bl','Psi_L','Psi_W','W_leaf','W_stem','Transp'

    write(fno3,'(9(a6,","),45(a8,","))')               &  ! Daily cohort
         'year','doy','hour','cID','PFT','layer',      &
         'Pheno','ndm','ncd','density','flayer','LAI', &
         'gpp','resp','transp','NPPL','NPPR','NPPW',   &
         'W_LF','W_SW','W_HW',                         &
         'NSC','seedC','leafC','rootC','SW-C','HW-C',  &
         'NSN','seedN','leafN','rootN','SW-N','HW-N',  &
         'GDD','ALT'

    write(fno4,'(2(a5,","),55(a10,","))')  'year','doy',    &  ! Daily tile
         'Tc','Prcp', 'totWs',  'Trsp', 'Evap','Runoff',    &
         'ws1','ws2','ws3', 'LAI','GPP', 'Rauto', 'Rh',     &
         'W_LF','W_SW','W_HW',                              &
         'NSC','seedC','leafC','rootC','SW-C','HW-C',       &
         'NSN','seedN','leafN','rootN','SW-N','HW-N',       &
         'fineL', 'strucL', 'McrbC', 'fastSOC', 'slowSOC',  &
         'fineN', 'strucN', 'McrbN', 'fastSON', 'slowSON',  &
         'mineralN', 'N_uptk','Kappa'

    write(fno5,'(3(a5,","),40(a9,","))')            &    ! Yearly cohort
         'yr','cID','PFT','layer','density','f_L',  &
         'dDBH','dbh','height','Acrown',            &
         'wood','nsc', 'NSN','NPPtr','seed',        &
         'NPPL','NPPR','NPPW','GPP','NPP','Transp', &
         'N_uptk','N_fix','mu','Atrunk','Asap',     &
         'Ktree','treeHU','treeW0',                 &
         'farea1','farea2','farea3','farea4','farea5'

    write(fno6,'(1(a5,","),80(a12,","))')  'year',         &  ! Yearly tile
         'CAI','LAI','treecover', 'grasscover',            &
         'GPP', 'Rauto', 'Rh', 'burned',                   &
         'rain','SiolWater','Transp','Evap','Runoff',      &
         'plantC', 'soilC', 'plantN', 'soilN', 'totN',     &
         'NSC', 'SeedC', 'leafC', 'rootC', 'swC', 'hwC',   &
         'NSN', 'SeedN', 'leafN', 'rootN', 'swN', 'hwN',   &
         'fineL', 'strucL', 'McrbC', 'fastSOC', 'slowSOC', &
         'fineN', 'strucN', 'McrbN', 'fastSON', 'slowSON', &
         'mineralN', 'N_fxed','N_uptk','N_yrMin','N_P2S',  &
         'N_loss','seedC','seedN','Seedling-C','Seedling-N'

end subroutine set_up_output_files

!===========for netcdf IO ============================
#ifdef USE_NETCDF

 !=====================================================
   subroutine nc_read_3D(FILE_NAME,field_idx,NX,NY,Ntime,DA)
   ! This is the name of the data file we will create.
   character (len = *), intent(in) :: FILE_NAME,field_idx
   integer, intent(in) :: NX, NY, Ntime
   real, intent(inout) :: DA(:,:,:)
   !------local vars ----------------
   ! When we create netCDF files, variables and dimensions, we get back
   ! an ID for each one.
   integer :: ncid, varid

   ! Loop indexes, and error handling.
   integer :: x, y, t

   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
   print *, 'ncid=',ncid
   ! Get the varid of the data variable, based on its name.
   call check( nf90_inq_varid(ncid, trim(field_idx), varid) )
   print *, 'varid=',varid
   ! Read the data.
   call check( nf90_get_var(ncid, varid, DA) )

   ! Check the data.
   !do x = 1, NX
   !   do y = 1, NY
   !      print *, "DA(", x, ", ", y, ") = ", DA(x,y,1)
   !   end do
   !end do

   ! Close the file, freeing all resources.
   call check( nf90_close(ncid) )

   print *,"*** SUCCESS reading example file ", FILE_NAME, "! "

   end subroutine nc_read_3D

 !=================================================================
   subroutine check(status)
     integer, intent ( in) :: status

     if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
     end if
   end subroutine check
 !=================================================================
#endif

!================================================
end module io_mod
