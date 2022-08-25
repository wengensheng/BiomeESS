module global_vars
  integer, parameter :: maxLayers = 3
  integer, parameter :: maxPFTs  = 5
  integer, parameter :: max_cc = 150
  integer, parameter :: n_cc = 15
  integer, parameter :: maxyears= 3600
  integer, parameter :: bins=21 ! 6
  integer, parameter :: bins2 = bins + 1
  real, parameter    ::  PI = 3.1415926
end module

!===========Main============
 program growthana

 implicit none

 integer :: istat1,istat2,istat3,istat4
 character(len=10) :: runID
 character(len=90) :: filepath
 character(len=90) :: filein1,filein2,fileout1,fileout2,fileout3
 integer :: fin_1,fin_2,fout_1,fout_2,fout_3

 runID    = 'BCI_hydro' !'OR_LMA' ! 'OR_Nfix' !'OR_GAPLUE'
 filepath = 'output/' ! 'testruns/'

 ! Yearly cohort
 filein1  = trim(filepath)//'Cohort_yearly'//'_'//trim(runID)//'.csv'
 filein2  = trim(filepath)//'Patch_yearly'//'_'//trim(runID)//'.csv'
 fileout1 = trim(filepath)//'FirstCohortYearly'//'_'//trim(runID)//'.csv'
 fileout2 = trim(filepath)//'Summary_'//trim(runID)//'.csv'

 fin_1=211; fin_2=212
 fout_1=15; fout_2=14

 open(fin_1,file=filein1,status='old',ACTION = 'read', IOSTAT=istat1)
 open(fin_2,file=filein2,status='old',ACTION = 'read', IOSTAT=istat2)

 ! output to files
 open(fout_1,file=fileout1)
 open(fout_2,file=fileout2)

 if(istat1==0 .and. istat2==0)then
   write(*,*)"Yearly Files opened!"
   call cohort_yearly(fin_1,fin_2,fout_1,fout_2)
 else
   write(*,*)"Files open failed",istat1, istat2
 endif

 close (fout_1); close (fout_2)
 close (fin_1);  close (fin_2)

 ! -------Hourly cohort----------
 filein1  = trim(filepath)//'Cohort_hourly'//'_'//trim(runID)//'.csv'
 fileout1 = trim(filepath)//'FirstCohortHourly'//'_'//trim(runID)//'.csv'
 fin_1=213; fout_1=16
 open(fin_1,file=filein1,status='old',ACTION = 'read', IOSTAT=istat1)
 ! output to files
 open(fout_1,file=fileout1)
 if(istat1==0 )then
   write(*,*)"Hourly Files opened!"
   call cohort_hourly(fin_1,fout_1)
 else
   write(*,*)"Files open failed",istat1, istat2
 endif
 close (fout_1)
 close (fin_1)

 ! -------Daily cohort----------
 filein1  = trim(filepath)//'Cohort_daily'//'_'//trim(runID)//'.csv'
 fileout1 = trim(filepath)//'FirstCohortDaily'//'_'//trim(runID)//'.csv'
 fin_1=213; fout_1=16
 open(fin_1,file=filein1,status='old',ACTION = 'read', IOSTAT=istat1)
 ! output to files
 open(fout_1,file=fileout1)
 if(istat1==0 )then
   write(*,*)"Daily Files opened!"
   call cohort_daily(fin_1,fout_1)
 else
   write(*,*)"Files open failed",istat1, istat2
 endif
 close (fout_1)
 close (fin_1)

end

!========================================================================

 subroutine cohort_yearly(fin_1,fin_2,fout_1,fout_2)
   use global_vars
   implicit none
   integer, intent(in):: fin_1,fin_2,fout_1,fout_2

   !------local vars -----------
   integer, parameter :: Columns = 28 + 210 ! 21
   integer, parameter :: N_veg_out = 22

   character(len=180) :: comments,header1
   character(len=8 ) :: BA_Header(maxPFTs), dD_Header(maxPFTs),farea_Header(210)
   integer,dimension(2,maxyears):: yearch
   real,   dimension(Columns,max_cc,maxyears):: dataarray,dataarray2
   real,   dimension(N_veg_out,maxPFTs,maxyears):: firstCC,L2firstCC
   real,   dimension(maxPFTs,maxyears):: BApft,dDBH,density
   real,   dimension(45,maxyears):: ecodata
   real,   dimension(maxyears):: HTstar,Layer1BM,Layer1Den,CApft,totalBM,basalA
   real,   dimension(maxyears):: BMmort,totCAI
   real,   dimension(N_veg_out,maxPFTs,maxLayers,maxyears):: meanvalues
   real,   dimension(bins2) :: DBHbins
   real,   dimension(bins,maxyears) :: DBHclasses,BMclasses

   integer :: totyears,m,n,i,j,k,iPFT,iLayer,layertype,yr
   integer :: commentlines
   integer :: istat2,istat3,cc
   real :: DBH, GPP, NPP, plantC, soilC,plantN, soilN, mineralN, Nmin

   ! Set file header
   do i=1,maxPFTs
     write (BA_Header(i), "(A6,I1)") "BA_PFT", i
     write (dD_Header(i), "(A6,I1)") "dD_PFT", i
   enddo
   header1 = 'yr,ID,PFT,L,dst,f_L,dD,dbh,HT,CA,wood,nsc,NSN,NPP,' &
           // 'seed,NPPL,NPPR,NPPW,GPP,NPP,Trsp,N_up,N_fix,mu,'  &
           // 'Atrunk,Asap,Ktree,treeHU,treeW0,farea1,farea2'
   !do i=1,Columns - 24
   !   write (farea_Header(i), "(A6,I3)") "farea_", i
   !enddo

   write(fout_1,'(a150,",")') header1 !, (farea_Header(i),i=1,Columns - 24)
   write(fout_2,'(360(a12,","))') 'year',   &
          (BA_Header(i),i=1,maxPFTs), (dD_Header(i),i=1,maxPFTs), &
          'CAI','critH','GPP', 'NPP','plantC','soilC','plantN','soilN', &
          'mineralN','Nmin','DBH1','DBH2','DBH3','DBH4','DBH5','DBH6'


  ! initiate bins
  DBHbins(1) = 0.0 ! 0.05
  do i=2,bins
      DBHbins(i) = DBHbins(i-1)+0.1 !0.1
  enddo
  DBHbins(1) = 0.025
  DBHbins(bins+1) = 9999.
  DBHclasses = 0.0
  BMclasses = 0.0

  ! read in file
  commentlines = 1 ! 91
  do i=1,commentlines !91
     read(fin_1,*,IOSTAT=istat3)comments ! the first line, header
  enddo

  firstCC = -9999.
  L2firstCC = -9999.
  meanvalues = 0.
  HTstar     = 0.
  Layer1BM   = 0.
  Layer1Den  = 0.
  CApft      = 0.
  totCAI     = 0.
  totalBM    = 0.
  BMmort     = 0.
  basalA     = 0.
  BApft      = 0.0
  dDBH       = 0.
  density    = 0.
  m=0
  do
     m=m+1 ! Year
     read(fin_1,*,IOSTAT=istat3)yearch(1,m),yearch(2,m)
     write(*,*)yearch(1,m),yearch(2,m)
     if(istat3<0 .or. yearch(2,m)<1)exit
     cc = yearch(2,m) ! total cohorts in this year
     do j=1,cc
         read(fin_1,*,IOSTAT=istat3)yr,(dataarray(i,j,m),i=1,Columns)
         if(istat3<0 .or. yr < 1)exit
         if(cc>4 .and. i==cc)cycle ! skip the last cohort

!             CA of plants above 5 meters
         if(dataarray(8,j,m) > 0.5) then ! height
             CApft(m) = CApft(m) + dataarray(9,j,m)* dataarray(4,j,m)/10000.0 ! crown area
         endif
!             first layer boimass and HT*
         if(dataarray(3,j,m) < 2.) then
             HTstar(m)    = dataarray(8,j,m)
             Layer1Den(m) = Layer1Den(m) + dataarray(4,j,m)
             Layer1BM(m)  = Layer1BM(m)  + dataarray(4,j,m) * &
                           (dataarray(10,j,m)+dataarray(11,j,m))/10000
         endif
!             total biomass
         totalBM(m) = totalBM(m) + (dataarray(11,j,m)+dataarray(10,j,m))  &
                    * dataarray(4,j,m)/10000
         totCAI(m)  = totCAI(m) + dataarray(9,j,m) * dataarray(4,j,m)/10000.0
!             woody biomass residence time
         DBH = dataarray(7,j,m)
         iLayer=dataarray(3,j,m)
         iPFT  =dataarray(2,j,m)-1

!             calculate means or sums
         dataarray2(:,j,m) = dataarray(:,j,m)
         dataarray2(6:Columns,j,m) = dataarray(6:Columns,j,m) * dataarray(4,j,m)

!             size classes distribution, including total individuals; biomass distribution
         do i=1,bins
             if(dataarray(7,j,m)>DBHbins(i) .and. dataarray(7,j,m)<DBHbins(i+1) ) then
                DBHclasses(i,m) = DBHclasses(i,m) + dataarray(4,j,m)
                BMclasses(i,m)  = BMclasses(i,m)  + dataarray(4,j,m)*         &
                                  (dataarray(11,j,m)+dataarray(10,j,m)) /10000
              endif
         enddo

         iLayer = dataarray(3,j,m)
         iPFT   = dataarray(2,j,m)  ! dataarray(2,j,m)-2
         ! Basal area with DBH > 0.1 m
         if(dataarray(7,j,m)> 0.0)then
             basalA(m)   = basalA(m)   + 0.25*PI*dataarray(7,j,m)**2 * dataarray(4,j,m)
             BApft(iPFT,m) = BApft(iPFT,m) + 0.25*PI*dataarray(7,j,m)**2 * dataarray(4,j,m)
         endif
         ! mean dDBH of the first layer
         if(dataarray(3,j,m)<2.0)then
             dDBH(iPFT,m) = dDBH(iPFT,m) + dataarray(6,j,m) * dataarray(4,j,m)
             density(iPFT,m) = density(iPFT,m) + dataarray(4,j,m)
         endif
     enddo ! j, cohorts
     do i=1,maxPFTs
        if(density(i,m)>0.0) dDBH(i,m)=dDBH(i,m)/density(i,m)
     enddo

     ! Output the first cohort
     write(fout_1,105)yearch(1,m),dataarray(:,1,m)

  enddo ! End of calculating
  totyears=m-1

   ! read in another file for combining data
   read(fin_2,*,IOSTAT=istat3)comments
   do m = 1, totyears
       read(fin_2,*,IOSTAT=istat3)n,(ecodata(i,m),i=1,45)
       GPP      = ecodata(5,m)
       NPP      = ecodata(5,m)-ecodata(6,m)
       plantC   = ecodata(14,m)
       soilC    = ecodata(15,m)
       plantN   = ecodata(16,m)
       soilN    = ecodata(17,m)
       mineralN = ecodata(41,m)
       Nmin     = ecodata(44,m)
       !if(m>500) &
       write(fout_2,140)m,(BApft(i,m),i=1,maxPFTs), (dDBH(i,m), i=1,maxPFTs), &
                  totCAI(m), HTstar(m), &
                  GPP, NPP, plantC, soilC,plantN, soilN, mineralN, Nmin, &
                  (DBHclasses(i,m),i=1,bins)
   enddo
   write(*,*)"Done!"

   140 format(1(I8,','),10(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
   105 format(1(I8,','),3(f6.1,','),300(E15.4,','))
   103 format(1(I8,','),11(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
   102 format(1(I8,','),4(f15.4,','),6(f15.4,','), 2(I8,',',10(f15.4,',')), 12(2(I8,','),10(f15.4,',')) )
   101 format(1(I8,','),4(f15.4,','),6(f15.4,','), 1(I8,',',10(f15.4,',')),  6(2(I8,','),10(f15.4,',')) )

end subroutine cohort_yearly

!================Cohort hourly ===========
subroutine cohort_hourly(fin_1,fout_1)
  use global_vars
  implicit none
  integer, intent(in):: fin_1,fout_1

  !------local vars -----------
  integer, parameter :: Columns = 19 ! Hourly output
  character(len=180) :: comments,header1
  character(len=8 ) :: BA_Header(maxPFTs), dD_Header(maxPFTs)
  real,dimension(Columns,max_cc):: dataarray

  integer :: m, j, i, yr
  integer :: cc, istat3, commentlines

  header1 ='yr,doy,h,cID,sp,layer,density,dbh,height,Acrown,bl,LAI,GPP,NPP,' &
          // 'Transp,Psi_L,Psi_W,W_leaf,W_stem'
  write(fout_1,'(a100,",")') trim(header1)
  ! read in file
  commentlines = 1
  do i=1,commentlines
     read(fin_1,*,IOSTAT=istat3)comments ! the first line, header
  enddo

  m=0
  do while(istat3==0)
     m=m+1 ! Year
     read(fin_1,*,IOSTAT=istat3)cc
     write(*,*)"total cohorts:",cc
     if(istat3<0 .or. cc<1)exit
     do j=1,cc
         read(fin_1,*,IOSTAT=istat3)(dataarray(i,j),i=1,Columns)
         if(istat3<0)exit
     enddo
     ! Output the first cohort
     write(fout_1,105)dataarray(:,1)
     write(*,*)"Writing hourly first cohort",m

  enddo ! End of file reading

  105 format(6(f8.0,','),30(f15.4,','))
end subroutine cohort_hourly

!================Cohort daily ===========
subroutine cohort_daily(fin_1,fout_1)
  use global_vars
  implicit none
  integer, intent(in):: fin_1,fout_1

  !------local vars -----------
  integer, parameter :: Columns = 33 ! Daily output
  character(len=180) :: comments,header1
  character(len=8 ) :: BA_Header(maxPFTs), dD_Header(maxPFTs)
  real,dimension(Columns,max_cc):: dataarray

  integer :: m, j, i, yr,doy
  integer :: cc, istat3, commentlines

  header1 ='yr,doy,hour,cID,PFT,layer,Pheno,ndm,ncd,density,f_L,LAI,' &
        // 'gpp,resp,transp,NPPL,NPPR,NPPW,W_LF,W_SW,W_HW,NSC,seedC,' &
        // 'leafC,rootC,SW-C,HW-C,NSN,seedN,leafN,rootN,SW-N,HW-N'
  write(fout_1,'(a170,",")') trim(header1)
  ! read in file
  commentlines = 1
  do i=1,commentlines
     read(fin_1,*,IOSTAT=istat3)comments ! the first line, header
  enddo

  m=0
  do while(istat3==0)
     m=m+1 ! Year
     read(fin_1,*,IOSTAT=istat3)yr,doy,cc
     write(*,*)"yr,doy,cc:",yr,doy,cc
     if(istat3<0 .or. cc<1)exit
     do j=1,cc
         read(fin_1,*,IOSTAT=istat3)(dataarray(i,j),i=1,Columns)
         if(istat3<0)exit
     enddo
     ! Output the first cohort
     write(fout_1,105)dataarray(:,1)
     write(*,*)"Writing daily first cohort",m

  enddo ! End of file reading

  105 format(8(f8.0,','),30(f15.4,','))
end subroutine cohort_daily

!=====================
