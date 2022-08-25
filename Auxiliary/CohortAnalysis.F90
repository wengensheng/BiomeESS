 program growthana

 implicit none
 integer, parameter :: no_file=1
 integer, parameter :: items = 22
 integer, parameter :: planttraits = 18
 integer, parameter :: maxLayers =2
 integer, parameter :: maxPFTs  =10
 integer, parameter :: id_format = 103
 integer, parameter :: cohorts = 120
 integer, parameter :: maxyears=1800
 integer, parameter :: bins=10 ! 6
 integer, parameter :: bins2 = bins + 1 
 real,    parameter :: PI = 3.1415926
 character(len=90) :: filein(no_file)
 character(len=60) :: filepath,filesuffix,chaSOM(12)
 character(len=50) :: parameterfile(10), runtag
 character(len=160) :: filein1,filein2,fileout1,fileout2,fileout3
 character(len=360) :: comments,format_string
 character(len=12) :: head(10,20),numPFTs
 integer,dimension(2,maxyears):: yearch
 real,   dimension(items,cohorts,maxyears):: dataarray
 real,   dimension(planttraits,maxPFTs,maxyears):: firstCC,L2firstCC
 real,   dimension(maxPFTs,maxyears):: BApft,CApft,dDBH,density,densityL1,maxLAI
 real,   dimension(maxPFTs,maxyears):: NPPlf,NPPwd,NPPfr,NPPsd
 real,   dimension(45,maxyears):: ecodata
 real,   dimension(maxyears):: HTstar,Layer1BM,Layer1Den,CAabv5m,totalBM,basalA
 real,   dimension(maxyears):: BMmort,totCAI
 real,   dimension(maxyears,10):: siteGPP,siteNPP,siteNPPlf,siteNPPwd,siteNPPfr,siteNPPsd
 real,   dimension(maxyears,10):: siteBA,siteCAI,siteHTs,siteNSC,siteNSN
 real,   dimension(maxyears,10):: siteplantC, sitesoilC,siteplantN,siteSoilN,siteMineralN
 real,   dimension(maxyears,10,10)::compBA
 real,   dimension(planttraits,maxPFTs,maxLayers,maxyears):: meanvalues 

 real,   dimension(bins2) :: DBHbins
 real,   dimension(bins,maxyears) :: DBHclasses,BMclasses

 integer :: totyears,m,n,i,j,k,iPFT,iLayer,layertype
 integer :: siteplus, iSOM, iTests, totTests, initialPFTs, SOMlevels
 integer :: commentlines
 integer :: istat2,istat3,cc
 real :: DBH, GPP, NPP, plantC, soilC,plantN, soilN, mineralN, Nmin

 ! Input files
 !  filepath = 'output/savedAllocation6PFTs/' ! 'testruns/'
   filepath = 'output/rerun0607/' ! 'output/lowDCRate/' ! 'testruns/'
   SOMlevels = 8
   initialPFTs   =  8 ! 8
   siteplus = 9 ! 9 ! 7
   totTests = siteplus - initialPFTs

   runtag = '_aCO2_' ! '_eCO2_' !   '_fN2FixA_' ! '_LoosedA_' !  fNSNmax =2, LMA =0.14
   chaSOM = (/'SC04','SC06','SC08','SC10','SC12','SC14', &
              'SC16','SC18','SC20','SC22','SC24','SC26'/)
   !chaSOM = (/'SC06','SC09','SC12','SC15','SC18', 'SC21', &
   !           'SC24','SC27','SC30','SC33','SC36','SC39'/)

   if(initialPFTs == 3)then
      parameterfile = (/'phiRL1-3','phiRL2-4','phiRL3-5','phiRL4-6','phiRL5-7', &
                        'phiRL6-8','phiRL7-9','phiRL8-1','phiRL9-1','_decdu33' /)
      numPFTs = '_Three'

   elseif(initialPFTs == 2)then
      parameterfile = (/'phiRL1-2','phiRL2-3','phiRL3-4','phiRL4-5','phiRL5-6', &
                        'phiRL6-7','phiRL7-8','phiRL8-9','phiRL9-1','_decdu33' /)
      numPFTs = '_Two'
   elseif(initialPFTs == 1)then
      parameterfile = (/'phiRL1','phiRL2','phiRL3','phiRL4','phiRL5', &
                        'phiRL6','phiRL7','phiRL8','phiRL9','decdu3' /)
      numPFTs = '_One'
    elseif(initialPFTs > 4)then
      parameterfile = 'All'
      numPFTs = '_8PFTs'
   endif


!!  for column names
!    do j=1, 13
!    j=1
!    do i=1, 10
!        if (i < 10) then
!            format_string = "(A5,I1)"
!        else
!            format_string = "(A5,I2)"
!       endif
!        write (clm_name(i,j),format_string) "BA-", i
!        print *, trim(clm_name(i,j))
!    enddo ! i, different tests
!    enddo ! j, different variables




do iSOM = 1,SOMlevels ! 2, 12, 3
       siteGPP      = 0.0
       siteNPP      = 0.0
       siteNPPlf    = 0.0
       siteNPPwd    = 0.0
       siteNPPfr    = 0.0
       siteNPPsd    = 0.0
       siteBA       = 0.0
       siteplantC   = 0.0
       sitesoilC    = 0.0
       siteplantN   = 0.0
       siteSoilN    = 0.0
       siteMineralN = 0.0
       compBA = 0.0

do iTests =1, totTests
 write(*,*)iSOM,iTests
 filesuffix  = trim(chaSOM(iSOM))//trim(runtag)//trim(parameterfile(iTests))//'.csv'

 ! Analysis starts here
 commentlines = 1 ! 91 ! 91
 filein1  = trim(filepath)//'Annual_cohorts'//trim(filesuffix)
 filein2  = trim(filepath)//'Ecosystem_yearly'//trim(filesuffix)
 write(*,*)filein1
 write(*,*)filein2

 fileout2 = trim(filepath)//'Sum'//trim(runtag)//trim(chaSOM(iSOM))//trim(numPFTs)//'.csv'

 ! Open input files
 open(211,file=filein1,status='old',ACTION='read',IOSTAT=istat2)
 open(212,file=filein2,status='old',ACTION='read',IOSTAT=istat3) ! for combining data

! Open output files
 open(14,file=fileout2)
 !open(13,file=fileout1)


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
 if(istat2==0)then
       write(*,*)" open successfully!"
       do i=1,commentlines !91
          read(211,*,IOSTAT=istat3)comments
       enddo

       firstCC = -9999.
       L2firstCC = -9999.
       meanvalues = 0.
       HTstar     = 0. 
       Layer1BM   = 0.
       Layer1Den  = 0.
       CAabv5m    = 0.
       totCAI     = 0.
       totalBM    = 0.
       BMmort     = 0.
       basalA     = 0.
       BApft      = 0.0
       CApft      = 0.0
       NPPlf   = 0.0
       NPPwd   = 0.0
       NPPfr   = 0.0
       NPPsd   = 0.0
       maxLAI     = 0.
       dDBH       = 0.
       density    = 0.
       densityL1  = 0.
       m=0
       do 
          m=m+1

          read(211,*,IOSTAT=istat3)yearch(1,m),yearch(2,m)
          !write(*,*)yearch(1,m),yearch(2,m)
          if(istat3<0)exit
          cc = yearch(2,m)
          if(yearch(2,m) > 4) cc = yearch(2,m) -1 ! skip the last cohort
          do j=1,cc
              read(211,*,IOSTAT=istat3)(dataarray(i,j,m),i=1,items)
              if(istat3<0)exit
!             CA of plants above 5 meters
              if(dataarray(8,j,m) > 0.5) then ! height
                  CAabv5m(m) = CAabv5m(m) + dataarray(9,j,m)* dataarray(4,j,m)/10000.0 ! crown area
              endif
!             first layer boimass and HT*
              if(dataarray(3,j,m) < 2.) then
                  HTstar(m)    = dataarray(8,j,m)
                  Layer1Den(m) = Layer1Den(m) + dataarray(4,j,m)
                  Layer1BM(m)  = Layer1BM(m)  + dataarray(4,j,m) * &
                                (dataarray(10,j,m)+dataarray(11,j,m))/10000
              endif
!             total biomass
              totalBM(m)  = totalBM(m)  + (dataarray(11,j,m)+dataarray(10,j,m))  &
                            * dataarray(4,j,m)/10000
              totCAI(m)   = totCAI(m) + dataarray(9,j,m) * dataarray(4,j,m)/10000.0
!             woody biomass residence time
              DBH = dataarray(7,j,m)
              iLayer=dataarray(3,j,m)
              iPFT  =dataarray(2,j,m)


!             size classes distribution, including total individuals; biomass distribution
              do i=1,bins
                  if(dataarray(7,j,m)>DBHbins(i) .and. dataarray(7,j,m)<=DBHbins(i+1) ) then
                     DBHclasses(i,m) = DBHclasses(i,m) + dataarray(4,j,m)
                     BMclasses(i,m)  = BMclasses(i,m)  + dataarray(4,j,m)*         &
                                       (dataarray(11,j,m)+dataarray(10,j,m)) /10000
                   endif
              enddo

              iLayer=dataarray(3,j,m)
              iPFT  =dataarray(2,j,m)  ! dataarray(2,j,m)-2
!             Basal area with DBH > 0.1 m
              if(dataarray(7,j,m)> 0.0)then
                  basalA(m)   = basalA(m)   + 0.25*PI*dataarray(7,j,m)**2 * dataarray(4,j,m)
                  BApft(iPFT,m) = BApft(iPFT,m) + 0.25*PI*dataarray(7,j,m)**2 * dataarray(4,j,m)
                  CApft(iPFT,m) = CApft(iPFT,m) + dataarray(9,j,m) * dataarray(4,j,m)
                  ! for growth
                  NPPlf(iPFT,m) = NPPlf(iPFT,m) + dataarray(13,j,m) * dataarray(15,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  NPPwd(iPFT,m) = NPPwd(iPFT,m) + dataarray(13,j,m) * dataarray(17,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  NPPfr(iPFT,m) = NPPfr(iPFT,m) + dataarray(13,j,m) * dataarray(16,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  NPPsd(iPFT,m) = NPPsd(iPFT,m) + dataarray(13,j,m) * dataarray(14,j,m) * dataarray(4,j,m)/10.0 ! gC per m2


                  siteNPPlf(m,iTests) = siteNPPlf(m,iTests) + dataarray(13,j,m) * dataarray(15,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  siteNPPwd(m,iTests) = siteNPPwd(m,iTests) + dataarray(13,j,m) * dataarray(17,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  siteNPPfr(m,iTests) = siteNPPfr(m,iTests) + dataarray(13,j,m) * dataarray(16,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  siteNPPsd(m,iTests) = siteNPPsd(m,iTests) + dataarray(13,j,m) * dataarray(14,j,m) * dataarray(4,j,m)/10.0 ! gC per m2

                  siteGPP(m,iTests)   = siteGPP(m,iTests)   + dataarray(18,j,m) * dataarray(4,j,m)/10.0 ! gC per m2
                  siteNPP(m,iTests)   = siteNPP(m,iTests)   + dataarray(19,j,m) * dataarray(4,j,m)/10.0 ! gC per m2

                  density(iPFT,m) = density(iPFT,m) + dataarray(4,j,m)
              endif
              ! mean dDBH of the first layer
              if(dataarray(3,j,m)<2.0)then
                  dDBH(iPFT,m) = dDBH(iPFT,m) + dataarray(6,j,m) * dataarray(4,j,m)
                  densityL1(iPFT,m) = densityL1(iPFT,m) + dataarray(4,j,m)
              endif
              ! Record the potential LAI for iPFT
              if(maxLAI(iPFT,m)<0.001)maxLAI(iPFT,m) = dataarray(22,j,m)
              !if(dDBH(iPFT,m) == 0.0 .and. dataarray(3,j,m)<2.0)then
              !    dDBH(iPFT,m) = dataarray(6,j,m) * dataarray(4,j,m)
              !    density(iPFT,m) = density(iPFT,m) + dataarray(4,j,m)
              !endif
          enddo ! j, cohorts

          ! keep BA
          do iPFT=1, i
             compBA(m,iPFT,iTests)=BApft(iPFT,m)
          enddo
          !siteNPPlf(m,iTests)=siteNPPlf(m,iTests)/siteNPP(m,iTests)
          !siteNPPwd(m,iTests)=siteNPPwd(m,iTests)/siteNPP(m,iTests)
          !siteNPPfr(m,iTests)=siteNPPfr(m,iTests)/siteNPP(m,iTests)
          !siteNPPsd(m,iTests)=siteNPPsd(m,iTests)/siteNPP(m,iTests)
          siteBA(m,iTests)   = basalA(m)
          siteCAI(m,iTests)  = totCAI(m)
          siteHTs(m,iTests)  = HTstar(m)


          do i=1,maxPFTs
             if(densityL1(i,m)>0.0) then
                dDBH(i,m)=dDBH(i,m)/densityL1(i,m)
             else
                dDBH(i,m)= 0.0
             endif
          enddo
!         read the last cc and discard
          if(yearch(2,m) > cc)  &
              read(211,*,IOSTAT=istat3)(dataarray(i,cc+1,m),i=1,items)
          if(istat3<0)exit


       enddo ! End of calculating
       totyears=m-1
     else
           write(*,*)" open failed !"
           stop
     endif

 ! Output
! read in another file for combining data
    read(212,*,IOSTAT=istat3)comments
    Data head &
    /   'BA1','BA2','BA3','BA4', 'BA5','BA6','BA7','BA8','BA9','BA10',                               &    ! 1
        'NPP1','NPP2','NPP3','NPP4', 'NPP5','NPP6','NPP7','NPP8','NPP9','NPP10',                     &    ! 2
        'NPPlf1','NPPlf2','NPPlf3','NPPlf4', 'NPPlf5','NPPlf6','NPPlf7','NPPlf8','NPPlf9','NPPlf10', &    ! 3
        'NPPfr1','NPPfr2','NPPfr3','NPPfr4', 'NPPfr5','NPPfr6','NPPfr7','NPPfr8','NPPfr9','NPPfr10', &    ! 4
        'NPPwd1','NPPwd2','NPPwd3','NPPwd4', 'NPPwd5','NPPwd6','NPPwd7','NPPwd8','NPPwd9','NPPwd10', &    ! 5
        'NPPsd1','NPPsd2','NPPsd3','NPPsd4', 'NPPsd5','NPPsd6','NPPsd7','NPPsd8','NPPsd9','NPPsd10', &    ! 6
        'dDBH1','dDBH2','dDBH3','dDBH4', 'dDBH5','dDBH6','dDBH7','dDBH8','dDBH9','dDBH10',           &    ! 7
        'N_PFT1','N_PFT2','N_PFT3','N_PFT4', 'N_PFT5','N_PFT6','N_PFT7','N_PFT8','N_PFT9','N_PFT10', &    ! 8
        'CA1','CA2','CA3','CA4', 'CA5','CA6','CA7','CA8','CA9','CA10',           &    ! 9
        'DBH0.1','DBH0.2','DBH0.3','DBH0.4', 'DBH0.5','DBH0.6','DBH0.7','DBH0.8','DBH0.9','DBH1.0',  &    ! 10
        'CAI-1','CAI-2','CAI-3','CAI-4', 'CAI-5','CAI-6','CAI-7','CAI-8','CAI-9','CAI-10',           &    ! 11
        'HTcri1','HTcri2','HTcri3','HTcri4', 'HTcri5','HTcri6','HTcri7','HTcri8','HTcri9','HTcri10', &    ! 12
        'PlantC1','PlantC2','PlantC3','PlantC4', 'PlantC5','PlantC6','PlantC7','PlantC8','PlantC9','PlantC10', & ! 13
        'SoilC1','SoilC2','SoilC3','SoilC4', 'SoilC5','SoilC6','SoilC7','SoilC8','SoilC9','SoilC10',           & ! 14
        'PlantN1','PlantN2','PlantN3','PlantN4', 'PlantN5','PlantN6','PlantN7','PlantN8','PlantN9','PlantN10', & ! 15
        'SoilN1','SoilN2','SoilN3','SoilN4', 'SoilN5','SoilN6','SoilN7','SoilN8','SoilN9','SoilN10',           & ! 16
        'Nmnr1','Nmnr2','Nmnr3','Nmnr4', 'Nmnr5','Nmnr6','Nmnr7','Nmnr8','Nmnr9','Nmnr10',                     & ! 17
        'GPP1','GPP2','GPP3','GPP4', 'GPP5','GPP6','GPP7','GPP8','GPP9','GPP10',                               & ! 18
        'NSC1','NSC2','NSC3','NSC4', 'NSC5','NSC6','NSC7','NSC8','NSC9','NSC10',                     & ! 19
        'NSN1','NSN2','NSN3','NSN4', 'NSN5','NSN6','NSN7','NSN8','NSN9','NSN10'                     & ! 20
    /
    write(14,'(360(a12,","))') 'year', &
        (head(i,1),i=iTests,iTests+initialPFTs-1), & ! 'BA1',...,'BA10', &
        (head(i,3),i=iTests,iTests+initialPFTs-1), & ! 'NPPlf1',...,'NPPlf10', &
        (head(i,4),i=iTests,iTests+initialPFTs-1), & ! 'NPPfr1',...,'NPPfr10', &
        (head(i,5),i=iTests,iTests+initialPFTs-1), & ! 'NPPwd1','NPPwd2','NPPwd3','NPPwd4', 'NPPwd5','NPPwd6','NPPwd7','NPPwd8','NPPwd9','NPPwd10', &
        (head(i,6),i=iTests,iTests+initialPFTs-1), & ! 'NPPsd1','NPPsd2','NPPsd3','NPPsd4', 'NPPsd5','NPPsd6','NPPsd7','NPPsd8','NPPsd9','NPPsd10', &
        (head(i,7),i=iTests,iTests+initialPFTs-1), & ! 'dDBH1','dDBH2','dDBH3','dDBH4', 'dDBH5','dDBH6','dDBH7','dDBH8','dDBH9','dDBH10', &
        (head(i,8),i=iTests,iTests+initialPFTs-1), & ! 'N_PFT1','N_PFT2','N_PFT3','N_PFT4', 'N_PFT5','N_PFT6','N_PFT7','N_PFT8','N_PFT9','N_PFT10', &
        (head(i,9),i=iTests,iTests+initialPFTs-1), & ! 'CA1','CA2','CA3','CA4', 'CA5','CA6','CA7','CA8','CA9','CA10', &
        (head(i,10),i=1,bins), & ! 'DBH0.1','DBH0.2','DBH0.3','DBH0.4', 'DBH0.5','DBH0.6','DBH0.7','DBH0.8','DBH0.9','DBH1.0', &
        'siteNPP','siteNPPLF','siteNPPWD','siteNPPFR','siteNPPsd', &
        'totBA','critH',  &
        'CAI', 'LAI','GPP','Rauto','Rh', &
        'rain','SiolWater','Transp','Evap','Runoff', &
        'plantC','soilC','plantN','soilN','totN','NSC', &
        'SeedC','leafC','rootC','SapwoodC','WoodC','NSN','SeedN', &
        'leafN','rootN','SapwoodN','WoodN','McrbC','fastSOM','SlowSOM', &
        'McrbN','fastSoilN','slowSoilN','mineralN','N_fxed','N_uptk','N_yrMin'

 do m = 1, totyears
     read(212,*,IOSTAT=istat3)n,(ecodata(i,m),i=1,37)
     GPP      = ecodata(3,m)
     NPP      = ecodata(4,m)
     plantC   = ecodata(11,m)
     soilC    = ecodata(12,m)
     plantN   = ecodata(13,m)
     soilN    = ecodata(14,m)
     mineralN = ecodata(28,m)
     Nmin     = ecodata(37,m)

     !siteGPP(m,iTests)      = ecodata(3,m)*1000
     !siteNPP(m,iTests)      = (ecodata(3,m) - ecodata(4,m))*1000
     siteplantC(m,iTests)   = ecodata(11,m)
     sitesoilC(m,iTests)    = ecodata(12,m)
     siteplantN(m,iTests)   = ecodata(13,m)
     siteSoilN(m,iTests)    = ecodata(14,m)
     siteNSC(m,iTests)      = ecodata(16,m)
     siteNSN(m,iTests)      = ecodata(22,m)
     siteMineralN(m,iTests) = ecodata(37,m) !ecodata(28,m)

     !if(m>500) &
     write(14,140)m,(BApft(i,m),i=iTests,iTests+initialPFTs-1), &
                (NPPlf(i,m),i=iTests,iTests+initialPFTs-1),  &
                (NPPfr(i,m),i=iTests,iTests+initialPFTs-1),  &
                (NPPwd(i,m),i=iTests,iTests+initialPFTs-1),  &
                (NPPsd(i,m),i=iTests,iTests+initialPFTs-1),  &
                (dDBH(i,m),i=iTests,iTests+initialPFTs-1),  &
                (density(i,m),i=iTests,iTests+initialPFTs-1),  &
                !(maxLAI(i,m),i=iTests,iTests+initialPFTs-1),  &
                (CApft(i,m),i=iTests,iTests+initialPFTs-1),  &
                (DBHclasses(i,m),i=1,bins),                &
                siteNPP(m,iTests),siteNPPlf(m,iTests),siteNPPwd(m,iTests),siteNPPfr(m,iTests),siteNPPsd(m,iTests), &
                basalA(m), HTstar(m),  &
                (ecodata(i,m),i=1,37) ! &

 enddo

   close (12)
   close (13)
   close (14)
   close (212)

   enddo ! iTests
  !! Output NPP for all Tests at a given soil
  fileout3 = trim(filepath)//'Tests'//trim(runtag)//trim(chaSOM(iSOM))//trim(numPFTs)//'.csv'
  open(155,file=fileout3)

  write(155,'(960(a12,","))') 'year',             &
        (head(iTests,1),iTests=1,siteplus-initialPFTs), & ! 'BA1',...,'BA10', &
        (head(iTests,18),iTests=1,siteplus-initialPFTs), & ! 'GPP1',..., 'GPP10', &
        (head(iTests,2),iTests=1,siteplus-initialPFTs), & ! 'NPP1',..., 'NPP10', &
        (head(iTests,3),iTests=1,siteplus-initialPFTs), & ! 'NPPlf1',...,'NPPlf10', &
        (head(iTests,4),iTests=1,siteplus-initialPFTs), & ! 'NPPfr1',..., 'NPPfr10', &
        (head(iTests,5),iTests=1,siteplus-initialPFTs), & ! 'NPPwd1',...,'NPPwd10', &
        (head(iTests,6),iTests=1,siteplus-initialPFTs), & ! 'NPPsd1',...,'NPPsd10', &
        (head(iTests,19),iTests=1,siteplus-initialPFTs), & ! 'NSC1',...,'NSC10', &
        (head(iTests,20),iTests=1,siteplus-initialPFTs), & ! 'NSN', &

        (head(iTests,11),iTests=1,siteplus-initialPFTs), & ! 'CAI-1',...,'CAI-10', &
        (head(iTests,12),iTests=1,siteplus-initialPFTs), & ! 'HTcri1',...,'HTcri10', &
        (head(iTests,13),iTests=1,siteplus-initialPFTs), & ! 'PlantC1',...,'PlantC10', &
        (head(iTests,14),iTests=1,siteplus-initialPFTs), & ! 'SoilC1',...,'SoilC10', &
        (head(iTests,15),iTests=1,siteplus-initialPFTs), & ! 'PlantN1','...,'PlantN10', &
        (head(iTests,16),iTests=1,siteplus-initialPFTs), & ! 'SoilN1',...,'SoilN10', &
        (head(iTests,17),iTests=1,siteplus-initialPFTs) ! , & ! 'mineralN1', ..., 'mineralN10', &
        ! 'BA1','BA2','BA2','BA3','BA3','BA4','BA4','BA5','BA5','BA6','BA6','BA7','BA7','BA8','BA8','BA9'
        ! 'BA1','BA2','BA3','BA2','BA3','BA4', 'BA3','BA4','BA5','BA4','BA5','BA6', &
        ! 'BA5','BA6','BA7','BA6','BA7','BA8','BA7','BA8','BA9'


  do m = 1, totyears
      write(155,140)m,(siteBA(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteGPP(m,iTests),   iTests=1,siteplus-initialPFTs), &
                    (siteNPP(m,iTests),   iTests=1,siteplus-initialPFTs), &
                    (siteNPPlf(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteNPPfr(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteNPPwd(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteNPPsd(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteNSC(m,iTests), iTests=1,siteplus-initialPFTs),   &
                    (siteNSN(m,iTests), iTests=1,siteplus-initialPFTs),   &
                    (siteCAI(m,iTests), iTests=1,siteplus-initialPFTs),   &
                    (siteHTs(m,iTests), iTests =1,siteplus-initialPFTs),  &
                    (sitePlantC(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteSoilC(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (sitePlantN(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteSoilN(m,iTests), iTests=1,siteplus-initialPFTs), &
                    (siteMineralN(m,iTests), iTests=1,siteplus-initialPFTs) !,&
                   ! (compBA(m,iTests,iTests),compBA(m,iTests+1,iTests),iTests=1,totTests)
                   ! (compBA(m,iTests,iTests),compBA(m,iTests+1,iTests),compBA(m,iTests+2,iTests),iTests=1,totTests)
 enddo !
 close(155)

enddo ! iSOM



140 format(1(I8,','),10(f15.4,','),65(f15.4,','), 120(f15.4,','))
150 format(13(I8,',',10(f15.4,',')),1(I8,','),120(f15.4,','))
103 format(1(I8,','),11(f15.4,','),42(f15.4,','), 6(I8,',',13(f15.4,',')), 18(2(I8,','),13(f15.4,',')) )
102 format(1(I8,','),4(f15.4,','),6(f15.4,','), 2(I8,',',10(f15.4,',')), 12(2(I8,','),10(f15.4,',')) )
101 format(1(I8,','),4(f15.4,','),6(f15.4,','), 1(I8,',',10(f15.4,',')),  6(2(I8,','),10(f15.4,',')) )

end program
