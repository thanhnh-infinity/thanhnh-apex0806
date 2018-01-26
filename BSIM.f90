      SUBROUTINE BSIM                                                                
!     APEX0806                                                                       
!     THIS SUBPROGRAM DRIVES THE DAILY SIMULATION FOR THE ENTIRE WATER-              
!     SHED FOR USER SPECIFIED NUMBER OF YEARS. IT CALLS SUBAREA, ROUTING             
!     , ADDING, AND RESERVOIR COMPONENTS IN PROPER SEQUENCE. 
      USE PARM                                                                       
      CHARACTER(4)::CMDN 
      DIMENSION CMDN(5)                       
      DIMENSION KTP(5),MNST(5)                                                       
      DIMENSION VARW(NSM),ZTX(33),TDFP(MPS),TDRP(MPS),TDSP(MPS),&
      TLHP(MPS),TRFP(MPS),TSSP(MPS),XTX(14),YTP(13),ZTZ(6)
      DIMENSION XTMP(MPS,MHY)                   
      DATA MNST/1,2,3,4,5/,IIP/1/
      DATA ALMN/.02083/,YTP/13*0./
      CMDN=(/'HYDV','ROUT','ADD ','RSRT','PDRT'/)
      RWSX=RWSA(NCMD)
      IF(NPSO>0)THEN
          DO I=1,NPSO
              II=I+LND
              DO J=1,6
                  READ(KR(II),'(A80)')ADUM
              END DO
          END DO
      END IF
      IPSF=0.
      T2=0.
      DO 940 IY=1,NBYR                                                               
      JDA=IBD
      JJ=IBD-1                                                                       
      IF(ICO2==1)THEN                                                                 
          IF(IYX<25)THEN                                                             
              CO2=280.33                                                
          ELSE                                                                       
              X1=IYX                                                                 
              CO2=280.33-X1*(.1879-X1*.0077)                                         
          END IF                                                                     
      END IF                                                                         
      IF(IPD>2)THEN                                                                  
          IF(NOP/=0)THEN                                                             
              CALL APAGE(0)                                                          
              IF(KFL(1)>0)WRITE(KW(1),'(//T11,A,A)')'SA(# ID)   Y M D',&             
              ' OPERATION'                                                           
          END IF                                                                     
      ELSE                                                                           
          CALL APAGE(0)                                                              
      END IF                                                                         
      ND=366-NYD                                                                     
      DO I=1,LC                                                                      
     	    WA(I)=100.*CO2/(CO2+EXP(WAC2(1,I)-CO2*WAC2(2,I)))                          
      END DO                                                                         
      WRITE(KW(1),'(T10,A,F5.0)')'ATMOS CO2 = ',CO2                                  
      IPC=1                                                                          
      KOMP=0                                                                         
      DO ISA=1,MSA                                                                   
          KTMX(ISA)=1                                                                
          TFLG(ISA)=0.                                                               
          AFLG(ISA)=0.                                                               
          IRO(ISA)=IRO(ISA)+1                                                        
          IF(IRO(ISA)>NRO(ISA))IRO(ISA)=1                                            
          DO I=1,LC                                                                  
              JE(I,ISA)=12                                                           
              IF(KGO(I,ISA)==0)CYCLE                                                 
              DO J=1,NCP(IRO(ISA),ISA)                                               
                  IF(I==LY(IRO(ISA),J,ISA))EXIT                                      
              END DO                                                                 
              JE(J,ISA)=I                                                            
          END DO
          IF(IPAT>0)THEN
		      IF(PDPLC(ISA)<20.)THEN
			      APMU=2.25*(30.-PDPLC(ISA))
			      JJK=LY(IRO(ISA),1,ISA)
	    	      IF(APMU>45.)CALL NFERT(APMU,7,IAUF(ISA),JRT)
		      END IF
	      END IF                                                                     
          LD1=LID(1,ISA)                                                             
          ANA(IRO(ISA),ISA)=0.                                                       
          HSM(ISA)=0.                                                                
          X1=.2+.3*EXP(-.0256*SAN(LD1,ISA)*(1.-.01*SIL(LD1,ISA)))                    
          X2=(SIL(LD1,ISA)/(CLA(LD1,ISA)+SIL(LD1,ISA)))**.3                          
          X5=.1*WOC(LD1,ISA)/WT(LD1,ISA)                                             
          X3=1.-.25*X5/(X5+EXP(3.718-2.947*X5))                                      
          X4=1.-.01*SAN(LD1,ISA)                                                     
          X5=1.-.7*X4/(X4+EXP(-5.509+22.899*X4))                                     
          EK(ISA)=X1*X2*X3*X5                                                        
          USL(ISA)=EK(ISA)*SLF(ISA)*PEC(ISA)                                         
	      SUM=(SAN(LD1,ISA)*.0247-SIL(LD1,ISA)*3.65-CLA(LD1,ISA)*6.908)/&                 
          100.                                                                       
          DG=EXP(SUM)                                                                
          REK=9.811*(.0034+.0405*EXP(-.5*((LOG10(DG)+1.659)/.7101)**2))              
          RSK(ISA)=REK*PEC(ISA)*RSF(ISA)                                             
          RSLK(ISA)=RSK(ISA)*RLF(ISA)                                                
          CALL EWIK
          SAMA(ISA)=0.                                                               
          XMAP(ISA)=0.
          SELECT CASE(MNUL)                                                               
	          CASE(1)                                                                          
                  IF(PDPLC(ISA)<62.)THEN                                                 
                      XMAP(ISA)=2.*UPR                                                   
                      CYCLE                                                              
                  END IF                                                                 
                  IF(PDPLC(ISA)<120.)THEN                                                
                      XMAP(ISA)=1.5*UPR                                                  
                      CYCLE                                                              
                  END IF                                                                 
                  IF(PDPLC(ISA)<200.)THEN                                                
                      XMAP(ISA)=UPR                                                      
                      CYCLE                                                              
                  END IF
              CASE(2)                                                                
	              IF(PDPLC(ISA)<200.)CYCLE                                               
  	              XMAP(ISA)=UNR                                                             
	              CYCLE                                                                  
              CASE(3)
                  IF(IAPL(ISA)<=0)CYCLE
                  IF(PDPLC(ISA)<200.)XMAP(ISA)=2.*UPR
	          CASE DEFAULT     
	      END SELECT                                                                
      END DO                                                                              
      NT1=0                                                                          
      WRITE(*,50)IY,NBYR                                                             
!     BEGIN DAILY SIMULATION.                                                        
      DO IDA=IBD,ND
          NTX=0                                                                          
          ISA=1                                                                          
          NOFL=40 
          NHY=0                                                                       
          IOF=1
    !     SLAP=0.                                                                          
          CALL AICL                                                                      
          XDA=KDA                                                                        
          IPC=MAX(IPC,IDA)                                                               
          XDA1=31.-XDA                                                                   
          TSFN=0.;QRP=0.;SST=0.;YP=0.;YN=0.;YC=0.;QC=0.;YCOU=0.;YPOU=0.                       
          YNOU=0.;QN=0.;QP=0.;QPU=0.;QURB=0.;YMNU=0.;RQRB=0.;QVOL=0.                             
	      YSD=0.;QPST=0.;YPST=0.                                                              
  	      RFV=0.;CVF=0.                                                                     
  	      IF(NPSO>0)THEN                                                                      
	          DO I=1,NPSO
	              IF(IPSF(I)>0)THEN
				      PSOQ(I)=0.
				      PSOY(I)=0.
				      PSON(I)=0.
				      PSOP(I)=0.
				      PSO3(I)=0.
				      PSSP(I)=0.
				      PQPS(I)=0.
				      PYPS(I)=0.
				      CYCLE
			      END IF                                                                     
	              II=I+LND                                                                    
                  !     POINT SOURCE INPUT
                  !             I1   = DAY OF YEAR
                  !             I2   = YEAR
                  !             PSOQ = FLOW(M3/D)
                  !             PSOY = SEDIMENT LOAD(T/D)
                  !             PSON = ORGANIC N(KG/D)
                  !             PSOP = ORGANIC P(KG/D) 
                  !             PSO3 = NITRATE N LOAD(KG/D)
                  !             PSH3 = AMMONIA N LOAD(KG/D)
                  !             PSO2 = NITRITE N LOAD(KG/D)
                  !             PSSP = SOLUBLE P LOAD(KG/D)
                  !             PBOD = BOD LOAD(KG/D)
                  !             PSDO = DISOLVED OXYGEN LOAD(KG/D)
                  !             PSCA = CHLORAFIL LOAD(KG/D)
                  !             PQPS = SOLUBLE PESTICIDE LOAD(G/D)
                  !             PYPS = ADSORBED PESTICIDE LOAD(G/D)
                  !             KPSN = PESTICIDE # FROM PESTCOM.DAT
			      READ(KR(II),*,IOSTAT=NFL)I1,I2,PSOQ(I),PSOY(I),PSON(I),&
                  PSOP(I),PSO3(I),PSH3,PSO2,PSSP(I),PBOD,PSDO,PSCA,PQPS(I),&
                  PYPS(I),X1,X2,X3,X4,X5,KPSN(I)
			      IF(NFL/=0)THEN
			          WRITE(KW(36),'(5X,A,I8,A,3I4)')'POINT SOURCE #',I,&
	                  ' ENDS (Y M D)',IYR,MO,KDA                      
	                  IPSF(I)=1                           
	              END IF
	              IF(NDP>0)THEN
				      DO L=1,NDP
					      IF(KDP(L)==KPSN(I))EXIT
				      END DO
			          IF(L<=NDP)CYCLE
			      END IF
		          NDP=NDP+1
			      KDP(NDP)=KPSN(I)
			      KDP1(KPSN(I))=NDP
                  ! READ PESTICIDE TABLE
                  !	1  PSTN = PESTICIDE NAME
                  ! 2  PSOL = PESTICIDE SOLUBILITY (ppm)
                  ! 3  PHLS = PESTICIDE HALF LIFE IN SOIL (d)
                  ! 4  PHLF = PESTICIDE HALF LIFE ON FOLIAGE (d)
                  ! 5  PWOF = PESTICIDE WASH OFF FRACTION
                  ! 6  PKOC = PESTICIDE ORGANIC C ADSORPTION COEF
                  ! 7  PCST = PESTICIDE COST ($/KG)
                  J1=-1
			      DO WHILE(J1/=KPSN(I))
				      READ(KR(8),123,IOSTAT=NFL)J1,PSTN(NDP),(XTP(L),L=2,7)
				      IF(NFL/=0)THEN
					      WRITE(*,*)'PEST NO = ',KPSN(I),' NOT IN PEST FILE'
					      PAUSE
					      STOP
				      END IF
                  END DO
			      PSOL(NDP)=XTP(2)
			      PHLS(NDP)=XTP(3)
			      PHLF(NDP)=XTP(4)
			      PWOF(NDP)=XTP(5)
			      PKOC(NDP)=XTP(6)
			      PCST(NDP)=XTP(7)
			      REWIND KR(8)                                     
              END DO                                                                     
          END IF                                                                         
          REP=0.                                                                         
          IF(IHY>0)THEN
              NRF=0                                                                          
              RFDT=0.
              YHY=0.
          END IF
          IF(INFL==5)THEN
              CALL HRFIN
		      RFV0(1)=RFDT(NRF)
		      GO TO 141
          END IF
          IF(NGN<=0)GO TO 140                                                            
          DO I=1,NWTH                                                                    
              II=I+KND                                                                   
              !     READ DAILY WEATHER IF NOT GENERATED                                            
              !  1  SRAD = SOLAR RADIAION(MJ/m2 OR LY)(BLANK TO GENERATE)                          
              !  2  TMX  = MAX TEMP(c)                                                             
              !  3  TMN  = MIN TEMP(c)                                                             
              !  4  RFV0 = RAINFALL(mm)(999. TO GENERATE OCCURRENCE & AMOUNT;                      
              !            -1. TO GENERATE AMOUNT GIVEN OCCURRENCE)                                
              !  5  RHD  = RELATIVE HUMIDITY(FRACTION)(BLANK TO GENERATE)                          
              !  6  U10  = WIND VELOCITY(m/s)(BLANK TO GENERATE)                                   
              !  7  CO2I = ATMOSPHERIC CO2 CONC(ppm)                                                 
              !  8  REP  = PEAK RAINFALL RATE(mm/h)
              !  9  ORSD = OBSERVED SOIL SURFACE CROP RESIDUE(t/ha) 
              READ(KR(II),1070,IOSTAT=NFL)SRAD(I),TMX(I),TMN(I),RFV0(I),RHD&             
              (I),U10(I),CO2I,REP,ORSD(I)                                                             
              IF(NFL/=0)THEN                                                             
                  NGN=-1                                                                 
                  IRF=1
                  KGN=0                                                                  
                  WRITE(KW(1),'(5X,A,A,3I4)')'GENERATED WEATHER STARTS',&               
                  ' (Y M D)',IYR,MO,KDA                                                  
                  GO TO 140                                                              
              END IF
              SRAD(I)=SRAD(I)*RUNT                                                            
              RHD(I)=RHD(I)/RHCF                                                         
              III=0                                                                      
              IF(RFV0(I)>900..OR.RFV0(I)<0.)THEN                                                         
                  CALL WRWD(I,0)                                                           
                  RFV0(I)=RFV0(I)*RNCF(MO)                                               
              END IF                                                                     
              IF(KGN(2)>0.AND.ABS(TMX(I)-TMN(I))>0.)THEN
                  IF(TMX(I)>100..OR.TMN(I)>100.)THEN
                      CALL WGN                                                                   
                      CALL WTAIX(I)                                                              
                      GO TO 111
                  ELSE            
                      IF(TMX(I)>TMN(I))GO TO 110
                  END IF
              END IF
              CALL WGN                                                                   
              CALL WTAIR(I)                                                                 
      111     X1=TMX(I)-TMN(I)                                                           
              TMX(I)=TMX(I)+TMXF(MO)                                                     
              TMN(I)=TMX(I)-TMNF(MO)*X1                                                  
      110     IF(RHD(I)==0..OR.KGN(5)==0)THEN                                            
                  IF(III==0)CALL WGN                                                         
                  CALL WRLHUM(I)                                                                
              END IF                                                                     
              IF(U10(I)>0..AND.KGN(4)>0)CYCLE                                            
              CALL WNSPD                                                                 
          END DO                                                                         
          IWI=1                                                                          
          GO TO 160
      140 U10(1)=0.                                                                      
          RHD(1)=0.
          IF(NWP>1)THEN                                                                  
              XX=AUNIF(IDG(12))                                                          
              DO IWI=1,NWP                                                               
                  IF(XX<=AWXP(IWI))EXIT                                                    
              END DO                                                                     
          ELSE                                                                           
              IWI=1                                                                      
          END IF                                                                         
          CALL WRWD(1,0)                                                                   
          NWPD(IWI)=NWPD(IWI)+1                                                          
      141 RFV0(1)=RFV0(1)*RNCF(MO)                                                       
          CALL WGN                                                                       
          CALL WTAIR(1)                                                                     
          X1=TMX(1)-TMN(1)                                                               
          TMX(1)=TMX(1)+TMXF(MO)                                                         
          TMN(1)=TMX(1)-TMNF(MO)*X1                                                      
          CALL WRLHUM(1)                                                                    
          CALL WNSPD                                                                     
          TMX=TMX(1)                                                                     
          TMN=TMN(1)                                                                     
          SRAD=SRAD(1)                                                                   
          RHD=RHD(1)                                                                     
          U10=U10(1)                                                                     
      160 IXP=NXP(90)
          IF(ICO2==2)THEN
              IF(CO2I>0.)CO2=CO2I                                                                    
          END IF
          IF(RFV0(1)>0.)THEN
              SDRF=SDRF+RFV0(1)*RFV0(1)                                                      
              IF(RFV0(1)>RMX0)RMX0=RFV0(1)                                                   
              IF(ITYP==0)THEN                                                                
                  AJP=1.-EXP(-125./(RFV0(1)+5.))                                             
                  AL5=ATRI(ALMN,WI(IWI,MO),AJP,4)                                            
              ELSE                                                                           
                  AL5=WI(IWI,MO)                                                             
              END IF                                                                         
              PRFF=-2.*LOG(1.-AL5)                                                           
              DUR=MIN(24.,4.605/PRFF)                                                        
              IF(NGN<0)THEN
                  DO ISA=1,MSA                                                                   
                      RFV(ISA)=RFV0(1)                                                           
                  END DO
              END IF
              IF(NGN==0)THEN                             
                  XX=AUNIF(IDG(10))                                                              
                  YY=AUNIF(IDG(11))                                                              
                  XX=XX*XSL+XCS                                                                  
                  YY=YY*YSL+YCS                                                                  
                  SUM=0.                                                                         
                  X1=DUR**(-.1478)                                                               
                  DO ISA=1,MSA                                                                   
                      XA=ABS(XCT(ISA))                                                           
                      YA=ABS(YCT(ISA))                                                           
                      D2=(XX-XA)**2+(YY-YA)**2                                                   
                      D=SQRT(D2)                                                                 
                      X2=MAX(.001,1.-X1*D/(D+EXP(SCRP(21,1)-SCRP(21,2)*D)))                      
                      ZTP(ISA)=X2*(1.+BXCT*XA)*(1.+BYCT*YA)                                      
                      SUM=SUM+ZTP(ISA)                                                           
                  END DO                                                                         
                  RTO=XSA*RFV0(1)/SUM                                                            
                  RFU=RFV0(1)                                                                    
                  RFL=RFV0(1)                                                                    
                  DO ISA=1,MSA                                                                   
                      X1=ATRI(.8,1.,1.2,6)                                                       
                      RFV(ISA)=RTO*ZTP(ISA)*X1                                                   
                      IF(RFV(ISA)>RFU)THEN
                          RFU=RFV(ISA)
                      ELSE   
                          IF(RFV(ISA)<RFL)RFL=RFV(ISA)
                      END IF                                                                             
                  END DO                                                                         
              END IF
          END IF
          DO IOW=1,NBON                                                                  
              DO IHD=1,NHRD(IOW)                                                         
	              II=IHBS(IHD,IOW)                                                            
	              IYMD=10000*IYHO(IHD,IOW)+100*MO+KDA                                         
	              IF(IHDT(II,IHD,IOW)/=IYMD)CYCLE                                             
	              NCOW(IHD,IOW)=NBSX(II,IHD,IOW)                                              
	              IF(KFL(1)>0)WRITE(KW(1),892)IY,MO,KDA,IOW,IHD,NCOW(IHD,&                    
                  IOW)                                                                   
  	              IHBS(IHD,IOW)=IHBS(IHD,IOW)+1                                             
	              IF(IHBS(IHD,IOW)>NHBS(IHD,IOW))IHBS(IHD,IOW)=1                              
	          END DO                                                                          
              IF(NSAS(IOW)>=2.AND.MNUL<3)THEN                                           
                  X5=PDPL(IDSS(1,IOW))                                                   
                  ISAS(IOW)=IDSS(1,IOW)                                                  
                  DO J=2,NSAS(IOW)                                                       
                      IF(X5<PDPL(IDSS(J,IOW)))CYCLE                                      
                      X5=PDPL(IDSS(J,IOW))                                               
                      ISAS(IOW)=IDSS(J,IOW)                                              
                  END DO                                                                 
              END IF                                                                     
              IF(NSAL(IOW)==0)CYCLE                                                      
              X5=PDPL(IDSL(1,IOW))                                                       
              ISAL(IOW)=IDSL(1,IOW)                                                      
              IF(NSAL(IOW)<2)CYCLE                                                       
              DO J=2,NSAL(IOW)                                                           
                  IF(X5<PDPL(IDSL(J,IOW)))CYCLE                                          
                  X5=PDPL(IDSL(J,IOW))                                                   
                  ISAL(IOW)=IDSL(J,IOW)                                                  
              END DO                                                                     
          END DO                                                                         
          SYW=0.                                                                         
          SDEP=0.                                                                        
	      SDEG=0.                                                                             
	      SRYF=0.                                                                             
	      VAR=0.
	      VARH=0.
	      TAC=0.                                                                              
          DO ICMD=1,NCMD                                                             
              SELECT CASE(ICDT(ICMD))                                                        
                  CASE(1)                                                                    
                      CALL BSUB(IXP)                                                         
                      WSAX=WSA(NISA(IDNB(IDO)))                                              
	                  X1=WSAX*(QN(IDO)+RSFN(IDO)+TSFN(IDO)+QRFN(IDO)+QDRN(IDO))
	                  VARH(13,IDO)=X1                                                             
	                  SMMH(13,MO,IDO)=SMMH(13,MO,IDO)+X1                                          
                      X1=WSAX*QP(IDO)                                                        
	                  VARH(19,IDO)=X1                                                             
	                  SMMH(19,MO,IDO)=SMMH(19,MO,IDO)+X1                                          
	                  SST(IDO)=SST(IDO)+QRP(IDO)                                                  
	              CASE(2)
	                  TDEG=0.
                      TDEP=0.                                                                         
                      CALL ROUTE                                                             
                  CASE(3)                                                                    
                      CALL RTADD                                                             
                  CASE(4)                                                                    
                      CALL RESRT                                                             
                  CASE(5)                                                                    
                      CALL RESPOND                                                           
              END SELECT                                                                     
              SRCH(7,IDO)=SRCH(7,IDO)+QVOL(IDO)                                              
              SRCH(8,IDO)=SRCH(8,IDO)+YSD(NDRV,IDO)                                          
              SRCH(9,IDO)=SRCH(9,IDO)+YN(IDO)+YNOU(IDO)                                                
              SRCH(10,IDO)=SRCH(10,IDO)+YP(IDO)                                              
              SRCH(11,IDO)=SRCH(11,IDO)+QN(IDO)                                              
              SRCH(12,IDO)=SRCH(12,IDO)+QP(IDO)                                              
              SRCH(15,IDO)=SRCH(15,IDO)+SST(IDO)                                             
              SRCH(16,IDO)=SRCH(16,IDO)+RSSF(IDO)                                            
              SRCH(17,IDO)=SRCH(17,IDO)+TSFN(IDO)                                            
              SRCH(21,IDO)=SRCH(21,IDO)+QRFN(IDO)                                            
              SRCH(22,IDO)=SRCH(22,IDO)+QDRN(IDO)                                            
              SRCH(18,IDO)=SRCH(18,IDO)+RSFN(IDO)                                            
              SRCH(19,IDO)=SRCH(19,IDO)+QRF(IDO)                                             
              SRCH(20,IDO)=SRCH(20,IDO)+QDR(IDO)                                             
              SRCH(23,IDO)=SRCH(23,IDO)+CPVH(IDO)                                            
              SRCH(24,IDO)=SRCH(24,IDO)+YMNU(IDO)                                            
              SRCH(25,IDO)=SRCH(25,IDO)+.001*YC(IDO)                                         
	          SRCH(26,IDO)=SRCH(26,IDO)+QPU(IDO)
	          VARW(13)=QVOL(IDO)
              VARW(15)=SST(IDO)
              VARW(NDVSS)=YSD(NDRV,IDO)
              VARW(37)=YN(IDO)
              VARW(38)=QN(IDO)
              VARW(39)=TSFN(IDO)
              VARW(72)=RSSF(IDO)
              VARW(80)=RSFN(IDO)
              VARW(17)=QDR(IDO)
              VARW(47)=QDRN(IDO)
              VARW(48)=YP(IDO)
              VARW(49)=QP(IDO)
              VARW(76)=QC(IDO)
              VARW(77)=YC(IDO)
              VARW(83)=QRF(IDO)
              VARW(84)=QRFN(IDO)
              VARW(88)=YMNU(IDO)
              VARW(108)=QPU(IDO)
              VARW(117)=WYLD(IDO)
              IF(IHY>0)THEN
                  WRITE(KW(26),'(1X,3I4,1X,A4)')IY,MO,KDA,CMDN(ICDT(ICMD))
                  T1=0.
                  AD1=0.
                  X2=0.
                  AD1=-.5*(QHY(1,IDO,IHX(1))+QHY(NPD,IDO,IHX(1)))
                  DO K=1,NPD
                      X1=QHY(K,IDO,IHX(1)) 
                      AD1=AD1+X1
                      IF(X1>X2)THEN
                          X2=X1
                          X3=T1
                      END IF
                      T1=T1+DTHY
                  END DO
                  AD1=AD1*DTHY*360./RWSA(IDO)
	              HYDV(IDO)=AD1
                  IF(AD1>=QTH)THEN
                      T1=0.
	                  DO K=1,NPD
	                      WRITE(KW(26),26)T1,(QHY(K,IDO,IHX(J)),J=1,MHX)
	                      T1=T1+DTHY
	                  END DO
                      IF(IDO==NCMD.AND.KFL(12)>0)THEN
	                      T1=0.
	                      DO K=2,NPD
	                          T1=T1+DTHY
	                          T2=T2+DTHY
	                          WRITE(KW(12),'(5X,3I4,F8.3,2F10.2)')IY,MO,KDA,T1,&
	                          T2,QHY(K,IDO,IHX(1))
	                          ADHY=ADHY+QHY(K,IDO,IHX(1))
                          END DO
	                  END IF    
	              END IF
	              SQVL(IDO)=SQVL(IDO)+QVOL(IDO)
                  SHYD(IDO)=SHYD(IDO)+AD1
	              IF(ICDT(ICMD)==3)THEN
	                  II=0
	              ELSE
	                  II=IDNB(IDO)
	              END IF
	              WRITE(KW(26),8)IDO,II,AD1,X2,X3
                  WRITE(KW(25),12)CMDX(ICDT(IDO)),IDO,II,IY,MO,KDA,X2,&
                  X3,AD1,SQVL(IDO),SHYD(IDO)
	          END IF 
	          !IF(ICDT(ICMD)==2.AND.IDO==5)THEN
	              !DF=YSD(NDRV,4)+TDEG-TDEP-YSD(NDRV,IDO)
	              !IF(ABS(DF)>.001)WRITE(KW(1),35)YSD(NDRV,4),TDEG,TDEP,YSD(NDRV,IDO),DF
	              !DF=SRCH(8,4)+SRCH(14,IDO)-SRCH(13,IDO)-SRCH(8,IDO)
	              !IF(ABS(DF)>.1)WRITE(KW(1),35)SRCH(8,4),SRCH(14,IDO),&
	              !SRCH(13,IDO),SRCH(8,IDO),DF
              !END IF	                                                        
        !	SUM=0.                                                                             
        !     DO 908 L=1,NBSL(ISA)                                                           
        !     ISL=LID(L,ISA)                                                                 
        !     YTP(ISL)=QSF(ISL,ISA)+SSF(ISL,ISA)+CPFH(ISL,ISA)                               
        !     SUM=SUM+YTP(ISL)                                                               
        !     ZTP(ISL)=100.*PSSF(1,ISL,ISA)/(YTP(ISL)+.001)                                  
        ! 908 CONTINUE                                                                       
        !     WRITE(KW(1),154)IYR,MO,KDA,IDO,RFV(ISA),QVOL(IDO),TMX(IRF(ISA)),               
        !    1TMN(IRF(ISA)),SRAD(IRF(ISA)),EO,AET                                            
        !     WRITE(KW(1),906)(PO(LID(L,ISA),ISA),L=1,NBSL(ISA))                             
        !     WRITE(KW(1),906)(FC(LID(L,ISA),ISA),L=1,NBSL(ISA))                             
        !     WRITE(KW(1),906)(S15(LID(L,ISA),ISA),L=1,NBSL(ISA))                            
        !     WRITE(KW(1),906)(ST(LID(L,ISA),ISA),L=1,NBSL(ISA))                             
        !     WRITE(KW(1),906)(PKRZ(LID(L,ISA)),L=1,NBSL(ISA))                               
        !     WRITE(KW(1),906)(QSF(LID(L,ISA),ISA),L=1,NBSL(ISA))                            
        !     WRITE(KW(1),906)(SSF(LID(L,ISA),ISA),L=1,NBSL(ISA))                            
        !     WRITE(KW(1),906)(CPFH(LID(L,ISA),ISA),L=1,NBSL(ISA))                           
        !     WRITE(KW(1),906)(YTP(LID(L,ISA)),L=1,NBSL(ISA)),SUM                            
        !     WRITE(KW(1),906)(PSTZ(1,LID(L,ISA),ISA),L=1,NBSL(ISA))                         
        !     WRITE(KW(1),906)(PSSF(1,LID(L,ISA),ISA),L=1,NBSL(ISA))                         
        !     WRITE(KW(1),906)(ZTP(L),L=1,NBSL(ISA))                                         
          END DO
          TAC=TAC+YC(IDO)
          IF(NDP>0)THEN                                                                  
	          X1=SST(IDO)+QRF(IDO)+CPVH(IDO)+QDR(IDO)                                                      
	          IF(KFL(6)>0)THEN
	              DO K=1,NDP
	                  IF(K==1)THEN
	                      WRITE(KW(6),903)MSA,NBSA(MSA),IYR,MO,KDA,RFV(IRF(ISA)),&
	                      QVOL(IDO),X1,YSD(NDRV,IDO),PSTN(K),(VARP(J,K,IDO),J=1,12)             
                      ELSE
                          WRITE(KW(6),904)MSA,NBSA(MSA),IYR,MO,KDA,PSTN(K),&
                          (VARP(J,K,IDO),J=1,12)        
                      END IF
                  END DO
              END IF
  	          IF(KFL(29)>0)WRITE(KW(29),909)MSA,NBSA(MSA),IYR,MO,KDA,RFV(IRF&                     
              (ISA)),QVOL(IDO),X1,YSD(NDRV,IDO),(PSTN(K),QPST(K,IDO),TSPS(K,IDO)&            
              ,YPST(K,IDO),K=1,NDP)
    !         IF(MSA>1)THEN
    !             QQ=QVOL(IDO)+X1
    !	          YY=YSD(NDRV,IDO)
    !	          SSPZ=0.
    !             CALL PSTFRQ(QQ,YY,SSPZ,NCMD,IXP)                                                          
    !         END IF
          END IF                                                                         
          VARP=0.                                                                        
   	      IF(KFL(16)>0)WRITE(KW(16),894)IYR,MO,KDA,RFV0(IRF(1)),(VARW&
          (KD(K)),K=1,NKD)
          IF(IHRD==2)THEN                                                           
              DO I=1,NBON                                                                
                  DO IHD=1,NHRD(I)                                                           
                      I1=IGZX(IHD,I)                                                                 
                      IF(NGZA(IHD,I)<2)CYCLE
                      X1=AGPM(I1)                                                                    
                      IGZ(I1)=IGZ(I1)+1                                                              
                      IF(IGZ(I1)<=LGRZ)THEN
                          IF(X1>GZLM(IHD,I1).AND.NGZ(IHD,I1)>0)CYCLE
                      END IF
                      J1=IGZO(IHD,I)                                                                 
                      DO K1=1,NGZA(IHD,I)                                                            
                          J1=J1+1                                                                        
                          IF(J1>NGZA(IHD,I))J1=1                                                         
                          I2=NGIX(J1,IHD,I)                                                                  
                          IF(IGZ(I2)==0)EXIT                                                          
                      END DO                                                                         
                      IGZX(IHD,I)=I2                                                                 
                      IGZ(I2)=1                                                                      
                      IGZ(I1)=0                                                                      
                      IGZO(IHD,I)=J1                                                                 
                  END DO
              END DO         
          END IF 
          DO I=90,2,-1                                                                   
              I1=I-1                                                                         
              NXP(I)=NXP(I1)                                                                 
          END DO                                                                         
          NXP(1)=IXP                                                                     
          SMMR(11,MO)=SMMR(11,MO)+QVOL(IDO)                                              
          SMMR(12,MO)=SMMR(12,MO)+QVOL(IDO)+RSSF(IDO)+QRF(IDO)+SST(IDO)+&                
          QDR(IDO)                                                                       
          SMMR(13,MO)=SMMR(13,MO)+YSD(NDRV,IDO)                                          
          SMMR(14,MO)=SMMR(14,MO)+YN(IDO)                                                
          SMMR(15,MO)=SMMR(15,MO)+YP(IDO)                                                
          SMMR(16,MO)=SMMR(16,MO)+QN(IDO)+RSFN(IDO)+TSFN(IDO)+QRFN(IDO)+&                
          QDRN(IDO)                                                                      
          SMMR(17,MO)=SMMR(17,MO)+QP(IDO)                                                
          SMMR(18,MO)=SMMR(18,MO)+YMNU(IDO)                                              
          SMMR(19,MO)=SMMR(19,MO)+QPU(IDO)                                                    
	      SMMR(20,MO)=SMMR(20,MO)+YC(IDO)                                                     
	      SMMR(21,MO)=SMMR(21,MO)+QN(IDO)                                                     
          DO K=1,NDP                                                                     
              SMMRP(3,K,MO)=SMMRP(3,K,MO)+QPST(K,IDO)+TSPS(K,IDO)+RSPS(K,IDO)                                        
              SMMRP(4,K,MO)=SMMRP(4,K,MO)+YPST(K,IDO)                                        
          END DO                                                                         
          !CALL NCONT(3)                                                               
          IF(KFL(35)>0.AND.IPD>=6.AND.MSA>1)THEN
              DO I=1,MSA                                                                     
                  I1=NBSA(IBSA(I))                                                               
	              I2=NISA(I1)
	              IF(IEXT(I2)>0)THEN
			          II=IDOA(I2)
			          X2=WSA(I2)
		          ELSE
			          II=IDOA(I2)-1
			          X2=RWSA(II)+WSA(I2)
		          END IF		
		          IF(VARH(7,II)<1.)VARH(7,II)=VARH(6,II)/(VARH(2,II)*.0864+1.E-10)                    
	              IF(NTX(II)==0)WRITE(KW(35),472)I1,IYR,IDA,X2,VARH(1,II),&
	              VARH(2,II),VARH(33,II),VARH(34,II),(VARH(K,II),K=3,NSH-2),&
	              (PCTH(J,II),J=1,NSZ),(PCT(J,II),J=1,NSZ)                             
     	          NTX(II)=1                                                                      
	          END DO
          END IF	                                                                                    
          IF(KFL(34)>0.AND.IPD>5)THEN                                                    
	          DO I=1,MSA                                                                          
	              I1=NBSA(IBSA(I))                                                                    
	              I2=NISA(I1)                                                                         
	              II=IDOA(I2)
	              X2=.001*ZOC(I2)
	              WRITE(KW(34),471)I1,IYR,IDA,WSA(I2),(VAR(KY(K),I2),K=1,NKY),&
	              X2,(PCTH(J,II),J=1,NSZ),(PCT(J,II),J=1,NSZ),(YLD1(LY(IRO(I2),J,&
	              I2),I2),YLD2(LY(IRO(I2),J,I2),I2),(VARC(K,LY(IRO(I2),J,I2),&
	              I2),K=1,15),CPNM(LY(IRO(I2),J,I2)),J=1,NCP(IRO(I2),I2))                                    
  	          END DO                                                                            
          END IF 
          IF(IHY>0)THEN
              DO I=1,NCMD
                  DO K=1,NPD
                      QHY(K,I,IHX(1))=0. 
                  END DO
              END DO
              II=IHX(1)
              DO I=2,MHX
                  IHX(I-1)=IHX(I)
              END DO
              IHX(MHX)=II
          END IF
          IF(KFL(5)>0)THEN                                                               
              X2=10.*RWSX
              X3=MIN(.95,DRSW*PRMT(63))
              X4=MIN(.95,DRSW*PRMT(75))                                                                      
              ZTZ(1)=QVOL(IDO)+RSSF(IDO)+QRF(IDO)+SST(IDO)+QDR(IDO)                     
	          SMSQ=SMSQ+ZTZ(1)*X2                                                                    
              ZTZ(2)=YSD(NDRV,IDO)*DRSW                                                        
              ZTZ(3)=YN(IDO)*X3                                                              
              ZTZ(4)=YP(IDO)*X4                                                              
              ZTZ(5)=QN(IDO)+RSFN(IDO)+QRFN(IDO)+TSFN(IDO)+QDRN(IDO)                    
              ZTZ(6)=QP(IDO)+QPU(IDO)
              IF(NDP>0)THEN
                  DO K=1,NDP
                      XTMP(K,IDO)=QPST(K,IDO)+TSPS(K,IDO)+RSPS(K,IDO)
                  END DO
              END IF
              WRITE(KW(5),902)IDA,IYR,(ZTZ(I),I=1,6),(XTMP(K,IDO),YPST(K,IDO),&
              K=1,NDP)                                         
          END IF
          IF(IDA==NSTP)STOP
          DO                                                                         
              CALL AXMON                                                                     
              IF(MO==MO1)EXIT
              DO K=1,21                                                                      
                  SMYR(K)=SMYR(K)+SMMR(K,MO1)                                                    
              END DO                                                                         
              DO K=1,NDP                                                                     
                  DO J=1,4                                                                       
                      SMYRP(J,K)=SMYRP(J,K)+SMMRP(J,K,MO1)                                           
                  END DO                                                                         
              END DO                                                                         
              AVRF=0.                                                                        
              ZTZ=0.
              TWMP=0.
              TSMU=0.
              AVRF=0.
              XTX=0.
	          TPRK=0.
	          TSSP=0.
	          TLHP=0.
	          TDFP=0.
	          TDSP=0.
	          TDRP=0.
	          TRFP=0.
	          ZTX=0.                                                                         
              DO ISA=1,MSA                                                               
                  LD1=LID(1,ISA)
                  WSAX=WSA(ISA)                                                                 
                  PMOEO=SMM(10,MO1,ISA)                                                          
                  PMORF=SMM(4,MO1,ISA)-SMM(13,MO1,ISA)                                           
                  XX=IDA-JJ                                                                      
                  STV(6,MO1,ISA)=RZSW(ISA)                                                       
                  STV(7,MO1,ISA)=WTBL(ISA)                                                       
                  STV(2,MO1,ISA)=ZNMN(ISA)                                                       
                  STV(1,MO1,ISA)=ZNMA(ISA)                                                       
                  STV(3,MO1,ISA)=PDPLC(ISA)                                                      
                  STV(8,MO1,ISA)=GWST(ISA)                                                       
                  STV(14,MO1,ISA)=SWLT(ISA)                                                      
                  STV(15,MO1,ISA)=SNO(ISA)                                                       
                  STV(17,MO1,ISA)=GWSN(ISA)                                                      
                  SMM(120,MO1,ISA)=RZSW(ISA)                                                     
                  NN=NCP(IRO(ISA),ISA)                                                           
                  N2=0                                                                           
                  DO K=1,NN                                                                      
                      K1=LY(IRO(ISA),K,ISA)                                                          
                      NTP(K1)=0                                                                      
                  END DO                                                                         
                  DO K=1,NN                                                                  
                      K1=LY(IRO(ISA),K,ISA)                                                          
                      IF(NTP(K1)>0)THEN
                          N2=N2+1                                                                        
                          CYCLE
                      END IF                                                                      
                      NTP(K1)=1                                                                      
                      IF(KGO(K1,ISA)>0)THEN                                                          
                          SMMC(1,K1,MO1,ISA)=HUI(K1,ISA)                                                 
                          SMMC(2,K1,MO1,ISA)=SLAI(K1,ISA)                                                
                          SMMC(3,K1,MO1,ISA)=RD(K1,ISA)                                                  
                          SMMC(4,K1,MO1,ISA)=RW(K1,ISA)                                                  
                          SMMC(5,K1,MO1,ISA)=DM(K1,ISA)                                                  
                          SMMC(6,K1,MO1,ISA)=STL(K1,ISA)                                                 
                          SMMC(7,K1,MO1,ISA)=CPHT(K1,ISA)                                                
                      END IF                                                                         
                      ISM=0                                                                          
                      SMMC(8,K1,MO1,ISA)=STD(K1,ISA)                                                 
                      SMMC(9,K1,MO1,ISA)=STDL(K1,ISA)                                                
                      DO J=1,5                                                                       
                          KTP(J)=SFMO(J,K1,ISA)+.5                                                       
                          ISM=ISM+KTP(J)                                                                 
                          SFMO(J,K1,ISA)=0.                                                              
                      END DO                                                                         
                      IF(ISM==0)THEN
                          KDT(MO1,K1,ISA)=0                                                              
                          CYCLE                                                                      
                      ELSE
                          CALL ASORTI(KTP,MNST,5)                                                        
                          KDT(MO1,K1,ISA)=KTP(MNST(3))+100*MNST(3)+1000*KTP(MNST(4))+&                   
                          100000*MNST(4)+1000000*KTP(MNST(5))+100000000*MNST(5)
                      END IF                          
                  END DO         
                  NN=NN-N2                                                                       
                  STV(10,MO1,ISA)=RSD(LD1,ISA)                                                   
                  STV(16,MO1,ISA)=RSDM(LD1,ISA)                                                  
                  STV(9,MO1,ISA)=STDO(ISA)                                                       
                  SMM(2,MO1,ISA)=SMM(2,MO1,ISA)/XX                                               
                  SMM(1,MO1,ISA)=SMM(1,MO1,ISA)/XX                                               
                  SMM(3,MO1,ISA)=SMM(3,MO1,ISA)/XX                                               
                  SMM(9,MO1,ISA)=SMM(9,MO1,ISA)/XX                                               
                  SMM(59,MO1,ISA)=SMM(59,MO1,ISA)/XX                                             
                  SMM(7,MO1,ISA)=SMM(7,MO1,ISA)/XX                                               
                  SMM(8,MO1,ISA)=SMM(8,MO1,ISA)/XX                                               
                  SMM(33,MO1,ISA)=SMM(33,MO1,ISA)/XX                                             
                  SMM(34,MO1,ISA)=SMM(34,MO1,ISA)/XX                                             
                  SMM(35,MO1,ISA)=SMM(35,MO1,ISA)/XX                                             
                  SMM(60,MO1,ISA)=SMM(60,MO1,ISA)/XX                                             
                  SMM(32,MO1,ISA)=SMM(32,MO1,ISA)/XX                                             
                  SSW(ISA)=SSW(ISA)/XX                                                           
                  ASW(MO1,ISA)=ASW(MO1,ISA)+SSW(ISA)                                             
                  SSW(ISA)=0.                                                                    
                  TR(MO1,ISA)=TR(MO1,ISA)+SMM(4,MO1,ISA)                                         
                  TSN(MO1,ISA)=TSN(MO1,ISA)+SMM(16,MO1,ISA)                                      
                  TSY(MO1,ISA)=TSY(MO1,ISA)+SMM(NDVSS,MO1,ISA)                                   
                  TYON(MO1,ISA)=TYON(MO1,ISA)+SMM(37,MO1,ISA)                                    
                  TYTP(MO1,ISA)=TYTP(MO1,ISA)+SMM(48,MO1,ISA)                                    
                  TQP(MO1,ISA)=TQP(MO1,ISA)+SMM(49,MO1,ISA)                                      
	              TQPU(MO1,ISA)=TQPU(MO1,ISA)+SMM(108,MO1,ISA)                                        
                  TQN(MO1,ISA)=TQN(MO1,ISA)+SMM(38,MO1,ISA)                                      
                  TYW(MO1,ISA)=TYW(MO1,ISA)+SMM(36,MO1,ISA)                                      
                  TQ(MO1,ISA)=TQ(MO1,ISA)+SMM(13,MO1,ISA)                                        
                  SET(MO1,ISA)=SET(MO1,ISA)+SMM(10,MO1,ISA)                                      
                  TXMX(MO1,ISA)=TXMX(MO1,ISA)+SMM(1,MO1,ISA)                                     
                  TXMN(MO1,ISA)=TXMN(MO1,ISA)+SMM(2,MO1,ISA)                                     
                  TSR(MO1,ISA)=TSR(MO1,ISA)+SMM(3,MO1,ISA)                                       
                  TRHT(MO1,ISA)=TRHT(MO1,ISA)+SMM(33,MO1,ISA)                                    
                  TET(MO1,ISA)=TET(MO1,ISA)+SMM(11,MO1,ISA)                                      
                  IF(ISA==MSA)JJ=IDA                                                             
                  DO K=1,NSM                                                                     
                      SMY(K,ISA)=SMY(K,ISA)+SMM(K,MO1,ISA)                                           
                  END DO                                                                         
                  IF(NDP>0)THEN                                                                  
                      DO K=1,NDP                                                                 
                          DO K1=1,7                                                                      
                              SMYP(K1,K,IDOA(ISA))=SMYP(K1,K,IDOA(ISA))+SMMP(K1,K,&
                              MO1,IDOA(ISA))                               
                          END DO
                          DO K1=10,13                                                                         
                              SMYP(K1,K,IDOA(ISA))=SMYP(K1,K,IDOA(ISA))+SMMP(K1,K,&
                              MO1,IDOA(ISA))
                          END DO
                      END DO                              
                  END IF                                                                         
                  TCVF(MO1,ISA)=TCVF(MO1,ISA)+SMM(25,MO1,ISA)                                    
                  TEI(MO1,ISA)=TEI(MO1,ISA)+SMM(24,MO1,ISA)                                      
                  SMM(25,MO1,ISA)=SMM(25,MO1,ISA)/(SMM(24,MO1,ISA)+1.E-20)                       
                  X1=JCN(ISA)-JCN0(ISA)                                                          
                  SMM(14,MO1,ISA)=SMM(14,MO1,ISA)/(X1+1.E-20)                                    
	              TCN(MO1,ISA)=TCN(MO1,ISA)+SMM(14,MO1,ISA)                                           
                  JCN0(ISA)=JCN(ISA)                                                             
                  ! WRITE MO VALUES AND SUM YEARLY VALUES                                          
                  IF(MO>MO1)CYCLE
                  ACC=0.
                  ADD=0.
                  SUM=0.
                  TOT=0.
                  DO J=1,LC
                      ACC=ACC+STDP(J,ISA)
                      ADD=ADD+STDN(J,ISA)
                      TOT=TOT+UP1(J,ISA)
                      SUM=SUM+UN1(J,ISA)
                  END DO
                  FTN=ZNMN(ISA)+ZNMA(ISA)+ZON(ISA)+ADD+STDON(ISA)+SUM
	              FTP=ZPML(ISA)+ZPMS(ISA)+ZPMA(ISA)+ZPO(ISA)+ZFOP(ISA)+ACC+STDOP&
                  (ISA)+TOT
	              ZTX(23)=ZTX(23)+BTNX(ISA)*WSAX 
	              ZTX(24)=ZTX(24)+BTPX(ISA)*WSAX 
	              !BTNX(ISA)=FTN
	              !BTPX(ISA)=FTP
	              ZTX(25)=ZTX(25)+FTN*WSAX
	              ZTX(26)=ZTX(26)+FTP*WSAX 
	              ZTX(27)=ZTX(27)+BTCX(ISA)*WSAX
	              ZTX(28)=ZTX(28)+ZOC(ISA)*WSAX                                                           
                  SMY(1,ISA)=SMY(1,ISA)/12.                                                      
                  SMY(2,ISA)=SMY(2,ISA)/12.                                                      
                  SMY(3,ISA)=SMY(3,ISA)/12.                                                      
                  SMY(59,ISA)=SMY(59,ISA)/12.                                                    
                  SMY(7,ISA)=SMY(7,ISA)/12.                                                      
                  SMY(8,ISA)=SMY(8,ISA)/12.                                                      
                  SMY(9,ISA)=SMY(9,ISA)/12.                                                      
                  SMY(33,ISA)=SMY(33,ISA)/12.                                                    
                  SMY(34,ISA)=SMY(34,ISA)/12.                                                    
                  SMY(35,ISA)=SMY(35,ISA)/12.                                                    
                  AVRF=AVRF+SMY(4,ISA)                                                           
                  IF(KFL(21)>0)CALL SOCIOA(KDA)                                                  
                  IF(LM(ISA)==0.AND.LUN(ISA)/=35)THEN
                      CALL NLIME                                                                     
	                  IF(TLA>0.)THEN
                          X3=TLA*COL                                                                     
                          COST(ISA)=COST(ISA)+X3                                                         
                          X1=COTL(IAUL(ISA))                                                             
                          X2=X1-COOP(IAUL(ISA))                                                          
                          COST(ISA)=COST(ISA)+X1                                                         
                          CSFX=CSFX+X2                                                                   
                          IF(KFL(31)>0)THEN                                                              
                              WRITE(KW(31),667)ISA,NBSA(ISA),IYR,MO1,KDA,KDC(JJK),X3,X3,TLA                  
                              WRITE(KW(31),666)ISA,NBSA(ISA),IYR,MO1,KDA,TIL(IAUL(ISA)),&
                              KDC(JJK),IHC(IAUL(ISA)),NBE(IAUL(ISA)),NBT(IAUL(ISA)),X1,&
                              X2,FULU(IAUL(ISA))
                          END IF
                      END IF                                                                         
                      SMM(58,MO,ISA)=SMM(58,MO,ISA)+TLA                                              
                      VAR(58,ISA)=TLA                                                                
                      SMY(58,ISA)=TLA
                  END IF                                                                
                  DO K=1,NSM                                                                     
                      SM(K,ISA)=SM(K,ISA)+SMY(K,ISA)                                                 
                  END DO                                                                         
                  DO K=1,NSH                                                                     
                      SMH(K,ISA)=SMH(K,ISA)+SMYH(K,ISA)                                              
                  END DO                                                                         
                  IF(NDP>0)THEN                                                                  
                      DO K=1,NDP                                                                     
                          DO I=1,7                                                                     
                              SMAP(I,K,IDOA(ISA))=SMAP(I,K,IDOA(ISA))+SMYP(I,K,IDOA(ISA))
                              IF(MSA==1)CYCLE
                              SMAP(I,K,NCMD)=SMAP(I,K,NCMD)+WSAX*SMYP(I,K,IDOA(ISA))                                    
                          END DO
                          DO I=10,13                                                                       
                              SMAP(I,K,IDOA(ISA))=SMAP(I,K,IDOA(ISA))+SMYP(I,K,IDOA(ISA))
                              IF(MSA==1)CYCLE
                              SMAP(I,K,NCMD)=SMAP(I,K,NCMD)+WSAX*SMYP(I,K,IDOA(ISA))
                          END DO                                                     
                      END DO                                                                         
                  END IF                                                                         
                  DMX=MIN(PRMT(31),Z(LID(NBSL(ISA),ISA),ISA))                                    
                  IF(DMX>Z(LD1,ISA).AND.LUN(ISA)/=35)CALL TMIX(PRMT(29),DMX,1,1)                     
                  SMY(25,ISA)=SMY(25,ISA)/(SMY(24,ISA)+1.E-20)
                  IF(KFL(27)>0)THEN
	                  WRITE(KW(27),537)ISA,NBSA(ISA),IYR,IY,SMY(13,ISA),SMY(15,ISA),&
                      SMY(16,ISA),SMY(17,ISA),SMY(NDVSS,ISA),SMY(77,ISA),PSTN(1),&
                      (SMYP(J,1,IDOA(ISA)),J=1,7),SMYP(10,1,IDOA(ISA)),SMYP(11,1,&
                      IDOA(ISA))
	                  IF(NDP>1)THEN
	                      DO K=2,NDP 
                              WRITE(KW(27),538)PSTN(K),(SMYP(J,K,IDOA(ISA)),J=1,7),&
                              SMYP(10,K,IDOA(ISA)),SMYP(11,K,IDOA(ISA))
                          END DO
                      END IF
                  END IF                                   
                  II=0                                                                           
                  X1=JCN(ISA)-JCN1(ISA)                                                          
                  SMY(14,ISA)=SMY(14,ISA)/(X1+1.E-20)                                            
                  JCN1(ISA)=JCN(ISA)                                                             
                  II=0                                                                           
                  CALL SPRNT                                                                     
                  IF(KFL(7)>0)THEN
                      DO J=1,NKA                                                                     
                          K=KA(J)                                                                        
                          WRITE(KW(7),1081)ISA,NBSA(ISA),IYR,IY,HED(K),(SMM(K,I,ISA),&
                          I=1,12),SMY(K,ISA),HED(K)                                                             
                      END DO                                                                         
                      DO J=1,NKS                                                                     
                          K=KS(J)                                                                        
                          WRITE(KW(7),1084)ISA,NBSA(ISA),IYR,IY,HEDS(K),(STV(K,I,ISA),&
                          I=1,12),HEDS(K)                                                                      
                      END DO
                  END IF                                                                         
                  !IF(IY/=IPY.OR.KFL(1)==0)GO TO 770                                              
                  IF(IY/=IPY)GO TO 770                                                           
                  IIP=IPYI                                                                       
                  IF(IPD>2.OR.NBSA(ISA)==ISAP)THEN
                      IF(NDP>0)THEN
                          DO K=1,NDP                                                                 
                              IF(K==6.OR.K==1)THEN                                                           
                                  CALL APAGE(1)                                                                  
                                  WRITE(KW(1),1064)ISA,NBSA(ISA),IYR,IY                                          
                                  WRITE(KW(1),'(T48,A)')'PESTICIDE SIMULATION(g/ha)'                             
                              END IF                                                                         
                              WRITE(KW(1),'(T35,A,A16,A)')'-------------------------',&
                              PSTN(K),'-------------------------'                                                      
                              DO L=1,7                                                                       
                                  WRITE(KW(1),145)HEDP(L),(SMMP(L,K,J,IDOA(ISA)),J=1,12),&
                                  SMYP(L,K,IDOA(ISA)),HEDP(L)                                                                        
                              END DO                                                                         
                              DO L=8,9                                                                       
                                  WRITE(KW(1),148)HEDP(L),(SMMP(L,K,J,IDOA(ISA)),J=1,12),HEDP(L)                       
                              END DO
                              DO L=10,13                                                                         
                                  WRITE(KW(1),145)HEDP(L),(SMMP(L,K,J,IDOA(ISA)),J=1,12),&
                                  SMYP(L,K,IDOA(ISA)),HEDP(L)
                              END DO                                                                       
                          END DO
                      END IF    
                      CALL APAGE(1)                                                                  
                      WRITE(KW(1),1060)ISA,NBSA(ISA),IYR,IY                                          
	                  IPG=NKA                                                                             
                      IF(NKA>0)THEN
                          !     PRINTOUT MONTHLY                                                               
                          DO J=1,NKA                                                                     
                              K=KA(J)                                                                        
                              WRITE(KW(1),1130)HED(K),(SMM(K,I,ISA),I=1,12),SMY(K,ISA),HED(K)                
                          END DO
                      END IF
                      IF(IPG==50)THEN
                          CALL APAGE(1)                                                                  
                          WRITE(KW(1),1060)ISA,NBSA(ISA),IYR,IY                                          
                          IPG=1                                                                          
                      END IF          
                      IF(NKS>0)THEN
                          DO J=1,NKS                                                                     
                              K=KS(J)                                                                        
                              WRITE(KW(1),1000)HEDS(K),(STV(K,I,ISA),I=1,12),HEDS(K)                         
                              IPG=IPG+1                                                                      
                              IF(IPG<50)CYCLE                                                                
                              IF(J==NKS)CYCLE                                                                
                              CALL APAGE(1)                                                                  
                              IPG=1                                                                          
                          END DO                                                                         
                      END IF          
                      IF(NJC==0)GO TO 770                                                            
                        DO J=1,NJC                                                                     
                          K=JC(J)                                                                        
                          WRITE(KW(1),1130)HED(K),(SMM(K,I,ISA),I=1,12),SMY(K,ISA),HED(K)                
                          IPG=IPG+1                                                                      
                          IF(IPG<50)CYCLE                                                                
                          CALL APAGE(1)                                                                  
                          IPG=1                                                                          
                      END DO                                                                         
                      GO TO 760                                                                      
                  END IF
            !     PRINTOUT ANNUAL                                                                
                  IF(IPD==0)GO TO 770                                                            
                  WRITE(KW(1),1061)ISA,NBSA(ISA),IYR,IY                                          
                  WRITE(KW(1),1010)IYR,(HED(KA(K)),SMY(KA(K),ISA),K=1,NKA)                       
                  WRITE(KW(1),1011)(HED(JC(K)),SMY(JC(K),ISA),K=1,NJC)                           
                  IF(IPD<=1)GO TO 770                                                            
                  GO TO 751                                                                      
              760 IF(IPD/=4)GO TO 770                                                            
              751 WRITE(KW(1),'(T5,A)')'SOIL DATA'                                               
                  CALL SOLIOP                                                                    
                  CALL SOLIOC
              770 ZTX(1)=ZTX(1)+SMY(11,ISA)*WSAX
                  ZTX(2)=ZTX(2)+SMY(16,ISA)*WSAX
	              ZTX(3)=ZTX(3)+SMY(76,ISA)*WSAX
                  ZTX(4)=ZTX(4)+SMY(18,ISA)*WSAX
	              ZTX(5)=ZTX(5)+SMY(40,ISA)*WSAX
	              ZTX(6)=ZTX(6)+SMY(42,ISA)*WSAX
	              ZTX(7)=ZTX(7)+SMY(43,ISA)*WSAX
	              ZTX(8)=ZTX(8)+SMY(46,ISA)*WSAX
	              ZTX(10)=ZTX(10)+SMY(51,ISA)*WSAX
	              ZTX(11)=ZTX(11)+SMY(53,ISA)*WSAX
	              ZTX(12)=ZTX(12)+SMY(54,ISA)*WSAX
	              ZTX(13)=ZTX(13)+SMY(55,ISA)*WSAX
	              ZTX(14)=ZTX(14)+SMY(56,ISA)*WSAX
	              ZTX(15)=ZTX(15)+SMY(57,ISA)*WSAX
	              ZTX(16)=ZTX(16)+SMY(36,ISA)*WSAX
	              ZTX(17)=ZTX(17)+SMY(134,ISA)*WSAX
	              ZTX(18)=ZTX(18)+SMY(135,ISA)*WSAX
	              ZTX(19)=ZTX(19)+SMY(136,ISA)*WSAX
	              ZTX(20)=ZTX(20)+SMY(4,ISA)*WSAX*RFNC
	              ZTX(31)=ZTX(31)+SMY(17,ISA)*WSAX 
                  ZTX(32)=ZTX(32)+SMY(139,ISA)*WSAX
                  ZTX(33)=ZTX(33)+SMY(44,ISA)*WSAX 	                                                                     
                  IF(KFL(33)>0)WRITE(KW(33),842)ASTN,SITEFILE,TITSO(ISA),TITOP&                  
                  (ISA),IRUN,IRO(ISA),IGN,IYR,IYR,SMY(4,ISA),SMY(10,ISA),SMY(11,ISA)&            
                  ,SMY(13,ISA),SMY(15,ISA),SMY(16,ISA),SMY(17,ISA),SMY(18,ISA),SMY&              
                  (19,ISA),SMY(25,ISA),SMY(26,ISA),SMY(27,ISA),SMY(28,ISA),SMY(29,&              
                  ISA),SMY(30,ISA),SMY(31,ISA),SMY(36,ISA),SMY(37,ISA),SMY(38,ISA),&             
                  SMY(39,ISA),SMY(40,ISA),SMY(49,ISA),SMY(42,ISA),SMY(43,ISA),SMY(44&            
                  ,ISA),SMY(45,ISA),SMY(46,ISA),SMY(47,ISA),SMY(48,ISA),SMY(51,ISA),&            
                  SMY(53,ISA),SMY(54,ISA),SMY(55,ISA),SMY(56,ISA),SMY(57,ISA),SMY(58&            
                  ,ISA),SMY(60,ISA),SMY(66,ISA),OCPD(ISA),ZCO(ISA),ZCOB(ISA),PDPLC&              
                  (ISA),ZPML(ISA),ZNMN(ISA),(CPNM(LY(IRO(ISA),J,ISA)),STD(LY(IRO(&               
                  ISA),J,ISA),ISA),STL(LY(IRO(ISA),J,ISA),ISA),SLAI(LY(IRO(ISA),J,&              
                  ISA),ISA),YLD1(LY(IRO(ISA),J,ISA),ISA),YLD2(LY(IRO(ISA),J,ISA),&               
                  ISA),DMF(LY(IRO(ISA),J,ISA),ISA),YLNF(LY(IRO(ISA),J,ISA),ISA),YLPF&            
                  (LY(IRO(ISA),J,ISA),ISA),J=1,NN)                                               
                  DO K=1,NN                                                                  
                      IF(K==1.OR.K==6)THEN
                          IF(IPD>0)THEN
                              CALL APAGE(1)                                                                  
                              IF(KFL(1)>0)WRITE(KW(1),1060)ISA,NBSA(ISA),IYR,IY
                          END IF
                      END IF                              
                      J=LY(IRO(ISA),K,ISA)                                                           
                      IYH(J,ISA)=IYH(J,ISA)+1                                                        
                !     PRINTOUT CROP MONTHLY                                                          
                      DO I=1,15                                                                      
                          IF(IPD>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),1000)HEDC(I),(SMMC(&                   
                          I,J,K1,ISA),K1=1,12),HEDC(I)                                                   
                      END DO                                                                         
                      DO I=1,6                                                                       
                          TSFC(I,J,ISA)=TSFC(I,J,ISA)+SFCP(I,J,ISA)                                      
                      END DO                                                                         
                      IF(IPD>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),1120)'STRS',(KDT(I,J,&                 
                      ISA),I=1,12),'STRS'                                                            
                      IF(ETG(J,ISA)<1.E-5)THEN                                                       
                          ETG(J,ISA)=ACET(J,ISA)                                                         
                          ACET(J,ISA)=0.                                                                 
                      END IF                                                                         
                      VAL1=YLD1(J,ISA)*PRYG(J)                                                       
                      VAL2=YLD2(J,ISA)*PRYF(J)                                                       
                      IF(VIR(J,ISA)>0.)THEN                                                          
                          VIR(J,ISA)=VIR(J,ISA)-VIR0(ISA)                                                
                      ELSE                                                                           
                          VIR(J,ISA)=VIRT(ISA)                                                           
                          VIR0(ISA)=VIRT(ISA)                                                            
                      END IF                                                                         
                      IF(CAW(J,ISA)<1.E-10)THEN
                          CAW(J,ISA)=AWC(J,ISA)                                                          
                          AWC(J,ISA)=0.                                                                  
                          IF(IDC(J)==NDC(1).OR.IDC(J)==NDC(2).OR.IDC(J)==NDC(4).OR.IDC&                  
                          (J)==NDC(5).OR.IDC(J)==NDC(9))GO TO 820                                        
                          JP(J,ISA)=0                                                                    
                          VIRT(ISA)=0.
                      END IF                                                                   
                      TCAW(J,ISA)=TCAW(J,ISA)+CAW(J,ISA)                                             
                  820 VALF1(ISA)=VALF1(ISA)+VAL1+VAL2
                      IF(ETG(J,ISA)>0.)THEN                                                
                          XX=1000.*YLD1(J,ISA)/ETG(J,ISA)
                      ELSE
                          XX=0.
                      END IF                                       
                      IF(IY==IPY)THEN
                          IF(CSTF(J,ISA)<1.E-10)THEN
                              CSTF(J,ISA)=COST(ISA)                                                          
                              COST(ISA)=0.
                          END IF                                                                   
                    !     PRINTOUT CROP ANNUAL                                                           
                          IF(IPD>0.OR.NBSA(ISA)==ISAP)THEN
                              IF(IDC(J)==NDC(7).OR.IDC(J)==NDC(8).OR.IDC(J)==NDC(10))THEN
                                  X1=.0001*PPL0(J,ISA)
                              ELSE
                                  X1=PPL0(J,ISA)
                              END IF
                              WRITE(KW(1),1100)CPNM(J),YLD1(J,ISA),YLD2(J,ISA),DMF&
                              (J,ISA),XX,YLNF(J,ISA),YLPF(J,ISA),CAW(J,ISA),X1,CSTF&
                              (J,ISA),VAL1,VAL2,VIR(J,ISA),EK(ISA),WK(ISA),THK(ISA)                                       
                              WRITE(KW(1),1020)(SFCP(I,J,ISA),I=1,6)
                          END IF
                      END IF                                         
                      IF(KFL(24)>0)THEN
                          X2=.001*ZOC(ISA)
                          WRITE(KW(24),98)ISA,NBSA(ISA),IYR,IY,CPNM(J),YLD1(J,ISA),&
                          YLD2(J,ISA),DMF(J,ISA),(SFCP(I,J,ISA),I=1,6),ZNMN(ISA),&               
                          ZPML(ISA),PDPLC(ISA),X2,OCPD(ISA),RSD(LD1,ISA),ARSD(ISA),&
                          SMY(54,ISA),SMY(55,ISA),SMY(53,ISA),SMY(57,ISA),SMY(56,ISA)
                      END IF                         
                      TDM(J,ISA)=TDM(J,ISA)+DMF(J,ISA)                                               
                      TYL1(J,ISA)=TYL1(J,ISA)+YLD1(J,ISA)                                            
                      TYL2(J,ISA)=TYL2(J,ISA)+YLD2(J,ISA)                                            
                      TYLN(J,ISA)=TYLN(J,ISA)+YLNF(J,ISA)                                            
                      TYLP(J,ISA)=TYLP(J,ISA)+YLPF(J,ISA)
                      TFTN(J,ISA)=TFTN(J,ISA)+FRTN(J,ISA)
                      TFTP(J,ISA)=TFTP(J,ISA)+FRTP(J,ISA)
                      X1=YLD1(J,ISA)+YLD2(J,ISA)
                      IF(X1>0.)THEN                                            
                          TCPA(J)=TCPA(J)+WSAX
                          IF(IDC(J)==9)X1=YLD1(J,ISA)
                          TCPY(J)=TCPY(J)+WSAX*X1
                      END IF
                      ZTX(21)=ZTX(21)+YLNF(J,ISA)*WSAX
                      ZTX(22)=ZTX(22)+YLPF(J,ISA)*WSAX
                      IF(YLNF(J,ISA)>0.)THEN                                                         
                          NYLN(J,ISA)=NYLN(J,ISA)+1                                                      
	                      XX=NYLN(J,ISA)                                                                      
                          X1=XX+1.                                                                       
                          ! UNA(J,ISA)=(UNA(J,ISA)+PRMT(28)*TYLN(J,ISA))/XX                                
                          UNA(J,ISA)=(UNA(J,ISA)*XX+1000.*DMF(J,ISA)*BN(3,J))/X1                         
                      END IF                                                                         
                      TRD(J,ISA)=TRD(J,ISA)+RDF(J,ISA)                                               
                      THU(J,ISA)=THU(J,ISA)+HUF(J,ISA)                                               
                      TETG(J,ISA)=TETG(J,ISA)+ETG(J,ISA)                                             
                      CST1(ISA)=CST1(ISA)+CSTF(J,ISA)                                                
                      YLNF(J,ISA)=0.                                                                 
                      YLPF(J,ISA)=0.                                                                 
                      DMF(J,ISA)=0.
                      FRTN(J,ISA)=0.
                      FRTP(J,ISA)=0.                                                                  
                      VIR(J,ISA)=0.                                                                  
                      CAW(J,ISA)=0.                                                                  
                      RDF(J,ISA)=0.                                                                  
                      HUF(J,ISA)=0.                                                                  
                      CSTF(J,ISA)=0.                                                                 
                      ETG(J,ISA)=0.                                                                  
                      IF(IDC(J)/=NDC(2).AND.IDC(J)/=NDC(5))DM1(J,ISA)=0.                             
                      DO I=1,6                                                                       
                          SFCP(I,J,ISA)=0.                                                               
                      END DO
                  END DO                                                                         
                  IF(KFL(4)>0)WRITE(KW(4),99)ISA,NBSA(ISA),IYR,IY,(SMY(KY(J1),&                  
                  ISA),J1=1,NKY)                                                                 
                  IF(KFL(39)>0.AND.IEXT(ISA)==0)WRITE(KW(39),472)NBSA(ISA),IGC,&                 
                  IYR,RWSA(ISA),(SMYH(K,ISA),K=1,NSH)                                            
                  DO I=1,12                                                                      
                      XTX(I)=XTX(I)+SMM(4,I,ISA)                                                     
                  END DO                                                                         
                  XTX(13)=XTX(13)+SMY(4,ISA)
                  IF(ISA==1)THEN
                      ZTX(29)=SMY(31,1)
                      ZTX(30)=BTCX(1)
                  END IF
                  RSVF(ISA)=RSV(ISA)                                                             
                  ARSD(ISA)=0.                                                                   
                  RSYF(ISA)=STV(12,12,ISA)                                                       
                  DO K=1,NKS                                                                     
                      DO I=1,12                                                                    
                          STV(K,I,ISA)=0.                                                              
                      END DO                                                                       
                  END DO                                                                         
                  KT(ISA)=1                                                                      
                  KP1(ISA)=1                                                                     
                  KF1(ISA)=1                                                                     
                  KC(ISA)=0                                                                      
                  NDFA(ISA)=0                                                                    
                  IF(NDP>0)THEN
                      DO K=1,NDP
                          TLHP(K)=TLHP(K)+SMYP(3,K,IDOA(ISA))*WSAX
		                  TSSP(K)=TSSP(K)+SMYP(4,K,IDOA(ISA))*WSAX
		                  TDSP(K)=TDSP(K)+SMYP(7,K,IDOA(ISA))*WSAX
		                  TDFP(K)=TDFP(K)+SMYP(6,K,IDOA(ISA))*WSAX
		                  TDRP(K)=TDRP(K)+SMYP(10,K,IDOA(ISA))*WSAX
		                  TRFP(K)=TRFP(K)+SMYP(11,K,IDOA(ISA))*WSAX                                                                                                                                          
                          DO K1=1,13                                                                    
                              SMYP(K1,K,IDOA(ISA))=0.                                                            
                              DO I=1,12                                                                
                                  SMMP(K1,K,I,IDOA(ISA))=0.                                                      
                              END DO                                                                   
                          END DO                                                                       
                      END DO
                  END IF                                                                                   
              END DO                                                                       
              MO1=MO                                                                         
        !     CALL NCONT                                                                     
              IF(MO1==12)EXIT
          END DO
      END DO                                                                
      IF(NDP>0)THEN
          D15=SRCH(15,NCMD)-D150
	      TPRK=TPRK/RWSX
	      DO K=1,NDP
		      TLHP(K)=TLHP(K)/RWSX
		      TSSP(K)=TSSP(K)/RWSX
		      TDSP(K)=TDSP(K)/RWSX
  		      TDFP(K)=TDFP(K)/RWSX
		      TDRP(K)=TDRP(K)/RWSX
		      TRFP(K)=TRFP(K)/RWSX
	      END DO
          WRITE(KW(1),18)                                                                
          WRITE(KW(1),1062)IYR,IY                                                        
          DO K=1,NDP                                                                  
              WRITE(KW(1),3)PSTN(K)                                                          
              SUM=0.                                                                         
              DO I=1,12                                                                    
                  SMMRP(5,K,I)=SMMRP(5,K,I)/RWSX                                         
                  SUM=SUM+SMMRP(5,K,I)                                                         
                  SMRP(5,K,I)=SMRP(5,K,I)+SMMRP(5,K,I)                                         
              END DO                                                                       
              WRITE(KW(1),2)HEDP(1),(SMMRP(5,K,I),I=1,12),SUM,HEDP(1)                        
              DO I=1,12                                                                    
                  SMMRP(5,K,I)=0.                                                              
              END DO
              I1=0
              ZTX(31)=ZTX(31)/RWSX
              IF(KFL(27)>0)THEN
                  X1=SMYRP(3,K)
                  IF(K==1)THEN
	                  WRITE(KW(27),537)I1,I1,IYR,IY,SMYR(9),D15,TPRK,&
	                  ZTX(31),SMYR(11),TAC,PSTN(1),SUM,X1,TLHP(K),&
	                  TSSP(K),SMYRP(4,1),TDFP(K),TDSP(K),TDRP(K),TRFP(K)
                  ELSE
       	              WRITE(KW(27),538)PSTN(K),SUM,X1,TLHP(K),TSSP(K),&
       	              SMYRP(4,K),TDFP(K),TDSP(K),TDRP(K),TRFP(K)
                  END IF
              END IF
              TAC=0.                                                                       
              DO J=1,2                                                                     
                  J1=J+2                                                                       
                  SMYRP(J,K)=SMYRP(J,K)/RWSX                                             
                  DO I=1,12                                                                
                      SMMRP(J,K,I)=SMMRP(J,K,I)/RWSX                                      
                      SMRP(J,K,I)=SMRP(J,K,I)+SMMRP(J,K,I)                                     
                      SMRP(J1,K,I)=SMRP(J1,K,I)+SMMRP(J1,K,I)                                  
                  END DO                                                                   
                  WRITE(KW(1),2)HDRP(J),(SMMRP(J,K,I),I=1,12),SMYRP(J,K),HDRP(J)                 
                  WRITE(KW(1),2)HDRP(J1),(SMMRP(J1,K,I),I=1,12),SMYRP(J1,K),HDRP(J1)             
                  SMYRP(J,K)=0.                                                                  
                  SMYRP(J1,K)=0.                                                                 
              END DO
          END DO              
          D150=SRCH(15,NCMD)
      END IF              
      CALL APAGE(0)
      DO I=1,28
          ZTX(I)=ZTX(I)/RWSX
      END DO
      DO I=1,13
          XTX(I)=XTX(I)/XSA
      END DO
      DO J=1,10                                                                      
          SMYR(J)=SMYR(J)/RWSX                                                     
      END DO
      ZTX(32)=ZTX(32)/RWSX                                                                         
      ZTX(33)=ZTX(33)/RWSX                                                                         
      YTP(1)=SRCH(15,NCMD)-YTP(1)
      YTP(2)=SRCH(16,NCMD)-YTP(2)
      YTP(3)=SRCH(19,NCMD)-YTP(3)
      YTP(4)=SRCH(20,NCMD)-YTP(4)
      YTP(5)=SRCH(17,NCMD)-YTP(5)
      YTP(6)=SRCH(21,NCMD)-YTP(6)
      YTP(7)=SRCH(18,NCMD)-YTP(7)
      YTP(8)=SRCH(22,NCMD)-YTP(8)
      YTP(9)=SRCH(25,NCMD)-YTP(9)
      YTP(10)=SRCH(11,NCMD)-YTP(10)
      WRITE(KW(8),936)IYR,XTX(13),ZTX(1),SMYR(1),YTP(1),YTP(2),YTP(3),&
      YTP(4),ZTX(2),ZTX(4),SMYR(2),SMYR(13),ZTX(16),YTP(10),YTP(5),YTP&
      (6),YTP(7),SMYR(14),ZTX(17),YTP(8),ZTX(5),ZTX(6),ZTX(8),ZTX(7),ZTX&
      (11),ZTX(12),ZTX(13),SMYR(17),SMYR(15),ZTX(18),ZTX(10),ZTX(14),ZTX&
      (15),ZTX(3),YTP(9),(ZTX(I),I=19,28),ZTX(30),ZOC(1),ZTX(29),ZTX(32),&
      ZTX(33)
      YTP(1)=SRCH(15,NCMD)
      YTP(2)=SRCH(16,NCMD)
      YTP(3)=SRCH(19,NCMD)
      YTP(4)=SRCH(20,NCMD)
      YTP(5)=SRCH(17,NCMD)
      YTP(6)=SRCH(21,NCMD)
      YTP(7)=SRCH(18,NCMD)
      YTP(8)=SRCH(22,NCMD)
      YTP(9)=SRCH(25,NCMD)
      YTP(10)=SRCH(11,NCMD)
      IF(IPD>0)CALL APAGE(0)                                                         
      IF(KFL(28)>0)THEN
          DO I=1,12                                                                      
	          CY=.1*SMMR(13,I)/(SMMR(11,I)+.1)
	          X3=MIN(.95,DRSW*PRMT(63))                                                               
              X4=MIN(.95,DRSW*PRMT(75))                                                               
	          ZTZ(1)=SMMR(12,I)                                                                   
	          ZTZ(2)=SMMR(13,I)*DRSW                                                              
              ZTZ(3)=SMMR(14,I)*X3                                                           
              ZTZ(4)=SMMR(15,I)*X4                                                           
              ZTZ(5)=SMMR(16,I)                                                              
              ZTZ(6)=SMMR(17,I)
              DO K=1,NDP
                  SMMRP(4,K,I)=X4*SMMRP(4,K,I)
              END DO                                                              
  	          WRITE(KW(28),895)IYR,I,(ZTZ(J),J=1,6),(SMMRP(3,K,I),&
  	          SMMRP(4,K,I),K=1,NDP)
  	          DO K=1,6
  	              SMSW(K)=SMSW(K)+ZTZ(K)
  	          END DO
  	                                                  
          END DO                                                                         
    !     AVRF=AVRF/XSA                                                                  
    !     WRITE(KW(4),26)AVRF
      END IF
      SMMRP=0.
      AVRF=0.                                                                        
      IF(KFL(35)>0.AND.MSA>1)THEN
          NTX=0                                                                          
	      DO IDO=1,NCMD	                                                                      
	          IF(ICDT(IDO)==3)CYCLE                                                               
	          DO MO=1,12                                                                          
		          SMMH(7,MO,IDO)=1.E5*SMMH(6,MO,IDO)/(SMMH(2,MO,IDO)*RWSA&
		          (IDO)+1.E-10)                                                                       
	              IF(SMMH(35,MO,IDO)<1.E-10)SMMH(35,MO,IDO)=SMMH(35,MO,IDOA&
                  (NISA(IDNB(IDO))))                                                                 
              END DO                                                                         
              DO K=1,NSH                                                                 
                  DO MO=1,12                                                                   
                      IF(K==7)CYCLE                                                                
		              SMYH(K,IDO)=SMYH(K,IDO)+SMMH(K,MO,IDO)                                             
		              SMH(K,IDO)=SMH(K,IDO)+SMMH(K,MO,IDO)                                               
		              IF(K<5.OR.K>32)THEN                                                                
                          N2=NC(MO+1)                                                                  
                          IF(MO==2)N2=N2-NYD                                                           
		                  XM=N2-NC(MO)                                                                       
		                  SMMH(K,MO,IDO)=SMMH(K,MO,IDO)*RWSA(IDO)/(8640.*XM)                                 
		              END IF                                                                             
		          END DO
		      END DO                                                                           
              SMYH(7,IDO)=1.E5*SMYH(6,IDO)/(SMYH(2,IDO)*RWSA(IDO)+1.E-10)                    
		      DO K=1,NSH                                                                         
		          IF(K<5.OR.K>32)THEN                                                                
		              X1=366-NYD                                                                         
		              SMYH(K,IDO)=SMYH(K,IDO)*RWSA(IDO)/(8640.*X1)                                       
		          END IF                                                                             
		      END DO                                                                             
          END DO                                                                         
	      IF(IPD<3)THEN                                                                       
	          DO I=1,MSA                                                                          
	              I1=NBSA(IBSA(I))                                                                    
	              I2=NISA(I1)                                                                         
	              IF(IEXT(I2)>0)THEN
	                  II=IDOA(I2)
	                  X2=WSA(I2)
	              ELSE
	                  II=IDOA(I2)-1
	                  X2=RWSA(II)+WSA(I2)
	              END IF                                                                      
                  IF(NTX(II)>0)CYCLE                                                             
                  WRITE(KW(35),472)I1,II,IYR,X2,SMYH(1,II),SMYH(2,II),SMYH&                
                  (33,II),SMYH(34,II),(SMYH(K,II),K=3,NSH-2),(PCTH(J,II),&
                  J=1,NSZ),(PCT(J,II),J=1,NSZ)                                                            
	              NTX(II)=1                                                                           
	          END DO                                                                              
	      END IF                                                                              
	      IF(KFL(37)>0)THEN                                                                   
	          DO MO=1,12                                                                          
	              DO I=1,NCMO                                                                       
	                  II=ICMO(I)                                                                        
	                  WRITE(KW(37),472)II,IYR,MO,RWSA(II),SMMH(1,MO,II),&
	                  SMMH(2,MO,II),SMMH(33,MO,II),SMMH(34,MO,II),(SMMH&
	                  (K,MO,II),K=3,NSH-2)                                                                       
                  END DO                                                                       
              END DO                                                                         
          END IF                                                                         
	      IF(IPD>2.AND.IPD<6)THEN                                                             
	          DO MO=1,12                                                                          
                  NTX=0                                                                          
                  DO I=1,MSA                                                                         
	                  I1=NBSA(IBSA(I))                                                                   
	                  I2=NISA(I1)                                                                        
	                  IF(IEXT(I2)>0)THEN
	                      II=IDOA(I2)
	                      X2=WSA(I2)
	                  ELSE
	                      II=IDOA(I2)-1
	                      X2=RWSA(II)+WSA(I2)
	                  END IF                                                                     
	                  IF(NTX(II)>0)CYCLE                                                                 
	                  WRITE(KW(35),472)I1,IYR,MO,X2,SMMH(1,MO,II),SMMH(2,MO,II),&
	                  SMMH(33,MO,II),SMMH(34,MO,II),(SMMH(K,MO,II),K=3,NSH-2),&
	                  (PCTH(J,II),J=1,NSZ),(PCT(J,II),J=1,NSZ)                               
	                  NTX(II)=1                                                                          
	              END DO                                                                             
	          END DO                                                                              
	      END IF                                                                              
	      SMMH=0.                                                                             
	      SMYH=0.
      END IF	                                                                                   
      IF(KFL(34)>0)THEN
          IF(IPD>2.AND.IPD<6)THEN                                                        
	          DO MO=1,12                                                                          
	              DO I=1,MSA                                                                         
	                  I1=NBSA(IBSA(I))                                                                   
	                  I2=NISA(I1)                                                                        
	                  II=IDOA(I2)
	                  X1=.001*ZOC(I2)	                                                                       
	                  WRITE(KW(34),471)I1,IYR,MO,WSA(I2),(SMM(KY(K),MO,I2),K=1,&
	                  NKY),X1,(PCTH(J,II),J=1,NSZ),(PCT(J,II),J=1,NSZ),(YLD1&
	                  (LY(IRO(I2),J,I2),I2),YLD2(LY(IRO(I2),J,I2),I2),(SMMC&
	                  (K,LY(IRO(I2),J,I2),MO,I2),K=1,15),CPNM(LY(IRO(I2),J,I2))&
	                  ,J=1,NCP(IRO(I2),I2))
                  END DO                                                                        
              END DO      	                                                                  
          END IF                                                                         
          DO I=1,MSA                                                                     
	          I1=NBSA(IBSA(I))                                                                    
	          I2=NISA(I1)	                                                                        
	          II=IDOA(I2)
	          X1=.001*ZOC(I2)	                                                                        
     	      IF(IPD>0.AND.IPD<3)WRITE(KW(34),471)I1,IGC,IYR,WSA(I2),(SMY&
     	      (KY(K),I2),K=1,NKY),X1,(PCTH(J,II),J=1,NSZ),(PCT(J,II),J=1,NSZ)&
     	      ,(YLD1(LY(IRO(I2),J,I2),I2),YLD2(LY(IRO(I2),J,I2),I2),(SMMC&
     	      (K,LY(IRO(I2),J,I2),12,I2),K=1,15),CPNM(LY(IRO(I2),J,I2)),&
     	      J=1,NCP(IRO(I2),I2))                              
          END DO
      END IF                                                                                   
      SMM=0.                                                                         
      SMMC=0.                                                                        
      SMY=0.                                                                         
      YLD1=0.                                                                        
      YLD2=0.                                                                        
      WRITE(KW(1),18)                                                                
      WRITE(KW(1),1062)IYR,IY                                                        
      WRITE(KW(1),1130)HED(4),(XTX(I),I=1,13),HED(4)                                 
      IF(KFL(15)>0)WRITE(KW(15),6)IYR,(XTX(I),I=1,13)                                
      IF(KFL(17)>0)WRITE(KW(17),4083)IYR,XTX(13),(SMYR(J),SMYR(J+10),&               
      J=1,10)                                                                        
      DO J=1,10                                                                   
          J1=J+10                                                                        
          DO I=1,12                                                                      
              SMMR(J,I)=SMMR(J,I)/RWSX                                                 
              SMR(J1,I)=SMR(J1,I)+SMMR(J1,I)                                                 
              SMR(J,I)=SMR(J,I)+SMMR(J,I)                                                    
          END DO                                                                         
          WRITE(KW(1),1130)HEDR(J),(SMMR(J,I),I=1,12),SMYR(J),HEDR(J)                    
          WRITE(KW(1),1130)HEDR(J1),(SMMR(J1,I),I=1,12),SMYR(J1),HEDR(J1)                
          IF(KFL(15)>0)WRITE(KW(15),1130)HEDR(J1),(SMMR(J1,I),I=1,12),&                  
          SMYR(J1),HEDR(J1)                                                              
          SMYR(J)=0.                                                                     
          SMYR(J1)=0.                                                                    
          DO I=1,12                                                                      
              SMMR(J,I)=0.                                                                   
              SMMR(J1,I)=0.                                                                  
          END DO
      END DO                                                                         
      DO I=1,12                                                                      
          SMR(21,I)=SMR(21,I)+SMMR(21,I)                                                 
          SMMR(21,I)=0.                                                                  
      END DO                                                                         
	  SMYR(21)=0.                                                                         
      TMAF=TMAF+TMAP                                                                 
      TMAP=0.                                                                        
      IBD=1                                                                          
      MO=1                                                                           
      IYR=IYR+1                                                                      
      DO IOW=1,NBON                                                                  
	      IF(NHRD(IOW)==0)CYCLE                                                               
          DO IHD=1,NHRD(IOW)                                                             
	          IF(NYHO(IHD,IOW)==0)CYCLE                                                           
	          IYHO(IHD,IOW)=IYHO(IHD,IOW)+1                                                       
	          IF(IYHO(IHD,IOW)>NYHO(IHD,IOW))IYHO(IHD,IOW)=1                                      
	      END DO                                                                              
	  END DO                                                                              
  	  IYX=IYX+1                                                                         
      NYD=1                                                                          
      IPY=IPY+IIP                                                                    
	  IF(ISW<=1.OR.ISW==4)THEN                                                                      
	      DO ISA=1,MSA                                                                    
              XX=0.                                                                          
              DO I=1,NBSL(ISA)                                                             
                  J=LID(I,ISA)                                                                 
                  X1=1000.*(Z(J,ISA)-XX)                                                       
                  Y1=.1*WOC(J,ISA)/WT(J,ISA)
                  IF(ISW<4)THEN
                      CALL SWRTNR(CLA(J,ISA),SAN(J,ISA),Y1,X2,X3)
                  ELSE
                      CALL SWNN(CLA(J,ISA),SAN(J,ISA),Y1,X2,X3)
                  END IF                                                   
                  XY=1.-ROK(J,ISA)*.01                                                         
                  S15(J,ISA)=X2*X1*XY                                                          
                  FC(J,ISA)=X3*X1*XY                                                           
                  CALL SPOFC(J)                                                                
                  XX=Z(J,ISA)                                                                  
              END DO 
          END DO                                                                      
      END IF                                                                         
!     IF(IGSD<0.AND.IRO(ISA)>=NRO(ISA))THEN
!         NGN=NGN0                                                                       
!         IYR=IYR0                                                                       
!         IGSD=IGSD+NRO(ISA)
!     END IF                                                             
      CALL ALPYR(IYR,NYD,LPYR)                                                       
      IF(ISTA==0)GO TO 940                                                           
      DO ISA=1,MSA                                                                   
          DO L=1,NBSL(ISA)                                                             
              J=LID(ISL,ISA)                                                               
              WPMA(J,ISA)=SOL(1,J,ISA)                                                     
              WON(J,ISA)=SOL(2,J,ISA)                                                      
              WPO(J,ISA)=SOL(4,J,ISA)                                                      
              WPMS(J,ISA)=SOL(5,J,ISA)                                                     
              WOC(J,ISA)=SOL(7,J,ISA)                                                      
          END DO                                                                       
      END DO                                                                         
  940 CONTINUE                                                                       
      L=KND                                                                          
      DO I=1,NWTH                                                                    
          L=L+1                                                                          
          REWIND KR(L)                                                                   
      END DO                                                                         
      IY=NBYR+1                                                                      
      RETURN                                                                         
    2 FORMAT(1X,A4,13E12.4,2X,A4)                                                    
    3 FORMAT(6X,A8)                                                                  
    6 FORMAT(1X,I4,13F9.2)
    8 FORMAT(5X,'IDO= ',I8,2X,'IDSA= ',I8,2X,'HYD VOL= ',F8.3,' mm',2X,&
      'PEAK RATE= ',E13.5,' m3/s'/10X,'TP= ',F7.2,' h')
   12 FORMAT(5X,A2,2I8,3I4,5F10.2)                                                                    
   18 FORMAT(//1X,'______________ANNUAL WATERSHED TABLE_________________&            
      _'/T10,'SUM OF SUBAREA OUTFLOWS/TOTAL WATERSHED OUTFLOW'/)                     
   26 FORMAT(10X,100F10.2)
   35 FORMAT(1X,'!!!!!',5E16.6)                                                              
   50 FORMAT(/T5,'YEAR ',I4,' OF ',I4,/)                                             
   98 FORMAT(1X,2I8,1X,I4,1X,I4,1X,A4,60F10.2)                                       
   99 FORMAT(1X,2I8,1X,2I4,60F10.2)                                                  
  105 FORMAT(16X,3F8.0)                                                              
  123 FORMAT(I5,1X,A16,F12.0,3F8.0,F10.0,F8.0)    
  145 FORMAT(1X,A4,13F12.5,2X,A4)                                                    
  148 FORMAT(1X,A4,12F12.5,14X,A4)                                                   
  154 FORMAT(1X,2I8,1X,I4,2I2,100F10.4)
  471 FORMAT('BIGSUB',I4,I9,I5,1X,E9.3,9F10.2,F10.4,25F10.2,6F10.4,&
      5(17F10.2,A10))                                                                                   
  472 FORMAT('REACH',I5,I9,I6,50E12.4)
  537 FORMAT(1X,4I5,6F8.2,1X,A16,11E13.6)
  538 FORMAT(70X,A16,11E13.6) 
  666 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)                   
  667 FORMAT(1X,2I8,1X,I4,2I2,2X,'LIME',12X,I6,6X,'   9',8X,F10.2,10X,&              
      2F10.2)                                                                        
  842 FORMAT(1X,A20,2X,3A12,5I4,44F8.2,10(2X,A4,8F8.2))                               
  892 FORMAT(2X,'BUY/SELL',I4,2I2,2X,'IDON=',I4,2X,'HRD#=',I3,2X,'SIZE='&            
      ,I8,' HD')                                                                     
  894 FORMAT(1X,3I4,100F10.4)                                                        
  895 FORMAT(1X,I4,1X,I4,1X,100(1X,F10.3))                                            
  902 FORMAT(1X,I4,1X,I4,1X,100(1X,E16.10))                                           
  903 FORMAT(1X,2I8,I5,2I3,1X,4F10.2,12(1X,A16,12F10.4))
  904 FORMAT(1X,2I8,I5,2I3,41X,1X,A16,12F10.4)                                       
! 906 FORMAT(15X,10F10.3)                                                                                           
! 921 FORMAT(1X,'ANNUAL MANURE BAL'/I3,8D16.8)                                       
  909 FORMAT(1X,2I8,1X,I4,2I2,1X,4F8.2,15(3X,A16,3E12.5))                            
  936 FORMAT(5X,I4,38F10.3,8F10.0,2F10.2,F10.1)     
 1000 FORMAT(1X,A4,12E12.4,14X,A4)                                                   
 1010 FORMAT(//I5,9(2X,A4,F8.2)/(5X,9(2X,A4,F8.2)))                                  
 1011 FORMAT(5X,9(2X,A4,1X,E10.4))                                                       
 1020 FORMAT(T10,'STRESS DAYS(BIOMASS)--  WATER=',F5.1,2X,'N=',F5.1,2X,&
      'P=',F5.1,2X,'TEMP=',F5.1,2X,'AIR=',F5.1,2X,'SALT=',F5.1)                           
 1060 FORMAT(35X,'SUBAREA NO=',I8,' ID=',I8,' YR=',I4,' YR#=',I4/T11,&               
      'JAN',9X,'FEB',9X,'MAR',9X,'APR',9X,'MAY',9X,'JUN',9X,'JUL',9X,&               
      'AUG',9X,'SEP',9X,'OCT',9X,'NOV',9X,'DEC',9X,' YR')                            
 1061 FORMAT(35X,'SUBAREA NO=',I8,' ID=',I8,' YR=',I4,' YR#=',I4)                    
 1062 FORMAT(35X,' YR=',I4,' YR#=',I4/T11,'JAN',9X,'FEB',9X,'MAR',9X,&               
      'APR',9X,'MAY',9X,'JUN',9X,'JUL',9X,'AUG',9X,'SEP',9X,'OCT',9X,&               
      'NOV',9X,'DEC',9X,' YR')                                                       
 1064 FORMAT(35X,'SUBAREA NO=',I8,' ID=',I8,' YR=',I4,' YR#=',I4/T14,&               
      'JAN',9X,'FEB',9X,'MAR',9X,'APR',9X,'MAY',9X,'JUN',9X,'JUL',9X,&               
      'AUG',9X,'SEP',9X,'OCT',9X,'NOV',9X,'DEC',9X,' YR')                            
 1070 FORMAT(14X,15F6.0)                                                              
!1080 FORMAT(1X,A4,13F12.2,2X,A4)                                                    
 1081 FORMAT(1X,2I8,1X,2I4,1X,A4,13E13.5,2X,A4)                                      
 1084 FORMAT(1X,2I8,1X,2I4,1X,A4,12E13.5,11X,A4)                                     
 1100 FORMAT(2X,A4,1X,'YLD=',F5.1,'/',F5.1,'t/ha',2X,'BIOM=',F6.2,'t/ha'&            
      ,2X,'WUEF=',F6.2,'t/mm',2X,'YLN=',F5.0,'kg/ha',2X,'YLP=',F5.0,&                
      'kg/ha',2X,'CAW=',F7.0,' mm',2X,'POP=',F9.4,'P/m2'/4X,'COST=',&
      F7.0,'$/ha',2X,'RTRN=',F5.0,'/',F5.0,'$/ha',2X,'IRGA=',F5.0,' mm',&
      2X,'EK=',F5.2,2X,'WK=',F5.2,2X,'THK=',F5.0,' mm')                                  
 1120 FORMAT(1X,A4,12I9,11X,A4)                                                      
 1130 FORMAT(1X,A4,13E12.4,2X,A4)                                                    
 4083 FORMAT(1X,I4,F10.0,20F10.2)                                                    
      END                                                                            
