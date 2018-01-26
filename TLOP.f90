      SUBROUTINE TLOP(CSTX,COX,JRT)
!     APEX0806
!     THIS SUBPROGRAM CONTROLS ALL TILLAGE OPERATIONS INCLUDING PLANTING
!     & HARVESTING.
      USE PARM 
      FNPP(X)=DMLA(JJK)*X/(X+EXP(PPCF(1,JJK)-PPCF(2,JJK)*X))
      JRT=0
      II=IHC(JT1)
      LD1=LID(1,ISA)
      NN=NCP(IRO(ISA),ISA)
      N1=MAX(1,NN)
      X1=CND(IRO(ISA),KT(ISA),ISA)
      IF(ABS(X1-CN0(ISA))>0.)THEN
          X2=SMX(ISA)
          CALL HCNSLP(X1,X3)
          CN0(ISA)=X1
          CN2(ISA)=X1
          SCI(ISA)=SMX(ISA)*SCI(ISA)/X2
      END IF
      SELECT CASE(II)
          CASE(1,2,3,19,22)
              GO TO 10
          CASE(5)
              GO TO 61
          CASE(6)
              GO TO 53
          CASE(7,8)
              GO TO 57
          CASE(10)
              CSTX=-CSTX*YLD1(JJK,ISA)/(1.-WCY(JJK))
              COX=CSTX
              GO TO 57
          CASE(11)
              CSTX=-CSTX*YLD(JJK)/(1.-WCY(JJK))
              COX=CSTX
              GO TO 57
          CASE(12,13)
              IF(ICUS(JT1)==0)GO TO 57
              CSTX=-CSTX*YLD(JJK)/(1.-WCY(JJK))
              COX=CSTX
              GO TO 57
          CASE(14)
              GO TO 42
          CASE(21)
              GO TO 43
          CASE(23)
              GO TO 63
          CASE(24)
              GO TO 64
          CASE DEFAULT
              GO TO 6
      END SELECT
   53 IDRL(ISA)=1
   61 ISL=LID(2,ISA)
      DO K=1,NN
          I2=LY(IRO(ISA),K,ISA)
          IF(JH(IRO(ISA),KT(ISA),ISA)==KDC(I2))GO TO 27
      END DO
      GO TO 26
   27 IF(KGO(I2,ISA)>0)GO TO 26
      ZX=0.
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          Z1=Z(ISL,ISA)
          ZZ=.5*(ZX+Z1)
          IF(ZZ>=.075)EXIT
   	      ZX=Z1 
      END DO
      LRD(ISA)=I 
      IF(STMP(ISL,ISA)<TBSC(I2)+2.)THEN
          KOMP(KT(ISA),ISA)=0
          JRT=1
          RETURN 
      END IF
      AWC(JJK,ISA)=RZSW(ISA)
      IGO(ISA)=IGO(ISA)+1
!	  TCPA(I2)=TCPA(I2)+WSA(ISA)
	  KC(ISA)=1
	  DO WHILE(JE(KC(ISA),ISA)<12)
          KC(ISA)=KC(ISA)+1
          IF(KC(ISA)>NN)KC(ISA)=1
      END DO
      JE(KC(ISA),ISA)=I2
      JJK=I2
      KGO(JJK,ISA)=1
      JP(JJK,ISA)=0
	  IYH(JJK,ISA)=1
      SWH(JJK,ISA)=0.
      SWP(JJK,ISA)=0.
      ACET(JJK,ISA)=0.
      XDLAI(JJK)=DLAI(JJK)
	  IF(CPNM(JJK)=='FALW')NCR(JJK,ISA)=NCR(JJK,ISA)+1
      WCYD=.3
      STDO(ISA)=STDO(ISA)+STD(JJK,ISA)
      STDON(ISA)=STDON(ISA)+STDN(JJK,ISA)
      STDOP(ISA)=STDOP(ISA)+STDP(JJK,ISA)
      STD(JJK,ISA)=0.
      STDN(JJK,ISA)=0.
      STDP(JJK,ISA)=0.
      STDL(JJK,ISA)=0.
      RD(JJK,ISA)=TLD(JT1)
      HU(JJK,ISA)=0.
      DM(JJK,ISA)=SDW(JJK)*5.E-4
      DM1(JJK,ISA)=DM(JJK,ISA)
      RW(JJK,ISA)=.4*DM(JJK,ISA)
      RWT(ISL,JJK,ISA)=RW(JJK,ISA)
      ROSP(ISA)=RIN(JT1)
      PPL0(JJK,ISA)=POP(JJK,IHU(JJK,ISA),ISA)
      XLAI(JJK,ISA)=FNPP(PPL0(JJK,ISA))
      DMLX(JJK)=XLAI(JJK,ISA)
	  X1=SDW(JJK)*CSTS(JJK)
      COST(ISA)=COST(ISA)+X1
      LRD(ISA)=MAX(2,LRD(ISA))
      JPL(JJK,ISA)=1
      IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),32)ISA,NBSA(ISA),IYR,&
      MO,KDA,CPNM(JJK),CV(ISA)
      IF(KFL(31)>0)WRITE(KW(31),49)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JT1),&
      KDC(JJK),II,NBE(JT1),NBT(JT1),X1,X1,SDW(JJK)
    6 EE=EMX(JT1)
      IF(II==NHC(19))THEN
          DO IHD=1,NHRD(IOW)
              IF(IGZX(IHD,IOW)==ISA)EXIT
          END DO
          IF(IHD<=NHRD(IOW))EE=EE*GCOW(IHD,ISA)/WSA(ISA)
      END IF
      PPL0(JJK,ISA)=(1.-FPOP(JT1))*PPL0(JJK,ISA)
      XLAI(JJK,ISA)=FNPP(PPL0(JJK,ISA))
      DMLX(JJK)=XLAI(JJK,ISA)
      DMX=TLD(JT1)
	  CALL TMIX(EE,DMX,0,0)
      IF(DMX>BIG(ISA))TLD(JT1)=BIG(ISA)
      IF(II==NHC(22).OR.II==NHC(19))GO TO 26
	  IF(II==NHC(15))THEN
	      SATC(LID(2,ISA),ISA)=PRMT(39)
	  ELSE
	      IF(II==NHC(16))SATC(LID(2,ISA),ISA)=SATK(ISA)
      END IF
   57 IF(IDR(ISA)>0)THEN
          IF(II==NHC(25))THEN
              HCL(IDR(ISA),ISA)=HCLN(ISA)
          ELSE
              IF(II==NHC(26))HCL(IDR(ISA),ISA)=HCLD(ISA)
          END IF
      END IF
      IF(KFL(31)>0)WRITE(KW(31),50)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JT1),&
      KDC(JJK),II,NBE(JT1),NBT(JT1),CSTX,COX,FULU(JT1)
      SMFU(ISA)=SMFU(ISA)+FULU(JT1)
      SMST(ISA)=SMST(ISA)+STIR(JT1)
    7 XX=TLD(JT1)*1000.
      IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),28)ISA,NBSA(ISA),IYR,&
     &MO,KDA,TIL(JT1),XX,XHSM(ISA)
      IF(II/=NHC(17).AND.II/=NHC(18))GO TO 26
      IF(II/=NHC(18))THEN
          DHT(ISA)=DKH(JT1)
          DKHL(ISA)=DHT(ISA)
          DKIN(ISA)=DKI(JT1)
          IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),30)ISA,NBSA(ISA),IYR,&
          &MO,KDA,DHT(ISA),DKIN(ISA),XHSM(ISA)
          GO TO 26
      END IF
      DHT(ISA)=0.
      DKHL(ISA)=0.
      IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),30)ISA,NBSA(ISA),IYR,&
      MO,KDA,DHT(ISA),DKIN(ISA),XHSM(ISA)
      GO TO 26
   42 CALL TBURN
      GO TO 7
   63 ICV=1
      GO TO 57
   64 ICV=0
      GO TO 57
   43 KOMP(KT(ISA),ISA)=0
      ISCP(ISA)=ISCP(ISA)+1
      IF(ISCP(ISA)<MSCP)THEN
          JRT=1
          RETURN
      END IF
      ISCP(ISA)=0
      XX=1.-ORHI(JT1)
      X4=RSDM(LD1,ISA)*ORHI(JT1)
      RSDM(LD1,ISA)=RSDM(LD1,ISA)-X4
      X1=WCMU(LD1,ISA)*ORHI(JT1)
      WCMU(LD1,ISA)=WCMU(LD1,ISA)-X1
      X2=WCOU(LD1,ISA)*ORHI(JT1)
      WCOU(LD1,ISA)=WCOU(LD1,ISA)-X2
      WNOU(LD1,ISA)=WNOU(LD1,ISA)*XX      
      WNMU(LD1,ISA)=WNMU(LD1,ISA)*XX
      WLM(LD1,ISA)=WLM(LD1,ISA)*XX
      WLS(LD1,ISA)=WLS(LD1,ISA)*XX
	  WLSL(LD1,ISA)=WLSL(LD1,ISA)*XX
      X2=WLMC(LD1,ISA)*ORHI(JT1)
      WLMC(LD1,ISA)=WLMC(LD1,ISA)-X2
      X1=WLSC(LD1,ISA)*ORHI(JT1)
      WLSC(LD1,ISA)=WLSC(LD1,ISA)-X1
      SMM(101,MO,ISA)=SMM(101,MO,ISA)+X1+X2
      VAR(101,ISA)=X1+X2
      WLSLC(LD1,ISA)=WLSLC(LD1,ISA)*XX
      WLSLNC(LD1,ISA)=WLSLNC(LD1,ISA)*XX
      X2=ORHI(JT1)*WNMN(LD1,ISA)
      X3=ORHI(JT1)*WNMA(LD1,ISA)
      WNMN(LD1,ISA)=MAX(1.E-5,WNMN(LD1,ISA)-X2)
      WNMA(LD1,ISA)=WNMA(LD1,ISA)-X3
      X1=ORHI(JT1)*WLMN(LD1,ISA)
      X5=ORHI(JT1)*WLSN(LD1,ISA)
      WLMN(LD1,ISA)=WLMN(LD1,ISA)-X1
      WLSN(LD1,ISA)=WLSN(LD1,ISA)-X5
      SMM(89,MO,ISA)=SMM(89,MO,ISA)+X1+X2+X3+X5
      VAR(89,ISA)=X1+X2+X3+X5
      X3=WPOU(LD1,ISA)*ORHI(JT1)
      WPOU(LD1,ISA)=WPOU(LD1,ISA)-X3
      X1=ORHI(JT1)*FOP(LD1,ISA)
      FOP(LD1,ISA)=FOP(LD1,ISA)-X1
      X2=ORHI(JT1)*WPML(LD1,ISA)
      WPML(LD1,ISA)=WPML(LD1,ISA)-X2
      X5=WPMU(LD1,ISA)*ORHI(JT1)
      WPMU(LD1,ISA)=WPMU(LD1,ISA)-X5
      SMM(90,MO,ISA)=SMM(90,MO,ISA)+X1+X2+X3+X5
      VAR(90,ISA)=X1+X2+X3+X5
      SMNU(IDON(ISA))=SMNU(IDON(ISA))+X4*WSA(ISA)
      IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),44)ISA,NBSA(ISA),IYR,&
      MO,KDA,TIL(JT1),X4,ORHI(JT1),RSD(LD1,ISA),RSDM(LD1,ISA),SMNU(IDON(&
      ISA)),XHSM(ISA)
      RETURN
   10 CALL PESTF
      X6=PSTF(ISA)
      IF(II==NHC(19).AND.LGZ==0)THEN
          CALL TGRAZ(JRT)
          IF(JRT==0)GO TO 6
          RETURN
      END IF
      IF(ORHI(JT1)>0..OR.II==NHC(22))THEN
          CALL THVOR(X6,JRT)
	      IF(JRT==0)GO TO 6
	      RETURN
	  END IF
      DO K=1,NN
          IF(KGO(JE(K,ISA),ISA)==0)CYCLE
          IF(JH(IRO(ISA),KT(ISA),ISA)==KDC(JE(K,ISA)))GO TO 37
      END DO
      CSTX=0.
      COX=0.
      JRT=1
      RETURN 
   37 JJK=JE(K,ISA)
      IF(II==NHC(1))GO TO 22
      IF(JP(JJK,ISA)==0)THEN
          JP(JJK,ISA)=1
          NCR(JJK,ISA)=NCR(JJK,ISA)+1
      END IF
      HUF(JJK,ISA)=MAX(HUF(JJK,ISA),HU(JJK,ISA))
      DMF(JJK,ISA)=DM1(JJK,ISA)
      TRA(JJK,ISA)=SRA(JJK,ISA)+TRA(JJK,ISA)
      IF(RD(JJK,ISA)>RDF(JJK,ISA))RDF(JJK,ISA)=RD(JJK,ISA)
      XX=DM(JJK,ISA)+.001
      X2=UN1(JJK,ISA)/XX
      X7=X2*.001
      X3=UP1(JJK,ISA)/XX
      XX=STD(JJK,ISA)+1.E-10
      RNR=STDN(JJK,ISA)/XX
      RPR=STDP(JJK,ISA)/XX
      STDL(JJK,ISA)=CLG(ISA)*XX
      RLR=STDL(JJK,ISA)/XX
      XX=SWH(JJK,ISA)
      F=XX/(XX+EXP(SCRP(10,1)-SCRP(10,2)*XX))
      XX=MAX(AJHI(JJK,ISA)-WSYF(JJK),0.)
      FT=MAX(.1,1.+PRMT(81)*(IYR-2000))
      X1=MIN(F*XX+WSYF(JJK),.9*DM(JJK,ISA)/(STL(JJK,ISA)+1.E-10))*FT
      X1=MAX(X1,WSYF(JJK))
      X2=1000.*CNY(JJK)*(X7/BN(3,JJK))**.1
      X3=1000.*CPY(JJK)*(.001*X3/BP(3,JJK))**.1
      XZ=X1*STL(JJK,ISA)
      AJHI(JJK,ISA)=0.
      YZ=X1*STD(JJK,ISA)
      ZZ=MAX(.01,1.-X1)
      CPHT(JJK,ISA)=CPHT(JJK,ISA)*ZZ
      HU(JJK,ISA)=MAX(.1,HU(JJK,ISA)*ZZ)
      SLAI(JJK,ISA)=SLAI(JJK,ISA)*ZZ
      STL(JJK,ISA)=STL(JJK,ISA)*ZZ
      STD(JJK,ISA)=STD(JJK,ISA)*ZZ
      STDL(JJK,ISA)=STDL(JJK,ISA)*ZZ
      YLD(JJK)=XZ*HE(JT1)*X6
      YLSD=YZ*HE(JT1)
      Y4=YZ*RNR
      Y5=YZ*RPR
      STDN(JJK,ISA)=MAX(1.E-10,STDN(JJK,ISA)-Y4)
      STDP(JJK,ISA)=MAX(1.E-10,STDP(JJK,ISA)-Y5)
      STDL(JJK,ISA)=MAX(STDL(JJK,ISA)-YZ*RLR,.1*STD(JJK,ISA))
      X4=MIN(XZ*X2,UN1(JJK,ISA))
      X5=MIN(XZ*X3,UP1(JJK,ISA))
      X11=XZ-YLD(JJK)+YZ-YLSD
      Z2=YLSD*RNR
      Z3=YLSD*RPR
      YLN=MIN(.9*(UN1(JJK,ISA)+STDN(JJK,ISA)),YLD(JJK)*X2+Z2)
      YLP=MIN(.9*(UP1(JJK,ISA)+STDP(JJK,ISA)),YLD(JJK)*X3+Z3)
      X10=X4-YLN+Y4
      CALL NCNSTD(.05,X11,X10,LD1)
      FOP(LD1,ISA)=MAX(.01,FOP(LD1,ISA)+X5-YLP+Y5)
      YY=YLD(JJK)+YLSD
      YLD(JJK)=YY
	  IF(IDC(JJK)/=NDC(9))THEN
          YLD1(JJK,ISA)=YLD1(JJK,ISA)+YY
      ELSE
          YLD1(JJK,ISA)=YLD1(JJK,ISA)+FTO(JJK)*YY
          YLD2(JJK,ISA)=YLD2(JJK,ISA)+YLD1(JJK,ISA)*(1./FLT(JJK)-1.)
      END IF
      X2=DM(JJK,ISA)
      X3=RW(JJK,ISA)
      JD(ISA)=JJK
      SRA(JJK,ISA)=0.
      UN1(JJK,ISA)=UN1(JJK,ISA)-X4
      UP1(JJK,ISA)=UP1(JJK,ISA)-X5
      DM(JJK,ISA)=DM(JJK,ISA)-XZ
      YLNF(JJK,ISA)=YLNF(JJK,ISA)+YLN
      YLPF(JJK,ISA)=YLPF(JJK,ISA)+YLP
      TYN(ISA)=TYN(ISA)+YLN
      TYP(ISA)=TYP(ISA)+YLP
	  IF(ICUS(JT1)/=0.AND.CSTX<1.E-10)THEN
          CSTX=-CSTX*YLD(JJK)
          COX=CSTX
      END IF
      IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),29)ISA,NBSA(ISA),IYR,&
     &MO,KDA,TIL(JT1),CPNM(JD(ISA)),YY,X2,X3,X1,X6,X7,XHSM(ISA)
      GO TO 6
   22 IF(IPD==5)THEN
          CALL SPRNT
          WRITE(KW(1),31)
          CALL SOLIOP
          CALL SOLIOC
      END IF
      CALL TRDST
      X1=TRSD(ISA)+STD(JJK,ISA)+STDO(ISA)
      SRSD(ISA)=SRSD(ISA)+X1
      VIR(JJK,ISA)=VIRT(ISA)
      VIRT(ISA)=0.
      VIR0(ISA)=0.
      WS(ISA)=1.
      IGO(ISA)=MAX(0,IGO(ISA)-1)
      KGO(JJK,ISA)=0
	  JE(K,ISA)=12
      HU(JJK,ISA)=0.
      HUI(JJK,ISA)=0.
      HSM(ISA)=0.
      SLAI(JJK,ISA)=0.
      WLV(JJK,ISA)=0.
	  IYH(JJK,ISA)=0
      NII(ISA)=IRI(ISA)
      CSTF(JJK,ISA)=COST(ISA)
      COST(ISA)=0.
      IHU(JJK,ISA)=IHU(JJK,ISA)+1
      IF(IHU(JJK,ISA)>NHU(JJK,ISA))IHU(JJK,ISA)=1
      CAW(JJK,ISA)=AWC(JJK,ISA)
      ETG(JJK,ISA)=ACET(JJK,ISA)+ETG(JJK,ISA)
      PSTF(ISA)=0.
      IPST(ISA)=0
	  FGC(ISA)=0.
	  FGSL(ISA)=0.
      GO TO 6
   26 XX=1.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          BDU=1.6+.005*SAN(ISL,ISA)
          X2=100.*Z(ISL,ISA)
          IF(X2>20.)THEN
              X1=0.
          ELSE
              X1=EXP(-X2)
          END IF
          BDP(ISL,ISA)=BDP(ISL,ISA)+(BDU-BDP(ISL,ISA))*FRCP(JT1)*.5*(X1+XX)
          XX=X1
      END DO
      RETURN
   28 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,2X,'DPTH = ',F7.0,' mm',2X,'HUSC = ',F6.2)
   29 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,2X,A4,2X,'YLD=',F6.2,'t/ha',2X,&
      'BIOM=',F6.2,'t/ha',2X,'RW=',F5.2,'t/ha',2X,'HI=',F6.2,2X,'PSTF=',&
      F5.2,2X,'NCN=',F6.3,'G/G',2X,'HUSC=',F5.2)
   30 FORMAT(1X,2I8,1X,I4,2I2,2X,'DKH = ',F6.0,' mm',3X,'DKI = ',F7.2,&
      ' m',2X,'HUSC= ',F5.2)
   31 FORMAT(T5,'SOIL DATA')
   32 FORMAT(1X,2I8,1X,I4,2I2,2X,A4,2X,'RSD = ',F5.1,'t')
   44 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,2X,'MNU SCRP=',F8.2,'t',2X,'SCRP EF=&
      ',F6.3,2X,'RSD REMAIN=',E13.5,'t/ha',2X,'MNU REMAIN=',E13.5,'t/ha'&
      ,2X,'MNU STK PL=',E13.5,'t',2X,'HUSC=',F5.2)
   49 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,F10.2,10X,3F10.2)
   50 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
      END