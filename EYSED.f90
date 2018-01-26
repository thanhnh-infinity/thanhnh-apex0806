      SUBROUTINE EYSED(JRT)
!     APEX0806
!     THIS SUBPROGRAM PREDICTS DAILY SOIL LOSS CAUSED BY WATER EROSION
!     AND ESTIMATES THE NUTRIENT ENRICHMENT RATIO.
      USE PARM
      JRT=0
      LD1=LID(1,ISA)      
      CALL EYCC
      XX=CVF(ISA)*USL(ISA)
      YLM=.5*QVOL(IDO)
      IF(RFV(IRF(ISA))>12.7)THEN
          YSD(3,IDO)=MIN(YLM,EI*XX*1.292)
          SMM(25,MO,ISA)=SMM(25,MO,ISA)+CVF(ISA)*EI
          SMM(24,MO,ISA)=SMM(24,MO,ISA)+EI
          VAR(24,ISA)=EI
          X1=EI*CVF(ISA)
          YSD(7,IDO)=MIN(YLM,X1*RSLK(ISA))
	      RXM=BETA/(1.+BETA)
	      RLFX=UPSX(ISA)**RXM
	      YI=MIN(YLM,.5*X1*RSK(ISA))
	      SUM=PSZX(1)*SAN(LD1,ISA)
	      SUM=SUM+PSZX(2)*SIL(LD1,ISA)
	      SUM=SUM+PSZX(3)*CLA(LD1,ISA)
	      SUM=.01*PRMT(65)*SUM/(QPR(IDO)+1.E-5)
	      T2=PRMT(66)*QPR(IDO)*STP(ISA)
	      IF(T2>YI)THEN
	          YSD(8,IDO)=1.5*X1*RLFX*RSK(ISA)
	      ELSE
	          YSD(8,IDO)=MAX(0.,YI+SUM*(T2-YI))
	      END IF
      END IF
      IF(QVOL(IDO)<1.)THEN
          JRT=1
          RETURN
      END IF
      REPC=REP*(QVOL(IDO)/(RFV(IRF(ISA))+.1))**.1
      YSD(2,IDO)=MIN(YLM,(.646*EI+.45*QVOL(IDO)*RQRB(IDO)**.3333)*XX)
      QQ=QVOL(IDO)*RQRB(IDO)
      X3=1.+WSA(ISA)
      YSD(4,IDO)=MIN(YLM,.79*X3**.009*QQ**.65*XX)
      X1=SQRT(QQ)
      YSD(1,IDO)=MIN(YLM,PRMT(33)*X1*XX)
      YSD(6,IDO)=MIN(YLM,BUS(1)*X3**BUS(4)*QVOL(IDO)**BUS(2)*RQRB(IDO)&
      **BUS(3)*XX)
      CX(MO,ISA)=CX(MO,ISA)+1.
      YSD(5,IDO)=MIN(YLM,1.586*X3**.12*QQ**.56*XX)
      LD1=LID(1,ISA)      
      IF(RSDM(LD1,ISA)>0.)THEN
          X2=PRMT(71)*AGPM(ISA)
          IF(X2<10.)THEN
              X2=EXP(-X2)
          ELSE
              X2=.00005
          END IF
          YMNU(IDO)=MIN(.9*RSDM(LD1,ISA),PRMT(62)*X1*SLF(ISA)*&
          PEC(ISA)*RSDM(LD1,ISA)**PRMT(68))
          YMNU(IDO)=MAX(0.,YMNU(IDO))
          IF(IDFH(ISA)>0)SYMU=SYMU+YMNU(IDO)*WSA(ISA)
      END IF
      IF(URBF(ISA)>0.)THEN
          DO I=1,8
              YSD(I,IDO)=YSD(I,IDO)*(1.-URBF(ISA))+URBF(ISA)*QURB(IDO)*.0003
          END DO
      END IF
      DRTO=MIN(.99,(RQRB(IDO)/REPC)**.56)
      DRAV(IDO)=DRAV(IDO)+DRTO
      B1=LOG(DRTO)/4.47
      SUM=0.
      DO I=1,NSZ
          PCTH(I,IDO)=PCT(I,ISA)
          X1=MAX(-10.,B1*PSZ(I))
          PCT(I,IDO)=MAX(.0001,PCT(I,ISA)*EXP(X1))
          SUM=SUM+PCT(I,IDO)
      END DO
      PSZM(IDO)=0.
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDO)/(SUM+1.E-10)
          PSZM(IDO)=PSZM(IDO)+PCT(I,IDO)*PSZY(I)
      END DO
      CY=.1*YSD(NDRV,IDO)/QVOL(IDO)
      IF(IERT==0)THEN
          B2=LOG10(DRTO)/2.699
          B1=1./.1**B2
          ERTO=B1*(CY+1.E-4)**B2
          ERTP=ERTO
      ELSE
          ERTO=PRMT(54)/(CY+1.E-4)**PRMT(55)
          ERTP=PRMT(57)/(CY+1.E-4)**PRMT(58)
      END IF
      RETURN
      END