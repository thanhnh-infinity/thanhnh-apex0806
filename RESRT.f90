      SUBROUTINE RESRT
!     APEX0806
!     THIS SUBPROGRAM ROUTES WATER AND SEDIMENT THROUGH RESEVOIRS.
!     COMPUTES EVAPORATION AND SEEPAGE FROM THE RESERVOIR.
      USE PARM
      IDO=IDOT(ICMD)
      IDN1=IDN1T(ICMD)
      IDN2=IDN2T(ICMD)
      IDRO(IDN2)=IDO
	  IDX=IDOA(IDN2)
	  IDNB(IDO)=NBSA(IDN2)
	  LD1=LID(1,IDN2)
	  RWSA(IDO)=RWSA(IDN1)
	  YB=RSYS(IDN2)
	  DEP=0.
	  Y2=0.
      Y2O=0.
	  IF(IEXT(IDN2)>0)THEN
	      A1=WSA(IDN2)
	  ELSE
	      A1=RWSA(IDN1)
	  END IF
      A10=10.*A1
      RSSF(IDO)=RSSF(IDN1)
      V0=RSV(IDN2)
      XX=V0/A10
      STV(11,MO,IDN2)=XX
      VARS(11)=XX
      TC(IDO)=TC(IDN1)
      SST(IDO)=SST(IDN1)
	  XX=MAX(.1,V0+RVP0(IDN2)-RSVP(IDN2)+RVE0(IDN2)-RSVE(IDN2))
      X2=BR1(IDN2)*XX**BR2(IDN2)
   	  RSSA(IDN2)=MIN(X2,RSAE(IDN2))
      STV(13,MO,IDN2)=RSSA(IDN2)
      VARS(13)=RSSA(IDN2)
      EV=10.*EO*RSSA(IDN2)
      SP=RSHC(IDN2)*RSSA(IDN2)*240.
      RFRA=RFV(IRF(IDN2))*RSSA(IDN2)*10.
	  SALA=10.*(WSA(IDN2)-RSSA(IDN2))
	  X1=A10*QVOL(IDN1)
	  X2=10.*RSSA(IDN2)*QVOL(IDX)
	  SMM(125,MO,IDN2)=SMM(125,MO,IDN2)+X1
	  SMM(126,MO,IDN2)=SMM(126,MO,IDN2)+X2
	  SMM(127,MO,IDN2)=SMM(127,MO,IDN2)+RFRA
	  Q1=X1-X2+RFRA
	  Y1=YSD(NDRV,IDN1)*A1
	  YU1=YMNU(IDN1)*A1
	  RSM1=RSDM(LD1,IDN2)
	  XX=SP+EV
      X1=Q1/A10
      SMM(64,MO,IDN2)=SMM(64,MO,IDN2)+Q1
      VAR(64,IDN2)=Q1
      V0=V0+Q1
      IF(V0<=XX)THEN
          X1=V0/(XX+1.E-10)
          EV=EV*X1
          SP=SP*X1
          RSV(IDN2)=0.
          QVOL(IDO)=0.
          YSD(NDRV,IDO)=0.
          YN(IDO)=0.
          YP(IDO)=0.
          QN(IDO)=0.
          QP(IDO)=0.
          TSFN(IDO)=0.
          RSFN(IDO)=0.
	      RSYS(IDN2)=0.
	      RSOP(IDN2)=.0001
          RSON(IDN2)=.0001
	      RSSP(IDN2)=.0001
	      RSO3(IDN2)=.0001
	      DEP=YB
	      GO TO 7
      END IF
      V0=V0-XX
      VRR=V0-RSVE(IDN2)
      OFLO=0.
	  OFP=0.
	  I1=1
	  A3=A1
      IF(VRR>0.)OFLO=VRR
      VVR=V0-RSVP(IDN2)
      IF(VVR>0.)THEN
	      X1=MIN(RSRR(IDN2),RSVE(IDN2)-RSVP(IDN2),VVR)
	      IF(ISAO(IDN2)==0)THEN
	          OFLO=OFLO+X1
	      ELSE
	          OFP=X1
	          I2=NISA(ISAO(IDN2))
	          I1=IDOA(I2)
	          A3=WSA(I2)
              QRP(I1)=QRP(I1)+.1*OFP/A3
 	          SMM(111,MO,IDN2)=SMM(111,MO,IDN2)+OFP
	          VAR(111,IDN2)=VAR(111,IDN2)+OFP
	          V0=V0-OFP
	      END IF    
          V0=V0-OFLO
      END IF
      VV=MAX(1.E-5,V0+OFLO+OFP)
      SMM(65,MO,IDN2)=SMM(65,MO,IDN2)+OFLO
      VAR(65,IDN2)=OFLO
	  QVOL(IDO)=OFLO/A10
	  YY=RSYS(IDN2)+Y1
	  IF(YY<1.E-10)GO TO 7
      SMM(68,MO,IDN2)=SMM(68,MO,IDN2)+Y1
      VAR(68,IDN2)=Y1
      CY=YY/VV
      CD=MAX(RSYN(IDN2),(CY-RSYN(IDN2))*RSDP(IDN2)+RSYN(IDN2))
	  X1=MIN(CD,CY)
      Y2=OFLO*X1
	  Y2O=MIN(YY,OFP*X1)
	  YSD(NDRV,I1)=YSD(NDRV,I1)+Y2O/A3
	  SMM(112,MO,IDN2)=SMM(112,MO,IDN2)+Y2O
	  VAR(112,IDN2)=Y2O
	  YSD(NDRV,IDO)=Y2/A1
	  RSM1=WSA(IDN2)*RSM1
	  YYU=RSM1+YU1
	  CUN=YYU/VV
	  CDU=CD*CUN/CY
	  X1=MIN(CDU,CUN)
	  YU2=MIN(OFLO*X1,.5*YU1+.1*RSM1)
	  YMNU(IDO)=MAX(0.,YU2)/A1
	  RSDM(LD1,IDN2)=(YYU-YU2)/WSA(IDN2)
	  SMM(69,MO,IDN2)=SMM(69,MO,IDN2)+Y2
      VAR(69,IDN2)=Y2
      DEP=MAX(0.,V0*(CY-CD))
      X1=DEP/A1
      SRCH(13,IDO)=SRCH(13,IDO)+X1
      SMM(70,MO,IDN2)=SMM(70,MO,IDN2)+DEP
      VAR(70,IDN2)=DEP
      X1=YY-Y2-Y2O-DEP
      RSYS(IDN2)=MAX(1.E-10,X1)
      XX=X1/A1
      STV(12,MO,IDN2)=X1
      VARS(12)=XX
      TOT=0.
      DRTO=MIN(1.,(.001+Y2)/(Y1+.001))
      B1=LOG(DRTO)/4.47
      DO I=1,NSZ
          PCTH(I,IDO)=PCT(I,IDN1)
          X1=MAX(-10.,B1*PSZY(I))
          PCT(I,IDO)=PCT(I,IDN1)*EXP(X1)
          TOT=TOT+PCT(I,IDO)
      END DO
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDO)/(TOT+1.E-10)
      END DO
      YX=YY+1.E-5
	  RSON(IDN2)=RSON(IDN2)+YN(IDN1)*A1
      CON=RSON(IDN2)/YX
	  RSOP(IDN2)=RSOP(IDN2)+YP(IDN1)*A1
      COP=RSOP(IDN2)/YX
	  X2=Y2+DEP+Y2O
      X1=CON*X2
      RTO=X1/RSON(IDN2)
      IF(RTO>.1)THEN
          X1=.1*RSON(IDN2)
          CON=X1/X2
      END IF
      YN(IDO)=CON*Y2/A1
      DRTP=MIN(1.,(.001+Y2O)/(Y1+.001))
	  X3=MIN(RSON(IDN2),CON*Y2O)
	  SMM(113,MO,IDN2)=SMM(113,MO,IDN2)+X3
      YN(I1)=YN(I1)+X3/A3
      RSON(IDN2)=MAX(.01,RSON(IDN2)-X1)
      X1=COP*X2
      RTO=X1/RSOP(IDN2)
      IF(RTO>.1)THEN
          X1=.1*RSOP(IDN2)
          COP=X1/X2
      END IF
      YP(IDO)=COP*Y2/A1
	  X3=MIN(RSOP(IDN2),COP*Y2O)
	  SMM(114,MO,IDN2)=SMM(114,MO,IDN2)+X3
      YP(I1)=YP(I1)+X3/A3
      RSOP(IDN2)=MAX(.01,RSOP(IDN2)-X1)
      RSO3(IDN2)=RSO3(IDN2)+A1*QN(IDN1)
      CO3=RSO3(IDN2)/VV
      RSSP(IDN2)=RSSP(IDN2)+A1*QP(IDN1)
      CSP=RSSP(IDN2)/VV
      QN(IDO)=OFLO*CO3
	  X3=OFP*CO3
	  SMM(115,MO,IDN2)=SMM(115,MO,IDN2)+X3
	  QN(I1)=QN(I1)+X3/A3
      RSO3(IDN2)=MAX(.01,RSO3(IDN2)-QN(IDO)-X3)
      QP(IDO)=OFLO*CSP
	  X3=OFP*CSP
	  SMM(116,MO,IDN2)=SMM(116,MO,IDN2)+X3
	  QP(I1)=QP(I1)+X3/A3
      RSSP(IDN2)=MAX(.01,RSSP(IDN2)-QP(IDO)-X3)
      RSV(IDN2)=V0
      QN(IDO)=QN(IDO)/A1
      QP(IDO)=QP(IDO)/A1
    7 SMM(66,MO,IDN2)=SMM(66,MO,IDN2)+EV
      VAR(66,IDN2)=EV
      TREV=TREV+EV
      SAET=SALA*AET
      TSAE=TSAE+SAET
      SMM(67,MO,IDN2)=SMM(67,MO,IDN2)+SP
      VAR(67,IDN2)=SP
      GWST(IDN2)=GWST(IDN2)+.1*SP/WSA(IDN2)
	  X1=DEP/RSBD(IDN2)
      RSVP(IDN2)=MAX(0.,RSVP(IDN2)-X1)
	  RSVE(IDN2)=MAX(0.,RSVE(IDN2)-X1)
	  DF=YB+Y1-Y2-Y2O-DEP-RSYS(IDN2)
	  IF(KFL(13)==0)RETURN
	  IF(ABS(DF)>1.E-5)WRITE(KW(13),13)DF
	  WRITE(KW(13),12)IDN2,NBSA(IDN2),IY,MO,KDA,RFV(IRF(IDN2)),Q1,EV,SP,&
      OFLO,RSV(IDN2),RSVP(IDN2),RSVE(IDN2),Y1,Y2,Y2O,DEP,RSYS(IDN2),RSSA&
      (IDN2)
      RETURN     
   12 FORMAT(1X,2I4,1X,I4,2I2,2X,8F10.0,20F10.2)
   13 FORMAT(1X,'!!!!!',E16.6)	 
	  END 
      